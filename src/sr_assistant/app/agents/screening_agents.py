# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Screening prompts and agents/chains. Also Streamlit chat UI agent(s).

These are not agents, just chains.

Considering we're using OpenAI for the prototype, it doesn't make sense to use two
different models to emulate human reviewers. Instead, we'll use the same model but
with two different prompts to better simulate human reviewers.

One is "conservative" and the other is "comprehensive"..

The key differences are in:

1. How they interpret missing information
2. What constitutes sufficient evidence
3. Confidence scoring interpretation
4. Threshold for uncertainty

Both get the same criteria and be asked for the same structured output, but their
interpretation frameworks differ. This should give us meaningful differences in
screening decisions while maintaining high-quality assessment.

## Concurrency

As we run with Streamlit for the prototype, concurrency could be handled either with
a threadpool executor, TaskGroup, or a queue with, e.g., Supabase edge function workers.
For this initial prototype, we'll stay within Streamlit and leverage LangChain
runnables' parallel and batching capabilities, which use a threadpool executor
under the hood.

## TODO

- Either move types here or from here to core.types/schemas. Not both.
"""

from __future__ import annotations

import typing as t
import uuid
from copy import deepcopy
from datetime import datetime, timezone

import streamlit as st
from langchain_community.callbacks.manager import get_openai_callback
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnableConfig, RunnableParallel
from langchain_google_genai.chat_models import ChatGoogleGenerativeAI
from langchain_openai.chat_models import ChatOpenAI
from loguru import logger
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    model_validator,
)

import sr_assistant.app.utils as ut
from sr_assistant.app.config import get_settings
from sr_assistant.core import models, schemas
from sr_assistant.core.schemas import (
    ResolverOutputSchema,
    ScreeningResponse,
    ScreeningResult,
)
from sr_assistant.core.types import ScreeningStrategyType

if t.TYPE_CHECKING:
    from langchain_community.callbacks.openai_info import OpenAICallbackHandler
    from langchain_core.tracers.schemas import Run


# NOTE: IMPORTANT - the order of with_structured_output and with_retry matters.
#       with_structured_output must be before with_retry.
resolver_model = (
    ChatGoogleGenerativeAI(
        model="gemini-2.5-pro-preview-05-06",
        temperature=0,
        max_tokens=None,
        # thinking_budget=24576,  # TODO: This is a langchain-google-genai v2.1.4 featere, but we can't upgrade due to LangChain minor version upgrade breaking the way on_end listener works (its modifications to RunTree are not present in chain output in later versions. We will file a bug report for this as it's an undocumented breaking change.) For now we've pinned LangChain and Pydantic versions to known working versions. We may need to refactor screening logic later on, for now we stick with this setup. Gemini uses thinking by default, without a budget, how deeply it thinks is dependent on the prompt, so it must encourage deep analysis!)  # noqa: W505
        timeout=None,
        max_retries=5,
        api_key=get_settings().GOOGLE_API_KEY,
        convert_system_message_to_human=True,  # Gemini doesn't support system messages
    )
    .with_structured_output(ResolverOutputSchema)
    .with_retry(
        stop_after_attempt=5,
        wait_exponential_jitter=True,
        retry_if_exception_type=(Exception,),
    )
)


RESOLVER_SYSTEM_PROMPT = """You are an expert systematic review screening assistant, specializing in conflict resolution. You are an analytical, critical, and deep thinker.

Your primary function is to resolve discrepancies in screening decisions between two independent reviewers (one 'conservative', one 'comprehensive') or to provide a definitive decision when both reviewers are 'uncertain'. Your goal is to arrive at a final, well-reasoned screening decision for a given research article based on the provided context.

To achieve this, you must think step-by-step. Analyze all provided information meticulously. Consider the perspectives and rationales of both reviewers, delving into their decision-making processes by examining their provided decision, confidence, rationale, extracted quotes, and any exclusion reasons. Explain your reasoning comprehensively before stating your final decision and confidence. Show your work clearly.

Your final output must conform to the specified structured format. Ensure your reasoning is thorough and directly supports your final decision.
"""

RESOLVER_HUMAN_PROMPT_TEMPLATE = """Please resolve the screening decision for the following research article based on the following detailed information.

<search_result_context>
    <title>{article_title}</title>
    <abstract>
{article_abstract}
    </abstract>
    <journal>{article_journal}</journal>
    <year>{article_year}</year>
    <source_id>{article_source_id}</source_id>
    <keywords>{article_keywords}</keywords>
    <mesh_terms>{article_mesh_terms}</mesh_terms>
</search_result_context>

<systematic_review_protocol>
    <background>
{review_background}
    </background>
    <research_question>
{review_research_question}
    </research_question>
    <criteria_framework>{review_criteria_framework}</criteria_framework>
    <inclusion_criteria>
{review_inclusion_criteria}
    </inclusion_criteria>
    <exclusion_criteria>
{review_exclusion_criteria}
    </exclusion_criteria>
</systematic_review_protocol>

<conservative_reviewer_assessment>
    <decision>{conservative_reviewer_decision}</decision>
    <confidence_score>{conservative_reviewer_confidence}</confidence_score>
    <rationale>
{conservative_reviewer_rationale}
    </rationale>
    <extracted_quotes>
{conservative_reviewer_quotes}
    </extracted_quotes>
    <exclusion_reasons>
{conservative_reviewer_exclusion_reasons}
    </exclusion_reasons>
</conservative_reviewer_assessment>

<comprehensive_reviewer_assessment>
    <decision>{comprehensive_reviewer_decision}</decision>
    <confidence_score>{comprehensive_reviewer_confidence}</confidence_score>
    <rationale>
{comprehensive_reviewer_rationale}
    </rationale>
    <extracted_quotes>
{comprehensive_reviewer_quotes}
    </extracted_quotes>
    <exclusion_reasons>
{comprehensive_reviewer_exclusion_reasons}
    </exclusion_reasons>
</comprehensive_reviewer_assessment>

<instructions_for_resolver>
Based on all the information above, please provide your final resolution.

First, articulate your detailed thought process within <thinking>...</thinking> XML tags. This thinking block is for your internal reasoning and will be reviewed for quality, but it is not part of the final structured output fields like 'resolver_reasoning'. Explain how you evaluated the reviewers' inputs against the protocol and the search result. Analyze any disagreements or uncertainties in depth.

After your thinking process, provide the structured output including:
- `resolver_decision`: Your final decision (e.g., 'include', 'exclude', 'uncertain').
- `resolver_reasoning`: A concise summary of your detailed thought process that justifies your decision. This should be a self-contained explanation.
- `resolver_confidence_score`: Your confidence in this final decision (0.0 to 1.0).
- `contributing_strategies`: A list of original screening strategies (e.g., ['conservative', 'comprehensive']) that you found to have contributed meaningfully to your final decision or that your decision aligns with. This can be an empty list if neither was particularly influential or if your reasoning diverged significantly.
</instructions_for_resolver>
"""

resolver_prompt = ChatPromptTemplate.from_messages(
    [
        ("system", RESOLVER_SYSTEM_PROMPT),
        ("human", RESOLVER_HUMAN_PROMPT_TEMPLATE),
    ]
)


class ResolverPromptInput(BaseModel):
    model_config = ConfigDict(frozen=True)

    article_title: str
    article_abstract: str
    article_journal: str | None
    article_year: str | None
    article_source_id: str
    article_keywords: str
    article_mesh_terms: str

    review_background: str | None
    review_research_question: str
    review_criteria_framework: str | None
    review_inclusion_criteria: str | None
    review_exclusion_criteria: str

    conservative_reviewer_decision: str
    conservative_reviewer_confidence: float
    conservative_reviewer_rationale: str
    conservative_reviewer_quotes: str
    conservative_reviewer_exclusion_reasons: str

    comprehensive_reviewer_decision: str
    comprehensive_reviewer_confidence: float
    comprehensive_reviewer_rationale: str
    comprehensive_reviewer_quotes: str
    comprehensive_reviewer_exclusion_reasons: str


def _format_extracted_quotes(quotes: list[str] | None) -> str:
    if not quotes:
        return "N/A"
    return "\n".join(f"- {q}" for q in quotes)


def _format_exclusion_reasons(reasons: schemas.ExclusionReasons | None) -> str:
    if not reasons:
        return "N/A"

    formatted_parts = []
    if reasons.population_exclusion_reasons:
        formatted_parts.append(
            f"Population: {', '.join(reason for reason in reasons.population_exclusion_reasons)}"
        )
    if reasons.intervention_exclusion_reasons:
        formatted_parts.append(
            f"Intervention: {', '.join(reason for reason in reasons.intervention_exclusion_reasons)}"
        )
    if reasons.comparison_exclusion_reasons:
        formatted_parts.append(
            f"Comparison: {', '.join(reason for reason in reasons.comparison_exclusion_reasons)}"
        )
    if reasons.outcome_exclusion_reasons:
        formatted_parts.append(
            f"Outcome: {', '.join(reason for reason in reasons.outcome_exclusion_reasons)}"
        )
    if reasons.reporting_exclusion_reasons:
        formatted_parts.append(
            f"Reporting: {', '.join(reason for reason in reasons.reporting_exclusion_reasons)}"
        )
    if reasons.study_design_exclusion_reasons:
        formatted_parts.append(
            f"Study Design: {', '.join(reason for reason in reasons.study_design_exclusion_reasons)}"
        )

    if not formatted_parts:
        return "N/A"
    return "\n".join(formatted_parts)


def _format_raw_data_list(raw_data_list: t.Any) -> str:
    if isinstance(raw_data_list, list):
        return ", ".join(str(item) for item in raw_data_list)
    if isinstance(raw_data_list, str):
        return raw_data_list
    return "N/A"


def prepare_resolver_inputs_for_prompt(
    search_result: models.SearchResult,
    review: models.SystematicReview,
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
) -> ResolverPromptInput:
    """Prepares the input dictionary for the resolver_prompt.

    Formats raw data from SearchResult, SystematicReview, and ScreeningResult
    instances into a flat dictionary suitable for the resolver prompt template,
    including formatting lists and structured objects into strings.
    """
    return ResolverPromptInput(
        article_title=search_result.title or "N/A",
        article_abstract=search_result.abstract or "N/A",
        article_journal=search_result.journal or "N/A",
        article_year=search_result.year or "N/A",
        article_source_id=search_result.source_id or "N/A",
        article_keywords=", ".join(search_result.keywords)
        if search_result.keywords
        else "N/A",
        article_mesh_terms=_format_raw_data_list(
            search_result.raw_data.get("mesh_terms")
            or search_result.raw_data.get("MeshHeadings")
        ),
        review_background=review.background or "N/A",
        review_research_question=review.research_question or "N/A",
        review_criteria_framework=review.criteria_framework.value
        if review.criteria_framework
        else "N/A",
        review_inclusion_criteria=review.inclusion_criteria or "N/A",
        review_exclusion_criteria=review.exclusion_criteria or "N/A",
        conservative_reviewer_decision=conservative_result.decision.value,
        conservative_reviewer_confidence=conservative_result.confidence_score,
        conservative_reviewer_rationale=conservative_result.rationale or "N/A",
        conservative_reviewer_quotes=_format_extracted_quotes(
            conservative_result.extracted_quotes
        ),
        conservative_reviewer_exclusion_reasons=_format_exclusion_reasons(
            conservative_result.exclusion_reason_categories
        ),
        comprehensive_reviewer_decision=comprehensive_result.decision.value,
        comprehensive_reviewer_confidence=comprehensive_result.confidence_score,
        comprehensive_reviewer_rationale=comprehensive_result.rationale or "N/A",
        comprehensive_reviewer_quotes=_format_extracted_quotes(
            comprehensive_result.extracted_quotes
        ),
        comprehensive_reviewer_exclusion_reasons=_format_exclusion_reasons(
            comprehensive_result.exclusion_reason_categories
        ),
    )


class ResolverChainInputDict(t.TypedDict):
    """Represents resolver chain prompt template input variables.

    Attributes:
        article_title (str): Title of the related article.
        article_abstract (str): Abstract of the related article.
        article_journal (str): Journal of the related article.
        article_year (str): Year of the related article.
        article_mesh_terms (list[str]): Mesh terms used in the search query. Set to "N/A" if not available.
        article_keywords (list[str]): Keywords used in the search query. Set to "N/A" if not available.
        review_research_question (str): Research question from the systematic review protocol (`SystematicReview.research_question`).
        review_inclusion_criteria (str): Inclusion criteria from the systematic review protocol (`SystematicReview.criteria_framework_answers).
        review_exclusion_criteria (str): Exclusion criteria from the systematic review protocol (`SystematicReview.exclusion_criteria`).
        review_background (str): Background from the systematic review protocol (`SystematicReview.background`).
        comprehensive_reviewer_decision (ScreeningResult): Comprehensive screening result to be used for context.
        comprehensive_reviewer_reasoning (str): Reasoning for the comprehensive screening decision.
        comprehensive_reviewer_confidence (float): Confidence score for the comprehensive screening decision. 1.0 is 100% confidence.
        conservative_reviewer_decision (ScreeningResult): Conservative reviewr model's decision: `include` | `exclude` | `uncertain` (`Scree)
        conservative_reviewer_reasoning (str): Reasoning for the conservative screening decision.
        conservative_reviewer_confidence (float): Confidence score for the conservative screening decision.
        current_date (str): Current date in ISO 8601 format. Partially bound to template, datetime.now(timezone.utc).
    """


class ResolverChainBatchInputDict(t.TypedDict):
    """ressolve_screening_result_batch.batch(**this_input_dict).

    Attributes:
        inputs (list[ResolverChainInputDict])
        config (list[RunnableConfig])
    """

    inputs: list[ResolverChainInputDict]
    config: list[RunnableConfig]


# TBD: class ResolverError (BaseModel or Exception?)


class ResolverChainOutputTuple(t.NamedTuple):
    """Tuple of search result and its resolver chain output.

    Attributes:
        search_result (SearchResult): Search result associated with this resolver
        comprehensive_result_id (uuid.UUID): ID of the comprehensive screening result
        conservative_result_id (uuid.UUID): ID of the conservative screening result
        resolver_result (ResolverChainOutput | ResolverError): Resolver result for conservative strategy
    """

    search_result: models.SearchResult


# TODO: This wsa simply copied from abstract screening chain, must be refactored for resolver chain
#       The reason for the listener is two-fold: 1) we don't want the model to see fields not related to the task in the output JSONSSchema (function call schema),
#       2) We want to grab the LangSmith trace_id and run id and other metadata from the LangSmith RunTree (imported from langchain as `Run`)
# @logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
# def resolver_chain_on_end_cb(run_obj: Run) -> None:
#     """Listener for resolver chain run on_end hook.


#     See `.with_listeners <https://python.langchain.com/api_reference/_modules/langchain_core/runnables/base.html#RunnableBinding.with_listeners>`_

#     Modifications to the `Run`|`RunTree` persist and no need to return anything. It's a
#
#     Reads the ScreeningResolutionSchema from the LLM structured output (whihch is handled by  LangChain)
#     and maps it to (to be created) ResolverChainResult BaseModel (name TBD)

#     Populates the model fields:
#     - start_time from `run_obj.start_time`
#     - end_time from `run_obj.end_time`
#     - trace_id from `run_obj.trace_id`
#     - tags from `run_obj.tags`
#     - review_id from `run_obj.metadata.get("review_id")`
#     - search_result_id from `run_obj.metadata.get("search_result_id")`
#     - comprehensive_result from `run_obj.metadata.get("comprehensive_result_id")` (db lookup)
#     - conservative_result from `run_obj.metadata.get("conservative_result_id")` (db lookup)
#     - resolver_result from `run_obj.metadata.get("resolver_result")`
#     - model_name from `run_obj` either from API response or if passed in during invocation then from run_obj.metadata
#     - response_metadata from `run_obj.metadata.get("response_metadata")`
#     """
#     cr_logger = logger.bind(run_id=run_obj.id, run_name=run_obj.name)
#     cr_logger.debug(f"screen_abstracts_chain_on_end_cb got RunTree: {run_obj!r}")
#     if not run_obj.outputs:
#         cr_logger.warning(f"No outputs for run: {run_obj!r}")
#         return
#     # Each child run is one screening chain invocation within RunnableParallel. This
#     # parent RunTree is the RunnableParallel itself.
#     for cr in run_obj.child_runs:
#         cr_logger = logger.bind(
#             run_id=cr.id,
#             run_name=cr.name,
#             parent_run_id=run_obj.id,
#             parent_run_name=run_obj.name,
#         )
#         cr_logger.debug(f"screen_abstracts_chain_on_end_cb got child run: {cr!r}")
#         if cr.outputs is None:
#             cr_logger.warning(f"No outputs for child run: {cr!r}")
#             continue
#         # Can't pass the strategy in the invocation metadata as its parallel.
#         # So run tags populated by LC/LS is the only way to get it. `map:key:...` is
#         # a kwarg key in RunnableParallel input dict.
#         output_key: ScreeningStrategyType
#         if isinstance(cr.tags, list):
#             maybe_output_key = next(
#                 (
#                     t.removeprefix("map:key:")
#                     for t in cr.tags
#                     if t.startswith("map:key:")
#                 ),
#                 None,
#             )
#             if maybe_output_key in ScreeningStrategyType:
#                 output_key = ScreeningStrategyType(maybe_output_key)
#                 cr_logger.debug(f"Output key for run: {cr.id!r} is {output_key}")
#             else:
#                 cr_logger.error("Unknown output key for run")
#                 continue
#         else:
#             cr_logger.error(f"No tags for run: {cr.id!r}, tags: {cr.tags!r}")
#             continue
#         orig_resp_model = cr.outputs.get("output")
#         if not isinstance(orig_resp_model, ScreeningResponse):
#             # We simply return the output as is and let caller make it a ScreeningError
#             cr_logger.error(f"Unknown response type: {orig_resp_model!r}")
#             continue
#         review_id = cr.metadata.get("review_id")
#         if not review_id:
#             cr_logger.error(f"Missing review_id in metadata: {cr.metadata!r}")
#             continue
#         review_id = uuid.UUID(review_id)
#         search_result_id = cr.metadata.get("search_result_id")
#         if not search_result_id:
#             cr_logger.error(f"Missing search_result_id in metadata: {cr.metadata!r}")
#             continue
#         search_result_id = uuid.UUID(search_result_id)
#         result_id = cr.id
#         trace_id = cr.trace_id  # this is shared between invocations with same input
#         inputs = cr.inputs
#         start_time = cr.start_time
#         end_time = cr.end_time
#         if not end_time:
#             cr_logger.warning("No end time for run: {cr!r}, setting to now")
#             end_time = datetime.now(tz=timezone.utc)
#         resp_metadata = {}  # FIXME: these don't work
#         invocation_params = {}
#         for ccr in cr.child_runs:
#             # TODO: parse model_id from openai metadata
#             if ccr.run_type == "chat_model":
#                 resp_metadata = ccr.extra.get("metadata", {})
#                 if not resp_metadata:
#                     cr_logger.warning(
#                         "No response metadata for chat_model child child run: {ccr!r}",
#                         ccr=ccr,
#                     )
#                 invocation_params = ccr.extra.get("invocation_params", {})
#                 if not invocation_params:
#                     cr_logger.warning(
#                         "No invocation_params for chat_model child child run:  {ccr!r}",
#                         ccr=ccr,
#                     )
#                 cr_logger.debug(f"chat_model child child run: {ccr!r}")
#         model_name = resp_metadata.get(
#             "ls_model_name", invocation_params.get("model_name", "not_found")
#         )
#         screening_result = ScreeningResult(
#             id=result_id,
#             trace_id=trace_id,  # this is shared between invocations with same input
#             review_id=review_id,
#             search_result_id=search_result_id,
#             start_time=start_time,
#             end_time=end_time,
#             model_name=model_name,
#             screening_strategy=output_key,
#             response_metadata=dict(
#                 inputs=inputs,
#                 run_name=run_obj.name,
#                 **resp_metadata,
#                 **invocation_params,
#             ),
#             **orig_resp_model.model_dump(),
#         )
#         cr.metadata["screening_strategy"] = output_key

#         strategy_tag = f"sra:screening_strategy:{output_key}"
#         cr.tags.append(strategy_tag)
#         if run_obj.tags and strategy_tag not in run_obj.tags:
#             run_obj.tags.append(strategy_tag)

#         if output_key not in run_obj.outputs:
#             logger.warning(
#                 f"Missing output key: {output_key!r} in {run_obj.outputs!r}, assigning as new key ..."
#             )
#         run_obj.outputs[output_key] = screening_result


class ScreenAbstractsChainOutputDict(t.TypedDict):
    """Final output from screen_abstracts_chain RunnableParallel call.

    Attributes:
        conservative (ScreeningResult|ScreeningResponse|t.Any): The conservative screening
            strategy result. Should be a `ScreeningResult`
            but may be `ScreeningResponse` if ``on_end`` listener failed parsing
            :class:`langsmith.RunTree`s. Or `typing.Any` depending on listener input
            RunTree. Caller should handle this.
        comprehensive (ScreeningResult|ScreeningResponse|t.Any): Same as above, but for
            comprehensive strategy.
    """

    conservative: ScreeningResult | ScreeningResponse | t.Any
    comprehensive: ScreeningResult | ScreeningResponse | t.Any


class ScreenAbstractsChainInput(BaseModel):
    """screen_abstracts_chain() prompt input variables.

    Attributes:
        background (str | None), optional: Background of the systematic review.
        research_question (str): Research question of the systematic review.
        inclusion_criteria (str): Inclusion criteria of the systematic review.
        exclusion_criteria (str): Exclusion criteria of the systematic review.
        title (str): Title of the PubMed search result.
        journal (str): Journal of the PubMed search result.
        year (str): Year of the PubMed search result.
        abstract (str): Abstract of the PubMed search result.

    Todo:
        - Use prompt .get_input_schema():
            ``ScreenAbstractsChainInput = prompt.get_input_schema()``
            or sublass as we want config (can it be passed?)

    Note:
        Update when prompt changes
    """

    model_config = ConfigDict(
        from_attributes=True,
        validate_assignment=True,
        extra="ignore",
    )
    # Protocol
    background: str | None = Field(default="<no background provided>")
    research_question: str
    inclusion_criteria: str
    exclusion_criteria: str
    # Search result
    title: str
    journal: str
    year: str
    abstract: str


class ScreenAbstractsChainInputDict(t.TypedDict):
    """Represents ScreenAbstractsChainInput.model_dump(). For type checking reasons.

    Attributes:
        background (str | None), optional: Background of the systematic review.
        research_question (str): Research question of the systematic review.
        inclusion_criteria (str): Inclusion criteria of the systematic review.
        exclusion_criteria (str): Exclusion criteria of the systematic review.
        title (str): Title of the PubMed search result.
        journal (str): Journal of the PubMed search result.
        year (str): Year of the PubMed search result.
        abstract (str): Abstract of the PubMed search result.

    Todo:
        - Generate from schema

    Note:
        Update when prompt changes

    """

    background: t.NotRequired[str]
    research_question: str
    inclusion_criteria: str
    exclusion_criteria: str
    title: str
    journal: str
    year: str
    abstract: str


screen_abstracts_chain_input_dict_adapter = TypeAdapter(ScreenAbstractsChainInputDict)
"""TODO: .validate_python() in screening wrapper"""


class ScreenAbstractsChainBatchInputDict(t.TypedDict):
    """chain.batch(**this_input).

    Attributes:
        inputs (list[ScreenAbstractsChainInput])
        config (list[RunnableConfig])
    """

    inputs: list[ScreenAbstractsChainInputDict]
    config: list[RunnableConfig]


class ScreeningError(BaseModel):
    """Screening error.

    Attributes:
        search_result (SearchResult): Search result associated with this screening error
        error (t.Any): Error encountered during screening
        message (str | None): Human-readable error message. Default is None.
    """

    search_result: models.SearchResult = Field(
        description="Search result associated with this screening error"
    )
    error: t.Any = Field(description="Error encountered during screening")
    message: str | None = Field(
        default=None, description="Human-readable error message"
    )

    @model_validator(mode="after")
    def _validate_error_and_message(self) -> t.Self:  # this is an instance method
        if isinstance(self.error, ScreeningResult) and self.message is None:
            msg = "ScreeningResult should always have a message"
            logger.warning(msg)
            self.message = msg
        if isinstance(self.error, ScreeningResponse) and self.message is None:
            msg = "Chain on_end failed, this is nested chain's screening output"
            logger.warning(msg)
            self.message = msg
        return self


class ScreenAbstractResultTuple(t.NamedTuple):
    """Tuple of search result and its screening results.

    Attributes:
        search_result (SearchResult): Search result associated with this screening
        conservative_result (ScreeningResult | ScreeningError): Screening result for conservative strategy
        comprehensive_result (ScreeningResult | ScreeningError): Screening result for comprehensive strategy
    """

    search_result: models.SearchResult
    conservative_result: ScreeningResult | ScreeningError
    comprehensive_result: ScreeningResult | ScreeningError


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def chain_on_error_listener_cb(run_obj: Run) -> None:
    # Get the error from the run object
    error = run_obj.error
    cb_logger = logger.bind(run_id=run_obj.id, run_name=run_obj.name)

    if error:
        # Log the exception if an error exists
        cb_logger.exception(f"Error during chain execution: {error!r}", error=error)
    else:
        # Log a general error message if run_obj.error is None for some reason
        cb_logger.error(
            "Chain execution failed, but error details not found in run object."
        )

    # Log the associated run object details
    cb_logger.error("Associated run object: {run_obj!r}", run_obj=run_obj)
    for cr in run_obj.child_runs:
        cb_cr_logger = logger.bind(
            run_id=cr.id,
            run_name=cr.name,
            parent_run_id=run_obj.id,
            parent_run_name=run_obj.name,
        )
        # Log child run errors if they exist
        child_error = cr.error
        if child_error:
            cb_cr_logger.exception(
                f"Associated child run failed: {child_error!r}",
                error=child_error,
                cr=cr,
            )
        else:
            cb_cr_logger.error("Associated child run details: {cr!r}", cr=cr)


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def screen_abstracts_chain_on_end_cb(run_obj: Run) -> None:  # noqa: C901
    """Listener for screening batch run on_end hook.

    See `.with_listeners <https://python.langchain.com/api_reference/_modules/langchain_core/runnables/base.html#RunnableBinding.with_listeners>`_

    Modifications to the `Run`|`RunTree` persist and no need to return anything. It's a
    Pydantic v1 model.

    Todo:
        - More metadata extraction
    """
    cr_logger = logger.bind(run_id=run_obj.id, run_name=run_obj.name)
    cr_logger.debug(f"screen_abstracts_chain_on_end_cb got RunTree: {run_obj!r}")
    if not run_obj.outputs:
        cr_logger.warning(f"No outputs for run: {run_obj!r}")
        return
    # Each child run is one screening chain invocation within RunnableParallel. This
    # parent RunTree is the RunnableParallel itself.
    for cr in run_obj.child_runs:
        cr_logger = logger.bind(
            run_id=cr.id,
            run_name=cr.name,
            parent_run_id=run_obj.id,
            parent_run_name=run_obj.name,
        )
        cr_logger.debug(f"screen_abstracts_chain_on_end_cb got child run: {cr!r}")
        if cr.outputs is None:
            cr_logger.warning(f"No outputs for child run: {cr!r}")
            continue
        # Can't pass the strategy in the invocation metadata as its parallel.
        # So run tags populated by LC/LS is the only way to get it. `map:key:...` is
        # a kwarg key in RunnableParallel input dict.
        output_key: ScreeningStrategyType
        if isinstance(cr.tags, list):
            maybe_output_key = next(
                (
                    t.removeprefix("map:key:")
                    for t in cr.tags
                    if t.startswith("map:key:")
                ),
                None,
            )
            if maybe_output_key in ScreeningStrategyType:
                output_key = ScreeningStrategyType(maybe_output_key)
                cr_logger.debug(f"Output key for run: {cr.id!r} is {output_key}")
            else:
                cr_logger.error("Unknown output key for run")
                continue
        else:
            cr_logger.error(f"No tags for run: {cr.id!r}, tags: {cr.tags!r}")
            continue
        orig_resp_model = cr.outputs.get("output")
        if not isinstance(orig_resp_model, ScreeningResponse):
            # We simply return the output as is and let caller make it a ScreeningError
            cr_logger.error(f"Unknown response type: {orig_resp_model!r}")
            continue
        review_id = cr.metadata.get("review_id")
        if not review_id:
            cr_logger.error(f"Missing review_id in metadata: {cr.metadata!r}")
            continue
        review_id = uuid.UUID(review_id)
        search_result_id = cr.metadata.get("search_result_id")
        if not search_result_id:
            cr_logger.error(f"Missing search_result_id in metadata: {cr.metadata!r}")
            continue
        search_result_id = uuid.UUID(search_result_id)
        result_id = cr.id
        trace_id = cr.trace_id  # this is shared between invocations with same input
        inputs = cr.inputs
        start_time = cr.start_time
        end_time = cr.end_time
        if not end_time:
            cr_logger.warning("No end time for run: {cr!r}, setting to now")
            end_time = datetime.now(tz=timezone.utc)
        resp_metadata = {}  # FIXME: these don't work
        invocation_params = {}
        for ccr in cr.child_runs:
            # TODO: parse model_id from openai metadata
            if ccr.run_type == "chat_model":
                resp_metadata = ccr.extra.get("metadata", {})
                if not resp_metadata:
                    cr_logger.warning(
                        "No response metadata for chat_model child child run: {ccr!r}",
                        ccr=ccr,
                    )
                invocation_params = ccr.extra.get("invocation_params", {})
                if not invocation_params:
                    cr_logger.warning(
                        "No invocation_params for chat_model child child run:  {ccr!r}",
                        ccr=ccr,
                    )
                cr_logger.debug(f"chat_model child child run: {ccr!r}")
        model_name = resp_metadata.get(
            "ls_model_name", invocation_params.get("model_name", "not_found")
        )
        screening_result = ScreeningResult(
            id=result_id,
            trace_id=trace_id,  # this is shared between invocations with same input
            review_id=review_id,
            search_result_id=search_result_id,
            start_time=start_time,
            end_time=end_time,
            model_name=model_name,
            screening_strategy=output_key,
            response_metadata=dict(
                inputs=inputs,
                run_name=run_obj.name,
                **resp_metadata,
                **invocation_params,
            ),
            **orig_resp_model.model_dump(),
        )
        cr.metadata["screening_strategy"] = output_key

        strategy_tag = f"sra:screening_strategy:{output_key}"
        cr.tags.append(strategy_tag)
        if run_obj.tags and strategy_tag not in run_obj.tags:
            run_obj.tags.append(strategy_tag)

        if output_key not in run_obj.outputs:
            logger.warning(
                f"Missing output key: {output_key!r} in {run_obj.outputs!r}, assigning as new key ..."
            )
        run_obj.outputs[output_key] = screening_result


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def make_screen_abstracts_chain_input(
    search_results_batch: list[models.SearchResult],
    review: models.SystematicReview,
) -> ScreenAbstractsChainBatchInputDict:
    """Create input for a batch to be passed to abstracts screening chain.

    chain.batch(**input) input as returned from this function.

    TODO: We pass the same input and config to both reviewers, is it possible
          to bind configs to the sub chains such that metadata from the parallel
          chain doesn't overwrite anything?

          Currently reviewer/strategy must be looked up in the listner as it's
          populated by the chain in the RunTree.outputs keys and RunTree.tags.

    Args:
        search_results_batch (list[models.SearchResult]): list of search results to be screened
        review (models.SystematicReview): systematic review associated with these search results

    Returns:
        ScreenAbstractsChainBatchInputDict: dict of list of prompt inputs and list of RunnableConfigs.
            Can be passed directly to chain.batch(**this)
    """
    logger.info(f"Creating batch input for {len(search_results_batch)} search results")
    chain_inputs: list[
        ScreenAbstractsChainInputDict
    ] = []  # input is same regardless of strategy, same prompt variables
    configs: list[RunnableConfig] = []
    for i, search_result in enumerate(search_results_batch):
        configs.append(
            RunnableConfig(
                run_name="screen_abstracts_chain",
                metadata={
                    "review_id": str(review.id),
                    "search_result_id": str(search_result.id),
                },
                tags=[
                    "sra:prototype",
                    "sra:streamlit",
                    "sra:step3",
                    f"sra:review_id:{review.id}",
                    f"sra:search_result_id:{search_result.id}",
                    f"sra:screen_abstracts_chain:i:{i}",
                ],  # sra:screening_strategy:comprehensive|conservative # added by listener
            )
        )
        chain_inputs.append(
            ScreenAbstractsChainInputDict(
                **ScreenAbstractsChainInput(
                    background=review.background,
                    research_question=review.research_question,
                    inclusion_criteria=review.inclusion_criteria or "",
                    exclusion_criteria=review.exclusion_criteria or "",
                    title=search_result.title or "",
                    journal=search_result.journal or "",
                    year=search_result.year or "",
                    abstract=search_result.abstract or "",
                ).model_dump()
            )
        )
        logger.info(
            f"Created chain input for search result {i}: config: {configs[-1]!r}, input: {chain_inputs[-1]!r}"
        )
    return ScreenAbstractsChainBatchInputDict(inputs=chain_inputs, config=configs)


# Prompts for screening

conservative_reviewer_prompt_text = """\
You are a highly conservative systematic reviewer, focusing on methodological \
rigor and strict interpretation of inclusion criteria. Your priority is to avoid \
including studies that might not fully meet the criteria.

When assessing abstracts:
- Require explicit statements matching criteria
- Flag any methodological ambiguity
- Consider unclear reporting as potential exclusion
- Demand high specificity in study characteristics
- Interpret missing information as a reason for uncertainty

Confidence Scoring for Conservative Review:
- Score < 0.7: Mark as "uncertain" when any required information is implicit or missing
- Score 0.7-0.85: Clear criteria match but some details could be more explicit
- Score > 0.85: Only when all criteria are explicitly and unambiguously met"""

comprehensive_reviewer_prompt_text = """\
You are a comprehensive systematic reviewer, focusing on potential relevance and \
broader interpretation of inclusion criteria. Your priority is to avoid excluding \
potentially relevant studies.

When assessing abstracts:
- Consider both explicit and implicit indicators
- Look for contextual clues about methodology
- Interpret typical field conventions
- Allow for variance in reporting styles
- Seek reasons to include when borderline

Confidence Scoring for Comprehensive Review:
- Score < 0.7: Mark as "uncertain" only when critical information is completely absent
- Score 0.7-0.85: Can infer required information from context
- Score > 0.85: Clear match with criteria, either explicit or strongly implied"""


task_prompt_text = """\
# Review protocol outputs

## Background to review:
{background}

## Research Question:
{research_question}

## Inclusion Criteria:
{inclusion_criteria}

## Exclusion Criteria:
{exclusion_criteria}

Please assess the following abstract for inclusion in the systematic review:

Your rationale must explicitly connect abstract content to specific criteria. \
For uncertain decisions, clearly state what additional information would be needed \
to make a confident decision.

# Study details:

## Title:
{title}

## Publishing Year:
{year}

## Journal:
{journal}

## Abstract:
{abstract}"""

conservative_reviewer_prompt = ChatPromptTemplate.from_messages(
    [("system", conservative_reviewer_prompt_text), ("human", task_prompt_text)]
)

comprehensive_reviewer_prompt = ChatPromptTemplate.from_messages(
    [("system", comprehensive_reviewer_prompt_text), ("human", task_prompt_text)]
)

# TODO: probably no need for two of these, can reuse one
llm1 = ChatOpenAI(model="gpt-4o", temperature=0)
llm2 = ChatOpenAI(model="gpt-4o", temperature=0)

# These return Pydantic models
llm1_with_structured_output = llm1.with_structured_output(schemas.ScreeningResponse)
llm2_with_structured_output = llm2.with_structured_output(schemas.ScreeningResponse)

screen_abstracts_chain = RunnableParallel(
    conservative=(
        conservative_reviewer_prompt | llm1_with_structured_output
    ).with_retry(
        stop_after_attempt=5,
        wait_exponential_jitter=True,
        retry_if_exception_type=(Exception,),
    ),
    comprehensive=(
        comprehensive_reviewer_prompt | llm2_with_structured_output
    ).with_retry(
        stop_after_attempt=5,
        wait_exponential_jitter=True,
        retry_if_exception_type=(Exception,),
    ),
).with_listeners(
    on_end=screen_abstracts_chain_on_end_cb,
    on_error=chain_on_error_listener_cb,
)  # .with_types(
#        output_type=ScreenAbstractsChainOutputDict # pyright: ignore [reportArgumentType]
# )

logger.info("chain: {!r}", screen_abstracts_chain)


class ScreenAbstractsBatchOutput(t.NamedTuple):
    """Output tuple of screen_abstracts_batch().

    Attributes:
        results (list[ScreenAbstractResultTuple]): list of screening results
        cb (OpenAICallbackHandler): callback handler

    Note:
        cb fields:
            cb.total_cost: float
            cb.total_tokens: int
            cb.prompt_tokens: int
            cb.prompt_tokens_cached: int
            cb.completion_tokens: int
            cb.reasoning_tokens: int
            cb.successful_requests: int

    """

    results: list[ScreenAbstractResultTuple]
    cb: OpenAICallbackHandler


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def screen_abstracts_batch(
    batch: list[models.SearchResult], batch_idx: int, review: models.SystematicReview
) -> ScreenAbstractsBatchOutput | None:
    """Invoke screen_abstracts_chain on a batch of PubMed results.

    Results and errors are also written to ``st.session_state.screen_abstracts_results``
    and ``st.session_state.screen_abstracts_errors``.

    Caller should set ``st.session_state.screen_abstracts_batch_idx`` though it will
    be set if missing.

    Args:
        batch (list[SearchResult]): list of search results to be screened
        batch_idx (int): index of the batch
        review (SystematicReview): systematic review associated with this batch

    Returns:
        tuple[list[ScreenAbstractResultTuple], OpenAICallbackHandler]:
            list of screened search results and OpenAI cb for the batch run.
            None if an uncaught exception occured (check the logs).
    """
    batch_logger = logger.bind(batch_idx=batch_idx)
    chain_inputs = make_screen_abstracts_chain_input(batch, review)
    chain_outputs: list[ScreenAbstractResultTuple] = []

    with get_openai_callback() as cb_openai:
        try:
            # logger.debug(f"Invoking chain with inputs: {chain_inputs!r}")
            # OLD: res = st.session_state.screen_abstracts_chain.batch(**chain_inputs)
            res = screen_abstracts_chain.batch(
                **chain_inputs  # type: ignore
            )  # Use globally defined chain

            # Ensure results are correctly formed into ScreenAbstractResultTuple
            # The on_end listener (screen_abstracts_chain_on_end_cb) is supposed to populate
            for i, parallel_invocation in enumerate(res):
                search_result = batch[i]
                conservative = parallel_invocation.get(
                    ScreeningStrategyType.CONSERVATIVE, parallel_invocation
                )
                comprehensive = parallel_invocation.get(
                    ScreeningStrategyType.COMPREHENSIVE, parallel_invocation
                )

                if not isinstance(conservative, ScreeningResult):
                    conservative = ScreeningError(
                        search_result=search_result,
                        error=conservative,
                    )
                    batch_logger.error(f"Conservative reviewer error: {conservative!r}")
                else:
                    search_result.conservative_result_id = conservative.id

                if not isinstance(comprehensive, ScreeningResult):
                    comprehensive = ScreeningError(
                        search_result=search_result,
                        error=comprehensive,
                    )
                    batch_logger.error(
                        f"Comprehensive reviewer error: {comprehensive!r}"
                    )
                else:
                    search_result.comprehensive_result_id = comprehensive.id

                res_tuple = ScreenAbstractResultTuple(
                    search_result=search_result,
                    conservative_result=conservative,
                    comprehensive_result=comprehensive,
                )
                chain_outputs.append(res_tuple)

            return ScreenAbstractsBatchOutput(
                results=chain_outputs,
                cb=deepcopy(cb_openai),
            )
        except Exception:
            batch_logger.exception("Exception occurred during screen_abstracts_batch")
            return None


# TODO: update
# @logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None) # pyright: ignore [reportArgumentType]
# async def ascreen_abstracts_batch(batch: list[models.PubMedResult], batch_idx: int) -> ScreenAbstractsBatchOutput | None:
#    """Invoke screen_abstracts_chain on a batch of PubMed results.
#
#    Results and errors are also written to ``st.session_state.screen_abstracts_results``
#    and ``st.session_state.screen_abstracts_errors``.
#
#    Caller should set ``st.session_state.screen_abstracts_batch_idx`` though it will
#    be set if missing.
#
#    Args:
#        batch (list[models.PubMedResult]): list of pubmed results to be screened
#        batch_idx (int): index of the batch
#
#    Returns:
#        tuple[list[ScreenAbstractResultTuple], OpenAICallbackHandler]:
#            list of screened pubmed results and OpenAI cb for the batch run.
#            None if an uncaught exception occured (check the logs).
#    """
#    ut.init_state_key("screen_abstracts_batch_idx", batch_idx)
#    ut.init_state_key("screen_abstracts_results", [])
#    ut.init_state_key("screen_abstracts_errors", [])
#    ut.init_state_key("screen_abstracts_chain", screen_abstracts_chain)
#    #ut.init_state_key("screen_abstracts_total_tokens", 0)
#    #ut.init_state_key("screen_abstracts_successful_requests", 0)
#    #ut.init_state_key("screen_abstracts_total_cost", 0.0)
#
#    batch_logger = logger.bind(batch_idx=batch_idx)
#    chain_inputs = make_screen_abstracts_chain_input(batch)
#    chain_outputs: list[ScreenAbstractResultTuple] = []
#
#    with get_openai_callback() as cb:
#        res =  await st.session_state.screen_abstracts_chain.abatch(**chain_inputs)
#    if len(res) != len(batch):
#        pad_size = len(batch) - len(res)
#        res += [None] * pad_size
#        batch_logger.warning(f"screen_abstracts_chain returned {len(res)} results, expected {len(batch)}")
#    for i, parallel_invocation in enumerate(res):
#        if isinstance(parallel_invocation, dict):
#            conservative = parallel_invocation.get(ScreeningStrategyType.CONSERVATIVE, parallel_invocation)
#            comprehensive = parallel_invocation.get(ScreeningStrategyType.COMPREHENSIVE, parallel_invocation)
#        else:
#            conservative = parallel_invocation
#            comprehensive = parallel_invocation
#        if not isinstance(conservative, ScreeningResult):
#            conservative = ScreeningError(
#                search_result=batch[i],
#                error=conservative,
#            )
#            if ut.is_streamlit():
#                st.session_state.screen_abstracts_errors.append(conservative)
#            batch_logger.warning(f"Conservative reviewer error: {conservative!r}")
#        else:
#            if ut.is_streamlit():
#                st.session_state.screen_abstracts_results.append(conservative)
#
#        if not isinstance(comprehensive, ScreeningResult):
#            comprehensive = ScreeningError(
#                search_result=batch[i],
#                error=comprehensive,
#            )
#            if ut.is_streamlit():
#                st.session_state.screen_abstracts_errors.append(comprehensive)
#            batch_logger.warning(f"Comprehensive reviewer error: {comprehensive!r}")
#        else:
#            if ut.is_streamlit():
#                st.session_state.screen_abstracts_results.append(comprehensive)
#
#        res_tuple = ScreenAbstractResultTuple(
#            search_result=batch[i],
#            conservative_result=conservative,
#            comprehensive_result=comprehensive,
#        )
#        chain_outputs.append(res_tuple)
#    return ScreenAbstractsBatchOutput(
#        results=chain_outputs,
#        cb=cb,
#    )


resolver_chain = resolver_prompt | resolver_model


@logger.catch(
    onerror=lambda exc: (
        st.error(f"Resolver Chain Error: {exc!r}") if ut.in_streamlit() else None,
        None,
    )[-1]
)  # pyright: ignore [reportArgumentType]
def invoke_resolver_chain(
    search_result: models.SearchResult,
    review: models.SystematicReview,
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
    run_config: RunnableConfig | None = None,
) -> ResolverOutputSchema | None:
    """Invokes the resolver_chain with the given inputs and configuration.

    Args:
        search_result: The search result to be resolved.
        review: The systematic review protocol.
        conservative_result: The screening result from the conservative reviewer.
        comprehensive_result: The screening result from the comprehensive reviewer.
        run_config: Optional LangChain RunnableConfig for the invocation.

    Returns:
        ResolverOutputSchema if successful, None otherwise (error will be logged by @logger.catch).
    """
    chain_input_model = prepare_resolver_inputs_for_prompt(
        search_result=search_result,
        review=review,
        conservative_result=conservative_result,
        comprehensive_result=comprehensive_result,
    )
    chain_input_dict = chain_input_model.model_dump()

    logger.info(f"Invoking resolver chain for SearchResult ID: {search_result.id!r}")
    try:
        response = resolver_chain.invoke(chain_input_dict, config=run_config)
        if isinstance(response, ResolverOutputSchema):
            logger.success(
                f"Resolver chain completed for SearchResult ID: {search_result.id!r}"
            )
            return response

        logger.error(
            f"Resolver chain for SearchResult ID: {search_result.id!r} returned unexpected type: {type(response)!r}. Response: {response!r}"
        )
        return None
    except Exception:  # Exception details will be logged by @logger.catch
        # The @logger.catch decorator handles logging the exception.
        # This block ensures 'None' is returned to the caller on error.
        logger.exception(
            f"Resolver chain invocation failed internally for SearchResult ID: {search_result.id!r}. See @logger.catch details."
        )
        return None


# class ResolverChainBatchInputDict(t.TypedDict):
#     \"\"\"ressolve_screening_result_batch.batch(**this_input_dict).
# No more 'existing code' after this point, ensure the file ends cleanly.
