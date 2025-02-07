"""Screening prompts and agents/chains. Also Streamlit chat UI agent(s).

Considering we're using OpenAI for the prototype, it doesn't make sense to use two
different models to emulate human reviewers. Instead, we'll use the same model but
with different prompts to better simulate human reviewers.

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
a threadpool executor or a queue with, e.g., Supabase edge function workers. For this
initial prototype, we'll stay within Streamlit and leverage LangChain runnables'
parallel and batching capabilities, which use a threadpool executor under the hood.
"""

from __future__ import annotations

import typing as t
from collections.abc import AsyncGenerator, Sequence, Callable
from datetime import datetime, timezone
import uuid

import streamlit as st
from langchain_core.messages import (
    AIMessageChunk,
    BaseMessage,
    BaseMessageChunk,
    MessageLikeRepresentation,
    SystemMessage,
)
from pydantic import BaseModel, ConfigDict, Field
from langchain_core.runnables.schema import StreamEvent
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnableConfig, RunnableParallel, Runnable
from langchain_core.tools import BaseTool, tool
from langchain_core.tracers.schemas import Run
from langchain_core.language_models import LanguageModelLike
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.postgres.aio import AsyncPostgresSaver
from langgraph.graph.graph import CompiledGraph
from langgraph.prebuilt.chat_agent_executor import create_react_agent, AgentState
from loguru import logger
from psycopg_pool import AsyncConnectionPool

from sr_assistant.core.schemas import ScreeningResponse, ScreeningResult
from sr_assistant.core.types import ScreeningStrategyType
from sr_assistant.core import models

# These are improved LangGraph types as originals cause type errors
# Ref: https://github.com/langchain-ai/langgraph/blob/9fe61a42d2bf63c4ba4a94de71e6e8016729a951/libs/langgraph/langgraph/prebuilt/chat_agent_executor.py#L46C1-L53C36

# This is the graph state, i.e. what tracks node return values
# See the above link for the required keys, additional ones can be added
StateSchema_co = t.TypeVar("StateSchema_co", bound=AgentState, covariant=True)
type StateSchemaType[StateSchema_co] = type[StateSchema_co]

# LangGraph agent prompts are flexible
type GraphAgentPrompt[StateSchema_co] = (
    SystemMessage
    | str
    | Callable[[StateSchema_co], Sequence[BaseMessage]]
    | Runnable[StateSchema_co, Sequence[BaseMessage]]
)
"""Can take a few different forms:

- str: This is converted to a SystemMessage and added to the beginning of the list of messages in state["messages"].
- SystemMessage: this is added to the beginning of the list of messages in state["messages"].
- Callable: This function should take in full graph state and the output is then passed to the language model.
- Runnable: This runnable should take in full graph state and the output is then passed to the language model.
    E.g., a ChatPromptTemplate.

Example:
    >>> from typing import TypedDict
    >>>
    >>> from langgraph.managed import IsLastStep
    >>> from langgraph.prebuilt.chat_agent_executor import AgentState
    >>> prompt = ChatPromptTemplate.from_messages(
    ...     [
    ...         ("system", "Today is {today}"),
    ...         ("placeholder", "{messages}"),
    ...     ]
    ... )
    >>>
    >>> class CustomState(AgentState):
    ...     today: str
    >>>
    >>> graph = create_react_agent(model, tools, state_schema=CustomState, prompt=prompt)
    >>> inputs = {"messages": [("user", "What's today's date? And what's the weather in SF?")], "today": "July 16, 2004"}
    >>> for s in graph.stream(inputs, stream_mode="values"):
    ...     message = s["messages"][-1]
    ...     if isinstance(message, tuple):
    ...         print(message)
    ...     else:
    ...         message.pretty_print()
"""

class AgentGraphInput(t.TypedDict):
    """Input to a LangGraph agent."""

    messages: list[MessageLikeRepresentation]

class ScreenAbstractsChainOutputType(t.TypedDict):
    conservative: ScreeningResult
    comprehensive: ScreeningResult

class ScreenAbstractsChainInput(BaseModel):
    """screen_abstracts_chain prompt input variables."""

    model_config = ConfigDict(
        from_attributes=True,
        validate_assignment=True,
        extra="ignore",
    )
    background: str | None = Field(default="<no background provided>")
    research_question: str
    inclusion_criteria: str
    exclusion_criteria: str
    title: str
    journal: str
    year: str
    abstract: str


class ScreenAbstractsChainInputDict(t.TypedDict):
    """Represents ScreenAbstractsChainInput.model_dump(). For type checking reasons."""
    background: str
    research_question: str
    inclusion_criteria: str
    exclusion_criteria: str
    title: str
    journal: str
    year: str
    abstract: str

class ScreenAbstractsChainBatchInputType(t.TypedDict):
    inputs: list[ScreenAbstractsChainInputDict]
    config: list[RunnableConfig]

def screening_chain_on_end_listener(run_obj: Run) -> None:
    """Listener for screening batch run on_end hook.

    See `.with_listeners <https://python.langchain.com/api_reference/_modules/langchain_core/runnables/base.html#RunnableBinding.with_listeners>`_

    Modifications to the `Run`|`RunTree` persist and no need to return anything. It's a
    Pydantic v1 model.

    Todo:
        - More metadata extraction
    """
    if not run_obj.outputs:
        logger.warning(f"No outputs for run: {run_obj!r}")
        return
    # Each child run is one screening chain invocation within RunnableParallel. This
    # parent RunTree is the RunnableParallel itself.
    for cr in run_obj.child_runs:
        if cr.outputs is None:
            logger.warning(f"No outputs for child run: {cr!r}")
            continue
        orig_resp_model = cr.outputs.get("output")
        if not isinstance(orig_resp_model, ScreeningResponse):
            logger.error(f"Unknown response type: {orig_resp_model!r}")
            continue
        review_id = cr.metadata.get("review_id")
        if not review_id:
            logger.error(f"Missing review_id: {cr.metadata!r}")
            continue
        review_id = uuid.UUID(review_id)
        pubmed_result_id = cr.metadata.get("pubmed_result_id")
        if not pubmed_result_id:
            logger.error(f"Missing pubmed_result_id: {cr.metadata!r}")
            continue
        pubmed_result_id = uuid.UUID(pubmed_result_id)
        result_id = cr.id
        trace_id = cr.trace_id # this is shared between invocations with same input
        start_time = cr.start_time
        end_time = cr.end_time
        inputs = cr.inputs
        if not end_time:
            logger.warning("No end time for run: {cr!r}, setting to now")
            end_time = datetime.now(tz=timezone.utc)
        resp_metadata = {}
        invocation_params = {}
        for ccr in cr.child_runs:
            if ccr.run_type == "chat_model":
                resp_metadata = ccr.extra.get("metadata", {})
                invocation_params = ccr.extra.get("invocation_params", {})
        model_name = resp_metadata.get("ls_model_name", invocation_params.get("model_name", "not_found"))
        output_key: ScreeningStrategyType
        if isinstance(cr.tags, list):
            maybe_output_key = next(
                (t.removeprefix("map:key:") for t in cr.tags if t.startswith("map:key:")),
                None,
            )
            if maybe_output_key:
                output_key = ScreeningStrategyType(maybe_output_key) # make type checker happy
            else:
                output_key = cr.metadata["screening_strategy"]
        else:
            output_key = cr.metadata["screening_strategy"]
        screening_result = ScreeningResult(
            id=result_id,
            trace_id=trace_id, # this is shared between invocations with same input
            review_id=review_id,
            pubmed_result_id=pubmed_result_id,
            start_time=start_time,
            end_time=end_time,
            model_name=model_name,
            screening_strategy=output_key,
            response_metadata=dict(inputs=inputs, run_name=cr.run_name, **resp_metadata, **invocation_params),
            **orig_resp_model.model_dump(),
        )
        run_obj.outputs[output_key] = screening_result
        # cr.outputs["output"].additional_kwargs.update({"sr_metadata": cr.metadata})
        # print(f"metadata: {cr.metadata!r}\noutputs: {cr.outputs!r}")

# LangGraph ReAct agent we use takes a system message
# directly as a prompt. It does support templates but
# they're more work for little gain in our use case.
#
# It's chat history and full record of every message
# is persisted in Postgres and can be used to construct
# audit trails or to better understand how AI is used
# in the project.
#
# The graph state handles chat history automatically.
# This message will simply be the first in state["messages"].
chat_assistant_system_prompt_1 = SystemMessage(
    content="""\
You're' a helpful systematic review expert assisting in developing \
a cutting-edge AI-driven systematic review platform and operating it.

## Responsibilities

You're required to review and grade every step of the process \
ensuring the highest standards. For example, you'll be asked to \
assist in developing the review protoco including developing \
PICO, etc. using FINER.

Likewise you'll assist with the search strategy helping to \
leverage query expansion and block building. The system \
leverages an ensemble LLM approach to screening, so you \
won't be directly involved in screening but may be asked to \
evaluate prompts and grade LLM responses.

You're equipped with a permanent memory that lasts for the \
duration of a given systematic review: your memory is \
keyed by the review ID. As such expect to be assigned \
a wide variety of tasks which may change over time and often."""
)

# create_react_agent requires tools but is the simplest way
# to get started with a good graph quickly, so we create a
# simple tool.

@tool
def get_current_time() -> str:
    """Returns current time and date."""
    return datetime.now(tz=timezone.utc).isoformat()

#class AgentGraph:
#    """LangGraph agent wrapper."""
#
#    def __init__(
#        self,
#        model: LanguageModelLike | None = None,
#        prompt: GraphAgentPrompt[StateSchema_co] = chat_assistant_system_prompt_1,
#        *,
#        config: RunnableConfig | None = None,
#        tools: list[BaseTool] | None = None,
#        db_uri: str = st.session_state.config.DATABASE_URL,
#        default_thread_id: str | uuid.UUID | None = None,
#    ) -> None:
#        self.prompt = prompt
#        self.tools = tools
#        self.db_uri = db_uri
#        self.model = model or ChatOpenAI(
#            model="gpt-4o", temperature=0, tags=["sr_assistant:prototype"]
#        )
#        self.tools = tools
#        self.db_uri = db_uri
#        self.default_thread_id = default_thread_id
#        self.checkpoint = None
#        self.checkpointer = None
#        self._checkpointer_conn_kwargs = {
#            "autocommit": True,
#            "prepare_threshold": 0,
#        }

# This needs refactoring and isn't used currently
async def ainvoke_chat_agent(
    messages: AgentGraphInput,
    prompt: GraphAgentPrompt[StateSchema_co] = chat_assistant_system_prompt_1,
    *,
    tools: list[BaseTool] | None = None,
    db_uri: str = st.session_state.config.DATABASE_URL,
    thread_id: str = str(st.session_state.rewiew_id),  # type: ignore
) -> AsyncGenerator[BaseMessageChunk]:
    """Systematic Review chat assistant in Streamlit UI.

    TODO:
        - Turn to a class, simplify uses, helpers for steamlit chat
        - invocation should require only passing one humanmessage
        - How to make this proper async considering Streamllit?
          - loop executor for graph invocation? checkpointer already
            uses executor. concur futures would be better but this
            is all async for many reasons.
        - Shouold be split into reusable base class/components
          since reading all docs and configuring everythihng takes
          time and we may want to use this for screening and
          past the prototype (curreently simple chains for that)
        - methods for working with checkpointer, doing restores, etc.
        - Tools: access Streamlit UIs so it can write for users,
          e.g., PICO from question, etc.
        - Tools: Entrez pubmed query gen/expansion
        - Tools: MeSH etc. autocomplete with Entrez and LLM
        - Tools to query PG and help with e.g. PRISMA
          and data extraction tasks
        - Tools to help analyze screeening results
        - Tool to screen for duplicates, e.g., embeddings
        - Tools to analyze and assist with protocol and search strategy
    """
    if not tools:
        tools = [get_current_time]
    model = ChatOpenAI(model="gpt-4o", temperature=0, tags=["sr_assistant:prototype"])
    connection_kwargs = {
        "autocommit": True,
        "prepare_threshold": 0,
    }
    async with AsyncConnectionPool(
        conninfo=db_uri,
        max_size=20,
        kwargs=connection_kwargs,
    ) as pool:
        checkpointer = AsyncPostgresSaver(pool)  # type: ignore
        # Uses AgentState
        graph: CompiledGraph = create_react_agent(
            prompt=prompt, model=model, tools=tools, checkpointer=checkpointer
        )
        config = RunnableConfig(
            configurable={"thread_id": thread_id},
            recursion_limit=1000,
            run_name="sr_assistant:react-agent",
            tags=["sr_assistant:prototype", "step1", "step2", "step3", "streamlit"],
            max_concurrency=50,
        )
        # StreamEvent:
        #'event': 'on_chat_model_stream',
        # 'name': 'ChatOpenAI',
        # 'run_id': '3fdbf494-acce-402e-9b50-4eab46403859',
        # 'tags': ['seq:step:1'],
        # 'metadata': {'langgraph_step': 1,
        #  'langgraph_node': 'call_model',
        #  'langgraph_triggers': ['start:call_model'],
        #  'langgraph_task_idx': 0,
        #  'checkpoint_id': '1ef657a0-0f9d-61b8-bffe-0c39e4f9ad6c',
        #  'checkpoint_ns': 'call_model',
        #  'ls_provider': 'openai',
        #  'ls_model_name': 'gpt-4o-mini',
        #  'ls_model_type': 'chat',
        #  'ls_temperature': 0.7},
        # 'data': {'chunk': AIMessageChunk(content='Hello', id='run-3fdbf494-acce-402e-9b50-4eab46403859')},
        # 'parent_ids': []}
        async for event in graph.astream_events(messages, config=config, version="v2"):
            kind = event["event"]
            logger.info(f"{kind!r}: {event['name']!r}")
            logger.debug(f"StreamEvent: {event!r}")
            if kind == "on_chat_model_stream":
                yield event["data"].get(
                    "chunk", AIMessageChunk(content="")
                )  #  langchain_core.runnables.schema.StreamEvent

        # TODO: do something with this, accept a session state key as param?
        await checkpointer.aget(config)


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


def make_screen_abstracts_chain_input(
    pubmed_results_batch: list[models.PubMedResult],
) -> ScreenAbstractsChainBatchInputType:
    """Create input for a batch to be passed to abstracts screening chain.

    chain.batch(input) input as returned from this function.

    TODO: We pass the same input and config to both reviewers, is it possible
          to bind configs to the sub chains such that metadata from the parallel
          chain doesn't overwrite anything?

          Currently reviewer/strategy must be looked up in the listner as it's
          populated by the chain in the RunTree.outputs keys and RunTree.tags.

    Args:
        pubmed_results_batch (list[models.PubMedResult]): list of pubmed results to be screened

    Returns:
        ScreenAbstractsChainBatchInputType: dict of list of prompt inputs and list of RunnableConfigs
            Can be passed directly to chain.batch()
    """
    #batch_id = uuid.uuid4()
    #if "screening_batches" in st.session_state:
    #    st.session_state.screening_batches = {}
    #st.session_state.screening_batch_id.update(
    #        {f"{batch_id}": {"pubmed_result_batch": pubmed_results_batch}, "created_at": datetime.now(tz=timezone.utc)}
    #)
    logger.info(f"Creating batch input for {len(pubmed_results_batch)} pubmed results")
    chain_inputs: list[ScreenAbstractsChainInputDict] = [] # input is same regardless of strategy, same prompt variables
    configs: list[RunnableConfig] = []
    for i, pubmed_result in enumerate(pubmed_results_batch):
        configs.append(
            RunnableConfig(
                run_name=f"screen_abstracts_chain",
                metadata={
                    "review_id": pubmed_result.review_id,
                    "pubmed_result_id": pubmed_result.id,
                },
                tags=[
                    "sra:prototype",
                    "sra:streamlit",
                    "sra:step3",
                    f"sra:review_id:{pubmed_result.review_id}",
                    f"sra:pubmed_result_id:{pubmed_result.id}",
                    f"sra:screen_abstracts_chain:i:{i}",
                ],  # srastep3:screning_strategy:comprehensive|conservative # added by listener
            )
        )
        chain_inputs.append(
            ScreenAbstractsChainInputDict(
                **ScreenAbstractsChainInput(
                    background=st.session_state.background,
                    research_question=st.session_state.research_question,
                    inclusion_criteria=st.session_state.inclusion_criteria,
                    exclusion_criteria=st.session_state.exclusion_criteria,
                    title=pubmed_result.title,
                    journal=pubmed_result.journal,
                    year=pubmed_result.year,
                    abstract=pubmed_result.abstract,
                ).model_dump()
            )
        )
        logger.info(f"Created chain input for pubmed result {i}: config: {configs[-1]!r}, input: {chain_inputs[-1]!r}")
    return ScreenAbstractsChainBatchInputType(inputs=chain_inputs, config=configs)

llm1 = ChatOpenAI(model="gpt-4o", temperature=0).with_retry(
    stop_after_attempt=5,
    wait_exponential_jitter=True,
    retry_if_exception_type=(Exception,),
)
llm2 = ChatOpenAI(model="gpt-4o", temperature=0).with_retry(
    stop_after_attempt=5,
    wait_exponential_jitter=True,
    retry_if_exception_type=(Exception,),
)

# These return Pydantic models
llm1_with_structured_output = llm1.with_structured_output(ScreeningResponse)
llm2_with_structured_output = llm2.with_structured_output(ScreeningResponse)

conservative_reviewer_prompt = ChatPromptTemplate.from_messages(
    [("system", conservative_reviewer_prompt_text), ("human", task_prompt_text)]
)

comprehensive_reviewer_prompt = ChatPromptTemplate.from_messages(
    [("system", comprehensive_reviewer_prompt_text), ("human", task_prompt_text)]
)

screen_abstracts_chain = RunnableParallel(
    conservative=(conservative_reviewer_prompt | llm1_with_structured_output),
    comprehensive=(
        comprehensive_reviewer_prompt | llm2_with_structured_output
    ),
).with_listeners(on_end=screening_chain_on_end_listener).with_retry(
    stop_after_attempt=5,
    wait_exponential_jitter=True,
    retry_if_exception_type=(Exception,),
).with_types(,output_type=ScreenAbstractsChainOutputType)
