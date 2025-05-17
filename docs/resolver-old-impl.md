
## Resolver model

The resolver model is `gemini-2.5-pro-preview-05-06` and implemented using `ChatGoogleGenerativeAI` LangChain class. 
- We'll configure it with a coding tool from Google GenAI in case it needs to create a model to evaluate the situation, or perform calculations.

### Reasoning mode

Note: "Thinking is on by default in both the API and AI Studio because the 2.5 series models have the ability to automatically decide when and how much to think based on the prompt."
- So the resolver prompt **must** encourage deep thinking and analytical mindset! And MUST reference docs/prompt-engineering-guide.md.
- Note that we had an old resolver implementation in-progress that was abruptly stopped. It should serve as a reference for implementing the new chain and prompt. See the code below.
- The `thinking_budget` parameter gives the model guidance on the number of thinking tokens it can use when generating a response. A greater number of tokens is typically associated with more detailed thinking, which is needed for solving more complex tasks. thinkingBudget must be an integer in the range 0 to 24576. Setting the thinking budget to 0 disables thinking.
- `24576` is the `thinking_budget` we use (max). See @@https://ai.google.dev/gemini-api/docs/thinking#set-budget

### Model definition

```python
from pydantic import SecretStr
from langchain_google_genai import ChatGoogleGenerativeAI

from st_assistant.app.config import get_config

resolver_model = (
    ChatGoogleGenerativeAI(
        model="gemini-2.5-pro-preview-05-06", # TODO: this should come from settings via env, but for now we hardcode it. Can be refactored later.
        temperature=0,
        max_tokens=None,
        thinking_budget=24576  # TODO: also should come from settings, ideally tunable in UI. But for this story 2.2 we hardcode it.
        timeout=None,
        max_retries=5,
        api_key=get_settings().GOOGLE_API_KEY,  # SecretStr
        convert_system_message_to_human=True, # Gemini doesn't support system prompts
)
```

Ideally we'd read from env:

```python
resolver_model_params: dict[str, str | int | bool | SecretStr] = {
    "model": settings.RESOLVER_MODEL_NAME,
    "temperature": settings.RESOLVER_TEMPERATURE,
    "max_retries": settings.RESOLVER_MAX_RETRIES,
    "api_key": get_settings().GOOGLE_API_KEY,
    "safety_settings": RESOLVER_SAFETY_SETTINGS,
    "convert_system_message_to_human": True,
}
```

### Old prompt

```python
resolver_prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            """\
**Role**: You are an expert systematic review screening assistant, built to resolve disputes between reviewers with confidence, and make the final decision when woth reviewers are `uncertain`. You're an analytical critical deep thinker able to analyze the situation from all perspectives, follow reviewer's reasoning and use the quotes they extracted and exclusion reason categories if applicable.
- 

Two reviewers either disagreeing on whether a given Database search result (article) should be included in a systematic review, or both have flagged the article as "uncertain".

One reviewer is using a "comprehensive" screening strategy, and the other is using a "conservative" screening strategy. These prioritize sensitivity and accuracy as different human reviewers might.

You will be given the search result, the systematic review protocol, and the two reviewers' screening results.

Task: You must deeply analyze and think about the situation. Let's do this step by step, analyzing every perspective. You should assume the reviewers persona to better understand their reasoning. You've access to each reviewers rationale, decision, confidence score, extracted quotes supporting their decision, and if applicable, exclusion reasons mapping to PRISMA.

When both are "uncertain": Follow their train of thought, fully understand their reasoning and perspective. You may assign "uncertain" if you as well are not sure. Carefully analyze both reviewers reasoning, supporting extracted quotes, etc. and fill the `resolver_j
""",
        ),
        (
            "human",
            """\
<search_result>
{search_result}
</search_result>

<systematic_review>
{systematic_review}
</systematic_review>

<comprehensive_reviewer_result>
{comprehensive_reviewer_result}
</comprehensive_reviewer_result>

<conservative_reviewer_result>
{conservative_reviewer_result}
</conservative_reviewer_result>

<resolver_output>
{resolver_output}
</resolver_output>""",
        ),
    ]
)

# Create resolver chain with structured output
resolver_model_struct_output = RESOLVER_MODEL.with_structured_output(
    ResolverOutputSchema
)
resolver_chain = resolver_prompt | resolver_model_struct_output


def make_resolver_chain_input(
    search_result: models.SearchResult,
    review: models.SystematicReview,
    conservative_result: ScreeningResult,
    comprehensive_result: ScreeningResult,
) -> dict[str, str]:
    """Create input dict for screening conflict resolver chain.

    Args:
        search_result: The Database search result with conflicting screening decisions
        review: The systematic review
        conservative_result: The screening result from the conservative reviewer
        comprehensive_result: The screening result from the comprehensive reviewer

    Returns:
        dict: Input variables for the resolver prompt
    """
    search_result_formatted = f"""\
Title: {search_result.title}
Source ID: {search_result.source_id}
Journal: {search_result.journal or "<No Journal>"}
Year: {search_result.year or "<No Year>"}
Abstract: {search_result.abstract or "<No Abstract>"}
Keyworks: {search_result.keywords}
MeSH terms: {search_result.mesh_terms}
"""
    systematic_review_formatted = f"""\
Background: {review.background or "<No Background>"}
Research Question: {review.research_question}
Inclusion Criteria: {review.inclusion_criteria or "<Not Specified>"}
Exclusion Criteria: {review.exclusion_criteria or "<Not Specified>"}
"""
    conservative_reviewer_result_formatted = f"""\
Decision: {conservative_result.decision.value}
Confidence Score: {conservative_result.confidence_score}
Rationale: {conservative_result.rationale}
Extracted quotes: {conservative_result.extracted_quotes}
"""
    if conservative_result.extracted_quotes:
        conservative_reviewer_result_formatted += (
            f"Extractea Quotes: {'; '.join(conservative_result.extracted_quotes)}\n"
        )

    # Format comprehensive result for the prompt
    comprehensive_reviewer_result_formatted = f"""\
Decision: {comprehensive_result.decision.value}
Confidence Score: {comprehensive_result.confidence_score}
Rationale: {comprehensive_result.rationale}
"""
    if comprehensive_result.extracted_quotes:
        comprehensive_reviewer_result_formatted += (
            f"Extracted Quotes: {'; '.join(comprehensive_result.extracted_quotes)}\n"
        )

    return {
        "search_result": search_result_formatted,
        "systematic_review": systematic_review_formatted,
        "conservative_reviewer_result": conservative_reviewer_result_formatted,
        "comprehensive_reviewer_result": comprehensive_reviewer_result_formatted,
    }


def resolve_screening_conflict(
    search_result: models.SearchResult,
    review: models.SystematicReview,
    conservative_result: ScreeningResult,
    comprehensive_result: ScreeningResult,
) -> models.ScreeningResolution:
    """Resolve screening conflict or ambiguity.

    Args:
        search_result: Database search result with conflicting screening decisions
        review: Systematic review
        conservative_result: Screening result from the conservative reviewer
        comprehensive_result: Screening result from the comprehensive reviewer

    Returns:
        models.ScreeningResolution: The resolution decision

    Raises:
        Exception: If the resolver chain fails
    """
    logger.info(
        f"Resolving screening conflict for Source ID: {search_result.source_id}, conservative: {conservative_result.decision}, comprehensive: {comprehensive_result.decision}"
    )

    chain_input = make_resolver_chain_input(
        search_result=search_result,
        review=review,
        conservative_result=conservative_result,
        comprehensive_result=comprehensive_result,
    )
    try:
        # Assuming resolver_chain is defined and returns ResolverOutputSchema
        response: ResolverOutputSchema = resolver_chain.invoke(
            chain_input,
            config=RunnableConfig(
                run_name="resolve_screening_conflict",
                tags=[
                    "sra:resolver",
                    f"sra:review_id:{review.id}",
                    f"sra:search_result_id:{search_result.id}",
                ],
                metadata={
                    "review_id": str(review.id),
                    "search_result_id": str(search_result.id),
                    "conservative_result_id": str(conservative_result.id),
                    "comprehensive_result_id": str(comprehensive_result.id),
                },
            ),
        )
        logger.info(f"Resolver response: {response!r}")
        resolution = models.ScreeningResolution(
            # id=uuid.uuid4(), # Let DB generate default ID
            review_id=review.id,
            search_result_id=search_result.id,
            conservative_result_id=conservative_result.id,
            comprehensive_result_id=comprehensive_result.id,
            resolver_decision=response.resolver_decision,
            resolver_reasoning=response.resolver_reasoning,
            resolver_confidence_score=response.resolver_confidence_score,
            response_metadata={"raw_output": response.model_dump()},
        )
        return resolution
    except Exception as e:
        logger.exception(
            f"Error resolving screening conflict for Source ID: {search_result.source_id}",
            exc_info=e,
        )
        return models.ScreeningResolution(
            # id=uuid.uuid4(), # Let DB generate default ID
            review_id=review.id,
            search_result_id=search_result.id,
            conservative_result_id=conservative_result.id,
            comprehensive_result_id=comprehensive_result.id,
            resolver_decision=ScreeningDecisionType.UNCERTAIN,
            resolver_reasoning=f"Error during resolution: {e}",
            resolver_confidence_score=0.0,
            response_metadata={"error": str(e)},
        )