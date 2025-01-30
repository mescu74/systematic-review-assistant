"""Screening prompts and agents/chains.

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

from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnableParallel
from langchain_openai import ChatOpenAI

from sr_assistant.core.schemas.screening import ScreeningResponse

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
Research Question: {research_question}

Inclusion Criteria:
{inclusion_criteria}

Exclusion Criteria:
{exclusion_criteria}

Please assess the following abstract for inclusion in the systematic review:

Your rationale must explicitly connect abstract content to specific criteria. \
For uncertain decisions, clearly state what additional information would be needed \
to make a confident decision.

Title: {title}
Year: {year}
Abstract: {abstract}"""


llm1 = ChatOpenAI(model="gpt-4o", temperature=0)
llm2 = ChatOpenAI(model="gpt-4o", temperature=0)

# These return Pydantic models
llm1_with_structured_output = llm1.with_structured_output(ScreeningResponse)
llm2_with_structured_output = llm2.with_structured_output(ScreeningResponse)

conservative_reviewer_prompt = ChatPromptTemplate.from_messages([
    ("system", conservative_reviewer_prompt_text),
    ("human", task_prompt_text)
])

comprehensive_reviewer_prompt = ChatPromptTemplate.from_messages([
    ("system", comprehensive_reviewer_prompt_text),
    ("human", task_prompt_text)
])

# In [2]: b = [{"research_question": "test question", "inclusion_criteria": "test criteria", "exclusion_criteria": "test criteria", "title": "test title", "year": "test year", "abstract": "test abstract"}]
# In [3]: res = abstract_screening_chain.batch(b)
# In [4]: res
# Out[4]:
# {'llm1_response': ScreeningResponse(decision='uncertain', confidence_score=0.5, rationale='The abstract provided does not contain any specific information that can be directly matched to the inclusion or exclusion criteria. Without explicit details on the population, intervention, comparison, outcomes, or study design, it is impossible to determine if the study meets the criteria. Additional information on these aspects would be necessary to make a confident decision.', exclusion_reason_categories=None),
# 'llm2_response': ScreeningResponse(decision='uncertain', confidence_score=0.5, rationale="The abstract provided is a placeholder ('test abstract') and does not contain any specific information about the study's population, intervention, comparison, outcomes, or study design. Without this information, it is impossible to determine whether the study meets the inclusion or exclusion criteria. Additional information about the study's methodology, results, and relevance to the research question is needed to make a confident decision.", exclusion_reason_categories=None)}]

abstract_screening_chain = RunnableParallel(
    llm1_response=(conservative_reviewer_prompt | llm1_with_structured_output),
    llm2_response=(comprehensive_reviewer_prompt | llm2_with_structured_output)
)
