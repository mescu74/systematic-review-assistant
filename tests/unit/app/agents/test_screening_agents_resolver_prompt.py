# pyright: reportPrivateUsage=false
#
import typing as t
import uuid
from datetime import datetime, timezone

import pytest
from pydantic import ValidationError

from sr_assistant.app.agents.screening_agents import (
    RESOLVER_SYSTEM_PROMPT,
    ResolverOutputSchema,
    ResolverPromptInput,
    _format_exclusion_reasons,
    _format_extracted_quotes,
    _format_raw_data_list,
    prepare_resolver_inputs_for_prompt,
    resolver_prompt,
)
from sr_assistant.core import models, schemas
from sr_assistant.core.types import (
    CriteriaFramework,
    ScreeningDecisionType,
    ScreeningStrategyType,
    SearchDatabaseSource,
)

# --- Tests for _format_extracted_quotes ---


def test_format_extracted_quotes_none():
    assert _format_extracted_quotes(None) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_extracted_quotes_empty_list():
    assert _format_extracted_quotes([]) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_extracted_quotes_with_strings():
    quotes = ["Quote 1.", "Another quote."]
    expected = "- Quote 1.\n- Another quote."
    assert _format_extracted_quotes(quotes) == expected  # pyright: ignore[reportPrivateUsage]


# --- Tests for _format_raw_data_list ---


def test_format_raw_data_list_none():
    assert _format_raw_data_list(None) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_empty_list():
    assert _format_raw_data_list([]) == ""  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_list_of_strings():
    data = ["term1", "term2"]
    assert _format_raw_data_list(data) == "term1, term2"  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_list_of_mixed():
    data = ["term1", 123, "term3"]
    assert _format_raw_data_list(data) == "term1, 123, term3"  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_single_string():
    assert _format_raw_data_list("single_term") == "single_term"  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_integer():
    assert _format_raw_data_list(123) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_raw_data_list_dict():
    assert _format_raw_data_list({"key": "value"}) == "N/A"  # pyright: ignore[reportPrivateUsage]


# --- Tests for _format_exclusion_reasons ---


def test_format_exclusion_reasons_none():
    assert _format_exclusion_reasons(None) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_exclusion_reasons_all_empty_or_none():
    reasons = schemas.ExclusionReasons(
        population_exclusion_reasons=None, intervention_exclusion_reasons=[]
    )
    assert _format_exclusion_reasons(reasons) == "N/A"  # pyright: ignore[reportPrivateUsage]


def test_format_exclusion_reasons_populated():
    reasons = schemas.ExclusionReasons(
        population_exclusion_reasons=[
            "Wrong age range for inclusion criteria",
            "Non-human subjects",
        ],
        intervention_exclusion_reasons=["Dosage outside specified range"],
        comparison_exclusion_reasons=None,
        outcome_exclusion_reasons=[],
        reporting_exclusion_reasons=None,
        study_design_exclusion_reasons=None,
    )
    expected = (
        "Population: Wrong age range for inclusion criteria, Non-human subjects\n"
        "Intervention: Dosage outside specified range"
    )
    assert _format_exclusion_reasons(reasons) == expected  # pyright: ignore[reportPrivateUsage]


def test_format_exclusion_reasons_only_one_category():
    reasons = schemas.ExclusionReasons(
        intervention_exclusion_reasons=["Dosage outside specified range"]
    )
    expected = "Intervention: Dosage outside specified range"
    assert _format_exclusion_reasons(reasons) == expected  # pyright: ignore[reportPrivateUsage]


# --- Tests for prepare_resolver_inputs_for_prompt ---


def test_prepare_resolver_inputs_basic():
    search_result_id = uuid.uuid4()
    review_id = uuid.uuid4()

    search_result = models.SearchResult(
        id=search_result_id,
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="pmid123",
        title="Test Title",
        abstract="Test Abstract.",
        journal="Test Journal",
        year="2023",
        keywords=["test", "mock"],
        raw_data={"MeshHeadings": ["mesh1", "mesh2"]},
    )
    review = models.SystematicReview(
        id=review_id,
        research_question="RQ?",
        inclusion_criteria="Incl Crit.",
        exclusion_criteria="Excl Crit.",
        background="Review Background.",
        criteria_framework=CriteriaFramework.PICO,
    )
    conservative_res = schemas.ScreeningResult(
        id=uuid.uuid4(),
        review_id=review_id,
        search_result_id=search_result_id,
        trace_id=uuid.uuid4(),
        model_name="test_model_conservative",
        screening_strategy=ScreeningStrategyType.CONSERVATIVE,
        start_time=datetime.now(timezone.utc),
        end_time=datetime.now(timezone.utc),
        decision=ScreeningDecisionType.INCLUDE,
        confidence_score=0.9,
        rationale="Con Rationale",
        extracted_quotes=["c_quote1"],
        exclusion_reason_categories=None,
    )
    comprehensive_res = schemas.ScreeningResult(
        id=uuid.uuid4(),
        review_id=review_id,
        search_result_id=search_result_id,
        trace_id=uuid.uuid4(),
        model_name="test_model_comprehensive",
        screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
        start_time=datetime.now(timezone.utc),
        end_time=datetime.now(timezone.utc),
        decision=ScreeningDecisionType.EXCLUDE,
        confidence_score=0.8,
        rationale="Comp Rationale",
        extracted_quotes=["comp_quote1"],
        exclusion_reason_categories=schemas.ExclusionReasons(
            population_exclusion_reasons=["Wrong age range for inclusion criteria"]
        ),
    )

    expected_input = ResolverPromptInput(
        article_title="Test Title",
        article_abstract="Test Abstract.",
        article_journal="Test Journal",
        article_year="2023",
        article_source_id="pmid123",
        article_keywords="test, mock",
        article_mesh_terms="mesh1, mesh2",
        review_background="Review Background.",
        review_research_question="RQ?",
        review_criteria_framework="PICO",
        review_inclusion_criteria="Incl Crit.",
        review_exclusion_criteria="Excl Crit.",
        conservative_reviewer_decision="include",
        conservative_reviewer_confidence=0.9,
        conservative_reviewer_rationale="Con Rationale",
        conservative_reviewer_quotes="- c_quote1",
        conservative_reviewer_exclusion_reasons="N/A",
        comprehensive_reviewer_decision="exclude",
        comprehensive_reviewer_confidence=0.8,
        comprehensive_reviewer_rationale="Comp Rationale",
        comprehensive_reviewer_quotes="- comp_quote1",
        comprehensive_reviewer_exclusion_reasons="Population: Wrong age range for inclusion criteria",
    )

    actual_input = prepare_resolver_inputs_for_prompt(
        search_result, review, conservative_res, comprehensive_res
    )
    assert actual_input == expected_input


def test_prepare_resolver_inputs_minimal_data():
    search_result_id = uuid.uuid4()
    review_id = uuid.uuid4()

    search_result = models.SearchResult(
        id=search_result_id,
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="pmid123_min",
        title="",
        abstract=None,
        journal=None,
        year=None,
        keywords=None,
        raw_data={},
    )
    review = models.SystematicReview(
        id=review_id,
        research_question="RQ? Min",
        exclusion_criteria="Excl Crit. Min",
        background="",
        inclusion_criteria=None,
        criteria_framework=None,
    )
    conservative_res = schemas.ScreeningResult(
        id=uuid.uuid4(),
        review_id=review_id,
        search_result_id=search_result_id,
        trace_id=uuid.uuid4(),
        model_name="test_model_conservative_min",
        screening_strategy=ScreeningStrategyType.CONSERVATIVE,
        start_time=datetime.now(timezone.utc),
        end_time=datetime.now(timezone.utc),
        decision=ScreeningDecisionType.UNCERTAIN,
        confidence_score=0.5,
        rationale="Rationale Min Con.",
        extracted_quotes=None,
        exclusion_reason_categories=None,
    )
    comprehensive_res = schemas.ScreeningResult(
        id=uuid.uuid4(),
        review_id=review_id,
        search_result_id=search_result_id,
        trace_id=uuid.uuid4(),
        model_name="test_model_comprehensive_min",
        screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
        start_time=datetime.now(timezone.utc),
        end_time=datetime.now(timezone.utc),
        decision=ScreeningDecisionType.UNCERTAIN,
        confidence_score=0.5,
        rationale="Rationale Min Comp.",
        extracted_quotes=None,
        exclusion_reason_categories=None,
    )

    expected_input = ResolverPromptInput(
        article_title="N/A",
        article_abstract="N/A",
        article_journal="N/A",
        article_year="N/A",
        article_source_id="pmid123_min",
        article_keywords="N/A",
        article_mesh_terms="N/A",
        review_background="N/A",
        review_research_question="RQ? Min",
        review_criteria_framework="N/A",
        review_inclusion_criteria="N/A",
        review_exclusion_criteria="Excl Crit. Min",
        conservative_reviewer_decision="uncertain",
        conservative_reviewer_confidence=0.5,
        conservative_reviewer_rationale="Rationale Min Con.",
        conservative_reviewer_quotes="N/A",
        conservative_reviewer_exclusion_reasons="N/A",
        comprehensive_reviewer_decision="uncertain",
        comprehensive_reviewer_confidence=0.5,
        comprehensive_reviewer_rationale="Rationale Min Comp.",
        comprehensive_reviewer_quotes="N/A",
        comprehensive_reviewer_exclusion_reasons="N/A",
    )
    actual_input = prepare_resolver_inputs_for_prompt(
        search_result, review, conservative_res, comprehensive_res
    )
    assert actual_input == expected_input


# --- Tests for resolver_prompt formatting ---


def test_resolver_prompt_formatting():
    input_data = ResolverPromptInput(
        article_title="Test Title",
        article_abstract="Test Abstract Content.",
        article_journal="Journal of Tests",
        article_year="2024",
        article_source_id="test_id_123",
        article_keywords="unit, test, keywords",
        article_mesh_terms="mesh_A, mesh_B",
        review_background="Background of the systematic review being tested.",
        review_research_question="What is the effect of testing on software quality?",
        review_criteria_framework="PICO",
        review_inclusion_criteria="Studies involving software testing.",
        review_exclusion_criteria="Studies not in English.",
        conservative_reviewer_decision="include",
        conservative_reviewer_confidence=0.95,
        conservative_reviewer_rationale="Conservative rationale: Looks good.",
        conservative_reviewer_quotes="- Relevant quote 1 from abstract by conservative.",
        conservative_reviewer_exclusion_reasons="N/A",
        comprehensive_reviewer_decision="include",
        comprehensive_reviewer_confidence=0.88,
        comprehensive_reviewer_rationale="Comprehensive rationale: Seems okay.",
        comprehensive_reviewer_quotes="- Relevant quote A from abstract by comprehensive.",
        comprehensive_reviewer_exclusion_reasons="N/A",
    )

    prompt_value = resolver_prompt.format_prompt(**input_data.model_dump())
    messages = prompt_value.to_messages()

    assert len(messages) == 2
    system_message_content = messages[0].content
    human_message_content = messages[1].content

    assert system_message_content == RESOLVER_SYSTEM_PROMPT

    assert "<title>Test Title</title>" in human_message_content
    assert "<abstract>\nTest Abstract Content.\n  </abstract>" in human_message_content
    assert "<criteria_framework>PICO</criteria_framework>" in human_message_content
    assert "<decision>include</decision>" in human_message_content
    assert (
        "<rationale>\nConservative rationale: Looks good.\n  </rationale>"
        in human_message_content
    )
    assert "- Relevant quote A from abstract by comprehensive." in human_message_content
    assert (
        "`resolver_decision`: Your final decision (e.g., 'include', 'exclude', 'uncertain')."
        in human_message_content
    )


# --- Tests for ResolverOutputSchema parsing ---


def test_resolver_output_schema_valid():
    valid_data = {
        "resolver_decision": ScreeningDecisionType.INCLUDE,
        "resolver_reasoning": "The study meets all criteria.",
        "resolver_confidence_score": 0.9,
        "contributing_strategies": [
            ScreeningStrategyType.CONSERVATIVE,
            ScreeningStrategyType.COMPREHENSIVE,
        ],
    }
    parsed = ResolverOutputSchema.model_validate(valid_data)
    assert parsed.resolver_decision == ScreeningDecisionType.INCLUDE
    assert parsed.resolver_reasoning == "The study meets all criteria."
    assert parsed.resolver_confidence_score == 0.9
    assert parsed.contributing_strategies == [
        ScreeningStrategyType.CONSERVATIVE,
        ScreeningStrategyType.COMPREHENSIVE,
    ]


def test_resolver_output_schema_invalid_confidence():
    invalid_data = {
        "resolver_decision": ScreeningDecisionType.EXCLUDE,
        "resolver_reasoning": "Reasoning here.",
        "resolver_confidence_score": 1.5,
        "contributing_strategies": [],
    }
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(invalid_data)


def test_resolver_output_schema_missing_field():
    # resolver_reasoning is missing
    invalid_data = {
        "resolver_decision": ScreeningDecisionType.UNCERTAIN,
        "resolver_confidence_score": 0.5,
        "contributing_strategies": [],
    }
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(invalid_data)


def test_resolver_output_schema_wrong_type():
    # Test invalid string for ScreeningDecisionType
    invalid_decision_data = {
        "resolver_decision": t.cast(
            "t.Any", "maybe"
        ),  # Cast to Any for Pyright with invalid enum string
        "resolver_reasoning": "A reason.",
        "resolver_confidence_score": 0.7,
        "contributing_strategies": [],
    }
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(invalid_decision_data)

    # Test invalid string in contributing_strategies list
    invalid_strategies_data = {
        "resolver_decision": ScreeningDecisionType.INCLUDE,
        "resolver_reasoning": "A reason.",
        "resolver_confidence_score": 0.7,
        "contributing_strategies": [
            ScreeningStrategyType.CONSERVATIVE,
            t.cast("t.Any", "invalid_strategy"),
        ],
    }
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(invalid_strategies_data)
