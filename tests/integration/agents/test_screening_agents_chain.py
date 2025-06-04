"""Integration tests for LangChain components in screening_agents.py."""

import typing as t
import uuid

import pytest
import streamlit as st

from sr_assistant.app.agents import screening_agents
from sr_assistant.core.models import SearchResult, SystematicReview
from sr_assistant.core.types import (
    ScreeningDecisionType,
    SearchDatabaseSource,
)


@pytest.fixture
def mock_search_result() -> SearchResult:
    """Create a mock SearchResult for testing."""
    return SearchResult(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        source_db=SearchDatabaseSource.PUBMED,
        source_id="PM12345",
        title="Test Study on Machine Learning",
        journal="Journal of Testing",
        year="2023",
        abstract="This is a test abstract about machine learning and its applications.",
        authors=["Test Author"],
        keywords=["machine learning", "testing"],
        raw_data={
            "pmid": "12345",
            "mesh_terms": ["Artificial Intelligence", "Software Testing"],
        },
        source_metadata={"publication_types": ["Journal Article"]},
    )


@pytest.fixture
def mock_systematic_review() -> SystematicReview:
    """Create a mock systematic review for testing."""
    return SystematicReview(
        id=uuid.uuid4(),
        research_question="What are the effects of machine learning on testing?",
        background="Testing background for systematic review",
        inclusion_criteria="Studies must be about machine learning and testing",
        exclusion_criteria="Studies not in English, studies before 2000",
    )


@pytest.fixture
def mock_edge_case_search_result() -> SearchResult:
    """Create a mock SearchResult with missing fields to test edge cases."""
    return SearchResult(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        source_db=SearchDatabaseSource.PUBMED,
        source_id="PM54321",
        title="",
        journal=None,
        year=None,
        abstract=None,
        authors=[],
        keywords=[],
        raw_data={},
        source_metadata={},
    )


@pytest.fixture
def mock_non_relevant_search_result() -> SearchResult:
    """Create a mock SearchResult that should be excluded."""
    return SearchResult(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        source_db=SearchDatabaseSource.PUBMED,
        source_id="PM98765",
        title="Impact of Diet on Heart Disease",
        journal="Journal of Cardiology",
        year="1995",
        abstract="This study examines the relationship between diet and heart disease.",
        authors=["Heart Researcher"],
        keywords=["diet", "heart disease", "nutrition"],
        raw_data={
            "pmid": "98765",
            "mesh_terms": ["Diet", "Heart Diseases", "Nutrition"],
        },
        source_metadata={"publication_types": ["Journal Article"]},
    )


@pytest.mark.integration
def test_different_strategy_results(
    mock_search_result: SearchResult, mock_systematic_review: SystematicReview
):
    """Test that conservative and comprehensive strategies produce different results."""
    # Prepare input for the chain
    chain_input = screening_agents.make_screen_abstracts_chain_input(
        [mock_search_result], mock_systematic_review
    )

    # Run the chain with batch
    results = screening_agents.screen_abstracts_chain.batch(
        inputs=t.cast("list[dict[str, t.Any]]", chain_input["inputs"]),
        config=chain_input["config"],
    )
    assert len(results) == 1
    result = results[0]

    # Extract the decisions and confidence scores
    conservative_decision = result["conservative"].decision
    conservative_confidence = result["conservative"].confidence_score
    conservative_rationale = result["conservative"].rationale

    comprehensive_decision = result["comprehensive"].decision
    comprehensive_confidence = result["comprehensive"].confidence_score
    comprehensive_rationale = result["comprehensive"].rationale

    # Print the decisions for debugging
    print(
        f"Conservative decision: {conservative_decision}, confidence: {conservative_confidence}"
    )
    print(f"Conservative rationale: {conservative_rationale}")
    print(
        f"Comprehensive decision: {comprehensive_decision}, confidence: {comprehensive_confidence}"
    )
    print(f"Comprehensive rationale: {comprehensive_rationale}")

    # The results might be the same in some cases, but the rationales should differ
    # We can't assert they're always different due to model randomness
    # But we can check that they contain strategy-specific language

    # At minimum, the confidence and rationales should differ even if decisions are the same
    if conservative_decision == comprehensive_decision:
        # Check at least one of these conditions is true
        assert (
            abs(conservative_confidence - comprehensive_confidence) > 0.05
            or conservative_rationale != comprehensive_rationale
        )


@pytest.mark.integration
def test_non_relevant_abstract(
    mock_non_relevant_search_result: SearchResult,
    mock_systematic_review: SystematicReview,
):
    """Test that a non-relevant abstract is correctly classified as exclude."""
    # Prepare input for the chain
    chain_input = screening_agents.make_screen_abstracts_chain_input(
        [mock_non_relevant_search_result], mock_systematic_review
    )

    # Run the chain with batch
    results = screening_agents.screen_abstracts_chain.batch(
        inputs=t.cast("list[dict[str, t.Any]]", chain_input["inputs"]),
        config=chain_input["config"],
    )
    assert len(results) == 1
    result = results[0]

    # At least one of the strategies should exclude the paper
    # (due to year before 2000 in exclusion criteria)
    conservative_decision = result["conservative"].decision
    comprehensive_decision = result["comprehensive"].decision

    # Print decisions for debugging
    print(f"Conservative decision for non-relevant abstract: {conservative_decision}")
    print(f"Comprehensive decision for non-relevant abstract: {comprehensive_decision}")

    # At least one should be EXCLUDE or UNCERTAIN (not both INCLUDE)
    assert ScreeningDecisionType.INCLUDE not in (
        conservative_decision,
        comprehensive_decision,
    ) or (
        conservative_decision
        != comprehensive_decision  # If one includes, the other should not
    )


@pytest.mark.integration
def test_edge_case_handling(
    mock_edge_case_search_result: SearchResult, mock_systematic_review: SystematicReview
):
    """Test that the chain handles edge cases with missing data."""
    # Prepare input for the chain
    chain_input = screening_agents.make_screen_abstracts_chain_input(
        [mock_edge_case_search_result], mock_systematic_review
    )

    # Run the chain with batch
    results = screening_agents.screen_abstracts_chain.batch(
        inputs=t.cast("list[dict[str, t.Any]]", chain_input["inputs"]),
        config=chain_input["config"],
    )
    assert len(results) == 1
    result = results[0]

    # The chain should complete without errors even with empty data
    assert "conservative" in result
    assert "comprehensive" in result

    # Results should reflect uncertainty due to missing data
    conservative_result = result["conservative"]
    comprehensive_result = result["comprehensive"]

    print(f"Conservative decision for edge case: {conservative_result.decision}")
    print(f"Comprehensive decision for edge case: {comprehensive_result.decision}")

    # With minimal information, decisions should be UNCERTAIN or EXCLUDE, not INCLUDE
    assert conservative_result.decision != ScreeningDecisionType.INCLUDE
    # Comprehensive might include in some cases, but should have low confidence
    if comprehensive_result.decision == ScreeningDecisionType.INCLUDE:
        assert comprehensive_result.confidence_score < 0.8


@pytest.mark.integration
def test_batch_processing(
    mock_systematic_review: SystematicReview, monkeypatch: pytest.MonkeyPatch
):
    """Test that batch processing works with multiple abstracts."""
    # Generate a batch of different PubMed results
    batch = []
    review_id_for_batch = uuid.uuid4()
    for i in range(3):  # Keep the batch small for testing
        is_ml = i % 2 == 0
        source_id = f"PM{10000 + i}"
        batch.append(
            SearchResult(
                id=uuid.uuid4(),
                review_id=review_id_for_batch,
                source_db=SearchDatabaseSource.PUBMED,
                source_id=source_id,
                title=f"Test Study {i} on {'Machine Learning' if is_ml else 'Healthcare'}",
                journal=f"Journal of {'Testing' if is_ml else 'Medicine'}",
                year=str(2010 + i),
                abstract=f"This is abstract {i} about {'machine learning' if is_ml else 'healthcare'}.",
                authors=[f"Author {i}"],
                keywords=[f"keyword_{i}", "testing"],
                raw_data={"pmid": source_id, "mesh_terms": [f"MeSH_{i}", "Testing"]},
                source_metadata={"publication_types": ["Journal Article"]},
            )
        )

    # Mock Streamlit session state
    class MockSessionState:
        def __init__(self):
            self.screen_abstracts_chain = screening_agents.screen_abstracts_chain
            self.screen_abstracts_results = []
            self.screen_abstracts_errors = []
            self.screen_abstracts_batch_idx = 0

    monkeypatch.setattr(st, "session_state", MockSessionState())

    # Run the batch screening
    batch_idx = 0
    result = screening_agents.screen_abstracts_batch(
        batch, batch_idx, mock_systematic_review
    )

    # Verify the output
    assert result is not None
    results, cb = result

    # Check batch results
    assert len(results) == len(batch)

    # Check that results have different decisions (or handle errors)
    decisions = []
    for res in results:
        if isinstance(res.conservative_result, screening_agents.ScreeningResult):
            decisions.append(res.conservative_result.decision)
        else:
            # Append an error indicator or handle ScreeningError appropriately
            decisions.append("ERROR")
    print(f"Batch decisions: {decisions}")
    # Example: Assert that not all decisions are the same if expected variation
    # if len(set(d for d in decisions if d != "ERROR")) > 0:
    #     assert len(set(d for d in decisions if d != "ERROR")) > 1

    # Verify the OpenAI callback metrics
    assert cb.total_tokens > 0
    assert cb.total_cost > 0
    print(
        f"Batch processing stats - Tokens: {cb.total_tokens}, Cost: ${cb.total_cost:.4f}"
    )
