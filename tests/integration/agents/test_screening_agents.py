"""Integration tests for screening_agents.py."""

import pytest
from langchain.schema.runnable import RunnableConfig
from langchain_core.outputs import RunInfo
import streamlit as st

from sr_assistant.app.agents import screening_agents
from sr_assistant.core.models import PubMedResult, SystematicReview
from sr_assistant.core.schemas import ScreeningResponse, ScreeningResult


@pytest.fixture
def mock_pubmed_result() -> PubMedResult:
    """Create a mock PubMed result for testing."""
    return PubMedResult(
        pmid="12345",
        title="Test Study on Machine Learning",
        journal="Journal of Testing",
        year="2023",
        abstract="This is a test abstract about machine learning and its applications.",
        authors="Test Author",
        keywords=["machine learning", "testing"],
        mesh_terms=["Artificial Intelligence", "Software Testing"],
        publication_types=["Journal Article"],
    )


@pytest.fixture
def mock_systematic_review() -> SystematicReview:
    """Create a mock systematic review for testing."""
    return SystematicReview(
        title="Test Review",
        background="Testing background for systematic review",
        research_question="What are the effects of machine learning on testing?",
        inclusion_criteria="Studies must be about machine learning and testing",
        exclusion_criteria="Studies not in English, studies before 2000",
    )


def test_screen_abstracts_chain(mock_pubmed_result, mock_systematic_review):
    """Test that screen_abstracts_chain correctly processes a single abstract."""
    # Prepare input for the chain
    chain_input = screening_agents.make_screen_abstracts_chain_input(
        [mock_pubmed_result], mock_systematic_review
    )

    # Run the chain with batch
    results = screening_agents.screen_abstracts_chain.batch(
        inputs=chain_input["inputs"],
        config=chain_input["config"],
    )
    assert len(results) == 1
    result = results[0]

    # Verify the output structure
    assert isinstance(result, dict)
    assert "conservative" in result
    assert "comprehensive" in result

    # Check conservative reviewer result
    assert isinstance(result["conservative"], ScreeningResult)
    assert hasattr(result["conservative"], "decision")
    assert hasattr(result["conservative"], "confidence_score")
    assert hasattr(result["conservative"], "rationale")
    assert 0 <= result["conservative"].confidence_score <= 1
    assert hasattr(result["conservative"], "id")
    assert hasattr(result["conservative"], "review_id")
    assert hasattr(result["conservative"], "pubmed_result_id")
    assert hasattr(result["conservative"], "trace_id")
    assert hasattr(result["conservative"], "model_name")
    assert hasattr(result["conservative"], "screening_strategy")
    assert hasattr(result["conservative"], "start_time")
    assert hasattr(result["conservative"], "end_time")
    assert hasattr(result["conservative"], "response_metadata")

    # Check comprehensive reviewer result
    assert isinstance(result["comprehensive"], ScreeningResult)
    assert hasattr(result["comprehensive"], "decision")
    assert hasattr(result["comprehensive"], "confidence_score")
    assert hasattr(result["comprehensive"], "rationale")
    assert 0 <= result["comprehensive"].confidence_score <= 1
    assert hasattr(result["comprehensive"], "id")
    assert hasattr(result["comprehensive"], "review_id")
    assert hasattr(result["comprehensive"], "pubmed_result_id")
    assert hasattr(result["comprehensive"], "trace_id")
    assert hasattr(result["comprehensive"], "model_name")
    assert hasattr(result["comprehensive"], "screening_strategy")
    assert hasattr(result["comprehensive"], "start_time")
    assert hasattr(result["comprehensive"], "end_time")
    assert hasattr(result["comprehensive"], "response_metadata")


def test_screen_abstracts_batch(
    mock_pubmed_result, mock_systematic_review, monkeypatch
):
    """Test that screen_abstracts_batch correctly processes a batch of abstracts."""

    # Mock Streamlit session state
    class MockSessionState:
        def __init__(self):
            self.screen_abstracts_chain = screening_agents.screen_abstracts_chain
            self.screen_abstracts_results = []
            self.screen_abstracts_errors = []
            self.screen_abstracts_batch_idx = 0

    monkeypatch.setattr(st, "session_state", MockSessionState())

    # Create a batch of 2 identical results for testing
    batch = [mock_pubmed_result, mock_pubmed_result]
    batch_idx = 0

    # Run the batch screening
    result = screening_agents.screen_abstracts_batch(
        batch, batch_idx, mock_systematic_review
    )

    # Verify the output
    assert result is not None
    results, cb = result

    # Check batch results
    assert len(results) == 2
    for res in results:
        # Verify result structure
        assert hasattr(res, "search_result")
        assert hasattr(res, "conservative_result")
        assert hasattr(res, "comprehensive_result")

        # Check search result
        assert isinstance(res.search_result, PubMedResult)
        assert res.search_result.pmid == mock_pubmed_result.pmid

        # Check conservative result
        assert isinstance(res.conservative_result, ScreeningResult)
        assert hasattr(res.conservative_result, "decision")
        assert hasattr(res.conservative_result, "confidence_score")
        assert hasattr(res.conservative_result, "rationale")
        assert 0 <= res.conservative_result.confidence_score <= 1
        assert hasattr(res.conservative_result, "id")
        assert hasattr(res.conservative_result, "review_id")
        assert hasattr(res.conservative_result, "pubmed_result_id")
        assert hasattr(res.conservative_result, "trace_id")
        assert hasattr(res.conservative_result, "model_name")
        assert hasattr(res.conservative_result, "screening_strategy")
        assert hasattr(res.conservative_result, "start_time")
        assert hasattr(res.conservative_result, "end_time")
        assert hasattr(res.conservative_result, "response_metadata")

        # Check comprehensive result
        assert isinstance(res.comprehensive_result, ScreeningResult)
        assert hasattr(res.comprehensive_result, "decision")
        assert hasattr(res.comprehensive_result, "confidence_score")
        assert hasattr(res.comprehensive_result, "rationale")
        assert 0 <= res.comprehensive_result.confidence_score <= 1
        assert hasattr(res.comprehensive_result, "id")
        assert hasattr(res.comprehensive_result, "review_id")
        assert hasattr(res.comprehensive_result, "pubmed_result_id")
        assert hasattr(res.comprehensive_result, "trace_id")
        assert hasattr(res.comprehensive_result, "model_name")
        assert hasattr(res.comprehensive_result, "screening_strategy")
        assert hasattr(res.comprehensive_result, "start_time")
        assert hasattr(res.comprehensive_result, "end_time")
        assert hasattr(res.comprehensive_result, "response_metadata")

    # Verify OpenAI callback handler
    assert hasattr(cb, "total_tokens")
    assert hasattr(cb, "total_cost")
    assert cb.total_tokens > 0
