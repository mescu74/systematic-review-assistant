from __future__ import annotations

from uuid import UUID

import pytest
from sr_assistant.core.models import AbstractScreeningResult
from sr_assistant.core.repositories import ScreeningAbstractRepository
from sr_assistant.core.schemas import ScreeningResponse
from supabase import Client


@pytest.fixture
def repository(supabase_client: Client) -> ScreeningAbstractRepository:
    return ScreeningAbstractRepository(supabase_client)


@pytest.mark.integration
def test_screening_result_lifecycle(repository: ScreeningAbstractRepository):
    """Test the complete lifecycle of a screening result."""
    # Setup test data
    review_id = UUID("12345678-1234-5678-1234-567812345678")
    search_result_id = UUID("87654321-4321-8765-4321-876543210987")
    response = ScreeningResponse(
        decision="include",
        confidence_score=0.95,
        rationale="Integration test rationale",
        extracted_quotes=["quote1", "quote2"],
        exclusion_reason_categories=[],
    )

    # Test creation
    result = repository.create_screening_result(review_id, search_result_id, response)
    assert isinstance(result, AbstractScreeningResult)
    assert result.review_id == review_id
    assert result.search_result_id == search_result_id
    assert result.decision == response.decision

    # Test retrieval by ID
    retrieved = repository.get_screening_result(review_id, search_result_id)
    assert retrieved is not None
    assert retrieved.review_id == review_id
    assert retrieved.search_result_id == search_result_id
    assert retrieved.decision == response.decision

    # Test retrieval of all results
    all_results = repository.get_screening_results(review_id)
    assert len(all_results) >= 1
    assert any(r.search_result_id == search_result_id for r in all_results)

    # Test update
    updated_response = ScreeningResponse(
        decision="exclude",
        confidence_score=0.99,
        rationale="Updated integration test rationale",
        extracted_quotes=["updated_quote"],
        exclusion_reason_categories=["category1"],
    )
    updated = repository.update_screening_result(
        review_id, search_result_id, updated_response
    )
    assert updated.decision == "exclude"
    assert updated.rationale == "Updated integration test rationale"

    # Verify update
    final = repository.get_screening_result(review_id, search_result_id)
    assert final is not None
    assert final.decision == "exclude"
    assert final.rationale == "Updated integration test rationale"


@pytest.mark.integration
def test_get_nonexistent_screening_result(repository: ScreeningAbstractRepository):
    """Test retrieving a non-existent screening result."""
    result = repository.get_screening_result(
        UUID("00000000-0000-0000-0000-000000000000"),
        UUID("00000000-0000-0000-0000-000000000000"),
    )
    assert result is None


@pytest.mark.integration
def test_get_screening_results_empty(repository: ScreeningAbstractRepository):
    """Test retrieving screening results for a review with no results."""
    results = repository.get_screening_results(
        UUID("00000000-0000-0000-0000-000000000000")
    )
    assert isinstance(results, list)
    assert len(results) == 0
