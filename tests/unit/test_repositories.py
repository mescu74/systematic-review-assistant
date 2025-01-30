from __future__ import annotations

from uuid import UUID

import pytest
from pydantic import BaseModel
from sr_assistant.core.models import AbstractScreeningResult, ScreeningDecisionType
from sr_assistant.core.repositories import ScreeningAbstractRepository
from sr_assistant.core.schemas import ScreeningResponse


class MockSupabaseResponse(BaseModel):
    data: list[dict]


class MockSupabaseTable:
    def __init__(self, data=None):
        self.data = data or []
        self.last_insert = None
        self.last_update = None
        self.conditions = []

    def insert(self, data):
        if "decision" in data:
            data = dict(data)  # Create a copy to avoid modifying the original
            data["decision"] = ScreeningDecisionType(data["decision"])
        self.last_insert = data
        return self

    def update(self, data):
        if "decision" in data:
            data = dict(data)  # Create a copy to avoid modifying the original
            data["decision"] = ScreeningDecisionType(data["decision"])
        self.last_update = data
        return self

    def select(self, *args):
        return self

    def eq(self, field, value):
        self.conditions.append((field, value))
        return self

    def delete(self):
        return self

    def execute(self):
        try:
            if self.last_insert:
                self.data.append(self.last_insert)
                return MockSupabaseResponse(data=[self.last_insert])
            if self.last_update:
                filtered_data = self.data
                for field, value in self.conditions:
                    filtered_data = [
                        item for item in filtered_data if str(item.get(field)) == str(value)
                    ]
                if filtered_data:
                    for item in filtered_data:
                        for key, value in self.last_update.items():
                            item[key] = value
                    return MockSupabaseResponse(data=[filtered_data[0]])
                return MockSupabaseResponse(data=[])
            filtered_data = self.data
            for field, value in self.conditions:
                filtered_data = [
                    item for item in filtered_data if str(item.get(field)) == str(value)
                ]
            return MockSupabaseResponse(data=filtered_data)
        finally:
            # Reset state after each operation
            self.last_insert = None
            self.last_update = None
            self.conditions = []


class MockSupabase:
    def __init__(self, data=None):
        self.mock_table = MockSupabaseTable(data)

    def table(self, _):
        return self.mock_table


@pytest.fixture
def mock_supabase():
    return MockSupabase()


@pytest.fixture
def repository(mock_supabase):
    return ScreeningAbstractRepository(mock_supabase)


def test_create_screening_result(repository):
    review_id = UUID("12345678-1234-5678-1234-567812345678")
    search_result_id = UUID("87654321-4321-8765-4321-876543210987")
    response = ScreeningResponse(
        decision="include",
        confidence_score=0.95,
        rationale="Test rationale",
        extracted_quotes=["quote1", "quote2"],
        exclusion_reason_categories=None,
    )

    result = repository.create_screening_result(review_id, search_result_id, response)

    assert isinstance(result, AbstractScreeningResult)
    assert result.review_id == review_id
    assert result.search_result_id == search_result_id
    assert result.decision == response.decision
    assert result.confidence_score == response.confidence_score
    assert result.rationale == response.rationale
    assert result.extracted_quotes == response.extracted_quotes
    assert result.exclusion_reason_categories == response.exclusion_reason_categories


def test_get_screening_results(repository):
    review_id = UUID("12345678-1234-5678-1234-567812345678")
    search_result_id = UUID("87654321-4321-8765-4321-876543210987")
    response = ScreeningResponse(
        decision="include",
        confidence_score=0.95,
        rationale="Test rationale",
        extracted_quotes=["quote1"],
        exclusion_reason_categories=None,
    )

    # Create a test result first
    repository.create_screening_result(review_id, search_result_id, response)

    # Get all results
    results = repository.get_screening_results(review_id)

    assert len(results) == 1
    assert isinstance(results[0], AbstractScreeningResult)
    assert results[0].review_id == review_id
    assert results[0].search_result_id == search_result_id


def test_get_screening_result(repository):
    review_id = UUID("12345678-1234-5678-1234-567812345678")
    search_result_id = UUID("87654321-4321-8765-4321-876543210987")
    response = ScreeningResponse(
        decision="include",
        confidence_score=0.95,
        rationale="Test rationale",
        extracted_quotes=["quote1"],
        exclusion_reason_categories=None,
    )

    # Create a test result first
    repository.create_screening_result(review_id, search_result_id, response)

    # Get specific result
    result = repository.get_screening_result(review_id, search_result_id)

    assert isinstance(result, AbstractScreeningResult)
    assert result.review_id == review_id
    assert result.search_result_id == search_result_id

    # Test non-existent result
    non_existent = repository.get_screening_result(
        UUID("00000000-0000-0000-0000-000000000000"),
        UUID("00000000-0000-0000-0000-000000000000"),
    )
    assert non_existent is None


def test_update_screening_result(repository):
    review_id = UUID("12345678-1234-5678-1234-567812345678")
    search_result_id = UUID("87654321-4321-8765-4321-876543210987")
    initial_response = ScreeningResponse(
        decision="include",
        confidence_score=0.95,
        rationale="Initial rationale",
        extracted_quotes=["quote1"],
        exclusion_reason_categories=None,
    )

    # Create initial result
    repository.create_screening_result(review_id, search_result_id, initial_response)

    # Update with new response
    updated_response = ScreeningResponse(
        decision="exclude",
        confidence_score=0.99,
        rationale="Updated rationale",
        extracted_quotes=["quote2"],
        exclusion_reason_categories=["Wrong study design type"],
    )

    updated_result = repository.update_screening_result(
        review_id, search_result_id, updated_response
    )

    assert isinstance(updated_result, AbstractScreeningResult)
    assert updated_result.review_id == review_id
    assert updated_result.search_result_id == search_result_id
    assert updated_result.decision == updated_response.decision
    assert updated_result.confidence_score == updated_response.confidence_score
    assert updated_result.rationale == updated_response.rationale
    assert updated_result.extracted_quotes == updated_response.extracted_quotes
    assert updated_result.exclusion_reason_categories == updated_response.exclusion_reason_categories
