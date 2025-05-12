import uuid

import pytest
from pydantic import ValidationError

from sr_assistant.core.schemas import SearchResultFilter
from sr_assistant.core.types import SearchDatabaseSource


def test_search_result_filter_empty_is_valid():
    """Test that an empty SearchResultFilter is valid (all fields optional)."""
    try:
        _ = SearchResultFilter()
    except ValidationError as e:
        pytest.fail(f"Empty SearchResultFilter should be valid, but raised: {e}")


def test_search_result_filter_all_fields_valid():
    """Test SearchResultFilter with all fields populated with valid types."""
    review_id = uuid.uuid4()
    data = {
        "review_id": review_id,
        "source_db": SearchDatabaseSource.PUBMED,
        "source_id": "12345",
        "doi": "10.1000/xyz123",
        "title": "Test Title",
        "year": "2023",
    }
    try:
        filt = SearchResultFilter(**data)  # type: ignore[arg-type]
        assert filt.review_id == review_id
        assert filt.source_db == SearchDatabaseSource.PUBMED
        assert filt.source_id == "12345"
        assert filt.doi == "10.1000/xyz123"
        assert filt.title == "Test Title"
        assert filt.year == "2023"
    except ValidationError as e:
        pytest.fail(f"SearchResultFilter with valid data raised: {e}")


def test_search_result_filter_partial_fields_valid():
    """Test SearchResultFilter with a subset of fields populated."""
    data = {
        "title": "Another Test",
        "year": "2024",
    }
    try:
        filt = SearchResultFilter(**data)  # type: ignore[arg-type]
        assert filt.title == "Another Test"
        assert filt.year == "2024"
        assert filt.review_id is None
        assert filt.doi is None
    except ValidationError as e:
        pytest.fail(f"SearchResultFilter with partial data raised: {e}")


def test_search_result_filter_invalid_type_for_field():
    """Test SearchResultFilter with an invalid type for a field (e.g., year as int)."""
    data = {
        "year": 2023  # Year should be str | None
    }
    with pytest.raises(ValidationError) as excinfo:
        SearchResultFilter(**data)  # type: ignore[arg-type]

    assert "year" in str(excinfo.value).lower()
    assert "string_type" in str(excinfo.value).lower()

    # Example of a clear type error Pydantic won't easily coerce for enums or UUIDs
    invalid_data_source_db = {"source_db": "INVALID_SOURCE"}
    with pytest.raises(ValidationError):
        SearchResultFilter(**invalid_data_source_db)  # type: ignore[arg-type]

    invalid_data_review_id = {"review_id": "not-a-uuid"}
    with pytest.raises(ValidationError):
        SearchResultFilter(**invalid_data_review_id)  # type: ignore[arg-type]
