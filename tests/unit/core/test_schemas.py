import uuid
from datetime import datetime, timedelta, timezone
from typing import Any

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


# --- SystematicReview Schemas Tests ---

from sr_assistant.core.schemas import (
    SystematicReviewCreate,
    SystematicReviewRead,
    SystematicReviewUpdate,
)
from sr_assistant.core.types import CriteriaFramework


@pytest.fixture
def valid_review_data() -> dict[str, Any]:
    return {
        "background": "Test Background",
        "research_question": "Test Question?",
        "criteria_framework": CriteriaFramework.PICO,
        "criteria_framework_answers": {"P": "Population", "I": "Intervention"},
        "inclusion_criteria": "Test Inclusion",
        "exclusion_criteria": "Test Exclusion",
        "review_metadata": {"key": "value"},
    }


def test_systematic_review_create_valid(valid_review_data: dict[str, Any]):
    """Test SystematicReviewCreate with valid data."""
    try:
        review = SystematicReviewCreate.model_validate(valid_review_data)
        assert review.research_question == "Test Question?"
        assert review.exclusion_criteria == "Test Exclusion"
        assert review.criteria_framework == CriteriaFramework.PICO
        assert review.criteria_framework_answers["P"] == "Population"
    except ValidationError as e:
        pytest.fail(f"SystematicReviewCreate validation failed for valid data: {e}")


def test_systematic_review_create_missing_required(valid_review_data: dict[str, Any]):
    """Test SystematicReviewCreate fails when required fields are missing."""
    data = valid_review_data.copy()
    del data["research_question"]
    with pytest.raises(ValidationError):
        SystematicReviewCreate.model_validate(data)

    data = valid_review_data.copy()
    del data["exclusion_criteria"]
    with pytest.raises(ValidationError):
        SystematicReviewCreate.model_validate(data)


def test_systematic_review_update_valid_partial(valid_review_data: dict[str, Any]):
    """Test SystematicReviewUpdate with partial valid data."""
    update_data = {"background": "Updated Background"}
    try:
        review_update = SystematicReviewUpdate.model_validate(update_data)
        assert review_update.background == "Updated Background"
        assert review_update.research_question is None  # Other fields should be None
    except ValidationError as e:
        pytest.fail(f"SystematicReviewUpdate validation failed for partial data: {e}")


def test_systematic_review_update_all_fields(valid_review_data: dict[str, Any]):
    """Test SystematicReviewUpdate with all fields valid."""
    try:
        review_update = SystematicReviewUpdate.model_validate(valid_review_data)
        assert review_update.research_question == "Test Question?"
    except ValidationError as e:
        pytest.fail(f"SystematicReviewUpdate validation failed for full data: {e}")


def test_systematic_review_read_valid(valid_review_data: dict[str, Any]):
    """Test SystematicReviewRead instantiation (simulating data from DB)."""
    read_data = {
        "id": uuid.uuid4(),
        "background": valid_review_data["background"],
        "research_question": valid_review_data["research_question"],
        "criteria_framework": valid_review_data["criteria_framework"],
        "criteria_framework_answers": valid_review_data["criteria_framework_answers"],
        "inclusion_criteria": valid_review_data["inclusion_criteria"],
        "exclusion_criteria": valid_review_data["exclusion_criteria"],
        "review_metadata": valid_review_data["review_metadata"],
        "created_at": datetime.now(timezone.utc),
        "updated_at": datetime.now(timezone.utc),
    }
    try:
        review_read = SystematicReviewRead.model_validate(read_data)
        assert review_read.id is not None
        assert review_read.created_at is not None
        assert review_read.updated_at is not None
        assert review_read.research_question == valid_review_data["research_question"]
    except ValidationError as e:
        pytest.fail(f"SystematicReviewRead validation failed for valid data: {e}")


def test_systematic_review_read_missing_optional_db_fields(
    valid_review_data: dict[str, Any],
):
    """Test SystematicReviewRead works even if optional DB fields (created_at) are None."""
    read_data = {
        "id": uuid.uuid4(),
        "background": valid_review_data["background"],
        "research_question": valid_review_data["research_question"],
        "criteria_framework": valid_review_data["criteria_framework"],
        "criteria_framework_answers": valid_review_data["criteria_framework_answers"],
        "inclusion_criteria": valid_review_data["inclusion_criteria"],
        "exclusion_criteria": valid_review_data["exclusion_criteria"],
        "review_metadata": valid_review_data["review_metadata"],
        "created_at": None,
        "updated_at": None,
    }
    try:
        review_read = SystematicReviewRead.model_validate(read_data)
        assert review_read.id is not None
        assert review_read.created_at is None
        assert review_read.updated_at is None
    except ValidationError as e:
        pytest.fail(
            f"SystematicReviewRead failed when optional DB fields are None: {e}"
        )


# --- ScreeningResult Schemas Tests ---

from sr_assistant.core.schemas import (
    ExclusionReasons,
    ScreeningResultCreate,
    ScreeningResultRead,
    ScreeningResultUpdate,
)
from sr_assistant.core.types import (
    ScreeningDecisionType,
    ScreeningStrategyType,
)


@pytest.fixture
def valid_screening_response_data() -> dict[str, Any]:
    """Provides base data for a ScreeningResponse component of ScreeningResult."""
    return {
        "decision": ScreeningDecisionType.INCLUDE,
        "confidence_score": 0.9,
        "rationale": "Clear evidence of inclusion.",
        "extracted_quotes": ["Quote 1", "Quote 2"],
        "exclusion_reason_categories": ExclusionReasons(
            population_exclusion_reasons=["Wrong age range for inclusion criteria"]
        ).model_dump(exclude_none=True),
    }


@pytest.fixture
def valid_screening_result_context_data() -> dict[str, Any]:
    """Provides contextual data for a ScreeningResult."""
    return {
        "id": uuid.uuid4(),
        "review_id": uuid.uuid4(),
        "trace_id": uuid.uuid4(),
        "model_name": "test_model",
        "screening_strategy": ScreeningStrategyType.CONSERVATIVE,
        "start_time": datetime.now(timezone.utc),
        "end_time": datetime.now(timezone.utc),
        "response_metadata": {"tokens": 100},
    }


def test_screening_result_create_valid(
    valid_screening_response_data: dict[str, Any],
    valid_screening_result_context_data: dict[str, Any],
):
    """Test ScreeningResultCreate with valid data."""
    create_data = {
        **valid_screening_response_data,
        **valid_screening_result_context_data,
    }
    try:
        result = ScreeningResultCreate.model_validate(create_data)
        assert result.decision == ScreeningDecisionType.INCLUDE
        assert result.id == valid_screening_result_context_data["id"]
    except ValidationError as e:
        pytest.fail(f"ScreeningResultCreate failed for valid data: {e}")


def test_screening_result_create_missing_required(
    valid_screening_response_data: dict[str, Any],
    valid_screening_result_context_data: dict[str, Any],
):
    """Test ScreeningResultCreate fails if required fields are missing."""
    full_data = {**valid_screening_response_data, **valid_screening_result_context_data}

    required_fields = [
        "decision",
        "confidence_score",
        "rationale",
        "id",
        "review_id",
        "model_name",
        "screening_strategy",
    ]
    for field in required_fields:
        data = full_data.copy()
        del data[field]
        with pytest.raises(ValidationError, match=field):
            ScreeningResultCreate.model_validate(data)


def test_screening_result_update_valid_partial():
    """Test ScreeningResultUpdate with partial valid data."""
    update_data = {"rationale": "Updated rationale."}
    try:
        result_update = ScreeningResultUpdate.model_validate(update_data)
        assert result_update.rationale == "Updated rationale."
        assert result_update.decision is None
    except ValidationError as e:
        pytest.fail(f"ScreeningResultUpdate failed for partial data: {e}")


def test_screening_result_read_valid(
    valid_screening_response_data: dict[str, Any],
    valid_screening_result_context_data: dict[str, Any],
):
    """Test ScreeningResultRead instantiation."""
    db_timestamps = {
        "created_at": datetime.now(timezone.utc),
        "updated_at": datetime.now(timezone.utc),
    }
    response_data_for_read = valid_screening_response_data.copy()
    read_data = {
        **response_data_for_read,
        **valid_screening_result_context_data,
        **db_timestamps,
    }
    assert "id" in read_data

    try:
        result_read = ScreeningResultRead.model_validate(read_data)
        assert result_read.id == valid_screening_result_context_data["id"]
        assert result_read.decision == ScreeningDecisionType.INCLUDE
        assert result_read.created_at is not None
        assert isinstance(result_read.exclusion_reason_categories, dict)
    except ValidationError as e:
        pytest.fail(f"ScreeningResultRead failed for valid data: {e}")


# --- ScreeningResolution Schemas Tests ---

from sr_assistant.core.schemas import (
    ResolverOutputSchema,
    ScreeningResolutionCreate,
    ScreeningResolutionRead,  # For LLM output part
)

# ScreeningDecisionType already imported
# AwareDatetime from pydantic.types already imported globally or earlier


@pytest.fixture
def valid_resolver_output_schema_data() -> dict[str, Any]:
    """Provides base data simulating LLM output for ResolverOutputSchema."""
    return {
        "resolver_decision": ScreeningDecisionType.INCLUDE,
        "resolver_reasoning": "The study meets all inclusion criteria based on abstract review.",
        "resolver_confidence_score": 0.95,
        # "resolver_include" is legacy/optional in ResolverOutputSchema, can omit for base valid data
    }


@pytest.fixture
def valid_resolution_context_data() -> dict[str, Any]:
    """Provides contextual data for creating a ScreeningResolution."""
    return {
        "search_result_id": uuid.uuid4(),
        "review_id": uuid.uuid4(),
        "conservative_result_id": uuid.uuid4(),
        "comprehensive_result_id": uuid.uuid4(),
        "resolver_model_name": "resolver_v1",
        "response_metadata": {"call_params": "default"},
        "start_time": datetime.now(timezone.utc) - timedelta(seconds=10),
        "end_time": datetime.now(timezone.utc),
        "trace_id": uuid.uuid4(),
    }


def test_screening_resolution_create_valid(
    valid_resolver_output_schema_data: dict[str, Any],
    valid_resolution_context_data: dict[str, Any],
):
    """Test ScreeningResolutionCreate with valid data."""
    create_data = {**valid_resolver_output_schema_data, **valid_resolution_context_data}
    try:
        resolution = ScreeningResolutionCreate.model_validate(create_data)
        assert resolution.resolver_decision == ScreeningDecisionType.INCLUDE
        assert (
            resolution.search_result_id
            == valid_resolution_context_data["search_result_id"]
        )
    except ValidationError as e:
        pytest.fail(f"ScreeningResolutionCreate failed for valid data: {e}")


def test_screening_resolution_read_valid(
    valid_resolver_output_schema_data: dict[str, Any],
    valid_resolution_context_data: dict[str, Any],
):
    """Test ScreeningResolutionRead instantiation."""
    db_timestamps = {
        "id": uuid.uuid4(),
        "created_at": datetime.now(timezone.utc),
        "updated_at": datetime.now(timezone.utc),
    }
    read_data = {
        **valid_resolver_output_schema_data,
        **valid_resolution_context_data,
        **db_timestamps,
    }
    try:
        resolution_read = ScreeningResolutionRead.model_validate(read_data)
        assert resolution_read.id == db_timestamps["id"]
        assert resolution_read.resolver_decision == ScreeningDecisionType.INCLUDE
        assert resolution_read.created_at is not None
    except ValidationError as e:
        pytest.fail(f"ScreeningResolutionRead failed for valid data: {e}")


# --- PicosSuggestions Schema Tests ---

from sr_assistant.core.schemas import PicosSuggestions


@pytest.fixture
def valid_picos_data() -> dict[str, str]:
    """Provides valid data for PicosSuggestions."""
    return {
        "population": "Patients with type 2 diabetes",
        "intervention": "Metformin 1000mg daily",
        "comparison": "Placebo",
        "outcome": "HbA1c levels",
        "general_critique": "Well-defined PICO components.",
    }


def test_picos_suggestions_valid(valid_picos_data: dict[str, str]):
    """Test PicosSuggestions with valid data."""
    try:
        picos = PicosSuggestions.model_validate(valid_picos_data)
        assert picos.population == "Patients with type 2 diabetes"
        assert picos.general_critique == "Well-defined PICO components."
    except ValidationError as e:
        pytest.fail(f"PicosSuggestions validation failed for valid data: {e}")


def test_picos_suggestions_missing_required(valid_picos_data: dict[str, str]):
    """Test PicosSuggestions fails if a required field is missing."""
    required_fields = list(valid_picos_data.keys())
    for field in required_fields:
        data = valid_picos_data.copy()
        del data[field]
        with pytest.raises(ValidationError, match=field):
            PicosSuggestions.model_validate(data)


@pytest.mark.unit
def test_resolver_output_schema_invalid_decision(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test ScreeningResolutionSchema with an invalid decision type."""
    invalid_data = valid_resolver_output_schema_data.copy()
    invalid_data["resolver_decision"] = "MAYBE_NOT"
    with pytest.raises(ValidationError):
        ResolverOutputSchema(**invalid_data)


@pytest.mark.unit
def test_resolver_output_schema_invalid_confidence(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test ScreeningResolutionSchema with an out-of-range confidence score."""
    # Test too high
    invalid_data_high = valid_resolver_output_schema_data.copy()
    invalid_data_high["resolver_confidence_score"] = 1.5
    with pytest.raises(ValidationError):
        ResolverOutputSchema(**invalid_data_high)

    # Test too low
    invalid_data_low = valid_resolver_output_schema_data.copy()
    invalid_data_low["resolver_confidence_score"] = -0.5
    with pytest.raises(ValidationError):
        ResolverOutputSchema(**invalid_data_low)


@pytest.mark.unit
def test_resolver_output_schema_missing_reasoning(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test ScreeningResolutionSchema with missing mandatory resolver_reasoning."""
    invalid_data = valid_resolver_output_schema_data.copy()
    del invalid_data["resolver_reasoning"]
    with pytest.raises(ValidationError):
        ResolverOutputSchema(**invalid_data)


# NOTE: The following tests (`test_resolver_output_schema_resolver_include_valid`
#       and `test_resolver_output_schema_resolver_include_invalid_type`) are
#       commented out because they test the legacy `resolver_include` field
#       which has been renamed to `contributing_strategies`. The tests for
#       `contributing_strategies`
#       (`test_resolver_output_schema_contributing_strategies_valid` and
#       `test_resolver_output_schema_contributing_strategies_invalid_type`)
#       supersede these.
# @pytest.mark.unit
# def test_resolver_output_schema_resolver_include_valid(
#     valid_resolver_output_schema_data: dict[str, Any],
# ):
#     """Test ScreeningResolutionSchema with valid 'resolver_include' values."""
#     data = valid_resolver_output_schema_data.copy()
#     data["resolver_include"] = ["conservative", "comprehensive"]
#     loaded = ResolverOutputSchema(**data)
#     assert loaded.resolver_include == ["conservative", "comprehensive"]
#
#     data["resolver_include"] = []
#     loaded = ResolverOutputSchema(**data)
#     assert loaded.resolver_include == []
#
#
# @pytest.mark.unit
# def test_resolver_output_schema_resolver_include_invalid_type(
#     valid_resolver_output_schema_data: dict[str, Any],
# ):
#     """Test ScreeningResolutionSchema with invalid type for 'resolver_include' values."""
#     data = valid_resolver_output_schema_data.copy()
#     data["resolver_include"] = ["conservative", 123]  # 123 is not a string
#     with pytest.raises(ValidationError):
#         ResolverOutputSchema(**data)


@pytest.mark.unit
def test_resolver_output_schema_contributing_strategies_valid(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test ResolverOutputSchema with valid 'contributing_strategies' values."""
    data = valid_resolver_output_schema_data.copy()
    # Test with valid strategy strings (Pydantic will convert them to enums)
    data["contributing_strategies"] = ["conservative", "comprehensive"]
    loaded = ResolverOutputSchema.model_validate(data)
    assert loaded.contributing_strategies == [
        ScreeningStrategyType.CONSERVATIVE,
        ScreeningStrategyType.COMPREHENSIVE,
    ]

    # Test with empty list
    data["contributing_strategies"] = []
    loaded = ResolverOutputSchema.model_validate(data)
    assert loaded.contributing_strategies == []

    # Test with actual enum members
    data["contributing_strategies"] = [ScreeningStrategyType.CONSERVATIVE]
    loaded = ResolverOutputSchema.model_validate(data)
    assert loaded.contributing_strategies == [ScreeningStrategyType.CONSERVATIVE]


@pytest.mark.unit
def test_resolver_output_schema_contributing_strategies_invalid_type(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test ResolverOutputSchema with invalid type for 'contributing_strategies' values."""
    data = valid_resolver_output_schema_data.copy()
    # Test with a non-string, non-enum member in the list
    data["contributing_strategies"] = ["conservative", 123]
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(data)

    # Test with a string that is not a valid enum member value
    data["contributing_strategies"] = ["conservative", "not_a_strategy"]
    with pytest.raises(ValidationError):
        ResolverOutputSchema.model_validate(data)


# Test if fields intended to be populated by caller can be None (if defaults allow)
@pytest.mark.unit
def test_resolver_output_schema_caller_populated_fields_default(
    valid_resolver_output_schema_data: dict[str, Any],
):
    """Test that caller-populated fields in ScreeningResolutionSchema can default to None if allowed by schema."""
    # Minimal data, relying on defaults for caller-populated fields in ScreeningResolutionSchema
    schema = ResolverOutputSchema(**valid_resolver_output_schema_data)
    assert schema.review_id is None
    assert schema.search_result_id is None
