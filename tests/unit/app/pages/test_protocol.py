# tests/unit/app/pages/test_protocol.py
"""Unit tests for protocol.py page logic."""

import uuid
from collections.abc import Generator
from typing import Any
from unittest.mock import MagicMock

import pytest
import streamlit as st
from pytest import MonkeyPatch

from sr_assistant.app import services  # For services.ReviewService

# Import functions and classes to test
from sr_assistant.app.pages.protocol import (
    build_review_model_from_pico,
    persist_review,
)
from sr_assistant.core import (
    schemas as app_schemas,  # For app_schemas.SystematicReviewCreate etc.
)

# Import models, schemas, and services needed for mocking and type checking
from sr_assistant.core.models import CriteriaFramework, SystematicReview


@pytest.fixture
def mock_session_state_with_pico(
    monkeypatch: MonkeyPatch,
) -> Generator[dict[str, Any], None, None]:
    """Fixture to mock st.session_state with PICO data."""
    mock_state = {
        "review_id": uuid.uuid4(),
        "background": "Test background about disease X.",
        "research_question": "Does drug Y work for disease X?",
        "pico_population": "Adults with disease X",
        "pico_intervention": "Drug Y 50mg",
        "pico_comparison": "Placebo",
        "pico_outcome": "Symptom reduction",
        "exclusion_criteria": "Pregnant women",
    }
    for key, value in mock_state.items():
        monkeypatch.setattr(st.session_state, key, value, raising=False)
    mock_session_factory = MagicMock()
    monkeypatch.setattr(
        st.session_state, "session_factory", mock_session_factory, raising=False
    )
    yield mock_state
    for key in mock_state:
        monkeypatch.delattr(st.session_state, key, raising=False)
    monkeypatch.delattr(st.session_state, "session_factory", raising=False)


# --- Tests for build_review_model_from_pico ---


def test_build_review_model_from_pico_all_fields(
    mock_session_state_with_pico: dict[str, Any],
):
    """Test building model when all PICO fields are present."""
    mock_response_widget = MagicMock()
    review = build_review_model_from_pico(mock_response_widget)

    assert isinstance(review, SystematicReview)
    assert review.id == st.session_state.review_id
    assert review.background == "Test background about disease X."
    assert review.research_question == "Does drug Y work for disease X?"
    assert (
        "Population: Adults with disease X" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert (
        "Intervention: Drug Y 50mg" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert (
        "Comparison: Placebo" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert (
        "Outcome: Symptom reduction" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert review.exclusion_criteria == "Pregnant women"
    assert review.criteria_framework == CriteriaFramework.PICO
    assert review.criteria_framework_answers == {
        "population": "Adults with disease X",
        "intervention": "Drug Y 50mg",
        "comparison": "Placebo",
        "outcome": "Symptom reduction",
    }
    mock_response_widget.success.assert_called_once()


def test_build_review_model_from_pico_missing_fields(
    mock_session_state_with_pico: dict[str, Any],
):
    """Test building model when some PICO fields are missing."""
    st.session_state.pico_comparison = ""
    st.session_state.pico_outcome = ""
    st.session_state.exclusion_criteria = ""

    mock_response_widget = MagicMock()
    review = build_review_model_from_pico(mock_response_widget)

    assert isinstance(review, SystematicReview)
    assert (
        "Population: Adults with disease X" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert (
        "Intervention: Drug Y 50mg" in review.inclusion_criteria
        if review.inclusion_criteria
        else False
    )
    assert (
        "Comparison:" not in review.inclusion_criteria
        if review.inclusion_criteria is not None
        else True
    )
    assert (
        "Outcome:" not in review.inclusion_criteria
        if review.inclusion_criteria is not None
        else True
    )
    assert review.exclusion_criteria == "Not specified"
    assert review.criteria_framework == CriteriaFramework.PICO
    assert review.criteria_framework_answers == {
        "population": "Adults with disease X",
        "intervention": "Drug Y 50mg",
        "comparison": "",
        "outcome": "",
    }
    mock_response_widget.success.assert_called_once()


# --- Tests for persist_review ---


def test_persist_review_new_record(
    mock_session_state_with_pico: dict[str, Any], mocker: MagicMock
):
    """Test saving a new review with PICO data using ReviewService."""
    review_to_save = build_review_model_from_pico(MagicMock())

    mock_review_service_instance = MagicMock(spec=services.ReviewService)
    mock_review_service_instance.get_review.return_value = None
    mock_review_service_instance.create_review.return_value = review_to_save

    mocker.patch(
        "sr_assistant.app.pages.protocol.ReviewService",
        return_value=mock_review_service_instance,
    )

    persisted_review = persist_review(review_to_save)

    mock_review_service_instance.get_review.assert_called_once_with(review_to_save.id)
    mock_review_service_instance.create_review.assert_called_once()
    call_args = mock_review_service_instance.create_review.call_args[0][0]
    assert isinstance(call_args, app_schemas.SystematicReviewCreate)
    assert call_args.research_question == review_to_save.research_question
    assert persisted_review is review_to_save


def test_persist_review_update_record(
    mock_session_state_with_pico: dict[str, Any], mocker: MagicMock
):
    """Test updating an existing review with PICO data using ReviewService."""
    existing_id = st.session_state.review_id
    updated_model_from_pico = build_review_model_from_pico(MagicMock())
    updated_model_from_pico.id = existing_id

    mock_existing_review_in_db = SystematicReview(
        id=existing_id,
        research_question="Old question?",
        exclusion_criteria="Old exclusion",
    )

    mock_review_service_instance = MagicMock(spec=services.ReviewService)
    mock_review_service_instance.get_review.return_value = mock_existing_review_in_db
    mock_review_service_instance.update_review.return_value = updated_model_from_pico

    mocker.patch(
        "sr_assistant.app.pages.protocol.ReviewService",
        return_value=mock_review_service_instance,
    )

    persisted_review = persist_review(updated_model_from_pico)

    mock_review_service_instance.get_review.assert_called_once_with(existing_id)
    mock_review_service_instance.update_review.assert_called_once()

    update_call_args = mock_review_service_instance.update_review.call_args[1]
    assert update_call_args["review_id"] == existing_id
    update_payload = update_call_args["review_update_data"]
    assert isinstance(update_payload, app_schemas.SystematicReviewUpdate)
    assert update_payload.research_question == updated_model_from_pico.research_question

    assert persisted_review is updated_model_from_pico
