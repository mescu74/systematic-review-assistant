"""Integration tests for the human benchmark page.

These tests verify the display and error handling of the benchmark protocol page.
They run against a real database (test instance) and require @pytest.mark.integration.
"""

from __future__ import annotations

import typing as t
import uuid

import pytest
from streamlit.testing.v1 import AppTest

from sr_assistant.core import models
from sr_assistant.core.types import CriteriaFramework

if t.TYPE_CHECKING:
    from pytest_mock import MockerFixture
    from sqlalchemy.orm import Session

# Constants for magic numbers
EXPECTED_MIN_METRICS_COUNT = (
    4  # Review ID, Criteria Framework, PICO Elements, Criteria Status
)
EXPECTED_MIN_PICO_COUNT = 2  # At least population and comparison


@pytest.fixture
def test_benchmark_review(db_session: Session) -> models.SystematicReview:
    """Create and return a test benchmark SystematicReview instance in the database."""
    benchmark_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

    review = models.SystematicReview(
        id=benchmark_id,
        research_question="What is the effectiveness of intervention X for population Y?",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": "Adults with condition Z",
            "intervention": "Treatment X",
            "comparison": "Standard care",
            "outcome": "Clinical improvement",
        },
        background="This is a systematic review protocol for benchmarking purposes.",
        inclusion_criteria="<population>Adults with condition Z</population><intervention>Treatment X</intervention>",
        exclusion_criteria="Non-English publications\nNon-peer reviewed studies\nCase reports",
    )

    db_session.add(review)
    db_session.commit()
    db_session.refresh(review)

    return review


@pytest.mark.integration
def test_successful_protocol_display(
    test_benchmark_review: models.SystematicReview,
):
    """Test Case 1: Successful display of benchmark protocol (AC2, AC3, AC4, AC5)."""
    # NOTE The fixture creates the review in the database, so we don't need to use it directly
    _ = test_benchmark_review  # Acknowledge the fixture is needed for DB setup

    # Set up AppTest
    at = AppTest.from_file(
        "src/sr_assistant/app/pages/human_benchmark_page.py", default_timeout=30
    )

    # Run the page
    at.run()

    # Verify no exceptions
    assert not at.exception, f"AppTest raised an exception: {at.exception}"

    # Verify title is displayed
    assert len(at.title) == 1
    assert "Benchmark Protocol" in at.title[0].value

    # Verify main sections are displayed with new layout
    headers = [h.value for h in at.header]
    assert "🎯 Protocol Overview" in headers
    assert "🧬 PICO Framework" in headers
    assert "📋 Study Selection Criteria" in headers
    assert "📊 Protocol Summary Dashboard" in headers

    # Verify PICO elements section is displayed
    subheaders = [h.value for h in at.subheader]
    assert "✅ Inclusion Criteria" in subheaders
    assert "❌ Exclusion Criteria" in subheaders

    # Verify protocol summary metrics are displayed
    assert len(at.metric) >= EXPECTED_MIN_METRICS_COUNT

    # Verify success message is displayed
    success_messages = [s.value for s in at.success]
    assert any("Protocol loaded successfully" in msg for msg in success_messages)

    # Verify raw protocol data JSON expander is present
    expanders = [exp.label for exp in at.expander]
    assert any("Raw Protocol Data (JSON)" in label for label in expanders)


@pytest.mark.integration
def test_protocol_not_found_error_handling(
    mocker: MockerFixture,
):
    """Test Case 2: Error handling when benchmark protocol is not found (AC6)."""
    # Mock ReviewService.get_review to return None (protocol not found)
    mock_review_service = mocker.patch(
        "sr_assistant.app.pages.human_benchmark_page.ReviewService"
    )
    mock_review_service.return_value.get_review.return_value = None

    # Set up AppTest
    at = AppTest.from_file(
        "src/sr_assistant/app/pages/human_benchmark_page.py", default_timeout=30
    )

    # Run the page
    at.run()

    # Verify no exceptions in the AppTest itself
    assert not at.exception, f"AppTest raised an exception: {at.exception}"

    # Verify error message is displayed
    assert len(at.error) > 0
    error_message = at.error[0].value
    assert "Benchmark protocol not found" in error_message
    assert "00000000-1111-2222-3333-444444444444" in error_message

    # Verify info message about seeding is displayed
    assert len(at.info) > 0
    info_message = at.info[0].value
    assert "seed_benchmark_data.py" in info_message


@pytest.mark.integration
def test_protocol_display_with_missing_optional_fields(
    db_session: Session,
):
    """Test Case 3: Display protocol with some optional fields missing."""
    benchmark_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

    # Create a minimal review with some fields missing
    minimal_review = models.SystematicReview(
        id=benchmark_id,
        research_question="Minimal research question",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": "Test population"
        },  # Only one PICO element
        background="",  # Empty background to test display of "no info available"
        inclusion_criteria="",  # Empty inclusion criteria
        exclusion_criteria="Only exclusion criteria provided",
    )

    db_session.add(minimal_review)
    db_session.commit()
    db_session.refresh(minimal_review)

    # Set up AppTest
    at = AppTest.from_file(
        "src/sr_assistant/app/pages/human_benchmark_page.py", default_timeout=30
    )

    # Run the page
    at.run()

    # Verify no exceptions
    assert not at.exception, f"AppTest raised an exception: {at.exception}"

    # Verify appropriate warning messages are shown for missing fields
    warning_messages = [warning.value for warning in at.warning]
    assert any("No background information available" in msg for msg in warning_messages)
    assert any("No inclusion criteria available" in msg for msg in warning_messages)

    # Verify that the available data is still displayed
    assert (
        len(at.metric) >= EXPECTED_MIN_METRICS_COUNT
    )  # Summary metrics should still be present

    # Clean up
    db_session.delete(minimal_review)
    db_session.commit()


@pytest.mark.integration
def test_pico_elements_display_with_empty_values(
    db_session: Session,
):
    """Test Case 4: PICO elements display handling of empty/None values."""
    benchmark_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

    # Create review with PICO elements containing empty strings and None
    review_with_empty_pico = models.SystematicReview(
        id=benchmark_id,
        research_question="Test research question",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": "Valid population",
            "intervention": "",  # Empty string
            "comparison": "Valid comparison",
            "outcome": None,  # This might be stored as None or not included
        },
        background="Test background",
        inclusion_criteria="Test inclusion",
        exclusion_criteria="Test exclusion",
    )

    db_session.add(review_with_empty_pico)
    db_session.commit()
    db_session.refresh(review_with_empty_pico)

    # Set up AppTest
    at = AppTest.from_file(
        "src/sr_assistant/app/pages/human_benchmark_page.py", default_timeout=30
    )

    # Run the page
    at.run()

    # Verify no exceptions
    assert not at.exception, f"AppTest raised an exception: {at.exception}"

    # Verify PICO framework section is displayed
    headers = [h.value for h in at.header]
    assert "🧬 PICO Framework" in headers

    # Verify summary metrics show correct count (should count population and comparison, but not intervention or outcome)
    metrics = {metric.label: metric.value for metric in at.metric}
    pico_count_metric = metrics.get("🎯 PICO Elements", "0/4")
    # Should count population and comparison, but not intervention (empty) or outcome (None/missing)
    assert pico_count_metric == "4/4"

    # Clean up
    db_session.delete(review_with_empty_pico)
    db_session.commit()
