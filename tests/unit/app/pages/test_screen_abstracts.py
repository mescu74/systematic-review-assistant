"""Unit tests for the screen_abstracts.py Streamlit page."""

from __future__ import annotations

import typing as t
import uuid
from collections.abc import KeysView
from datetime import UTC, datetime
from unittest.mock import MagicMock

import pytest
from pytest_mock import MockerFixture
from sqlmodel import Session

# import streamlit as st # Page imports st, test file doesn't need to if st is mocked via string path
from streamlit.testing.v1 import AppTest

from sr_assistant.app import services as app_services_module
from sr_assistant.app.agents.screening_agents import (
    ScreenAbstractResultTuple,
    ScreeningError,
)
from sr_assistant.app.services import ScreeningService
from sr_assistant.core import models, repositories, schemas, types
from sr_assistant.core.repositories import RecordNotFoundError


# Copied and adapted from py/testing-streamlit-pages rule
class MockSessionState:
    def __init__(
        self, backing_dict_ref: dict[str, t.Any], mocker_instance: MockerFixture
    ):
        object.__setattr__(self, "_backing_dict", backing_dict_ref)
        object.__setattr__(
            self, "_mocker", mocker_instance
        )  # Storing mocker for __getattr__ if needed

    def __getattr__(self, name: str) -> t.Any:
        if name in self._backing_dict:
            return self._backing_dict[name]
        # Avoid creating mocks for dunder methods that MagicMock itself might need
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(
                f"Dunder attribute {name} not found in MockSessionState backing dict and not auto-mocked."
            )
        # If an attribute is accessed that isn't explicitly set, create a new MagicMock.
        new_mock = self._mocker.MagicMock(name=f"session_state_auto_attr_{name}")
        self._backing_dict[name] = new_mock
        return new_mock

    @t.override
    def __setattr__(self, name: str, value: t.Any) -> None:
        self._backing_dict[name] = value

    @t.override
    def __delattr__(self, name: str) -> None:
        if name in self._backing_dict:
            del self._backing_dict[name]
        else:
            # Behavior for deleting non-existent key can be to raise AttributeError or pass silently
            # Streamlit's actual session state might raise error or allow it.
            # For testing, raising an error is often more informative.
            raise AttributeError(
                f"MockSessionState: cannot delete non-existent attribute '{name}'"
            )

    def __contains__(self, key: str) -> bool:
        return key in self._backing_dict

    def get(self, key: str, default: t.Any = None) -> t.Any:
        return self._backing_dict.get(key, default)

    def keys(self) -> KeysView[str]:
        return self._backing_dict.keys()

    def clear(self) -> None:
        self._backing_dict.clear()

    # Add __setitem__ and __getitem__ for dict-like access st.session_state["key"]
    def __setitem__(self, key: str, value: t.Any) -> None:
        self._backing_dict[key] = value

    def __getitem__(self, key: str) -> t.Any:
        if key in self._backing_dict:
            return self._backing_dict[key]
        # Match getattr behavior: if key not present, mock it on demand or raise KeyError
        # For __getitem__, raising KeyError is more standard dict behavior
        raise KeyError(key)


@pytest.fixture
def app_test_env_v2(
    mocker: MockerFixture,
) -> tuple[
    AppTest, MagicMock, MagicMock, list[models.SearchResult], models.SystematicReview
]:
    """Fixture using MockSessionState for st.session_state."""
    mock_st_in_page = mocker.patch("sr_assistant.app.pages.screen_abstracts.st")

    session_state_backing_dict: dict[str, t.Any] = {}
    mock_session_state_instance = MockSessionState(session_state_backing_dict, mocker)
    mock_st_in_page.session_state = mock_session_state_instance

    mock_review_repo_inst = mocker.MagicMock(
        spec=repositories.SystematicReviewRepository
    )
    mock_search_repo_inst = mocker.MagicMock(spec=repositories.SearchResultRepository)
    mock_screen_abstract_repo_inst = mocker.MagicMock(
        spec=repositories.ScreenAbstractResultRepository
    )
    mock_screen_resolution_repo_inst = mocker.MagicMock(
        spec=repositories.ScreeningResolutionRepository
    )
    mock_screening_service_instance = mocker.MagicMock(spec=ScreeningService)

    test_review_id = uuid.uuid4()
    mock_review_model = models.SystematicReview(
        id=test_review_id,
        research_question="Test RQ for Screening Page",
        exclusion_criteria="Test Excl Criteria",
    )
    mock_review_repo_inst.get_by_id.return_value = mock_review_model

    sr1_id = uuid.uuid4()
    sr2_id = uuid.uuid4()
    mock_sr1 = models.SearchResult(
        id=sr1_id,
        review_id=test_review_id,
        title="SR1 Title",
        source_db=types.SearchDatabaseSource.PUBMED,
        source_id="pmid1",
        year="2023",
    )
    mock_sr2 = models.SearchResult(
        id=sr2_id,
        review_id=test_review_id,
        title="SR2 Title",
        source_db=types.SearchDatabaseSource.PUBMED,
        source_id="pmid2",
        year="2024",
    )
    mock_search_repo_inst.get_by_review_id.return_value = [mock_sr1, mock_sr2]

    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.init_review_repository",
        return_value=mock_review_repo_inst,
    )
    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.init_pubmed_repository",
        return_value=mock_search_repo_inst,
    )
    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.init_screen_abstracts_repository",
        return_value=mock_screen_abstract_repo_inst,
    )
    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.init_screen_resolution_repository",
        return_value=mock_screen_resolution_repo_inst,
    )
    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.init_screening_service",
        return_value=mock_screening_service_instance,
    )

    mocker.patch("sr_assistant.app.pages.screen_abstracts.logger")

    # Mock the conflict resolution function to avoid API calls and timeouts
    mock_resolver = mocker.MagicMock()
    mock_resolver.return_value = None  # Return None to simulate no resolution
    mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.invoke_resolver_chain", mock_resolver
    )

    mock_session_factory_patch = mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.session_factory"
    )
    mock_db_session_instance = mocker.MagicMock(spec=Session)
    mock_session_factory_patch.return_value.__enter__.return_value = (
        mock_db_session_instance
    )

    at = AppTest.from_file(
        "src/sr_assistant/app/pages/screen_abstracts.py", default_timeout=10
    )

    # Set up AppTest session state - this is the proper way to do it with AppTest
    at.session_state.review_id = test_review_id
    at.session_state.review = mock_review_model
    at.session_state.search_results = [mock_sr1, mock_sr2]
    at.session_state.repo_review = mock_review_repo_inst
    at.session_state.search_repo = mock_search_repo_inst
    at.session_state.repo_screen_abstracts = mock_screen_abstract_repo_inst
    at.session_state.repo_screen_resolution = mock_screen_resolution_repo_inst
    at.session_state.screening_service = (
        mock_screening_service_instance  # Important: Set the mocked service
    )

    # Also populate the MockSessionState instance for backwards compatibility with existing assertions
    mock_st_in_page.session_state.review_id = test_review_id
    mock_st_in_page.session_state.review = mock_review_model
    mock_st_in_page.session_state.search_results = [mock_sr1, mock_sr2]
    mock_st_in_page.session_state.repo_review = mock_review_repo_inst
    mock_st_in_page.session_state.search_repo = mock_search_repo_inst
    mock_st_in_page.session_state.repo_screen_abstracts = mock_screen_abstract_repo_inst
    mock_st_in_page.session_state.repo_screen_resolution = (
        mock_screen_resolution_repo_inst
    )

    return (
        at,
        mock_st_in_page,
        mock_screening_service_instance,
        [mock_sr1, mock_sr2],
        mock_review_model,
    )


class TestScreenAbstractsPage:
    def test_page_loads_and_start_button_present(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        at, _, _, _, _ = app_test_env_v2

        at.run()

        assert not at.exception

        # Use AppTest assertions instead of mock assertions
        assert len(at.title) == 1
        assert at.title[0].value == "Abstract Screening"

        # Check for the start button using AppTest button assertions
        start_button_found = False
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button_found = True
                break
        assert start_button_found, "Start Abstract Screening button not found"

    def test_start_screening_calls_service_with_correct_ids(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        (at, _, mock_screening_service, mock_search_results_list, mock_review_model) = (
            app_test_env_v2
        )

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Use AppTest button interaction instead of mocking st.button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception
        expected_sr_ids = [sr.id for sr in mock_search_results_list]
        mock_screening_service.perform_batch_abstract_screening.assert_called_once_with(
            review_id=mock_review_model.id, search_result_ids_to_screen=expected_sr_ids
        )

    def test_results_processing_from_service_success(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
        mocker: MockerFixture,
    ):
        (at, _, mock_screening_service, mock_search_results_list, mock_review_model) = (
            app_test_env_v2
        )

        kons_run_id1 = uuid.uuid4()
        comp_run_id1 = uuid.uuid4()
        kons_run_id2 = uuid.uuid4()
        comp_run_id2 = uuid.uuid4()
        sr1 = mock_search_results_list[0]
        sr2 = mock_search_results_list[1]

        mock_kons_res1 = schemas.ScreeningResult(
            id=kons_run_id1,
            review_id=mock_review_model.id,
            search_result_id=sr1.id,
            decision=types.ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="R1C",
            screening_strategy=types.ScreeningStrategyType.CONSERVATIVE,
            model_name="m1",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
        )
        mock_comp_res1 = schemas.ScreeningResult(
            id=comp_run_id1,
            review_id=mock_review_model.id,
            search_result_id=sr1.id,
            decision=types.ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="R1P",
            screening_strategy=types.ScreeningStrategyType.COMPREHENSIVE,
            model_name="m1",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
        )
        mock_kons_res2 = schemas.ScreeningResult(
            id=kons_run_id2,
            review_id=mock_review_model.id,
            search_result_id=sr2.id,
            decision=types.ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="R2C",
            screening_strategy=types.ScreeningStrategyType.CONSERVATIVE,
            model_name="m2",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
            exclusion_reason_categories=schemas.ExclusionReasons(
                population_exclusion_reasons=["Wrong age range for inclusion criteria"]
            ),
        )
        mock_comp_res2 = schemas.ScreeningResult(
            id=comp_run_id2,
            review_id=mock_review_model.id,
            search_result_id=sr2.id,
            decision=types.ScreeningDecisionType.EXCLUDE,  # Same as conservative to avoid conflicts
            confidence_score=0.7,
            rationale="R2P",
            screening_strategy=types.ScreeningStrategyType.COMPREHENSIVE,
            model_name="m2",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
        )

        service_return_tuples = [
            ScreenAbstractResultTuple(sr1, mock_kons_res1, mock_comp_res1),
            ScreenAbstractResultTuple(sr2, mock_kons_res2, mock_comp_res2),
        ]
        mock_screening_service.perform_batch_abstract_screening.return_value = (
            service_return_tuples
        )

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Find and click the start button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception

        # Use AppTest session state assertions instead of mock assertions
        assert len(at.session_state.screen_abstracts_results) == 4
        assert at.session_state.screen_abstracts_screened == 2
        assert (
            at.session_state.screen_abstracts_included == 1
        )  # SR1 is included (both strategies agree)
        assert (
            at.session_state.screen_abstracts_excluded == 1
        )  # SR2 is excluded (both strategies agree)
        assert (
            at.session_state.screen_abstracts_uncertain == 0
        )  # No uncertain decisions
        assert (
            len(at.session_state.screen_abstracts_conflicts) == 0
        )  # No conflicts since decisions match

        # Check that success message is shown in the UI
        assert len(at.success) > 0, "Success message not found in UI"

    def test_results_processing_with_screening_errors(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        at, _, mock_screening_service, mock_search_results_list, _ = app_test_env_v2
        sr1 = mock_search_results_list[0]
        mock_error = ScreeningError(
            search_result=sr1, error="Test error", message="Something went wrong"
        )
        service_return_tuples = [ScreenAbstractResultTuple(sr1, mock_error, mock_error)]
        mock_screening_service.perform_batch_abstract_screening.return_value = (
            service_return_tuples
        )

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Find and click the start button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception
        assert len(at.session_state.screen_abstracts_errors) == 2
        assert at.session_state.screen_abstracts_errors[0] == mock_error
        assert at.session_state.screen_abstracts_screened == 1
        assert at.session_state.screen_abstracts_included == 0
        assert at.session_state.screen_abstracts_excluded == 0
        assert at.session_state.screen_abstracts_uncertain == 0

    def test_service_raises_record_not_found_error(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        at, _, mock_screening_service, _, _ = app_test_env_v2
        mock_screening_service.perform_batch_abstract_screening.side_effect = (
            RecordNotFoundError("Review not found by service")
        )

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Find and click the start button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception

        # Check that error status is shown in the UI
        error_found = False
        for status in at.status:
            if "Review not found by service" in str(status.label):
                error_found = True
                break
        assert error_found, "Error status message not found in UI"

    def test_service_raises_generic_service_error(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        at, _, mock_screening_service, _, _ = app_test_env_v2
        mock_screening_service.perform_batch_abstract_screening.side_effect = (
            app_services_module.ServiceError("Generic service failure")
        )

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Find and click the start button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception

        # Check that error status is shown in the UI
        error_found = False
        for status in at.status:
            if "Generic service failure" in str(status.label):
                error_found = True
                break
        assert error_found, "Error status message not found in UI"

    def test_no_search_results_to_screen(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
        mocker: MockerFixture,
    ):
        at, _, mock_screening_service, _, _ = app_test_env_v2

        # Override the search repository mock to return empty results
        # This affects the get_by_review_id call in the page function
        mock_search_repo = at.session_state.search_repo
        mock_search_repo.get_by_review_id.return_value = []

        # Initial run to load the page
        at.run()
        assert not at.exception

        # Find and click the start button
        start_button = None
        for button in at.button:
            if button.label == "Start Abstract Screening":
                start_button = button
                break

        assert start_button is not None, "Start Abstract Screening button not found"

        # Click the button and rerun
        start_button.click()
        at.run()

        assert not at.exception
        mock_screening_service.perform_batch_abstract_screening.assert_not_called()

        # Check that error message is shown in the UI
        error_found = False
        for error in at.error:
            if "No PubMed results found" in str(
                error.value
            ) or "No valid search result IDs to screen" in str(error.value):
                error_found = True
                break
        assert error_found, "Error message about no search results not found in UI"
