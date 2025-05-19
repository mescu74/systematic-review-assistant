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

    mock_session_factory_patch = mocker.patch(
        "sr_assistant.app.pages.screen_abstracts.session_factory"
    )
    mock_db_session_instance = mocker.MagicMock(spec=Session)
    mock_session_factory_patch.return_value.__enter__.return_value = (
        mock_db_session_instance
    )

    # Populate the MockSessionState instance (which is mock_st_in_page.session_state)
    mock_st_in_page.session_state.review_id = test_review_id
    mock_st_in_page.session_state.review = mock_review_model
    mock_st_in_page.session_state.search_results = [mock_sr1, mock_sr2]
    mock_st_in_page.session_state.repo_review = mock_review_repo_inst
    mock_st_in_page.session_state.search_repo = mock_search_repo_inst
    mock_st_in_page.session_state.repo_screen_abstracts = mock_screen_abstract_repo_inst
    mock_st_in_page.session_state.repo_screen_resolution = (
        mock_screen_resolution_repo_inst
    )

    at = AppTest.from_file(
        "src/sr_assistant/app/pages/screen_abstracts.py", default_timeout=10
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
        at, mock_st, _, _, _ = app_test_env_v2

        at.run()

        assert not at.exception
        mock_st.title.assert_called_once_with("Abstract Screening")

        start_button_call_found = False
        for call in mock_st.button.call_args_list:
            if len(call.args) > 0 and call.args[0] == "Start Abstract Screening":
                start_button_call_found = True
                break
            if call.kwargs.get("key") == "Start Abstract Screening":
                start_button_call_found = True
                break
        assert start_button_call_found, (
            "'Start Abstract Screening' button not defined via st.button"
        )

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
        (
            at,
            mock_st,
            mock_screening_service,
            mock_search_results_list,
            mock_review_model,
        ) = app_test_env_v2

        def button_side_effect(
            label: str, key: str | None = None, **kwargs: t.Any
        ) -> bool:
            if key == "Start Abstract Screening" or label == "Start Abstract Screening":
                return True
            return False

        mock_st.button.side_effect = button_side_effect

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
        (
            at,
            mock_st,
            mock_screening_service,
            mock_search_results_list,
            mock_review_model,
        ) = app_test_env_v2

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
            decision=types.ScreeningDecisionType.UNCERTAIN,
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

        mock_st.button.side_effect = (
            lambda label, key=None, **kwargs: key == "Start Abstract Screening"
        )
        at.run()

        assert not at.exception
        assert len(mock_st.session_state.screen_abstracts_results) == 4
        assert mock_st.session_state.screen_abstracts_screened == 2
        assert mock_st.session_state.screen_abstracts_included == 2
        assert mock_st.session_state.screen_abstracts_excluded == 1
        assert mock_st.session_state.screen_abstracts_uncertain == 1
        assert len(mock_st.session_state.screen_abstracts_conflicts) == 1
        assert (
            mock_st.success.called
        )  # Check if st.success was called by the page script

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
        at, mock_st, mock_screening_service, mock_search_results_list, _ = (
            app_test_env_v2
        )
        sr1 = mock_search_results_list[0]
        mock_error = ScreeningError(
            search_result=sr1, error="Test error", message="Something went wrong"
        )
        service_return_tuples = [ScreenAbstractResultTuple(sr1, mock_error, mock_error)]
        mock_screening_service.perform_batch_abstract_screening.return_value = (
            service_return_tuples
        )

        mock_st.button.side_effect = (
            lambda label, key=None, **kwargs: key == "Start Abstract Screening"
        )
        at.run()

        assert not at.exception
        assert len(mock_st.session_state.screen_abstracts_errors) == 2
        assert mock_st.session_state.screen_abstracts_errors[0] == mock_error
        assert mock_st.session_state.screen_abstracts_screened == 1
        assert mock_st.error.called or mock_st.write.called  # Check calls on mock_st

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
        at, mock_st, mock_screening_service, _, _ = app_test_env_v2
        mock_screening_service.perform_batch_abstract_screening.side_effect = (
            RecordNotFoundError("Review not found by service")
        )

        mock_st.button.side_effect = (
            lambda label, key=None, **kwargs: key == "Start Abstract Screening"
        )
        at.run()

        assert not at.exception
        assert mock_st.error.called
        mock_st.error.assert_any_call(
            "Error during screening: Review not found by service"
        )

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
        at, mock_st, mock_screening_service, _, _ = app_test_env_v2
        mock_screening_service.perform_batch_abstract_screening.side_effect = (
            app_services_module.ServiceError("Generic service failure")
        )

        mock_st.button.side_effect = (
            lambda label, key=None, **kwargs: key == "Start Abstract Screening"
        )
        at.run()

        assert not at.exception
        assert mock_st.error.called
        mock_st.error.assert_any_call(
            "An unexpected error occurred during screening: Generic service failure"
        )

    def test_no_search_results_to_screen(
        self,
        app_test_env_v2: tuple[
            AppTest,
            MagicMock,
            MagicMock,
            list[models.SearchResult],
            models.SystematicReview,
        ],
    ):
        at, mock_st, mock_screening_service, _, _ = app_test_env_v2
        mock_st.session_state.search_results = []

        mock_st.button.side_effect = (
            lambda label, key=None, **kwargs: key == "Start Abstract Screening"
        )
        at.run()

        assert not at.exception
        mock_screening_service.perform_batch_abstract_screening.assert_not_called()
        assert mock_st.error.called

        error_calls = [
            str(call.args[0]) for call in mock_st.error.call_args_list if call.args
        ]
        assert any(
            "No PubMed results found" in err_msg for err_msg in error_calls
        ) or any(
            "No valid search result IDs to screen" in err_msg for err_msg in error_calls
        )
