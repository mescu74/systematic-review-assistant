# pyright: ignore [reportPrivateUsage]

from __future__ import annotations

import typing
import uuid
from collections.abc import KeysView
from unittest.mock import MagicMock, patch

import pytest
import streamlit as st
from pytest_mock import MockerFixture

from sr_assistant.app.pages.search import gen_query_cb, search_page
from sr_assistant.core.models import SystematicReview
from sr_assistant.core.schemas import SearchResultRead
from sr_assistant.core.types import SearchDatabaseSource


# Helper to create a mock review
def mock_systematic_review(review_id: uuid.UUID | None = None) -> SystematicReview:
    return SystematicReview(
        id=review_id or uuid.uuid4(),
        research_question="Test question?",
        exclusion_criteria="Test exclusion",
        # Add other mandatory fields if any
    )


# Helper to create mock search results
def mock_search_result_read(count: int = 1) -> list[SearchResultRead]:
    results = []
    for i in range(count):
        results.append(
            SearchResultRead(
                id=uuid.uuid4(),
                review_id=uuid.uuid4(),  # This will be the review_id of the parent review in real scenarios
                source_db=SearchDatabaseSource.PUBMED,
                source_id=f"PMID{i + 1}",
                doi=f"doi_{i + 1}",
                title=f"Test Title {i + 1}",
                abstract=f"Abstract {i + 1} content.",
                journal=f"Test Journal {i + 1}",
                year=str(2020 + i),
                authors=[f"Author {A}" for A in "ABC"],
                keywords=["test", f"keyword{i + 1}"],
                raw_data={},
                source_metadata={},
                created_at=None,  # Explicitly set to None
                updated_at=None,  # Explicitly set to None
            )
        )
    return results


@pytest.fixture(autouse=True)
def clear_st_session_state(mocker: MockerFixture):
    """Clears the *actual* streamlit.session_state before each test.
    This is important if any code under test (or other tests) might modify it directly.
    For tests fully mocking 'st', this fixture primarily ensures a clean global state
    if st is imported elsewhere without being mocked.
    """
    # Use streamlit.session_state for clearing, as this is the global object.
    # The test-specific mock_st.session_state will be a fresh MagicMock attribute.
    global_st_keys = list(st.session_state.keys())
    for key in global_st_keys:
        del st.session_state[key]
    # If config is needed globally, set it on actual st.session_state
    # st.session_state.config = mocker.MagicMock()
    # st.session_state.config.NCBI_EMAIL = "test@example.com"
    # st.session_state.config.NCBI_API_KEY.get_secret_value.return_value = "test_api_key"
    # However, for tests where 'st' is fully mocked, interaction should be via the mock.


@pytest.fixture
def mock_session_state_manager(
    mocker: MockerFixture,
) -> tuple[MagicMock, dict[str, typing.Any]]:
    """Mocks sr_assistant.app.pages.search.st.session_state with a custom MockSessionState object
    that uses a dictionary as a backing store.
    Returns the mock for 'st' module and the backing dictionary for direct assertions if needed.
    """
    _backing_dict = {}

    class MockSessionState:
        def __init__(
            self,
            backing_dict_ref: dict[str, typing.Any],
            mocker_instance: MockerFixture,
        ):
            object.__setattr__(self, "_backing_dict", backing_dict_ref)
            object.__setattr__(self, "_mocker", mocker_instance)

        def __getattr__(self, name: str) -> typing.Any:
            if name in self._backing_dict:
                return self._backing_dict[name]
            new_mock = self._mocker.MagicMock(name=f"session_state_auto_attr_{name}")
            self._backing_dict[name] = new_mock
            return new_mock

        def __setattr__(self, name: str, value: typing.Any) -> None:
            self._backing_dict[name] = value

        def __delattr__(self, name: str) -> None:
            if name in self._backing_dict:
                del self._backing_dict[name]
            else:
                pass

        def __contains__(self, key: str) -> bool:
            return key in self._backing_dict

        def get(self, key: str, default: typing.Any = None) -> typing.Any:
            return self._backing_dict.get(key, default)

        def keys(self) -> KeysView[str]:
            return self._backing_dict.keys()

    patched_st_object = mocker.patch("sr_assistant.app.pages.search.st")
    patched_st_object.session_state = MockSessionState(_backing_dict, mocker)
    return patched_st_object, _backing_dict


def test_search_page_initial_load_and_search(
    mocker: MockerFixture,
    mock_session_state_manager: tuple[MagicMock, dict[str, typing.Any]],
):
    """Test initial page load, query generation, search execution, and result display."""
    mock_st_object, _ = mock_session_state_manager

    MockSearchService = mocker.patch("sr_assistant.app.pages.search.SearchService")
    mock_init_review_repo = mocker.patch(
        "sr_assistant.app.pages.search.init_review_repository"
    )
    mock_init_query_chain = mocker.patch(
        "sr_assistant.app.pages.search.init_query_chain"
    )

    expected_query_value = "initial_query_from_get_query"
    expected_max_results = 20

    mock_get_query = mocker.patch(
        "sr_assistant.app.pages.search.get_query",
        return_value=expected_query_value,
    )

    mock_st_object.columns.return_value = (
        mocker.MagicMock(),
        mocker.MagicMock(),
        mocker.MagicMock(),
    )
    mock_st_object.text_area.return_value = expected_query_value
    mock_st_object.slider.return_value = expected_max_results

    review_id = uuid.uuid4()
    mock_review_obj = mock_systematic_review(review_id)

    mock_repo_instance = mocker.MagicMock()
    mock_repo_instance.get_by_id.return_value = mock_review_obj

    def mock_init_review_repository_side_effect():
        mock_st_object.session_state.repo_review = mock_repo_instance
        return mock_repo_instance

    mock_init_review_repo.side_effect = mock_init_review_repository_side_effect

    mock_search_service_instance = MockSearchService.return_value
    mock_results = mock_search_result_read(2)
    mock_search_service_instance.search_pubmed_and_store_results.return_value = (
        mock_results
    )

    mock_st_object.session_state.review_id = review_id
    mock_st_object.session_state.query_value = None
    mock_st_object.session_state.max_results = expected_max_results

    if "review" in mock_st_object.session_state:
        del mock_st_object.session_state.review

    def button_side_effect(
        label: str, *args: typing.Any, **kwargs: typing.Any
    ) -> bool | typing.Any:
        if label == "Search":
            return True
        return False

    mock_st_object.button.side_effect = button_side_effect

    search_page(review_id)

    mock_init_review_repo.assert_called_once()
    mock_repo_instance.get_by_id.assert_called_with(review_id)
    mock_init_query_chain.assert_called_once()
    mock_get_query.assert_called_once_with(mock_review_obj)

    mock_search_service_instance.search_pubmed_and_store_results.assert_called_once_with(
        review_id=review_id,
        query=expected_query_value,
        max_results=expected_max_results,
    )

    expected_df_data = [
        {
            "Source ID": r.source_id,
            "DOI": r.doi,
            "Title": r.title,
            "Journal": r.journal,
            "Year": r.year,
        }
        for r in mock_results
    ]

    mock_st_object.dataframe.assert_called_once_with(
        expected_df_data, use_container_width=True
    )

    assert mock_st_object.session_state.search_results == mock_results
    mock_st_object.success.assert_any_call(
        f"Search complete. Found {len(mock_results)} articles."
    )


def test_search_page_search_service_error(
    mocker: MockerFixture,
    mock_session_state_manager: tuple[MagicMock, dict[str, typing.Any]],
):
    """Test error handling when SearchService raises an exception."""
    mock_st_object, _ = mock_session_state_manager

    MockSearchService = mocker.patch("sr_assistant.app.pages.search.SearchService")
    mock_init_review_repo = mocker.patch(
        "sr_assistant.app.pages.search.init_review_repository"
    )
    mocker.patch("sr_assistant.app.pages.search.init_query_chain")

    initial_setup_query = "error_query_for_page_load"
    mocker.patch(
        "sr_assistant.app.pages.search.get_query", return_value=initial_setup_query
    )

    mock_st_object.columns.return_value = (
        mocker.MagicMock(),
        mocker.MagicMock(),
        mocker.MagicMock(),
    )

    actual_query_for_service = "actual_error_query_from_ui"
    actual_max_results_for_service = 10
    mock_st_object.text_area.return_value = actual_query_for_service
    mock_st_object.slider.return_value = actual_max_results_for_service

    review_id = uuid.uuid4()
    mock_review_obj = mock_systematic_review(review_id)

    mock_repo_instance = mocker.MagicMock()
    mock_repo_instance.get_by_id.return_value = mock_review_obj

    def mock_init_review_repository_side_effect():
        mock_st_object.session_state.repo_review = mock_repo_instance
        return mock_repo_instance

    mock_init_review_repo.side_effect = mock_init_review_repository_side_effect

    mock_search_service_instance = MockSearchService.return_value
    mock_search_service_instance.search_pubmed_and_store_results.side_effect = (
        Exception("Service Error")
    )

    mock_st_object.session_state.review_id = review_id
    mock_st_object.session_state.query_value = None
    mock_st_object.session_state.max_results = actual_max_results_for_service

    if "review" in mock_st_object.session_state:
        del mock_st_object.session_state.review

    def button_side_effect(
        label: str, *args: typing.Any, **kwargs: typing.Any
    ) -> bool | typing.Any:
        if label == "Search":
            return True
        return False

    mock_st_object.button.side_effect = button_side_effect

    search_page(review_id)

    mock_search_service_instance.search_pubmed_and_store_results.assert_called_once_with(
        review_id=review_id,
        query=actual_query_for_service,
        max_results=actual_max_results_for_service,
    )
    mock_st_object.error.assert_called_with("Service Error")


def test_search_page_generate_query_button(
    mocker: MockerFixture,
    mock_session_state_manager: tuple[MagicMock, dict[str, typing.Any]],
):
    """Test the 'Generate query' button callback."""
    mock_st_object, _ = mock_session_state_manager  # Use _ for unused backing_dict

    mocker.patch("sr_assistant.app.pages.search.SearchService")
    mock_init_review_repo = mocker.patch(
        "sr_assistant.app.pages.search.init_review_repository"
    )
    mocker.patch("sr_assistant.app.pages.search.init_query_chain")
    mock_get_query = mocker.patch("sr_assistant.app.pages.search.get_query")

    mock_st_object.columns.return_value = (
        mocker.MagicMock(),
        mocker.MagicMock(),
        mocker.MagicMock(),
    )
    mock_st_object.text_area.return_value = "initial_text_area_query"
    mock_st_object.slider.return_value = 50

    review_id = uuid.uuid4()
    mock_review_obj = mock_systematic_review(review_id)

    mock_repo_instance = mocker.MagicMock()
    mock_repo_instance.get_by_id.return_value = mock_review_obj

    def mock_init_review_repository_side_effect():
        mock_st_object.session_state.repo_review = mock_repo_instance
        return mock_repo_instance

    mock_init_review_repo.side_effect = mock_init_review_repository_side_effect

    mock_st_object.session_state.review_id = review_id
    mock_st_object.session_state.review = mock_review_obj
    mock_st_object.session_state.query_value = "query_before_gen_click"
    mock_st_object.session_state.search_results = []

    def temp_button_side_effect(
        label: str, *args: typing.Any, on_click: typing.Any = None, **kwargs: typing.Any
    ) -> bool:
        return False

    mock_st_object.button.side_effect = temp_button_side_effect
    # Call search_page to set up initial state, including st.session_state.review
    # If query_value is already set, get_query won't be called here by search_page itself.
    mock_get_query.return_value = mock_st_object.session_state.query_value
    search_page(review_id)

    expected_query_after_gen = "newly_generated_query_by_button"
    mock_get_query.return_value = expected_query_after_gen
    mock_get_query.reset_mock()

    gen_query_cb()

    mock_get_query.assert_called_once_with(mock_review_obj)
    assert mock_st_object.session_state.query_value == expected_query_after_gen


@patch("sr_assistant.app.pages.search.st")
@patch("sr_assistant.app.pages.search.SearchService")
@patch("sr_assistant.app.pages.search.init_review_repository")
@patch("sr_assistant.app.pages.search.init_query_chain")
def test_search_page_clear_results_button(
    mock_init_query_chain: MagicMock,
    mock_init_review_repo: MagicMock,
    MockSearchService: MagicMock,
    mock_st: MagicMock,
):
    """Test the 'Clear search results' button."""
    mock_st.columns.return_value = (MagicMock(), MagicMock(), MagicMock())
    review_id = uuid.uuid4()
    mock_review = mock_systematic_review(review_id)

    mock_repo_instance = MagicMock()
    mock_repo_instance.get_by_id.return_value = mock_review

    def mock_init_side_effect():
        mock_st.session_state.repo_review = mock_repo_instance
        return mock_repo_instance

    mock_init_review_repo.side_effect = mock_init_side_effect

    mock_st.session_state.review_id = review_id
    mock_st.session_state.search_results = mock_search_result_read(3)
    mock_st.session_state.review = mock_review
    mock_st.session_state.query_value = "some_query_for_this_test"

    def button_side_effect(label: str, *args: typing.Any, **kwargs: typing.Any) -> bool:
        if label == "Clear search results":
            return True
        return False

    mock_st.button.side_effect = button_side_effect
    if "review" in mock_st.session_state:
        del mock_st.session_state.review
    search_page(review_id)

    assert mock_st.session_state.search_results == []
    mock_st.success.assert_called_with(
        "Search results cleared from display. Re-search to fetch again."
    )
    mock_st.rerun.assert_called_once()


def test_search_page_display_selected_article_details(
    mocker: MockerFixture, mock_session_state_manager
):
    """Test displaying details of a selected article."""
    mock_st_object, _ = mock_session_state_manager

    mocker.patch("sr_assistant.app.pages.search.SearchService")
    mock_init_review_repo = mocker.patch(
        "sr_assistant.app.pages.search.init_review_repository"
    )
    mocker.patch("sr_assistant.app.pages.search.init_query_chain")
    mocker.patch(
        "sr_assistant.app.pages.search.get_query", return_value="display_detail_query"
    )

    mock_st_object.columns.return_value = (
        mocker.MagicMock(),
        mocker.MagicMock(),
        mocker.MagicMock(),
    )
    mock_st_object.text_area.return_value = "display_detail_query"
    mock_st_object.slider.return_value = 50

    review_id = uuid.uuid4()
    mock_review_obj = mock_systematic_review(review_id)

    mock_repo_instance = mocker.MagicMock()
    mock_repo_instance.get_by_id.return_value = mock_review_obj

    def mock_init_review_repository_side_effect():
        mock_st_object.session_state.repo_review = mock_repo_instance
        return mock_repo_instance

    mock_init_review_repo.side_effect = mock_init_review_repository_side_effect

    mock_results_list = mock_search_result_read(2)
    selected_article_source_id = "SELECTED_PMID_XYZ"
    mock_results_list[0].source_id = selected_article_source_id
    selected_article = mock_results_list[0]

    mock_st_object.session_state.review_id = review_id
    mock_st_object.session_state.search_results = mock_results_list
    mock_st_object.session_state.review = mock_review_obj
    mock_st_object.session_state.query_value = "display_detail_query"

    mock_st_object.selectbox.return_value = selected_article_source_id

    def button_side_effect(label: str, *args: typing.Any, **kwargs: typing.Any) -> bool:
        return False

    mock_st_object.button.side_effect = button_side_effect

    search_page(review_id)

    mock_st_object.subheader.assert_called_with(selected_article.title)
    mock_st_object.text.assert_called_with(
        f"{selected_article.journal} ({selected_article.year})"
    )
    mock_st_object.write.assert_any_call(selected_article.abstract)
    mock_st_object.json.assert_called_once_with(
        selected_article.model_dump(mode="json"), expanded=True
    )
    mock_st_object.selectbox.assert_called_with(
        "Select article by Source ID:", [r.source_id for r in mock_results_list]
    )


# More tests can be added for:
# - Initial query generation if st.session_state.query_value is None
# - Behavior when no review_id is provided or review is not found
# - Different states of st.session_state (e.g. no search_results yet)
