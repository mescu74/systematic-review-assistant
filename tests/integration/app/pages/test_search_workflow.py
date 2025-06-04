"""Integration tests for the search page workflow.

These tests verify the interaction between the search page UI and the Search service.
They run against a real database (test instance) and require @pytest.mark.integration.
"""

from __future__ import annotations

import uuid

import pytest
from pytest_mock import MockerFixture
from sqlalchemy.orm import Session
from sqlmodel import select
from streamlit.testing.v1 import AppTest

from sr_assistant.core import models
from sr_assistant.core.types import CriteriaFramework


@pytest.fixture
def mockable_streamlit_page_link(mocker: MockerFixture) -> None:
    """Mocks st.page_link to prevent KeyError: 'url_pathname' in test environment."""
    return mocker.patch("streamlit.page_link", return_value=None)


@pytest.fixture
def test_review(db_session: Session) -> models.SystematicReview:
    """Create and return a test SystematicReview instance in the database."""
    review = models.SystematicReview(
        id=uuid.uuid4(),
        research_question="Test Research Question",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={"population": "test pop"},
        background="Test Background",
        exclusion_criteria="Test Exclusion",
    )

    db_session.add(review)
    db_session.commit()
    db_session.refresh(review)

    return review


@pytest.mark.integration
def test_successful_search_and_display(
    mocker: MockerFixture,
    test_review: models.SystematicReview,
    db_session: Session,
    mockable_streamlit_page_link: None,
):
    """Test Case 1: Successful PubMed Search and Display (AC1, AC2, AC4)."""
    # Mock Entrez to avoid actual PubMed API calls
    mock_entrez = mocker.patch("sr_assistant.app.services.Entrez")
    mock_entrez.email = "test@example.com"
    mock_entrez.api_key = "fake_api_key"

    # Mock environment variables for NCBI credentials
    mocker.patch(
        "sr_assistant.app.services.os.getenv",
        side_effect=lambda key, default=None: {
            "NCBI_EMAIL": "test@example.com",
            "NCBI_API_KEY": "fake_api_key",
        }.get(key, default),
    )

    # Mock the query generation
    mock_get_query = mocker.patch("sr_assistant.app.pages.search.get_query")
    generated_query = "test query for PubMed"
    mock_get_query.return_value = generated_query

    # Mock Entrez.read for esearch results (PMIDs)
    mock_entrez.read.side_effect = [
        {"IdList": ["12345", "67890"]},  # esearch results - PMIDs
        [  # efetch results - PubMed records
            {
                "MedlineCitation": {
                    "PMID": "12345",
                    "Article": {
                        "ArticleTitle": "Test Article 1",
                        "Abstract": {"AbstractText": ["This is a test abstract."]},
                        "Journal": {
                            "Title": "Test Journal",
                            "JournalIssue": {"PubDate": {"Year": "2023"}},
                        },
                        "AuthorList": [{"LastName": "Smith", "ForeName": "John"}],
                    },
                },
                "PubmedData": {
                    "ArticleIdList": [{"IdType": "doi", "$": "10.1234/test"}],
                    "PublicationStatus": "epublish",
                },
            },
            {
                "MedlineCitation": {
                    "PMID": "67890",
                    "Article": {
                        "ArticleTitle": "Test Article 2",
                        "Abstract": {"AbstractText": ["Another test abstract."]},
                        "Journal": {
                            "Title": "Journal 2",
                            "JournalIssue": {"PubDate": {"Year": "2024"}},
                        },
                        "AuthorList": [{"LastName": "Jones", "ForeName": "Alice"}],
                    },
                },
                "PubmedData": {"ArticleIdList": [], "PublicationStatus": "epublish"},
            },
        ],
    ]

    # Set up AppTest
    at = AppTest.from_file("src/sr_assistant/app/pages/search.py", default_timeout=30)

    # Initialize session state with the real review instance
    at.session_state.review = test_review
    at.session_state.review_id = test_review.id
    at.session_state.query_value = generated_query

    # Initial run
    at.run()

    # Verify initial state
    assert at.session_state.query_value == generated_query

    # Click search button
    at.button(key="search_button").click().run()

    # Verify no exceptions
    assert not at.exception, f"AppTest raised an exception: {at.exception}"

    # Verify Entrez was called correctly
    mock_entrez.esearch.assert_called_once()
    mock_entrez.efetch.assert_called_once()

    # Verify results are displayed in a dataframe
    assert len(at.dataframe) == 1

    # Verify success message
    assert len(at.success) > 0
    assert "Found" in at.success[0].value

    # Verify records were stored in the database
    statement = select(models.SearchResult).where(
        models.SearchResult.review_id == test_review.id
    )
    records = db_session.execute(statement).scalars().all()
    assert len(records) == 2
    assert records[0].source_id == "12345"
    assert records[1].source_id == "67890"


@pytest.mark.integration
def test_duplicate_search_ui_behavior(
    mocker: MockerFixture,
    test_review: models.SystematicReview,
    db_session: Session,
    mockable_streamlit_page_link: None,
):
    """Test Case 2: UI behavior when a search yielding duplicates is performed (AC3, AC5)."""
    # Mock Entrez to avoid actual PubMed API calls
    mock_entrez = mocker.patch("sr_assistant.app.services.Entrez")
    mock_entrez.email = "test@example.com"
    mock_entrez.api_key = "fake_api_key"

    # Mock environment variables for NCBI credentials
    mocker.patch(
        "sr_assistant.app.services.os.getenv",
        side_effect=lambda key, default=None: {
            "NCBI_EMAIL": "test@example.com",
            "NCBI_API_KEY": "fake_api_key",
        }.get(key, default),
    )

    # Mock the query generation
    mock_get_query = mocker.patch("sr_assistant.app.pages.search.get_query")
    test_query = "duplicate test query"
    mock_get_query.return_value = test_query

    # First search will return one result
    first_search_response = {"IdList": ["12345"]}
    second_search_response = {"IdList": ["12345"]}  # Same PMID

    test_record = {
        "MedlineCitation": {
            "PMID": "12345",
            "Article": {
                "ArticleTitle": "Duplicate Article",
                "Abstract": {"AbstractText": ["Abstract text"]},
                "Journal": {
                    "Title": "Journal C",
                    "JournalIssue": {"PubDate": {"Year": "2023"}},
                },
                "AuthorList": [{"LastName": "Author", "ForeName": "C"}],
            },
        },
        "PubmedData": {"ArticleIdList": [], "PublicationStatus": "epublish"},
    }

    # Configure mock for esearch and efetch calls
    mock_entrez.read.side_effect = [
        first_search_response,  # First esearch
        [test_record],  # First efetch
        second_search_response,  # Second esearch
        [test_record],  # Second efetch (should be filtered by repo)
    ]

    # Set up AppTest
    at = AppTest.from_file("src/sr_assistant/app/pages/search.py", default_timeout=30)
    at.session_state.review = test_review
    at.session_state.review_id = test_review.id
    at.session_state.query_value = test_query
    at.run()

    # First search
    at.button(key="search_button").click().run()

    # Verify no exceptions for first search
    assert not at.exception, f"AppTest raised an exception (1st search): {at.exception}"

    # Perform second search with same query
    at.button(key="search_button").click().run()

    # Verify no exceptions for second search
    assert not at.exception, f"AppTest raised an exception (2nd search): {at.exception}"

    # Verify Entrez was called twice for esearch
    assert mock_entrez.esearch.call_count == 2

    # For duplicates, the UI shows message in a status component or warning
    assert len(at.status) > 0 or len(at.warning) > 0, (
        "Expected status or warning message about search results"
    )

    # Verify duplicate handling - the record should only be in the database once
    statement = select(models.SearchResult).where(
        models.SearchResult.review_id == test_review.id,
        models.SearchResult.source_id == "12345",
    )
    records = db_session.execute(statement).scalars().all()
    assert len(records) == 1


@pytest.mark.integration
def test_ui_error_handling_on_search_failure(
    mocker: MockerFixture,
    test_review: models.SystematicReview,
    db_session: Session,
    mockable_streamlit_page_link: None,
):
    """Test Case 3: UI error handling when SearchService raises an exception."""
    # Mock query generation
    mock_get_query = mocker.patch("sr_assistant.app.pages.search.get_query")
    test_query = "error case query"
    mock_get_query.return_value = test_query

    # Force an error in the service by making Entrez raise an exception
    mock_entrez = mocker.patch("sr_assistant.app.services.Entrez")
    mock_entrez.esearch.side_effect = Exception("Simulated PubMed API failure")

    # Set up AppTest
    at = AppTest.from_file("src/sr_assistant/app/pages/search.py", default_timeout=30)
    at.session_state.review = test_review
    at.session_state.review_id = test_review.id
    at.session_state.query_value = test_query
    at.run()

    # Try search (will error)
    at.button(key="search_button").click().run()

    # Verify AppTest itself didn't crash
    assert not at.exception, f"AppTest itself raised an exception: {at.exception}"

    # There should be an error message visible in the UI
    assert len(at.error) > 0
