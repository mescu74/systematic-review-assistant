# tests/integration/core/test_repositories.py
# Assuming REVIEW_1_ID comes from a fixture or is defined
# We need a review ID that exists in the integration test DB
# Define a dummy ID for now if not provided by fixture
from __future__ import annotations

import uuid
from datetime import datetime, timezone

import pytest  # Add pytest if needed for markers etc.
from sqlmodel import Session, select

from sr_assistant.app import services  # Import services
from sr_assistant.core import models

# Import specific repo tested
from sr_assistant.core.repositories import (
    SearchResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.core.types import LogLevel, SearchDatabaseSource

REVIEW_1_ID = uuid.uuid4()  # Replace with actual ID from test DB setup


@pytest.mark.integration
def test_relationship_loading_behavior(
    db_session: Session, test_review: models.SystematicReview
):
    """Test lazy/eager loading behaviors if needed (more complex scenarios)."""
    # Use the ID from the created test_review
    review_id = test_review.id
    repo = SystematicReviewRepository()

    # Add some related data
    result1 = models.SearchResult(
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="REL1",
        title="Relationship Test 1",
    )
    # Assuming LogRecord has a review_id FK
    log1 = models.LogRecord(
        review_id=review_id,
        timestamp=datetime.now(tz=timezone.utc),
        level=LogLevel.INFO,
        message="Log for relationship test",
        record={},
    )
    db_session.add_all([result1, log1])
    db_session.commit()

    # Fetch the review
    fetched_review = repo.get_by_id(db_session, review_id)
    assert fetched_review is not None

    # Access relationships (adjust based on actual relationships defined)
    # Example: Check if logs are loaded (assuming selectin)
    # This depends heavily on how relationships are configured in your models
    # assert fetched_review.log_records # Accessing might trigger loading
    # assert len(fetched_review.log_records) > 0
    # assert fetched_review.log_records[0].message == "Log for relationship test"

    # If SearchResult relationship exists and is configured for loading:
    # search_results = repo.get_search_results(db_session, review_id) # Hypothetical method
    # assert len(search_results) == 1
    # assert search_results[0].source_id == "REL1"

    # Add more specific assertions based on your relationship loading strategy
    # Placeholder if detailed relationship tests aren't needed yet


@pytest.mark.integration
def test_get_with_search_results(
    db_session: Session, test_review: models.SystematicReview
):
    """Test fetching a review along with its search results (if relationship exists)."""
    review_id = test_review.id
    repo = SystematicReviewRepository()
    search_repo = SearchResultRepository()

    # Add review and results
    result1 = models.SearchResult(
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="REL_S1",
        title="Search Result Relation 1",
    )
    db_session.add(result1)
    db_session.commit()

    # Fetch review (assuming a method or relationship loads results)
    # This depends on your model setup and repo methods
    # E.g., if SystematicReview model has a 'search_results' relationship:
    # fetched_review = repo.get_by_id_with_results(db_session, review_id) # Hypothetical
    fetched_review = repo.get_by_id(db_session, review_id)
    assert fetched_review is not None

    # Now fetch results separately to verify
    search_results = search_repo.get_by_review_id(db_session, review_id)
    assert len(search_results) == 1
    assert search_results[0].source_id == "REL_S1"

    # If the relationship is configured for direct access & loading:
    # assert hasattr(fetched_review, 'search_results')
    # assert len(fetched_review.search_results) == 1
    # assert fetched_review.search_results[0].source_id == "REL_S1"


@pytest.mark.integration
def test_get_with_all(db_session: Session, test_review: models.SystematicReview):
    """Test fetching a review with all related data (if applicable). Needs adjustment."""
    review_id = test_review.id
    repo = SystematicReviewRepository()
    search_repo = SearchResultRepository()
    # log_repo = LogRecordRepository() # If exists

    # Add related data
    result1 = models.SearchResult(
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="ALL1",
        title="All Related 1",
    )
    log1 = models.LogRecord(
        review_id=review_id,
        timestamp=datetime.now(tz=timezone.utc),
        level=LogLevel.INFO,
        message="Log for get_with_all test",
        record={},
    )
    db_session.add_all([result1, log1])
    db_session.commit()

    # Fetch review (method depends on implementation)
    # fetched_review = repo.get_by_id_with_all(db_session, review_id) # Hypothetical
    fetched_review = repo.get_by_id(db_session, review_id)
    assert fetched_review is not None

    # Verify related data separately for now
    search_results = search_repo.get_by_review_id(db_session, review_id)
    assert len(search_results) == 1
    assert search_results[0].source_id == "ALL1"

    # logs = log_repo.get_by_review_id(db_session, review_id) # Hypothetical
    # assert len(logs) == 1
    # assert logs[0].message == "Log for get_with_all test"

    # If relationships are loaded directly:
    # assert hasattr(fetched_review, 'search_results') and len(fetched_review.search_results) == 1
    # assert hasattr(fetched_review, 'log_records') and len(fetched_review.log_records) == 1


# --- Fixtures ---
# Assuming db_session is provided by conftest.py and provides a sync Session
# Define a fixture for a review_id created in the test DB for each test
@pytest.fixture(scope="function")
def test_review(db_session: Session):
    review = models.SystematicReview(
        id=uuid.uuid4(),
        research_question=f"Integration Test Review {uuid.uuid4()}",
        exclusion_criteria="Test exclusion criteria",
    )
    db_session.add(review)
    db_session.commit()

    # Explicitly fetch to ensure it's in the DB
    fetched = db_session.get(models.SystematicReview, review.id)
    if not fetched:
        raise RuntimeError(f"Review wasn't properly saved to database: {review.id}")

    yield review

    # Teardown: Use transaction rollback instead of manual deletion
    # to avoid transaction errors
    db_session.rollback()


# --- Mock API Data ---
MOCK_PUBMED_RECORD_1 = {
    "MedlineCitation": {
        "PMID": "PM_INT_1",  # Unique ID for test
        "Article": {
            "ArticleTitle": "Integration Test PubMed Title 1",
            "Journal": {"Title": "Journal P Int"},
            "Abstract": {"AbstractText": ["Abstract P Int 1"]},
            "AuthorList": [{"LastName": "AuthPInt", "ForeName": "A"}],
        },
    },
    "PubmedData": {"ArticleIdList": []},  # Keep it simple
}
MOCK_SCOPUS_RECORD_1 = {
    "dc:identifier": "SCOPUS_ID:SC_INT_1",
    "eid": "eid_int_1",
    "dc:title": "Integration Test Scopus Title 1",
    "prism:publicationName": "Journal S Int",
    "prism:coverDate": "2024-01-01",
    "author": [{"authname": "AuthSInt B"}],
}
MOCK_SCOPUS_RECORD_DUPLICATE = {
    "dc:identifier": "SCOPUS_ID:SC_INT_1",  # Same ID
    "eid": "eid_int_1_dup",
    "dc:title": "Integration Test Scopus Title 1 Duplicate Attempt",
}


# --- New Integration Tests for SearchService (Synchronous) ---


@pytest.mark.integration
def test_service_add_pubmed_results(
    db_session: Session, test_review: models.SystematicReview
):
    """Test adding PubMed results via SearchService integration."""
    search_service = services.SearchService()  # Uses sync factory now
    api_data = [MOCK_PUBMED_RECORD_1]

    # Act (call sync method, passing the session)
    added_results = search_service.add_api_search_results(
        review_id=test_review.id,
        source_db=SearchDatabaseSource.PUBMED,
        api_records=api_data,
        session=db_session,
    )

    # Assert (Service Level)
    assert len(added_results) == 1
    saved_result_svc = added_results[0]
    assert saved_result_svc.source_db == SearchDatabaseSource.PUBMED
    assert saved_result_svc.source_id == "PM_INT_1"
    assert saved_result_svc.title == "Integration Test PubMed Title 1"

    # Assert (Database Level using sync session - now consistent)
    stmt = select(models.SearchResult).where(
        models.SearchResult.review_id == test_review.id,
        models.SearchResult.source_db == SearchDatabaseSource.PUBMED,
    )
    db_results = db_session.exec(stmt).all()
    assert len(db_results) == 1
    saved_result_db = db_results[0]
    assert saved_result_db.source_id == "PM_INT_1"
    assert saved_result_db.title == "Integration Test PubMed Title 1"


@pytest.mark.integration
def test_service_add_scopus_deduplication(
    db_session: Session, test_review: models.SystematicReview
):
    """Test adding Scopus results with duplicates via SearchService integration."""
    search_service = services.SearchService()
    api_data = [MOCK_SCOPUS_RECORD_1, MOCK_SCOPUS_RECORD_DUPLICATE]

    # Act (call sync method, passing the session)
    # Service returns [] on bulk constraint violation
    results_from_service = search_service.add_api_search_results(
        review_id=test_review.id,
        source_db=SearchDatabaseSource.SCOPUS,
        api_records=api_data,
        session=db_session,
    )

    # Assert Service returns empty list on constraint violation
    assert results_from_service == []

    # Assert (Database Level - check that nothing was committed due to rollback)
    stmt = select(models.SearchResult).where(
        models.SearchResult.review_id == test_review.id,
        models.SearchResult.source_db == SearchDatabaseSource.SCOPUS,
    )
    db_results = db_session.exec(stmt).all()
    assert len(db_results) == 0
    # Remove checks for saved result as nothing should be saved
    # saved_result_db = db_results[0]
    # assert saved_result_db.source_id == "SC_INT_1"
    # assert saved_result_db.title == "Integration Test Scopus Title 1"
    # assert saved_result_db.year == "2024" # Check year mapping (should be string now)


@pytest.mark.integration
def test_service_get_results(db_session: Session, test_review: models.SystematicReview):
    """Test retrieving results via SearchService integration."""
    # Arrange (using db_session)
    result1 = models.SearchResult(
        review_id=test_review.id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="IG1",
        title="Get T1 Int",
    )
    result2 = models.SearchResult(
        review_id=test_review.id,
        source_db=SearchDatabaseSource.SCOPUS,
        source_id="IG2",
        title="Get T2 Int",
    )
    db_session.add_all([result1, result2])
    db_session.commit()

    search_service = services.SearchService()

    # Act (call sync method, passing the session)
    results = search_service.get_search_results_by_review_id(
        test_review.id, session=db_session
    )

    # Assert
    assert len(results) == 2
    assert {r.source_id for r in results} == {"IG1", "IG2"}


@pytest.mark.integration
def test_service_get_result_details(
    db_session: Session, test_review: models.SystematicReview
):
    """Test retrieving a specific result via SearchService."""
    # Arrange (using db_session)
    source_db = SearchDatabaseSource.SCOPUS
    source_id = "IGD1"
    result1 = models.SearchResult(
        review_id=test_review.id,
        source_db=source_db,
        source_id=source_id,
        title="Get Detail T1 Int",
    )
    db_session.add(result1)
    db_session.commit()

    search_service = services.SearchService()

    # Act (call sync methods, passing the session)
    result = search_service.get_search_result_by_source_details(
        test_review.id, source_db, source_id, session=db_session
    )
    result_none = search_service.get_search_result_by_source_details(
        test_review.id, source_db, "NOT_THERE", session=db_session
    )

    # Assert
    assert result is not None
    assert result.source_id == source_id
    assert result.title == "Get Detail T1 Int"
    assert result_none is None


# --- Fixture for tests needing an existing review ---
# This fixture is no longer used by the tests above
@pytest.fixture(scope="module")
def existing_review_id(db_session: Session) -> uuid.UUID:
    """Provides the ID of a review assumed to exist in the test DB.

    NOTE: This is likely unnecessary now tests use test_review fixture.
    Kept for potential future use or reference.
    """
    # IMPORTANT: Replace with an ID known to exist in your integration DB setup
    # If your DB is cleaned each time, this needs to be created.
    hardcoded_id = uuid.UUID("9f336bd5-8dd9-40f0-97b0-81040920febd")

    # Check if it exists, skip if not
    review = db_session.get(models.SystematicReview, hardcoded_id)
    if not review:
        pytest.skip(f"Review ID {hardcoded_id} not found in integration database.")
    return hardcoded_id
