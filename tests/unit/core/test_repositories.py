"""Unit tests for repository classes."""

import uuid
from datetime import datetime, timezone
from unittest.mock import MagicMock, create_autospec

import pytest
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlmodel import Session as SQLModelSession

from sr_assistant.core.models import (
    LogRecord,
    ScreenAbstractResult,
    ScreeningDecisionType,
    ScreeningResolution,
    ScreeningStrategyType,
    SearchResult,
    SystematicReview,
)
from sr_assistant.core.repositories import (
    ConstraintViolationError,
    LogRepository,
    RecordNotFoundError,
    RepositoryError,
    ScreenAbstractResultRepository,
    ScreeningResolutionRepository,
    SearchResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.core.schemas import SearchResultFilter
from sr_assistant.core.types import LogLevel, SearchDatabaseSource


@pytest.fixture
def mock_session() -> MagicMock:
    session = create_autospec(SQLModelSession, instance=True)
    mock_exec_result = MagicMock()
    session.exec.return_value = mock_exec_result
    mock_exec_result.first.return_value = None
    mock_exec_result.all.return_value = []
    session.get.return_value = None
    return session


@pytest.fixture
def review_repo() -> SystematicReviewRepository:
    return SystematicReviewRepository()


@pytest.fixture
def search_repo() -> SearchResultRepository:
    return SearchResultRepository()


@pytest.fixture
def sample_search_result() -> SearchResult:
    return SearchResult(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        source_db=SearchDatabaseSource.PUBMED,
        source_id="test_source_id",
        title="Test Title",
        abstract="Test abstract.",
        doi="10.1234/test.doi",
        year="2023",
        authors=["Author A", "Author B"],
    )


def test_create_review(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    review_id = uuid.uuid4()
    review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
        inclusion_criteria=None,
    )
    mock_session.add.return_value = None
    mock_session.flush.return_value = None

    result = review_repo.add(mock_session, review)

    mock_session.add.assert_called_once_with(review)
    mock_session.flush.assert_called_once_with([review])
    assert result == review


def test_get_review(
    review_repo: SystematicReviewRepository, mock_session: MagicMock, mocker: MagicMock
):
    review_id = uuid.uuid4()
    expected_review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="Test question",
        exclusion_criteria="Test exclusion",
        inclusion_criteria=None,
    )
    mock_get_by_id = mocker.patch.object(
        review_repo, "get_by_id", return_value=expected_review
    )

    result = review_repo.get_by_id(mock_session, review_id)

    mock_get_by_id.assert_called_once_with(mock_session, review_id)
    assert result == expected_review


def test_update_review(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    review_id = uuid.uuid4()
    existing_review = SystematicReview(
        id=review_id,
        research_question="old q",
        exclusion_criteria="old exc",
        inclusion_criteria=None,
    )
    review_update_data = SystematicReview(
        id=review_id,
        research_question="new q",
        exclusion_criteria="new exc",
        inclusion_criteria=None,
    )
    mock_session.get.return_value = existing_review
    mock_session.merge.return_value = review_update_data

    result = review_repo.update(mock_session, review_update_data)

    mock_session.get.assert_called_once_with(SystematicReview, review_id)
    mock_session.merge.assert_called_once_with(review_update_data)
    mock_session.flush.assert_called_once_with([review_update_data])
    assert result.research_question == "new q"
    assert result.id == review_id


def test_delete_review(
    review_repo: SystematicReviewRepository, mock_session: MagicMock, mocker: MagicMock
):
    review_id = uuid.uuid4()
    review_to_delete = SystematicReview(
        id=review_id,
        research_question="To delete",
        exclusion_criteria="Del",
        inclusion_criteria=None,
    )
    mock_get_by_id = mocker.patch.object(
        review_repo, "get_by_id", return_value=review_to_delete
    )

    review_repo.delete(mock_session, review_id)

    mock_get_by_id.assert_called_once_with(mock_session, review_id)
    mock_session.delete.assert_called_once_with(review_to_delete)
    mock_session.flush.assert_called_once_with()


def test_get_all_reviews(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    reviews_data = [
        SystematicReview(
            id=uuid.uuid4(),
            research_question="Q1",
            exclusion_criteria="E1",
            inclusion_criteria=None,
        ),
        SystematicReview(
            id=uuid.uuid4(),
            research_question="Q2",
            exclusion_criteria="E2",
            inclusion_criteria=None,
        ),
    ]
    mock_session.exec.return_value.all.return_value = reviews_data

    result = review_repo.get_all(mock_session)

    mock_session.exec.assert_called_once()
    assert result == reviews_data


def test_list_reviews(review_repo: SystematicReviewRepository, mock_session: MagicMock):
    repo = SystematicReviewRepository()
    expected = [
        SystematicReview(
            id=uuid.uuid4(),
            research_question="test q 1",
            exclusion_criteria="exc 1",
            inclusion_criteria=None,
        ),
        SystematicReview(
            id=uuid.uuid4(),
            research_question="test q 2",
            exclusion_criteria="exc 2",
            inclusion_criteria=None,
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.list(mock_session)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_get_review_with_search_results(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
        inclusion_criteria=None,
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_with_search_results(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_get_review_with_all(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
        inclusion_criteria=None,
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_with_all(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_log_repository_get_by_review_id(mock_session: MagicMock) -> None:
    repo = LogRepository()
    review_id = uuid.uuid4()
    expected = [
        LogRecord(
            level=LogLevel.INFO,
            message="Test log 1",
            timestamp=datetime.now(tz=timezone.utc),
        ),
        LogRecord(
            level=LogLevel.WARNING,
            message="Test log 2",
            timestamp=datetime.now(tz=timezone.utc),
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_by_review_id(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_log_repository_add(mock_session: MagicMock) -> None:
    repo = LogRepository()
    log_record = LogRecord(
        level=LogLevel.INFO, message="Test", timestamp=datetime.now(tz=timezone.utc)
    )
    repo.add(mock_session, log_record)
    mock_session.add.assert_called_once_with(log_record)
    mock_session.flush.assert_called_once_with([log_record])


def test_log_repository_add_all(mock_session: MagicMock) -> None:
    repo = LogRepository()
    logs = [
        LogRecord(
            level=LogLevel.INFO,
            message="Test log 1",
            timestamp=datetime.now(tz=timezone.utc),
        ),
        LogRecord(
            level=LogLevel.WARNING,
            message="Test log 2",
            timestamp=datetime.now(tz=timezone.utc),
        ),
    ]
    repo.add_all(mock_session, logs)
    mock_session.add_all.assert_called_once_with(logs)
    mock_session.flush.assert_called_once_with(logs)


def test_search_result_repo_get_by_review_id(mock_session: MagicMock) -> None:
    repo = SearchResultRepository()
    review_id = uuid.uuid4()
    expected_results = [
        SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            source_db=SearchDatabaseSource.PUBMED,
            source_id="1",
            title="Title 1",
        ),
        SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            source_db=SearchDatabaseSource.SCOPUS,
            source_id="2",
            title="Title 2",
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected_results
    results = repo.get_by_review_id(mock_session, review_id)
    assert results == expected_results
    mock_session.exec.assert_called_once()


def test_search_result_repo_get_by_source_details(mock_session: MagicMock) -> None:
    repo = SearchResultRepository()
    review_id = uuid.uuid4()
    source_db = SearchDatabaseSource.PUBMED
    source_id = "12345"
    expected_result = SearchResult(
        id=uuid.uuid4(),
        review_id=review_id,
        source_db=source_db,
        source_id=source_id,
        title="Specific Title",
    )
    mock_session.exec.return_value.first.return_value = expected_result
    result = repo.get_by_source_details(mock_session, review_id, source_db, source_id)
    assert result == expected_result
    mock_session.exec.assert_called_once()


def test_screen_repo_get_by_review_id(mock_session: MagicMock) -> None:
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale 1",
            model_name="test_model",
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale 2",
            model_name="test_model",
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_by_review_id(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repo_get_by_strategy(mock_session: MagicMock) -> None:
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    strategy = ScreeningStrategyType.CONSERVATIVE
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            screening_strategy=strategy,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale 1",
            model_name="test_model",
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_by_strategy(mock_session, review_id, strategy)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repo_get_by_exclusion_category(mock_session: MagicMock) -> None:
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    exclusion_reason = "Non-human subjects"
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
            model_name="test_model",
            exclusion_reason_categories={
                "population_exclusion_reasons": [exclusion_reason]
            },
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_by_exclusion_category(mock_session, review_id, exclusion_reason)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_resolution_repo_get_by_search_result_id(mock_session: MagicMock) -> None:
    repo = ScreeningResolutionRepository()
    search_result_id_val = uuid.uuid4()
    expected = ScreeningResolution(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        search_result_id=search_result_id_val,
        resolver_decision=ScreeningDecisionType.INCLUDE,
        resolver_reasoning="Conflict resolved",
        resolver_confidence_score=1.0,
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_by_search_result_id(mock_session, search_result_id_val)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_systematic_review_repository_error_handling(mock_session: MagicMock) -> None:
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    mock_session.exec.side_effect = SQLAlchemyError("DB error")
    with pytest.raises(RepositoryError, match="get_by_id.*DB error"):
        repo.get_with_search_results(mock_session, review_id)

    mock_session.exec.side_effect = None
    mock_session.exec.return_value.first.return_value = None
    mock_session.exec.side_effect = SQLAlchemyError("DB error")

    with pytest.raises(RepositoryError, match="get_by_id.*DB error"):
        repo.get_with_all(mock_session, review_id)


def test_delete_record_not_found(mock_session: MagicMock) -> None:
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    mock_session.exec.return_value.first.return_value = None
    with pytest.raises(RecordNotFoundError, match="not found for deletion"):
        repo.delete(mock_session, review_id)
    mock_session.delete.assert_not_called()
    mock_session.flush.assert_not_called()


def test_base_repo_add_constraint_error(
    review_repo: SystematicReviewRepository, mock_session: MagicMock
):
    review = SystematicReview(
        id=uuid.uuid4(),
        research_question="Constrain me",
        exclusion_criteria="No",
        inclusion_criteria=None,
    )
    mock_session.flush.side_effect = IntegrityError("test", "params", "orig")

    with pytest.raises(ConstraintViolationError):
        review_repo.add(mock_session, review)
    mock_session.add.assert_called_once_with(review)


def test_search_result_repo_get_by_doi_found(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
    sample_search_result: SearchResult,
) -> None:
    mock_session.exec.return_value.first.return_value = sample_search_result
    assert sample_search_result.doi is not None
    result = search_repo.get_by_doi(mock_session, doi=sample_search_result.doi)
    assert result == sample_search_result


def test_search_result_repo_get_by_title_and_year_found(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
    sample_search_result: SearchResult,
) -> None:
    mock_session.exec.return_value.first.return_value = sample_search_result
    assert sample_search_result.title is not None
    assert sample_search_result.year is not None
    result = search_repo.get_by_title_and_year(
        mock_session, title=sample_search_result.title, year=sample_search_result.year
    )
    assert result == sample_search_result


def test_search_result_repo_advanced_search_no_filters(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
    sample_search_result: SearchResult,
) -> None:
    """Test advanced_search with no filters, returning all (mocked) results."""
    another_result = sample_search_result.model_copy(deep=True)
    another_result.id = uuid.uuid4()
    another_result.title = "Another Test Title"
    expected_results = [sample_search_result, another_result]
    mock_session.exec.return_value.all.return_value = expected_results

    results = search_repo.advanced_search(
        mock_session, search_params=SearchResultFilter(), skip=0, limit=10
    )

    assert results == expected_results
    mock_session.exec.assert_called_once()
    # Check that the select statement was called without a where clause initially, then offset and limit
    # This is a bit fragile as it depends on the exact construction of the query
    args, _ = mock_session.exec.call_args
    query_str = str(args[0]).upper()
    assert "WHERE" not in query_str  # No WHERE clause from filters
    assert "OFFSET" in query_str  # Check for OFFSET keyword
    assert "LIMIT" in query_str  # Check for LIMIT keyword


def test_search_result_repo_advanced_search_with_filters(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
    sample_search_result: SearchResult,
) -> None:
    """Test advanced_search with various filters."""
    mock_session.exec.return_value.all.return_value = [sample_search_result]
    review_id_filter = uuid.uuid4()

    search_filters = SearchResultFilter(
        review_id=review_id_filter,
        source_db=SearchDatabaseSource.SCOPUS,
        title="Specific Title",
        year="2024",
        doi="10.specific/doi",
    )

    results = search_repo.advanced_search(
        mock_session, search_params=search_filters, skip=0, limit=10
    )

    assert results == [sample_search_result]
    mock_session.exec.assert_called_once()
    args, _ = mock_session.exec.call_args
    query_str = str(args[0]).upper()

    assert "search_results.review_id = :REVIEW_ID_1".upper() in query_str
    assert "search_results.source_db = :SOURCE_DB_1".upper() in query_str
    # The actual SQL might use functions like LOWER for ILIKE
    assert "LOWER(SEARCH_RESULTS.TITLE) LIKE LOWER(:TITLE_1)" in query_str
    assert "search_results.year = :YEAR_1".upper() in query_str
    assert "search_results.doi = :DOI_1".upper() in query_str
    assert "OFFSET" in query_str
    assert "LIMIT" in query_str


def test_search_result_repo_advanced_search_pagination(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
) -> None:
    """Test advanced_search pagination (skip and limit)."""
    mock_session.exec.return_value.all.return_value = []  # Don't care about results for this

    search_repo.advanced_search(
        mock_session, search_params=SearchResultFilter(), skip=10, limit=50
    )

    mock_session.exec.assert_called_once()
    args, _ = mock_session.exec.call_args
    query_str = str(args[0]).upper()
    assert "OFFSET" in query_str  # Check for OFFSET keyword
    assert "LIMIT" in query_str  # Check for LIMIT keyword


def test_search_result_repo_advanced_search_error_handling(
    mock_session: MagicMock, search_repo: SearchResultRepository
) -> None:
    """Test error handling in advanced_search."""
    mock_session.exec.side_effect = SQLAlchemyError("DB search error")

    with pytest.raises(RepositoryError, match="DB search error"):
        search_repo.advanced_search(mock_session, search_params=SearchResultFilter())


def test_search_result_repo_count_no_filters(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
) -> None:
    """Test count with no filters."""
    expected_count = 5
    mock_session.exec.return_value.scalar_one.return_value = expected_count

    count = search_repo.count(mock_session, search_params=None)

    assert count == expected_count
    mock_session.exec.assert_called_once()
    args, _ = mock_session.exec.call_args
    query_str = str(args[0]).upper()
    assert "WHERE" not in query_str
    assert "COUNT(*)" in query_str


def test_search_result_repo_count_with_filters(
    mock_session: MagicMock,
    search_repo: SearchResultRepository,
) -> None:
    """Test count with various filters."""
    expected_count = 1
    mock_session.exec.return_value.scalar_one.return_value = expected_count
    review_id_filter = uuid.uuid4()

    search_filters = SearchResultFilter(
        review_id=review_id_filter, title="Filtered Title"
    )

    count = search_repo.count(mock_session, search_params=search_filters)

    assert count == expected_count
    mock_session.exec.assert_called_once()
    args, _ = mock_session.exec.call_args
    query_str = str(args[0]).upper()

    assert "search_results.review_id = :REVIEW_ID_1".upper() in query_str
    assert "LOWER(SEARCH_RESULTS.TITLE) LIKE LOWER(:TITLE_1)" in query_str
    assert "COUNT(*)" in query_str


def test_search_result_repo_count_error_handling(
    mock_session: MagicMock, search_repo: SearchResultRepository
) -> None:
    """Test error handling in count."""
    mock_session.exec.side_effect = SQLAlchemyError("DB count error")

    with pytest.raises(RepositoryError, match="DB count error"):
        search_repo.count(mock_session, search_params=SearchResultFilter())
