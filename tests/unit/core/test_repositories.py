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


def test_create_review(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    review_to_add = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
    )
    repo.add(mock_session, review_to_add)
    mock_session.add.assert_called_once_with(review_to_add)
    mock_session.flush.assert_called_once_with([review_to_add])


def test_get_review(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_by_id(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_update_review(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    existing_review = SystematicReview(
        id=review_id, research_question="old q", exclusion_criteria="old exc"
    )
    review_update_data = SystematicReview(
        id=review_id, research_question="new q", exclusion_criteria="new exc"
    )
    mock_session.get.return_value = existing_review
    mock_session.merge.return_value = review_update_data
    result = repo.update(mock_session, review_update_data)
    mock_session.get.assert_called_once_with(SystematicReview, review_id)
    mock_session.merge.assert_called_once_with(review_update_data)
    mock_session.flush.assert_called_once_with([review_update_data])
    assert result.research_question == "new q"


def test_delete_review(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    review_to_delete = SystematicReview(
        id=review_id, research_question="test q", exclusion_criteria="test exc"
    )
    mock_session.exec.return_value.first.return_value = review_to_delete
    repo.delete(mock_session, review_id)
    mock_session.delete.assert_called_once_with(review_to_delete)
    mock_session.flush.assert_called_once()


def test_get_all_reviews(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    expected = [
        SystematicReview(
            id=uuid.uuid4(), research_question="test q 1", exclusion_criteria="exc 1"
        ),
        SystematicReview(
            id=uuid.uuid4(), research_question="test q 2", exclusion_criteria="exc 2"
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_all(mock_session)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_list_reviews(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    expected = [
        SystematicReview(
            id=uuid.uuid4(), research_question="test q 1", exclusion_criteria="exc 1"
        ),
        SystematicReview(
            id=uuid.uuid4(), research_question="test q 2", exclusion_criteria="exc 2"
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.list(mock_session)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_get_review_with_search_results(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
    )
    search_results_data = [
        SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            source_db=SearchDatabaseSource.PUBMED,
            source_id="12345",
            title="Test Article 1",
        ),
    ]
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_with_search_results(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_get_review_with_all(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        exclusion_criteria="test exc",
    )
    search_result_1 = SearchResult(
        id=uuid.uuid4(),
        review_id=review_id,
        source_db=SearchDatabaseSource.PUBMED,
        source_id="12345",
        title="Test Article 1",
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_with_all(mock_session, review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_log_repository_get_by_review_id(mock_session: MagicMock):
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


def test_log_repository_add(mock_session: MagicMock):
    repo = LogRepository()
    log_record = LogRecord(
        level=LogLevel.INFO, message="Test", timestamp=datetime.now(tz=timezone.utc)
    )
    repo.add(mock_session, log_record)
    mock_session.add.assert_called_once_with(log_record)
    mock_session.flush.assert_called_once_with([log_record])


def test_log_repository_add_all(mock_session: MagicMock):
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


def test_search_result_repo_get_by_review_id(mock_session: MagicMock):
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


def test_search_result_repo_get_by_source_details(mock_session: MagicMock):
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


def test_screen_repo_get_by_review_id(mock_session: MagicMock):
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


def test_screen_repo_get_by_strategy(mock_session: MagicMock):
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


def test_screen_repo_get_by_exclusion_category(mock_session: MagicMock):
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


def test_resolution_repo_get_by_search_result_id(mock_session: MagicMock):
    repo = ScreeningResolutionRepository()
    search_result_id = uuid.uuid4()
    expected = ScreeningResolution(
        id=uuid.uuid4(),
        review_id=uuid.uuid4(),
        resolver_decision=ScreeningDecisionType.INCLUDE,
        resolver_reasoning="Conflict resolved",
        resolver_confidence_score=1.0,
    )
    mock_session.exec.return_value.first.return_value = expected
    result = repo.get_by_search_result_id(mock_session, search_result_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_systematic_review_repository_error_handling(mock_session: MagicMock):
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


def test_delete_record_not_found(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    mock_session.exec.return_value.first.return_value = None
    with pytest.raises(RecordNotFoundError, match="not found for deletion"):
        repo.delete(mock_session, review_id)
    mock_session.delete.assert_not_called()
    mock_session.flush.assert_not_called()


def test_base_repo_add_constraint_error(mock_session: MagicMock):
    repo = SystematicReviewRepository()
    review_to_add = SystematicReview(
        id=uuid.uuid4(), research_question="test q", exclusion_criteria="test exc"
    )
    mock_session.flush.side_effect = IntegrityError(
        "duplicate key", params=None, orig=None
    )
    with pytest.raises(ConstraintViolationError, match="duplicate key"):
        repo.add(mock_session, review_to_add)
    mock_session.add.assert_called_once_with(review_to_add)
