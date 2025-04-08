"""Unit tests for repository classes."""

import uuid
from datetime import datetime, timezone
from unittest.mock import MagicMock, create_autospec

import pytest
from sqlalchemy.exc import SQLAlchemyError
from sqlmodel import Session as SQLModelSession

from sr_assistant.core.models import (
    LogRecord,
    ScreenAbstractResult,
    ScreeningDecisionType,
    ScreeningStrategyType,
    SystematicReview,
)
from sr_assistant.core.repositories import (
    LogRepository,
    RecordNotFoundError,
    RepositoryError,
    ScreenAbstractResultRepository,
    SystematicReviewRepository,
)


@pytest.fixture
def mock_session() -> MagicMock:
    session = create_autospec(SQLModelSession, instance=True)
    mock_exec_result = MagicMock()
    session.exec.return_value = mock_exec_result
    mock_exec_result.first.return_value = None
    mock_exec_result.all.return_value = []
    session.get.return_value = None
    return session


def test_create_review(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    repo.add(
        SystematicReview(
            id=review_id,
            background="test bg",
            research_question="test q",
            inclusion_criteria="test inc",
            exclusion_criteria="test exc",
        )
    )

    mock_session.add.assert_called_once()
    mock_session.commit.assert_called_once()


def test_get_review(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_by_id(review_id)
    assert result == expected


def test_update_review(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    mock_session.get.return_value = review
    mock_session.merge.return_value = review

    result = repo.update(review)
    assert result == review

    mock_session.get.assert_called_once()
    mock_session.merge.assert_called_once_with(review)


def test_delete_review(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    mock_session.exec.return_value.first.return_value = review

    repo.delete(review_id)

    mock_session.delete.assert_called_once_with(review)


def test_get_all_reviews(mock_session):
    repo = SystematicReviewRepository()

    expected = [
        SystematicReview(
            id=uuid.uuid4(),
            background="test bg 1",
            research_question="test q 1",
            inclusion_criteria="test inc 1",
            exclusion_criteria="test exc 1",
        ),
        SystematicReview(
            id=uuid.uuid4(),
            background="test bg 2",
            research_question="test q 2",
            inclusion_criteria="test inc 2",
            exclusion_criteria="test exc 2",
        ),
    ]

    mock_session.exec.return_value = expected

    result = repo.get_all()
    assert result == expected


def test_list_reviews(mock_session):
    repo = SystematicReviewRepository()

    expected = [
        SystematicReview(
            id=uuid.uuid4(),
            background="test bg 1",
            research_question="test q 1",
            inclusion_criteria="test inc 1",
            exclusion_criteria="test exc 1",
        ),
        SystematicReview(
            id=uuid.uuid4(),
            background="test bg 2",
            research_question="test q 2",
            inclusion_criteria="test inc 2",
            exclusion_criteria="test exc 2",
        ),
    ]

    mock_session.exec.return_value.all.return_value = expected

    result = repo.list()
    assert result == expected

    mock_session.exec.assert_called_once()


def test_get_review_with_pubmed_results(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()

    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )
    expected.pubmed_results = [
        PubMedResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_id="12345",
            title="Test Article 1",
            abstract="Test abstract 1",
        ),
        PubMedResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_id="67890",
            title="Test Article 2",
            abstract="Test abstract 2",
        ),
    ]

    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_with_pubmed_results(review_id)
    assert result == expected
    assert len(result.pubmed_results) == 2

    mock_session.exec.assert_called_once()


def test_get_review_with_all(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()

    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    pubmed_result_1 = PubMedResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_id="12345",
        title="Test Article 1",
        abstract="Test abstract 1",
    )
    pubmed_result_2 = PubMedResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_id="67890",
        title="Test Article 2",
        abstract="Test abstract 2",
    )
    expected.pubmed_results = [pubmed_result_1, pubmed_result_2]

    pubmed_result_1.conservative_result = ScreenAbstractResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_result_id=pubmed_result_1.id,
        strategy=ScreeningStrategyType.CONSERVATIVE,
        decision=ScreeningDecisionType.INCLUDE,
        confidence_score=0.9,
        rationale="Test rationale 1",
    )
    pubmed_result_2.comprehensive_result = ScreenAbstractResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_result_id=pubmed_result_2.id,
        strategy=ScreeningStrategyType.COMPREHENSIVE,
        decision=ScreeningDecisionType.EXCLUDE,
        confidence_score=0.8,
        rationale="Test rationale 2",
    )

    expected.log_records = [
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="INFO",
            message="Test log 1",
            timestamp=datetime.now(tz=timezone.utc),
        ),
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="WARNING",
            message="Test log 2",
            timestamp=datetime.now(tz=timezone.utc),
        ),
    ]

    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_with_all(review_id)
    assert result == expected
    assert len(result.pubmed_results) == 2
    assert len(result.log_records) == 2
    assert result.pubmed_results[0].conservative_result is not None
    assert result.pubmed_results[1].comprehensive_result is not None

    mock_session.exec.assert_called_once()


def test_log_repository_get_by_review_id(mock_session):
    repo = LogRepository()
    review_id = uuid.uuid4()

    expected = [
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="INFO",
            message="Test log 1",
            timestamp=datetime.now(tz=timezone.utc),
        ),
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="WARNING",
            message="Test log 2",
            timestamp=datetime.now(tz=timezone.utc),
        ),
    ]

    mock_session.exec.return_value = expected

    result = repo.get_by_review_id(review_id)
    assert result == expected

    mock_session.exec.assert_called_once()


def test_log_repository_add(mock_session):
    repo = LogRepository()
    review_id = uuid.uuid4()

    log = LogRecord(
        id=uuid.uuid4(),
        review_id=review_id,
        level="INFO",
        message="Test log",
        timestamp=datetime.now(tz=timezone.utc),
    )

    mock_session.add.return_value = log
    mock_session.commit.return_value = log
    mock_session.refresh.return_value = log

    result = repo.add(log)

    mock_session.add.assert_called_once_with(log)
    mock_session.commit.assert_called_once()
    mock_session.refresh.assert_called_once_with(log)
    assert result == log


def test_log_repository_add_all(mock_session):
    repo = LogRepository()
    review_id = uuid.uuid4()

    logs = [
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="INFO",
            message="Test log 1",
            timestamp=datetime.now(tz=timezone.utc),
        ),
        LogRecord(
            id=uuid.uuid4(),
            review_id=review_id,
            level="WARNING",
            message="Test log 2",
            timestamp=datetime.now(tz=timezone.utc),
        ),
    ]

    mock_session.add_all.return_value = logs
    mock_session.commit.return_value = logs

    result = repo.add_all(logs)

    mock_session.add_all.assert_called_once_with(logs)
    mock_session.commit.assert_called_once()
    assert result == logs
    for log in logs:
        mock_session.refresh.assert_any_call(log)


def test_screen_repo_get_by_review_id(mock_session):
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    search_result_id_1 = uuid.uuid4()
    search_result_id_2 = uuid.uuid4()

    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=search_result_id_1,
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale 1",
            model_name="test_model",
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=search_result_id_2,
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


def test_screen_repo_get_by_search_result(mock_session):
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    search_result_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=search_result_id,
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
            model_name="test_model",
        )
    ]
    mock_session.exec.return_value.all.return_value = expected
    results = repo.get_by_search_result(mock_session, review_id, search_result_id)
    assert results == expected
    mock_session.exec.assert_called_once()


def test_screen_repo_get_by_strategy(mock_session):
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    strategy = ScreeningStrategyType.CONSERVATIVE
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=uuid.uuid4(),
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


def test_screen_repo_get_by_exclusion_category(mock_session):
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    exclusion_reason = "Non-human subjects"
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=uuid.uuid4(),
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


def test_screen_repo_get_by_response_metadata(mock_session):
    repo = ScreenAbstractResultRepository()
    review_id = uuid.uuid4()
    metadata = {"model": "gpt-4", "temperature": 0.0}
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            search_result_id=uuid.uuid4(),
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
            model_name="test_model",
            response_metadata=metadata,
        ),
    ]
    mock_session.exec.return_value.all.return_value = expected
    result = repo.get_by_response_metadata(mock_session, review_id, metadata)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_systematic_review_repository_error_handling(mock_session):
    repo = SystematicReviewRepository()
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Database error in get_with_pubmed_results for {repo.model_cls.__name__}: Database error",
    ):
        repo.get_with_pubmed_results(uuid.uuid4())

    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Database error in get_with_all for {repo.model_cls.__name__}: Database error",
    ):
        repo.get_with_all(uuid.uuid4())


def test_delete_record_not_found(mock_session):
    repo = SystematicReviewRepository()
    review_id = uuid.uuid4()
    mock_session.exec.return_value.first.return_value = None

    with pytest.raises(RecordNotFoundError, match="not found for deletion"):
        repo.delete(mock_session, review_id)

    mock_session.delete.assert_not_called()
    mock_session.flush.assert_not_called()
