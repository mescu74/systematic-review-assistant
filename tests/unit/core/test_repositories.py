"""Unit tests for repository classes."""

import uuid
from datetime import datetime, timezone
from unittest.mock import MagicMock, create_autospec

import pytest
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlmodel import Session

from sr_assistant.core.models import (
    LogRecord,
    PubMedResult,
    ScreenAbstractResult,
    ScreeningDecisionType,
    ScreeningStrategyType,
    SystematicReview,
)
from sr_assistant.core.repositories import (
    ConstraintViolationError,
    LogRepository,
    PubMedResultRepository,
    RepositoryError,
    ScreenAbstractResultRepository,
    SystematicReviewRepository,
)


@pytest.fixture
def mock_session():
    session = create_autospec(Session, instance=True)
    session.begin.return_value = session
    session.__enter__.return_value = session
    session.__exit__.return_value = None
    return session


@pytest.fixture
def mock_session_factory(mock_session):
    factory = MagicMock()
    factory.begin.return_value = mock_session
    factory.return_value = mock_session
    return factory


@pytest.fixture
def mock_pubmed_result_repository(mock_session_factory):
    return PubMedResultRepository(mock_session_factory)


@pytest.fixture
def mock_screen_result_repository(mock_session_factory):
    return ScreenAbstractResultRepository(mock_session_factory)


def test_create_review(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
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

    mock_session = mock_session_factory()
    mock_session.add.assert_called_once()
    mock_session.commit.assert_called_once()


def test_get_review(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    review_id = uuid.uuid4()
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_by_id(review_id)
    assert result == expected


def test_update_review(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    review_id = uuid.uuid4()
    review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    mock_session = mock_session_factory.begin.return_value
    mock_session.get.return_value = review  # Simulate existing record
    mock_session.merge.return_value = review

    result = repo.update(review)
    assert result == review

    mock_session.get.assert_called_once()
    mock_session.merge.assert_called_once_with(review)


def test_delete_review(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    review_id = uuid.uuid4()
    review = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    # Mock get_by_id to return our review
    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value.first.return_value = review

    repo.delete(review_id)

    # Check that delete was called
    mock_session.delete.assert_called_once_with(review)


def test_get_all_reviews(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)

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

    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = repo.get_all()
    assert result == expected


def test_list_reviews(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)

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

    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value.all.return_value = expected

    result = repo.list()
    assert result == expected

    # Verify that exec was called with a select statement
    mock_session.exec.assert_called_once()


def test_get_review_with_pubmed_results(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    review_id = uuid.uuid4()

    # Create a review with PubMed results
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

    # Mock the session and query
    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_with_pubmed_results(review_id)
    assert result == expected
    assert len(result.pubmed_results) == 2

    # Verify the query was constructed correctly
    mock_session.exec.assert_called_once()


def test_get_review_with_all(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    review_id = uuid.uuid4()

    # Create a review with all relationships
    expected = SystematicReview(
        id=review_id,
        background="test bg",
        research_question="test q",
        inclusion_criteria="test inc",
        exclusion_criteria="test exc",
    )

    # Add PubMed results
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

    # Add screening results
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

    # Add log records
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

    # Mock the session and query
    mock_session = mock_session_factory.begin.return_value
    mock_session.exec.return_value.first.return_value = expected

    result = repo.get_with_all(review_id)
    assert result == expected
    assert len(result.pubmed_results) == 2
    assert len(result.log_records) == 2
    assert result.pubmed_results[0].conservative_result is not None
    assert result.pubmed_results[1].comprehensive_result is not None

    # Verify the query was constructed correctly
    mock_session.exec.assert_called_once()


def test_log_repository_get_by_review_id(mock_session_factory):
    repo = LogRepository(mock_session_factory)
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

    mock_session = mock_session_factory.begin.return_value
    # Mock the query chain
    mock_session.exec.return_value = expected

    result = repo.get_by_review_id(review_id)
    assert result == expected

    # Verify the query was constructed correctly
    mock_session.exec.assert_called_once()


def test_log_repository_add(mock_session_factory):
    repo = LogRepository(mock_session_factory)
    review_id = uuid.uuid4()

    log = LogRecord(
        id=uuid.uuid4(),
        review_id=review_id,
        level="INFO",
        message="Test log",
        timestamp=datetime.now(tz=timezone.utc),
    )

    mock_session = mock_session_factory.return_value
    result = repo.add(log)

    # Verify add was called
    mock_session.add.assert_called_once_with(log)
    mock_session.commit.assert_called_once()
    mock_session.refresh.assert_called_once_with(log)
    assert result == log


def test_log_repository_add_all(mock_session_factory):
    repo = LogRepository(mock_session_factory)
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

    mock_session = mock_session_factory.return_value
    result = repo.add_all(logs)

    # Verify add_all was called
    mock_session.add_all.assert_called_once_with(logs)
    mock_session.commit.assert_called_once()
    assert result == logs
    for log in logs:
        mock_session.refresh.assert_any_call(log)


def test_pubmed_repository_get_by_review_id(mock_pubmed_result_repository):
    review_id = uuid.uuid4()
    expected = [
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

    # Mock the session and query
    mock_session = mock_pubmed_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_pubmed_result_repository.get_by_review_id(review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_pubmed_repository_get_with_screening_results(mock_pubmed_result_repository):
    pubmed_id = uuid.uuid4()
    review_id = uuid.uuid4()
    expected = PubMedResult(
        id=pubmed_id,
        review_id=review_id,
        pubmed_id="12345",
        title="Test Article",
        abstract="Test abstract",
    )
    expected.conservative_result = ScreenAbstractResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_result_id=pubmed_id,
        strategy=ScreeningStrategyType.CONSERVATIVE,
        decision=ScreeningDecisionType.INCLUDE,
        confidence_score=0.9,
        rationale="Test rationale",
    )
    expected.comprehensive_result = ScreenAbstractResult(
        id=uuid.uuid4(),
        review_id=review_id,
        pubmed_result_id=pubmed_id,
        strategy=ScreeningStrategyType.COMPREHENSIVE,
        decision=ScreeningDecisionType.EXCLUDE,
        confidence_score=0.8,
        rationale="Test rationale",
    )

    # Mock the session and query
    mock_session = mock_pubmed_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value.first.return_value = expected

    result = mock_pubmed_result_repository.get_with_screening_results(pubmed_id)
    assert result == expected
    assert result.conservative_result is not None
    assert result.comprehensive_result is not None
    mock_session.exec.assert_called_once()


def test_pubmed_repository_get_by_screening_strategy(mock_pubmed_result_repository):
    review_id = uuid.uuid4()
    expected = [
        PubMedResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_id="12345",
            title="Test Article 1",
            abstract="Test abstract 1",
            conservative_result=ScreenAbstractResult(
                id=uuid.uuid4(),
                review_id=review_id,
                pubmed_result_id=uuid.uuid4(),
                strategy=ScreeningStrategyType.CONSERVATIVE,
                decision=ScreeningDecisionType.INCLUDE,
                confidence_score=0.9,
                rationale="Test rationale",
            ),
        ),
        PubMedResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_id="67890",
            title="Test Article 2",
            abstract="Test abstract 2",
            conservative_result=ScreenAbstractResult(
                id=uuid.uuid4(),
                review_id=review_id,
                pubmed_result_id=uuid.uuid4(),
                strategy=ScreeningStrategyType.CONSERVATIVE,
                decision=ScreeningDecisionType.EXCLUDE,
                confidence_score=0.8,
                rationale="Test rationale",
            ),
        ),
    ]

    # Mock the session and query
    mock_session = mock_pubmed_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_pubmed_result_repository.get_by_screening_strategy(
        review_id, ScreeningStrategyType.CONSERVATIVE
    )
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repository_get_by_review_id(mock_screen_result_repository):
    review_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale 1",
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.COMPREHENSIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale 2",
        ),
    ]

    # Mock the session and query
    mock_session = mock_screen_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_screen_result_repository.get_by_review_id(review_id)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repository_get_by_pubmed_result(mock_screen_result_repository):
    review_id = uuid.uuid4()
    pubmed_result_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=pubmed_result_id,
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=pubmed_result_id,
            strategy=ScreeningStrategyType.COMPREHENSIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale",
        ),
    ]

    # Mock the session and query
    mock_session = mock_screen_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_screen_result_repository.get_by_pubmed_result(
        review_id, pubmed_result_id
    )
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repository_get_by_strategy(mock_screen_result_repository):
    review_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale 1",
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale 2",
        ),
    ]

    # Mock the session and query
    mock_session = mock_screen_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_screen_result_repository.get_by_strategy(
        review_id, ScreeningStrategyType.CONSERVATIVE
    )
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repository_get_by_exclusion_category(mock_screen_result_repository):
    review_id = uuid.uuid4()
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
            population_exclusion_reasons=["Non-human subjects"],
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.COMPREHENSIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale",
            population_exclusion_reasons=["Non-human subjects", "Wrong age group"],
        ),
    ]

    # Mock the session and query
    mock_session = mock_screen_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_screen_result_repository.get_by_exclusion_category(
        review_id, "Non-human subjects"
    )
    assert result == expected
    mock_session.exec.assert_called_once()


def test_screen_repository_get_by_response_metadata(mock_screen_result_repository):
    review_id = uuid.uuid4()
    metadata = {"model": "gpt-4", "temperature": 0.0}
    expected = [
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.CONSERVATIVE,
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
            response_metadata={"model": "gpt-4", "temperature": 0.0, "max_tokens": 100},
        ),
        ScreenAbstractResult(
            id=uuid.uuid4(),
            review_id=review_id,
            pubmed_result_id=uuid.uuid4(),
            strategy=ScreeningStrategyType.COMPREHENSIVE,
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.8,
            rationale="Test rationale",
            response_metadata={"model": "gpt-4", "temperature": 0.0, "max_tokens": 200},
        ),
    ]

    # Mock the session and query
    mock_session = mock_screen_result_repository.session_factory.begin.return_value
    mock_session.exec.return_value = expected

    result = mock_screen_result_repository.get_by_response_metadata(review_id, metadata)
    assert result == expected
    mock_session.exec.assert_called_once()


def test_systematic_review_repository_error_handling(mock_session_factory):
    repo = SystematicReviewRepository(mock_session_factory)
    mock_session = mock_session_factory()
    mock_session.__enter__.return_value = mock_session

    # Test get_with_pubmed_results error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Database error in get_with_pubmed_results for {repo.model_cls.__name__}: Database error",
    ):
        repo.get_with_pubmed_results(uuid.uuid4())

    # Test get_with_all error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Database error in get_with_all for {repo.model_cls.__name__}: Database error",
    ):
        repo.get_with_all(uuid.uuid4())


def test_pubmed_result_repository_error_handling(mock_pubmed_result_repository):
    mock_session = mock_pubmed_result_repository.session_factory()
    mock_session.__enter__.return_value = mock_session

    # Test get_by_review_id error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch PubMed results for review {mock_pubmed_result_repository.model_cls.__name__}: Database error",
    ):
        mock_pubmed_result_repository.get_by_review_id(uuid.uuid4())

    # Test get_with_screening_results error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch PubMed result with screening results for {mock_pubmed_result_repository.model_cls.__name__}: Database error",
    ):
        mock_pubmed_result_repository.get_with_screening_results(uuid.uuid4())

    # Test get_by_screening_strategy error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Database error in get_by_screening_strategy for {mock_pubmed_result_repository.model_cls.__name__}: Database error",
    ):
        mock_pubmed_result_repository.get_by_screening_strategy(
            uuid.uuid4(), ScreeningStrategyType.CONSERVATIVE
        )


def test_screen_result_repository_error_handling(mock_screen_result_repository):
    mock_session = mock_screen_result_repository.session_factory()
    mock_session.__enter__.return_value = mock_session

    # Test get_by_review_id error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch screening results for review {mock_screen_result_repository.model_cls.__name__}: Database error",
    ):
        mock_screen_result_repository.get_by_review_id(uuid.uuid4())

    # Test get_by_pubmed_result error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch screening results for PubMed result {mock_screen_result_repository.model_cls.__name__}: Database error",
    ):
        mock_screen_result_repository.get_by_pubmed_result(uuid.uuid4(), uuid.uuid4())

    # Test get_by_strategy error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch screening results by strategy for {mock_screen_result_repository.model_cls.__name__}: Database error",
    ):
        mock_screen_result_repository.get_by_strategy(
            uuid.uuid4(), ScreeningStrategyType.CONSERVATIVE
        )

    # Test get_by_exclusion_category error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch records by exclusion category for {mock_screen_result_repository.model_cls.__name__}: Database error",
    ):
        mock_screen_result_repository.get_by_exclusion_category(uuid.uuid4(), "test")

    # Test get_by_response_metadata error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch records by response metadata for {mock_screen_result_repository.model_cls.__name__}: Database error",
    ):
        mock_screen_result_repository.get_by_response_metadata(
            uuid.uuid4(), {"test": "test"}
        )


def test_log_repository_error_handling(mock_session_factory):
    repo = LogRepository(mock_session_factory)
    mock_session = mock_session_factory()
    mock_session.__enter__.return_value = mock_session

    # Test get_by_review_id error
    mock_session.exec.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to fetch records by review ID for {repo.model_cls.__name__}: Database error",
    ):
        repo.get_by_review_id(uuid.uuid4())

    # Test add error
    mock_session.add.side_effect = SQLAlchemyError("Database error")
    with pytest.raises(
        RepositoryError,
        match=f"Failed to add {repo.model_cls.__name__}: Database error",
    ):
        repo.add(
            LogRecord(
                id=uuid.uuid4(),
                timestamp=datetime.now(timezone.utc),
                review_id=uuid.uuid4(),
                level="INFO",
                message="test",
            )
        )

    # Test add_all error
    mock_session.add_all.side_effect = IntegrityError("statement", "params", "orig")
    with pytest.raises(
        ConstraintViolationError,
        match=f"Failed to add multiple {repo.model_cls.__name__} records: .*orig.*",
    ):
        repo.add_all(
            [
                LogRecord(
                    id=uuid.uuid4(),
                    timestamp=datetime.now(timezone.utc),
                    review_id=uuid.uuid4(),
                    level="INFO",
                    message="test",
                )
            ]
        )
