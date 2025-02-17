"""Synchronous repository implementations for systematic review models.

This module provides sync repositories for database operations using SQLModel.
Includes examples for common operations and complex queries including JSONB.

Examples:
    Basic repository usage:
    ```python
    # Initialize repository
    review_repo = SystematicReviewRepository()

    # Create new review
    review = SystematicReview(
        research_question="What is the effectiveness of...",
        inclusion_criteria="RCTs published after 2010...",
        exclusion_criteria="Case studies, non-English...",
    )
    created = review_repo.create(review)

    # Get with relationships
    full_review = review_repo.get_with_pubmed_results(created.id)
    ```
"""

from __future__ import annotations

import typing as t
import uuid

from pydantic.types import JsonValue
from loguru import logger
from sqlmodel import Session, select
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import selectinload

from sr_assistant.app.database import SQLModelSession, session_factory, sessionmaker
from sr_assistant.core.models import (
    PubMedResult,
    ScreenAbstractResultModel,
    SystematicReview,
)
from sr_assistant.core.types import ScreeningStrategyType


type RepositoryType = SystematicReview | PubMedResult | ScreenAbstractResultModel


class BaseRepository[T: RepositoryType]:
    # class BaseRepository[T: SystematicReview | PubMedResult | ScreenAbstractResultModel]:
    """Base repository implementing common database operations.

    Examples:
        ```python
        repo = BaseRepository[SystematicReview]()

        # Get by ID
        review = repo.get_by_id(uuid.UUID("..."))

        # Get all with limit
        reviews = repo.get_all(limit=10)

        # Create
        new_review = SystematicReview(...)
        created = repo.create(new_review)

        # Update
        updated = repo.update(review.id, review)

        # Delete
        deleted = repo.delete(review.id)
        ```
    """

    def __init__(
        self,
        session_factory: sessionmaker[SQLModelSession] | None = None,
    ) -> None:
        """Initialize repository with session factory."""
        self.session_factory = session_factory

    model_class: t.ClassVar[t.type[T]]

    def __init__(self, session_factory: SQLModelSession | None = None) -> None:
        """Initialize repository with session factory."""
        self.get_session = get_session or Session

    def get_by_id(self, id: uuid.UUID) -> T | None:
        """Get a record by its ID."""
        try:
            with self.session_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                return session.exec(query).first()
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_id: {e}")
            raise

    def get_all(self, limit: int | None = None) -> list[T]:
        """Get all records with optional limit."""
        try:
            with self.session_factory.begin() as session:
                query = select(self.model_class)
                if limit is not None:
                    query = query.limit(limit)
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_all: {e}")
            raise

    def create(self, model: T) -> T:
        """Create a new record."""
        try:
            with self.session_factory.begin() as session:
                session.add(model)
                session.commit()
                session.refresh(model)
                return model
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in create: {e}")
            raise

    def update(self, id: uuid.UUID, model: T) -> T:
        """Update an existing record."""
        try:
            with self.session_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                db_model = session.exec(query).first()

                if db_model is None:
                    msg = f"{self.model_class.__name__} with ID {id} not found"
                    raise ValueError(msg)

                # Update fields excluding id and timestamps
                model_data = model.model_dump(
                    exclude={"id", "created_at", "updated_at"}
                )
                for key, value in model_data.items():
                    setattr(db_model, key, value)

                session.add(db_model)
                session.commit()
                session.refresh(db_model)
                return db_model
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in update: {e}")
            raise

    def delete(self, id: uuid.UUID) -> bool:
        """Delete a record by ID."""
        try:
            with self.session_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                model = session.exec(query).first()

                if model is None:
                    return False

                session.delete(model)
                session.commit()
                return True
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in delete: {e}")
            raise


class SystematicReviewRepository(BaseRepository[SystematicReview]):
    """Repository for SystematicReview model operations.

    Examples:
        ```python
        repo = SystematicReviewRepository()

        # Get review with all relationships
        review = repo.get_with_pubmed_results(review_id)

        # Access relationships
        for result in review.pubmed_results:
            print(f"Title: {result.title}")
            if result.conservative_result:
                print(f"Conservative: {result.conservative_result.decision}")
            if result.comprehensive_result:
                print(f"Comprehensive: {result.comprehensive_result.decision}")
        ```
    """

    model_class = SystematicReview

    def get_with_pubmed_results(self, id: uuid.UUID) -> SystematicReview | None:
        """Get a review with its PubMed results and screening results loaded."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(SystematicReview)
                    .where(SystematicReview.id == id)
                    .options(
                        selectinload(SystematicReview.pubmed_results),
                        selectinload(SystematicReview.screen_abstract_results),
                    )
                )
                return session.exec(query).first()
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_with_pubmed_results: {e}")
            raise


class PubMedResultRepository(BaseRepository[PubMedResult]):
    """Repository for PubMedResult model operations.

    Examples:
        ```python
        repo = PubMedResultRepository()

        # Get all results for a review
        results = repo.get_by_review_id(review_id)

        # Get with screening results
        result = repo.get_with_screening_results(pubmed_id)

        # Get results by screening strategy
        conservative = repo.get_by_screening_strategy(
            review_id, ScreeningStrategyType.CONSERVATIVE
        )
        ```
    """

    model_class = PubMedResult

    def get_by_review_id(self, review_id: uuid.UUID) -> list[PubMedResult]:
        """Get all PubMed results for a review."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(PubMedResult)
                    .where(PubMedResult.review_id == review_id)
                    .options(
                        selectinload(PubMedResult.conservative_result),
                        selectinload(PubMedResult.comprehensive_result),
                    )
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_review_id: {e}")
            raise

    def get_with_screening_results(self, id: uuid.UUID) -> PubMedResult | None:
        """Get a PubMed result with both screening results loaded."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(PubMedResult)
                    .where(PubMedResult.id == id)
                    .options(
                        selectinload(PubMedResult.conservative_result),
                        selectinload(PubMedResult.comprehensive_result),
                    )
                )
                return session.exec(query).first()
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_with_screening_results: {e}")
            raise

    def get_by_screening_strategy(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> list[PubMedResult]:
        """Get PubMed results that have a specific screening strategy result."""
        try:
            with self.session_factory.begin() as session:
                query = select(PubMedResult).where(PubMedResult.review_id == review_id)

                if strategy == ScreeningStrategyType.CONSERVATIVE:
                    query = query.where(
                        PubMedResult.conservative_result_id.is_not(None)
                    )
                    query = query.options(
                        selectinload(PubMedResult.conservative_result)
                    )
                else:
                    query = query.where(
                        PubMedResult.comprehensive_result_id.is_not(None)
                    )
                    query = query.options(
                        selectinload(PubMedResult.comprehensive_result)
                    )

                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_screening_strategy: {e}")
            raise


class ScreenAbstractResultRepository(BaseRepository[ScreenAbstractResultModel]):
    """Repository for ScreenAbstractResult model operations.

    Examples:
        ```python
        repo = ScreenAbstractResultRepository()

        # Get all results for a review
        results = repo.get_by_review_id(review_id)

        # Get results by strategy
        conservative = repo.get_by_strategy(review_id, ScreeningStrategyType.CONSERVATIVE)

        # Query by JSONB fields (exclusion reasons)
        wrong_population = repo.get_by_exclusion_category(review_id, "wrong_population")

        # Complex JSONB query
        results = repo.get_by_response_metadata(
            review_id, {"model": "gpt-4", "temperature": 0.0}
        )
        ```
    """

    model_class = ScreenAbstractResultModel

    def get_by_review_id(self, review_id: uuid.UUID) -> list[ScreenAbstractResultModel]:
        """Get all screening results for a review."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(ScreenAbstractResultModel.review_id == review_id)
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_review_id: {e}")
            raise

    def get_by_pubmed_result(
        self,
        review_id: uuid.UUID,
        pubmed_result_id: uuid.UUID,
    ) -> list[ScreenAbstractResultModel]:
        """Get screening results for a specific PubMed result."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(
                        ScreenAbstractResultModel.review_id == review_id,
                        ScreenAbstractResultModel.pubmed_result_id == pubmed_result_id,
                    )
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_pubmed_result: {e}")
            raise

    def get_by_strategy(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> list[ScreenAbstractResultModel]:
        """Get screening results for a specific strategy."""
        try:
            with self.session_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(
                        ScreenAbstractResultModel.review_id == review_id,
                        ScreenAbstractResultModel.screening_strategy == strategy,
                    )
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_strategy: {e}")
            raise

    def get_by_exclusion_category(
        self,
        review_id: uuid.UUID,
        category: str,
    ) -> list[ScreenAbstractResultModel]:
        """Get results where exclusion_reason_categories JSONB contains category.

        Example:
            ```python
            # Find all results excluded due to wrong population
            results = repo.get_by_exclusion_category(review_id, "wrong_population")
            ```
        """
        try:
            with self.session_factory.begin() as session:
                query = select(ScreenAbstractResultModel).where(
                    ScreenAbstractResultModel.review_id == review_id,
                    ScreenAbstractResultModel.exclusion_reason_categories[
                        category
                    ].is_not(None),
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_exclusion_category: {e}")
            raise

    def get_by_response_metadata(
        self,
        review_id: uuid.UUID,
        metadata: dict[str, JsonValue],
    ) -> list[ScreenAbstractResultModel]:
        """Get results by matching response_metadata JSONB fields.

        Example:
            ```python
            # Find all results from a specific model configuration
            results = repo.get_by_response_metadata(
                review_id, {"model": "gpt-4", "temperature": 0.0}
            )
            ```
        """
        try:
            with self.session_factory.begin() as session:
                query = select(ScreenAbstractResultModel).where(
                    ScreenAbstractResultModel.review_id == review_id,
                    *[
                        ScreenAbstractResultModel.response_metadata[key] == value
                        for key, value in metadata.items()
                    ],
                )
                return list(session.exec(query))
        except SQLAlchemyError as e:
            loggger.exception(f"Database error in get_by_response_metadata: {e}")
            raise
