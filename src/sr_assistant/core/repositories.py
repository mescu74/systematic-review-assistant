"""Repository implementations for systematic review models.

This module provides async repositories for database operations using SQLModel.
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
    created = await review_repo.create(review)

    # Get with relationships
    full_review = await review_repo.get_with_pubmed_results(created.id)

    # Update
    review.inclusion_criteria = "Updated criteria..."
    updated = await review_repo.update(review.id, review)
    ```

    Working with screening results:
    ```python
    # Get both conservative and comprehensive results
    pubmed_repo = PubMedResultRepository()
    result = await pubmed_repo.get_with_screening_results(pubmed_id)

    if result.conservative_result:
        print(f"Conservative decision: {result.conservative_result.decision}")
    if result.comprehensive_result:
        print(f"Comprehensive decision: {result.comprehensive_result.decision}")

    # Query by JSONB fields
    screen_repo = ScreenAbstractResultRepository()
    results = await screen_repo.get_by_exclusion_category(
        review_id, category="wrong_population"
    )
    ```
"""

from __future__ import annotations

import typing as t
import uuid

from loguru import logger
from pydantic.types import JsonValue
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.ext.asyncio import async_sessionmaker
from sqlalchemy.orm import selectinload
from sqlmodel import select

from sr_assistant.app.database import AsyncSQLModelSession, asession_factory
from sr_assistant.core.models import (
    PubMedResult,
    ScreenAbstractResultModel,
    SystematicReview,
)
from sr_assistant.core.types import ScreeningStrategyType

type RepositoryType = SystematicReview | PubMedResult | ScreenAbstractResultModel


class BaseRepository[T: RepositoryType]:
    """Base repository implementing common database operations.

    Examples:
        ```python
        repo = BaseRepository[SystematicReview]()

        # Get by ID
        review = await repo.get_by_id(uuid.UUID("..."))

        # Get all with limit
        reviews = await repo.get_all(limit=10)

        # Create
        new_review = SystematicReview(...)
        created = await repo.create(new_review)

        # Update
        updated = await repo.update(review.id, review)

        # Delete
        deleted = await repo.delete(review.id)
        ```
    """

    model_class: t.ClassVar[t.type[T]]

    def __init__(
        self,
        asession_factory: async_sessionmaker[AsyncSQLModelSession] = asession_factory,
    ) -> None:
        """Initialize repository with session factory."""
        self.asession_factory = asession_factory

    async def get_by_id(self, id: uuid.UUID) -> T | None:
        """Get a record by its ID."""
        try:
            async with self.asession_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                return await session.exec(query).first()
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_id: {e}")
            raise

    async def get_all(self, limit: int | None = None) -> list[T]:
        """Get all records with optional limit."""
        try:
            async with self.asession_factory.begin() as session:
                query = select(self.model_class)
                if limit is not None:
                    query = query.limit(limit)
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_exclusion_category: {e}")
            raise

    async def get_by_response_metadata(
        self,
        review_id: uuid.UUID,
        metadata: dict[str, JsonValue],
    ) -> list[ScreenAbstractResultModel]:
        """Get results by matching response_metadata JSONB fields.

        Example:
            ```python
            # Find all results from a specific model configuration
            results = await repo.get_by_response_metadata(
                review_id, {"model": "gpt-4", "temperature": 0.0}
            )
            ```
        """
        try:
            async with self.asession_factory.begin() as session:
                query = select(ScreenAbstractResultModel).where(
                    ScreenAbstractResultModel.review_id == review_id,
                    *[
                        ScreenAbstractResultModel.response_metadata[key] == value
                        for key, value in metadata.items()
                    ],
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_response_metadata: {e}")
            raise

    async def create(self, model: T) -> T:
        """Create a new record."""
        try:
            async with self.asession_factory.begin() as session:
                session.add(model)
                await session.commit()
                await session.refresh(model)
                return model
        except SQLAlchemyError as e:
            logger.exception(f"Database error in create: {e}")
            raise

    async def update(self, id: uuid.UUID, model: T) -> T:
        """Update an existing record."""
        try:
            async with self.asession_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                db_model = await session.exec(query).first()

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
                await session.commit()
                await session.refresh(db_model)
                return db_model
        except SQLAlchemyError as e:
            logger.exception(f"Database error in update: {e}")
            raise

    async def delete(self, id: uuid.UUID) -> bool:
        """Delete a record by ID."""
        try:
            async with self.asession_factory.begin() as session:
                query = select(self.model_class).where(self.model_class.id == id)
                model = await session.exec(query).first()

                if model is None:
                    return False

                await session.delete(model)
                await session.commit()
                return True
        except SQLAlchemyError as e:
            logger.exception(f"Database error in delete: {e}")
            raise


class SystematicReviewRepository(BaseRepository[SystematicReview]):
    """Repository for SystematicReview model operations.

    Examples:
        ```python
        repo = SystematicReviewRepository()

        # Get review with all relationships
        review = await repo.get_with_pubmed_results(review_id)

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

    async def get_with_pubmed_results(self, id: uuid.UUID) -> SystematicReview | None:
        """Get a review with its PubMed results and screening results loaded."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(SystematicReview)
                    .where(SystematicReview.id == id)
                    .options(
                        selectinload(SystematicReview.pubmed_results),
                        selectinload(SystematicReview.screen_abstract_results),
                    )
                )
                return await session.exec(query).first()
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_with_pubmed_results: {e}")
            raise


class PubMedResultRepository(BaseRepository[PubMedResult]):
    """Repository for PubMedResult model operations.

    Examples:
        ```python
        repo = PubMedResultRepository()

        # Get all results for a review
        results = await repo.get_by_review_id(review_id)

        # Get with screening results
        result = await repo.get_with_screening_results(pubmed_id)

        # Get results by screening strategy
        conservative = await repo.get_by_screening_strategy(
            review_id, ScreeningStrategyType.CONSERVATIVE
        )
        comprehensive = await repo.get_by_screening_strategy(
            review_id, ScreeningStrategyType.COMPREHENSIVE
        )
        ```
    """

    model_class = PubMedResult

    async def get_by_review_id(self, review_id: uuid.UUID) -> list[PubMedResult]:
        """Get all PubMed results for a review."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(PubMedResult)
                    .where(PubMedResult.review_id == review_id)
                    .options(
                        selectinload(PubMedResult.conservative_result),
                        selectinload(PubMedResult.comprehensive_result),
                    )
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_review_id: {e}")
            raise

    async def get_with_screening_results(self, id: uuid.UUID) -> PubMedResult | None:
        """Get a PubMed result with both screening results loaded."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(PubMedResult)
                    .where(PubMedResult.id == id)
                    .options(
                        selectinload(PubMedResult.conservative_result),
                        selectinload(PubMedResult.comprehensive_result),
                    )
                )
                return await session.exec(query).first()
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_with_screening_results: {e}")
            raise

    async def get_by_screening_strategy(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> list[PubMedResult]:
        """Get PubMed results that have a specific screening strategy result."""
        try:
            async with self.asession_factory.begin() as session:
                # Build base query
                query = select(PubMedResult).where(PubMedResult.review_id == review_id)

                # Add join based on strategy
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

                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_screening_strategy: {e}")
            raise


class ScreenAbstractResultRepository(BaseRepository[ScreenAbstractResultModel]):
    """Repository for ScreenAbstractResult model operations.

    Examples:
        ```python
        repo = ScreenAbstractResultRepository()

        # Get all results for a review
        results = await repo.get_by_review_id(review_id)

        # Get results by strategy
        conservative = await repo.get_by_strategy(
            review_id, ScreeningStrategyType.CONSERVATIVE
        )

        # Query by JSONB fields (exclusion reasons)
        wrong_population = await repo.get_by_exclusion_category(
            review_id, "wrong_population"
        )

        # Complex JSONB query
        results = await repo.get_by_response_metadata(
            review_id, {"model": "gpt-4", "temperature": 0.0}
        )
        ```
    """

    model_class = ScreenAbstractResultModel

    async def get_by_review_id(
        self, review_id: uuid.UUID
    ) -> list[ScreenAbstractResultModel]:
        """Get all screening results for a review."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(ScreenAbstractResultModel.review_id == review_id)
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_review_id: {e}")
            raise

    async def get_by_pubmed_result(
        self,
        review_id: uuid.UUID,
        pubmed_result_id: uuid.UUID,
    ) -> list[ScreenAbstractResultModel]:
        """Get screening results for a specific PubMed result."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(
                        ScreenAbstractResultModel.review_id == review_id,
                        ScreenAbstractResultModel.pubmed_result_id == pubmed_result_id,
                    )
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_pubmed_result: {e}")
            raise

    async def get_by_strategy(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> list[ScreenAbstractResultModel]:
        """Get screening results for a specific strategy."""
        try:
            async with self.asession_factory.begin() as session:
                query = (
                    select(ScreenAbstractResultModel)
                    .where(
                        ScreenAbstractResultModel.review_id == review_id,
                        ScreenAbstractResultModel.screening_strategy == strategy,
                    )
                    .options(selectinload(ScreenAbstractResultModel.pubmed_result))
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_strategy: {e}")
            raise

    async def get_by_exclusion_category(
        self,
        review_id: uuid.UUID,
        category: str,
    ) -> list[ScreenAbstractResultModel]:
        """Get results where exclusion_reason_categories JSONB contains category.

        Example:
            ```python
            # Find all results excluded due to wrong population
            results = await repo.get_by_exclusion_category(review_id, "wrong_population")
            ```
        """
        try:
            async with self.asession_factory.begin() as session:
                query = select(ScreenAbstractResultModel).where(
                    ScreenAbstractResultModel.review_id == review_id,
                    ScreenAbstractResultModel.exclusion_reason_categories[
                        category
                    ].is_not(None),
                )

                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_exclusion_category: {e}")
            raise

    async def get_by_response_metadata(
        self,
        review_id: uuid.UUID,
        metadata: dict[str, JsonValue],
    ) -> list[ScreenAbstractResultModel]:
        """Get results by matching response_metadata JSONB fields.

        Example:
            ```python
            # Find all results from a specific model configuration
            results = await repo.get_by_response_metadata(
                review_id, {"model": "gpt-4", "temperature": 0.0}
            )
            ```
        """
        try:
            async with self.asession_factory.begin() as session:
                query = select(ScreenAbstractResultModel).where(
                    ScreenAbstractResultModel.review_id == review_id,
                    *[
                        ScreenAbstractResultModel.response_metadata[key] == value
                        for key, value in metadata.items()
                    ],
                )
                return list(await session.exec(query))
        except SQLAlchemyError as e:
            logger.exception(f"Database error in get_by_response_metadata: {e}")
            raise
