"""Synchronous repository implementations for systematic review models.

This module provides sync repositories for database operations using SQLModel.
Includes examples for common operations and complex queries including JSONB.

Note:
    - All repositories use a sessionmaker factory assigned as self.session_factory.

      ``with self.session_factory.begin() as session``:
            - Executes the query within a session.
            - Commits the session.
            - Rollbacks the session on error.
            - DOES NOT SUPPORT MANUAL ``session.commit()``
            - Cannot be used when ``session.refresh()`` is needed.

      ``with self.session_factory() as session``:
            - No automatic commit/rollback
            - Can use ``session.commit()`` and ``session.refresh()`` manually.

    - Async repositories are in `sr_assistant.core.repositories_async`.

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
    created = review_repo.add(review)

    # Get with relationships
    review_with_search_results = review_repo.get_with_search_results(created.id)
    full_review = review_repo.get_with_all(created.id) # search, logs, screening results
    full_review.search_results ...
    full_review.screen_abstract_results ...
    full_review.log_records ...
    ```
"""

from __future__ import annotations

import types
import typing as t
import uuid

from loguru import logger
from pydantic.types import JsonValue
from sqlalchemy import func, text
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlmodel import Session, and_, col, or_, select
from sqlmodel.sql.expression import SelectOfScalar

from sr_assistant.core.models import (
    Base,
    BenchmarkResultItem,
    BenchmarkRun,
    LogRecord,
    ScreenAbstractResult,
    ScreeningResolution,
    SearchResult,
    SystematicReview,
)
from sr_assistant.core.schemas import ExclusionReasons, SearchResultFilter
from sr_assistant.core.types import (
    LogLevel,
    ScreeningStrategyType,
    SearchDatabaseSource,
)

if t.TYPE_CHECKING:
    from collections.abc import Sequence


# Define a protocol for models with an ID
class ModelWithID(t.Protocol):
    id: uuid.UUID


class RepositoryError(Exception):
    """Base exception for persistence layer failures."""


class ConstraintViolationError(RepositoryError):
    """Database constraint violation (unique constraint, foreign key, etc.)."""


class RecordNotFoundError(RepositoryError):
    """Requested record was not found in the database."""


class BaseRepository[T: Base]:
    """Base repository implementing common sync database operations.

    Methods accept a Session provided by the calling service.
    Repositories do not manage transactions (commit/rollback).
    """

    @property
    def model_cls(self) -> type[T]:
        """Get the model class associated with the repository."""
        # Assumes direct inheritance like: class SubRepo(BaseRepository[ActualModel]):
        # Determine the generic base type (BaseRepository[ActualModel])
        generic_base = next(
            (
                base
                for base in types.get_original_bases(type(self))
                if t.get_origin(base) is BaseRepository
            ),
            None,
        )
        if generic_base is None:
            raise TypeError(
                f"Could not determine the generic base for {type(self).__name__}"
            )

        # Extract the type argument (ActualModel)
        model_arg = t.get_args(generic_base)[0]
        if not isinstance(model_arg, type):
            raise TypeError(
                f"Expected a type argument for BaseRepository, got {model_arg}"
            )
        return model_arg

    def _construct_get_stmt(self, id: uuid.UUID) -> SelectOfScalar[T]:
        Model = self.model_cls
        # Now T is bound to Base, which implicitly has 'id' via SQLModel
        return select(Model).where(Model.id == id)

    def get_by_id(self, session: Session, id: uuid.UUID) -> T | None:
        try:
            stmt = self._construct_get_stmt(id)
            return session.exec(stmt).first()
        except SQLAlchemyError as exc:
            msg = f"Database error in get_by_id for {self.model_cls.__name__}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_all(self, session: Session, limit: int | None = None) -> Sequence[T]:
        try:
            query = select(self.model_cls)
            if limit is not None:
                query = query.limit(limit)
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch all records for {self.model_cls.__name__}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def _construct_list_stmt(self, **filters: t.Any) -> SelectOfScalar[T]:
        Model = self.model_cls
        stmt = select(Model)
        where_clauses = []
        for c, v in filters.items():
            if not hasattr(Model, c):
                msg = f"Invalid column name {c} for model {Model.__name__}"
                logger.warning(msg)
                raise ValueError(msg)
            where_clauses.append(getattr(Model, c) == v)

        if where_clauses:
            stmt = stmt.where(and_(*where_clauses))
        return stmt

    def list(self, session: Session, **filters: t.Any) -> Sequence[T]:
        try:
            stmt = self._construct_list_stmt(**filters)
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to list {self.model_cls.__name__} with filters {filters}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def add(self, session: Session, record: T) -> T:
        try:
            session.add(record)
            session.flush([record])
            return record
        except IntegrityError as exc:
            msg = f"Constraint violation adding {self.model_cls.__name__}: {exc}"
            logger.exception(msg)
            raise ConstraintViolationError(msg) from exc
        except SQLAlchemyError as exc:
            msg = f"Database error adding {self.model_cls.__name__}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def add_all(self, session: Session, objs: Sequence[T]) -> Sequence[T]:
        try:
            session.add_all(objs)
            session.flush(objs)
            return objs
        except IntegrityError as exc:
            msg = f"Constraint violation adding multiple {self.model_cls.__name__} records: {exc}"
            logger.exception(msg)
            raise ConstraintViolationError(msg) from exc
        except SQLAlchemyError as exc:
            msg = f"Database error adding multiple {self.model_cls.__name__} records: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def update(self, session: Session, record: T) -> T:
        Model = self.model_cls
        record_id = getattr(record, "id", None)
        if not record_id:
            msg = f"{Model.__name__} object has no id, cannot update."
            raise ValueError(msg)

        try:
            existing = session.get(Model, record_id)
            if not existing:
                msg = f"{Model.__name__} with id {record_id} not found for update."
                logger.warning(msg)
                raise RecordNotFoundError(msg)

            merged_record = session.merge(record)
            session.flush([merged_record])
            return merged_record
        except IntegrityError as exc:
            msg = f"Constraint violation updating {self.model_cls.__name__} with id {record_id}: {exc}"
            logger.exception(msg)
            raise ConstraintViolationError(msg) from exc
        except SQLAlchemyError as exc:
            msg = f"Database error updating {self.model_cls.__name__} with id {record_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def delete(self, session: Session, id: uuid.UUID) -> None:
        try:
            record = self.get_by_id(session, id)
            if record is None:
                msg = f"{self.model_cls.__name__} with id {id} not found for deletion."
                logger.warning(msg)
                raise RecordNotFoundError(msg)

            session.delete(record)
            session.flush()
        except IntegrityError as exc:
            msg = f"Constraint violation deleting {self.model_cls.__name__} with id {id}: {exc}"
            logger.exception(msg)
            raise ConstraintViolationError(msg) from exc
        except SQLAlchemyError as exc:
            msg = (
                f"Database error deleting {self.model_cls.__name__} with id {id}: {exc}"
            )
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class SystematicReviewRepository(BaseRepository[SystematicReview]):
    """Repository for SystematicReview model operations.

    Examples:
        ```python
        repo = SystematicReviewRepository()

        # Get review with all relationships
        review = repo.get_with_search_results(review_id)

        # Access relationships
        for result in review.search_results:
            print(f"Title: {result.title}")
            if result.conservative_result:
                print(f"Conservative: {result.conservative_result.decision}")
            if result.comprehensive_result:
                print(f"Comprehensive: {result.comprehensive_result.decision}")
        ```
    """

    def get_with_search_results(
        self, session: Session, id: uuid.UUID
    ) -> SystematicReview | None:
        """Gets a review by ID. Relationship loading is handled by model defaults or service layer."""
        logger.debug("Calling get_by_id via get_with_search_results alias.")
        return self.get_by_id(session, id)

    def get_with_all(self, session: Session, id: uuid.UUID) -> SystematicReview | None:
        """Gets a review by ID. Relationship loading is handled by model defaults or service layer."""
        logger.debug("Calling get_by_id via get_with_all alias.")
        return self.get_by_id(session, id)


class SearchResultRepository(BaseRepository[SearchResult]):
    """Repository for SearchResult model operations."""

    def get_by_doi(self, session: Session, *, doi: str) -> SearchResult | None:
        """Retrieves a SearchResult by its DOI.

        Args:
            session: The database session.
            doi: The Digital Object Identifier to search for.

        Returns:
            The SearchResult if found, otherwise None.

        Raises:
            RepositoryError: If a database error occurs.
        """
        try:
            stmt = select(self.model_cls).where(self.model_cls.doi == doi)
            return session.exec(stmt).first()
        except SQLAlchemyError as exc:
            msg = f"Database error in get_by_doi for DOI '{doi}': {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_title_and_year(
        self, session: Session, *, title: str, year: str
    ) -> SearchResult | None:
        """Retrieves a SearchResult by its title and publication year.

        Args:
            session: The database session.
            title: The title to search for (case-insensitive exact match).
            year: The publication year (string) to search for.

        Returns:
            The SearchResult if found, otherwise None.

        Raises:
            RepositoryError: If a database error occurs.
        """
        try:
            # Using ilike for case-insensitive title matching, assuming PostgreSQL
            # For other DBs, lower(column) == lower(value) might be needed
            stmt = select(self.model_cls).where(
                # self.model_cls.title.ilike(title), # Hold off on ilike for now, depends on DB dialect features
                self.model_cls.title == title,  # Exact match for now
                self.model_cls.year == year,
            )
            return session.exec(stmt).first()
        except SQLAlchemyError as exc:
            msg = f"Database error in get_by_title_and_year for title '{title}' and year '{year}': {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[SearchResult]:
        """Get all SearchResults for a specific review."""
        try:
            query = (
                select(self.model_cls).where(self.model_cls.review_id == review_id)
                # Optionally load the review relationship if needed often
                # .options(selectinload(self.model_cls.review))
                # Remove loading of screening results
            )
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch SearchResults for review {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_source_details(
        self,
        session: Session,
        review_id: uuid.UUID,
        source_db: SearchDatabaseSource,
        source_id: str,
    ) -> SearchResult | None:
        """Get a specific SearchResult by its unique source identifiers within a review."""
        try:
            stmt = select(self.model_cls).where(
                self.model_cls.review_id == review_id,
                self.model_cls.source_db == source_db,
                self.model_cls.source_id == source_id,
            )
            return session.exec(stmt).first()
        except SQLAlchemyError as exc:
            msg = (
                f"Failed to fetch SearchResult for review {review_id}, "
                f"source {source_db.name}, id {source_id}: {exc}"
            )
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_existing_source_ids(
        self,
        session: Session,
        review_id: uuid.UUID,
        source_db: SearchDatabaseSource,
        source_ids: list[str],
    ) -> set[str]:
        """Checks a list of source_ids against the database for a given review_id and source_db,
        and returns a set of those that already exist.

        Args:
            session: The database session.
            review_id: The ID of the systematic review.
            source_db: The source database (e.g., PubMed, Scopus).
            source_ids: A list of source identifiers (e.g., PMIDs) to check.

        Returns:
            A set of source_ids from the input list that already exist in the database
            for the given review_id and source_db.

        Raises:
            RepositoryError: If a database error occurs.
        """
        if not source_ids:
            return set()
        try:
            stmt = (
                select(self.model_cls.source_id)
                .where(self.model_cls.review_id == review_id)
                .where(self.model_cls.source_db == source_db)
                .where(self.model_cls.source_id.in_(source_ids))  # type: ignore[attr-defined]
            )
            existing_ids = session.exec(stmt).all()
            return set(existing_ids)
        except SQLAlchemyError as exc:
            msg = (
                f"Database error in get_existing_source_ids for review {review_id}, "
                f"source {source_db.name}, with {len(source_ids)} source_ids: {exc}"
            )
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def advanced_search(
        self,
        session: Session,
        *,
        search_params: SearchResultFilter,
        skip: int = 0,
        limit: int = 100,
    ) -> list[SearchResult]:  # type: ignore
        """Retrieves SearchResults based on complex criteria with pagination.

        Args:
            session: The database session.
            search_params: A SearchResultFilter object with attributes to filter on.
            skip: Number of records to skip for pagination.
            limit: Maximum number of records to return.

        Returns:
            A list of SearchResult objects.

        Raises:
            RepositoryError: If a database error occurs.
        """
        try:
            query = select(self.model_cls)
            conditions = []

            if search_params.review_id is not None:
                conditions.append(self.model_cls.review_id == search_params.review_id)
            if search_params.source_db is not None:
                conditions.append(self.model_cls.source_db == search_params.source_db)
            if search_params.source_id is not None:
                conditions.append(self.model_cls.source_id == search_params.source_id)
            if search_params.doi is not None:
                conditions.append(self.model_cls.doi == search_params.doi)
            if search_params.title is not None:
                # Using `ilike` for case-insensitive partial matching (e.g., PostgreSQL).
                # For broader compatibility without specific DB features, one might use
                # from sqlalchemy import func
                # conditions.append(func.lower(self.model_cls.title).contains(search_params.title.lower()))
                conditions.append(
                    self.model_cls.title.ilike(f"%{search_params.title}%")
                )  # type: ignore[attr-defined]
            if search_params.year is not None:
                conditions.append(self.model_cls.year == search_params.year)

            if conditions:
                query = query.where(and_(*conditions))

            query = query.offset(skip).limit(limit)
            results = session.exec(query).all()
            return list(results)  # Ensure it's a list
        except SQLAlchemyError as exc:
            msg = f"Database error in advanced_search with params {search_params!r}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def count(
        self, session: Session, *, search_params: SearchResultFilter | None = None
    ) -> int:
        """Counts SearchResults, optionally filtered by search_params.

        Args:
            session: The database session.
            search_params: An optional SearchResultFilter object with attributes to filter on.

        Returns:
            The total count of matching SearchResult objects.

        Raises:
            RepositoryError: If a database error occurs.
        """
        try:
            stmt = select(func.count(text("*"))).select_from(self.model_cls)
            conditions = []

            if search_params:
                if search_params.review_id is not None:
                    conditions.append(
                        self.model_cls.review_id == search_params.review_id
                    )
                if search_params.source_db is not None:
                    conditions.append(
                        self.model_cls.source_db == search_params.source_db
                    )
                if search_params.source_id is not None:
                    conditions.append(
                        self.model_cls.source_id == search_params.source_id
                    )
                if search_params.doi is not None:
                    conditions.append(self.model_cls.doi == search_params.doi)
                if search_params.title is not None:
                    conditions.append(
                        self.model_cls.title.ilike(f"%{search_params.title}%")
                    )  # type: ignore[attr-defined]
                if search_params.year is not None:
                    conditions.append(self.model_cls.year == search_params.year)

            if conditions:
                stmt = stmt.where(and_(*conditions))

            count = session.exec(stmt).scalar_one()  # type: ignore[attr-defined]
            return count
        except SQLAlchemyError as exc:
            msg = f"Database error in count with params {search_params!r}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class ScreenAbstractResultRepository(BaseRepository[ScreenAbstractResult]):
    """Repository for ScreenAbstractResult model operations.

    Examples:
        ```python
        repo = ScreenAbstractResultRepository()

        # Get all results for a review
        results = repo.get_by_review_id(review_id)

        # Get results by strategy
        conservative = repo.get_by_strategy(
            review_id, ScreeningStrategyType.CONSERVATIVE
        )

        # Query by JSONB fields (exclusion reasons)
        wrong_population = repo.get_by_exclusion_category(
            review_id, "population_exclusion_reasons"
        )

        # Complex JSONB query
        results = repo.get_by_response_metadata(
            review_id, {"model": "gpt-4", "temperature": 0.0}
        )
        ```
    """

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[ScreenAbstractResult]:
        """Get all screening results for a review."""
        try:
            query = (
                select(self.model_cls).where(self.model_cls.review_id == review_id)
                # Remove loading - relationship commented out in model
                # .options(selectinload(self.model_cls.search_result))
            )
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch screening results for review {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_strategy(
        self,
        session: Session,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> Sequence[ScreenAbstractResult]:
        """Get screening results for a specific strategy."""
        try:
            query = (
                select(self.model_cls).where(
                    self.model_cls.review_id == review_id,
                    self.model_cls.screening_strategy == strategy,
                )
                # Remove loading - relationship commented out in model
                # .options(selectinload(self.model_cls.search_result))
            )
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch screening results by strategy for {review_id}, {strategy}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_exclusion_category(
        self,
        session: Session,
        review_id: uuid.UUID,
        exclusion_reason: str,
    ) -> Sequence[ScreenAbstractResult]:
        """Get results where any exclusion category list contains the given reason."""
        try:
            category_fields = [
                field_name
                for field_name in ExclusionReasons.model_fields.keys()
                if field_name.endswith("_exclusion_reasons")
            ]
            category_conditions = [
                col(self.model_cls.exclusion_reason_categories)[field].op("@>")(
                    f'["{exclusion_reason}"]'
                )
                for field in category_fields
            ]
            query = (
                select(self.model_cls).where(
                    and_(
                        self.model_cls.review_id == review_id,
                        or_(*category_conditions),
                    )
                )
                # Remove loading - relationship commented out in model
                # .options(selectinload(self.model_cls.search_result))
            )
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch records by exclusion category for {review_id}, {exclusion_reason}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_response_metadata(
        self,
        session: Session,
        review_id: uuid.UUID,
        metadata: dict[str, JsonValue],
    ) -> Sequence[ScreenAbstractResult]:
        """Get results by matching response_metadata JSONB fields."""
        try:
            query = select(self.model_cls).where(
                self.model_cls.review_id == review_id,
                *[
                    col(self.model_cls.response_metadata)[key] == value
                    for key, value in metadata.items()
                ],
            )
            # No relationship needed usually for this query
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch records by response metadata for {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_review_strategy(
        self,
        session: Session,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> Sequence[ScreenAbstractResult]:
        """Gets screening results for a specific review and strategy."""
        try:
            stmt = select(self.model_cls).where(
                self.model_cls.review_id == review_id,
                self.model_cls.screening_strategy == strategy,
            )
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch ScreenAbstractResult by review/strategy: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class LogRepository(BaseRepository[LogRecord]):
    """Repository for log records."""

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[LogRecord]:
        """Get log records by review ID."""
        try:
            query = select(self.model_cls).where(self.model_cls.review_id == review_id)
            # Remove loading - relationship commented out in model
            # return session.exec(
            #     query.options(selectinload(self.model_cls.review))
            # ).all()
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch records by review ID for {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_level(
        self,
        session: Session,
        level: LogLevel,
        review_id: uuid.UUID | None = None,
    ) -> Sequence[LogRecord]:
        """Get log records by level and optional review ID."""
        try:
            if review_id is None:
                query = select(self.model_cls).where(self.model_cls.level == level)
            else:
                query = select(self.model_cls).where(
                    self.model_cls.review_id == review_id,
                    self.model_cls.level == level,
                )
            # Remove loading - relationship commented out in model
            # return session.exec(
            #     query.options(selectinload(self.model_cls.review))
            # ).all()
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch records by level for {review_id}, {level}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class ScreeningResolutionRepository(BaseRepository[ScreeningResolution]):
    """Repository for ScreeningResolution model operations.

    Handles saving and retrieving resolution decisions.
    """

    def get_by_search_result_id(
        self, session: Session, search_result_id: uuid.UUID
    ) -> ScreeningResolution | None:
        """Get a ScreeningResolution by SearchResult ID."""
        try:
            stmt = select(self.model_cls).where(
                self.model_cls.search_result_id == search_result_id
            )
            return session.exec(stmt).first()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch ScreeningResolution for search result {search_result_id}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[ScreeningResolution]:
        """Get all ScreeningResolutions for a specific review."""
        try:
            stmt = select(self.model_cls).where(self.model_cls.review_id == review_id)
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch ScreeningResolutions for review {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class BenchmarkRunRepository(BaseRepository["BenchmarkRun"]):  # type: ignore[name-defined]
    """Repository for BenchmarkRun model operations."""

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[BenchmarkRun]:  # type: ignore[name-defined]
        """Get all BenchmarkRuns for a specific review."""
        try:
            stmt = select(self.model_cls).where(self.model_cls.review_id == review_id)
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch BenchmarkRuns for review {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc


class BenchmarkResultItemRepository(BaseRepository["BenchmarkResultItem"]):  # type: ignore[name-defined]
    """Repository for BenchmarkResultItem model operations."""

    def get_by_benchmark_run_id(
        self, session: Session, benchmark_run_id: uuid.UUID
    ) -> Sequence[BenchmarkResultItem]:  # type: ignore[name-defined]
        """Get all BenchmarkResultItems for a specific benchmark run."""
        try:
            stmt = select(self.model_cls).where(
                self.model_cls.benchmark_run_id == benchmark_run_id
            )
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch BenchmarkResultItems for benchmark run {benchmark_run_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_search_result_id(
        self, session: Session, search_result_id: uuid.UUID
    ) -> Sequence[BenchmarkResultItem]:  # type: ignore[name-defined]
        """Get all BenchmarkResultItems for a specific search result."""
        try:
            stmt = select(self.model_cls).where(
                self.model_cls.search_result_id == search_result_id
            )
            return session.exec(stmt).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch BenchmarkResultItems for search result {search_result_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc
