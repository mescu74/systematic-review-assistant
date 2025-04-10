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

import typing as t

from loguru import logger
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlmodel import Session, and_, col, or_, select

from sr_assistant.core import models
from sr_assistant.core.models import (
    Base,
    LogRecord,
    ScreenAbstractResult,
    ScreeningResolution,
    SystematicReview,
)
from sr_assistant.core.schemas import ExclusionReasons
from sr_assistant.core.types import (
    LogLevel,
    ScreeningStrategyType,
)

if t.TYPE_CHECKING:
    import uuid
    from collections.abc import Sequence

    from pydantic.types import JsonValue
    from sqlmodel.sql.expression import SelectOfScalar


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

    _model_cls: type[T] | None = None  # Cache for resolved model class

    @property
    def model_cls(self) -> type[T]:
        """Get the model class associated with the repository, resolving forward refs."""
        if self._model_cls is None:
            orig_bases = getattr(type(self), "__orig_bases__", [])
            model_type_arg = None
            for base in orig_bases:
                if getattr(base, "__origin__", None) is BaseRepository:
                    args = getattr(base, "__args__", [])
                    if args:
                        model_type_arg = args[0]
                        break

            if model_type_arg is None:
                raise TypeError(
                    "Could not determine model class argument for repository"
                )

            # Resolve the type argument (which might be a string/ForwardRef)
            resolved_type = None
            if isinstance(model_type_arg, type) and issubclass(model_type_arg, Base):
                resolved_type = model_type_arg  # It was already a resolved type
            elif isinstance(model_type_arg, (str, t.ForwardRef)):
                # Need to resolve string or ForwardRef
                try:
                    # Use get_type_hints on a temporary namespace
                    temp_globals = models.__dict__.copy()
                    temp_locals = (
                        locals()
                    )  # Include local scope if needed, though usually not for models

                    # Directly evaluate the ForwardRef if possible, or use get_type_hints
                    if isinstance(model_type_arg, t.ForwardRef):
                        # Pass recursive_guard as an empty set for older Python versions compatibility if needed
                        # For Python 3.9+, frozenset() is correct.
                        resolved_type = model_type_arg._evaluate(
                            temp_globals, temp_locals, frozenset()
                        )
                    elif isinstance(model_type_arg, str):
                        # Use eval if it's just a string name
                        resolved_type = eval(model_type_arg, temp_globals, temp_locals)

                    if not (
                        isinstance(resolved_type, type)
                        and issubclass(resolved_type, Base)
                    ):
                        raise TypeError(
                            f"Resolved type {resolved_type!r} is not a valid Base model subclass."
                        )

                except NameError as e:
                    raise TypeError(
                        f"Could not resolve type name '{model_type_arg}': {e}"
                    ) from e
                except Exception as e:
                    # Catch potential errors during resolution
                    raise TypeError(
                        f"Error resolving type argument {model_type_arg!r}: {e}"
                    ) from e
            else:
                raise TypeError(
                    f"Unexpected type argument for repository: {model_type_arg!r}"
                )

            if resolved_type is None:
                raise TypeError(
                    f"Failed to resolve model class from {model_type_arg!r}"
                )

            self._model_cls = resolved_type

        return self._model_cls

    def _construct_get_stmt(self, id: uuid.UUID) -> SelectOfScalar[T]:
        Model = self.model_cls
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
            msg = f"Constraint violation updating {Model.__name__} with id {record_id}: {exc}"
            logger.exception(msg)
            raise ConstraintViolationError(msg) from exc
        except SQLAlchemyError as exc:
            msg = f"Database error updating {Model.__name__} with id {record_id}: {exc}"
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


class SearchResultRepository(BaseRepository["SearchResult"]):
    """Repository for SearchResult model operations."""

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


class ScreenAbstractResultRepository(BaseRepository[ScreenAbstractResult]):
    """Repository for ScreenAbstractResult model operations.

    Examples:
        ```python
        repo = ScreenAbstractResultRepository()

        # Get all results for a review
        results = repo.get_by_review_id(review_id)

        # Get results by strategy
        conservative = repo.get_by_strategy(review_id, ScreeningStrategyType.CONSERVATIVE)

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
        """Get resolution by SearchResult ID."""
        try:
            query = (
                select(self.model_cls).where(
                    self.model_cls.search_result_id == search_result_id
                )
                # .options(selectinload(self.model_cls.search_result))
            )
            return session.exec(query).first()
        except SQLAlchemyError as exc:
            msg = f"Database error in get_by_search_result_id for {search_result_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc

    def get_by_review_id(
        self, session: Session, review_id: uuid.UUID
    ) -> Sequence[ScreeningResolution]:
        """Get all resolutions for a specific review."""
        try:
            query = (
                select(self.model_cls).where(self.model_cls.review_id == review_id)
                # .options(
                #     selectinload(self.model_cls.review),
                #     selectinload(self.model_cls.search_result),
                # )
            )
            return session.exec(query).all()
        except SQLAlchemyError as exc:
            msg = f"Failed to fetch screening resolutions for review {review_id}: {exc}"
            logger.exception(msg)
            raise RepositoryError(msg) from exc
