"""Application Service Layer - Synchronous Implementation."""

from __future__ import annotations

import contextlib
import re
import typing as t
import uuid
from collections.abc import Iterator
from datetime import datetime

from loguru import logger
from sqlalchemy.orm import sessionmaker
from sqlmodel import Session

from sr_assistant.app.database import session_factory
from sr_assistant.core import models, repositories, schemas
from sr_assistant.core.types import (
    ScreeningStrategyType,
    SearchDatabaseSource,
)

if t.TYPE_CHECKING:
    from collections.abc import Sequence


# Define potential service-level errors
class ServiceError(Exception):
    """Base exception for service layer errors."""


class RecordMappingError(ServiceError):
    """Error during mapping of external data (e.g., API record) to internal model."""


class BaseService:
    """Base service providing session management (sync)."""

    def __init__(self, factory: sessionmaker[Session] = session_factory):
        self.session_factory = factory


class SearchService(BaseService):
    """Service for managing SearchResult data (sync)."""

    def __init__(
        self,
        factory: sessionmaker[Session] = session_factory,
        search_repo: repositories.SearchResultRepository | None = None,
    ):
        super().__init__(factory)
        self.search_repo = search_repo or repositories.SearchResultRepository()

    # --- Helper context manager for session handling ---
    @contextlib.contextmanager
    def _get_session(self, session: Session | None = None) -> Iterator[Session]:
        """Provides a session, either the one passed in or a new one from the factory."""
        if session:
            # If session is provided, yield it directly without context management
            yield session
        else:
            # If no session provided, create one using the factory's context manager
            with self.session_factory() as new_session:
                yield new_session

    # --- Mapping functions (to be implemented) ---
    def _map_pubmed_to_search_result(
        self, review_id: uuid.UUID, api_record: dict
    ) -> models.SearchResult | None:
        """Maps a raw PubMed API record (parsed XML dict) to a SearchResult model."""
        pmid: str = ""
        try:
            # Navigate through the nested dictionary structure
            medline = api_record.get("MedlineCitation", {})
            if not medline:
                return None  # Essential part missing

            pubmed_data = api_record.get("PubmedData", {})
            article_info = medline.get("Article", {})
            if not article_info:
                return None  # Essential part missing

            journal_info = article_info.get("Journal", {})
            journal_issue = journal_info.get("JournalIssue", {})
            pub_date = journal_issue.get("PubDate", {})

            # --- Extract Core IDs ---
            pmid_elem = medline.get("PMID")
            # PMID element might have attributes, get text content
            pmid = str(pmid_elem) if pmid_elem is not None else ""
            if not pmid:
                logger.warning("Skipping record: Missing PMID")
                return None

            pmc = ""
            doi = ""
            # Handle ArticleIdList structure variations
            article_id_list = pubmed_data.get("ArticleIdList", [])
            if isinstance(article_id_list, list):
                for v in article_id_list:
                    # Ensure v is a dict-like object with attributes
                    if hasattr(v, "attributes") and hasattr(v, "__str__"):
                        try:
                            id_type = v.attributes.get("IdType", "").lower()
                            id_value = str(v)
                            if id_type == "pmc":
                                pmc = id_value
                            elif id_type == "doi":
                                doi = id_value
                        except Exception as e:
                            logger.warning(f"Skipping malformed ArticleId {v}: {e}")
                            continue
                    else:
                        logger.warning(
                            f"Skipping unexpected item in ArticleIdList: {v}"
                        )

            # --- Extract Title ---
            title_elem = article_info.get("ArticleTitle")
            # Title element might have attributes or be complex, get text content
            title = str(title_elem) if title_elem is not None else ""
            if not title:
                logger.warning(f"Skipping record {pmid}: Missing Title")
                return None

            # --- Extract Abstract ---
            abstract_text_parts = []
            abstract_section = article_info.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_section, list):
                for part in abstract_section:
                    # Handle StringElement with attributes or plain strings
                    if hasattr(part, "attributes"):
                        label = part.attributes.get("Label", None)
                        text = str(part) if part is not None else ""
                        if label:
                            abstract_text_parts.append(f"{label}: {text}")
                        else:
                            abstract_text_parts.append(text)
                    elif isinstance(part, str):
                        abstract_text_parts.append(part)
            elif abstract_section:  # Handle single StringElement or string
                abstract_text_parts.append(str(abstract_section))

            abstract = " ".join(filter(None, abstract_text_parts)) or None

            # --- Extract Journal and Year ---
            journal_title = (
                str(journal_info.get("Title")) if journal_info.get("Title") else None
            )
            year_str: str | None = None
            if pub_date.get("Year"):
                year_str = str(pub_date.get("Year"))
            elif pub_date.get("MedlineDate"):  # Fallback
                medline_date_str = str(pub_date.get("MedlineDate"))
                match = re.search(r"^(\d{4})", medline_date_str)
                if match:
                    year_str = match.group(1)

            year = int(year_str) if year_str and year_str.isdigit() else None

            # --- Extract Authors ---
            authors_list = article_info.get("AuthorList", [])
            authors = []
            if isinstance(authors_list, list):
                for author_entry in authors_list:
                    # Check if author_entry is dict-like
                    if isinstance(author_entry, dict):
                        last_name = author_entry.get("LastName", "")
                        fore_name = author_entry.get("ForeName", "")
                        initials = author_entry.get("Initials", "")
                        # Construct name: Prefer ForeName + LastName, fallback to LastName + Initials, etc.
                        name_parts = [
                            n
                            for n in [fore_name, last_name]
                            if n and isinstance(n, str)
                        ]
                        if name_parts:
                            authors.append(" ".join(name_parts))
                        elif (
                            last_name
                            and initials
                            and isinstance(last_name, str)
                            and isinstance(initials, str)
                        ):
                            authors.append(f"{last_name} {initials}")
                        elif initials and isinstance(
                            initials, str
                        ):  # Fallback if only initials are present
                            authors.append(initials)
                        # Handle CollectiveName if present and no individual authors
                        collective_name = author_entry.get("CollectiveName")
                        if (
                            collective_name
                            and isinstance(collective_name, str)
                            and not authors
                        ):
                            authors.append(collective_name)
                    else:
                        logger.warning(
                            f"Skipping non-dict item in AuthorList: {author_entry}"
                        )
            authors = authors or None

            # --- Extract Keywords ---
            keyword_list_data = medline.get("KeywordList", [])
            keywords = []
            if isinstance(keyword_list_data, list):
                for kw_group in keyword_list_data:
                    if hasattr(kw_group, "__iter__") and not isinstance(
                        kw_group, (str, bytes, dict)
                    ):
                        # Handles [['keyword1'], ['keyword2']] structure sometimes seen
                        for inner_kw_list in kw_group:
                            if hasattr(inner_kw_list, "__iter__") and not isinstance(
                                inner_kw_list, (str, bytes)
                            ):
                                keywords.extend(str(kw) for kw in inner_kw_list if kw)
                            elif (
                                inner_kw_list
                            ):  # Handle plain string keyword in inner list
                                keywords.append(str(inner_kw_list))
                    elif isinstance(
                        kw_group, dict
                    ):  # Handles {'Keyword': ['kw1', 'kw2']} or {'Keyword': 'kw1'} etc.
                        kws = kw_group.get("Keyword", [])
                        if isinstance(kws, list):
                            keywords.extend(str(kw) for kw in kws if kw)
                        elif kws:
                            keywords.append(str(kws))
                    elif kw_group:  # Handle plain string keyword in top list
                        keywords.append(str(kw_group))

            keywords = list(set(keywords))  # Deduplicate
            keywords = keywords or None

            # --- Create SearchResult ---
            search_result = models.SearchResult(
                review_id=review_id,
                source_db=SearchDatabaseSource.PUBMED,
                source_id=pmid,
                doi=doi or None,
                title=title,
                abstract=abstract,
                journal=journal_title,
                year=year,
                authors=authors,
                keywords=keywords,
                raw_data=api_record,  # Store the original record
                source_metadata={  # Store potentially useful extra IDs/info
                    "pmc": pmc or None,
                    "publication_status": pubmed_data.get("PublicationStatus", None),
                    "mesh_headings": medline.get(
                        "MeshHeadingList", None
                    ),  # Example: include MeSH
                },
            )
            return search_result

        except Exception as e:
            pmid_for_log = pmid if pmid else "UNKNOWN"
            logger.opt(exception=True).error(
                f"Error mapping PubMed record {pmid_for_log}: {e}. Record snippet: {str(api_record)[:500]}"
            )
            return None  # Return None on any mapping error

    def _map_scopus_to_search_result(
        self, review_id: uuid.UUID, api_record: dict
    ) -> models.SearchResult | None:
        """Maps a raw Scopus API record (dict) to a SearchResult model."""
        source_id: str = ""
        try:
            # --- Extract Core IDs ---
            scopus_id_full = api_record.get(
                "dc:identifier", ""
            )  # e.g., "SCOPUS_ID:85100353187"
            if isinstance(scopus_id_full, str) and scopus_id_full.startswith(
                "SCOPUS_ID:"
            ):
                source_id = scopus_id_full.split(":")[1]
            else:
                # Fallback or alternative ID field if structure differs
                source_id = api_record.get("eid", "")  # Another common identifier

            if not source_id:
                logger.warning(
                    "Skipping Scopus record: Missing dc:identifier (SCOPUS_ID) or eid"
                )
                return None

            doi = api_record.get("prism:doi")

            # --- Extract Title ---
            title = api_record.get("dc:title")
            if not title:
                logger.warning(f"Skipping Scopus record {source_id}: Missing dc:title")
                return None

            # --- Extract Abstract ---
            # Note: Abstract might require specific API view ('STANDARD' vs 'COMPLETE')
            abstract = api_record.get("dc:description")

            # --- Extract Journal and Year ---
            journal_title = api_record.get("prism:publicationName")
            cover_date = api_record.get("prism:coverDate")  # Often YYYY-MM-DD
            year_value = cover_date.split("-")[0] if cover_date else None
            # Convert to string if it's an integer or other type
            if year_value is not None and not isinstance(year_value, str):
                year_value = str(year_value)

            # --- Extract Authors ---
            # Authors often in an 'author' list of dicts
            authors_list = api_record.get("author", [])
            authors = []
            if isinstance(authors_list, list):
                for author_entry in authors_list:
                    if isinstance(author_entry, dict):
                        # Scopus format varies, prefer 'authname', fallback '$'
                        name = author_entry.get("authname", author_entry.get("$", None))
                        if name and isinstance(name, str):
                            authors.append(name)
            authors = authors or None

            # --- Extract Keywords ---
            # Keywords often in 'authkeywords' as a |-separated string
            keywords_str = api_record.get("authkeywords")
            keywords = None
            if isinstance(keywords_str, str) and keywords_str:
                keywords = [kw.strip() for kw in keywords_str.split("|") if kw.strip()]
                keywords = keywords or None  # Ensure None if list becomes empty

            # --- Create SearchResult ---
            search_result = models.SearchResult(
                review_id=review_id,
                source_db=SearchDatabaseSource.SCOPUS,
                source_id=source_id,
                doi=doi,
                title=title,
                abstract=abstract,
                journal=journal_title,
                year=year_value,
                authors=authors,
                keywords=keywords,
                raw_data=api_record,  # Store the original record
                source_metadata={  # Store other potentially useful fields
                    "eid": api_record.get("eid"),
                    "subtype": api_record.get("subtype"),
                    "subtype_description": api_record.get("subtypeDescription"),
                    "aggregation_type": api_record.get("prism:aggregationType"),
                    # Add other relevant Scopus-specific fields here
                },
            )
            return search_result

        except Exception as e:
            sid_for_log = source_id if source_id else "UNKNOWN"
            logger.opt(exception=True).error(
                f"Error mapping Scopus record {sid_for_log}: {e}. Record snippet: {str(api_record)[:500]}"
            )
            return None  # Return None on any mapping error

    # --- Core Service Methods (Synchronous) ---
    def add_api_search_results(
        self,
        review_id: uuid.UUID,
        source_db: SearchDatabaseSource,
        api_records: list[dict[str, t.Any]],
        session: Session | None = None,
    ) -> Sequence[models.SearchResult]:
        """Adds search results from an external API (sync).

        Args:
            review_id: The ID of the review to add results to.
            source_db: The source database of the results.
            api_records: A list of raw API record dictionaries.
            session: An optional existing DB session to use.

        Returns:
            A sequence of the added SearchResult objects (refreshed).
        """
        added_results_output: list[models.SearchResult] = []
        if not api_records:
            return added_results_output

        logger.info(
            f"Processing {len(api_records)} results from {source_db.name} for review {review_id}"
        )

        # Use the helper context manager
        with self._get_session(session) as active_session:
            try:
                results_to_add: list[models.SearchResult] = []
                for record in api_records:
                    mapped_result: models.SearchResult | None = None
                    if source_db == SearchDatabaseSource.PUBMED:
                        mapped_result = self._map_pubmed_to_search_result(
                            review_id, record
                        )
                    elif source_db == SearchDatabaseSource.SCOPUS:
                        mapped_result = self._map_scopus_to_search_result(
                            review_id, record
                        )
                    else:
                        logger.warning(f"Unsupported source database: {source_db.name}")
                        continue

                    if mapped_result:
                        if not mapped_result.source_id or not mapped_result.title:
                            logger.warning(
                                f"Skipping record due to missing source_id or title: {mapped_result}"
                            )
                            continue
                        results_to_add.append(mapped_result)
                    else:
                        logger.warning(f"Failed to map record from {source_db.name}")

                if not results_to_add:
                    logger.info("No valid results mapped for addition.")
                    return []

                try:
                    # Pass the active_session to the repository method
                    added_results_in_session = self.search_repo.add_all(
                        active_session, results_to_add
                    )
                    logger.info(
                        f"Added/flushed {len(results_to_add)} results to session."
                    )
                except repositories.ConstraintViolationError as e:
                    logger.warning(
                        f"Constraint violation during add_all, likely duplicates: {e}"
                    )
                    # Always rollback the active session state after a constraint violation
                    # regardless of whether it was passed in or created internally.
                    active_session.rollback()
                    logger.info(
                        "Rolled back due to constraint violation. Returning empty list."
                    )
                    return []  # Return empty list as per current logic

                # Only commit if we created the session
                if not session:
                    active_session.commit()
                logger.info(
                    f"Committed {len(added_results_in_session)} search results for review {review_id}"
                )

                # Refresh results after commit/flush using the active session
                refreshed_results = []
                for (
                    result
                ) in added_results_in_session:  # Use the result from repo.add_all
                    try:
                        active_session.refresh(result)
                        refreshed_results.append(result)
                    except Exception as refresh_exc:
                        logger.warning(
                            f"Could not refresh result {getattr(result, 'id', '?')}: {refresh_exc}"
                        )
                        # Still append the potentially stale object if refresh fails
                        refreshed_results.append(result)
                added_results_output = refreshed_results

            except Exception as e:
                logger.exception(
                    f"Error adding API search results for review {review_id}: {e}"
                )
                # Always rollback the active session state after an error
                active_session.rollback()
                raise ServiceError(f"Failed to add search results: {e}") from e

        return added_results_output

    def get_search_results_by_review_id(
        self, review_id: uuid.UUID, session: Session | None = None
    ) -> Sequence[models.SearchResult]:
        """Retrieves all search results associated with a given review ID (sync)."""
        logger.debug(f"Getting search results for review_id: {review_id}")
        # Use the helper context manager
        with self._get_session(session) as active_session:
            try:
                # Pass the active_session to the repository method
                results = self.search_repo.get_by_review_id(active_session, review_id)
                logger.debug(
                    f"Found {len(results)} search results for review {review_id}"
                )
                return results
            except Exception as e:
                logger.exception(
                    f"Error getting search results for review {review_id}: {e}"
                )
                raise ServiceError(f"Failed to get search results: {e}") from e

    def get_search_result_by_source_details(
        self,
        review_id: uuid.UUID,
        source_db: SearchDatabaseSource,
        source_id: str,
        session: Session | None = None,
    ) -> models.SearchResult | None:
        """Retrieves a specific search result using its source database and ID (sync)."""
        logger.debug(
            f"Getting search result for review {review_id}, source {source_db.name}, id {source_id}"
        )
        # Use the helper context manager
        with self._get_session(session) as active_session:
            try:
                # Pass the active_session to the repository method
                result = self.search_repo.get_by_source_details(
                    active_session, review_id, source_db, source_id
                )
                if result:
                    logger.debug(f"Found search result {result.id}")
                else:
                    logger.debug("Search result not found.")
                return result
            except Exception as e:
                # Fix implicit string concatenation in f-string
                log_message = (
                    f"Error getting search result for review {review_id}, "
                    f"source {source_db.name}, id {source_id}: {e}"
                )
                logger.exception(log_message)
                raise ServiceError(
                    f"Failed to get search result by source details: {e}"
                ) from e

    # Implement sync update method
    def update_search_result(
        self, result_update: models.SearchResult, session: Session | None = None
    ) -> models.SearchResult:
        """Updates an existing search result."""
        if not result_update.id:
            raise ValueError("Cannot update SearchResult without an ID.")
        logger.debug(f"Updating search result id: {result_update.id}")
        # Use the helper context manager
        with self._get_session(session) as active_session:
            try:
                # Pass the active_session to the repository method
                updated_result = self.search_repo.update(active_session, result_update)
                # Only commit if we created the session
                if not session:
                    active_session.commit()
                # Refresh using the active session
                active_session.refresh(updated_result)
                logger.info(f"Successfully updated search result {updated_result.id}")
                return updated_result
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Update failed: {e}")
                # Always rollback the active session state after an error
                active_session.rollback()
                raise
            except Exception as e:
                logger.exception(
                    f"Error updating search result {result_update.id}: {e}"
                )
                # Always rollback the active session state after an error
                active_session.rollback()
                raise ServiceError(f"Failed to update search result: {e}") from e

    # Implement sync delete method
    def delete_search_result(
        self, result_id: uuid.UUID, session: Session | None = None
    ):
        """Deletes a search result by its ID."""
        logger.debug(f"Deleting search result id: {result_id}")
        # Use the helper context manager
        with self._get_session(session) as active_session:
            try:
                # Pass the active_session to the repository method
                self.search_repo.delete(active_session, result_id)
                # Only commit if we created the session
                if not session:
                    active_session.commit()
                logger.info(f"Successfully deleted search result {result_id}")
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Delete failed: {e}")
                # Always rollback the active session state after an error
                active_session.rollback()
            except Exception as e:
                logger.exception(f"Error deleting search result {result_id}: {e}")
                # Always rollback the active session state after an error
                active_session.rollback()
                raise ServiceError(f"Failed to delete search result: {e}") from e


# --- Review Service ---


class ReviewService(BaseService):
    """Service for managing SystematicReview data (sync)."""

    def __init__(
        self,
        factory: sessionmaker[Session] = session_factory,
        review_repo: repositories.SystematicReviewRepository | None = None,
    ):
        super().__init__(factory)
        self.review_repo = review_repo or repositories.SystematicReviewRepository()

    def create_review(
        self, review_data: schemas.SystematicReviewCreate
    ) -> models.SystematicReview:
        """Creates a new systematic review."""
        logger.debug(f"Creating review: {review_data.research_question[:50]}...")
        # Create the SQLModel instance from the Pydantic schema
        # Ensure all required fields are present in SystematicReviewCreate or have defaults
        review = models.SystematicReview.model_validate(review_data)
        with self.session_factory() as session:
            try:
                added_review = self.review_repo.add(session, review)
                session.commit()
                session.refresh(added_review)
                logger.info(f"Successfully created review {added_review.id}")
                return added_review
            except Exception as e:
                logger.exception("Error creating review")
                session.rollback()
                raise ServiceError(f"Failed to create review: {e}") from e

    def get_review(self, review_id: uuid.UUID) -> models.SystematicReview | None:
        """Gets a review by its ID."""
        logger.debug(f"Getting review id: {review_id}")
        with self.session_factory() as session:
            try:
                review = self.review_repo.get_by_id(session, review_id)
                if not review:
                    logger.warning(f"Review {review_id} not found.")
                    return None
                logger.debug(f"Found review {review_id}")
                return review
            except Exception as e:
                logger.exception(f"Error getting review {review_id}: {e}")
                raise ServiceError(f"Failed to get review: {e}") from e

    def get_all_reviews(self) -> Sequence[models.SystematicReview]:
        """Gets all reviews."""
        logger.debug("Getting all reviews")
        with self.session_factory() as session:
            try:
                reviews = self.review_repo.get_all(session)
                logger.debug(f"Found {len(reviews)} reviews")
                return reviews
            except Exception as e:
                logger.exception("Error getting all reviews")
                raise ServiceError(f"Failed to get all reviews: {e}") from e

    def update_review(
        self, review_update_data: schemas.SystematicReviewUpdate, review_id: uuid.UUID
    ) -> models.SystematicReview:
        """Updates an existing review using partial update schema."""
        logger.debug(f"Updating review id: {review_id}")
        with self.session_factory() as session:
            try:
                # Get the existing review
                db_review = self.review_repo.get_by_id(session, review_id)
                if not db_review:
                    raise repositories.RecordNotFoundError(
                        f"Review {review_id} not found for update."
                    )

                # Create dictionary of fields to update, excluding unset fields from input
                update_data = review_update_data.model_dump(exclude_unset=True)

                # Update the db_review object with new data
                db_review.sqlmodel_update(update_data)

                # Add to session (marks it as dirty) and commit
                session.add(db_review)
                session.commit()
                session.refresh(db_review)
                logger.info(f"Successfully updated review {db_review.id}")
                return db_review
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Update failed: {e}")
                session.rollback()
                raise  # Re-raise specific error
            except Exception as e:
                logger.exception(f"Error updating review {review_id}: {e}")
                session.rollback()
                raise ServiceError(f"Failed to update review: {e}") from e

    def delete_review(self, review_id: uuid.UUID) -> None:
        """Deletes a review. Assumes DB cascade handles related data or it's handled separately."""
        logger.debug(f"Deleting review id: {review_id}")
        with self.session_factory() as session:
            try:
                self.review_repo.delete(session, review_id)
                session.commit()
                logger.info(f"Successfully deleted review {review_id}")
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Delete failed: {e}")
                session.rollback()
            except Exception as e:
                logger.exception(f"Error deleting review {review_id}: {e}")
                session.rollback()
                raise ServiceError(f"Failed to delete review: {e}") from e


# --- Screening Service ---


class ScreeningService(BaseService):
    """Service for managing screening results and resolutions (sync)."""

    def __init__(
        self,
        factory: sessionmaker[Session] = session_factory,
        screen_repo: repositories.ScreenAbstractResultRepository | None = None,
        resolution_repo: repositories.ScreeningResolutionRepository | None = None,
        search_repo: repositories.SearchResultRepository | None = None,
    ):
        super().__init__(factory)
        self.screen_repo = screen_repo or repositories.ScreenAbstractResultRepository()
        self.resolution_repo = (
            resolution_repo or repositories.ScreeningResolutionRepository()
        )
        self.search_repo = search_repo or repositories.SearchResultRepository()

    def add_screening_result(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
        screening_response: schemas.ScreeningResponse,
        run_id: uuid.UUID | None = None,
        trace_id: uuid.UUID | None = None,
        start_time: datetime | None = None,
        end_time: datetime | None = None,
        model_name: str | None = "unknown",
        response_metadata: dict | None = None,
    ) -> models.ScreenAbstractResult:
        logger.debug(f"Adding {strategy.name} screening result for review {review_id}")
        result_id = run_id or uuid.uuid4()
        result_model = models.ScreenAbstractResult(
            id=result_id,
            review_id=review_id,
            decision=screening_response.decision,
            confidence_score=screening_response.confidence_score,
            rationale=screening_response.rationale,
            extracted_quotes=screening_response.extracted_quotes,
            exclusion_reason_categories=(
                screening_response.exclusion_reason_categories.model_dump()
                if screening_response.exclusion_reason_categories
                else {}
            ),
            screening_strategy=strategy,
            trace_id=trace_id,
            start_time=start_time,
            end_time=end_time,
            model_name=model_name or "unknown",
            response_metadata=response_metadata or {},
        )
        with self.session_factory() as session:
            try:
                added_result = self.screen_repo.add(session, result_model)
                session.commit()
                try:
                    session.refresh(added_result)
                except Exception:
                    logger.warning(f"Could not refresh screening result {result_id}")
                logger.info(f"Added screening result {added_result.id}")
                return added_result
            except Exception as e:
                logger.exception("Error adding screening result")
                session.rollback()
                raise ServiceError(f"Failed to add screening result: {e}") from e

    # --- New Methods ---
    def get_screening_result_by_strategy(
        self,
        session: Session,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
    ) -> models.ScreenAbstractResult | None:
        """Gets a screening result for a review and strategy. Returns first found."""
        try:
            results = self.screen_repo.get_by_review_strategy(
                session, review_id, strategy
            )
            if results:
                if len(results) > 1:
                    logger.warning(
                        f"Found multiple screening results for review {review_id} and strategy {strategy.name}. Returning first."
                    )
                return results[0]
            return None
        except Exception as e:
            logger.exception("Error getting specific screening result by strategy")
            return None

    def add_or_update_screening_result(
        self,
        review_id: uuid.UUID,
        strategy: ScreeningStrategyType,
        screening_response: schemas.ScreeningResponse,
        run_id: uuid.UUID | None = None,
        trace_id: uuid.UUID | None = None,
        start_time: datetime | None = None,
        end_time: datetime | None = None,
        model_name: str | None = "unknown",
        response_metadata: dict | None = None,
    ) -> models.ScreenAbstractResult:
        """Adds a new screening result or updates it based on run_id or review/strategy match."""
        logger.debug(
            f"Adding/Updating {strategy.name} screening result for review {review_id}"
        )

        with self.session_factory() as session:
            try:
                existing_result = None
                if run_id:
                    # Try finding by primary key first if provided
                    existing_result = self.screen_repo.get_by_id(session, run_id)

                if not existing_result:
                    # If no run_id or not found by run_id, try finding by review/strategy
                    # This assumes strategy should be unique per review (might need adjustment)
                    existing_result = self.get_screening_result_by_strategy(
                        session, review_id, strategy
                    )

                result_id = run_id or (
                    existing_result.id if existing_result else uuid.uuid4()
                )

                update_data = {
                    "review_id": review_id,
                    "decision": screening_response.decision,
                    "confidence_score": screening_response.confidence_score,
                    "rationale": screening_response.rationale,
                    "extracted_quotes": screening_response.extracted_quotes,
                    "exclusion_reason_categories": (
                        screening_response.exclusion_reason_categories.model_dump()
                        if screening_response.exclusion_reason_categories
                        else {}
                    ),
                    "screening_strategy": strategy,
                    "trace_id": trace_id,
                    "start_time": start_time,
                    "end_time": end_time,
                    "model_name": model_name or "unknown",
                    "response_metadata": response_metadata or {},
                }

                if existing_result:
                    logger.info(
                        f"Updating existing screening result {existing_result.id}"
                    )
                    updated_result_model = existing_result.sqlmodel_update(update_data)
                    # Pass the specific model instance to update
                    updated_result = self.screen_repo.update(
                        session, updated_result_model
                    )
                else:
                    logger.info(f"Creating new screening result with id {result_id}")
                    # Create new model instance including the ID
                    new_result_model = models.ScreenAbstractResult(
                        id=result_id, **update_data
                    )
                    updated_result = self.screen_repo.add(session, new_result_model)

                session.commit()
                session.refresh(updated_result)
                logger.info(
                    f"Successfully added/updated screening result {updated_result.id}"
                )
                return updated_result

            except Exception as e:
                logger.exception(
                    f"Error adding/updating screening result for review {review_id}"
                )
                session.rollback()
                raise ServiceError(f"Failed to add/update screening result: {e}") from e

    def get_screening_results_for_review(
        self, review_id: uuid.UUID
    ) -> Sequence[models.ScreenAbstractResult]:
        """Gets all screening results for a specific review."""
        logger.debug(f"Getting all screening results for review {review_id}")
        with self.session_factory() as session:
            try:
                results = self.screen_repo.get_by_review_id(session, review_id)
                logger.debug(f"Found {len(results)} screening results.")
                return results
            except Exception as e:
                logger.exception(
                    f"Error getting screening results for review {review_id}"
                )
                raise ServiceError(f"Failed to get screening results: {e}") from e

    # TODO: get_conflicting_results(...)
    # TODO: add_resolution(...)


# TODO: Define other services (ReviewService, ScreeningService, LogService) following the same pattern.
# Example:
# class ReviewService(BaseService):
#     def __init__(self, ...):
#         self.review_repo = repositories.SystematicReviewRepository()

#     async def create_review(self, review_data: schemas.SystematicReviewCreate) -> models.SystematicReview:
#         async with self.session_factory() as session:
#             try:
#                 review = models.SystematicReview.model_validate(review_data)
#                 added_review = await self.review_repo.add(session, review)
#                 await session.commit()
#                 await session.refresh(added_review)
#                 return added_review
#             except Exception as e:
#                 await session.rollback()
#                 raise ServiceError(f"Failed to create review: {e}") from e
