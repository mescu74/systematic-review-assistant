# pyright: reportPrivateUsage=false
"""Application Service Layer - Synchronous Implementation."""

from __future__ import annotations

import os  # Import os for getenv
import re
import typing as t
import uuid
from collections.abc import Mapping
from datetime import datetime

# Import BioPython Entrez for PubMed API interaction
# Assuming BioPython is installed and configured (email, api_key)
from Bio import Entrez
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

    # --- PubMed Data Cleaning and Parsing Helpers ---
    def _recursive_clean(self, data: t.Any) -> t.Any:
        """Recursively convert BioPython Entrez parser elements into plain Python types."""
        if data is None:  # Explicitly handle None first
            return None
        if isinstance(data, list):
            return [self._recursive_clean(item) for item in data]
        if isinstance(data, dict):
            return {k: self._recursive_clean(v) for k, v in data.items()}
        # Handle BioPython's StringElement and other similar elements that might have attributes
        # or a special _value attribute, but can often be treated as strings directly.
        # The key is to get to the actual string, int, etc. value.
        if hasattr(data, "__str__") and not isinstance(data, (str, int, float, bool)):
            # This is a broad catch; specific Bio.Entrez.Element types might need tailored handling
            # if str(data) isn't sufficient or if attributes need to be preserved differently.
            # For now, assume str() is a reasonable default for unknown Entrez elements.
            if hasattr(
                data, "_value"
            ):  # Some BioPython elements store main content in _value
                return data._value  # noqa: SLF001
            return str(data)  # Fallback to string representation
        return data

    def _extract_text_from_element(self, element: t.Any, default: str = "") -> str:
        """Safely extracts text, handling None or BioPython elements."""
        if element is None:
            return default
        if hasattr(element, "_value"):  # Common for Bio.Entrez.Element.StringElement
            return str(element._value)  # noqa: SLF001
        return str(element)

    def _parse_pubmed_ids(
        self, cleaned_article_data: dict[str, t.Any]
    ) -> tuple[str | None, str | None, str | None]:
        """Extracts PMID, DOI, and PMC from cleaned PubMed data."""
        pmid = self._extract_text_from_element(
            cleaned_article_data.get("MedlineCitation", {}).get("PMID")
        )

        doi: str | None = None
        pmc: str | None = None

        # DOI and PMC can be in PubmedData/ArticleIdList or Article/ELocationID
        pubmed_data_ids = cleaned_article_data.get("PubmedData", {}).get(
            "ArticleIdList", []
        )
        if isinstance(pubmed_data_ids, list):
            for item in pubmed_data_ids:
                if isinstance(item, dict):
                    id_type = self._extract_text_from_element(
                        item.get("IdType")
                    ).lower()
                    id_value = self._extract_text_from_element(
                        item.get("$")
                    )  # Entrez often puts value in '$'
                    if not id_value:  # Sometimes it's in '#text' or direct from element after cleaning
                        id_value = self._extract_text_from_element(item)

                    if id_type == "doi":
                        doi = id_value
                    elif id_type == "pmc":
                        pmc = id_value

        if not doi:
            article_elocation_ids = (
                cleaned_article_data.get("MedlineCitation", {})
                .get("Article", {})
                .get("ELocationID", [])
            )
            if isinstance(article_elocation_ids, list):
                for item in article_elocation_ids:
                    if isinstance(item, dict):
                        id_type = self._extract_text_from_element(
                            item.get("EIdType")
                        ).lower()
                        if id_type == "doi":
                            doi = self._extract_text_from_element(item.get("$"))
                            if not doi:
                                doi = self._extract_text_from_element(item)
                            break
        return pmid or None, doi or None, pmc or None

    def _parse_pubmed_title_abstract(
        self, cleaned_article_data: dict[str, t.Any]
    ) -> tuple[str | None, str | None]:
        """Extracts title and abstract from cleaned PubMed data."""
        article_info = cleaned_article_data.get("MedlineCitation", {}).get(
            "Article", {}
        )
        title = self._extract_text_from_element(article_info.get("ArticleTitle"))

        abstract_text_parts = []
        abstract_section = article_info.get("Abstract", {}).get("AbstractText", [])
        if isinstance(abstract_section, list):
            for part in abstract_section:
                if isinstance(
                    part, dict
                ):  # Handle sections like { 'Label': 'BACKGROUND', '$: 'text...'}
                    label = self._extract_text_from_element(part.get("Label"))
                    text = self._extract_text_from_element(part.get("$"))
                    if not text:  # Fallback if value not in '$'
                        text = self._extract_text_from_element(
                            part
                        )  # Check if part itself is the text after cleaning

                    if (
                        label and text and label.lower() != text.lower()
                    ):  # Avoid label being the text
                        abstract_text_parts.append(f"{label.upper()}: {text}")
                    elif text:
                        abstract_text_parts.append(text)
                else:  # Plain string part
                    abstract_text_parts.append(self._extract_text_from_element(part))
        elif isinstance(abstract_section, str):  # Single string abstract
            abstract_text_parts.append(abstract_section)

        abstract = (
            " ".join(filter(None, abstract_text_parts)) if abstract_text_parts else None
        )
        return title or None, abstract or None

    def _parse_pubmed_journal_year(
        self, cleaned_article_data: dict[str, t.Any]
    ) -> tuple[str | None, str | None]:
        """Extracts journal and year from cleaned PubMed data."""
        article_info = cleaned_article_data.get("MedlineCitation", {}).get(
            "Article", {}
        )
        journal_info = article_info.get("Journal", {})
        journal_title = self._extract_text_from_element(journal_info.get("Title"))

        year_str: str | None = None
        pub_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
        if isinstance(pub_date, dict):
            year_str = self._extract_text_from_element(pub_date.get("Year"))
            if not year_str:  # Fallback for MedlineDate e.g. "2023 Spring"
                medline_date_str = self._extract_text_from_element(
                    pub_date.get("MedlineDate")
                )
                if medline_date_str:
                    match = re.search(r"^(\d{4})", medline_date_str)
                    if match:
                        year_str = match.group(1)
        return journal_title or None, year_str or None

    def _parse_pubmed_authors(
        self, cleaned_article_data: dict[str, t.Any]
    ) -> list[str] | None:
        """Extracts authors from cleaned PubMed data."""
        authors_list_data = (
            cleaned_article_data.get("MedlineCitation", {})
            .get("Article", {})
            .get("AuthorList", [])
        )
        authors: list[str] = []
        if isinstance(authors_list_data, list):
            for author_entry in authors_list_data:
                if isinstance(author_entry, dict):
                    last_name = self._extract_text_from_element(
                        author_entry.get("LastName")
                    )
                    fore_name = self._extract_text_from_element(
                        author_entry.get("ForeName")
                    )
                    initials = self._extract_text_from_element(
                        author_entry.get("Initials")
                    )
                    collective_name = self._extract_text_from_element(
                        author_entry.get("CollectiveName")
                    )

                    if fore_name and last_name:
                        authors.append(f"{fore_name} {last_name}")
                    elif last_name and initials:
                        authors.append(f"{last_name} {initials}")
                    elif last_name:
                        authors.append(last_name)
                    elif collective_name:
                        authors.append(collective_name)
                    elif (
                        initials
                    ):  # Fallback if only initials are present and no other name part
                        authors.append(initials)
        return authors if authors else None

    def _parse_pubmed_keywords(
        self, cleaned_article_data: dict[str, t.Any]
    ) -> list[str] | None:
        """Extracts keywords from cleaned PubMed data."""
        keyword_list_data = (
            cleaned_article_data.get("MedlineCitation", {})
            .get("Article", {})
            .get("KeywordList", [])
        )
        keywords: list[str] = []
        if isinstance(keyword_list_data, list):
            for (
                kw_group
            ) in keyword_list_data:  # This can be a list of lists or list of dicts
                if isinstance(kw_group, list):  # e.g., [[kw1, kw2], [kw3]]
                    for inner_kw_list in kw_group:
                        if isinstance(inner_kw_list, list):
                            keywords.extend(
                                self._extract_text_from_element(kw)
                                for kw in inner_kw_list
                                if kw
                            )
                        else:  # inner_kw_list is a direct keyword string/element
                            keywords.append(
                                self._extract_text_from_element(inner_kw_list)
                            )
                elif isinstance(kw_group, dict):  # e.g. { Keyword: ['kw1', 'kw2']}
                    actual_keywords = kw_group.get(
                        "Keyword", []
                    )  # 'Keyword' is a common key
                    if isinstance(actual_keywords, list):
                        keywords.extend(
                            self._extract_text_from_element(kw)
                            for kw in actual_keywords
                            if kw
                        )
                    else:
                        keywords.append(
                            self._extract_text_from_element(actual_keywords)
                        )
                else:  # Direct keyword string/element in the outer list
                    keywords.append(self._extract_text_from_element(kw_group))

        return (
            list({k for k in keywords if k}) if keywords else None
        )  # Deduplicate and remove empty strings

    # --- Original Mapping function, now a coordinator ---
    def _map_pubmed_to_search_result(
        self,
        review_id: uuid.UUID,
        api_record: dict[str, t.Any],  # Will receive a CLEANED dict
    ) -> models.SearchResult | None:
        """Maps a CLEANED PubMed API record dictionary to a SearchResult model."""
        try:
            # The api_record is already cleaned by _recursive_clean
            pmid, doi, pmc = self._parse_pubmed_ids(api_record)
            title, abstract = self._parse_pubmed_title_abstract(api_record)
            journal, year = self._parse_pubmed_journal_year(api_record)
            authors = self._parse_pubmed_authors(api_record)
            keywords = self._parse_pubmed_keywords(api_record)

            if not pmid or not title:
                logger.warning(
                    f"Skipping record due to missing PMID or title after parsing. PMID: {pmid}, Title: {title}."
                )
                return None

            # Construct SearchResult
            search_result = models.SearchResult(
                review_id=review_id,
                source_db=SearchDatabaseSource.PUBMED,
                source_id=pmid,
                doi=doi,
                title=title,
                abstract=abstract,
                journal=journal,
                year=year,
                authors=authors,
                keywords=keywords,
                raw_data=api_record,
                source_metadata={
                    "pmc": pmc,
                    "publication_status": self._extract_text_from_element(
                        api_record.get("PubmedData", {}).get("PublicationStatus")
                    ),
                    # Add MeSH Headings back
                    "mesh_headings": api_record.get("MedlineCitation", {}).get(
                        "MeshHeadingList", []
                    ),
                },
            )
            return search_result

        except Exception as e:
            pmid_for_log = api_record.get("MedlineCitation", {}).get("PMID", "UNKNOWN")
            # Ensure pmid_for_log is a string for the log message
            if not isinstance(pmid_for_log, str):
                pmid_for_log = str(pmid_for_log)

            logger.opt(exception=True).error(
                f"Error mapping cleaned PubMed record {pmid_for_log}: {e!r}. Record snippet: {str(api_record)[:500]}"
            )
            return None

    def _map_scopus_to_search_result(
        self,
        review_id: uuid.UUID,
        api_record: dict[str, t.Any],  # Updated type hint
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
    def search_pubmed_and_store_results(
        self,
        review_id: uuid.UUID,
        query: str,
        max_results: int = 100,
    ) -> Sequence[schemas.SearchResultRead]:
        """Performs a PubMed search, maps results, stores them, and handles sessions.

        Args:
            review_id: The ID of the review to associate results with.
            query: The search query string for PubMed.
            max_results: The maximum number of results to fetch and store.

        Returns:
            A sequence of the added/stored SearchResult objects, converted to SearchResultRead schemas.

        Raises:
            ServiceError: If there is an issue with the API call or database operation,
                        or if NCBI credentials are not set in environment variables.
        """
        logger.info(
            f"Starting PubMed search for review {review_id!r} with query {query!r} (max: {max_results})"
        )

        entrez_email = os.getenv("NCBI_EMAIL")
        entrez_api_key = os.getenv("NCBI_API_KEY")

        if not entrez_email:
            logger.error(
                "NCBI_EMAIL environment variable not set. PubMed search cannot proceed."
            )
            raise ServiceError("NCBI_EMAIL environment variable not set.")

        Entrez.email = entrez_email
        if entrez_api_key:
            Entrez.api_key = entrez_api_key
        else:
            logger.warning("NCBI_API_KEY not set. PubMed searches may be rate-limited.")

        fetched_pmids: list[str] = []
        new_pmids_to_fetch_details: list[str] = []
        mapped_model_results: list[models.SearchResult] = []
        added_search_result_models_for_conversion: list[models.SearchResult] = []

        try:
            logger.debug("Executing Entrez.esearch for PMIDs...")
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=max_results, sort="relevance"
            )
            search_results_raw = Entrez.read(handle)
            handle.close()

            if isinstance(search_results_raw, Mapping):
                fetched_pmids = search_results_raw.get("IdList", [])
            else:
                logger.error(
                    f"Entrez.esearch returned unexpected type: {type(search_results_raw)}"
                )

            if not fetched_pmids:
                logger.info("No PMIDs found by Entrez.esearch for query.")
                return []

            logger.debug(
                f"Found {len(fetched_pmids)} PMIDs from PubMed initial search."
            )

            # Check for existing PMIDs in the database for this review
            with (
                self.session_factory() as session
            ):  # Use a new session for this read operation
                existing_pmids_in_db = self.search_repo.get_existing_source_ids(
                    session,
                    review_id,
                    SearchDatabaseSource.PUBMED,
                    fetched_pmids,
                )

            new_pmids_to_fetch_details = [
                pmid for pmid in fetched_pmids if pmid not in existing_pmids_in_db
            ]

            if not new_pmids_to_fetch_details:
                logger.info(
                    f"All {len(fetched_pmids)} PMIDs found already exist in the database for review {review_id!r}. No new articles to fetch."
                )
                # To provide feedback on existing articles, one might consider fetching and returning them here,
                # but current story implies returning newly added ones. For now, return empty.
                return []

            logger.debug(
                f"{len(new_pmids_to_fetch_details)} new PMIDs to fetch details for with Entrez.efetch..."
            )
            handle = Entrez.efetch(
                db="pubmed", id=new_pmids_to_fetch_details, rettype="xml", retmode="xml"
            )
            papers_raw_data = Entrez.read(handle)
            handle.close()

            cleaned_papers_input = self._recursive_clean(papers_raw_data)
            actual_articles_data: list[dict[str, t.Any]] = []

            if isinstance(cleaned_papers_input, list):
                actual_articles_data = [
                    item for item in cleaned_papers_input if isinstance(item, dict)
                ]
            elif isinstance(cleaned_papers_input, dict):
                if "PubmedArticle" in cleaned_papers_input and isinstance(
                    cleaned_papers_input["PubmedArticle"], list
                ):
                    actual_articles_data = [
                        item
                        for item in cleaned_papers_input["PubmedArticle"]
                        if isinstance(item, dict)
                    ]
                elif "PubmedArticleSet" in cleaned_papers_input and isinstance(
                    cleaned_papers_input["PubmedArticleSet"], list
                ):
                    actual_articles_data = [
                        item
                        for item in cleaned_papers_input["PubmedArticleSet"]
                        if isinstance(item, dict)
                    ]
                else:
                    actual_articles_data = [cleaned_papers_input]
            else:
                logger.warning(
                    f"Unexpected data structure after cleaning PubMed efetch results: {type(cleaned_papers_input)}"
                )

            if not actual_articles_data:
                logger.info(
                    "No article data extracted after cleaning PubMed efetch results for new PMIDs."
                )
                return []

            for article_dict in actual_articles_data:
                mapped_model = self._map_pubmed_to_search_result(
                    review_id, article_dict
                )
                if mapped_model:
                    mapped_model_results.append(mapped_model)

            if not mapped_model_results:
                logger.info(
                    f"No new results successfully mapped from PubMed for review {review_id!r}."
                )
                return []

        except Exception as e:
            logger.opt(exception=True).error(
                f"Error during PubMed API interaction or mapping for query '{query}': {e!r}"
            )
            raise ServiceError(
                f"PubMed API interaction or mapping failed for query '{query}'."
            ) from e

        # This block handles database storage and transaction for NEW articles
        if (
            not mapped_model_results
        ):  # Should be redundant if logic above is correct, but safeguard
            logger.info("No new mapped models to store.")
            return []

        try:
            with self.session_factory.begin() as session:
                persisted_models = self.search_repo.add_all(
                    session,
                    mapped_model_results,  # Only try to add new, mapped models
                )

                refreshed_models: list[models.SearchResult] = []
                for model_instance in persisted_models:
                    session.refresh(model_instance)
                    refreshed_models.append(model_instance)

                added_search_result_models_for_conversion = refreshed_models

                logger.info(
                    f"Successfully stored and refreshed {len(added_search_result_models_for_conversion)} new results from PubMed for review {review_id!r}"
                )

        except (
            repositories.ConstraintViolationError
        ) as e:  # Should be less likely with proactive check, but good to keep
            logger.warning(
                f"Constraint violation storing new PubMed results for review {review_id!r}, this is unexpected after proactive check: {e!r}"
            )
            # Depending on desired behavior, could try to return already existing ones or just fail.
            # For now, returning empty, signaling no *new* articles were successfully added and persisted in this transaction.
            return []
        except Exception as e:
            logger.opt(exception=True).error(
                f"Database error storing new PubMed results for review {review_id!r}: {e!r}"
            )
            raise ServiceError(
                f"Failed to store new PubMed results for review {review_id!r}."
            ) from e

        # Final conversion to schema and return
        # return [
        #     schemas.SearchResultRead.model_validate(model_res, from_attributes=True)
        #     for model_res in added_search_result_models_for_conversion
        # ]
        # Explicitly select fields for SearchResultRead to avoid validation errors with extra model attributes
        search_result_read_fields = schemas.SearchResultRead.model_fields.keys()
        return [
            schemas.SearchResultRead.model_validate(
                {
                    field: getattr(model_res, field)
                    for field in search_result_read_fields
                    if hasattr(model_res, field)
                }
            )
            for model_res in added_search_result_models_for_conversion
        ]

    def get_search_results_by_review_id(
        self, review_id: uuid.UUID
    ) -> Sequence[schemas.SearchResultRead]:  # MODIFIED return type
        """Retrieves all search results associated with a given review ID (sync), converted to schemas."""  # MODIFIED comment
        logger.debug(f"Getting search results for review_id: {review_id!r}")
        with self.session_factory.begin() as session:
            try:
                results_models = self.search_repo.get_by_review_id(session, review_id)
                logger.debug(
                    f"Found {len(results_models)} search results for review {review_id!r}"
                )
                search_result_read_fields = schemas.SearchResultRead.model_fields.keys()
                return [
                    schemas.SearchResultRead.model_validate(
                        {
                            field: getattr(res, field)
                            for field in search_result_read_fields
                            if hasattr(res, field)
                        }
                    )
                    for res in results_models
                ]
            except Exception as e_ex:
                logger.exception(
                    f"Error getting search results for review {review_id!r}"
                )
                raise ServiceError("Failed to get search results") from e_ex

    def get_search_result_by_source_details(
        self,
        review_id: uuid.UUID,
        source_db: SearchDatabaseSource,
        source_id: str,
    ) -> schemas.SearchResultRead | None:  # MODIFIED return type
        """Retrieves a specific search result using its source database and ID (sync), converted to schema."""  # MODIFIED comment
        logger.debug(
            f"Getting search result for review {review_id!r}, source {source_db.name!r}, id {source_id!r}"
        )
        with self.session_factory.begin() as session:
            try:
                result_model = self.search_repo.get_by_source_details(
                    session, review_id, source_db, source_id
                )
                if result_model:
                    logger.debug(f"Found search result {result_model.id!r}")
                    search_result_read_fields = (
                        schemas.SearchResultRead.model_fields.keys()
                    )
                    return schemas.SearchResultRead.model_validate(
                        {
                            field: getattr(result_model, field)
                            for field in search_result_read_fields
                            if hasattr(result_model, field)
                        }
                    )
                logger.debug("Search result not found by source details.")
                return None
            except Exception as e_ex:
                logger.exception(
                    f"Error getting search result for review {review_id!r}, source {source_db.name!r}, id {source_id!r}"
                )
                raise ServiceError(
                    "Failed to get search result by source details"
                ) from e_ex

    # Implement sync update method
    def update_search_result(
        self,
        result_id: uuid.UUID,
        update_data: schemas.SearchResultUpdate,  # Changed signature
    ) -> schemas.SearchResultRead:  # MODIFIED return type
        """Updates an existing search result using Pydantic schema for update data, returns updated schema."""  # MODIFIED comment
        logger.debug(f"Updating search result id: {result_id!r}")

        with self.session_factory.begin() as session:
            try:
                db_search_result_model = self.search_repo.get_by_id(session, result_id)
                if not db_search_result_model:
                    raise repositories.RecordNotFoundError(
                        f"SearchResult with id {result_id} not found for update."
                    )

                update_dict = update_data.model_dump(exclude_unset=True)
                for key, value in update_dict.items():
                    if key == "screening_decision":
                        db_search_result_model.final_decision = value
                    else:
                        setattr(db_search_result_model, key, value)

                updated_result_from_repo_model = self.search_repo.update(
                    session, db_search_result_model
                )
                session.refresh(updated_result_from_repo_model)

                logger.info(
                    f"Successfully updated search result {updated_result_from_repo_model.id!r}"
                )
                search_result_read_fields = schemas.SearchResultRead.model_fields.keys()
                return schemas.SearchResultRead.model_validate(
                    {
                        field: getattr(updated_result_from_repo_model, field)
                        for field in search_result_read_fields
                        if hasattr(updated_result_from_repo_model, field)
                    }
                )

            except repositories.RecordNotFoundError:
                logger.warning(
                    f"Update failed for SearchResult {result_id!r}: Record not found."
                )
                raise
            except Exception as e_ex:
                logger.exception(f"Error updating search result {result_id!r}")
                raise ServiceError(
                    f"Failed to update search result: {e_ex!r}"
                ) from e_ex

    # Implement sync delete method
    def delete_search_result(self, result_id: uuid.UUID) -> None:
        """Deletes a search result by its ID."""
        logger.debug(f"Deleting search result id: {result_id!r}")
        # Use the internal session factory with .begin() for transaction management
        with self.session_factory.begin() as session:
            try:
                # Pass the internally managed session to the repository method
                self.search_repo.delete(session, result_id)
                # Commit is handled automatically by .begin() context manager on successful exit
                logger.info(f"Successfully deleted search result {result_id!r}")
            except repositories.RecordNotFoundError as e:
                logger.warning(
                    f"Delete failed for search result {result_id!r}: Record not found."
                )
                # Rollback is handled automatically by .begin() context manager on exception
            except Exception as e:
                logger.exception(f"Error deleting search result {result_id!r}")
                # Rollback is handled automatically by .begin() context manager on exception
                raise ServiceError("Failed to delete search result") from e


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
        review_dict = review_data.model_dump(exclude_none=True)

        # Ensure the ID from review_data (if provided) is used, otherwise generate a new one.
        if review_data.id is not None:
            review_dict["id"] = review_data.id
        elif "id" not in review_dict:  # If not in review_data and not in dump, generate
            review_dict["id"] = uuid.uuid4()

        review = models.SystematicReview.model_validate(review_dict)

        with self.session_factory.begin() as session:
            try:
                added_review = self.review_repo.add(session, review)
                session.refresh(added_review)
                logger.info(f"Successfully created review {added_review.id!r}")
                return added_review
            except Exception as e:
                logger.exception("Error creating review")
                raise ServiceError(f"Failed to create review: {e!r}") from e

    def get_review(self, review_id: uuid.UUID) -> models.SystematicReview | None:
        """Gets a review by its ID."""
        logger.debug(f"Getting review id: {review_id!r}")
        with self.session_factory.begin() as session:
            try:
                review = self.review_repo.get_by_id(session, review_id)
                if not review:
                    logger.warning(f"Review {review_id!r} not found.")
                    return None
                logger.debug(f"Found review {review_id!r}")
                return review
            except Exception as e:
                logger.exception(f"Error getting review {review_id!r}")
                raise ServiceError(f"Failed to get review: {e!r}") from e

    def get_all_reviews(self) -> Sequence[models.SystematicReview]:
        """Gets all reviews."""
        logger.debug("Getting all reviews")
        with self.session_factory.begin() as session:
            try:
                reviews = self.review_repo.get_all(session)
                logger.debug(f"Found {len(reviews)} reviews")
                return reviews
            except Exception as e:
                logger.exception("Error getting all reviews")
                raise ServiceError(f"Failed to get all reviews: {e!r}") from e

    def update_review(
        self, review_update_data: schemas.SystematicReviewUpdate, review_id: uuid.UUID
    ) -> models.SystematicReview:
        """Updates an existing review using partial update schema."""
        logger.debug(f"Updating review id: {review_id!r}")
        with self.session_factory.begin() as session:
            try:
                db_review = self.review_repo.get_by_id(session, review_id)
                if not db_review:
                    raise repositories.RecordNotFoundError(
                        f"Review {review_id!r} not found for update."
                    )
                update_dict = review_update_data.model_dump(exclude_unset=True)

                db_review.sqlmodel_update(update_dict)

                updated_review = self.review_repo.update(session, db_review)
                session.refresh(updated_review)

                logger.info(f"Successfully updated review {updated_review.id!r}")
                return updated_review
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Update failed: {e!r}")
                raise  # Re-raise specific error
            except Exception as e:
                logger.exception(f"Error updating review {review_id!r}")
                raise ServiceError(f"Failed to update review: {e!r}") from e

    def delete_review(self, review_id: uuid.UUID) -> None:
        """Deletes a review. Assumes DB cascade handles related data or it's handled separately."""
        logger.debug(f"Deleting review id: {review_id!r}")
        with self.session_factory.begin() as session:
            try:
                self.review_repo.delete(session, review_id)
                logger.info(f"Successfully deleted review {review_id!r}")
            except repositories.RecordNotFoundError as e:
                logger.warning(f"Delete failed: {e!r}")
            except Exception as e:
                logger.exception(f"Error deleting review {review_id!r}")
                raise ServiceError(f"Failed to delete review: {e!r}") from e


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
        response_metadata: dict[str, t.Any] | None = None,
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
                except Exception as e_refresh:
                    logger.warning(
                        f"Could not refresh screening result {result_id}: {e_refresh!r}"
                    )
                logger.info(f"Added screening result {added_result.id}")
                return added_result
            except Exception as e_add:
                logger.exception("Error adding screening result")
                session.rollback()
                raise ServiceError(
                    f"Failed to add screening result: {e_add!r}"
                ) from e_add

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
        except Exception as e_ex:  # Changed e to e_ex
            logger.exception(
                f"Error getting specific screening result by strategy: {e_ex!r}"
            )  # Use e_ex
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
        response_metadata: dict[str, t.Any] | None = None,
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

                result_id_to_use = run_id or (
                    existing_result.id if existing_result else uuid.uuid4()
                )

                # Prepare a dictionary with correctly typed values for model creation/update
                update_payload: dict[str, t.Any] = {
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
                    # SQLModel's sqlmodel_update is good for applying partial updates
                    # It updates the instance in place and returns None.
                    existing_result.sqlmodel_update(update_payload)
                    updated_result_model = (
                        existing_result  # Use the instance that was updated in place
                    )

                    updated_result = self.screen_repo.update(
                        session, updated_result_model
                    )
                else:
                    logger.info(
                        f"Creating new screening result with id {result_id_to_use}"
                    )
                    new_result_model = models.ScreenAbstractResult(
                        id=result_id_to_use, **update_payload
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
