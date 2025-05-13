"""Unit tests for the Service layer in sr_assistant.app.services."""

from __future__ import annotations

import typing as t  # Add import for typing as t
import uuid
from unittest.mock import (  # Import call for checking multiple calls
    MagicMock,
    patch,
)

import pytest
from sqlalchemy.orm import sessionmaker
from sqlmodel import Session

from sr_assistant.app import services
from sr_assistant.core import models, repositories, schemas  # Import schemas
from sr_assistant.core.types import SearchDatabaseSource


# --- Test Data ---
# Minimal mock for Bio.Entrez.Element.StringElement or similar simple elements
class MockEntrezStringElement:
    def __init__(self, value: t.Any, attributes: dict[t.Any, t.Any] | None = None):
        self._value = str(value)
        self.attributes = attributes if attributes is not None else {}

    def __str__(self) -> str:
        return self._value


# Mock for more complex Bio.Entrez.Element.DictElement or ListElement if needed
# For now, recursive_clean handles dicts and lists directly.

SAMPLE_RAW_PUBMED_RECORD_SIMPLE = {
    "MedlineCitation": {
        "PMID": MockEntrezStringElement("12345"),
        "Article": {
            "ArticleTitle": MockEntrezStringElement("Test Title Simple"),
            "Journal": {"Title": MockEntrezStringElement("Test Journal Simple")},
            "Abstract": {"AbstractText": [MockEntrezStringElement("Abstract simple.")]},
            "AuthorList": [  # List of dict-like structures for authors
                {
                    "LastName": MockEntrezStringElement("Doe"),
                    "ForeName": MockEntrezStringElement("John"),
                }
            ],
            "KeywordList": [
                MockEntrezStringElement("kw1"),
                MockEntrezStringElement("kw2"),
            ],  # Simple list of keywords
            "Journal": {
                "JournalIssue": {"PubDate": {"Year": MockEntrezStringElement("2023")}}
            },
        },
    },
    "PubmedData": {
        "ArticleIdList": [
            {
                "IdType": MockEntrezStringElement("doi"),
                "$": MockEntrezStringElement("10.123/test"),
            }
        ],
        "PublicationStatus": MockEntrezStringElement("epublish"),
    },
}

SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE = {
    "MedlineCitation": {
        "PMID": "12345",
        "Article": {
            "ArticleTitle": "Test Title Simple",
            "Journal": {
                "Title": "Test Journal Simple",
                "JournalIssue": {"PubDate": {"Year": "2023"}},
            },
            "Abstract": {"AbstractText": ["Abstract simple."]},
            "AuthorList": [{"LastName": "Doe", "ForeName": "John"}],
            "KeywordList": ["kw1", "kw2"],
        },
    },
    "PubmedData": {
        "ArticleIdList": [{"IdType": "doi", "$": "10.123/test"}],
        "PublicationStatus": "epublish",
    },
}


# --- Test SearchService Helper Methods ---


@pytest.fixture
def search_service_no_deps() -> services.SearchService:
    """Returns a SearchService instance with no real dependencies for helper testing."""
    mock_session_factory = MagicMock(spec=sessionmaker)
    # For helper methods that don't interact with the repo, repo can be None or MagicMock
    mock_repo = MagicMock(spec=repositories.SearchResultRepository)
    return services.SearchService(factory=mock_session_factory, search_repo=mock_repo)


class TestSearchServiceRecursiveClean:
    def test_recursive_clean_simple_types(
        self, search_service_no_deps: services.SearchService
    ):
        assert search_service_no_deps._recursive_clean("test") == "test"
        assert search_service_no_deps._recursive_clean(123) == 123
        assert search_service_no_deps._recursive_clean(None) is None

    def test_recursive_clean_list(self, search_service_no_deps: services.SearchService):
        data = ["a", MockEntrezStringElement("b"), [MockEntrezStringElement("c"), 1]]
        expected = ["a", "b", ["c", 1]]
        assert search_service_no_deps._recursive_clean(data) == expected

    def test_recursive_clean_dict(self, search_service_no_deps: services.SearchService):
        data = {
            "key1": MockEntrezStringElement("val1"),
            "key2": [
                MockEntrezStringElement("v2a"),
                {"subkey": MockEntrezStringElement("v2b")},
            ],
        }
        expected = {"key1": "val1", "key2": ["v2a", {"subkey": "v2b"}]}
        assert search_service_no_deps._recursive_clean(data) == expected

    def test_recursive_clean_with_entrez_elements(
        self, search_service_no_deps: services.SearchService
    ):
        # Test with a structure similar to what Entrez might return, using MockEntrezStringElement
        data = {
            "PMID": MockEntrezStringElement("123"),
            "Article": {
                "Title": MockEntrezStringElement("A Title"),
                "Authors": [
                    MockEntrezStringElement("Auth1"),
                    MockEntrezStringElement("Auth2"),
                ],
            },
        }
        expected = {
            "PMID": "123",
            "Article": {"Title": "A Title", "Authors": ["Auth1", "Auth2"]},
        }
        assert search_service_no_deps._recursive_clean(data) == expected


class TestSearchServiceExtractText:
    def test_extract_text_from_none(
        self, search_service_no_deps: services.SearchService
    ):
        assert search_service_no_deps._extract_text_from_element(None) == ""
        assert (
            search_service_no_deps._extract_text_from_element(None, default="def")
            == "def"
        )

    def test_extract_text_from_string(
        self, search_service_no_deps: services.SearchService
    ):
        assert search_service_no_deps._extract_text_from_element("hello") == "hello"

    def test_extract_text_from_entrez_element(
        self, search_service_no_deps: services.SearchService
    ):
        element = MockEntrezStringElement("value_from_entrez")
        assert (
            search_service_no_deps._extract_text_from_element(element)
            == "value_from_entrez"
        )

    def test_extract_text_from_other_types(
        self, search_service_no_deps: services.SearchService
    ):
        assert search_service_no_deps._extract_text_from_element(123) == "123"
        assert search_service_no_deps._extract_text_from_element([1, 2]) == "[1, 2]"


# TODO: Add tests for _parse_pubmed_ids, _parse_pubmed_title_abstract, etc.
# TODO: Add tests for _map_pubmed_to_search_result (using cleaned data)
# TODO: Add tests for search_pubmed_and_store_results (mocking Entrez and repo)


class TestSearchServiceParsePubmedIds:
    def test_parse_ids_basic(self, search_service_no_deps: services.SearchService):
        pmid, doi, pmc = search_service_no_deps._parse_pubmed_ids(
            SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )
        assert pmid == "12345"
        assert doi == "10.123/test"
        assert pmc is None  # Not in simple sample

    def test_parse_ids_missing_data(
        self, search_service_no_deps: services.SearchService
    ):
        empty_record = {}
        pmid, doi, pmc = search_service_no_deps._parse_pubmed_ids(empty_record)
        assert pmid is None
        assert doi is None
        assert pmc is None

    # Add more specific cases, e.g., DOI in ELocationID, PMC present etc.


class TestSearchServiceParsePubmedTitleAbstract:
    def test_parse_title_abstract_basic(
        self, search_service_no_deps: services.SearchService
    ):
        title, abstract = search_service_no_deps._parse_pubmed_title_abstract(
            SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )
        assert title == "Test Title Simple"
        assert abstract == "Abstract simple."

    def test_parse_title_abstract_missing(
        self, search_service_no_deps: services.SearchService
    ):
        record = {
            "MedlineCitation": {"Article": {}}
        }  # Missing Title and Abstract sections
        title, abstract = search_service_no_deps._parse_pubmed_title_abstract(record)
        assert title is None
        assert abstract is None

    def test_parse_title_abstract_complex_abstract(
        self, search_service_no_deps: services.SearchService
    ):
        record = {
            "MedlineCitation": {
                "Article": {
                    "ArticleTitle": "Complex Abstract Test",
                    "Abstract": {
                        "AbstractText": [
                            {"Label": "BACKGROUND", "$": "Background text."},
                            "Followed by plain text.",
                            {"$": "Text without label."},
                        ]
                    },
                }
            }
        }
        title, abstract = search_service_no_deps._parse_pubmed_title_abstract(record)
        assert title == "Complex Abstract Test"
        assert (
            abstract
            == "BACKGROUND: Background text. Followed by plain text. Text without label."
        )


class TestSearchServiceParsePubmedJournalYear:
    def test_parse_journal_year_basic(
        self, search_service_no_deps: services.SearchService
    ):
        journal, year = search_service_no_deps._parse_pubmed_journal_year(
            SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )
        assert journal == "Test Journal Simple"
        assert year == "2023"

    def test_parse_journal_year_medline_date(
        self, search_service_no_deps: services.SearchService
    ):
        record = {
            "MedlineCitation": {
                "Article": {
                    "Journal": {
                        "JournalIssue": {"PubDate": {"MedlineDate": "2022 Spring"}}
                    }
                }
            }
        }
        journal, year = search_service_no_deps._parse_pubmed_journal_year(record)
        assert year == "2022"


class TestSearchServiceParsePubmedAuthors:
    def test_parse_authors_basic(self, search_service_no_deps: services.SearchService):
        authors = search_service_no_deps._parse_pubmed_authors(
            SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )
        assert authors == ["John Doe"]

    def test_parse_authors_multiple_and_collective(
        self, search_service_no_deps: services.SearchService
    ):
        record = {
            "MedlineCitation": {
                "Article": {
                    "AuthorList": [
                        {"LastName": "Smith", "Initials": "J"},
                        {"ForeName": "Alice", "LastName": "Wonder"},
                        {"CollectiveName": "A Study Group"},
                    ]
                }
            }
        }
        authors = search_service_no_deps._parse_pubmed_authors(record)
        assert authors == ["Smith J", "Alice Wonder", "A Study Group"]

    def test_parse_authors_empty(self, search_service_no_deps: services.SearchService):
        record = {"MedlineCitation": {"Article": {"AuthorList": []}}}
        authors = search_service_no_deps._parse_pubmed_authors(record)
        assert authors is None


class TestSearchServiceParsePubmedKeywords:
    def test_parse_keywords_basic(self, search_service_no_deps: services.SearchService):
        keywords = search_service_no_deps._parse_pubmed_keywords(
            SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )
        assert sorted(keywords) == sorted(["kw1", "kw2"])  # Compare sorted lists

    def test_parse_keywords_complex_structure(
        self, search_service_no_deps: services.SearchService
    ):
        record = {
            "MedlineCitation": {
                "Article": {
                    "KeywordList": [
                        ["topic A", "topic B"],  # list of lists
                        {"Keyword": ["topic C", "topic D"]},  # dict with list
                        "topic E",  # direct string
                    ]
                }
            }
        }
        keywords = search_service_no_deps._parse_pubmed_keywords(record)
        assert isinstance(keywords, list)
        # Order might not be guaranteed due to set for deduplication, so check content
        assert sorted(keywords) == sorted(
            ["topic A", "topic B", "topic C", "topic D", "topic E"]
        )

    def test_parse_keywords_empty(self, search_service_no_deps: services.SearchService):
        record = {"MedlineCitation": {"KeywordList": []}}
        keywords = search_service_no_deps._parse_pubmed_keywords(record)
        assert keywords is None


# TODO: Add tests for _map_pubmed_to_search_result (using cleaned data)
# TODO: Add tests for search_pubmed_and_store_results (mocking Entrez and repo)


class TestSearchServiceMapPubmedToSearchResult:
    def test_map_pubmed_basic(self, search_service_no_deps: services.SearchService):
        review_id = uuid.uuid4()
        mapped_result = search_service_no_deps._map_pubmed_to_search_result(
            review_id, SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        )

        assert isinstance(mapped_result, models.SearchResult)
        assert mapped_result.review_id == review_id
        assert mapped_result.source_db == SearchDatabaseSource.PUBMED
        assert mapped_result.source_id == "12345"
        assert mapped_result.doi == "10.123/test"
        assert mapped_result.title == "Test Title Simple"
        assert mapped_result.abstract == "Abstract simple."
        assert mapped_result.journal == "Test Journal Simple"
        assert mapped_result.year == "2023"
        assert mapped_result.authors == ["John Doe"]
        assert sorted(mapped_result.keywords) == sorted(
            ["kw1", "kw2"]
        )  # Compare sorted lists
        assert mapped_result.raw_data == SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE
        assert mapped_result.source_metadata.get("pmc") is None
        assert mapped_result.source_metadata.get("publication_status") == "epublish"
        assert (
            mapped_result.source_metadata.get("mesh_headings") == []
        )  # From _map_pubmed_to_search_result defaults

    def test_map_pubmed_missing_pmid_or_title(
        self, search_service_no_deps: services.SearchService
    ):
        review_id = uuid.uuid4()
        record_no_pmid = SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE.copy()
        record_no_pmid["MedlineCitation"] = record_no_pmid["MedlineCitation"].copy()
        record_no_pmid["MedlineCitation"]["PMID"] = ""  # Empty PMID
        assert (
            search_service_no_deps._map_pubmed_to_search_result(
                review_id, record_no_pmid
            )
            is None
        )

        record_no_title = SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE.copy()
        record_no_title["MedlineCitation"] = record_no_title["MedlineCitation"].copy()
        record_no_title["MedlineCitation"]["Article"] = record_no_title[
            "MedlineCitation"
        ]["Article"].copy()
        record_no_title["MedlineCitation"]["Article"]["ArticleTitle"] = (
            ""  # Empty Title
        )
        assert (
            search_service_no_deps._map_pubmed_to_search_result(
                review_id, record_no_title
            )
            is None
        )

    def test_map_pubmed_with_mesh_headings(
        self, search_service_no_deps: services.SearchService
    ):
        review_id = uuid.uuid4()
        record_with_mesh = SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE.copy()
        record_with_mesh["MedlineCitation"] = record_with_mesh["MedlineCitation"].copy()
        mesh_data = [{"DescriptorName": "COVID-19"}, {"DescriptorName": "Vaccines"}]
        record_with_mesh["MedlineCitation"]["MeshHeadingList"] = mesh_data

        # We need to ensure the parsing helpers for MeSH are robust or this data is passed correctly.
        # The _map_pubmed_to_search_result currently gets MeshHeadingList directly from api_record["MedlineCitation"]
        # So, the cleaned version should still have this structure.

        mapped_result = search_service_no_deps._map_pubmed_to_search_result(
            review_id, record_with_mesh
        )
        assert isinstance(mapped_result, models.SearchResult)
        assert mapped_result.source_metadata.get("mesh_headings") == mesh_data


# TODO: Add tests for search_pubmed_and_store_results (mocking Entrez and repo)


# This test does not need the Entrez class-level patch, move it outside the class
# and use search_service_no_deps for a clean service instance.
def test_search_pubmed_no_ncbi_email_standalone(
    search_service_no_deps: services.SearchService,
):
    # This test uses search_service_no_deps which is fine, no Entrez mocking needed here
    with patch("sr_assistant.app.services.os.getenv") as mock_getenv:
        mock_getenv.side_effect = lambda key, default=None: {
            "NCBI_API_KEY": "testkey123"
        }.get(key, default)  # Email is missing

        with pytest.raises(
            services.ServiceError, match="NCBI_EMAIL environment variable not set"
        ):
            search_service_no_deps.search_pubmed_and_store_results(
                uuid.uuid4(), "query"
            )


# Remove class-level patch, use mocker fixture instead
# @patch("sr_assistant.app.services.Entrez", new_callable=MagicMock)
class TestSearchServiceSearchPubmedAndStoreResults:
    @pytest.fixture
    def search_service_and_mocks(
        self, mocker: MagicMock
    ) -> tuple[services.SearchService, MagicMock, MagicMock]:
        """Provides a SearchService instance with mocked Entrez, os.getenv, session factory, and repo."""
        # 1. Mock os.getenv for NCBI credentials
        mock_getenv = mocker.patch("sr_assistant.app.services.os.getenv")
        mock_getenv.side_effect = lambda key, default=None: {
            "NCBI_EMAIL": "test@example.com",
            "NCBI_API_KEY": "testkey123",
        }.get(key, default)

        # 2. Mock Entrez module
        mock_entrez = mocker.patch("sr_assistant.app.services.Entrez")
        mock_entrez.email = ""
        mock_entrez.api_key = None
        mock_esearch_handle = MagicMock()
        mock_entrez.esearch.return_value = mock_esearch_handle
        # Default read sequence: first for esearch, second for efetch
        mock_entrez.read.side_effect = [
            {"IdList": ["pmid1", "pmid2"]},
            [
                SAMPLE_CLEANED_PUBMED_RECORD_SIMPLE,
                {
                    "MedlineCitation": {
                        "PMID": "pmid2",
                        "Article": {"ArticleTitle": "Title 2"},
                    },
                    "PubmedData": {},
                },
            ],
        ]

        # 3. Mock session factory and session
        mock_session_factory = MagicMock(spec=sessionmaker)
        mock_session = MagicMock(spec=Session)
        mock_session_factory.begin.return_value.__enter__.return_value = mock_session

        # 4. Mock repository
        mock_repo = MagicMock(spec=repositories.SearchResultRepository)
        mock_repo.add_all.side_effect = lambda session, items: items

        service_instance = services.SearchService(
            factory=mock_session_factory, search_repo=mock_repo
        )
        return (
            service_instance,
            mock_entrez,
            mock_repo,
        )  # Return repo mock too for assertions

    def test_search_pubmed_success(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, mock_repo = search_service_and_mocks

        review_id = uuid.uuid4()
        query = "test query"
        max_r = 2

        results = search_service.search_pubmed_and_store_results(
            review_id, query, max_results=max_r
        )

        assert len(results) == 2
        assert results[0].source_id == "12345"
        assert results[1].source_id == "pmid2"

        mock_entrez.esearch.assert_called_once_with(
            db="pubmed", term=query, retmax=max_r, sort="relevance"
        )
        mock_entrez.efetch.assert_called_once_with(
            db="pubmed", id=["pmid1", "pmid2"], rettype="medline", retmode="xml"
        )
        mock_repo.add_all.assert_called_once()
        search_service.session_factory.begin.return_value.__exit__.assert_called_once()

    def test_search_pubmed_esearch_fails(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, _ = search_service_and_mocks
        mock_entrez.esearch.side_effect = Exception("ESearch Boom!")
        with pytest.raises(
            services.ServiceError, match=r"PubMed search \(esearch\) failed"
        ):
            search_service.search_pubmed_and_store_results(uuid.uuid4(), "query")

    def test_search_pubmed_efetch_fails(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, _ = search_service_and_mocks
        # Ensure esearch is successful, then efetch (via Entrez.read for efetch) fails
        mock_entrez.read.side_effect = [{"IdList": ["1"]}, Exception("EFetch Boom!")]
        with pytest.raises(
            services.ServiceError, match=r"PubMed fetch \(efetch\) failed"
        ):
            search_service.search_pubmed_and_store_results(uuid.uuid4(), "query")

    def test_search_pubmed_no_pmids_found(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, _ = search_service_and_mocks
        mock_entrez.read.side_effect = [
            {"IdList": []}
        ]  # Mock for esearch returning no PMIDs
        results = search_service.search_pubmed_and_store_results(uuid.uuid4(), "query")
        assert len(results) == 0
        mock_entrez.efetch.assert_not_called()

    def test_search_pubmed_repo_add_all_fails_constraint(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, mock_repo = search_service_and_mocks
        # Default Entrez mock is fine (returns 2 PMIDs and records)
        mock_repo.add_all.side_effect = repositories.ConstraintViolationError(
            "DB constraint fail"
        )
        results = search_service.search_pubmed_and_store_results(uuid.uuid4(), "query")
        assert results == []

    def test_search_pubmed_mapping_returns_none(
        self,
        search_service_and_mocks: tuple[services.SearchService, MagicMock, MagicMock],
    ):
        search_service, mock_entrez, _ = search_service_and_mocks
        # Default Entrez mock returns data for 2 PMIDs / records
        with patch.object(
            search_service, "_map_pubmed_to_search_result", return_value=None
        ) as mock_mapper:
            results = search_service.search_pubmed_and_store_results(
                uuid.uuid4(), "query"
            )
            assert len(results) == 0
            assert (
                mock_mapper.call_count == 2
            )  # Called for each of the 2 records from mock Entrez.read

    # The test_search_pubmed_no_ncbi_email is now a standalone function test_search_pubmed_no_ncbi_email_standalone


# --- Test ReviewService Methods ---


@pytest.fixture
def review_service_with_mocks(
    mocker: MagicMock,
) -> tuple[services.ReviewService, MagicMock, MagicMock]:
    """Provides a ReviewService instance with mocked session factory and repository."""
    mock_session_factory = MagicMock(spec=sessionmaker)
    mock_session = MagicMock(spec=Session)
    mock_session_factory.begin.return_value.__enter__.return_value = mock_session

    mock_review_repo = MagicMock(spec=repositories.SystematicReviewRepository)

    service = services.ReviewService(
        factory=mock_session_factory, review_repo=mock_review_repo
    )
    return service, mock_session_factory, mock_review_repo


class TestReviewService:
    def test_create_review(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_data_create = schemas.SystematicReviewCreate(
            research_question="Test RQ", exclusion_criteria="Test Exclusion"
        )
        # The model_validate part happens inside the service method
        # Mock the repository's add method
        mock_created_review = models.SystematicReview(
            id=uuid.uuid4(),
            research_question="Test RQ",  # Use correct spelling
            exclusion_criteria="Test Exclusion",
        )
        mock_review_repo.add.return_value = mock_created_review

        created_review = service.create_review(review_data_create)

        mock_review_repo.add.assert_called_once()
        args, _ = mock_review_repo.add.call_args
        assert isinstance(args[1], models.SystematicReview)
        assert (
            args[1].research_question == review_data_create.research_question
        )  # Use correct spelling

        mock_session_factory.begin.return_value.__enter__.return_value.refresh.assert_called_with(
            mock_created_review
        )
        assert created_review == mock_created_review
        mock_session_factory.begin.return_value.__exit__.assert_called_once()  # Verifies commit/rollback context

    def test_get_review_found(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        mock_existing_review = models.SystematicReview(
            id=review_id, research_question="Existing RQ", exclusion_criteria="Excl"
        )
        mock_review_repo.get_by_id.return_value = mock_existing_review

        retrieved_review = service.get_review(review_id)

        mock_review_repo.get_by_id.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        assert retrieved_review == mock_existing_review
        mock_session_factory.begin.return_value.__exit__.assert_called_once()

    def test_get_review_not_found(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        mock_review_repo.get_by_id.return_value = None

        retrieved_review = service.get_review(review_id)

        mock_review_repo.get_by_id.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        assert retrieved_review is None
        mock_session_factory.begin.return_value.__exit__.assert_called_once()

    # TODO: Add tests for get_all_reviews, update_review, delete_review
    def test_get_all_reviews(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        mock_reviews_list = [
            models.SystematicReview(
                id=uuid.uuid4(), research_question="RQ1", exclusion_criteria="E1"
            ),
            models.SystematicReview(
                id=uuid.uuid4(), research_question="RQ2", exclusion_criteria="E2"
            ),
        ]
        mock_review_repo.get_all.return_value = mock_reviews_list

        reviews = service.get_all_reviews()

        mock_review_repo.get_all.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value
        )
        assert reviews == mock_reviews_list
        mock_session_factory.begin.return_value.__exit__.assert_called_once()

    def test_update_review_success(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        update_data_schema = schemas.SystematicReviewUpdate(
            research_question="Updated RQ"
        )

        mock_existing_review = models.SystematicReview(
            id=review_id,
            research_question="Original RQ",
            exclusion_criteria="Orig Excl",  # Use correct spelling
        )
        mock_review_repo.get_by_id.return_value = mock_existing_review

        # Assume repo.update returns the updated model (which it should, or the same instance)
        # and that it handles the actual update logic based on sqlmodel_update called in service.
        # The service calls repo.update then refreshes.
        mock_review_repo.update.return_value = (
            mock_existing_review  # Simulate update returning the instance
        )

        updated_review = service.update_review(update_data_schema, review_id)

        mock_review_repo.get_by_id.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        # Check that sqlmodel_update was effectively called (by checking a changed attribute)
        # This relies on the mock_existing_review object being modified.
        # The service uses update_data_schema.model_dump(exclude_unset=True), which will have 'research_question'.
        # SQLModel's sqlmodel_update should handle setting the correct field.
        assert (
            mock_existing_review.research_question == "Updated RQ"
        )  # Assert correctly spelled field
        mock_review_repo.update.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value,
            mock_existing_review,
        )
        mock_session_factory.begin.return_value.__enter__.return_value.refresh.assert_called_with(
            mock_existing_review
        )
        assert updated_review == mock_existing_review
        mock_session_factory.begin.return_value.__exit__.assert_called_once()

    def test_update_review_not_found(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        update_data_schema = schemas.SystematicReviewUpdate(
            research_question="Updated RQ"
        )
        mock_review_repo.get_by_id.return_value = None

        with pytest.raises(repositories.RecordNotFoundError):
            service.update_review(update_data_schema, review_id)

        mock_review_repo.get_by_id.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        mock_review_repo.update.assert_not_called()
        mock_session_factory.begin.return_value.__exit__.assert_called_once()  # For rollback

    def test_delete_review_success(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        mock_review_repo.delete.return_value = (
            None  # Delete doesn't typically return a value
        )

        service.delete_review(review_id)

        mock_review_repo.delete.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        mock_session_factory.begin.return_value.__exit__.assert_called_once()

    def test_delete_review_not_found(
        self,
        review_service_with_mocks: tuple[services.ReviewService, MagicMock, MagicMock],
    ):
        service, mock_session_factory, mock_review_repo = review_service_with_mocks
        review_id = uuid.uuid4()
        mock_review_repo.delete.side_effect = repositories.RecordNotFoundError(
            "Not found"
        )

        # The service currently suppresses RecordNotFoundError for delete and logs it.
        # So, we expect no exception here.
        service.delete_review(review_id)

        mock_review_repo.delete.assert_called_once_with(
            mock_session_factory.begin.return_value.__enter__.return_value, review_id
        )
        mock_session_factory.begin.return_value.__exit__.assert_called_once()  # For rollback due to error

    # TODO: Add tests for get_all_reviews, update_review, delete_review
