"""Unit tests for PubMed integration functions."""

import typing as t

import pytest
from sr_assistant.app.pubmed_integration import (
    extract_text,
    parse_pubmed_systematic_review_data,
    recursive_clean,
)


class TestExtractText:
    """Test cases for the extract_text function."""

    def test_dict_with_value(self) -> None:
        """Test extracting text from a dict with _value key."""
        element: dict[str, t.Any] = {"_value": "test value", "other": "data"}
        assert extract_text(element) == "test value"

    def test_dict_without_value(self) -> None:
        """Test extracting text from a dict without _value key."""
        element: dict[str, str] = {"data": "test"}
        assert extract_text(element) == "{'data': 'test'}"

    def test_object_with_value(self) -> None:
        """Test extracting text from an object with _value attribute."""

        class TestElement:
            _value = "object value"

        assert extract_text(TestElement()) == "object value"

    def test_plain_string(self) -> None:
        """Test extracting text from a plain string."""
        assert extract_text("plain string") == "plain string"

    def test_none_value(self) -> None:
        """Test extracting text from None."""
        assert extract_text(None) == "None"


class TestRecursiveClean:
    """Test cases for the recursive_clean function."""

    def test_list_cleaning(self) -> None:
        """Test cleaning a list of mixed elements."""
        input_data: list[t.Any] = [
            {"_value": "value1"},
            "plain string",
            {"nested": {"_value": "nested value"}},
        ]
        result = recursive_clean(input_data)
        assert result[0] == {"_value": "value1"}
        assert result[1] == "plain string"
        assert result[2] == {"nested": {"_value": "nested value"}}

    def test_dict_cleaning(self) -> None:
        """Test cleaning a dictionary with nested elements."""
        input_data: dict[str, t.Any] = {
            "key1": {"_value": "value1"},
            "key2": "plain string",
            "key3": {"nested": {"_value": "nested value"}},
        }
        result = recursive_clean(input_data)
        assert result["key1"] == {"_value": "value1"}
        assert result["key2"] == "plain string"
        assert result["key3"] == {"nested": {"_value": "nested value"}}

    def test_object_with_attributes(self) -> None:
        """Test cleaning an object with attributes."""

        class TestElement:
            attributes: dict[str, t.Any] = {}
            _value = "test value"

        element = TestElement()
        assert recursive_clean(element) == "test value"

    def test_deep_nested_structure(self) -> None:
        """Test cleaning a deeply nested structure."""
        input_data: dict[str, list[t.Any]] = {
            "level1": [
                {"_value": "v1"},
                {
                    "level2": {
                        "array": [{"_value": "v2"}, {"_value": "v3"}],
                        "_value": "v4",
                    }
                },
            ]
        }
        result = recursive_clean(input_data)
        assert result["level1"][0] == {"_value": "v1"}
        assert "_value" in result["level1"][1]["level2"]
        assert result["level1"][1]["level2"]["_value"] == "v4"


class MockArticleId:
    """Mock class for article IDs with attributes."""

    def __init__(self, value: str, id_type: str) -> None:
        """Initialize with value and ID type."""
        self._value = value
        self.attributes = {"IdType": id_type}


@pytest.fixture
def sample_pubmed_article() -> dict[str, list[dict[str, t.Any]]]:
    """Fixture providing a sample PubMed article record."""
    return {
        "PubmedArticle": [
            {
                "MedlineCitation": {
                    "PMID": {"_value": "12345678"},
                    "Article": {
                        "ArticleTitle": "Test Article Title",
                        "Abstract": {
                            "AbstractText": [
                                {
                                    "_value": "Background section",
                                    "attributes": {"Label": "BACKGROUND"},
                                },
                                {
                                    "_value": "Methods section",
                                    "attributes": {"Label": "METHODS"},
                                },
                            ]
                        },
                        "Journal": {
                            "Title": "Test Journal",
                            "JournalIssue": {
                                "Volume": "10",
                                "Issue": "2",
                                "PubDate": {
                                    "Year": "2025",
                                    "Month": "Feb",
                                    "Day": "15",
                                },
                            },
                        },
                        "AuthorList": [
                            {
                                "LastName": {"_value": "Smith"},
                                "ForeName": {"_value": "John"},
                                "Initials": {"_value": "J"},
                                "AffiliationInfo": [
                                    {"Affiliation": {"_value": "Test University"}}
                                ],
                            }
                        ],
                        "PublicationTypeList": [
                            {"_value": "Journal Article"},
                            {"_value": "Review"},
                        ],
                    },
                    "KeywordList": [
                        [
                            {"_value": "Keyword1"},
                            {"_value": "Keyword2"},
                        ]
                    ],
                    "MeshHeadingList": [
                        {
                            "DescriptorName": {"_value": "Test Mesh Term"},
                            "QualifierName": [{"_value": "methods"}],
                        }
                    ],
                },
                "PubmedData": {
                    "ArticleIdList": [
                        MockArticleId("10.1234/test", "doi"),
                        MockArticleId("PMC123456", "pmc"),
                    ]
                },
            }
        ]
    }


class TestParsePubmedSystematicReviewData:
    """Test cases for parse_pubmed_systematic_review_data function."""

    def test_complete_article_parsing(
        self, sample_pubmed_article: dict[str, list[dict[str, t.Any]]]
    ) -> None:
        """Test parsing a complete PubMed article with all fields."""
        result = parse_pubmed_systematic_review_data(sample_pubmed_article)
        pmid = extract_text(
            sample_pubmed_article["PubmedArticle"][0]["MedlineCitation"]["PMID"]
        )

        assert result["pmid"] == pmid
        assert result["article_title"] == "Test Article Title"
        assert "Background section" in result["abstract"]
        assert "Methods section" in result["abstract"]
        assert result["journal_title"] == "Test Journal"
        assert result["journal_pub_date"]["year"] == "2025"
        assert result["journal_pub_date"]["month"] == "Feb"
        assert result["journal_pub_date"]["day"] == "15"

        # Check author parsing
        assert len(result["authors"]) == 1
        author = result["authors"][0]
        assert author["last_name"] == "Smith"
        assert author["fore_name"] == "John"
        assert author["initials"] == "J"
        assert author["affiliations"] == ["Test University"]

        # Check IDs
        assert result["doi"] == "10.1234/test"
        assert result["pmcid"] == "PMC123456"

        # Check publication types
        assert "Journal Article" in result["publication_types"]
        assert "Review" in result["publication_types"]

        # Check keywords
        assert "Keyword1" in result["keywords"]
        assert "Keyword2" in result["keywords"]

        # Check MeSH terms
        assert len(result["mesh_headings"]) == 1
        assert result["mesh_headings"][0] == "Test Mesh Term"

    def test_empty_record(self) -> None:
        """Test parsing an empty record."""
        assert parse_pubmed_systematic_review_data({}) == {}

    def test_minimal_record(self) -> None:
        """Test parsing a minimal record with only required fields."""
        minimal_record: dict[str, list[dict[str, t.Any]]] = {
            "PubmedArticle": [
                {
                    "MedlineCitation": {
                        "PMID": {"_value": "12345"},
                        "Article": {
                            "ArticleTitle": "Minimal Title",
                            "Journal": {"Title": "Test Journal"},
                        },
                    },
                    "PubmedData": {},
                }
            ]
        }
        result = parse_pubmed_systematic_review_data(minimal_record)
        assert result["pmid"] == "12345"
        assert result["article_title"] == "Minimal Title"
        assert result["journal_title"] == "Test Journal"
        assert result["abstract"] == ""  # Empty when not provided
        assert result["authors"] == []  # Empty when not provided
        assert result["keywords"] == []  # Empty when not provided

    def test_multiple_abstract_sections(
        self, sample_pubmed_article: dict[str, list[dict[str, t.Any]]]
    ) -> None:
        """Test parsing multiple abstract sections."""
        result = parse_pubmed_systematic_review_data(sample_pubmed_article)
        assert "Background section" in result["abstract"]
        assert "Methods section" in result["abstract"]
        assert "\n\n" in result["abstract"]  # Sections should be separated

    @pytest.mark.parametrize(
        "id_type,id_value,expected_key",
        [
            ("doi", "10.1234/test", "doi"),
            ("pmc", "PMC123456", "pmcid"),
        ],
    )
    def test_article_id_extraction(
        self,
        sample_pubmed_article: dict[str, list[dict[str, t.Any]]],
        id_type: str,
        id_value: str,
        expected_key: str,
    ) -> None:
        """Test extraction of different article IDs."""
        sample_pubmed_article["PubmedArticle"][0]["PubmedData"]["ArticleIdList"] = [
            MockArticleId(id_value, id_type)
        ]
        result = parse_pubmed_systematic_review_data(sample_pubmed_article)
        assert result[expected_key] == id_value
