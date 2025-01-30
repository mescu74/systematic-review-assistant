"""Tests for pubmed_integration module."""

from __future__ import annotations

import pytest

from step2.pubmed_integration import (
    extract_article_info,
    pubmed_fetch_details,
    pubmed_search,
)


@pytest.mark.integration
def test_pubmed_search() -> None:
    """Test PubMed search functionality."""
    query = "systematic review AND machine learning"
    pmids = pubmed_search(query, max_results=5)

    assert isinstance(pmids, list)
    assert len(pmids) <= 5
    assert all(isinstance(pmid, str) for pmid in pmids)


@pytest.mark.integration
def test_pubmed_fetch() -> None:
    """Test fetching article details."""
    # Use a small sample of PMIDs
    pmids = ["123456", "789012"]
    records = pubmed_fetch_details(pmids)

    assert isinstance(records, dict)
    assert "PubmedArticle" in records


@pytest.mark.integration
def test_extract_article_info() -> None:
    """Test article info extraction."""
    # Mock article data
    article = {
         "MedlineCitation": {
             "PMID": "39826015",
             "Article": {
                 "ArticleTitle": "Progress report on multiple endocrine neoplasia type 1.",
                 "Abstract": {
                     "AbstractText": [
                         "Multiple endocrine neoplasia type 1 (MEN1) syndrome is an ",
                         "autosomal dominant disorder caused by a germline pathogenic ",
                         "variant in the MEN1 tumor suppressor gene. Patients with MEN1 ",
                         "have a high risk for primary hyperparathyroidism (PHPT) with a ",
                         "penetrance of nearly 100%, pituitary adenomas (PitAd) in 40% of ",
                         "patients, and neuroendocrine neoplasms (NEN) of the pancreas ",
                         "with a mainline mortality of 40%...",
                     ]
                 },
                 "Journal": {
                     "Title": "Familial cancer",
                 },
                 "JournalIssue": {
                     "PubDate": {
                         "Year": "2024",
                     },
                 },
             },
         },
         "PubmedData": {
             "ArticleIdList": [ ],
         },
     }
    info = extract_article_info(article)

    assert isinstance(info, dict)
    assert info["pmid"] == article["MedlineCitation"]["PMID"]
    assert info["title"] == article["MedlineCitation"]["Article"]["ArticleTitle"]
    assert info["abstract"] == "".join(article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"])
    assert info["journal"] == article["MedlineCitation"]["Article"]["Journal"]["Title"]
    assert info["year"] == article["MedlineCitation"]["Article"]["JournalIssue"]["PubDate"]["Year"]


@pytest.mark.integration
def test_pubmed_error_handling() -> None:
    """Test error handling in PubMed functions."""
    with pytest.raises(ValueError):
        pubmed_fetch_details([])  # Empty PMID list should raise ValueError
