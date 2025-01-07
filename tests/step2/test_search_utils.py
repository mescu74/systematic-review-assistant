from __future__ import annotations

import pytest

from step2.search_utils import build_pubmed_query


def test_basic_keyword_search() -> None:
    keywords = ["diabetes", "obesity"]
    result = build_pubmed_query(keywords=keywords, inclusion_criteria={}, exclusion_criteria={})
    expected = "(diabetes) OR (obesity)"
    assert result == expected


def test_with_inclusion_criteria() -> None:
    keywords = ["cancer"]
    inclusion_criteria = {"language": "english", "study_type": "clinical trial"}
    result = build_pubmed_query(
        keywords=keywords,
        inclusion_criteria=inclusion_criteria,
        exclusion_criteria={},
    )
    expected = "(cancer) AND english[lang] AND clinical trial[pt]"
    assert result == expected


def test_with_exclusion_criteria() -> None:
    keywords = ["alzheimer"]
    exclusion_criteria = {"exclude_study_type": "review"}
    result = build_pubmed_query(
        keywords=keywords,
        inclusion_criteria={},
        exclusion_criteria=exclusion_criteria,
    )
    expected = "(alzheimer) AND NOT review[pt]"
    assert result == expected


def test_with_additional_filters() -> None:
    keywords = ["covid"]
    additional_filters = {"date_range": "2020:2023"}
    result = build_pubmed_query(
        keywords=keywords,
        inclusion_criteria={},
        exclusion_criteria={},
        additional_filters=additional_filters,
    )
    expected = "(covid) AND 2020:2023[dp]"
    assert result == expected


def test_empty_inputs() -> None:
    with pytest.raises(ValueError, match="At least one keyword is required"):
        build_pubmed_query(keywords=[], inclusion_criteria={}, exclusion_criteria={})
