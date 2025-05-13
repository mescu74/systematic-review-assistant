"""Integration tests for the Service layer."""

from __future__ import annotations

import uuid  # Needed for test_review_with_criteria
from collections.abc import Sequence

import pytest
from sqlalchemy.orm import sessionmaker  # Import sessionmaker
from sqlmodel import Session  # Keep Session from sqlmodel

# Imports used by test and fixtures
from sr_assistant.app import services
from sr_assistant.core import models, repositories, schemas
from sr_assistant.core.types import CriteriaFramework

# Local test_db_session fixture is removed, conftest.py provides db_session


@pytest.fixture(scope="function")
def search_service_test_session_factory(db_session: Session) -> sessionmaker[Session]:
    """Creates a sessionmaker bound to the integration test DB engine."""
    test_engine = (
        db_session.get_bind()
    )  # Engine is obtained from the db_session fixture
    return sessionmaker(
        bind=test_engine,
        class_=Session,  # Use SQLModel Session
        expire_on_commit=False,
    )


@pytest.fixture(scope="function")
def test_review_with_criteria(db_session: Session) -> models.SystematicReview:
    """Creates and persists a SystematicReview with PICO criteria for testing."""
    review_id_val = uuid.uuid4()  # Use a local var for clarity
    pico_population = "Adults with Type 2 Diabetes"
    pico_intervention = "Metformin 1000mg"
    pico_comparison = "Placebo"
    pico_outcome = "HbA1c reduction"

    inclusion_criteria_parts = [
        f"Population: {pico_population}",
        f"Intervention: {pico_intervention}",
        f"Comparison: {pico_comparison}",
        f"Outcome: {pico_outcome}",
    ]
    inclusion_criteria_str = "; ".join(inclusion_criteria_parts)
    exclusion_criteria_str = "Studies on Type 1 Diabetes; Non-English studies"

    # Note: SystematicReview model has 'reearch_question' (typo)
    # Assuming the model field is 'reearch_question' as per models.py
    # If it's corrected to 'research_question' in models.py, this needs to match.
    # For now, using the typo 'reearch_question' from the provided models.py
    review_create_schema = schemas.SystematicReviewCreate(
        research_question=f"Integration Test RQ {review_id_val}",  # This goes to SystematicReviewCreate.research_question
        background="Test background for review with PICO criteria.",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": pico_population,
            "intervention": pico_intervention,
            "comparison": pico_comparison,
            "outcome": pico_outcome,
        },
        inclusion_criteria=inclusion_criteria_str,
        exclusion_criteria=exclusion_criteria_str,
    )

    # SQLModel.model_validate needs to be called on the SQLModel class
    # And the SystematicReview SQLModel expects 'reearch_question'
    # We need to map SystematicReviewCreate to SystematicReview, handling the typo

    review_model_data = review_create_schema.model_dump()

    # Fields not in SystematicReview model (if any from SystematicReviewCreate that don't map)
    # would ideally be handled by excluding them from model_dump or ensuring model_validate can ignore extras.
    # For now, we assume model_validate will work if all required model fields are present.

    review_model_instance = models.SystematicReview.model_validate(review_model_data)
    review_model_instance.id = review_id_val

    db_session.add(review_model_instance)
    db_session.commit()
    db_session.refresh(review_model_instance)
    return review_model_instance


@pytest.mark.integration
def test_search_pubmed_and_store_results(
    search_service_test_session_factory: sessionmaker[Session],
    test_review_with_criteria: models.SystematicReview,
):
    """Verify that search_pubmed_and_store_results fetches, maps, and stores data correctly."""
    # Arrange
    review_id = test_review_with_criteria.id
    assert review_id is not None, "Test review should have an ID"

    search_service = services.SearchService(factory=search_service_test_session_factory)

    # Using a very specific query for a known small set of results is better for stability.
    # This query might still vary over time or return 0 results.
    test_query = "(diabetes mellitus type 2[MeSH Major Topic]) AND (metformin[Title/Abstract]) AND (randomized controlled trial[Publication Type])"
    max_results_to_fetch = 1  # Keep very low for stable integration test

    # Act
    stored_results: Sequence[models.SearchResult] = []
    try:
        # NCBI_EMAIL and NCBI_API_KEY must be set in .env.test (loaded by conftest.py db_engine fixture)
        stored_results = search_service.search_pubmed_and_store_results(
            review_id=review_id,
            query=test_query,
            max_results=max_results_to_fetch,
            # Credentials are now handled by the service from env vars
        )
    except services.ServiceError as e:
        if "NCBI_EMAIL environment variable not set" in str(e):
            pytest.skip(
                "NCBI_EMAIL not set in environment, skipping PubMed API call test."
            )
        else:
            pytest.fail(
                f"SearchService.search_pubmed_and_store_results call failed: {e}"
            )
    except Exception as e:
        pytest.fail(f"An unexpected error occurred during PubMed search: {e}")

    # Assert
    assert isinstance(stored_results, list), "Service should return a list of results"
    assert len(stored_results) <= max_results_to_fetch, (
        f"Expected at most {max_results_to_fetch} results, got {len(stored_results)}"
    )

    # Verification of DB state using a new session from the same factory
    # The service uses factory.begin() which handles commit/rollback
    with search_service_test_session_factory.begin() as verification_session:
        search_repo_for_verification = repositories.SearchResultRepository()
        results_in_db = search_repo_for_verification.get_by_review_id(
            verification_session, review_id
        )

    assert len(results_in_db) == len(stored_results), (
        f"Mismatch in result count: service returned {len(stored_results)}, DB query found {len(results_in_db)}"
    )

    if not stored_results:
        print(
            f"Warning: PubMed query '{test_query}' returned no results. Detailed assertions for stored data will be skipped."
        )

    if stored_results:
        assert len(results_in_db) > 0, (
            "Should have stored results in DB if service returned any"
        )
        for result in results_in_db:
            assert isinstance(result, models.SearchResult)
            assert result.review_id == review_id
            assert result.source_db == models.SearchDatabaseSource.PUBMED
            assert result.source_id is not None
            assert result.title is not None
            assert result.year is None or isinstance(result.year, str), (
                f"Year should be str or None, got {type(result.year)}"
            )
            assert result.raw_data, (
                "raw_data (cleaned record) should not be empty"
            )  # Check if dict is not empty
            assert isinstance(result.raw_data, dict), "raw_data should be a dict"
            assert "mesh_headings" in result.source_metadata, (
                "MeSH headings should be in source_metadata"
            )
            if result.source_metadata.get(
                "mesh_headings"
            ):  # Check if mesh_headings list is not empty if present
                assert isinstance(result.source_metadata["mesh_headings"], list)
