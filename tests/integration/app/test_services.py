"""Integration tests for the Service layer."""

from __future__ import annotations

import uuid  # Needed for test_review_with_criteria

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
        id=review_id_val,
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
        review_metadata={},
    )

    # SQLModel.model_validate needs to be called on the SQLModel class
    # And the SystematicReview SQLModel expects 'reearch_question'
    # We need to map SystematicReviewCreate to SystematicReview, handling the typo

    review_model_data = review_create_schema.model_dump()

    # Fields not in SystematicReview model (if any from SystematicReviewCreate that don't map)
    # would ideally be handled by excluding them from model_dump or ensuring model_validate can ignore extras.
    # For now, we assume model_validate will work if all required model fields are present.

    review_model_instance = models.SystematicReview.model_validate(review_model_data)

    db_session.add(review_model_instance)
    db_session.commit()
    db_session.refresh(review_model_instance)
    return review_model_instance


@pytest.mark.integration
def test_search_pubmed_and_store_results(
    search_service_test_session_factory: sessionmaker[Session],
    test_review_with_criteria: models.SystematicReview,
):
    """Verify that search_pubmed_and_store_results fetches, maps, stores, and handles duplicates."""
    # Arrange
    review_id = test_review_with_criteria.id
    assert review_id is not None, "Test review should have an ID"

    search_service = services.SearchService(factory=search_service_test_session_factory)

    test_query = "(diabetes mellitus type 2[MeSH Major Topic]) AND (metformin[Title/Abstract]) AND (randomized controlled trial[Publication Type])"
    max_results_to_fetch = 1

    # Act - First Call
    stored_results_first_call = []
    try:
        stored_results_first_call = search_service.search_pubmed_and_store_results(
            review_id=review_id,
            query=test_query,
            max_results=max_results_to_fetch,
        )
    except services.ServiceError as e:
        if "NCBI_EMAIL environment variable not set" in str(e):
            pytest.skip(
                "NCBI_EMAIL not set in environment, skipping PubMed API call test."
            )
        else:
            pytest.fail(
                f"SearchService.search_pubmed_and_store_results call (1st) failed: {e}"
            )
    except Exception as e:
        pytest.fail(f"An unexpected error occurred during 1st PubMed search: {e}")

    # Assert - First Call
    assert isinstance(stored_results_first_call, list), (
        "Service should return a list (1st call)"
    )
    assert len(stored_results_first_call) <= max_results_to_fetch, (
        f"Expected at most {max_results_to_fetch} results (1st call), got {len(stored_results_first_call)}"
    )

    for result in stored_results_first_call:
        assert isinstance(result, schemas.SearchResultRead), (
            "Result should be a SearchResultRead schema (1st call)"
        )
        assert result.review_id == review_id, (
            "Result should have correct review_id (1st call)"
        )
        assert result.source_db == schemas.SearchDatabaseSource.PUBMED, (
            "Source DB should be PubMed (1st call)"
        )

    # Verification of DB state after first call
    with search_service_test_session_factory.begin() as verification_session_1:
        search_repo_for_verification = repositories.SearchResultRepository()
        results_in_db_after_first_call = search_repo_for_verification.get_by_review_id(
            verification_session_1, review_id
        )

    count_after_first_call = len(results_in_db_after_first_call)
    assert count_after_first_call == len(stored_results_first_call), (
        f"Mismatch in result count after 1st call: service returned {len(stored_results_first_call)}, DB has {count_after_first_call}"
    )

    if not stored_results_first_call:
        print(
            f"Warning: PubMed query '{test_query}' returned no results on 1st call. Duplicate check might not be effective."
        )
        # Depending on strictness, one might pytest.skip here if no results means the test cannot proceed meaningfully for duplicates.

    # Act - Second Call (to test duplicate handling)
    stored_results_second_call = []
    try:
        # Calling again with the same parameters
        stored_results_second_call = search_service.search_pubmed_and_store_results(
            review_id=review_id,
            query=test_query,
            max_results=max_results_to_fetch,
        )
    except services.ServiceError as e:
        # Assuming NCBI creds were okay for the first call if we reached here
        pytest.fail(
            f"SearchService.search_pubmed_and_store_results call (2nd) failed: {e}"
        )
    except Exception as e:
        pytest.fail(f"An unexpected error occurred during 2nd PubMed search: {e}")

    # Assert - Second Call & Duplicate Check
    # If the first call successfully stored results, the second call processing the same PMIDs
    # should not identify them as new additions to be returned by the service in this list.
    if stored_results_first_call:  # If there were results from the first call
        assert len(stored_results_second_call) == 0, (
            f"Expected 0 results from second call when duplicates were processed, got {len(stored_results_second_call)}. "
            f"Service should ideally return only newly added items."
        )
    else:  # If first call had no results, second call's length is less constrained by duplication logic
        assert len(stored_results_second_call) <= max_results_to_fetch, (
            f"Expected at most {max_results_to_fetch} results (2nd call, 1st was empty), got {len(stored_results_second_call)}"
        )

    with search_service_test_session_factory.begin() as verification_session_2:
        results_in_db_after_second_call = search_repo_for_verification.get_by_review_id(
            verification_session_2, review_id
        )

    count_after_second_call = len(results_in_db_after_second_call)
    assert count_after_second_call == count_after_first_call, (
        f"Duplicate check failed: DB count changed after 2nd call. "
        f"Before: {count_after_first_call}, After: {count_after_second_call}. "
        f"Expected no new unique items matching the exact same query to be added."
    )

    # If results were found, verify their content (can reuse some logic from first call)
    if stored_results_first_call:  # Check based on first call's success
        # Ensure the items in the DB after the second call are consistent with the first successful fetch
        assert len(results_in_db_after_second_call) == len(
            results_in_db_after_first_call
        )
        if results_in_db_after_first_call:  # if there was something to compare
            db_result_original = results_in_db_after_first_call[0]
            db_result_after_dupe_check = next(
                (
                    r
                    for r in results_in_db_after_second_call
                    if r.id == db_result_original.id
                ),
                None,
            )
            assert db_result_after_dupe_check is not None, (
                "Original result not found after duplicate call."
            )
            assert db_result_original.title == db_result_after_dupe_check.title
            assert db_result_original.source_id == db_result_after_dupe_check.source_id
    else:
        print(
            f"Warning: PubMed query '{test_query}' returned no results. Detailed assertions for stored data will be skipped."
        )
