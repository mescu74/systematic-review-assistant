"""Integration tests for the Service layer."""

from __future__ import annotations

import uuid  # Needed for test_review_with_criteria
from datetime import UTC, datetime

import pytest
from pytest_mock import MockerFixture
from sqlalchemy.orm import sessionmaker
from sqlmodel import (
    Session,  # NOTE: This is the SQLModel Session, NEVER use SQLAlchemy Session!
)

# Imports used by test and fixtures
from sr_assistant.app import services
from sr_assistant.core import models, repositories, schemas, types
from sr_assistant.core.types import CriteriaFramework


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
            review_id=review_id, query=test_query, max_results=max_results_to_fetch
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
            review_id=review_id, query=test_query, max_results=max_results_to_fetch
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


# ============================================================================
# ScreeningService Integration Tests (AC8)
# ============================================================================


@pytest.fixture(scope="function")
def test_review_with_search_results(
    db_session: Session,
) -> tuple[models.SystematicReview, list[models.SearchResult]]:
    """Creates a test review with search results for screening tests."""
    # Create review
    review_id = uuid.uuid4()
    review = models.SystematicReview(
        id=review_id,
        research_question="Test screening integration",
        background="Test background",
        inclusion_criteria="Include diabetes studies",
        exclusion_criteria="Exclude non-English studies",
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": "Adults with diabetes",
            "intervention": "Metformin",
            "comparison": "Placebo",
            "outcome": "HbA1c reduction",
        },
    )
    db_session.add(review)
    db_session.commit()  # Commit review first
    db_session.refresh(review)

    # Create search results after review is committed
    search_results = []
    for i in range(3):
        sr = models.SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            title=f"Test Study {i + 1}: Metformin Effects",
            abstract=f"This study examines metformin effects in diabetes patients. Study {i + 1} shows significant results.",
            source_db=types.SearchDatabaseSource.PUBMED,
            source_id=f"pmid{12345678 + i}",
            year="2023",
            authors=["Smith J", "Doe A"],
            journal="Diabetes Care",
        )
        search_results.append(sr)
        db_session.add(sr)

    db_session.commit()  # Commit search results
    for sr in search_results:
        db_session.refresh(sr)

    return review, search_results


@pytest.fixture(scope="function")
def screening_service_test_session_factory(
    db_session: Session,
) -> sessionmaker[Session]:
    """Creates a sessionmaker bound to the integration test DB engine for ScreeningService."""
    test_engine = db_session.get_bind()
    return sessionmaker(bind=test_engine, class_=Session, expire_on_commit=False)


@pytest.mark.integration
def test_screening_service_perform_batch_abstract_screening_success(
    screening_service_test_session_factory: sessionmaker[Session],
    test_review_with_search_results: tuple[
        models.SystematicReview, list[models.SearchResult]
    ],
    mocker: MockerFixture,
):
    """Test ScreeningService.perform_batch_abstract_screening with mocked LLM calls."""
    review, search_results = test_review_with_search_results
    search_result_ids = [sr.id for sr in search_results]

    # Mock the screening agent to avoid actual LLM calls
    from datetime import UTC, datetime

    from sr_assistant.app.agents.screening_agents import ScreenAbstractResultTuple
    from sr_assistant.core import schemas

    # Create mock screening results
    mock_screening_results = []
    for i, sr in enumerate(search_results):
        # Conservative result
        conservative_result = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=review.id,
            search_result_id=sr.id,
            decision=types.ScreeningDecisionType.INCLUDE
            if i % 2 == 0
            else types.ScreeningDecisionType.EXCLUDE,
            confidence_score=0.85,
            rationale=f"Conservative screening rationale for study {i + 1}",
            screening_strategy=types.ScreeningStrategyType.CONSERVATIVE,
            model_name="gpt-4",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
        )

        # Comprehensive result
        comprehensive_result = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=review.id,
            search_result_id=sr.id,
            decision=types.ScreeningDecisionType.INCLUDE
            if i % 2 == 0
            else types.ScreeningDecisionType.EXCLUDE,
            confidence_score=0.90,
            rationale=f"Comprehensive screening rationale for study {i + 1}",
            screening_strategy=types.ScreeningStrategyType.COMPREHENSIVE,
            model_name="gpt-4",
            start_time=datetime.now(UTC),
            end_time=datetime.now(UTC),
            trace_id=uuid.uuid4(),
        )

        # Add exclusion reasons for excluded studies
        if i % 2 == 1:
            conservative_result.exclusion_reason_categories = schemas.ExclusionReasons(
                population_exclusion_reasons=["Wrong age range for inclusion criteria"]
            )
            comprehensive_result.exclusion_reason_categories = schemas.ExclusionReasons(
                population_exclusion_reasons=["Wrong age range for inclusion criteria"]
            )

        mock_screening_results.append(
            ScreenAbstractResultTuple(sr, conservative_result, comprehensive_result)
        )

    # Create mock callback handler
    from langchain_community.callbacks import OpenAICallbackHandler

    from sr_assistant.app.agents.screening_agents import ScreenAbstractsBatchOutput

    mock_cb = OpenAICallbackHandler()
    mock_cb.total_cost = 0.05
    mock_cb.total_tokens = 1000

    # Create proper ScreenAbstractsBatchOutput
    mock_batch_output = ScreenAbstractsBatchOutput(
        results=mock_screening_results, cb=mock_cb
    )

    # Mock the screen_abstracts_batch function
    mock_screen_abstracts_batch = mocker.patch(
        "sr_assistant.app.services.screen_abstracts_batch",
        return_value=mock_batch_output,
    )

    # Create ScreeningService instance
    screening_service = services.ScreeningService(
        factory=screening_service_test_session_factory
    )

    # Act
    result_tuples = screening_service.perform_batch_abstract_screening(
        review_id=review.id, search_result_ids_to_screen=search_result_ids
    )

    # Assert service call
    assert len(result_tuples) == 3
    mock_screen_abstracts_batch.assert_called_once()

    # Verify database persistence
    with screening_service_test_session_factory.begin() as verification_session:
        # Check that ScreenAbstractResult records were created
        screen_repo = repositories.ScreenAbstractResultRepository()
        screen_results = screen_repo.get_by_review_id(verification_session, review.id)

        # Should have 6 results (2 strategies × 3 search results)
        assert len(screen_results) == 6

        # Verify conservative and comprehensive results exist
        conservative_results = [
            r
            for r in screen_results
            if r.screening_strategy == types.ScreeningStrategyType.CONSERVATIVE
        ]
        comprehensive_results = [
            r
            for r in screen_results
            if r.screening_strategy == types.ScreeningStrategyType.COMPREHENSIVE
        ]

        assert len(conservative_results) == 3
        assert len(comprehensive_results) == 3

        # Check that SearchResult instances were updated with result IDs
        search_repo = repositories.SearchResultRepository()
        updated_search_results = search_repo.get_by_review_id(
            verification_session, review.id
        )

        for updated_sr in updated_search_results:
            assert updated_sr.conservative_result_id is not None
            assert updated_sr.comprehensive_result_id is not None

            # Verify the linkage is correct
            conservative_result = next(
                (
                    r
                    for r in conservative_results
                    if r.id == updated_sr.conservative_result_id
                ),
                None,
            )
            comprehensive_result = next(
                (
                    r
                    for r in comprehensive_results
                    if r.id == updated_sr.comprehensive_result_id
                ),
                None,
            )

            assert conservative_result is not None
            assert comprehensive_result is not None
            # NOTE: ScreenAbstractResult doesn't have search_result_id field, it's linked via foreign key
            # The linkage is verified by the fact that the IDs match the SearchResult foreign key references


@pytest.mark.integration
def test_screening_service_perform_batch_abstract_screening_with_errors(
    screening_service_test_session_factory: sessionmaker[Session],
    test_review_with_search_results: tuple[
        models.SystematicReview, list[models.SearchResult]
    ],
    mocker: MockerFixture,
):
    """Test ScreeningService handles screening errors properly."""
    review, search_results = test_review_with_search_results
    search_result_ids = [sr.id for sr in search_results]

    from langchain_community.callbacks import OpenAICallbackHandler

    from sr_assistant.app.agents.screening_agents import (
        ScreenAbstractResultTuple,
        ScreenAbstractsBatchOutput,
        ScreeningError,
    )

    # Create mixed results with some errors
    mock_screening_results = []
    for i, sr in enumerate(search_results):
        if i == 1:  # Second result has an error
            error = ScreeningError(
                search_result=sr,
                error="LLM timeout",
                message="The screening agent timed out",
            )
            mock_screening_results.append(ScreenAbstractResultTuple(sr, error, error))
        else:
            # Normal successful results
            conservative_result = schemas.ScreeningResult(
                id=uuid.uuid4(),
                review_id=review.id,
                search_result_id=sr.id,
                decision=types.ScreeningDecisionType.INCLUDE,
                confidence_score=0.85,
                rationale=f"Rationale for study {i + 1}",
                screening_strategy=types.ScreeningStrategyType.CONSERVATIVE,
                model_name="gpt-4",
                start_time=datetime.now(UTC),
                end_time=datetime.now(UTC),
                trace_id=uuid.uuid4(),
            )

            comprehensive_result = schemas.ScreeningResult(
                id=uuid.uuid4(),
                review_id=review.id,
                search_result_id=sr.id,
                decision=types.ScreeningDecisionType.INCLUDE,
                confidence_score=0.90,
                rationale=f"Rationale for study {i + 1}",
                screening_strategy=types.ScreeningStrategyType.COMPREHENSIVE,
                model_name="gpt-4",
                start_time=datetime.now(UTC),
                end_time=datetime.now(UTC),
                trace_id=uuid.uuid4(),
            )

            mock_screening_results.append(
                ScreenAbstractResultTuple(sr, conservative_result, comprehensive_result)
            )

    # Create mock callback handler and batch output
    mock_cb = OpenAICallbackHandler()
    mock_cb.total_cost = 0.03
    mock_cb.total_tokens = 800

    mock_batch_output = ScreenAbstractsBatchOutput(
        results=mock_screening_results, cb=mock_cb
    )

    # Mock the screen_abstracts_batch function
    mocker.patch(
        "sr_assistant.app.services.screen_abstracts_batch",
        return_value=mock_batch_output,
    )

    # Create ScreeningService instance
    screening_service = services.ScreeningService(
        factory=screening_service_test_session_factory
    )

    # Act
    result_tuples = screening_service.perform_batch_abstract_screening(
        review_id=review.id, search_result_ids_to_screen=search_result_ids
    )

    # Assert
    assert len(result_tuples) == 3

    # Verify database state - only successful results should be persisted
    with screening_service_test_session_factory.begin() as verification_session:
        screen_repo = repositories.ScreenAbstractResultRepository()
        screen_results = screen_repo.get_by_review_id(verification_session, review.id)

        # Should have 4 results (2 strategies × 2 successful search results)
        assert len(screen_results) == 4

        # Check SearchResult updates - only successful ones should have result IDs
        search_repo = repositories.SearchResultRepository()
        updated_search_results = search_repo.get_by_review_id(
            verification_session, review.id
        )

        successful_count = 0
        error_count = 0

        for updated_sr in updated_search_results:
            if (
                updated_sr.conservative_result_id is not None
                and updated_sr.comprehensive_result_id is not None
            ):
                successful_count += 1
            else:
                error_count += 1
                # The errored search result should not have been updated
                assert updated_sr.conservative_result_id is None
                assert updated_sr.comprehensive_result_id is None

        assert successful_count == 2  # First and third search results
        assert error_count == 1  # Second search result had errors


@pytest.mark.integration
def test_screening_service_review_not_found(
    screening_service_test_session_factory: sessionmaker[Session],
):
    """Test ScreeningService raises RecordNotFoundError for non-existent review."""
    screening_service = services.ScreeningService(
        factory=screening_service_test_session_factory
    )

    non_existent_review_id = uuid.uuid4()
    search_result_ids = [uuid.uuid4()]

    with pytest.raises(
        repositories.RecordNotFoundError, match="SystematicReview with ID .* not found"
    ):
        screening_service.perform_batch_abstract_screening(
            review_id=non_existent_review_id,
            search_result_ids_to_screen=search_result_ids,
        )


@pytest.mark.integration
def test_screening_service_no_search_results_found(
    screening_service_test_session_factory: sessionmaker[Session],
    test_review_with_search_results: tuple[
        models.SystematicReview, list[models.SearchResult]
    ],
):
    """Test ScreeningService handles case where no search results are found for given IDs."""
    review, _ = test_review_with_search_results

    screening_service = services.ScreeningService(
        factory=screening_service_test_session_factory
    )

    # Use non-existent search result IDs
    non_existent_ids = [uuid.uuid4(), uuid.uuid4()]

    # Should not raise an error, but return empty results
    result_tuples = screening_service.perform_batch_abstract_screening(
        review_id=review.id, search_result_ids_to_screen=non_existent_ids
    )

    assert len(result_tuples) == 0

    # Verify no database changes
    with screening_service_test_session_factory.begin() as verification_session:
        screen_repo = repositories.ScreenAbstractResultRepository()
        screen_results = screen_repo.get_by_review_id(verification_session, review.id)
        assert len(screen_results) == 0
