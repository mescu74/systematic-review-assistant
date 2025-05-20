import csv
import uuid
from pathlib import Path
import typing as t 
import html # Added for XML escaping
import sys # <--- Add this import

from loguru import logger
from sqlmodel import Session # Removed select
from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from sqlalchemy import delete as sqlalchemy_delete # <--- Add this import
from psycopg.errors import UniqueViolation # Import specific DB error

from sr_assistant.app.database import session_factory
from sr_assistant.core import models 
from sr_assistant.core.types import CriteriaFramework, SearchDatabaseSource

# Fixed UUID for the benchmark review
BENCHMARK_REVIEW_ID = uuid.UUID("00000000-1111-2222-3333-444444444444")
BENCHMARK_PROTOCOL_PATH = Path("docs/benchmark/bechmark-protocol.md") # Not used for PICO/exclusion per story
BENCHMARK_CSV_PATH = Path("docs/benchmark/human-reviewer-results-to-bench-against.csv")

def parse_protocol_and_create_review(db_session: Session) -> models.SystematicReview:
    """Uses the user-provided PICO and exclusion criteria to create or update a SystematicReview model."""
    logger.info(f"Attempting to create/update benchmark SystematicReview object with ID: {BENCHMARK_REVIEW_ID}")

    # --- Restored PICO Elements from Story 4.1 / previous correct script version ---
    pico_population_text = (
        "Individuals experiencing homelessness. Studies must include data collected in the "
        "Republic of Ireland. The review focuses on the health of homeless individuals "
        "themselves, not on reports from key informants about their needs. Studies that "
        "only include international/European datasets without specific outcomes for the "
        "Republic of Ireland should be excluded."
    )
    pico_intervention_text = ( # Note: Story 4.1 refers to this as 'Phenomenon of Interest' for non-intervention studies
        "The review focuses on studies that generate empirical data (quantitative or "
        "qualitative) on the following health-related topics for the homeless population:\n"
        "- Overall health status.\n"
        "- Health care access, utilisation, and quality.\n"
        "- Specific health conditions (e.g., addiction, diabetes, cancer, "
        "communicable/non-communicable diseases, STIs, pregnancy and childbirth).\n"
        "- Health behaviours (e.g., nutrition, child development, tobacco use, "
        "vaccination).\n"
        "- Social determinants of health (e.g., social and community context, education, "
        "economic stability)."
    )
    pico_comparison_text = (
        "Studies that include a comparison/control group comprising the general, housed "
        "population are of interest.\n"
        "Such studies should contain a method for comparing health indicator(s) between the "
        "homeless (exposed) group and the general housed (control) group (e.g., using "
        "relative risk, absolute difference, slope/relative index of inequality)."
    )
    pico_outcome_text = (
        "- Empirical indicators of health status.\n"
        "- Empirical indicators of health care access.\n"
        "- Empirical indicators of health care quality.\n"
        "- Empirical indicators of health care utilisation."
    )

    inclusion_criteria_xml_parts = [
        f"<population>{html.escape(pico_population_text)}</population>",
        f"<intervention>{html.escape(pico_intervention_text)}</intervention>", # Or PhenomenonOfInterest
        f"<comparison>{html.escape(pico_comparison_text)}</comparison>",
        f"<outcome>{html.escape(pico_outcome_text)}</outcome>",
    ]
    inclusion_criteria_str = "\n".join(inclusion_criteria_xml_parts)

    # --- Restored Exclusion Criteria from Story 4.1 / previous correct script version ---
    exclusion_criteria_str = ("""Population-related:
- Studies focusing on key informants reporting on the needs of the homeless, rather than the homeless population directly.
- Studies with no data from the Republic of Ireland.
- Studies using international/European datasets that include data from Ireland but do not report outcomes specific to the Republic of Ireland.
Study Design & Publication Type-related:
- Studies that do not generate empirical primary or secondary data on a health topic. This includes:
    - Modelling studies.
    - Commentaries/Letters.
    - Individual case reports.
- Conference Abstracts.
- Policy papers.
- Guidelines.
- Grey literature (e.g., government documents and reports, pre-print articles, research reports, statistical reports) - due to resource limitations for a thorough search.
Topic-related:
- Animal studies.
- Studies on economic, health care, or housing policy that do not relate directly to health outcomes or access for the homeless population.
Language-related:
- Studies not published in English.
Date-related:
- Studies published before January 1, 2012.
Full-Text Screening Specific Exclusions:
- Studies that do not contain empirical health indicators for the general, housed population (when a comparison is implied or attempted).
- Studies that do not provide a method for comparing health indicator(s) between the exposed (homeless) and control (general housed) groups (e.g., missing denominators for calculating rates or risks).""").strip()

    review_data_dict = {
        "title": "Benchmark Systematic Review on Homelessness and Healthcare Access in Dublin", # Kept concise title
        "research_question": "Benchmark: What is the health status, healthcare access/utilization/quality, and what are the health conditions, health behaviours, and social determinants of health for individuals experiencing homelessness in the Republic of Ireland, and how do these compare to the general housed population where data allows?",
        "background":("This is a benchmark systematic review protocol based on the pre-defined PICO and exclusion criteria "
                     "from Story 4.1 for Epic 4 (SRA Benchmarking Module Implementation). "
                     "It uses the dataset from 'docs/benchmark/human-reviewer-results-to-bench-against.csv' "
                     "and is intended for testing and evaluation purposes of the SRA screening module."),
        "criteria_framework": CriteriaFramework.PICO, # Set the enum member
        "criteria_framework_answers": { # Directly create the dict
            "population": pico_population_text,
            "intervention": pico_intervention_text,
            "comparison": pico_comparison_text,
            "outcome": pico_outcome_text,
        },
        "inclusion_criteria": inclusion_criteria_str,
        "exclusion_criteria": exclusion_criteria_str, # Already stripped
        "review_metadata": {
            "benchmark_source_protocol_version": "Story 4.1 - Pre-defined PICO and Exclusion Criteria",
            "benchmark_data_source_csv": str(BENCHMARK_CSV_PATH)
        },
    }

    review_to_save: models.SystematicReview | None = None 
    existing_review = db_session.get(models.SystematicReview, BENCHMARK_REVIEW_ID)

    if existing_review:
        logger.info(f"SystematicReview with ID {BENCHMARK_REVIEW_ID} already exists. Updating.")
        # Directly pass review_data_dict to sqlmodel_update
        existing_review.sqlmodel_update(review_data_dict) # New: Pass the original dict
        review_to_save = existing_review
    else:
        logger.info(f"Creating new SystematicReview with ID {BENCHMARK_REVIEW_ID}.")
        review_to_save = models.SystematicReview.model_validate({"id": BENCHMARK_REVIEW_ID, **review_data_dict})

    if not review_to_save:
        logger.error("Review object to save is unexpectedly None. Aborting review creation/update.")
        # This condition should ideally not be hit if the logic above is sound.
        # If it is, it indicates a flaw in how existing_review or new review_to_save is assigned.
        raise ValueError("Failed to prepare review object for saving due to unexpected None.")

    db_session.add(review_to_save)
    try:
        db_session.commit()
        db_session.refresh(review_to_save)
        logger.info(f"SystematicReview {BENCHMARK_REVIEW_ID} successfully saved/updated.")
        return review_to_save
    except IntegrityError as e:
        db_session.rollback()
        logger.warning(f"IntegrityError during commit for review {BENCHMARK_REVIEW_ID}: {e}. Attempting to re-fetch and update.")
        if isinstance(e.orig, UniqueViolation) and "systematic_reviews_pkey" in str(e.orig).lower():
            logger.info(f"UniqueViolation for primary key. Review {BENCHMARK_REVIEW_ID} was likely created by another process. Re-fetching and updating.")
            freshly_fetched_review = db_session.get(models.SystematicReview, BENCHMARK_REVIEW_ID)
            if freshly_fetched_review:
                logger.info(f"Successfully re-fetched review {BENCHMARK_REVIEW_ID}. Proceeding with update.")
                for key, value in review_data_dict.items():
                    setattr(freshly_fetched_review, key, value)
                db_session.add(freshly_fetched_review)
                db_session.commit()
                db_session.refresh(freshly_fetched_review)
                logger.info(f"SystematicReview {BENCHMARK_REVIEW_ID} successfully updated after re-fetch.")
                return freshly_fetched_review
            else:
                logger.error(f"Failed to re-fetch review {BENCHMARK_REVIEW_ID} after UniqueViolation. This is unexpected.")
            raise 
        else:
            logger.error(f"Unhandled IntegrityError for review {BENCHMARK_REVIEW_ID}. Original error: {e.orig}")
            raise 
    except Exception as e:
        db_session.rollback()
        logger.exception(f"Unexpected error during commit for review {BENCHMARK_REVIEW_ID}: {e}")
        raise

def infer_source_db(key: str) -> SearchDatabaseSource:
    """Infers the source database from the benchmark key."""
    if key.startswith("rayyan-"):
        return SearchDatabaseSource.OTHER 
    elif key.startswith("pmid_") or key.isdigit(): 
        return SearchDatabaseSource.PUBMED
    # Add more specific checks if CSV 'key' or other fields (like 'url') indicate Scopus/Embase
    # For example, if URL contains 'scopus.com' or 'embase.com'
    # Placeholder for now, defaulting to OTHER for non-PubMed like keys if not Rayyan
    return SearchDatabaseSource.OTHER

def parse_csv_and_create_search_results(review_id: uuid.UUID) -> list[models.SearchResult]:
    """Parses the benchmark CSV and creates SearchResult models."""
    logger.info(f"Parsing benchmark CSV data from {BENCHMARK_CSV_PATH}")
    parsed_search_results: list[models.SearchResult] = []
    current_raw_data: dict[str, t.Any] = {} 
    rows_yielded_by_reader = 0 # Counter for rows from DictReader
    
    try:
        with open(BENCHMARK_CSV_PATH, mode='r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            logger.info("CSV DictReader initialized.")
            for row_num, row_dict in enumerate(reader):
                rows_yielded_by_reader += 1
                logger.debug(f"Processing CSV file row number (0-indexed from data): {row_num}, Total rows from reader so far: {rows_yielded_by_reader}")
                current_raw_data = dict(row_dict) 
                try:
                    key = current_raw_data.get("key", f"benchmark_row_{row_num+2}")
                    logger.debug(f"Attempting to parse row with key: {key}")
                    source_db_enum = infer_source_db(key)
                    source_id_val = key 

                    authors_list = [author.strip() for author in current_raw_data.get("authors", "").split(';') if author.strip()] if current_raw_data.get("authors") else []
                    keywords_list = [kw.strip() for kw in current_raw_data.get("keywords", "").split(';') if kw.strip()] if current_raw_data.get("keywords") else []
                    
                    included_r2_val = current_raw_data.get("included_round2", "").strip().upper()
                    human_decision: bool | None = None
                    if included_r2_val == 'Y':
                        human_decision = True
                    elif included_r2_val == 'N':
                        human_decision = False
                    else:
                        included_r1_val = current_raw_data.get("included_round1", "").strip().upper()
                        if included_r1_val == 'Y':
                            human_decision = True
                        elif included_r1_val == 'N':
                            human_decision = False
                        else:
                            logger.warning(f"Row {row_num+2}: Could not determine benchmark decision for key {key}. Defaulting to None.")

                    result = models.SearchResult(
                        review_id=review_id,
                        source_db=source_db_enum,
                        source_id=source_id_val,
                        doi=current_raw_data.get("doi") if current_raw_data.get("doi") else None,
                        title=current_raw_data.get("title", "N/A"),
                        abstract=current_raw_data.get("abstract") if current_raw_data.get("abstract") else None,
                        journal=current_raw_data.get("journal") if current_raw_data.get("journal") else None,
                        year=str(current_raw_data.get("year")) if current_raw_data.get("year") else None,
                        authors=authors_list,
                        keywords=keywords_list,
                        raw_data=current_raw_data, 
                        source_metadata={
                            "benchmark_human_decision": human_decision, 
                            "original_key": key, 
                            "included_round1": current_raw_data.get("included_round1", ""),
                            "included_round2": current_raw_data.get("included_round2", ""),
                            "exclusion_stage_round1": current_raw_data.get("exclusion_stage_round1", "")
                        }
                    )
                    parsed_search_results.append(result)
                except Exception:
                    log_key = current_raw_data.get('key', f'row_{row_num+2}') 
                    logger.exception(f"Error processing CSV row {row_num+2} (key: {log_key!r})") 
                    logger.debug(f"Problematic row data: {current_raw_data!r}")
                    continue 
            logger.info(f"Finished iterating through DictReader. Total rows yielded: {rows_yielded_by_reader}")
                    
    except FileNotFoundError:
        logger.exception(f"Benchmark CSV file not found at {BENCHMARK_CSV_PATH}") 
        return [] 
    except Exception:
        logger.exception(f"Outer exception: Failed to read or parse CSV {BENCHMARK_CSV_PATH}") 
        return [] 
        
    logger.info(f"Successfully parsed {len(parsed_search_results)} search results from CSV. Expected around 585 for full benchmark CSV.")
    return parsed_search_results

def seed_data_to_db(db: Session, review: models.SystematicReview, results: list[models.SearchResult]) -> None:
    """Seeds the review and search results into the database."""
    try:
        # The 'review' object is assumed to be committed and refreshed by parse_protocol_and_create_review.
        
        # Step 1: Delete old SearchResult records for this review using a bulk delete.
        logger.info(f"Deleting existing SearchResult records for review ID: {review.id} using bulk delete.")
        delete_stmt = sqlalchemy_delete(models.SearchResult).where(
            models.SearchResult.review_id == review.id # pyright: ignore
        )
        result_proxy = db.execute(delete_stmt) # pyright: ignore[reportDeprecated] # Use SQLAlchemy's execute for broader compatibility
        db.commit()
        deleted_count = result_proxy.rowcount # pyright: ignore[reportUnknownMemberType, reportAttributeAccessIssue]
        logger.info(f"Bulk deleted {deleted_count} old SearchResult records.")

        # Step 2: Add new SearchResult records.
        if results:
            logger.info(f"Adding {len(results)} new benchmark search results to database for review ID: {review.id}.")
            for res in results: 
                res.review_id = review.id # Ensure review_id is correctly set on each result
            db.add_all(results)
            db.commit()
            logger.info(f"{len(results)} SearchResult records successfully added.")
        else:
            logger.info("No new SearchResult records to add.")

        logger.info(f"Benchmark data successfully seeded for review ID: {review.id}.")

    except SQLAlchemyError: # Catch specific SQLAlchemy errors
        db.rollback()
        logger.exception("Database error during seeding") 
        raise
    except Exception:
        db.rollback()
        logger.exception("Unexpected error during seeding") 
        raise

if __name__ == "__main__":
    # Setup logger to also print to console for direct script runs
    logger.remove()
    logger.add(sys.stderr, level="DEBUG") # Print DEBUG and above to console
    logger.add("logs/seed_benchmark_data.log", rotation="5 MB", level="INFO") # INFO and above to file
    logger.info("Starting benchmark data seeding process...")

    with session_factory() as main_session: # Renamed from 'session' to 'main_session' for clarity
        try:
            logger.info("Parsing benchmark protocol to create/update SystematicReview object...")
            benchmark_review_model = parse_protocol_and_create_review(main_session) # Pass main_session
            
            logger.info("Parsing benchmark CSV to create SearchResult objects...")
            benchmark_search_results = parse_csv_and_create_search_results(benchmark_review_model.id) # Use the ID from the returned model

            logger.info("Seeding data to the database...")
            seed_data_to_db(main_session, benchmark_review_model, benchmark_search_results) # Pass main_session

            logger.info("Benchmark data seeding process completed successfully.")
        except Exception as e:
            logger.exception("An error occurred during the seeding process. Rolling back changes.")
            if 'main_session' in locals() and main_session.is_active:
                main_session.rollback()
            sys.exit(1) # Exit with error code 