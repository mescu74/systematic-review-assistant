import uuid
from pathlib import Path
import sys # <--- Add this import

import pandas as pd  # Add pandas for Excel reading
from loguru import logger
from sqlmodel import Session # Removed select
from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from sqlalchemy import delete as sqlalchemy_delete # <--- Add this import
from psycopg.errors import UniqueViolation # Import specific DB error
from sqlalchemy import text # Added for SQL text queries

from sr_assistant.app.database import session_factory
from sr_assistant.core import models 
from sr_assistant.core.types import SearchDatabaseSource

# Fixed UUID for the benchmark review
BENCHMARK_REVIEW_ID = uuid.UUID("00000000-1111-2222-3333-444444444444")
BENCHMARK_PROTOCOL_PATH = Path("docs/benchmark/bechmark-protocol.md") # Not used for PICO/exclusion per story
BENCHMARK_EXCEL_PATH = Path("docs/benchmark/benchmark_human_ground_truth.xlsx")  # Changed from CSV to Excel

def parse_protocol_and_create_review(db_session: Session) -> models.SystematicReview:
    """Uses the user-provided PICO and exclusion criteria to create or update a SystematicReview model."""
    logger.info(f"Attempting to create/update benchmark SystematicReview object with ID: {BENCHMARK_REVIEW_ID}")

    inclusion_criteria_xml_parts = """
  <InclusionCriteria>
    <Population>
      <Item>Homeless</Item>
      <Item>Key informants reporting on needs of Homeless</Item>
    </Population>
    <Setting>
      <Item>Data collected in Republic of Ireland</Item>
    </Setting>
    <StudyDesign>
      <Item>Empirical primary or secondary data on a health topic</Item>
      <Item>Quantitative studies</Item>
      <Item>Qualitative studies</Item>
    </StudyDesign>
    <PublicationType>
      <Item>Peer-reviewed publications</Item>
      <Item>Conference Abstracts</Item>
      <Item>Systematic reviews examined for individual studies meeting inclusion criteria</Item>
    </PublicationType>
    <Topic>
      <Item>Health conditions (e.g., addiction, diabetes, cancer, communicable/non-communicable disease, STI, pregnancy and childbirth, etc.)</Item>
      <Item>Health behaviours (e.g., nutrition, child development, tobacco use, vaccination, etc.)</Item>
      <Item>Health care access, utilisation, quality</Item>
      <Item>Social determinants of health (e.g., social and community context, education, economic stability)</Item>
    </Topic>
    <Language>English</Language>
    <Date>Published in 2012 or later</Date>
  </InclusionCriteria>"""

    # --- Restored Exclusion Criteria from Story 4.1 / previous correct script version ---
    exclusion_criteria_str = """\
  <ExclusionCriteria>
    <Population></Population>
    <Setting>
      <Item>No data from the Republic of Ireland</Item>
      <Item>Studies using international/European datasets without specific outcomes for Republic of Ireland</Item>
    </Setting>
    <StudyDesign>
      <Item>No empirical primary or secondary data on a health topic</Item>
      <Item>Modelling studies</Item>
      <Item>Commentaries/Letters</Item>
      <Item>Individual case reports</Item>
    </StudyDesign>
    <PublicationType>
      <Item>Policy papers</Item>
      <Item>Guidelines</Item>
      <Item>Systematic reviews containing studies meeting inclusion criteria</Item>
      <Item>Grey literature (government documents and reports, pre-print articles, research reports, statistical reports)</Item>
    </PublicationType>
    <Topic>
      <Item>Animal study</Item>
      <Item>Economic/health care/housing policy not relating to health</Item>
    </Topic>
    <Language>Any language that is not English</Language>
    <Date>Published before 2012</Date>
  </ExclusionCriteria>"""
    review_data_dict = {
        "title": "Benchmark Systematic Review on Homelessness and Healthcare Access in Dublin", # Kept concise title
        "research_question": "Benchmark: What is the health status, healthcare access/utilization/quality, and what are the health conditions, health behaviours, and social determinants of health for individuals experiencing homelessness in the Republic of Ireland, and how do these compare to the general housed population where data allows?",
        "background":("This is a benchmark systematic review protocol based on the pre-defined PICO and exclusion criteria "
                     "from Story 4.1 for Epic 4 (SRA Benchmarking Module Implementation). "
                     "It uses the dataset from 'docs/benchmark/human-reviewer-results-to-bench-against.csv' "
                     "and is intended for testing and evaluation purposes of the SRA screening module."),
        "inclusion_criteria": inclusion_criteria_xml_parts,
        "exclusion_criteria": exclusion_criteria_str, # Already stripped
        "review_metadata": {
            "benchmark_source_protocol_version": "Story 4.1 - Pre-defined PICO and Exclusion Criteria",
            "benchmark_data_source_csv": str(BENCHMARK_EXCEL_PATH)
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

def parse_excel_and_create_search_results(review_id: uuid.UUID) -> list[models.SearchResult]:
    """Parses the benchmark Excel file and creates SearchResult models."""
    logger.info(f"Parsing benchmark Excel data from {BENCHMARK_EXCEL_PATH}")
    parsed_search_results: list[models.SearchResult] = []
    
    try:
        # Read Excel file using pandas
        df = pd.read_excel(BENCHMARK_EXCEL_PATH)
        logger.info(f"Excel file loaded. Found {len(df)} rows.")
        
        # Process each row
        for idx, row in df.iterrows():
            # idx is index from pandas, cast to int to avoid type issues
            if isinstance(idx, int):
                row_idx = idx
            else:
                # If pandas changes and idx is not int, convert it
                row_idx = int(str(idx))  # Convert to string first then int for safety
            
            # Check if key is NaN - if so, we've reached summary rows, stop processing
            key_val = row.get("key", row.get(df.columns[0], pd.NA))
            if pd.isna(key_val):
                logger.info(f"Reached end of data rows at index {row_idx} (summary row)")
                break
                
            key = str(key_val) if not pd.isna(key_val) else f"benchmark_row_{row_idx+2}"
            logger.debug(f"Processing row {row_idx} with key: {key}")
            
            try:
                source_db_enum = infer_source_db(key)
                source_id_val = key
                
                # Parse authors and keywords
                authors_str = str(row.get("authors", "")) if pd.notna(row.get("authors")) else ""
                authors_list = [author.strip() for author in authors_str.split(';') if author.strip()]
                
                keywords_str = str(row.get("keywords", "")) if pd.notna(row.get("keywords")) else ""
                keywords_list = [kw.strip() for kw in keywords_str.split(';') if kw.strip()]
                
                # Determine ground truth based on "Included after T&A screen" column
                included_ta_screen = str(row.get("Included after T&A screen", "")) if pd.notna(row.get("Included after T&A screen")) else ""
                
                # If "Included after T&A screen" is "Y", they were included (passed title/abstract screening)
                # Otherwise (including "N" or any other value), they were excluded
                human_decision = True if included_ta_screen == "Y" else False
                
                logger.debug(f"Row {row_idx}: Included after T&A screen='{included_ta_screen}', human_decision={human_decision}")
                
                # Clean up year - strip .0 from float values
                year_val = row.get("year")
                if pd.notna(year_val):
                    year_str = str(year_val)
                    # Remove .0 from float values like 2022.0
                    if year_str.endswith('.0'):
                        year_str = year_str[:-2]
                else:
                    year_str = None
                
                # Convert row to dict and replace NaN values with None for JSON compatibility
                raw_data_dict = {}
                for col, val in row.items():
                    if pd.isna(val):
                        raw_data_dict[col] = None
                    else:
                        raw_data_dict[col] = val
                
                # Create the result object
                result = models.SearchResult(
                    review_id=review_id,
                    source_db=source_db_enum,
                    source_id=source_id_val,
                    doi=str(row.get("doi")) if pd.notna(row.get("doi")) else None,
                    title=str(row.get("title", "N/A")) if pd.notna(row.get("title")) else "N/A",
                    abstract=str(row.get("abstract")) if pd.notna(row.get("abstract")) else None,
                    journal=str(row.get("journal")) if pd.notna(row.get("journal")) else None,
                    year=year_str,
                    authors=authors_list,
                    keywords=keywords_list,
                    raw_data=raw_data_dict,  # Use the cleaned dict
                    source_metadata={
                        "benchmark_human_decision": human_decision,
                        "original_key": key,
                        "included_after_ta_screen": included_ta_screen
                    }
                )
                parsed_search_results.append(result)
                
            except Exception:
                logger.exception(f"Error processing Excel row {row_idx} (key: {key!r})")
                continue
                
    except FileNotFoundError:
        logger.exception(f"Benchmark Excel file not found at {BENCHMARK_EXCEL_PATH}")
        return []
    except Exception:
        logger.exception(f"Failed to read or parse Excel file {BENCHMARK_EXCEL_PATH}")
        return []
    
    # Log summary statistics
    total_included = sum(1 for r in parsed_search_results if r.source_metadata.get("benchmark_human_decision") is True)
    total_excluded = sum(1 for r in parsed_search_results if r.source_metadata.get("benchmark_human_decision") is False)
    logger.info(f"Successfully parsed {len(parsed_search_results)} search results from Excel.")
    logger.info(f"Ground truth summary: {total_included} included (passed title/abstract), {total_excluded} excluded at title/abstract")
    
    return parsed_search_results

def seed_data_to_db(db: Session, review: models.SystematicReview, results: list[models.SearchResult]) -> None:
    """Seeds the review and search results into the database."""
    try:
        # The 'review' object is assumed to be committed and refreshed by parse_protocol_and_create_review.
        
        # Step 1: Delete related benchmark_result_items first to handle foreign key constraints
        logger.info(f"Deleting existing BenchmarkResultItem records for review ID: {review.id}")
        # Delete benchmark_result_items that reference search_results for this review
        delete_items_stmt = text("""
            DELETE FROM benchmark_result_items 
            WHERE search_result_id IN (
                SELECT id FROM search_results WHERE review_id = :review_id
            )
        """)
        db.execute(delete_items_stmt, {"review_id": review.id})  # pyright: ignore[reportDeprecated]
        db.commit()
        logger.info("Deleted related BenchmarkResultItem records.")
        
        # Step 2: Delete old SearchResult records for this review using a bulk delete.
        logger.info(f"Deleting existing SearchResult records for review ID: {review.id} using bulk delete.")
        delete_stmt = sqlalchemy_delete(models.SearchResult).where(
            models.SearchResult.review_id == review.id # pyright: ignore
        )
        result_proxy = db.execute(delete_stmt) # pyright: ignore[reportDeprecated] # Use SQLAlchemy's execute for broader compatibility
        db.commit()
        deleted_count = result_proxy.rowcount # pyright: ignore[reportUnknownMemberType, reportAttributeAccessIssue]
        logger.info(f"Bulk deleted {deleted_count} old SearchResult records.")

        # Step 3: Add new SearchResult records.
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
            
            logger.info("Parsing benchmark Excel to create SearchResult objects...")
            benchmark_search_results = parse_excel_and_create_search_results(benchmark_review_model.id) # Use the ID from the returned model

            logger.info("Seeding data to the database...")
            seed_data_to_db(main_session, benchmark_review_model, benchmark_search_results) # Pass main_session

            logger.info("Benchmark data seeding process completed successfully.")
        except Exception as e:
            logger.exception("An error occurred during the seeding process. Rolling back changes.")
            if 'main_session' in locals() and main_session.is_active:
                main_session.rollback()
            sys.exit(1) # Exit with error code 
