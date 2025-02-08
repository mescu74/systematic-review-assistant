"""Generate test data for systematic review database using SQLModel."""
from datetime import UTC, datetime
import uuid

from sqlmodel import Session, SQLModel, create_engine, select

from sr_assistant.core.models import (
    SystematicReview,
    PubMedResult,
    ScreenAbstractResultModel,
)
from sr_assistant.core.types import ScreeningDecisionType, ScreeningStrategyType

def create_test_data(db_url: str) -> None:
    """Create test data in database.

    Creates reviews, PubMed results, and screening results with proper relationships.
    Uses SQLModel for clean, type-safe database operations.
    """
    engine = create_engine(db_url, echo=True)
    SQLModel.metadata.create_all(engine)

    with Session(engine) as session:
        # Create reviews
        reviews = [
            SystematicReview(
                id=uuid.uuid4(),
                background="Test background for systematic review on cancer treatment",
                research_question="What is the efficacy of immunotherapy in treating melanoma?",
                inclusion_criteria="RCTs studying immunotherapy in melanoma patients",
                exclusion_criteria="Non-RCT studies, non-melanoma cancers",
            ),
            SystematicReview(
                id=uuid.uuid4(),
                background="COVID-19 vaccine effectiveness study",
                research_question="How effective are mRNA vaccines against COVID-19 variants?",
                inclusion_criteria="Clinical trials of mRNA vaccines",
                exclusion_criteria="Non-clinical studies, non-mRNA vaccines",
            ),
        ]
        # Add and commit reviews
        session.add_all(reviews)
        session.commit()

        # Create PubMed results for each review
        pubmed_results = []
        for review in reviews:
            for i in range(5):
                result = PubMedResult(
                    id=uuid.uuid4(),
                    review_id=review.id,
                    query=f"Test query for {review.research_question}",
                    pmid=f"PMD{i+1}",
                    pmc=f"PMC00{i+1}",
                    doi=f"10.1000/test.{i+1}",
                    title=f"Test Article {i+1} for Review {review.id}",
                    abstract=f"This is a test abstract for article {i+1}",
                    journal=f"Test Journal {i+1}",
                    year=str(2020 + i),
                )
                pubmed_results.append(result)

        # Add and commit PubMed results
        session.add_all(pubmed_results)
        session.commit()

        # Create screening results for each PubMed result
        screening_results = []
        screening_pairs = []  # Store pairs for later FK update
        for result in pubmed_results:
            trace_id = uuid.uuid4()

            # Conservative screening
            conservative = ScreenAbstractResultModel(
                id=uuid.uuid4(),
                review_id=result.review_id,
                pubmed_result_id=result.id,
                trace_id=trace_id,
                decision=ScreeningDecisionType.INCLUDE if result.year > "2022" else ScreeningDecisionType.EXCLUDE,
                confidence_score=0.85,
                rationale="Test rationale for conservative screening",
                extracted_quotes=["Test quote 1", "Test quote 2"],
                screening_strategy=ScreeningStrategyType.CONSERVATIVE,
                model_name="test-model",
                start_time=datetime.now(UTC),
                end_time=datetime.now(UTC),
                response_metadata={"test": "metadata"},
            )
            screening_results.append(conservative)

            # Comprehensive screening
            comprehensive = ScreenAbstractResultModel(
                id=uuid.uuid4(),
                review_id=result.review_id,
                pubmed_result_id=result.id,
                trace_id=trace_id,
                decision=ScreeningDecisionType.INCLUDE,
                confidence_score=0.75,
                rationale="Test rationale for comprehensive screening",
                extracted_quotes=["Test quote 3", "Test quote 4"],
                screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
                model_name="test-model",
                start_time=datetime.now(UTC),
                end_time=datetime.now(UTC),
                response_metadata={"test": "metadata"},
            )
            screening_results.append(comprehensive)

            # Store screening result IDs in PubMed result
            #result.conservative_result_id = conservative.id
            #result.comprehensive_result_id = comprehensive.id

            screening_results.extend([conservative, comprehensive])
            screening_pairs.append((result, conservative.id, comprehensive.id))

        # Add screening results first
        session.add_all(screening_results)
        session.commit()

         # Now update PubMed results with the committed screening result IDs
        for result, cons_id, comp_id in screening_pairs:
            result.conservative_result_id = cons_id
            result.comprehensive_result_id = comp_id
        session.commit()

        # Verify results using SQLModel select
        stmt = select(ScreenAbstractResultModel).order_by(ScreenAbstractResultModel.created_at)
        results = session.exec(stmt).all()
        print(f"Created {len(results)} screening results")
        for result in results:
            print(f"- {result.screening_strategy}: {result.decision} ({result.confidence_score})")

if __name__ == "__main__":
    import os
    from dotenv import load_dotenv

    load_dotenv(".env.test")
    dburl = os.getenv("SRA_DATABASE_URL", "")
    if "sra_integration_test" not in dburl:
        print("Error: SRA_DATABASE_URL must be set to integration test database")
        exit(1)
    create_test_data(dburl)
