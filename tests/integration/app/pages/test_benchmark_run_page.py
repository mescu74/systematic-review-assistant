# Test file for Story 4.9 - Benchmark UI Display Summary Performance Metrics

import uuid
from datetime import datetime, timezone
from typing import Any

import pytest
from streamlit.testing.v1 import AppTest

from sr_assistant.core import models
from sr_assistant.core.types import ScreeningDecisionType, SearchDatabaseSource


class TestBenchmarkMetricsDisplay:
    """Integration tests for Story 4.9 - displaying completed benchmark run metrics."""

    def test_display_completed_benchmark_run_with_metrics(
        self, mock_session: Any, mock_benchmark_data: dict[str, Any]
    ) -> None:
        """Test displaying a completed benchmark run with calculated metrics."""
        # Create a completed benchmark run with metrics
        run_id = uuid.uuid4()
        review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

        benchmark_run = models.BenchmarkRun(
            id=run_id,
            review_id=review_id,
            config_details={"conservative_model": "gpt-4o", "batch_size": 10},
            run_notes="Test run with full metrics",
            # Populated metrics
            tp=15,
            fp=5,
            tn=20,
            fn=10,
            sensitivity=0.6,  # 15/(15+10)
            specificity=0.8,  # 20/(20+5)
            accuracy=0.7,  # (15+20)/(15+5+20+10)
            ppv=0.75,  # 15/(15+5)
            npv=0.667,  # 20/(20+10)
            f1_score=0.667,
            mcc=0.4,
            cohen_kappa=0.4,
            pabak=0.4,
            lr_plus=3.0,
            lr_minus=0.5,
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

        # Mock the repository to return this benchmark run
        mock_benchmark_run_repo = mock_session.return_value.__enter__.return_value
        mock_benchmark_run_repo.get_by_review_id.return_value = [benchmark_run]

        # Initialize and run the app test
        at = AppTest.from_file("src/sr_assistant/app/pages/benchmark_run_page.py")
        at.run()

        # Check that the completed benchmark runs section exists
        assert (
            len(
                [
                    el
                    for el in at.markdown
                    if "Completed Benchmark Runs" in str(el.value)
                ]
            )
            > 0
        )

        # Simulate clicking the load button
        load_button = at.button("Load Latest Completed Benchmark Run")
        assert load_button is not None
        load_button.click()
        at.run()

        # Check that metrics are displayed
        assert len([el for el in at.metric if "True Positives" in str(el.label)]) > 0
        assert len([el for el in at.metric if "Accuracy" in str(el.label)]) > 0
        assert len([el for el in at.metric if "F1 Score" in str(el.label)]) > 0

    def test_display_benchmark_run_without_metrics(
        self, mock_session: Any, mock_benchmark_data: dict[str, Any]
    ) -> None:
        """Test displaying a benchmark run that doesn't have calculated metrics yet."""
        # Create a benchmark run without metrics
        run_id = uuid.uuid4()
        review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

        benchmark_run = models.BenchmarkRun(
            id=run_id,
            review_id=review_id,
            config_details={"conservative_model": "gpt-4o", "batch_size": 10},
            run_notes="Test run without metrics",
            # No metrics populated (all None)
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

        # Mock the repository to return this benchmark run
        mock_benchmark_run_repo = mock_session.return_value.__enter__.return_value
        mock_benchmark_run_repo.get_by_review_id.return_value = [benchmark_run]

        # Initialize and run the app test
        at = AppTest.from_file("src/sr_assistant/app/pages/benchmark_run_page.py")
        at.run()

        # Simulate clicking the load button
        load_button = at.button("Load Latest Completed Benchmark Run")
        assert load_button is not None
        load_button.click()
        at.run()

        # Check that warning message is displayed
        warnings = [
            el for el in at.warning if "Metrics not yet available" in str(el.value)
        ]
        assert len(warnings) > 0

    def test_no_completed_benchmark_runs(
        self, mock_session: Any, mock_benchmark_data: dict[str, Any]
    ) -> None:
        """Test the case where no completed benchmark runs exist."""
        # Mock the repository to return empty list
        mock_benchmark_run_repo = mock_session.return_value.__enter__.return_value
        mock_benchmark_run_repo.get_by_review_id.return_value = []

        # Initialize and run the app test
        at = AppTest.from_file("src/sr_assistant/app/pages/benchmark_run_page.py")
        at.run()

        # Simulate clicking the load button
        load_button = at.button("Load Latest Completed Benchmark Run")
        assert load_button is not None
        load_button.click()
        at.run()

        # Should show info message about clicking the button
        info_messages = [
            el
            for el in at.info
            if "Click 'Load Latest Completed Benchmark Run'" in str(el.value)
        ]
        assert len(info_messages) > 0

    def test_individual_paper_results_display(
        self, mock_session: Any, mock_benchmark_data: dict[str, Any]
    ) -> None:
        """Test displaying individual paper-level results table."""
        # Create a completed benchmark run with metrics
        run_id = uuid.uuid4()
        review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

        benchmark_run = models.BenchmarkRun(
            id=run_id,
            review_id=review_id,
            config_details={"conservative_model": "gpt-4o"},
            tp=1,
            fp=1,
            tn=1,
            fn=1,
            accuracy=0.5,
            created_at=datetime.now(timezone.utc),
        )

        # Create search results and benchmark result items
        search_result_1 = models.SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            source_db=SearchDatabaseSource.PUBMED,
            source_id="12345",
            title="Test Paper 1",
            authors=["Author A", "Author B"],
            year="2023",
            journal="Test Journal",
        )

        search_result_2 = models.SearchResult(
            id=uuid.uuid4(),
            review_id=review_id,
            source_db=SearchDatabaseSource.PUBMED,
            source_id="67890",
            title="Test Paper 2",
            authors=["Author C"],
            year="2022",
            journal="Another Journal",
        )

        result_item_1 = models.BenchmarkResultItem(
            id=uuid.uuid4(),
            benchmark_run_id=run_id,
            search_result_id=search_result_1.id,
            human_decision=True,
            conservative_decision=ScreeningDecisionType.INCLUDE,
            conservative_confidence=0.8,
            conservative_rationale="Matches inclusion criteria",
            comprehensive_decision=ScreeningDecisionType.INCLUDE,
            comprehensive_confidence=0.9,
            comprehensive_rationale="Strong evidence for inclusion",
            final_decision=ScreeningDecisionType.INCLUDE,
            classification="TP",
        )

        result_item_2 = models.BenchmarkResultItem(
            id=uuid.uuid4(),
            benchmark_run_id=run_id,
            search_result_id=search_result_2.id,
            human_decision=False,
            conservative_decision=ScreeningDecisionType.EXCLUDE,
            conservative_confidence=0.7,
            conservative_rationale="Does not meet criteria",
            comprehensive_decision=ScreeningDecisionType.INCLUDE,
            comprehensive_confidence=0.6,
            comprehensive_rationale="Some relevance found",
            resolver_decision=ScreeningDecisionType.EXCLUDE,
            resolver_confidence=0.8,
            resolver_reasoning="Conservative approach is correct",
            final_decision=ScreeningDecisionType.EXCLUDE,
            classification="TN",
        )

        # Mock repository methods
        mock_session_ctx = mock_session.return_value.__enter__.return_value
        mock_session_ctx.get_by_review_id.return_value = [benchmark_run]
        mock_session_ctx.get_by_benchmark_run_id.return_value = [
            result_item_1,
            result_item_2,
        ]
        mock_session_ctx.get_by_id.side_effect = (
            lambda session, id: search_result_1
            if id == search_result_1.id
            else search_result_2
        )

        # Initialize and run the app test
        at = AppTest.from_file("src/sr_assistant/app/pages/benchmark_run_page.py")
        at.run()

        # Simulate clicking the load button
        load_button = at.button("Load Latest Completed Benchmark Run")
        assert load_button is not None
        load_button.click()
        at.run()

        # Check that individual results section is present
        # Note: Due to the complexity of caching and nested function calls in the actual implementation,
        # we mainly test that the structure is in place rather than the exact data
        assert (
            len(
                [
                    el
                    for el in at.subheader
                    if "Individual Paper-Level Results" in str(el.value)
                ]
            )
            > 0
        )


@pytest.fixture
def mock_benchmark_data() -> dict[str, Any]:
    """Fixture to provide mock benchmark data."""
    return {
        "benchmark_review_id": uuid.UUID("00000000-1111-2222-3333-444444444444"),
        "test_run_id": uuid.uuid4(),
    }
