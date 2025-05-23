# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Human Benchmark Page - Display benchmark protocol details."""

from __future__ import annotations

import typing as t
import uuid

import streamlit as st
from loguru import logger

from sr_assistant.app.agents.screening_agents import (
    invoke_resolver_chain,
    screen_abstracts_batch,
)
from sr_assistant.app.database import session_factory
from sr_assistant.app.services import ReviewService
from sr_assistant.core import models, schemas
from sr_assistant.core.repositories import (
    BenchmarkResultItemRepository,
    BenchmarkRunRepository,
    SearchResultRepository,
)
from sr_assistant.core.types import ScreeningDecisionType

if t.TYPE_CHECKING:
    from collections.abc import MutableMapping

    from sr_assistant.core import models

# Constants
MIN_PICO_ELEMENTS_FOR_READINESS = 3
READINESS_SCORE_PER_SECTION = 25
READY_THRESHOLD = 75
NEEDS_ATTENTION_THRESHOLD = 50

# Import the benchmark review ID from the seeding tool
try:
    from tools.seed_benchmark_data import BENCHMARK_REVIEW_ID

    _benchmark_review_id = BENCHMARK_REVIEW_ID
except ImportError:
    # Fallback if tools module not available in deployment
    _benchmark_review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")


def display_protocol_overview(benchmark_review: models.SystematicReview) -> None:
    """Display the protocol overview in a prominent card format."""
    st.header("🎯 Protocol Overview")

    # Create an attractive overview card
    with st.container():
        st.markdown("### Research Question")
        if benchmark_review.research_question:
            st.info(benchmark_review.research_question)
        else:
            st.warning("No research question available.")

        st.markdown("### Background")
        if benchmark_review.background:
            with st.expander("📖 View Background", expanded=False):
                st.markdown(benchmark_review.background)
        else:
            st.warning("No background information available.")


def display_pico_elements(
    criteria_framework_answers: MutableMapping[str, t.Any] | None,
) -> None:
    """Display PICO elements in an attractive card layout."""
    st.header("🧬 PICO Framework")

    if not criteria_framework_answers:
        st.warning("No PICO elements available.")
        return

    # Create a 2x2 grid for PICO elements
    col1, col2 = st.columns(2)

    pico_config = [
        ("population", "👥 Population", col1),
        ("intervention", "💊 Intervention", col2),
        ("comparison", "⚖️ Comparison", col1),
        ("outcome", "📊 Outcome", col2),
    ]

    for key, label, column in pico_config:
        with column:
            if key in criteria_framework_answers:
                value = criteria_framework_answers[key]
                if value and str(value).strip():
                    st.markdown(f"**{label}**")
                    st.markdown(f"> {value}")
                else:
                    st.markdown(f"**{label}**")
                    st.caption("_Not specified_")
            else:
                st.markdown(f"**{label}**")
                st.caption("_Not specified_")
            st.markdown("---")


def _categorize_exclusion_criteria(exclusion_lines: list[str]) -> dict[str, list[str]]:
    """Categorize exclusion criteria into groups."""
    categories = {
        "language_related": [],
        "study_type_related": [],
        "publication_related": [],
        "other_criteria": [],
    }

    for line in exclusion_lines:
        line_lower = line.lower()
        if any(word in line_lower for word in ["language", "english", "non-english"]):
            categories["language_related"].append(line)
        elif any(
            word in line_lower
            for word in ["study", "design", "case report", "review", "meta-analysis"]
        ):
            categories["study_type_related"].append(line)
        elif any(
            word in line_lower
            for word in ["publication", "published", "peer", "preprint"]
        ):
            categories["publication_related"].append(line)
        else:
            categories["other_criteria"].append(line)

    return categories


def _display_categorized_exclusion_criteria(categories: dict[str, list[str]]) -> None:
    """Display categorized exclusion criteria."""
    category_config = [
        ("language_related", "🌐 Language Requirements:"),
        ("study_type_related", "🔬 Study Design:"),
        ("publication_related", "📰 Publication Type:"),
        ("other_criteria", "📝 Other Criteria:"),
    ]

    for key, title in category_config:
        if categories[key]:
            st.markdown(f"**{title}**")
            for item in categories[key]:
                st.markdown(f"• {item}")


def _display_inclusion_criteria(benchmark_review: models.SystematicReview) -> None:
    """Display inclusion criteria section."""
    st.subheader("✅ Inclusion Criteria")
    if benchmark_review.inclusion_criteria:
        with st.expander("View XML Criteria (for LLM processing)", expanded=False):
            st.code(benchmark_review.inclusion_criteria, language="xml")

        # Try to extract readable criteria from XML
        st.markdown("**Criteria Summary:**")
        criteria_terms = ["population", "intervention", "comparison", "outcome"]
        for term in criteria_terms:
            if term in benchmark_review.inclusion_criteria.lower():
                st.markdown(f"• {term.title()} criteria defined")
    else:
        st.warning("No inclusion criteria available.")


def _display_exclusion_criteria(benchmark_review: models.SystematicReview) -> None:
    """Display exclusion criteria section."""
    st.subheader("❌ Exclusion Criteria")
    if benchmark_review.exclusion_criteria:
        exclusion_lines = [
            line.strip()
            for line in benchmark_review.exclusion_criteria.split("\n")
            if line.strip()
        ]

        if exclusion_lines:
            categories = _categorize_exclusion_criteria(exclusion_lines)
            _display_categorized_exclusion_criteria(categories)
        else:
            st.markdown(benchmark_review.exclusion_criteria)
    else:
        st.warning("No exclusion criteria available.")


def display_criteria_sections(benchmark_review: models.SystematicReview) -> None:
    """Display inclusion and exclusion criteria in organized sections."""
    st.header("📋 Study Selection Criteria")

    # Create two columns for inclusion and exclusion criteria
    col1, col2 = st.columns([1, 1])

    with col1:
        _display_inclusion_criteria(benchmark_review)

    with col2:
        _display_exclusion_criteria(benchmark_review)


def _calculate_readiness_score(
    benchmark_review: models.SystematicReview,
) -> tuple[int, list[str]]:
    """Calculate protocol readiness score and factors."""
    readiness_score = 0
    readiness_factors = []

    # Check research question
    if benchmark_review.research_question:
        readiness_score += READINESS_SCORE_PER_SECTION
        readiness_factors.append("✅ Research question defined")
    else:
        readiness_factors.append("❌ Research question missing")

    # Check background
    if benchmark_review.background:
        readiness_score += READINESS_SCORE_PER_SECTION
        readiness_factors.append("✅ Background provided")
    else:
        readiness_factors.append("❌ Background missing")

    # Check PICO elements
    if (
        benchmark_review.criteria_framework_answers
        and len(benchmark_review.criteria_framework_answers)
        >= MIN_PICO_ELEMENTS_FOR_READINESS
    ):
        readiness_score += READINESS_SCORE_PER_SECTION
        readiness_factors.append("✅ PICO elements defined")
    else:
        readiness_factors.append("❌ Incomplete PICO elements")

    # Check selection criteria
    if benchmark_review.inclusion_criteria and benchmark_review.exclusion_criteria:
        readiness_score += READINESS_SCORE_PER_SECTION
        readiness_factors.append("✅ Selection criteria complete")
    else:
        readiness_factors.append("❌ Selection criteria incomplete")

    return readiness_score, readiness_factors


def _display_readiness_indicator(
    readiness_score: int, readiness_factors: list[str]
) -> None:
    """Display protocol readiness indicator."""
    col1, col2 = st.columns([1, 2])

    with col1:
        if readiness_score >= READY_THRESHOLD:
            st.success(f"**Protocol Readiness: {readiness_score}%**")
            st.success("🚀 Ready for benchmarking!")
        elif readiness_score >= NEEDS_ATTENTION_THRESHOLD:
            st.warning(f"**Protocol Readiness: {readiness_score}%**")
            st.warning("⚠️ Needs attention")
        else:
            st.error(f"**Protocol Readiness: {readiness_score}%**")
            st.error("❌ Not ready")

    with col2:
        st.markdown("**Readiness Checklist:**")
        for factor in readiness_factors:
            st.markdown(f"  {factor}")


def display_protocol_summary(benchmark_review: models.SystematicReview) -> None:
    """Display the protocol summary in an attractive dashboard format."""
    st.header("📊 Protocol Summary Dashboard")

    # Create metrics in a more prominent layout
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            label="📋 Review ID",
            value=str(benchmark_review.id)[:8] + "...",  # Shortened for display
            help=f"Full ID: {benchmark_review.id}",
        )

    with col2:
        framework_value = (
            benchmark_review.criteria_framework.value
            if benchmark_review.criteria_framework
            else "Not specified"
        )
        st.metric(label="🧬 Framework", value=framework_value)

    with col3:
        pico_count = (
            len(benchmark_review.criteria_framework_answers)
            if benchmark_review.criteria_framework_answers
            else 0
        )
        st.metric(label="🎯 PICO Elements", value=f"{pico_count}/4")

    with col4:
        criteria_status = (
            "Complete ✅"
            if benchmark_review.inclusion_criteria
            and benchmark_review.exclusion_criteria
            else "Incomplete ⚠️"
        )
        st.metric(label="📋 Criteria Status", value=criteria_status)

    # Add protocol readiness indicator
    st.markdown("---")
    readiness_score, readiness_factors = _calculate_readiness_score(benchmark_review)
    _display_readiness_indicator(readiness_score, readiness_factors)


def _needs_resolver(
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
) -> bool:
    """Determine if resolver is needed based on conservative and comprehensive results.

    Args:
        conservative_result: Result from conservative screening agent
        comprehensive_result: Result from comprehensive screening agent

    Returns:
        True if resolver should be invoked, False otherwise
    """
    # Check for disagreement between conservative and comprehensive
    if conservative_result.decision != comprehensive_result.decision:
        return True

    # Check if both are uncertain
    if (
        conservative_result.decision == ScreeningDecisionType.UNCERTAIN
        and comprehensive_result.decision == ScreeningDecisionType.UNCERTAIN
    ):
        return True

    # Check for low confidence even with agreement
    if (
        conservative_result.confidence_score < 0.8
        or comprehensive_result.confidence_score < 0.8
    ):
        return True

    return False


def _determine_final_decision(
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
    resolver_result: schemas.ResolverOutputSchema | None = None,
) -> ScreeningDecisionType:
    """Determine the final decision based on all available results.

    Args:
        conservative_result: Result from conservative screening agent
        comprehensive_result: Result from comprehensive screening agent
        resolver_result: Result from resolver agent (if invoked)

    Returns:
        The final screening decision
    """
    if resolver_result is not None:
        return resolver_result.resolver_decision

    # If no resolver was invoked, use the agreed-upon decision
    if conservative_result.decision == comprehensive_result.decision:
        return conservative_result.decision

    # This shouldn't happen if _needs_resolver logic is correct
    msg = f"No resolver result but decisions disagree: conservative={conservative_result.decision}, comprehensive={comprehensive_result.decision}"
    logger.warning(msg)
    return ScreeningDecisionType.UNCERTAIN


def _calculate_classification(
    final_decision: ScreeningDecisionType,
    human_decision: bool | None,
) -> str:
    """Calculate the classification (TP, FP, TN, FN) based on final decision vs human decision.

    Args:
        final_decision: The SRA's final decision
        human_decision: The human ground truth decision (True=include, False=exclude, None=unknown)

    Returns:
        Classification string: "TP", "FP", "TN", "FN", or "UNKNOWN"
    """
    if human_decision is None:
        return "UNKNOWN"

    # Convert final_decision to boolean for comparison
    ai_include = final_decision == ScreeningDecisionType.INCLUDE
    human_include = human_decision

    if ai_include and human_include:
        return "TP"  # True Positive
    if ai_include and not human_include:
        return "FP"  # False Positive
    if not ai_include and not human_include:
        return "TN"  # True Negative
    if not ai_include and human_include:
        return "FN"  # False Negative
    return "UNKNOWN"


def run_benchmark_screening(benchmark_review: models.SystematicReview) -> None:
    """Execute the benchmark screening process.

    Args:
        benchmark_review: The systematic review protocol to use for benchmarking
    """
    # Initialize repositories
    benchmark_run_repo = BenchmarkRunRepository()
    benchmark_result_item_repo = BenchmarkResultItemRepository()
    search_result_repo = SearchResultRepository()

    # Create progress indicators
    progress_bar = st.progress(0)
    status_text = st.empty()

    try:
        with session_factory() as session:
            # Step 1: Create BenchmarkRun record
            status_text.text("Creating benchmark run record...")

            benchmark_run = models.BenchmarkRun(
                review_id=benchmark_review.id,
                config_details={
                    "conservative_model": "gemini-2.5-pro-preview-05-06",
                    "comprehensive_model": "gemini-2.5-pro-preview-05-06",
                    "resolver_model": "gemini-2.5-pro-preview-05-06",
                    "run_type": "benchmark_screening",
                    "timestamp": str(uuid.uuid1()),  # Include timestamp for uniqueness
                },
                run_notes="Automated benchmark run from UI",
            )

            created_run = benchmark_run_repo.add(session, benchmark_run)
            session.commit()

            logger.info(f"Created benchmark run {created_run.id}")

            # Step 2: Fetch all SearchResults for the benchmark review
            status_text.text("Fetching benchmark dataset...")
            search_results = search_result_repo.get_by_review_id(
                session, benchmark_review.id
            )

            if not search_results:
                st.error("No search results found for the benchmark review.")
                return

            total_items = len(search_results)
            logger.info(f"Processing {total_items} search results for benchmark")

            # Step 3: Process each SearchResult
            processed_items = 0

            for i, search_result in enumerate(search_results):
                progress = (i + 1) / total_items
                progress_bar.progress(progress)
                title_preview = (
                    search_result.title[:50] if search_result.title else "Untitled"
                )
                status_text.text(
                    f"Processing item {i + 1} of {total_items}: {title_preview}..."
                )

                try:
                    # Step 3a: Call screen_abstracts_batch
                    batch_output = screen_abstracts_batch(
                        batch=[search_result], batch_idx=0, review=benchmark_review
                    )

                    if batch_output is None or not batch_output.results:
                        logger.error(
                            f"No screening results for search result {search_result.id}"
                        )
                        continue

                    result_tuple = batch_output.results[0]
                    conservative_result = result_tuple.conservative_result
                    comprehensive_result = result_tuple.comprehensive_result

                    # Check if results are errors and skip if so
                    if not isinstance(conservative_result, schemas.ScreeningResult):
                        logger.error(
                            f"Conservative screening error for search result {search_result.id}: {conservative_result!r}"
                        )
                        continue

                    if not isinstance(comprehensive_result, schemas.ScreeningResult):
                        logger.error(
                            f"Comprehensive screening error for search result {search_result.id}: {comprehensive_result!r}"
                        )
                        continue

                    # Step 3b: Check if resolver is needed
                    resolver_result = None
                    if _needs_resolver(conservative_result, comprehensive_result):
                        status_text.text(f"Resolving conflicts for item {i + 1}...")
                        resolver_result = invoke_resolver_chain(
                            search_result=search_result,
                            review=benchmark_review,
                            conservative_result=conservative_result,
                            comprehensive_result=comprehensive_result,
                        )

                    # Step 3c: Determine final decision
                    final_decision = _determine_final_decision(
                        conservative_result, comprehensive_result, resolver_result
                    )

                    # Step 3d: Get human decision from source_metadata
                    human_decision_raw = search_result.source_metadata.get(
                        "benchmark_human_decision"
                    )
                    human_decision: bool | None = None
                    if isinstance(human_decision_raw, bool):
                        human_decision = human_decision_raw
                    elif human_decision_raw is not None:
                        # Try to convert string representations to bool
                        if str(human_decision_raw).lower() in (
                            "true",
                            "1",
                            "yes",
                            "include",
                        ):
                            human_decision = True
                        elif str(human_decision_raw).lower() in (
                            "false",
                            "0",
                            "no",
                            "exclude",
                        ):
                            human_decision = False

                    # Step 3e: Calculate classification
                    classification = _calculate_classification(
                        final_decision, human_decision
                    )

                    # Step 3f: Create BenchmarkResultItem
                    result_item = models.BenchmarkResultItem(
                        benchmark_run_id=created_run.id,
                        search_result_id=search_result.id,
                        human_decision=human_decision,
                        conservative_decision=conservative_result.decision,
                        conservative_confidence=conservative_result.confidence_score,
                        conservative_rationale=conservative_result.rationale,
                        comprehensive_decision=comprehensive_result.decision,
                        comprehensive_confidence=comprehensive_result.confidence_score,
                        comprehensive_rationale=comprehensive_result.rationale,
                        resolver_decision=resolver_result.resolver_decision
                        if resolver_result
                        else None,
                        resolver_confidence=resolver_result.resolver_confidence_score
                        if resolver_result
                        else None,
                        resolver_reasoning=resolver_result.resolver_reasoning
                        if resolver_result
                        else None,
                        final_decision=final_decision,
                        classification=classification,
                    )

                    benchmark_result_item_repo.add(session, result_item)
                    processed_items += 1

                except Exception as exc:
                    logger.exception(
                        f"Error processing search result {search_result.id}: {exc!r}"
                    )
                    continue

            # Commit all result items
            session.commit()

            # Update progress to complete
            progress_bar.progress(1.0)
            status_text.text(
                f"Benchmark run completed! Processed {processed_items} of {total_items} items."
            )

            # Show success message
            st.success("✅ Benchmark run completed successfully!")
            st.info(f"**Run ID:** {created_run.id}")
            st.info(f"**Items processed:** {processed_items} of {total_items}")

            logger.info(
                f"Benchmark run {created_run.id} completed with {processed_items} items processed"
            )

    except Exception as exc:
        logger.exception(f"Error during benchmark run: {exc!r}")
        st.error(f"❌ Benchmark run failed: {exc}")
        raise


def main() -> None:
    """Main function for the human benchmark page."""
    # Main title with description
    st.title("🎯 Benchmark Protocol")
    st.markdown(
        """
        <div style='background-color: #f0f2f6; padding: 1rem; border-radius: 0.5rem; margin-bottom: 2rem;'>
            <p style='margin: 0; color: #262730;'>
                <strong>📖 Purpose:</strong> Review the systematic review protocol that will be used for benchmarking.
                This protocol defines the research question, selection criteria, and framework for evaluating AI screening performance.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Initialize ReviewService
    review_service = ReviewService()

    try:
        # Fetch the benchmark review
        benchmark_review: models.SystematicReview | None = review_service.get_review(
            review_id=_benchmark_review_id
        )

        if benchmark_review is None:
            st.error(
                f"🚫 **Benchmark protocol not found** (ID: `{_benchmark_review_id}`)"
            )
            st.info(
                """
                **📥 To seed the benchmark data:**
                ```bash
                uv run python tools/seed_benchmark_data.py
                ```
                """
            )
            return

        # Display protocol sections with improved layout
        display_protocol_overview(benchmark_review)
        st.markdown("---")

        display_pico_elements(benchmark_review.criteria_framework_answers)
        st.markdown("---")

        display_criteria_sections(benchmark_review)
        st.markdown("---")

        display_protocol_summary(benchmark_review)

        # Add benchmark execution section
        st.markdown("---")
        st.header("🚀 Benchmark Execution")

        st.markdown(
            """
            <div style='background-color: #e8f4fd; padding: 1rem; border-radius: 0.5rem; margin-bottom: 1rem;'>
                <p style='margin: 0; color: #1f4e79;'>
                    <strong>⚡ Run Benchmark:</strong> Execute the AI screening pipeline on the entire benchmark dataset.
                    This will process all search results using conservative and comprehensive agents, with resolver intervention when needed.
                </p>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Benchmark execution button
        if st.button(
            "🎯 Run New Benchmark Screening", type="primary", use_container_width=True
        ):
            st.markdown("### 🔄 Benchmark Execution in Progress")
            run_benchmark_screening(benchmark_review)

        # Add raw protocol data view
        st.markdown("---")
        with st.expander("🔍 Raw Protocol Data (JSON)", expanded=False):
            st.caption(
                "Technical view of the complete protocol record for debugging and verification purposes."
            )
            st.json(benchmark_review.model_dump(mode="json"), expanded=False)

        # Success footer
        st.success("✅ Protocol loaded successfully!")

        logger.info(
            f"Successfully displayed benchmark protocol for review {benchmark_review.id}"
        )

    except Exception:
        logger.exception("Error loading benchmark protocol")
        st.error(
            """
            🚨 **An error occurred while loading the benchmark protocol.**

            Please check the application logs for more details.
            """
        )


if __name__ == "__main__":
    main()
