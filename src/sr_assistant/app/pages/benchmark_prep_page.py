# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Human Benchmark Page - Display benchmark protocol details."""

from __future__ import annotations

import typing as t
import uuid

import streamlit as st
from loguru import logger

from sr_assistant.app.services import ReviewService
from sr_assistant.core import models

# Import the benchmark review ID from the seeding tool
try:
    from tools.seed_benchmark_data import BENCHMARK_REVIEW_ID as _benchmark_review_id
except ImportError:
    # Fallback if tools module not available in deployment
    _benchmark_review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")

if t.TYPE_CHECKING:
    from collections.abc import MutableMapping

# Constants
MIN_PICO_ELEMENTS_FOR_READINESS = 3
READINESS_SCORE_PER_SECTION = 25
READY_THRESHOLD = 75
NEEDS_ATTENTION_THRESHOLD = 50


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
