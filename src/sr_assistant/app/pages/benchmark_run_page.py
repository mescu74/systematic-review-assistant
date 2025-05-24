# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Benchmark Tool Page - Execute AI screening benchmarks with advanced batching and metrics."""

from __future__ import annotations

import io
import typing as t
import uuid
from datetime import datetime, timezone
from statistics import mean, median, stdev

import pandas as pd
import streamlit as st
from loguru import logger
from sklearn.metrics import (
    accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
)

from sr_assistant.app.agents.screening_agents import (
    ScreeningError,
    invoke_resolver_chain,
    screen_abstracts_batch,
)
from sr_assistant.app.database import session_factory
from sr_assistant.core import models, schemas
from sr_assistant.core.repositories import (
    BenchmarkResultItemRepository,
    BenchmarkRunRepository,
    SearchResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.core.types import ScreeningDecisionType

if t.TYPE_CHECKING:
    import collections.abc

# Import the benchmark review ID from the seeding tool
try:
    from tools.seed_benchmark_data import BENCHMARK_REVIEW_ID
except ImportError:
    # Fallback if tools module not available in deployment
    _benchmark_review_id = uuid.UUID("00000000-1111-2222-3333-444444444444")
    BENCHMARK_REVIEW_ID = _benchmark_review_id  # pyright: ignore[reportConstantRedefinition]
    # Fallback if not available - we'll implement a simplified version


def _needs_resolver(
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
) -> bool:
    """Determine if resolver is needed based on conservative and comprehensive results."""
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
    min_confidence = min(
        conservative_result.confidence_score, comprehensive_result.confidence_score
    )
    if min_confidence < 0.7:
        return True

    return False


def _determine_final_decision(
    conservative_result: schemas.ScreeningResult,
    comprehensive_result: schemas.ScreeningResult,
    resolver_result: schemas.ScreeningResult | None,
) -> ScreeningDecisionType:
    """Determine final decision based on all screening results."""
    # If resolver was invoked and gave a definitive decision, use it
    if resolver_result and resolver_result.decision != ScreeningDecisionType.UNCERTAIN:
        return resolver_result.decision

    # If resolver was invoked but still uncertain, return uncertain
    if resolver_result and resolver_result.decision == ScreeningDecisionType.UNCERTAIN:
        return ScreeningDecisionType.UNCERTAIN

    # If no resolver, check for agreement
    if conservative_result.decision == comprehensive_result.decision:
        return conservative_result.decision

    # Disagreement without resolver - should not happen with proper _needs_resolver logic
    logger.warning(
        f"No resolver result but decisions disagree: conservative={conservative_result.decision}, comprehensive={comprehensive_result.decision}"
    )
    return ScreeningDecisionType.UNCERTAIN


def _calculate_classification(
    final_decision: ScreeningDecisionType,
    human_decision: bool | None,
) -> str:
    """Calculate the classification (TP, FP, TN, FN) based on final decision vs human decision."""
    if human_decision is None:
        return "UNKNOWN"

    # Convert final_decision to boolean for comparison
    ai_include = final_decision == ScreeningDecisionType.INCLUDE
    human_include = human_decision

    if ai_include and human_include:
        return "TP"  # True Positive
    if not ai_include and not human_include:
        return "TN"  # True Negative
    if ai_include and not human_include:
        return "FP"  # False Positive
    if not ai_include and human_include:
        return "FN"  # False Negative
    return "UNKNOWN"


def calculate_metrics(
    y_true: list[bool | None],
    y_pred_decision: list[ScreeningDecisionType | None],
    y_pred_confidence: list[float | None],
) -> dict[str, t.Any]:
    """Calculate screening performance metrics."""
    y_pred_bool = []
    valid_indices = []

    for i, pred_decision in enumerate(y_pred_decision):
        if i < len(y_true) and y_true[i] is not None and pred_decision is not None:
            if pred_decision == ScreeningDecisionType.INCLUDE:
                y_pred_bool.append(True)
                valid_indices.append(i)
            elif pred_decision == ScreeningDecisionType.EXCLUDE:
                y_pred_bool.append(False)
                valid_indices.append(i)

    y_true_filtered = [y_true[i] for i in valid_indices]
    metrics: dict[str, t.Any] = {}
    tp, fp, fn, tn = 0, 0, 0, 0

    if not y_true_filtered or not y_pred_bool:
        logger.warning(
            "Not enough valid data points for CM (empty after filtering Nones/UNCERTAINs)."
        )
        metrics["error"] = "Not enough valid data points for Confusion Matrix."
        metrics["Total Compared"] = 0
        metrics["True Positives (TP)"] = 0
        metrics["True Negatives (TN)"] = 0
        metrics["False Positives (FP)"] = 0
        metrics["False Negatives (FN)"] = 0
        return metrics

    if len(y_true_filtered) != len(y_pred_bool):
        logger.error(
            f"Mismatch in length of y_true_filtered ({len(y_true_filtered)}) and y_pred_bool ({len(y_pred_bool)})"
        )
        metrics["error"] = (
            "Length mismatch in true and predicted values after filtering."
        )
        metrics["Total Compared"] = len(y_true_filtered)
        metrics["True Positives (TP)"] = 0
        metrics["True Negatives (TN)"] = 0
        metrics["False Positives (FP)"] = 0
        metrics["False Negatives (FN)"] = 0
        return metrics

    try:
        unique_true_classes = set(y_true_filtered)
        unique_pred_classes = set(y_pred_bool)

        if len(unique_true_classes) == 1 and len(unique_pred_classes) == 1:
            if unique_true_classes == unique_pred_classes:
                if True in unique_true_classes:
                    tp = len(y_true_filtered)
                    tn, fp, fn = 0, 0, 0
                else:
                    tn = len(y_true_filtered)
                    tp, fp, fn = 0, 0, 0
            elif True in unique_true_classes:
                fn = len(y_true_filtered)
                tp, tn, fp = 0, 0, 0
            else:
                fp = len(y_true_filtered)
                tp, tn, fn = 0, 0, 0
        else:
            cm_values = confusion_matrix(
                y_true_filtered, y_pred_bool, labels=[False, True]
            ).ravel()
            if len(cm_values) == 4:  # noqa: PLR2004
                tn, fp, fn, tp = map(int, cm_values)
            else:
                logger.error(
                    f"Confusion matrix did not produce 4 values: {cm_values!r}."
                )
                metrics["error"] = "Confusion matrix calculation failed (wrong shape)."
    except ValueError:
        logger.exception(
            f"ValueError computing CM: y_true: {y_true_filtered!r}, y_pred: {y_pred_bool!r}"
        )
        metrics["error"] = "Could not compute confusion matrix due to ValueError."

    metrics["Total Compared"] = len(y_true_filtered)
    metrics["True Positives (TP)"] = int(tp)
    metrics["True Negatives (TN)"] = int(tn)
    metrics["False Positives (FP)"] = int(fp)
    metrics["False Negatives (FN)"] = int(fn)

    if "error" in metrics:
        logger.warning(
            f"Returning early from metrics calculation due to CM error: {metrics['error']}"
        )
        metrics.setdefault("Sensitivity (Recall)", 0.0)
        metrics.setdefault("Specificity", 0.0)
        metrics.setdefault("Precision (PPV)", 0.0)
        metrics.setdefault("Negative Predictive Value (NPV)", 0.0)
        metrics.setdefault("Accuracy", 0.0)
        metrics.setdefault("F1 Score", 0.0)
        metrics.setdefault("Matthews Correlation Coefficient (MCC)", 0.0)
        metrics.setdefault("Cohen's Kappa", 0.0)
        metrics.setdefault("PABAK", -1.0)
        metrics.setdefault("Positive Likelihood Ratio (LR+)", float("inf"))
        metrics.setdefault("Negative Likelihood Ratio (LR-)", float("inf"))
        return metrics

    try:
        metrics["Sensitivity (Recall)"] = recall_score(
            y_true_filtered, y_pred_bool, zero_division="warn"
        )
        metrics["Specificity"] = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
        metrics["Precision (PPV)"] = precision_score(
            y_true_filtered, y_pred_bool, zero_division="warn"
        )
        metrics["Negative Predictive Value (NPV)"] = (
            float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0
        )
        metrics["Accuracy"] = accuracy_score(y_true_filtered, y_pred_bool)
        metrics["F1 Score"] = f1_score(
            y_true_filtered, y_pred_bool, zero_division="warn"
        )

        if (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) == 0:
            metrics["Matthews Correlation Coefficient (MCC)"] = 0.0
            logger.warning(
                "MCC set to 0.0 due to zero in a CM quadrant sum for denominator."
            )
        else:
            metrics["Matthews Correlation Coefficient (MCC)"] = matthews_corrcoef(
                y_true_filtered, y_pred_bool
            )
        metrics["Cohen's Kappa"] = cohen_kappa_score(y_true_filtered, y_pred_bool)
        metrics["PABAK"] = (2 * metrics["Accuracy"]) - 1

        sensitivity = metrics["Sensitivity (Recall)"]
        specificity = metrics["Specificity"]
        metrics["Positive Likelihood Ratio (LR+)"] = (
            sensitivity / (1 - specificity)
            if (1 - specificity) > 0 and specificity != 1.0
            else float("inf")
        )
        metrics["Negative Likelihood Ratio (LR-)"] = (
            (1 - sensitivity) / specificity
            if specificity > 0 and specificity != 0.0
            else float("inf")
        )

    except Exception:
        logger.exception("Error calculating some metrics (post-CM)")
        if "error" not in metrics or "calculation_error" not in metrics:
            metrics["calculation_error"] = "Failed during secondary metric calculation."

    valid_confidences = []
    for original_idx in valid_indices:
        if (
            original_idx < len(y_pred_confidence)
            and y_pred_confidence[original_idx] is not None
        ):
            valid_confidences.append(y_pred_confidence[original_idx])

    if valid_confidences:
        metrics["Mean AI Confidence"] = mean(valid_confidences)
        if len(valid_confidences) > 1:
            metrics["Median AI Confidence"] = median(valid_confidences)
            metrics["StdDev AI Confidence"] = stdev(valid_confidences)

    return metrics


# WARN: This can be only invoked once (in main.py): `st.set_page_config(layout="wide", page_title="SRA Benchmark Tool")`

st.title("SR Assistant - Benchmark Tool")

if "benchmark_review" not in st.session_state:
    st.session_state.benchmark_review = None
if "benchmark_search_results" not in st.session_state:
    st.session_state.benchmark_search_results = []
if "benchmark_ai_decisions" not in st.session_state:
    st.session_state.benchmark_ai_decisions = []
if "benchmark_metrics" not in st.session_state:
    st.session_state.benchmark_metrics = None
if "benchmark_comparison_data" not in st.session_state:
    st.session_state.benchmark_comparison_data = []


with st.sidebar:
    st.header("Controls")
    if st.button("Load Benchmark Dataset"):
        with session_factory() as session:
            review_repo = SystematicReviewRepository()
            search_repo = SearchResultRepository()
            loaded_review: models.SystematicReview | None = review_repo.get_by_id(
                session=session, id=BENCHMARK_REVIEW_ID
            )
            st.session_state.benchmark_review = loaded_review

            if st.session_state.benchmark_review:
                loaded_search_results: collections.abc.Sequence[models.SearchResult] = (
                    search_repo.get_by_review_id(
                        session=session, review_id=BENCHMARK_REVIEW_ID
                    )
                )
                st.session_state.benchmark_search_results = list(loaded_search_results)

                st.success(
                    f"Loaded benchmark review '({st.session_state.benchmark_review.id}) and {len(st.session_state.benchmark_search_results)} search results."
                )
                st.session_state.benchmark_ai_decisions = []
                st.session_state.benchmark_metrics = None
                st.session_state.benchmark_comparison_data = []
            else:
                st.error(
                    f"Benchmark review with ID {BENCHMARK_REVIEW_ID} not found. Please seed data first."
                )

if st.session_state.benchmark_review and st.session_state.benchmark_search_results:
    st.header("Benchmark Protocol")
    review: models.SystematicReview = st.session_state.benchmark_review
    st.markdown(f"**Research Question:** {review.research_question}")
    st.markdown(f"**Background:** {review.background}")
    st.expander("Inclusion/Exclusion Criteria").write(
        {
            "PICO": review.criteria_framework_answers,
            "Inclusion String": review.inclusion_criteria,
            "Exclusion String": review.exclusion_criteria,
        }
    )
    st.info(
        f"{len(st.session_state.benchmark_search_results)} abstracts loaded for benchmark."
    )

    if st.button("Run AI Screening on Benchmark", type="primary"):
        # Initialize benchmark execution in session state
        st.session_state.benchmark_running = True
        st.session_state.benchmark_progress = 0.0
        st.session_state.benchmark_status = "Initializing benchmark run..."
        st.session_state.benchmark_batch_num = 0
        st.session_state.benchmark_total_batches = 0
        st.session_state.benchmark_run_id = None
        st.session_state.benchmark_search_results = []
        st.session_state.benchmark_stats = {
            "total_processed": 0,
            "conflicts_detected": 0,
            "resolver_invoked": 0,
            "screening_errors": 0,
            "accumulated_results": [],  # Store results for real-time metrics
        }
        st.session_state.benchmark_current_metrics = {}  # Clear previous metrics
        st.session_state.benchmark_phase = "creating_run"
        st.rerun()

# Handle benchmark execution phases
if st.session_state.get("benchmark_running", False):
    # Display current progress
    st.subheader("🔄 Benchmark Execution in Progress")

    # Progress bar
    progress_value = st.session_state.get("benchmark_progress", 0.0)
    st.progress(progress_value)

    # Status display
    status_text = st.session_state.get("benchmark_status", "Processing...")
    st.info(f"**Status:** {status_text}")

    # Stats display
    stats = st.session_state.get("benchmark_stats", {})
    if stats.get("total_processed", 0) > 0:
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Items Processed", stats["total_processed"])
        with col2:
            st.metric("Conflicts Detected", stats["conflicts_detected"])
        with col3:
            st.metric("Resolver Invoked", stats["resolver_invoked"])
        with col4:
            st.metric("Screening Errors", stats["screening_errors"])

        # Display real-time performance metrics
        current_metrics = st.session_state.get("benchmark_current_metrics", {})
        if current_metrics and "error" not in current_metrics:
            st.subheader("📊 Live Performance Metrics")

            # Key metrics in columns
            metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)

            with metric_col1:
                accuracy = current_metrics.get("Accuracy", 0.0)
                st.metric("Accuracy", f"{accuracy:.3f}" if accuracy else "0.000")

                precision = current_metrics.get("Precision (PPV)", 0.0)
                st.metric("Precision", f"{precision:.3f}" if precision else "0.000")

            with metric_col2:
                sensitivity = current_metrics.get("Sensitivity (Recall)", 0.0)
                st.metric(
                    "Recall/Sensitivity",
                    f"{sensitivity:.3f}" if sensitivity else "0.000",
                )

                specificity = current_metrics.get("Specificity", 0.0)
                st.metric(
                    "Specificity", f"{specificity:.3f}" if specificity else "0.000"
                )

            with metric_col3:
                f1_score = current_metrics.get("F1 Score", 0.0)
                st.metric("F1 Score", f"{f1_score:.3f}" if f1_score else "0.000")

                mcc = current_metrics.get("Matthews Correlation Coefficient (MCC)", 0.0)
                st.metric("MCC", f"{mcc:.3f}" if mcc else "0.000")

            with metric_col4:
                tp = current_metrics.get("True Positives (TP)", 0)
                fp = current_metrics.get("False Positives (FP)", 0)
                fn = current_metrics.get("False Negatives (FN)", 0)
                tn = current_metrics.get("True Negatives (TN)", 0)

                st.metric("True Positives", tp)
                st.metric("False Positives", fp)

            # Confusion matrix summary
            total_compared = current_metrics.get("Total Compared", 0)
            if total_compared > 0:
                st.info(
                    f"📈 **Current Analysis**: {total_compared} items with ground truth | TP: {tp} | TN: {tn} | FP: {fp} | FN: {fn}"
                )

        elif current_metrics and "error" in current_metrics:
            st.warning(f"⚠️ Metrics calculation: {current_metrics['error']}")

        elif stats.get("total_processed", 0) >= 5:  # Show message after a few items
            st.info(
                "🔄 **Real-time metrics will appear once sufficient data with ground truth is available**"
            )

    # Phase-based execution
    phase = st.session_state.get("benchmark_phase", "creating_run")

    if phase == "creating_run":
        # Phase 1: Create benchmark run
        st.session_state.benchmark_status = "Creating benchmark run in database..."

        with session_factory() as session:
            try:
                benchmark_run_repo = BenchmarkRunRepository()

                config_details = {
                    "conservative_model": "gpt-4o",
                    "comprehensive_model": "gpt-4o",
                    "resolver_model": "gpt-4o",
                    "batch_size": 10,
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                }

                benchmark_run = models.BenchmarkRun(
                    review_id=BENCHMARK_REVIEW_ID, config_details=config_details
                )

                benchmark_run = benchmark_run_repo.add(session, benchmark_run)
                session.flush()
                session.commit()

                # Store run ID
                st.session_state.benchmark_run_id = benchmark_run.id

                # Get search results
                search_repo = SearchResultRepository()
                search_results = search_repo.get_by_review_id(
                    session, BENCHMARK_REVIEW_ID
                )

                if not search_results:
                    st.error("No search results found for benchmark.")
                    st.session_state.benchmark_running = False
                    st.rerun()

                # Store search results and calculate batches
                st.session_state.benchmark_search_results = list(search_results)
                total_items = len(search_results)
                batch_size = 10
                st.session_state.benchmark_total_batches = (
                    total_items + batch_size - 1
                ) // batch_size

                # Move to next phase
                st.session_state.benchmark_phase = "processing_batches"
                st.session_state.benchmark_status = f"Created run {benchmark_run.id}. Ready to process {total_items} items in {st.session_state.benchmark_total_batches} batches."

                logger.info(
                    f"Created BenchmarkRun {benchmark_run.id} with {total_items} items"
                )

            except Exception as e:
                logger.exception("Failed to create benchmark run")
                st.error(f"Failed to create benchmark run: {e}")
                st.session_state.benchmark_running = False
                st.rerun()

        # Force UI update to next phase
        st.rerun()

    elif phase == "processing_batches":
        # Phase 2: Process items incrementally (1-2 items at a time for frequent updates)

        # Check if we have items left to process
        total_items = len(st.session_state.benchmark_search_results)
        processed_items = st.session_state.benchmark_stats["total_processed"]

        if processed_items < total_items:
            # Process next 1-2 items for frequent updates
            items_per_cycle = (
                2  # Process 2 items per UI cycle for better responsiveness
            )
            start_idx = processed_items
            end_idx = min(start_idx + items_per_cycle, total_items)

            current_batch_items = st.session_state.benchmark_search_results[
                start_idx:end_idx
            ]

            # Update status
            batch_num = (processed_items // 10) + 1
            item_in_batch = (processed_items % 10) + 1
            st.session_state.benchmark_status = f"Processing items {start_idx + 1}-{end_idx} (batch {batch_num}, items {item_in_batch}-{item_in_batch + len(current_batch_items) - 1})..."

            # Process these items
            with st.spinner(f"Processing items {start_idx + 1}-{end_idx}..."):
                with session_factory() as session:
                    try:
                        benchmark_result_item_repo = BenchmarkResultItemRepository()

                        # Call screen_abstracts_batch for this small chunk
                        if not st.session_state.benchmark_review:
                            st.error("Benchmark review not found in session state.")
                            st.session_state.benchmark_running = False
                            st.rerun()

                        # Use batch index based on the first item
                        batch_idx = start_idx // 10

                        batch_output = screen_abstracts_batch(
                            batch=list(current_batch_items),
                            batch_idx=batch_idx,
                            review=st.session_state.benchmark_review,
                        )

                        if not batch_output or not batch_output.results:
                            logger.error(
                                f"No screening results for items {start_idx + 1}-{end_idx}"
                            )
                            st.session_state.benchmark_stats["screening_errors"] += len(
                                current_batch_items
                            )
                            # Skip to next items
                            st.session_state.benchmark_stats["total_processed"] += len(
                                current_batch_items
                            )
                            st.rerun()

                        # Process each result in this small chunk
                        for result_tuple in batch_output.results:
                            search_result = result_tuple.search_result
                            conservative_result = result_tuple.conservative_result
                            comprehensive_result = result_tuple.comprehensive_result

                            # Check if we have valid results
                            if isinstance(
                                conservative_result, ScreeningError
                            ) or isinstance(comprehensive_result, ScreeningError):
                                logger.warning(
                                    f"Screening error for {search_result.source_id}"
                                )
                                st.session_state.benchmark_stats[
                                    "screening_errors"
                                ] += 1
                                st.session_state.benchmark_stats["total_processed"] += 1
                                continue

                            # Check for conflicts that need resolver
                            needs_resolver = _needs_resolver(
                                conservative_result, comprehensive_result
                            )
                            resolver_result = None

                            if needs_resolver:
                                st.session_state.benchmark_stats[
                                    "conflicts_detected"
                                ] += 1

                                # Invoke resolver immediately for this item
                                st.session_state.benchmark_status = f"Resolving conflict for item {st.session_state.benchmark_stats['total_processed'] + 1}..."
                                logger.info(
                                    f"Invoking resolver for search result {search_result.id}"
                                )

                                try:
                                    if not st.session_state.benchmark_review:
                                        st.error(
                                            "Benchmark review not found in session state."
                                        )
                                        logger.error(
                                            "Benchmark review not found in session state"
                                        )
                                        break

                                    resolver_output = invoke_resolver_chain(
                                        search_result=search_result,
                                        conservative_result=conservative_result,
                                        comprehensive_result=comprehensive_result,
                                        review=st.session_state.benchmark_review,
                                    )

                                    logger.info(
                                        f"Resolver output type: {type(resolver_output)}, value: {resolver_output}"
                                    )

                                    if resolver_output and isinstance(
                                        resolver_output, schemas.ResolverOutputSchema
                                    ):  # pyright: ignore[reportUnnecessaryIsInstance]
                                        # Create a mock ScreeningResult from ResolverOutputSchema for consistency
                                        resolver_result = schemas.ScreeningResult(
                                            id=uuid.uuid4(),
                                            review_id=st.session_state.benchmark_review.id,
                                            search_result_id=search_result.id,
                                            trace_id=uuid.uuid4(),
                                            model_name="resolver",
                                            screening_strategy=schemas.ScreeningStrategyType.COMPREHENSIVE,
                                            start_time=datetime.now(timezone.utc),
                                            end_time=datetime.now(timezone.utc),
                                            decision=resolver_output.resolver_decision,
                                            confidence_score=resolver_output.resolver_confidence_score,
                                            rationale=resolver_output.resolver_reasoning,
                                        )

                                        st.session_state.benchmark_stats[
                                            "resolver_invoked"
                                        ] += 1
                                        logger.info(
                                            f"Resolver successfully invoked. Total invocations: {st.session_state.benchmark_stats['resolver_invoked']}"
                                        )
                                    else:
                                        logger.warning(
                                            f"Resolver returned invalid result: {resolver_output}"
                                        )
                                except Exception as e:
                                    logger.error(
                                        f"Resolver failed for {search_result.source_id}: {e}"
                                    )

                            # Determine final decision
                            final_decision = _determine_final_decision(
                                conservative_result,
                                comprehensive_result,
                                resolver_result,
                            )

                            # Get human decision from source_metadata
                            human_decision_raw = search_result.source_metadata.get(
                                "benchmark_human_decision"
                            )
                            human_decision: bool | None = None
                            if isinstance(human_decision_raw, bool):
                                human_decision = human_decision_raw

                            # Calculate classification
                            classification = _calculate_classification(
                                final_decision, human_decision
                            )

                            # Create BenchmarkResultItem
                            benchmark_result_item = models.BenchmarkResultItem(
                                benchmark_run_id=st.session_state.benchmark_run_id,
                                search_result_id=search_result.id,
                                human_decision=human_decision,
                                conservative_decision=conservative_result.decision,
                                conservative_confidence=conservative_result.confidence_score,
                                conservative_rationale=conservative_result.rationale,
                                comprehensive_decision=comprehensive_result.decision,
                                comprehensive_confidence=comprehensive_result.confidence_score,
                                comprehensive_rationale=comprehensive_result.rationale,
                                resolver_decision=resolver_result.decision
                                if resolver_result
                                else None,
                                resolver_confidence=resolver_result.confidence_score
                                if resolver_result
                                else None,
                                resolver_reasoning=resolver_result.rationale
                                if resolver_result
                                else None,
                                final_decision=final_decision,
                                classification=classification,
                            )

                            benchmark_result_item_repo.add(
                                session, benchmark_result_item
                            )
                            st.session_state.benchmark_stats["total_processed"] += 1

                            # Store result for real-time metrics calculation
                            result_data = {
                                "human_decision": human_decision,
                                "final_decision": final_decision,
                                "conservative_confidence": conservative_result.confidence_score,
                                "comprehensive_confidence": comprehensive_result.confidence_score,
                                "resolver_confidence": resolver_result.confidence_score
                                if resolver_result
                                else None,
                            }
                            st.session_state.benchmark_stats[
                                "accumulated_results"
                            ].append(result_data)

                        # Commit this chunk
                        session.commit()

                        # Calculate and update real-time metrics
                        accumulated_results = st.session_state.benchmark_stats[
                            "accumulated_results"
                        ]
                        if accumulated_results:
                            # Prepare data for metrics calculation
                            y_true = []
                            y_pred_decision = []
                            y_pred_confidence = []

                            for result in accumulated_results:
                                y_true.append(result["human_decision"])
                                y_pred_decision.append(result["final_decision"])

                                # Use the highest confidence score available
                                confidence = result["conservative_confidence"]
                                if result["comprehensive_confidence"] is not None:
                                    confidence = max(
                                        confidence, result["comprehensive_confidence"]
                                    )
                                if result["resolver_confidence"] is not None:
                                    confidence = result["resolver_confidence"]
                                y_pred_confidence.append(confidence)

                            # Calculate current metrics
                            current_metrics = calculate_metrics(
                                y_true, y_pred_decision, y_pred_confidence
                            )
                            st.session_state.benchmark_current_metrics = current_metrics

                        # Update progress
                        st.session_state.benchmark_progress = end_idx / total_items

                        logger.info(
                            f"Completed processing items {start_idx + 1}-{end_idx}"
                        )

                    except Exception as e:
                        session.rollback()
                        logger.exception(f"Items {start_idx + 1}-{end_idx} failed")
                        st.error(f"Items {start_idx + 1}-{end_idx} failed: {e}")
                        st.session_state.benchmark_running = False
                        st.rerun()

            # Force UI update for next items (more frequent updates)
            st.rerun()

        else:
            # All items completed - move to completion phase
            st.session_state.benchmark_phase = "completed"
            st.rerun()

    elif phase == "completed":
        # Phase 3: Benchmark completed
        st.session_state.benchmark_running = False
        st.session_state.benchmark_status = "✅ Benchmark run completed successfully!"
        st.session_state.benchmark_progress = 1.0

        st.success(f"""
        🎉 **Benchmark Run Complete!**
        
        **Run ID:** `{st.session_state.benchmark_run_id}`
        - **Total Items Processed:** {st.session_state.benchmark_stats["total_processed"]}
        - **Conflicts Detected:** {st.session_state.benchmark_stats["conflicts_detected"]}
        - **Resolver Invocations:** {st.session_state.benchmark_stats["resolver_invoked"]}
        - **Screening Errors:** {st.session_state.benchmark_stats["screening_errors"]}
        
        All results have been saved to the database.
        """)

        # Store run ID for potential metrics calculation
        st.session_state.current_benchmark_run_id = st.session_state.benchmark_run_id

        # Clean up session state
        for key in [
            "benchmark_running",
            "benchmark_progress",
            "benchmark_status",
            "benchmark_batch_num",
            "benchmark_total_batches",
            "benchmark_run_id",
            "benchmark_stats",
            "benchmark_search_results",
            "benchmark_phase",
            "benchmark_current_metrics",
        ]:
            if key in st.session_state:
                del st.session_state[key]


if st.session_state.benchmark_metrics:
    st.header("Benchmark Results")
    metrics_display = st.session_state.benchmark_metrics
    if "error" in metrics_display:
        st.error(f"Could not calculate metrics: {metrics_display['error']}")
    else:
        st.subheader("Performance Metrics")
        cols_metrics = st.columns(3)
        metric_keys_display = list(metrics_display.keys())
        for i, key in enumerate(metric_keys_display):
            with cols_metrics[i % 3]:
                if isinstance(metrics_display[key], float):
                    st.metric(label=key, value=f"{metrics_display[key]:.4f}")
                else:
                    st.metric(label=key, value=metrics_display[key])

        st.subheader("Confusion Matrix")
        cm_tn = metrics_display.get("True Negatives (TN)", 0)
        cm_fp = metrics_display.get("False Positives (FP)", 0)
        cm_fn = metrics_display.get("False Negatives (FN)", 0)
        cm_tp = metrics_display.get("True Positives (TP)", 0)
        cm_df = pd.DataFrame(
            {
                "Predicted: EXCLUDE": [cm_tn, cm_fn],
                "Predicted: INCLUDE": [cm_fp, cm_tp],
            },
            index=["Actual: EXCLUDE", "Actual: INCLUDE"],
        )
        st.table(cm_df)

if st.session_state.benchmark_comparison_data:
    st.subheader("Detailed Comparison")
    df_compare = pd.DataFrame(st.session_state.benchmark_comparison_data)

    def classify_row(row: pd.Series[t.Any]) -> str:
        ai = row["ai_decision"]
        human = row["human_decision"]
        if ai == ScreeningDecisionType.INCLUDE.value and human == "INCLUDE":
            return "TP"
        if ai == ScreeningDecisionType.EXCLUDE.value and human == "EXCLUDE":
            return "TN"
        if ai == ScreeningDecisionType.INCLUDE.value and human == "EXCLUDE":
            return "FP"
        if ai == ScreeningDecisionType.EXCLUDE.value and human == "INCLUDE":
            return "FN"
        return "UNCERTAIN/ERROR"

    df_compare["classification"] = df_compare.apply(classify_row, axis=1)
    st.dataframe(df_compare, use_container_width=True, hide_index=True)

    csv_buffer = io.StringIO()
    df_compare.to_csv(csv_buffer, index=False)
    csv_export_bytes = csv_buffer.getvalue().encode("utf-8")
    csv_buffer.close()

    st.download_button(
        label="Download Comparison Data as CSV",
        data=csv_export_bytes,
        file_name=f"benchmark_comparison_results_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.csv",
        mime="text/csv",
        key="download_comparison_csv",
    )

    if (
        st.session_state.benchmark_metrics
        and "error" not in st.session_state.benchmark_metrics
    ):
        metrics_summary_dict = st.session_state.benchmark_metrics
        metrics_export_str = "Metric,Value\n"
        for key, value in metrics_summary_dict.items():
            metrics_export_str += f"{key},{value}\n"

        metrics_buffer = io.StringIO()
        metrics_buffer.write(metrics_export_str)
        metrics_export_bytes = metrics_buffer.getvalue().encode("utf-8")
        metrics_buffer.close()

        st.download_button(
            label="Download Summary Metrics as CSV",
            data=metrics_export_bytes,
            file_name=f"benchmark_summary_metrics_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv",
            key="download_metrics_csv",
        )
