# pyright: reportUnnecessaryComparison=false
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
from sr_assistant.benchmark.logic.metrics_calculator import (
    calculate_and_update_benchmark_metrics,
)
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
    return min_confidence < 0.7


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
    final_decision: ScreeningDecisionType, human_decision: bool | None
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

# Initialize session state keys if they don't exist
EXPECTED_SESSION_STATE_KEYS = [
    "benchmark_review",
    "benchmark_search_results",
    "benchmark_ai_decisions",
    "benchmark_metrics",
    "benchmark_comparison_data",
    # "selected_benchmark_run", # Temporarily remove as we're disabling the section
]
for key in EXPECTED_SESSION_STATE_KEYS:
    if key not in st.session_state:
        st.session_state[key] = None

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
                # st.session_state.selected_benchmark_run = None # Clear any previously loaded run
            else:
                st.error(
                    f"Benchmark review with ID {BENCHMARK_REVIEW_ID} not found. Please seed data first."
                )
        st.rerun()

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
        # st.session_state.benchmark_search_results = [] # Keep existing loaded results
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
    st.subheader("🔄 Benchmark Execution in Progress")
    progress_value = st.session_state.get("benchmark_progress", 0.0)
    st.progress(progress_value)
    status_text = st.session_state.get("benchmark_status", "Processing...")
    st.info(f"**Status:** {status_text}")

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

        current_metrics = st.session_state.get("benchmark_current_metrics", {})
        if current_metrics and "error" not in current_metrics:
            st.subheader("📊 Live Performance Metrics")

            # Store previous metrics for delta calculation
            prev_metrics = st.session_state.get("benchmark_prev_metrics", {})

            metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)

            with metric_col1:
                accuracy_val = current_metrics.get("Accuracy", 0.0)
                accuracy_delta = (
                    accuracy_val - prev_metrics.get("Accuracy", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Accuracy",
                    f"{accuracy_val:.3f}" if accuracy_val else "0.000",
                    delta=f"{accuracy_delta:+.3f}"
                    if accuracy_delta is not None and abs(accuracy_delta) >= 0.001
                    else None,
                )

                sensitivity_val = current_metrics.get("Sensitivity (Recall)", 0.0)
                sensitivity_delta = (
                    sensitivity_val - prev_metrics.get("Sensitivity (Recall)", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Sensitivity",
                    f"{sensitivity_val:.3f}" if sensitivity_val else "0.000",
                    delta=f"{sensitivity_delta:+.3f}"
                    if sensitivity_delta is not None and abs(sensitivity_delta) >= 0.001
                    else None,
                )

            with metric_col2:
                precision_val = current_metrics.get("Precision (PPV)", 0.0)
                precision_delta = (
                    precision_val - prev_metrics.get("Precision (PPV)", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Precision",
                    f"{precision_val:.3f}" if precision_val else "0.000",
                    delta=f"{precision_delta:+.3f}"
                    if precision_delta is not None and abs(precision_delta) >= 0.001
                    else None,
                )

                specificity_val = current_metrics.get("Specificity", 0.0)
                specificity_delta = (
                    specificity_val - prev_metrics.get("Specificity", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Specificity",
                    f"{specificity_val:.3f}" if specificity_val else "0.000",
                    delta=f"{specificity_delta:+.3f}"
                    if specificity_delta is not None and abs(specificity_delta) >= 0.001
                    else None,
                )

            with metric_col3:
                f1_score_val = current_metrics.get("F1 Score", 0.0)
                f1_delta = (
                    f1_score_val - prev_metrics.get("F1 Score", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "F1 Score",
                    f"{f1_score_val:.3f}" if f1_score_val else "0.000",
                    delta=f"{f1_delta:+.3f}"
                    if f1_delta is not None and abs(f1_delta) >= 0.001
                    else None,
                )

                npv_val = current_metrics.get("Negative Predictive Value (NPV)", 0.0)
                npv_delta = (
                    npv_val - prev_metrics.get("Negative Predictive Value (NPV)", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "NPV",
                    f"{npv_val:.3f}" if npv_val else "0.000",
                    delta=f"{npv_delta:+.3f}"
                    if npv_delta is not None and abs(npv_delta) >= 0.001
                    else None,
                )

            with metric_col4:
                mcc_val = current_metrics.get(
                    "Matthews Correlation Coefficient (MCC)", 0.0
                )
                mcc_delta = (
                    mcc_val
                    - prev_metrics.get("Matthews Correlation Coefficient (MCC)", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "MCC",
                    f"{mcc_val:.3f}" if mcc_val else "0.000",
                    delta=f"{mcc_delta:+.3f}"
                    if mcc_delta is not None and abs(mcc_delta) >= 0.001
                    else None,
                )

                cohen_kappa_val = current_metrics.get("Cohen's Kappa", 0.0)
                cohen_delta = (
                    cohen_kappa_val - prev_metrics.get("Cohen's Kappa", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Cohen's κ",
                    f"{cohen_kappa_val:.3f}" if cohen_kappa_val else "0.000",
                    delta=f"{cohen_delta:+.3f}"
                    if cohen_delta is not None and abs(cohen_delta) >= 0.001
                    else None,
                )

            # Additional metrics row
            st.subheader("🎯 Advanced Metrics")
            adv_col1, adv_col2, adv_col3, adv_col4 = st.columns(4)

            with adv_col1:
                pabak_val = current_metrics.get("PABAK", 0.0)
                pabak_delta = (
                    pabak_val - prev_metrics.get("PABAK", 0.0) if prev_metrics else None
                )
                st.metric(
                    "PABAK",
                    f"{pabak_val:.3f}" if pabak_val else "0.000",
                    delta=f"{pabak_delta:+.3f}"
                    if pabak_delta is not None and abs(pabak_delta) >= 0.001
                    else None,
                )

                mean_conf_val = current_metrics.get("Mean AI Confidence", 0.0)
                mean_conf_delta = (
                    mean_conf_val - prev_metrics.get("Mean AI Confidence", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Mean Confidence",
                    f"{mean_conf_val:.3f}" if mean_conf_val else "0.000",
                    delta=f"{mean_conf_delta:+.3f}"
                    if mean_conf_delta is not None and abs(mean_conf_delta) >= 0.001
                    else None,
                )

            with adv_col2:
                lr_plus_val = current_metrics.get(
                    "Positive Likelihood Ratio (LR+)", 0.0
                )
                if lr_plus_val == float("inf"):
                    lr_plus_str = "∞"
                    lr_plus_delta_str = None
                else:
                    lr_plus_str = f"{lr_plus_val:.3f}" if lr_plus_val else "0.000"
                    lr_plus_delta = (
                        lr_plus_val
                        - prev_metrics.get("Positive Likelihood Ratio (LR+)", 0.0)
                        if prev_metrics
                        else None
                    )
                    lr_plus_delta_str = (
                        f"{lr_plus_delta:+.3f}"
                        if lr_plus_delta is not None and abs(lr_plus_delta) >= 0.001
                        else None
                    )
                st.metric("LR+", lr_plus_str, delta=lr_plus_delta_str)

                median_conf_val = current_metrics.get("Median AI Confidence", 0.0)
                median_conf_delta = (
                    median_conf_val - prev_metrics.get("Median AI Confidence", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Median Confidence",
                    f"{median_conf_val:.3f}" if median_conf_val else "0.000",
                    delta=f"{median_conf_delta:+.3f}"
                    if median_conf_delta is not None and abs(median_conf_delta) >= 0.001
                    else None,
                )

            with adv_col3:
                lr_minus_val = current_metrics.get(
                    "Negative Likelihood Ratio (LR-)", 0.0
                )
                if lr_minus_val == float("inf"):
                    lr_minus_str = "∞"
                    lr_minus_delta_str = None
                else:
                    lr_minus_str = f"{lr_minus_val:.3f}" if lr_minus_val else "0.000"
                    lr_minus_delta = (
                        lr_minus_val
                        - prev_metrics.get("Negative Likelihood Ratio (LR-)", 0.0)
                        if prev_metrics
                        else None
                    )
                    lr_minus_delta_str = (
                        f"{lr_minus_delta:+.3f}"
                        if lr_minus_delta is not None and abs(lr_minus_delta) >= 0.001
                        else None
                    )
                st.metric("LR-", lr_minus_str, delta=lr_minus_delta_str)

                stdev_conf_val = current_metrics.get("StdDev AI Confidence", 0.0)
                stdev_conf_delta = (
                    stdev_conf_val - prev_metrics.get("StdDev AI Confidence", 0.0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "StdDev Confidence",
                    f"{stdev_conf_val:.3f}" if stdev_conf_val else "0.000",
                    delta=f"{stdev_conf_delta:+.3f}"
                    if stdev_conf_delta is not None and abs(stdev_conf_delta) >= 0.001
                    else None,
                )

            with adv_col4:
                total_compared = current_metrics.get("Total Compared", 0)
                total_delta = (
                    total_compared - prev_metrics.get("Total Compared", 0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "Items Analyzed",
                    str(total_compared),
                    delta=f"+{total_delta}"
                    if total_delta and total_delta > 0
                    else None,
                )

            # Confusion Matrix Summary with deltas
            st.subheader("📈 Confusion Matrix")
            cm_col1, cm_col2, cm_col3, cm_col4 = st.columns(4)

            with cm_col1:
                tp_val = current_metrics.get("True Positives (TP)", 0)
                tp_delta = (
                    tp_val - prev_metrics.get("True Positives (TP)", 0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "True Positives",
                    str(tp_val),
                    delta=f"+{tp_delta}" if tp_delta and tp_delta > 0 else None,
                )

            with cm_col2:
                tn_val = current_metrics.get("True Negatives (TN)", 0)
                tn_delta = (
                    tn_val - prev_metrics.get("True Negatives (TN)", 0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "True Negatives",
                    str(tn_val),
                    delta=f"+{tn_delta}" if tn_delta and tn_delta > 0 else None,
                )

            with cm_col3:
                fp_val = current_metrics.get("False Positives (FP)", 0)
                fp_delta = (
                    fp_val - prev_metrics.get("False Positives (FP)", 0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "False Positives",
                    str(fp_val),
                    delta=f"+{fp_delta}" if fp_delta and fp_delta > 0 else None,
                )

            with cm_col4:
                fn_val = current_metrics.get("False Negatives (FN)", 0)
                fn_delta = (
                    fn_val - prev_metrics.get("False Negatives (FN)", 0)
                    if prev_metrics
                    else None
                )
                st.metric(
                    "False Negatives",
                    str(fn_val),
                    delta=f"+{fn_delta}" if fn_delta and fn_delta > 0 else None,
                )

            # Store current metrics as previous for next update
            st.session_state.benchmark_prev_metrics = current_metrics.copy()

            # Confusion matrix summary
            total_compared = current_metrics.get("Total Compared", 0)
            if total_compared > 0:
                st.info(
                    f"📈 **Current Analysis**: {total_compared} items with ground truth | TP: {tp_val} | TN: {tn_val} | FP: {fp_val} | FN: {fn_val}"
                )

        elif current_metrics and "error" in current_metrics:
            st.warning(f"⚠️ Metrics calculation: {current_metrics['error']}")

        elif stats.get("total_processed", 0) >= 5:
            st.info(
                "🔄 **Real-time metrics will appear once sufficient data with ground truth is available**"
            )

        # LIVE TABLE DISPLAY - Shows during execution
        st.subheader("🔎 Live Processing Results")

        accumulated_results = stats.get("accumulated_results", [])
        search_results = st.session_state.get("benchmark_search_results", [])

        if accumulated_results and search_results:
            # Create DataFrame for live results
            live_data = []
            for i, result_data in enumerate(accumulated_results):
                if i < len(search_results):
                    search_result = search_results[i]

                    row = {
                        "Title": search_result.title[:100]
                        + ("..." if len(search_result.title) > 100 else ""),
                        "Year": search_result.year or "N/A",
                        "Human Decision": "Include"
                        if result_data["human_decision"] is True
                        else "Exclude"
                        if result_data["human_decision"] is False
                        else "Unknown",
                        "Final Decision": result_data["final_decision"].value,
                        "Classification": result_data["classification"],
                        "Conservative Decision": result_data[
                            "conservative_decision"
                        ].value
                        if result_data["conservative_decision"]
                        else "N/A",
                        "Conservative Confidence": f"{result_data['conservative_confidence']:.3f}"
                        if result_data["conservative_confidence"] is not None
                        else "N/A",
                        "Conservative Rationale": result_data["conservative_rationale"][
                            :150
                        ]
                        + (
                            "..."
                            if result_data["conservative_rationale"]
                            and len(result_data["conservative_rationale"]) > 150
                            else ""
                        )
                        if result_data["conservative_rationale"]
                        else "N/A",
                        "Conservative Quotes": result_data["conservative_quotes"][:150]
                        + (
                            "..."
                            if result_data["conservative_quotes"]
                            and len(result_data["conservative_quotes"]) > 150
                            else ""
                        )
                        if result_data["conservative_quotes"]
                        and result_data["conservative_quotes"] != "N/A"
                        else "N/A",
                        "Comprehensive Decision": result_data[
                            "comprehensive_decision"
                        ].value
                        if result_data["comprehensive_decision"]
                        else "N/A",
                        "Comprehensive Confidence": f"{result_data['comprehensive_confidence']:.3f}"
                        if result_data["comprehensive_confidence"] is not None
                        else "N/A",
                        "Comprehensive Rationale": result_data[
                            "comprehensive_rationale"
                        ][:150]
                        + (
                            "..."
                            if result_data["comprehensive_rationale"]
                            and len(result_data["comprehensive_rationale"]) > 150
                            else ""
                        )
                        if result_data["comprehensive_rationale"]
                        else "N/A",
                        "Comprehensive Quotes": result_data["comprehensive_quotes"][
                            :150
                        ]
                        + (
                            "..."
                            if result_data["comprehensive_quotes"]
                            and len(result_data["comprehensive_quotes"]) > 150
                            else ""
                        )
                        if result_data["comprehensive_quotes"]
                        and result_data["comprehensive_quotes"] != "N/A"
                        else "N/A",
                        "Resolver Decision": result_data["resolver_decision"].value
                        if result_data["resolver_decision"]
                        else "N/A",
                        "Resolver Confidence": f"{result_data['resolver_confidence']:.3f}"
                        if result_data["resolver_confidence"] is not None
                        else "N/A",
                        "Resolver Reasoning": result_data["resolver_rationale"][:150]
                        + (
                            "..."
                            if result_data["resolver_rationale"]
                            and len(result_data["resolver_rationale"]) > 150
                            else ""
                        )
                        if result_data["resolver_rationale"]
                        else "N/A",
                        "Authors": (
                            ", ".join(search_result.authors[:2])
                            + (
                                "..."
                                if search_result.authors
                                and len(search_result.authors) > 2
                                else ""
                            )
                        )
                        if search_result.authors
                        else "N/A",
                        "Source ID": search_result.source_id,
                        "DOI": search_result.doi or "N/A",
                        "_conservative_rationale_full": result_data[
                            "conservative_rationale"
                        ]
                        or "",
                        "_comprehensive_rationale_full": result_data[
                            "comprehensive_rationale"
                        ]
                        or "",
                        "_resolver_rationale_full": result_data["resolver_rationale"]
                        or "",
                        "_conservative_quotes_full": result_data["conservative_quotes"]
                        or "",
                        "_comprehensive_quotes_full": result_data[
                            "comprehensive_quotes"
                        ]
                        or "",
                    }
                    live_data.append(row)

            if live_data:
                live_df = pd.DataFrame(live_data)

                # Filter controls for live results
                filter_col1, filter_col2, filter_col3 = st.columns(3)
                with filter_col1:
                    classification_options = ["All"] + sorted(
                        live_df["Classification"].unique().tolist()
                    )
                    selected_classification_live = st.selectbox(
                        "Filter by Classification:",
                        options=classification_options,
                        key="classification_filter_live_run",
                    )
                with filter_col2:
                    show_mismatches_only_live = st.checkbox(
                        "Show only mismatches (FP/FN)", key="mismatches_filter_live_run"
                    )
                with filter_col3:
                    # Export buttons
                    export_col1, export_col2 = st.columns(2)
                    with export_col1:
                        st.download_button(
                            label="📄 TSV",
                            data=live_df.to_csv(sep="\t", index=False),
                            file_name=f"live_benchmark_results_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.tsv",
                            mime="text/tab-separated-values",
                            key="export_live_results_tsv",
                        )
                    with export_col2:
                        # Convert DataFrame to Excel
                        excel_buffer = io.BytesIO()
                        with pd.ExcelWriter(excel_buffer, engine="openpyxl") as writer:
                            live_df.to_excel(
                                writer, index=False, sheet_name="Live Results"
                            )
                        excel_data = excel_buffer.getvalue()

                        st.download_button(
                            label="📊 XLSX",
                            data=excel_data,
                            file_name=f"live_benchmark_results_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            key="export_live_results_xlsx",
                        )

                # Apply filters
                df_to_display_live = live_df
                if selected_classification_live != "All":
                    df_to_display_live = df_to_display_live[
                        df_to_display_live["Classification"]
                        == selected_classification_live
                    ]
                if show_mismatches_only_live:
                    df_to_display_live = df_to_display_live[
                        df_to_display_live["Classification"].isin(["FP", "FN"])
                    ]

                st.info(
                    f"Showing {len(df_to_display_live)} of {len(live_df)} results from current run"
                )

                # Display table
                display_columns_live = [
                    col for col in df_to_display_live.columns if not col.startswith("_")
                ]

                # Interactive table with selection
                selected_indices = st.dataframe(
                    df_to_display_live[display_columns_live],
                    use_container_width=True,
                    hide_index=True,
                    height=400,
                    on_select="rerun",
                    selection_mode="single-row",
                )

                # Paper detail view for live results
                if (
                    selected_indices
                    and "selection" in selected_indices
                    and "rows" in selected_indices["selection"]
                ):
                    if selected_indices["selection"]["rows"]:
                        selected_row_idx = selected_indices["selection"]["rows"][0]
                        selected_paper_live = df_to_display_live.iloc[selected_row_idx]

                        st.subheader("🔍 Paper Detail View")

                        # Display details for selected paper
                        detail_col1, detail_col2 = st.columns(2)
                        with detail_col1:
                            st.markdown("**Conservative Agent:**")
                            st.markdown(
                                f"Decision: {selected_paper_live['Conservative Decision']}"
                            )
                            st.markdown(
                                f"Confidence: {selected_paper_live['Conservative Confidence']}"
                            )
                            if selected_paper_live["_conservative_rationale_full"]:
                                with st.expander("Rationale"):
                                    st.markdown(
                                        selected_paper_live[
                                            "_conservative_rationale_full"
                                        ]
                                    )

                        with detail_col2:
                            st.markdown("**Comprehensive Agent:**")
                            st.markdown(
                                f"Decision: {selected_paper_live['Comprehensive Decision']}"
                            )
                            st.markdown(
                                f"Confidence: {selected_paper_live['Comprehensive Confidence']}"
                            )
                            if selected_paper_live["_comprehensive_rationale_full"]:
                                with st.expander("Rationale"):
                                    st.markdown(
                                        selected_paper_live[
                                            "_comprehensive_rationale_full"
                                        ]
                                    )

                        if selected_paper_live["Resolver Decision"] != "N/A":
                            st.markdown("**Resolver Agent:**")
                            st.markdown(
                                f"Decision: {selected_paper_live['Resolver Decision']}"
                            )
                            st.markdown(
                                f"Confidence: {selected_paper_live['Resolver Confidence']}"
                            )
                            if selected_paper_live["_resolver_rationale_full"]:
                                with st.expander("Reasoning"):
                                    st.markdown(
                                        selected_paper_live["_resolver_rationale_full"]
                                    )

                        st.markdown("---")
                        st.markdown(
                            f"**Human Ground Truth:** {selected_paper_live['Human Decision']}"
                        )
                        st.markdown(
                            f"**SRA Final Decision:** {selected_paper_live['Final Decision']}"
                        )
                        st.markdown(
                            f"**Classification:** {selected_paper_live['Classification']}"
                        )
        else:
            st.info("No processed results available to display yet.")

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
                    "resolver_model": "gemini-2.5-pro-preview-05-06",
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
                batch_size = 10  # Consistent batch size
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
        BATCH_SIZE_CONFIG = 10  # True batch size for the agent

        total_items = len(st.session_state.benchmark_search_results)
        processed_items_count = st.session_state.benchmark_stats["total_processed"]

        if processed_items_count < total_items:
            batch_start_idx = (
                processed_items_count // BATCH_SIZE_CONFIG
            ) * BATCH_SIZE_CONFIG
            batch_end_idx = min(batch_start_idx + BATCH_SIZE_CONFIG, total_items)
            agent_batch_idx = batch_start_idx // BATCH_SIZE_CONFIG

            # Check if we need to fetch a new batch for the AI agent
            if processed_items_count == batch_start_idx and not st.session_state.get(
                "current_batch_screening_results"
            ):
                current_conceptual_batch_items = (
                    st.session_state.benchmark_search_results[
                        batch_start_idx:batch_end_idx
                    ]
                )

                # Update status to indicate AI batch screening is about to happen
                st.session_state.benchmark_status = f"AI processing batch {agent_batch_idx + 1} (items {batch_start_idx + 1}-{batch_end_idx})..."
                # Consider a spinner here if screen_abstracts_batch is very long, but avoid st.rerun()
                # For now, we assume screen_abstracts_batch is reasonably fast or its internal calls are async if needed

                if not st.session_state.benchmark_review:
                    st.error("Benchmark review not found in session state.")
                    st.session_state.benchmark_running = False
                    st.rerun()

                # This is the single call to the agent for the batch of 10
                batch_output = screen_abstracts_batch(
                    batch=list(current_conceptual_batch_items),
                    batch_idx=agent_batch_idx,
                    review=st.session_state.benchmark_review,
                )
                st.session_state.current_batch_screening_results = (
                    batch_output.results if batch_output else []
                )
                st.session_state.current_batch_item_offset = (
                    0  # Reset offset for the new batch results
                )

                # If after fetching, there are no results for this batch (e.g. agent error for whole batch)
                if not st.session_state.current_batch_screening_results:
                    logger.error(
                        f"No screening results returned from agent for batch {agent_batch_idx + 1}"
                    )
                    # Advance processed_items_count by the size of this failed batch to avoid getting stuck
                    st.session_state.benchmark_stats["total_processed"] += len(
                        current_conceptual_batch_items
                    )
                    st.session_state.benchmark_stats["screening_errors"] += len(
                        current_conceptual_batch_items
                    )
                    st.session_state.current_batch_screening_results = []  # Ensure it's empty
                    st.session_state.current_batch_item_offset = 0
                    if (
                        st.session_state.benchmark_stats["total_processed"]
                        >= total_items
                    ):
                        st.session_state.benchmark_phase = "completed"
                    st.rerun()  # Rerun to update UI and potentially move to next batch or complete

            # Now, process ONE item from the current_batch_screening_results for UI update
            batch_results = st.session_state.get("current_batch_screening_results", [])
            item_offset = st.session_state.get("current_batch_item_offset", 0)

            if batch_results and item_offset < len(batch_results):
                result_tuple = batch_results[item_offset]
                actual_item_index_in_full_list = batch_start_idx + item_offset
                current_item_for_display = st.session_state.benchmark_search_results[
                    actual_item_index_in_full_list
                ]

                # Update status for the item being processed *after* its batch was screened by AI
                st.session_state.benchmark_status = f"Processing & resolving item {actual_item_index_in_full_list + 1}/{total_items}: {current_item_for_display.title[:50]}..."

                search_result = result_tuple.search_result
                conservative_result = result_tuple.conservative_result
                comprehensive_result = result_tuple.comprehensive_result

                # Process AI logic without database session
                needs_resolver = _needs_resolver(
                    conservative_result, comprehensive_result
                )
                resolver_result_obj = None

                if needs_resolver:
                    st.session_state.benchmark_stats["conflicts_detected"] += 1
                    try:
                        resolver_output_schema = invoke_resolver_chain(
                            search_result=search_result,
                            conservative_result=conservative_result,
                            comprehensive_result=comprehensive_result,
                            review=st.session_state.benchmark_review,
                        )
                        if resolver_output_schema:
                            resolver_result_obj = schemas.ScreeningResult(
                                id=uuid.uuid4(),
                                review_id=st.session_state.benchmark_review.id,
                                search_result_id=search_result.id,
                                trace_id=uuid.uuid4(),
                                model_name="resolver",
                                screening_strategy=schemas.ScreeningStrategyType.COMPREHENSIVE,
                                start_time=datetime.now(timezone.utc),
                                end_time=datetime.now(timezone.utc),
                                decision=resolver_output_schema.resolver_decision,
                                confidence_score=resolver_output_schema.resolver_confidence_score,
                                rationale=resolver_output_schema.resolver_reasoning,
                            )
                            st.session_state.benchmark_stats["resolver_invoked"] += 1
                    except Exception as e_resolve:
                        logger.exception(
                            f"Resolver failed for {search_result.source_id}: {e_resolve}"
                        )

                # Calculate final decision and classification
                final_decision = _determine_final_decision(
                    conservative_result, comprehensive_result, resolver_result_obj
                )
                human_decision_raw = search_result.source_metadata.get(
                    "benchmark_human_decision"
                )
                human_decision = (
                    human_decision_raw if isinstance(human_decision_raw, bool) else None
                )
                classification = _calculate_classification(
                    final_decision, human_decision
                )

                # Now perform database operations in a short-lived session
                db_success = False
                try:
                    with session_factory() as session:
                        benchmark_result_item_repo = BenchmarkResultItemRepository()

                        if isinstance(
                            conservative_result, ScreeningError
                        ) or isinstance(comprehensive_result, ScreeningError):
                            logger.warning(
                                f"Screening error for {search_result.source_id} from agent output."
                            )
                            st.session_state.benchmark_stats["screening_errors"] += 1
                        else:
                            benchmark_result_item = models.BenchmarkResultItem(
                                benchmark_run_id=st.session_state.benchmark_run_id,
                                search_result_id=search_result.id,
                                human_decision=human_decision,
                                conservative_decision=conservative_result.decision,
                                conservative_confidence=conservative_result.confidence_score,
                                conservative_rationale=conservative_result.rationale,
                                conservative_run_id=conservative_result.id,
                                conservative_trace_id=conservative_result.trace_id,
                                comprehensive_decision=comprehensive_result.decision,
                                comprehensive_confidence=comprehensive_result.confidence_score,
                                comprehensive_rationale=comprehensive_result.rationale,
                                comprehensive_run_id=comprehensive_result.id,
                                comprehensive_trace_id=comprehensive_result.trace_id,
                                resolver_decision=resolver_result_obj.decision
                                if resolver_result_obj
                                else None,
                                resolver_confidence=resolver_result_obj.confidence_score
                                if resolver_result_obj
                                else None,
                                resolver_reasoning=resolver_result_obj.rationale
                                if resolver_result_obj
                                else None,
                                final_decision=final_decision,
                                classification=classification,
                            )
                            benchmark_result_item_repo.add(
                                session, benchmark_result_item
                            )

                            result_data = {
                                "human_decision": human_decision,
                                "final_decision": final_decision,
                                "classification": classification,
                                "conservative_decision": conservative_result.decision,
                                "conservative_confidence": conservative_result.confidence_score,
                                "conservative_rationale": conservative_result.rationale,
                                "conservative_quotes": "; ".join(
                                    conservative_result.extracted_quotes
                                )
                                if conservative_result.extracted_quotes
                                else "N/A",
                                "comprehensive_decision": comprehensive_result.decision,
                                "comprehensive_confidence": comprehensive_result.confidence_score,
                                "comprehensive_rationale": comprehensive_result.rationale,
                                "comprehensive_quotes": "; ".join(
                                    comprehensive_result.extracted_quotes
                                )
                                if comprehensive_result.extracted_quotes
                                else "N/A",
                                "resolver_decision": resolver_result_obj.decision
                                if resolver_result_obj
                                else None,
                                "resolver_confidence": resolver_result_obj.confidence_score
                                if resolver_result_obj
                                else None,
                                "resolver_rationale": resolver_result_obj.rationale
                                if resolver_result_obj
                                else None,
                            }
                            st.session_state.benchmark_stats[
                                "accumulated_results"
                            ].append(result_data)

                        session.commit()
                        db_success = True
                        logger.info(
                            f"DB record saved for item {actual_item_index_in_full_list + 1}"
                        )

                    # Update progress and metrics after successful database operation
                    if db_success:
                        st.session_state.benchmark_stats["total_processed"] += 1
                        st.session_state.current_batch_item_offset += 1

                        # Calculate and update real-time metrics
                        accumulated_results = st.session_state.benchmark_stats[
                            "accumulated_results"
                        ]
                        y_true = [r["human_decision"] for r in accumulated_results]
                        y_pred_decision = [
                            r["final_decision"] for r in accumulated_results
                        ]
                        y_pred_confidence = []
                        for r_data in accumulated_results:
                            conf = r_data.get("conservative_confidence", 0.0) or 0.0
                            if r_data.get("comprehensive_confidence") is not None:
                                conf = max(
                                    conf,
                                    r_data.get("comprehensive_confidence", 0.0) or 0.0,
                                )
                            if r_data.get("resolver_confidence") is not None:
                                conf = r_data.get("resolver_confidence", 0.0) or 0.0
                            y_pred_confidence.append(conf)

                        st.session_state.benchmark_current_metrics = calculate_metrics(
                            y_true, y_pred_decision, y_pred_confidence
                        )
                        st.session_state.benchmark_progress = (
                            st.session_state.benchmark_stats["total_processed"]
                            / total_items
                        )
                        logger.info(
                            f"Progress updated for item {actual_item_index_in_full_list + 1}"
                        )

                except Exception as e_item:
                    logger.exception(
                        f"DB/processing for item {actual_item_index_in_full_list + 1} failed: {e_item}"
                    )
                    st.error(
                        f"Processing item {actual_item_index_in_full_list + 1} failed: {e_item}"
                    )
                    st.session_state.benchmark_stats["total_processed"] += (
                        1  # Count as processed to move on
                    )
                    st.session_state.benchmark_stats["screening_errors"] += 1
                    st.session_state.current_batch_item_offset += (
                        1  # Move to next item in batch results
                    )

                # Rerun after each item from the batch is processed and its DB record saved
                st.rerun()

            elif (
                not batch_results
                and processed_items_count < total_items
                and st.session_state.get("current_batch_screening_results") is not None
            ):
                # This case means screen_abstracts_batch was called, returned [], and we already handled advancing past it.
                # We might get here if the previous rerun after handling empty batch_results leads here.
                # Simply rerun to let the main logic decide to fetch the next batch or complete.
                logger.info(
                    "Batch results were empty, rerunning to fetch next batch or complete."
                )
                st.rerun()

            else:  # current_batch_item_offset >= len(batch_results) - All items in current AI batch results processed
                st.session_state.current_batch_screening_results = []  # Clear to trigger fetch of next AI batch
                st.session_state.current_batch_item_offset = 0
                if st.session_state.benchmark_stats["total_processed"] >= total_items:
                    st.session_state.benchmark_phase = "completed"
                st.rerun()  # Rerun to fetch next batch or go to completed phase

        else:  # All items processed (processed_items_count >= total_items)
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

        # Calculate and persist final metrics to database
        with st.spinner("Calculating and persisting final performance metrics..."):
            try:
                with session_factory() as session:
                    benchmark_run_repo = BenchmarkRunRepository()
                    benchmark_result_item_repo = BenchmarkResultItemRepository()

                    updated_run = calculate_and_update_benchmark_metrics(
                        session=session,
                        benchmark_run_id=st.session_state.benchmark_run_id,
                        benchmark_run_repo=benchmark_run_repo,
                        benchmark_result_item_repo=benchmark_result_item_repo,
                    )

                    logger.info(
                        f"Successfully calculated and persisted metrics for benchmark run {updated_run.id}"
                    )
                    st.success(
                        "✅ Final performance metrics calculated and saved to database!"
                    )

            except Exception as e:
                logger.exception("Failed to calculate and persist final metrics")
                st.error(f"Failed to calculate final metrics: {e}")

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

# Display accumulated results for drill-down (this is for the LIVE RUN or if a run just completed)
if st.session_state.get("benchmark_stats", {}).get(
    "total_processed", 0
) > 0 and not st.session_state.get("benchmark_running", False):
    st.subheader("🔎 Last Completed Run Results")

    accumulated_results = st.session_state.benchmark_stats["accumulated_results"]
    search_results = st.session_state.get("benchmark_search_results", [])

    if accumulated_results and search_results:
        st.info(
            "💡 **Tip:** Use the 'View Past Completed Benchmark Runs' section below to access this data with full metrics and export options."
        )
    else:
        st.info("No processed results available to display.")

# --- Re-enable section for loading PAST Completed Benchmark Runs ---
st.markdown("---")
st.header("📊 View Past Completed Benchmark Runs")


# Helper function for fetching all completed runs
def get_all_completed_benchmark_runs() -> list[models.BenchmarkRun]:
    """Fetch all completed benchmark runs with calculated metrics, sorted by newest first."""
    with session_factory() as session:
        try:
            benchmark_run_repo = BenchmarkRunRepository()
            benchmark_result_item_repo = BenchmarkResultItemRepository()
            all_runs = benchmark_run_repo.get_by_review_id(session, BENCHMARK_REVIEW_ID)

            completed_runs = []
            for run in all_runs:
                # Consider a run completed if either:
                # 1. It has calculated metrics (tp or accuracy populated)
                # 2. It has processed items (BenchmarkResultItem records exist)
                has_metrics = run.tp is not None or run.accuracy is not None

                if not has_metrics:
                    # Check if run has processed items
                    result_items = benchmark_result_item_repo.get_by_benchmark_run_id(
                        session, run.id
                    )
                    has_processed_items = len(result_items) > 0
                else:
                    has_processed_items = (
                        False  # Don't need to check if already has metrics
                    )

                if has_metrics or has_processed_items:
                    completed_runs.append(run)

            # Sort by created_at descending (newest first)
            completed_runs.sort(
                key=lambda x: x.created_at or datetime.min, reverse=True
            )
            return completed_runs
        except Exception:
            logger.exception("Failed to fetch completed benchmark runs from DB")
            return []


# Helper function for fetching latest completed run
def get_latest_completed_benchmark_run() -> models.BenchmarkRun | None:
    """Fetch the latest completed benchmark run with calculated metrics."""
    completed_runs = get_all_completed_benchmark_runs()
    return completed_runs[0] if completed_runs else None


# Helper function for formatting metric values (ensure it's defined or imported)
def format_metric_value(
    value: float | None, is_percentage: bool = False, decimal_places: int = 3
) -> str:
    if value is None:
        return "N/A"
    if is_percentage:
        # Ensure value is float for percentage formatting if it's an int
        value_float = float(value)
        return f"{value_float * 100:.{decimal_places - 2}f}%"
    if isinstance(value, float):
        if value == float("inf"):
            return "∞"
        if value == float("-inf"):
            return "-∞"
        return f"{value:.{decimal_places}f}"
    return str(value)  # Handles int directly


# Helper function for loading result items for a specific PAST run
def load_benchmark_result_items_for_past_run(run_id: str) -> pd.DataFrame:
    """Load and format individual benchmark result items for a specific past run."""
    with session_factory() as session:
        try:
            benchmark_result_item_repo = BenchmarkResultItemRepository()
            search_result_repo = SearchResultRepository()
            result_items = benchmark_result_item_repo.get_by_benchmark_run_id(
                session, uuid.UUID(run_id)
            )
            logger.info(f"Loading {len(result_items)} result items for run {run_id}")
            if not result_items:
                logger.warning(f"No result items found for run {run_id}")
                return pd.DataFrame()
            data = []
            for item in result_items:
                search_result = search_result_repo.get_by_id(
                    session, item.search_result_id
                )
                if search_result:
                    row = {
                        "Title": search_result.title[:100]
                        + ("..." if len(search_result.title) > 100 else ""),
                        "Year": search_result.year or "N/A",
                        "Human Decision": "Include"
                        if item.human_decision is True
                        else "Exclude"
                        if item.human_decision is False
                        else "Unknown",
                        "Final Decision": item.final_decision.value,
                        "Classification": item.classification,
                        "Conservative Decision": item.conservative_decision.value
                        if item.conservative_decision
                        else "N/A",
                        "Conservative Confidence": f"{item.conservative_confidence:.3f}"
                        if item.conservative_confidence is not None
                        else "N/A",
                        "Conservative Rationale": item.conservative_rationale[:150]
                        + (
                            "..."
                            if item.conservative_rationale
                            and len(item.conservative_rationale) > 150
                            else ""
                        )
                        if item.conservative_rationale
                        else "N/A",
                        "Comprehensive Decision": item.comprehensive_decision.value
                        if item.comprehensive_decision
                        else "N/A",
                        "Comprehensive Confidence": f"{item.comprehensive_confidence:.3f}"
                        if item.comprehensive_confidence is not None
                        else "N/A",
                        "Comprehensive Rationale": item.comprehensive_rationale[:150]
                        + (
                            "..."
                            if item.comprehensive_rationale
                            and len(item.comprehensive_rationale) > 150
                            else ""
                        )
                        if item.comprehensive_rationale
                        else "N/A",
                        "Resolver Decision": item.resolver_decision.value
                        if item.resolver_decision
                        else "N/A",
                        "Resolver Confidence": f"{item.resolver_confidence:.3f}"
                        if item.resolver_confidence is not None
                        else "N/A",
                        "Resolver Reasoning": item.resolver_reasoning[:150]
                        + (
                            "..."
                            if item.resolver_reasoning
                            and len(item.resolver_reasoning) > 150
                            else ""
                        )
                        if item.resolver_reasoning
                        else "N/A",
                        "Authors": (
                            ", ".join(search_result.authors[:2])
                            + (
                                "..."
                                if search_result.authors
                                and len(search_result.authors) > 2
                                else ""
                            )
                        )
                        if search_result.authors
                        else "N/A",
                        "Source ID": search_result.source_id,
                        "DOI": search_result.doi or "N/A",
                        "_conservative_rationale_full": item.conservative_rationale
                        or "",
                        "_comprehensive_rationale_full": item.comprehensive_rationale
                        or "",
                        "_resolver_rationale_full": item.resolver_reasoning or "",
                        "_conservative_run_id": str(item.conservative_run_id)
                        if item.conservative_run_id
                        else "N/A",
                        "_conservative_trace_id": str(item.conservative_trace_id)
                        if item.conservative_trace_id
                        else "N/A",
                        "_comprehensive_run_id": str(item.comprehensive_run_id)
                        if item.comprehensive_run_id
                        else "N/A",
                        "_comprehensive_trace_id": str(item.comprehensive_trace_id)
                        if item.comprehensive_trace_id
                        else "N/A",
                    }
                    data.append(row)
            return pd.DataFrame(data)
        except Exception as e_load_items:
            logger.exception(
                f"Failed to load benchmark result items for past run {run_id}: {e_load_items}"
            )
            st.error(
                f"Failed to load individual results for run {run_id}: {e_load_items}"
            )
            return pd.DataFrame()


# Initialize session state for run selection
if "selected_past_benchmark_run" not in st.session_state:
    st.session_state.selected_past_benchmark_run = None
if "available_benchmark_runs" not in st.session_state:
    st.session_state.available_benchmark_runs = []

# Load available runs and setup dropdown
with st.spinner("Loading available benchmark runs..."):
    available_runs = get_all_completed_benchmark_runs()
    st.session_state.available_benchmark_runs = available_runs

if available_runs:
    st.subheader("📊 Select Benchmark Run to View")

    # Create options for selectbox with run information
    run_options = []
    for i, run in enumerate(available_runs):
        created_str = (
            run.created_at.strftime("%Y-%m-%d %H:%M:%S")
            if run.created_at
            else "Unknown"
        )

        # Check if metrics are available
        has_metrics = run.accuracy is not None

        if has_metrics:
            accuracy_str = f"Accuracy: {run.accuracy:.3f}"
            status_str = "✅"
        else:
            # For runs without metrics, show number of processed items
            with session_factory() as session:
                benchmark_result_item_repo = BenchmarkResultItemRepository()
                result_items = benchmark_result_item_repo.get_by_benchmark_run_id(
                    session, run.id
                )
                item_count = len(result_items)
            accuracy_str = f"{item_count} items processed"
            status_str = "⏳" if item_count > 0 else "❌"

        run_options.append(f"{status_str} Run #{i + 1}: {created_str} ({accuracy_str})")

    # Add "None" option at the beginning
    run_options.insert(0, "-- Select a benchmark run --")

    selected_run_index = st.selectbox(
        "Choose a completed benchmark run:",
        options=range(len(run_options)),
        format_func=lambda x: run_options[x],
        key="benchmark_run_selector",
        index=0,
    )

    # Update selected run based on dropdown selection
    if selected_run_index > 0:  # Skip the "-- Select --" option
        selected_run = available_runs[selected_run_index - 1]
        st.session_state.selected_past_benchmark_run = selected_run
        logger.info(
            f"User selected benchmark run {selected_run.id} with accuracy: {selected_run.accuracy}"
        )
    else:
        st.session_state.selected_past_benchmark_run = None
        logger.info("User deselected benchmark run")
else:
    st.warning("⚠️ No completed benchmark runs found in the database.")
    st.session_state.selected_past_benchmark_run = None

if st.session_state.selected_past_benchmark_run:
    logger.info(
        f"Displaying results for run {st.session_state.selected_past_benchmark_run.id}"
    )
    st.subheader(
        f"📈 Benchmark Run Results - ID: {st.session_state.selected_past_benchmark_run.id}"
    )
    st.markdown(
        f"**Created:** {st.session_state.selected_past_benchmark_run.created_at}"
    )
    st.markdown(
        f"**Updated:** {st.session_state.selected_past_benchmark_run.updated_at}"
    )

    if st.session_state.selected_past_benchmark_run.run_notes:
        st.markdown(
            f"**Notes:** {st.session_state.selected_past_benchmark_run.run_notes}"
        )
    if st.session_state.selected_past_benchmark_run.config_details:
        with st.expander("🔧 Run Configuration"):
            st.json(st.session_state.selected_past_benchmark_run.config_details)

    # Check if this run has calculated metrics
    has_metrics = st.session_state.selected_past_benchmark_run.accuracy is not None

    if not has_metrics:
        # Check if run has processed items
        with session_factory() as session:
            benchmark_result_item_repo = BenchmarkResultItemRepository()
            result_items = benchmark_result_item_repo.get_by_benchmark_run_id(
                session, st.session_state.selected_past_benchmark_run.id
            )
            item_count = len(result_items)

        if item_count > 0:
            st.warning(
                f"⚠️ This run processed {item_count} items but metrics haven't been calculated yet."
            )

            if st.button("🔄 Calculate Metrics for This Run", type="primary"):
                with st.spinner("Calculating and persisting performance metrics..."):
                    try:
                        logger.info(
                            f"Starting metrics calculation for benchmark run {st.session_state.selected_past_benchmark_run.id}"
                        )

                        with session_factory() as session:
                            benchmark_run_repo = BenchmarkRunRepository()
                            benchmark_result_item_repo = BenchmarkResultItemRepository()

                            # Debug: Check if run exists and has items
                            run_check = benchmark_run_repo.get_by_id(
                                session, st.session_state.selected_past_benchmark_run.id
                            )
                            items_check = (
                                benchmark_result_item_repo.get_by_benchmark_run_id(
                                    session,
                                    st.session_state.selected_past_benchmark_run.id,
                                )
                            )

                            logger.info(
                                f"Run check: {run_check is not None}, Items found: {len(items_check) if items_check else 0}"
                            )

                            if not run_check:
                                st.error(
                                    f"❌ Benchmark run {st.session_state.selected_past_benchmark_run.id} not found in database!"
                                )
                                logger.error(
                                    f"Benchmark run {st.session_state.selected_past_benchmark_run.id} not found"
                                )
                                st.stop()

                            if not items_check:
                                st.error("❌ No result items found for this run!")
                                logger.error(
                                    f"No result items found for run {st.session_state.selected_past_benchmark_run.id}"
                                )
                                st.stop()

                            st.info(
                                f"Found {len(items_check)} result items. Calculating metrics..."
                            )

                            updated_run = calculate_and_update_benchmark_metrics(
                                session=session,
                                benchmark_run_id=st.session_state.selected_past_benchmark_run.id,
                                benchmark_run_repo=benchmark_run_repo,
                                benchmark_result_item_repo=benchmark_result_item_repo,
                            )

                            if updated_run:
                                # Update the selected run in session state
                                st.session_state.selected_past_benchmark_run = (
                                    updated_run
                                )
                                has_metrics = True

                                logger.info(
                                    f"Successfully calculated metrics for benchmark run {updated_run.id}"
                                )
                                logger.info(
                                    f"Calculated metrics: Accuracy={updated_run.accuracy}, TP={updated_run.tp}, TN={updated_run.tn}, FP={updated_run.fp}, FN={updated_run.fn}"
                                )

                                st.success(
                                    "✅ Metrics calculated and saved successfully!"
                                )
                                st.info(
                                    f"📊 New metrics: Accuracy={updated_run.accuracy:.3f}, TP={updated_run.tp}, TN={updated_run.tn}, FP={updated_run.fp}, FN={updated_run.fn}"
                                )
                                st.rerun()
                            else:
                                st.error(
                                    "❌ Metrics calculation returned None - check logs for details"
                                )
                                logger.error(
                                    "calculate_and_update_benchmark_metrics returned None"
                                )

                    except Exception as e:
                        logger.exception(
                            f"Failed to calculate metrics for selected run {st.session_state.selected_past_benchmark_run.id}"
                        )
                        st.error(f"Failed to calculate metrics: {e!s}")
                        st.error("Check the logs for detailed error information.")

                        # Show the full traceback in debug mode
                        import traceback

                        st.code(traceback.format_exc(), language="python")
        else:
            st.error(
                "❌ This run has no processed items. It may have failed during execution."
            )

    if has_metrics:
        # Display all metrics for the past run
        st.subheader("🔢 Confusion Matrix Counts (Past Run)")
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        with m_col1:
            st.metric(
                "True Positives (TP)",
                format_metric_value(st.session_state.selected_past_benchmark_run.tp),
            )
        with m_col2:
            st.metric(
                "False Positives (FP)",
                format_metric_value(st.session_state.selected_past_benchmark_run.fp),
            )
        with m_col3:
            st.metric(
                "False Negatives (FN)",
                format_metric_value(st.session_state.selected_past_benchmark_run.fn),
            )
        with m_col4:
            st.metric(
                "True Negatives (TN)",
                format_metric_value(st.session_state.selected_past_benchmark_run.tn),
            )

        st.subheader("📊 Primary Performance Metrics (Past Run)")
        pm_col1, pm_col2, pm_col3 = st.columns(3)
        with pm_col1:
            st.metric(
                "Accuracy",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.accuracy,
                    is_percentage=True,
                ),
            )
            st.metric(
                "Sensitivity (Recall)",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.sensitivity,
                    is_percentage=True,
                ),
            )
        with pm_col2:
            st.metric(
                "Specificity",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.specificity,
                    is_percentage=True,
                ),
            )
            st.metric(
                "Precision (PPV)",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.ppv, is_percentage=True
                ),
            )
        with pm_col3:
            st.metric(
                "NPV",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.npv, is_percentage=True
                ),
            )
            st.metric(
                "F1 Score",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.f1_score,
                    is_percentage=True,
                ),
            )

        st.subheader("🎯 Advanced Performance Metrics (Past Run)")
        am_col1, am_col2, am_col3 = st.columns(3)
        with am_col1:
            st.metric(
                "Matthews Corr. Coeff. (MCC)",
                format_metric_value(st.session_state.selected_past_benchmark_run.mcc),
            )
            st.metric(
                "Cohen's Kappa",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.cohen_kappa
                ),
            )
        with am_col2:
            st.metric(
                "PABAK",
                format_metric_value(st.session_state.selected_past_benchmark_run.pabak),
            )
            st.metric(
                "Positive Likelihood Ratio (LR+)",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.lr_plus
                ),
            )
        with am_col3:
            st.metric(
                "Negative Likelihood Ratio (LR-)",
                format_metric_value(
                    st.session_state.selected_past_benchmark_run.lr_minus
                ),
            )

        st.subheader("📋 Individual Paper-Level Results (Past Run)")
        past_run_results_df = load_benchmark_result_items_for_past_run(
            str(st.session_state.selected_past_benchmark_run.id)
        )

        if not past_run_results_df.empty:
            filter_col1, filter_col2, filter_col3 = st.columns(3)
            with filter_col1:
                classification_options = ["All"] + sorted(
                    past_run_results_df["Classification"].unique().tolist()
                )
                selected_classification_past = st.selectbox(
                    "Filter by Classification:",
                    options=classification_options,
                    key="classification_filter_past_run",
                )
            with filter_col2:
                show_mismatches_only = st.checkbox(
                    "Show only mismatches (FP/FN)", key="mismatches_filter_past_run"
                )
            with filter_col3:
                # Export buttons for past run
                past_export_col1, past_export_col2 = st.columns(2)
                with past_export_col1:
                    st.download_button(
                        label="📄 TSV",
                        data=past_run_results_df.to_csv(sep="\t", index=False),
                        file_name=f"past_benchmark_results_{st.session_state.selected_past_benchmark_run.id}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.tsv",
                        mime="text/tab-separated-values",
                        key="export_past_results_tsv",
                    )
                with past_export_col2:
                    # Convert DataFrame to Excel
                    past_excel_buffer = io.BytesIO()
                    with pd.ExcelWriter(past_excel_buffer, engine="openpyxl") as writer:
                        past_run_results_df.to_excel(
                            writer, index=False, sheet_name="Past Run Results"
                        )
                    past_excel_data = past_excel_buffer.getvalue()

                    st.download_button(
                        label="📊 XLSX",
                        data=past_excel_data,
                        file_name=f"past_benchmark_results_{st.session_state.selected_past_benchmark_run.id}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        key="export_past_results_xlsx",
                    )

            # Apply filters
            df_to_display = past_run_results_df
            if selected_classification_past != "All":
                df_to_display = df_to_display[
                    df_to_display["Classification"] == selected_classification_past
                ]
            if show_mismatches_only:
                df_to_display = df_to_display[
                    df_to_display["Classification"].isin(["FP", "FN"])
                ]

            st.info(
                f"Showing {len(df_to_display)} of {len(past_run_results_df)} results from past run"
            )
            display_columns_past = [
                col for col in df_to_display.columns if not col.startswith("_")
            ]
            st.dataframe(
                df_to_display[display_columns_past],
                use_container_width=True,
                hide_index=True,
                height=400,
            )

            if not df_to_display.empty:
                st.subheader("🔍 Paper Detail View (Past Run)")
                # Ensure options are generated based on the *currently displayed* DataFrame (df_to_display)
                # and that its index is reset if filtering changes the original indices.
                df_for_selection = df_to_display.reset_index(drop=True)
                selected_paper_idx_past = st.selectbox(
                    "Select paper from displayed results to view details:",  # Clarified label
                    options=range(len(df_for_selection)),
                    format_func=lambda x: f"{df_for_selection.iloc[x]['Title'][:50]}... ({df_for_selection.iloc[x]['Source ID']})",
                    key="selected_paper_detail_past_run",
                    index=None,  # Default to no selection or 0 if you prefer first item selected
                )
                if selected_paper_idx_past is not None:
                    selected_paper_past = df_for_selection.iloc[selected_paper_idx_past]
                    # Display details for selected_paper_past
                    detail_col1, detail_col2 = st.columns(2)
                    with detail_col1:
                        st.markdown("**Conservative Agent:**")
                        st.markdown(
                            f"Decision: {selected_paper_past['Conservative Decision']}"
                        )
                        st.markdown(
                            f"Confidence: {selected_paper_past['Conservative Confidence']}"
                        )
                        # LangSmith IDs for Conservative Agent
                        if selected_paper_past["_conservative_rationale_full"]:
                            st.markdown(
                                f"**LangSmith Run ID:** `{selected_paper_past['_conservative_rationale_full']}`"
                            )
                            if selected_paper_past["_conservative_quotes_full"]:
                                st.markdown(
                                    f"**Quotes:** {selected_paper_past['_conservative_quotes_full']}"
                                )
                        if selected_paper_past["_conservative_rationale_full"]:
                            with st.expander("Rationale"):
                                st.markdown(
                                    selected_paper_past["_conservative_rationale_full"]
                                )
                        # LangSmith IDs for Conservative Agent
                        if selected_paper_past["_conservative_run_id"] != "N/A":
                            st.markdown(
                                f"**LangSmith Run ID:** `{selected_paper_past['_conservative_run_id']}`"
                            )
                            if selected_paper_past["_conservative_trace_id"] != "N/A":
                                langsmith_url = f"https://smith.langchain.com/o/00000000-0000-0000-0000-000000000000/projects/p/r/{selected_paper_past['_conservative_run_id']}"
                                st.markdown(
                                    f"**Trace ID:** [`{selected_paper_past['_conservative_trace_id']}`]({langsmith_url})"
                                )
                    with detail_col2:
                        st.markdown("**Comprehensive Agent:**")
                        st.markdown(
                            f"Decision: {selected_paper_past['Comprehensive Decision']}"
                        )
                        st.markdown(
                            f"Confidence: {selected_paper_past['Comprehensive Confidence']}"
                        )
                        # LangSmith IDs for Comprehensive Agent
                        if selected_paper_past["_comprehensive_rationale_full"]:
                            st.markdown(
                                f"**LangSmith Run ID:** `{selected_paper_past['_comprehensive_rationale_full']}`"
                            )
                            if selected_paper_past["_comprehensive_quotes_full"]:
                                st.markdown(
                                    f"**Quotes:** {selected_paper_past['_comprehensive_quotes_full']}"
                                )
                        if selected_paper_past["_comprehensive_rationale_full"]:
                            with st.expander("Rationale"):
                                st.markdown(
                                    selected_paper_past["_comprehensive_rationale_full"]
                                )
                        # LangSmith IDs for Comprehensive Agent
                        if selected_paper_past["_comprehensive_run_id"] != "N/A":
                            st.markdown(
                                f"**LangSmith Run ID:** `{selected_paper_past['_comprehensive_run_id']}`"
                            )
                            if selected_paper_past["_comprehensive_trace_id"] != "N/A":
                                langsmith_url = f"https://smith.langchain.com/o/00000000-0000-0000-0000-000000000000/projects/p/r/{selected_paper_past['_comprehensive_run_id']}"
                                st.markdown(
                                    f"**Trace ID:** [`{selected_paper_past['_comprehensive_trace_id']}`]({langsmith_url})"
                                )
                        if selected_paper_past["_comprehensive_rationale_full"]:
                            with st.expander("Rationale"):
                                st.markdown(
                                    selected_paper_past["_comprehensive_rationale_full"]
                                )
                    if selected_paper_past["Resolver Decision"] != "N/A":
                        st.markdown("**Resolver Agent:**")
                        st.markdown(
                            f"Decision: {selected_paper_past['Resolver Decision']}"
                        )
                        st.markdown(
                            f"Confidence: {selected_paper_past['Resolver Confidence']}"
                        )
                        if selected_paper_past["_resolver_rationale_full"]:
                            with st.expander("Reasoning"):
                                st.markdown(
                                    selected_paper_past["_resolver_rationale_full"]
                                )
                    st.markdown("---")
                    st.markdown(
                        f"**Human Ground Truth:** {selected_paper_past['Human Decision']}"
                    )
                    st.markdown(
                        f"**SRA Final Decision:** {selected_paper_past['Final Decision']}"
                    )
                    st.markdown(
                        f"**Classification:** {selected_paper_past['Classification']}"
                    )
        else:
            st.warning("No individual results found for this past benchmark run.")

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
