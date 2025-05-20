# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

import collections.abc
import io
import typing as t
import uuid
from datetime import datetime
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
    ScreeningResult,
    screen_abstracts_batch,
)
from sr_assistant.app.database import session_factory
from sr_assistant.core import models
from sr_assistant.core.repositories import (
    SearchResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.core.schemas import ScreeningDecisionType

# Fixed UUID for the benchmark review (must match the one in seed_benchmark_data.py)
BENCHMARK_REVIEW_ID = uuid.UUID("00000000-1111-2222-3333-444444444444")


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
            y_true_filtered, y_pred_bool, zero_division=0
        )
        metrics["Specificity"] = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
        metrics["Precision (PPV)"] = precision_score(
            y_true_filtered, y_pred_bool, zero_division=0
        )
        metrics["Negative Predictive Value (NPV)"] = (
            float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0
        )
        metrics["Accuracy"] = accuracy_score(y_true_filtered, y_pred_bool)
        metrics["F1 Score"] = f1_score(y_true_filtered, y_pred_bool, zero_division=0)

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


st.set_page_config(layout="wide", page_title="SRA Benchmark Tool")
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
        with st.spinner("Running AI screening on benchmark dataset..."):
            benchmark_review_model: models.SystematicReview = (
                st.session_state.benchmark_review
            )
            search_results_to_process: list[models.SearchResult] = list(
                st.session_state.benchmark_search_results
            )

            logger.info(
                f"Starting benchmark screening for review ID: {benchmark_review_model.id} with {len(search_results_to_process)} abstracts."
            )
            try:
                batch_output = screen_abstracts_batch(
                    batch=search_results_to_process,
                    batch_idx=0,
                    review=benchmark_review_model,
                )
            except Exception:
                st.error("An error occurred while calling screen_abstracts_batch.")
                logger.exception("Error during screen_abstracts_batch")
                st.stop()

            if not batch_output or not batch_output.results:
                st.error("AI screening did not return any output or results.")
                logger.error(
                    "screen_abstracts_batch returned None, empty, or no results list."
                )
                st.stop()

            ai_decisions = []
            comparison_data = []
            num_ai_errors = 0

            for i, result_tuple in enumerate(batch_output.results):
                search_result_model = result_tuple.search_result
                conservative_result = result_tuple.conservative_result

                ai_decision: ScreeningDecisionType | None = None
                ai_confidence: float | None = None
                ai_rationale = "N/A"

                if isinstance(conservative_result, ScreeningResult):
                    ai_decision = conservative_result.decision
                    ai_confidence = conservative_result.confidence_score
                    ai_rationale = conservative_result.rationale
                elif type(conservative_result) is ScreeningError:
                    logger.warning(
                        f"ScreeningError for {search_result_model.source_id}: {conservative_result.error}"
                    )
                    ai_decision = ScreeningDecisionType.UNCERTAIN
                    ai_rationale = f"Screening Error: {conservative_result.error}"
                    num_ai_errors += 1
                else:
                    logger.error(
                        f"Unexpected result type for {search_result_model.source_id}: {type(conservative_result)}"
                    )
                    ai_decision = ScreeningDecisionType.UNCERTAIN
                    ai_rationale = "Unexpected AI result type"
                    num_ai_errors += 1

                ai_decisions.append(
                    {
                        "search_result_id": search_result_model.id,
                        "ai_decision": ai_decision,
                        "ai_confidence": ai_confidence,
                        "ai_rationale": ai_rationale,
                    }
                )

                human_decision_raw = search_result_model.source_metadata.get(
                    "benchmark_human_decision"
                )
                human_decision_bool: bool | None = None
                if isinstance(human_decision_raw, bool):
                    human_decision_bool = human_decision_raw
                elif human_decision_raw is not None:
                    logger.warning(
                        f"Human decision for {search_result_model.source_id} is not bool or None: {human_decision_raw}. Treating as UNKNOWN/None."
                    )

                comparison_data.append(
                    {
                        "search_result_source_id": search_result_model.source_id,
                        "title": search_result_model.title,
                        "ai_decision": ai_decision.value if ai_decision else "ERROR",
                        "human_decision": "INCLUDE"
                        if human_decision_bool is True
                        else ("EXCLUDE" if human_decision_bool is False else "UNKNOWN"),
                        "classification": "N/A",
                    }
                )

            st.session_state.benchmark_ai_decisions = ai_decisions
            st.session_state.benchmark_comparison_data = comparison_data
            logger.info(
                f"AI screening complete. {len(ai_decisions)} decisions processed. {num_ai_errors} AI errors."
            )
            st.success(f"AI Screening complete. {num_ai_errors} AI errors encountered.")

            y_true_for_metrics: list[bool | None] = []
            for sr_model_item in search_results_to_process:
                human_decision_val = sr_model_item.source_metadata.get(
                    "benchmark_human_decision"
                )
                if isinstance(human_decision_val, bool):
                    y_true_for_metrics.append(human_decision_val)
                else:
                    y_true_for_metrics.append(None)

            y_pred_decision_for_metrics: list[ScreeningDecisionType | None] = [
                decision_info["ai_decision"] for decision_info in ai_decisions
            ]
            y_pred_confidence_for_metrics: list[float | None] = [
                decision_info["ai_confidence"] for decision_info in ai_decisions
            ]

            st.session_state.benchmark_metrics = calculate_metrics(
                y_true_for_metrics,
                y_pred_decision_for_metrics,
                y_pred_confidence_for_metrics,
            )

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
        if ai == "INCLUDE" and human == "INCLUDE":
            return "TP"
        if ai == "EXCLUDE" and human == "EXCLUDE":
            return "TN"
        if ai == "INCLUDE" and human == "EXCLUDE":
            return "FP"
        if ai == "EXCLUDE" and human == "INCLUDE":
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
        file_name=f"benchmark_comparison_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
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
            file_name=f"benchmark_summary_metrics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv",
            key="download_metrics_csv",
        )
