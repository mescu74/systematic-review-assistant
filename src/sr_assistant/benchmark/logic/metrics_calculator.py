"""Performance metrics calculator for systematic review benchmarking.

This module provides functions to calculate standard screening performance metrics
based on True Positives (TP), False Positives (FP), True Negatives (TN), and
False Negatives (FN) counts, following the formulas defined in docs/sr_metrics.md.
"""

from __future__ import annotations

import math
import typing as t

from loguru import logger

from sr_assistant.core import models, schemas
from sr_assistant.core.types import ScreeningDecisionType

if t.TYPE_CHECKING:
    import uuid
    from collections.abc import Sequence


def _classify_decision(
    final_decision: ScreeningDecisionType,
    human_decision: bool | None,
) -> str:
    """Classify a single decision as TP, FP, TN, FN, or UNKNOWN.

    Args:
        final_decision: The SRA's final decision for the item
        human_decision: The ground truth decision (True=include, False=exclude, None=unknown)

    Returns:
        Classification string: "TP", "FP", "TN", "FN", or "UNKNOWN"
    """
    if human_decision is None:
        return "UNKNOWN"

    # Convert final_decision to boolean for comparison
    # UNCERTAIN is treated as incorrect (either FP or FN depending on ground truth)
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


def calculate_confusion_matrix_counts(
    result_items: Sequence[models.BenchmarkResultItem],
) -> dict[str, int]:
    """Calculate TP, FP, TN, FN counts from benchmark result items.

    Args:
        result_items: List of BenchmarkResultItem records for a benchmark run

    Returns:
        Dictionary with counts: {"tp": int, "fp": int, "tn": int, "fn": int, "total_compared": int}
    """
    tp = fp = tn = fn = 0

    for item in result_items:
        classification = _classify_decision(item.final_decision, item.human_decision)

        if classification == "TP":
            tp += 1
        elif classification == "FP":
            fp += 1
        elif classification == "TN":
            tn += 1
        elif classification == "FN":
            fn += 1
        # UNKNOWN items are not counted in confusion matrix

    total_compared = tp + fp + tn + fn

    logger.debug(
        f"Confusion matrix counts: TP={tp}, FP={fp}, TN={tn}, FN={fn}, Total={total_compared}"
    )

    return {
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "total_compared": total_compared,
    }


def calculate_sensitivity(tp: int, fn: int) -> float | None:
    """Calculate Sensitivity (Recall).

    Formula: Sensitivity = TP / (TP + FN)

    Args:
        tp: True Positives count
        fn: False Negatives count

    Returns:
        Sensitivity value [0.0, 1.0] or None if undefined (TP + FN = 0)
    """
    denominator = tp + fn
    if denominator == 0:
        logger.warning("Sensitivity undefined: TP + FN = 0")
        return None
    return float(tp) / denominator


def calculate_specificity(tn: int, fp: int) -> float | None:
    """Calculate Specificity.

    Formula: Specificity = TN / (TN + FP)

    Args:
        tn: True Negatives count
        fp: False Positives count

    Returns:
        Specificity value [0.0, 1.0] or None if undefined (TN + FP = 0)
    """
    denominator = tn + fp
    if denominator == 0:
        logger.warning("Specificity undefined: TN + FP = 0")
        return None
    return float(tn) / denominator


def calculate_accuracy(tp: int, tn: int, fp: int, fn: int) -> float | None:
    """Calculate Accuracy.

    Formula: Accuracy = (TP + TN) / (TP + TN + FP + FN)

    Args:
        tp: True Positives count
        tn: True Negatives count
        fp: False Positives count
        fn: False Negatives count

    Returns:
        Accuracy value [0.0, 1.0] or None if undefined (total = 0)
    """
    total = tp + tn + fp + fn
    if total == 0:
        logger.warning("Accuracy undefined: total count = 0")
        return None
    return float(tp + tn) / total


def calculate_ppv(tp: int, fp: int) -> float | None:
    """Calculate Positive Predictive Value (Precision).

    Formula: PPV = TP / (TP + FP)

    Args:
        tp: True Positives count
        fp: False Positives count

    Returns:
        PPV value [0.0, 1.0] or None if undefined (TP + FP = 0)
    """
    denominator = tp + fp
    if denominator == 0:
        logger.warning("PPV undefined: TP + FP = 0")
        return None
    return float(tp) / denominator


def calculate_npv(tn: int, fn: int) -> float | None:
    """Calculate Negative Predictive Value.

    Formula: NPV = TN / (TN + FN)

    Args:
        tn: True Negatives count
        fn: False Negatives count

    Returns:
        NPV value [0.0, 1.0] or None if undefined (TN + FN = 0)
    """
    denominator = tn + fn
    if denominator == 0:
        logger.warning("NPV undefined: TN + FN = 0")
        return None
    return float(tn) / denominator


def calculate_f1_score(tp: int, fp: int, fn: int) -> float | None:
    """Calculate F1 Score.

    Formula: F1 = 2 * (Precision * Recall) / (Precision + Recall)
    Which simplifies to: F1 = 2 * TP / (2 * TP + FP + FN)

    Args:
        tp: True Positives count
        fp: False Positives count
        fn: False Negatives count

    Returns:
        F1 Score value [0.0, 1.0] or None if undefined (2*TP + FP + FN = 0)
    """
    denominator = 2 * tp + fp + fn
    if denominator == 0:
        logger.warning("F1 Score undefined: 2*TP + FP + FN = 0")
        return None
    return float(2 * tp) / denominator


def calculate_mcc(tp: int, tn: int, fp: int, fn: int) -> float | None:
    """Calculate Matthews Correlation Coefficient.

    Formula: MCC = (TP * TN - FP * FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))

    Args:
        tp: True Positives count
        tn: True Negatives count
        fp: False Positives count
        fn: False Negatives count

    Returns:
        MCC value [-1.0, 1.0] or None if undefined (denominator = 0)
    """
    numerator = tp * tn - fp * fn
    denominator_terms = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)

    if denominator_terms == 0:
        logger.warning("MCC undefined: denominator = 0")
        return None

    return float(numerator) / math.sqrt(denominator_terms)


def calculate_cohen_kappa(tp: int, tn: int, fp: int, fn: int) -> float | None:
    """Calculate Cohen's Kappa.

    Formula: κ = (p_o - p_e) / (1 - p_e)
    Where:
    - p_o = observed agreement = (TP + TN) / Total
    - p_e = expected agreement = ((TP+FP)(TP+FN) + (TN+FP)(TN+FN)) / Total²

    Args:
        tp: True Positives count
        tn: True Negatives count
        fp: False Positives count
        fn: False Negatives count

    Returns:
        Cohen's Kappa value [-1.0, 1.0] or None if undefined
    """
    total = tp + tn + fp + fn
    if total == 0:
        logger.warning("Cohen's Kappa undefined: total = 0")
        return None

    # Observed agreement
    p_o = (tp + tn) / total

    # Expected agreement by chance
    p_e = ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (total * total)

    if p_e == 1.0:
        logger.warning("Cohen's Kappa undefined: expected agreement = 1.0")
        return None

    return (p_o - p_e) / (1 - p_e)


def calculate_pabak(tp: int, tn: int, fp: int, fn: int) -> float | None:
    """Calculate Prevalence and Bias Adjusted Kappa (PABAK).

    Formula: PABAK = 2 * p_o - 1
    Where p_o = observed agreement = (TP + TN) / Total

    Args:
        tp: True Positives count
        tn: True Negatives count
        fp: False Positives count
        fn: False Negatives count

    Returns:
        PABAK value [-1.0, 1.0] or None if undefined (total = 0)
    """
    total = tp + tn + fp + fn
    if total == 0:
        logger.warning("PABAK undefined: total = 0")
        return None

    p_o = (tp + tn) / total
    return 2 * p_o - 1


def calculate_lr_plus(
    sensitivity: float | None, specificity: float | None
) -> float | None:
    """Calculate Positive Likelihood Ratio.

    Formula: LR+ = Sensitivity / (1 - Specificity)

    Args:
        sensitivity: Sensitivity value
        specificity: Specificity value

    Returns:
        LR+ value or None if undefined (specificity = 1.0 or inputs are None)
    """
    if sensitivity is None or specificity is None:
        logger.warning("LR+ undefined: sensitivity or specificity is None")
        return None

    if specificity == 1.0:
        logger.warning("LR+ undefined: specificity = 1.0 (division by zero)")
        return None

    return sensitivity / (1 - specificity)


def calculate_lr_minus(
    sensitivity: float | None, specificity: float | None
) -> float | None:
    """Calculate Negative Likelihood Ratio.

    Formula: LR- = (1 - Sensitivity) / Specificity

    Args:
        sensitivity: Sensitivity value
        specificity: Specificity value

    Returns:
        LR- value or None if undefined (specificity = 0.0 or inputs are None)
    """
    if sensitivity is None or specificity is None:
        logger.warning("LR- undefined: sensitivity or specificity is None")
        return None

    if specificity == 0.0:
        logger.warning("LR- undefined: specificity = 0.0 (division by zero)")
        return None

    return (1 - sensitivity) / specificity


def calculate_all_metrics(
    result_items: Sequence[models.BenchmarkResultItem],
) -> schemas.BenchmarkRunUpdate:
    """Calculate all performance metrics for a benchmark run.

    This is the main function that calculates all metrics and returns a
    BenchmarkRunUpdate schema that can be used to update the BenchmarkRun database record.

    Args:
        result_items: List of BenchmarkResultItem records for the benchmark run

    Returns:
        BenchmarkRunUpdate schema populated with all calculated metrics
    """
    logger.info(f"Calculating metrics for {len(result_items)} benchmark result items")

    # Calculate confusion matrix counts
    counts = calculate_confusion_matrix_counts(result_items)
    tp, fp, tn, fn = counts["tp"], counts["fp"], counts["tn"], counts["fn"]

    logger.info(f"Confusion matrix - TP: {tp}, FP: {fp}, TN: {tn}, FN: {fn}")

    # Calculate primary metrics
    sensitivity = calculate_sensitivity(tp, fn)
    specificity = calculate_specificity(tn, fp)
    accuracy = calculate_accuracy(tp, tn, fp, fn)
    ppv = calculate_ppv(tp, fp)
    npv = calculate_npv(tn, fn)

    # Calculate composite metrics
    f1_score = calculate_f1_score(tp, fp, fn)
    mcc = calculate_mcc(tp, tn, fp, fn)
    cohen_kappa = calculate_cohen_kappa(tp, tn, fp, fn)
    pabak = calculate_pabak(tp, tn, fp, fn)

    # Calculate likelihood ratios
    lr_plus = calculate_lr_plus(sensitivity, specificity)
    lr_minus = calculate_lr_minus(sensitivity, specificity)

    # Create the update schema
    metrics_update = schemas.BenchmarkRunUpdate(
        tp=tp,
        fp=fp,
        tn=tn,
        fn=fn,
        sensitivity=sensitivity,
        specificity=specificity,
        accuracy=accuracy,
        ppv=ppv,
        npv=npv,
        f1_score=f1_score,
        mcc=mcc,
        cohen_kappa=cohen_kappa,
        pabak=pabak,
        lr_plus=lr_plus,
        lr_minus=lr_minus,
    )

    logger.info("Successfully calculated all performance metrics")
    logger.debug(
        f"Metrics summary: Accuracy={accuracy or 'N/A'}, F1={f1_score or 'N/A'}, MCC={mcc or 'N/A'}"
    )

    return metrics_update


def calculate_and_update_benchmark_metrics(
    session: t.Any,  # SQLModel Session
    benchmark_run_id: uuid.UUID,
    benchmark_run_repo: t.Any,  # BenchmarkRunRepository
    benchmark_result_item_repo: t.Any,  # BenchmarkResultItemRepository
) -> models.BenchmarkRun:
    """Calculate metrics for a benchmark run and update the database record.

    This function fetches all BenchmarkResultItem records for a run, calculates metrics,
    and updates the BenchmarkRun record with the results.

    Args:
        session: Database session
        benchmark_run_id: ID of the benchmark run to calculate metrics for
        benchmark_run_repo: BenchmarkRunRepository instance
        benchmark_result_item_repo: BenchmarkResultItemRepository instance

    Returns:
        Updated BenchmarkRun model instance

    Raises:
        ValueError: If benchmark run is not found
        Exception: If database operations fail
    """
    logger.info(f"Starting metrics calculation for benchmark run {benchmark_run_id}")

    # Fetch the benchmark run
    benchmark_run = benchmark_run_repo.get_by_id(session, benchmark_run_id)
    if not benchmark_run:
        msg = f"BenchmarkRun with ID {benchmark_run_id} not found"
        logger.error(msg)
        raise ValueError(msg)

    # Fetch all result items for this run
    result_items = benchmark_result_item_repo.get_by_benchmark_run_id(
        session, benchmark_run_id
    )

    if not result_items:
        logger.warning(f"No result items found for benchmark run {benchmark_run_id}")
        # Still update with zero counts
        result_items = []

    # Calculate all metrics
    metrics_update = calculate_all_metrics(result_items)

    # Update the benchmark run fields
    for field_name, value in metrics_update.model_dump(exclude_unset=True).items():
        if hasattr(benchmark_run, field_name):
            setattr(benchmark_run, field_name, value)

    # Save the updated benchmark run
    updated_run = benchmark_run_repo.update(session, benchmark_run)

    logger.info(
        f"Successfully updated benchmark run {benchmark_run_id} with calculated metrics"
    )

    return updated_run
