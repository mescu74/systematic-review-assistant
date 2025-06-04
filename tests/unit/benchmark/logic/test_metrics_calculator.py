# tests/unit/benchmark/logic/test_metrics_calculator.py
"""Unit tests for metrics calculator module."""

from __future__ import annotations

import math
import uuid
from unittest.mock import MagicMock

import pytest
from pytest_mock import MockerFixture

from sr_assistant.benchmark.logic.metrics_calculator import (
    calculate_accuracy,
    calculate_all_metrics,
    calculate_and_update_benchmark_metrics,
    calculate_cohen_kappa,
    calculate_confusion_matrix_counts,
    calculate_f1_score,
    calculate_lr_minus,
    calculate_lr_plus,
    calculate_mcc,
    calculate_npv,
    calculate_pabak,
    calculate_ppv,
    calculate_sensitivity,
    calculate_specificity,
)
from sr_assistant.core import schemas
from sr_assistant.core.types import ScreeningDecisionType


class TestConfusionMatrixCounts:
    """Test the calculate_confusion_matrix_counts function."""

    def test_confusion_matrix_counts_mixed(self) -> None:
        """Test confusion matrix calculation with mixed results."""
        # Create mock BenchmarkResultItem objects
        items = [
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=True
            ),  # TP
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=False
            ),  # FP
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=True
            ),  # FN
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=False
            ),  # TN
            MagicMock(
                final_decision=ScreeningDecisionType.UNCERTAIN, human_decision=True
            ),  # FN
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=None
            ),  # UNKNOWN
        ]

        result = calculate_confusion_matrix_counts(items)

        assert result["tp"] == 1
        assert result["fp"] == 1
        assert result["tn"] == 1
        assert result["fn"] == 2  # One EXCLUDE with True, one UNCERTAIN with True
        assert result["total_compared"] == 5  # UNKNOWN not counted

    def test_confusion_matrix_counts_empty(self) -> None:
        """Test confusion matrix calculation with empty list."""
        result = calculate_confusion_matrix_counts([])

        assert result["tp"] == 0
        assert result["fp"] == 0
        assert result["tn"] == 0
        assert result["fn"] == 0
        assert result["total_compared"] == 0

    def test_confusion_matrix_counts_all_unknown(self) -> None:
        """Test confusion matrix calculation with all unknown human decisions."""
        items = [
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=None
            ),
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=None
            ),
        ]

        result = calculate_confusion_matrix_counts(items)

        assert result["tp"] == 0
        assert result["fp"] == 0
        assert result["tn"] == 0
        assert result["fn"] == 0
        assert result["total_compared"] == 0


class TestSensitivity:
    """Test the calculate_sensitivity function."""

    def test_sensitivity_normal(self) -> None:
        """Test sensitivity calculation with normal values."""
        result = calculate_sensitivity(tp=8, fn=2)
        assert result == pytest.approx(0.8)

    def test_sensitivity_perfect(self) -> None:
        """Test sensitivity calculation with perfect score."""
        result = calculate_sensitivity(tp=10, fn=0)
        assert result == pytest.approx(1.0)

    def test_sensitivity_zero_tp(self) -> None:
        """Test sensitivity calculation with zero TP."""
        result = calculate_sensitivity(tp=0, fn=5)
        assert result == pytest.approx(0.0)

    def test_sensitivity_undefined(self) -> None:
        """Test sensitivity calculation when undefined (TP + FN = 0)."""
        result = calculate_sensitivity(tp=0, fn=0)
        assert result is None


class TestSpecificity:
    """Test the calculate_specificity function."""

    def test_specificity_normal(self) -> None:
        """Test specificity calculation with normal values."""
        result = calculate_specificity(tn=7, fp=3)
        assert result == pytest.approx(0.7)

    def test_specificity_perfect(self) -> None:
        """Test specificity calculation with perfect score."""
        result = calculate_specificity(tn=10, fp=0)
        assert result == pytest.approx(1.0)

    def test_specificity_zero_tn(self) -> None:
        """Test specificity calculation with zero TN."""
        result = calculate_specificity(tn=0, fp=5)
        assert result == pytest.approx(0.0)

    def test_specificity_undefined(self) -> None:
        """Test specificity calculation when undefined (TN + FP = 0)."""
        result = calculate_specificity(tn=0, fp=0)
        assert result is None


class TestAccuracy:
    """Test the calculate_accuracy function."""

    def test_accuracy_normal(self) -> None:
        """Test accuracy calculation with normal values."""
        result = calculate_accuracy(tp=8, tn=7, fp=3, fn=2)
        assert result == pytest.approx(0.75)  # (8+7)/(8+7+3+2) = 15/20

    def test_accuracy_perfect(self) -> None:
        """Test accuracy calculation with perfect score."""
        result = calculate_accuracy(tp=10, tn=10, fp=0, fn=0)
        assert result == pytest.approx(1.0)

    def test_accuracy_zero(self) -> None:
        """Test accuracy calculation with zero correct predictions."""
        result = calculate_accuracy(tp=0, tn=0, fp=5, fn=5)
        assert result == pytest.approx(0.0)

    def test_accuracy_undefined(self) -> None:
        """Test accuracy calculation when undefined (total = 0)."""
        result = calculate_accuracy(tp=0, tn=0, fp=0, fn=0)
        assert result is None


class TestPPV:
    """Test the calculate_ppv function."""

    def test_ppv_normal(self) -> None:
        """Test PPV calculation with normal values."""
        result = calculate_ppv(tp=8, fp=2)
        assert result == pytest.approx(0.8)

    def test_ppv_perfect(self) -> None:
        """Test PPV calculation with perfect score."""
        result = calculate_ppv(tp=10, fp=0)
        assert result == pytest.approx(1.0)

    def test_ppv_zero_tp(self) -> None:
        """Test PPV calculation with zero TP."""
        result = calculate_ppv(tp=0, fp=5)
        assert result == pytest.approx(0.0)

    def test_ppv_undefined(self) -> None:
        """Test PPV calculation when undefined (TP + FP = 0)."""
        result = calculate_ppv(tp=0, fp=0)
        assert result is None


class TestNPV:
    """Test the calculate_npv function."""

    def test_npv_normal(self) -> None:
        """Test NPV calculation with normal values."""
        result = calculate_npv(tn=7, fn=3)
        assert result == pytest.approx(0.7)

    def test_npv_perfect(self) -> None:
        """Test NPV calculation with perfect score."""
        result = calculate_npv(tn=10, fn=0)
        assert result == pytest.approx(1.0)

    def test_npv_zero_tn(self) -> None:
        """Test NPV calculation with zero TN."""
        result = calculate_npv(tn=0, fn=5)
        assert result == pytest.approx(0.0)

    def test_npv_undefined(self) -> None:
        """Test NPV calculation when undefined (TN + FN = 0)."""
        result = calculate_npv(tn=0, fn=0)
        assert result is None


class TestF1Score:
    """Test the calculate_f1_score function."""

    def test_f1_score_normal(self) -> None:
        """Test F1 score calculation with normal values."""
        # F1 = 2 * TP / (2 * TP + FP + FN) = 2 * 8 / (2 * 8 + 2 + 2) = 16 / 20 = 0.8
        result = calculate_f1_score(tp=8, fp=2, fn=2)
        assert result == pytest.approx(0.8)

    def test_f1_score_perfect(self) -> None:
        """Test F1 score calculation with perfect score."""
        result = calculate_f1_score(tp=10, fp=0, fn=0)
        assert result == pytest.approx(1.0)

    def test_f1_score_zero_tp(self) -> None:
        """Test F1 score calculation with zero TP."""
        result = calculate_f1_score(tp=0, fp=5, fn=5)
        assert result == pytest.approx(0.0)

    def test_f1_score_undefined(self) -> None:
        """Test F1 score calculation when undefined (2*TP + FP + FN = 0)."""
        result = calculate_f1_score(tp=0, fp=0, fn=0)
        assert result is None


class TestMCC:
    """Test the calculate_mcc function."""

    def test_mcc_normal(self) -> None:
        """Test MCC calculation with normal values."""
        # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
        # MCC = (8*7 - 2*3) / sqrt((8+2)(8+3)(7+2)(7+3)) = (56-6) / sqrt(10*11*9*10) = 50 / sqrt(9900)
        result = calculate_mcc(tp=8, tn=7, fp=2, fn=3)
        expected = 50 / math.sqrt(9900)
        assert result == pytest.approx(expected)

    def test_mcc_perfect_positive(self) -> None:
        """Test MCC calculation with perfect positive correlation."""
        result = calculate_mcc(tp=10, tn=10, fp=0, fn=0)
        assert result == pytest.approx(1.0)

    def test_mcc_perfect_negative(self) -> None:
        """Test MCC calculation with perfect negative correlation."""
        result = calculate_mcc(tp=0, tn=0, fp=10, fn=10)
        assert result == pytest.approx(-1.0)

    def test_mcc_undefined(self) -> None:
        """Test MCC calculation when undefined (denominator = 0)."""
        result = calculate_mcc(tp=0, tn=0, fp=0, fn=0)
        assert result is None

    def test_mcc_undefined_partial_zero(self) -> None:
        """Test MCC calculation when one denominator term is zero."""
        result = calculate_mcc(tp=0, tn=5, fp=0, fn=5)  # (TP+FP) = 0
        assert result is None


class TestCohenKappa:
    """Test the calculate_cohen_kappa function."""

    def test_cohen_kappa_normal(self) -> None:
        """Test Cohen's Kappa calculation with normal values."""
        tp, tn, fp, fn = 8, 7, 2, 3
        total = tp + tn + fp + fn  # 20
        p_o = (tp + tn) / total  # 15/20 = 0.75
        p_e = ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (total * total)
        # p_e = ((8+2)*(8+3) + (7+2)*(7+3)) / (20*20) = (10*11 + 9*10) / 400 = (110+90)/400 = 0.5
        expected = (p_o - p_e) / (
            1 - p_e
        )  # (0.75 - 0.5) / (1 - 0.5) = 0.25 / 0.5 = 0.5

        result = calculate_cohen_kappa(tp=8, tn=7, fp=2, fn=3)
        assert result == pytest.approx(expected)

    def test_cohen_kappa_perfect_agreement(self) -> None:
        """Test Cohen's Kappa calculation with perfect agreement."""
        result = calculate_cohen_kappa(tp=10, tn=10, fp=0, fn=0)
        assert result == pytest.approx(1.0)

    def test_cohen_kappa_undefined_total_zero(self) -> None:
        """Test Cohen's Kappa calculation when total = 0."""
        result = calculate_cohen_kappa(tp=0, tn=0, fp=0, fn=0)
        assert result is None

    def test_cohen_kappa_undefined_pe_one(self) -> None:
        """Test Cohen's Kappa calculation when expected agreement = 1.0."""
        # This is a theoretical edge case that's hard to construct naturally
        # We'll test the logic by mocking the scenario
        result = calculate_cohen_kappa(tp=10, tn=0, fp=0, fn=0)
        # In this case p_e = (10*10 + 0*0) / 100 = 1.0, so kappa is undefined
        assert result is None


class TestPABAK:
    """Test the calculate_pabak function."""

    def test_pabak_normal(self) -> None:
        """Test PABAK calculation with normal values."""
        tp, tn, fp, fn = 8, 7, 2, 3
        total = tp + tn + fp + fn  # 20
        p_o = (tp + tn) / total  # 15/20 = 0.75
        expected = 2 * p_o - 1  # 2 * 0.75 - 1 = 0.5

        result = calculate_pabak(tp=8, tn=7, fp=2, fn=3)
        assert result == pytest.approx(expected)

    def test_pabak_perfect_agreement(self) -> None:
        """Test PABAK calculation with perfect agreement."""
        result = calculate_pabak(tp=10, tn=10, fp=0, fn=0)
        assert result == pytest.approx(1.0)

    def test_pabak_no_agreement(self) -> None:
        """Test PABAK calculation with no agreement."""
        result = calculate_pabak(tp=0, tn=0, fp=10, fn=10)
        assert result == pytest.approx(-1.0)

    def test_pabak_undefined(self) -> None:
        """Test PABAK calculation when undefined (total = 0)."""
        result = calculate_pabak(tp=0, tn=0, fp=0, fn=0)
        assert result is None


class TestLikelihoodRatios:
    """Test the likelihood ratio functions."""

    def test_lr_plus_normal(self) -> None:
        """Test LR+ calculation with normal values."""
        sensitivity = 0.8
        specificity = 0.7
        expected = sensitivity / (1 - specificity)  # 0.8 / 0.3 = 2.667

        result = calculate_lr_plus(sensitivity, specificity)
        assert result == pytest.approx(expected)

    def test_lr_plus_perfect_specificity(self) -> None:
        """Test LR+ calculation when specificity = 1.0 (undefined)."""
        result = calculate_lr_plus(0.8, 1.0)
        assert result is None

    def test_lr_plus_none_inputs(self) -> None:
        """Test LR+ calculation with None inputs."""
        assert calculate_lr_plus(None, 0.7) is None
        assert calculate_lr_plus(0.8, None) is None
        assert calculate_lr_plus(None, None) is None

    def test_lr_minus_normal(self) -> None:
        """Test LR- calculation with normal values."""
        sensitivity = 0.8
        specificity = 0.7
        expected = (1 - sensitivity) / specificity  # 0.2 / 0.7 = 0.286

        result = calculate_lr_minus(sensitivity, specificity)
        assert result == pytest.approx(expected)

    def test_lr_minus_zero_specificity(self) -> None:
        """Test LR- calculation when specificity = 0.0 (undefined)."""
        result = calculate_lr_minus(0.8, 0.0)
        assert result is None

    def test_lr_minus_none_inputs(self) -> None:
        """Test LR- calculation with None inputs."""
        assert calculate_lr_minus(None, 0.7) is None
        assert calculate_lr_minus(0.8, None) is None
        assert calculate_lr_minus(None, None) is None


class TestCalculateAllMetrics:
    """Test the calculate_all_metrics function."""

    def test_calculate_all_metrics_normal(self) -> None:
        """Test calculate_all_metrics with normal benchmark result items."""
        # Create mock BenchmarkResultItem objects
        items = [
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=True
            ),  # TP
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=True
            ),  # TP
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=False
            ),  # FP
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=False
            ),  # TN
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=False
            ),  # TN
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=True
            ),  # FN
        ]

        result = calculate_all_metrics(items)

        # Verify it returns a BenchmarkRunUpdate schema
        assert isinstance(result, schemas.BenchmarkRunUpdate)

        # Check confusion matrix counts
        assert result.tp == 2
        assert result.fp == 1
        assert result.tn == 2
        assert result.fn == 1

        # Check calculated metrics (spot check a few)
        assert result.sensitivity == pytest.approx(2 / 3)  # TP/(TP+FN) = 2/3
        assert result.specificity == pytest.approx(2 / 3)  # TN/(TN+FP) = 2/3
        assert result.accuracy == pytest.approx(4 / 6)  # (TP+TN)/Total = 4/6

    def test_calculate_all_metrics_empty(self) -> None:
        """Test calculate_all_metrics with empty list."""
        result = calculate_all_metrics([])

        assert isinstance(result, schemas.BenchmarkRunUpdate)
        assert result.tp == 0
        assert result.fp == 0
        assert result.tn == 0
        assert result.fn == 0

        # All metrics should be None when no data
        assert result.sensitivity is None
        assert result.specificity is None
        assert result.accuracy is None


class TestCalculateAndUpdateBenchmarkMetrics:
    """Test the calculate_and_update_benchmark_metrics function."""

    def test_calculate_and_update_benchmark_metrics_success(
        self, mocker: MockerFixture
    ) -> None:
        """Test successful metrics calculation and database update."""
        # Mock dependencies
        mock_session = mocker.MagicMock()
        mock_benchmark_run_repo = mocker.MagicMock()
        mock_benchmark_result_item_repo = mocker.MagicMock()

        benchmark_run_id = uuid.uuid4()

        # Mock benchmark run
        mock_benchmark_run = mocker.MagicMock()
        mock_benchmark_run_repo.get_by_id.return_value = mock_benchmark_run

        # Mock result items
        mock_result_items = [
            MagicMock(
                final_decision=ScreeningDecisionType.INCLUDE, human_decision=True
            ),  # TP
            MagicMock(
                final_decision=ScreeningDecisionType.EXCLUDE, human_decision=False
            ),  # TN
        ]
        mock_benchmark_result_item_repo.get_by_benchmark_run_id.return_value = (
            mock_result_items
        )

        # Mock updated run
        mock_updated_run = mocker.MagicMock()
        mock_benchmark_run_repo.update.return_value = mock_updated_run

        # Call the function
        result = calculate_and_update_benchmark_metrics(
            session=mock_session,
            benchmark_run_id=benchmark_run_id,
            benchmark_run_repo=mock_benchmark_run_repo,
            benchmark_result_item_repo=mock_benchmark_result_item_repo,
        )

        # Verify calls
        mock_benchmark_run_repo.get_by_id.assert_called_once_with(
            mock_session, benchmark_run_id
        )
        mock_benchmark_result_item_repo.get_by_benchmark_run_id.assert_called_once_with(
            mock_session, benchmark_run_id
        )
        mock_benchmark_run_repo.update.assert_called_once_with(
            mock_session, mock_benchmark_run
        )

        # Verify result
        assert result == mock_updated_run

    def test_calculate_and_update_benchmark_metrics_run_not_found(
        self, mocker: MockerFixture
    ) -> None:
        """Test error when benchmark run is not found."""
        # Mock dependencies
        mock_session = mocker.MagicMock()
        mock_benchmark_run_repo = mocker.MagicMock()
        mock_benchmark_result_item_repo = mocker.MagicMock()

        benchmark_run_id = uuid.uuid4()

        # Mock benchmark run not found
        mock_benchmark_run_repo.get_by_id.return_value = None

        # Call the function and expect ValueError
        with pytest.raises(
            ValueError, match=f"BenchmarkRun with ID {benchmark_run_id} not found"
        ):
            calculate_and_update_benchmark_metrics(
                session=mock_session,
                benchmark_run_id=benchmark_run_id,
                benchmark_run_repo=mock_benchmark_run_repo,
                benchmark_result_item_repo=mock_benchmark_result_item_repo,
            )

    def test_calculate_and_update_benchmark_metrics_no_result_items(
        self, mocker: MockerFixture
    ) -> None:
        """Test metrics calculation with no result items."""
        # Mock dependencies
        mock_session = mocker.MagicMock()
        mock_benchmark_run_repo = mocker.MagicMock()
        mock_benchmark_result_item_repo = mocker.MagicMock()

        benchmark_run_id = uuid.uuid4()

        # Mock benchmark run
        mock_benchmark_run = mocker.MagicMock()
        mock_benchmark_run_repo.get_by_id.return_value = mock_benchmark_run

        # Mock no result items
        mock_benchmark_result_item_repo.get_by_benchmark_run_id.return_value = []

        # Mock updated run
        mock_updated_run = mocker.MagicMock()
        mock_benchmark_run_repo.update.return_value = mock_updated_run

        # Call the function
        result = calculate_and_update_benchmark_metrics(
            session=mock_session,
            benchmark_run_id=benchmark_run_id,
            benchmark_run_repo=mock_benchmark_run_repo,
            benchmark_result_item_repo=mock_benchmark_result_item_repo,
        )

        # Should still work with empty list
        assert result == mock_updated_run
        mock_benchmark_run_repo.update.assert_called_once()
