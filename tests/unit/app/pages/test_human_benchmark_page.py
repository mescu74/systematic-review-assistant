# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Unit tests for human benchmark page."""
# pyright: reportPrivateUsage=false

from __future__ import annotations

import uuid
from datetime import datetime, timezone

from sr_assistant.app.pages.human_benchmark_page import (
    _calculate_classification,
    _determine_final_decision,
    _needs_resolver,
)
from sr_assistant.core import schemas
from sr_assistant.core.types import ScreeningDecisionType, ScreeningStrategyType


class TestBenchmarkLogic:
    """Test the benchmark logic functions."""

    def test_needs_resolver_disagreement(self) -> None:
        """Test that resolver is needed when conservative and comprehensive disagree."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        assert _needs_resolver(conservative, comprehensive) is True

    def test_needs_resolver_both_uncertain(self) -> None:
        """Test that resolver is needed when both are uncertain."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.UNCERTAIN,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.UNCERTAIN,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        assert _needs_resolver(conservative, comprehensive) is True

    def test_needs_resolver_low_confidence(self) -> None:
        """Test that resolver is needed when confidence is low."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.7,  # Low confidence
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        assert _needs_resolver(conservative, comprehensive) is True

    def test_needs_resolver_agreement_high_confidence(self) -> None:
        """Test that resolver is not needed when both agree with high confidence."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        assert _needs_resolver(conservative, comprehensive) is False

    def test_determine_final_decision_with_resolver(self) -> None:
        """Test final decision determination when resolver is used."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.EXCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        resolver = schemas.ResolverOutputSchema(
            resolver_decision=ScreeningDecisionType.EXCLUDE,
            resolver_reasoning="Resolver reasoning",
            resolver_confidence_score=0.95,
        )

        final_decision = _determine_final_decision(
            conservative, comprehensive, resolver
        )
        assert final_decision == ScreeningDecisionType.EXCLUDE

    def test_determine_final_decision_without_resolver(self) -> None:
        """Test final decision determination when no resolver is used."""
        conservative = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.CONSERVATIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        comprehensive = schemas.ScreeningResult(
            id=uuid.uuid4(),
            review_id=uuid.uuid4(),
            search_result_id=uuid.uuid4(),
            trace_id=uuid.uuid4(),
            model_name="test-model",
            screening_strategy=ScreeningStrategyType.COMPREHENSIVE,
            start_time=datetime.now(timezone.utc),
            end_time=datetime.now(timezone.utc),
            decision=ScreeningDecisionType.INCLUDE,
            confidence_score=0.9,
            rationale="Test rationale",
        )

        final_decision = _determine_final_decision(conservative, comprehensive, None)
        assert final_decision == ScreeningDecisionType.INCLUDE

    def test_calculate_classification_tp(self) -> None:
        """Test classification calculation for True Positive."""
        classification = _calculate_classification(ScreeningDecisionType.INCLUDE, True)
        assert classification == "TP"

    def test_calculate_classification_fp(self) -> None:
        """Test classification calculation for False Positive."""
        classification = _calculate_classification(ScreeningDecisionType.INCLUDE, False)
        assert classification == "FP"

    def test_calculate_classification_tn(self) -> None:
        """Test classification calculation for True Negative."""
        classification = _calculate_classification(ScreeningDecisionType.EXCLUDE, False)
        assert classification == "TN"

    def test_calculate_classification_fn(self) -> None:
        """Test classification calculation for False Negative."""
        classification = _calculate_classification(ScreeningDecisionType.EXCLUDE, True)
        assert classification == "FN"

    def test_calculate_classification_unknown_human_decision(self) -> None:
        """Test classification calculation when human decision is unknown."""
        classification = _calculate_classification(ScreeningDecisionType.INCLUDE, None)
        assert classification == "UNKNOWN"

    def test_calculate_classification_uncertain_ai_decision(self) -> None:
        """Test classification calculation when AI decision is uncertain."""
        classification = _calculate_classification(
            ScreeningDecisionType.UNCERTAIN, True
        )
        assert classification == "FN"  # Uncertain is treated as exclude

        classification = _calculate_classification(
            ScreeningDecisionType.UNCERTAIN, False
        )
        assert classification == "TN"  # Uncertain is treated as exclude
