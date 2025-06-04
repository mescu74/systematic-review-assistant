# Copyright 2025 Gareth Morgan
# SPDX-License-Identifier: MIT

"""Unit tests for benchmark run page."""

from __future__ import annotations

from sr_assistant.core.types import ScreeningDecisionType


class TestScreeningDecisionTypeValues:
    """Test that ScreeningDecisionType uses correct lowercase values."""

    def test_enum_values_are_lowercase(self) -> None:
        """Test that ScreeningDecisionType enum values are lowercase as expected."""
        assert ScreeningDecisionType.INCLUDE.value == "include"
        assert ScreeningDecisionType.EXCLUDE.value == "exclude"
        assert ScreeningDecisionType.UNCERTAIN.value == "uncertain"

    def test_enum_comparison_with_strings(self) -> None:
        """Test that enum values can be compared with lowercase strings."""
        decision = ScreeningDecisionType.INCLUDE
        assert decision.value == "include"
        assert decision.value != "INCLUDE"  # Uppercase should not match

        decision = ScreeningDecisionType.EXCLUDE
        assert decision.value == "exclude"
        assert decision.value != "EXCLUDE"  # Uppercase should not match

        decision = ScreeningDecisionType.UNCERTAIN
        assert decision.value == "uncertain"
        assert decision.value != "UNCERTAIN"  # Uppercase should not match
