"""Global pytest configuration."""
from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _setup_test_env() -> None:
    """Set up test environment for all tests."""
    pass  # Add any global test setup here if needed
