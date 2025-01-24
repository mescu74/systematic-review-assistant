"""Utility functions for working with UUIDv7 timestamps and other common operations."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Literal

import uuid6


def id_to_timestamp(id_str: str, suffix: Literal["Z", "+00:00"] = "+00:00") -> str:
    """Convert a UUIDv7 string to an ISO 8601 timestamp.

    Args:
        id_str: UUIDv7 string to convert
        suffix: Timezone suffix format, either 'Z' or '+00:00' (default: '+00:00')

    Returns:
        ISO 8601 formatted timestamp string with the specified timezone suffix

    Examples:
        >>> id_to_timestamp("018f668c-1d7b-7cea-b078-2bf357352967")
        '2024-05-11T07:26:49.723000+00:00'
        >>> id_to_timestamp("018f668c-1d7b-7cea-b078-2bf357352967", suffix="Z")
        '2024-05-11T07:26:49.723000Z'
    """
    uuid_obj = uuid6.UUID(id_str)
    timestamp = datetime.fromtimestamp(uuid_obj.time / 1000, tz=timezone.utc)
    return timestamp.replace(tzinfo=None).isoformat() + suffix


def id_to_unix_ms(id_str: str) -> int:
    """Get Unix timestamp in milliseconds from a UUIDv7 string.

    Args:
        id_str: UUIDv7 string to convert

    Returns:
        Unix timestamp in milliseconds

    Examples:
        >>> id_to_unix_ms("018f668c-1d7b-7cea-b078-2bf357352967")
        1715412409723
    """
    return uuid6.UUID(id_str).time


def id_to_unix_s(id_str: str) -> float:
    """Get Unix timestamp in seconds from a UUIDv7 string.

    Args:
        id_str: UUIDv7 string to convert

    Returns:
        Unix timestamp in seconds with millisecond precision

    Examples:
        >>> id_to_unix_s("018f668c-1d7b-7cea-b078-2bf357352967")
        1715412410.723
    """
    return uuid6.UUID(id_str).time / 1000


def validate_uuid7(id_str: str) -> None:
    """Validate that a string is a valid UUIDv7.

    Args:
        id_str: String to validate

    Raises:
        ValueError: If the string is not a valid UUIDv7

    Examples:
        >>> validate_uuid7("018f668c-1d7b-7cea-b078-2bf357352967")  # Valid
        >>> validate_uuid7("123e4567-e89b-12d3-a456-426614174000")  # Invalid UUIDv4
        Traceback (most recent call last):
            ...
        ValueError: Invalid UUID version, want 7, got 4
    """
    if not id_str:
        msg = "UUID string is required"
        raise ValueError(msg)

    try:
        uuid_obj = uuid6.UUID(id_str)
    except ValueError as e:
        msg = f"Invalid UUID format: {e}"
        raise ValueError(msg) from e

    if uuid_obj.version != 7:
        msg = f"Invalid UUID version, want 7, got {uuid_obj.version}"
        raise ValueError(msg)
