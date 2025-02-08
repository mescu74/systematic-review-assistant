"""Utility functions for working with UUIDv7 timestamps and other common operations.

The UUID7 type is defined in `sr_assistant.core.types`.
"""

# ruff: noqa: F401 # unused import

from __future__ import annotations

from datetime import datetime, timezone
from zoneinfo import ZoneInfo
import typing as t
import uuid

import regex as re
from pydantic import GetPydanticSchema
from pydantic_core import core_schema
from loguru import logger
import uuid6
from langchain_core.runnables.config import (
    ContextThreadPoolExecutor,
    acall_func_with_variable_args,
    call_func_with_variable_args,
    ensure_config,
    get_async_callback_manager_for_config,
    get_callback_manager_for_config,
    get_executor_for_config,
    merge_configs,
    patch_config,
    run_in_executor,
)
import asyncio
from concurrent import futures as cf


"""ContextThreadPoolExecutor.

ThreadPoolExecutor that copies the context to the child thread.

Methods:
    submit(func: Callable[..., T], *args: Any, **kwargs: Any) -> Future[T]:
        Add a function to the executor.

        Args:
            func (Callable[..., T]): The function to submit.
            *args (Any): The positional arguments to the function.
            **kwargs (Any): The keyword arguments to the function.

        Returns:
            Future[T]: The future for the function.

    map(fn: Callable[..., T], *iterables: Iterable[Any], timeout: float | None = None, chunksize: int = 1) -> Iterator[T]:
        Map a function to multiple iterables.

        Args:
            fn (Callable[..., T]): The function to map.
            *iterables (Iterable[Any]): The iterables to map over.
            timeout (float | None, optional): The timeout for the map.
                Defaults to None.
            chunksize (int, optional): The chunksize for the map. Defaults to 1.

        Returns:
            Iterator[T]: The iterator for the mapped function.

See `LangChain docs <https://python.langchain.com/api_reference/_modules/langchain_core/runnables/config.html#ContextThreadPoolExecutor>`_

This is especially useful in Streamlit where the script thread is the only one supposed
to invoke st APIs.
"""

def is_pmid(v: str) -> bool:
    return bool(re.match(r"^\d{1,8}$", v))

def is_pmcid(v: str) -> bool:
    return bool(re.match(r"^PMC\d{1,8}$", v))

def is_utc_datetime(dt: datetime) -> bool:
    return dt.tzinfo in (timezone.utc, ZoneInfo("UTC"))

def utcnow() -> datetime:
    return datetime.now(timezone.utc)

def is_uuid(u: str | uuid.UUID | uuid6.UUID) -> bool:
    """Check if input is a valid UUID."""
    for c in (uuid.UUID, uuid6.UUID):
        try:
            c(u)
            return True
        except ValueError:
            continue
    return False


def validate_uuid(val: t.Any, version: t.Literal[1, 3, 4, 5, 6, 7, 8]) -> uuid6.UUID|uuid.UUID:
    """Validate a UUID against given version: 1, 3, 4, 5, 6, 7, 8.

    Args:
        val: UUID to validate
        version: UUID version to validate against

    Returns:
        UUID object

    Raises:
        TypeError: if val is not a uuid.UUID or uuid6.UUID
        ValueError: if val is not a UUID{version}
    """
    if isinstance(val, str):
        try:
            return uuid.UUID(val)
        except Exception as e:
            logger.opt(exception=True).warning(
                f"Failed to convert {val!r} to uuid.UUID: {e!r}"
            )
            try:
                return uuid6.UUID(val)
            except Exception as e:
                logger.opt(exception=True).warning(
                    f"Failed to convert {val!r} to uuid6.UUID: {e!r}"
                )
    if not isinstance(val, uuid6.UUID) and not isinstance(val, uuid.UUID):
        msg = f"Expected a uuid.UUID or uuid6.UUID, got {type(val)}"
        raise TypeError(msg)
    if val.version != version:
        msg = f"Expected a UUID{version}, got UUID{val.version}"
        raise ValueError(msg)
    return val


def id_to_timestamp(id_str: str, suffix: t.Literal["Z", "+00:00"] = "+00:00") -> str:
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
