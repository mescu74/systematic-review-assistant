"""Utility functions for working with UUIDv7 timestamps and other common operations.

The UUID7 type is defined in `sr_assistant.core.types`.
"""

# ruff: noqa: F401 # unused import

from __future__ import annotations

import asyncio
import typing as t
import uuid
from concurrent.futures.thread import ThreadPoolExecutor
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

import regex as re
import streamlit as st
import uuid6
from langchain_core.runnables.config import (
    ContextThreadPoolExecutor,
    acall_func_with_variable_args,
    call_func_with_variable_args,
    ensure_config,
    get_async_callback_manager_for_config,
    get_callback_manager_for_config,
    get_executor_for_config,  # TODO: reimplement with stconfig executor
    merge_configs,
    patch_config,
    run_in_executor,
)
from loguru import logger
from pydantic import GetPydanticSchema
from pydantic_core import core_schema
from streamlit.runtime.scriptrunner_utils.script_run_context import (
    ScriptRunContext,
    add_script_run_ctx,
    get_script_run_ctx,
)


class StContextThreadPoolExecutor(ContextThreadPoolExecutor):
    """ThreadPoolExecutor that copies both contextvars and Streamlit context.

    Wrapper around LangChain's ContextThreadPoolExecutor to add Streamlit context.

    See Also:
        - `LangChain ContextThreadPoolExecutor <https://python.langchain.com/api_reference/_modules/langchain_core/runnables/config.html#ContextThreadPoolExecutor>`_
        - `ThreadPoolExecutor https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor`_
        - `Executor <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.Executor>`_
    """

    def __init__(
        self,
        max_workers: int | None = None,
        thread_name_prefix: str = "",
        *,
        st_ctx: ScriptRunContext | None = None,
    ) -> None:
        """Initialize ThreadPoolExecutor with Streamlit context and contextvars.

        Args:
            max_workers (int | None, optional): Maximum number of worker threads. Defaults to None.
                If None, for >=3.8,<3.13 defaults to ``min(32, os.cpu_count() + 4)``,
                for >=3.13 defaults to ``min(32, (os.process_cpu_count() or 1) + 4)``.
            thread_name_prefix (str, optional): Prefix for thread names. Defaults to "".
            st_ctx (st.ScriptRunContext | None, optional): Streamlit context.
                Defaults to None. If None, uses the current thread's st context
                (one running the executor) by invoking st helper get_script_run_ctx().
        """
        _st_ctx = st_ctx or get_script_run_ctx()
        if _st_ctx is None:
            logger.warning("No Streamlit context found")
        super().__init__(
            max_workers=max_workers,
            thread_name_prefix=thread_name_prefix,
            initializer=lambda: add_script_run_ctx(None, _st_ctx),
        )


class StThreadPoolExecutor(ThreadPoolExecutor):
    """ThreadPoolExecutor that copies Streamlit context to threads.

    To copy context vars as well, use `StContextThreadPoolExecutor`.

    See Also:
        - `ThreadPoolExecutor https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor`_
        - `Executor <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.Executor>`_
    """

    def __init__(
        self,
        max_workers: int | None = None,
        thread_name_prefix: str = "",
        *,
        st_ctx: ScriptRunContext | None = None,
    ) -> None:
        """Initialize ThreadPoolExecutor with Streamlit context.

        Args:
            max_workers (int | None, optional): Maximum number of worker threads.
                Defaults to None.
                If None, for >=3.8,<3.13 defaults to ``min(32, os.cpu_count() + 4)``,
                for >=3.13 defaults to ``min(32, (os.process_cpu_count() or 1) + 4)``.
            thread_name_prefix (str, optional): Prefix for thread names. Defaults to "".
            st_ctx (st.ScriptRunContext | None, optional): Streamlit context.
                Defaults to None. If None, uses the current thread's st context
                (one running the executor) by invoking st helper get_script_run_ctx().
        """
        _st_ctx = st_ctx or get_script_run_ctx()
        if _st_ctx is None:
            logger.warning("No Streamlit context found")
        super().__init__(
            max_workers=max_workers,
            thread_name_prefix=thread_name_prefix,
            initializer=lambda: add_script_run_ctx(None, _st_ctx),
        )


def is_pmid(v: str) -> bool:
    return bool(re.match(r"^\d{1,8}$", v))


def is_pmcid(v: str) -> bool:
    return bool(re.match(r"^PMC\d{1,8}$", v))


def is_utc_datetime(dt: datetime) -> bool:
    return dt.tzinfo in (timezone.utc, ZoneInfo("UTC"))


def utcnow() -> datetime:
    return datetime.now(tz=timezone.utc)


def is_uuid(u: t.Any) -> uuid.UUID | uuid6.UUID | None:
    """Check if input `u` is a valid UUID.

    Supports v1-v8 UUIDs.

    Args:
        u (Any): The input to check.

    Returns:
        uuid.UUID | uuid6.UUID | None: UUID object if `u` is a valid UUID, None otherwise.
    """
    if isinstance(u, (uuid.UUID | uuid6.UUID)):
        return u
    for c in (uuid.UUID, uuid6.UUID):
        if isinstance(u, str):
            try:
                return c(hex=u)
            except (ValueError, TypeError):
                continue
    return None


def in_streamlit() -> bool:
    return bool(get_script_run_ctx(suppress_warning=True))


def init_state_key(key: str, value: t.Any) -> t.Any:
    """Initialize a session state key with the given value if it doesn't exist.

    Args:
        key: The key to initialize.
        value: The value to initialize the key with.

    Returns:
        The value of st.session_state[key] regardless if it was initialized or not.
    """
    if key not in st.session_state:
        st.session_state[key] = value
    return st.session_state[key]


def validate_uuid(
    val: t.Any, version: t.Literal[1, 3, 4, 5, 6, 7, 8]
) -> uuid6.UUID | uuid.UUID:
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
