from __future__ import annotations

import typing as t
from functools import cache

import streamlit as st
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
)
from pydantic.networks import PostgresDsn  # noqa: TC002
from pydantic_settings import BaseSettings

import sr_assistant.app.utils as ut
from sr_assistant.core.constants import ModelId


class ScreeningTemperature(BaseModel):
    """Temperature for first, second, and third screeners.

    If using different models, keep at 0.
    """

    first: float = Field(default=0.0, ge=0.0, le=1.0)
    second: float = Field(default=0.0, ge=0.0, le=1.0)
    third: float = Field(default=0.0, ge=0.0, le=1.0)


class Settings(BaseSettings):
    """Application settings.

    These are loaded from environment variables. Can be prefixed with SRA_ or sra_.

    # TODO: add types and defaults
    Attributes:
        OPENAI_API_KEY (str): OpenAI API key, from OPENAI_API_KEY env var or APIKEY
        ANTHROPIC_API_KEY (str): Anthropic API key
        SUPABASE_URL (str): Supabase URL
        SUPABASE_KEY (str): Supabase service role key
        DATABASE_URL (PostgresDsn): Database URL
        IN_STREAMLIT (bool, optional): Whether we're executed via streamlit run or not.
            Defaults to checking st context of main thread.
        NCBI_EMAIL (str): NCBI email
        NCBI_API_KEY (str): NCBI API key
        LOG_LEVEL (str): Log level
        env: (Literal["local", "test", "prototype"], optional): Environment to run in.
            ``test`` is for integration tests and uses ``sra_integration_test``
            SupaBase DB. Defaults to ``prototype``.
            Loaded from env vars ``env`` or ``environment`` also with ``sra_`` prefix,
            case insensitive.
        debug (bool): Enable debug mode
    """

    model_config: ConfigDict = {  # pyright: ignore
        "env_prefix": "sra_",  # pyright: ignore
        "env_file": (
            ".env",
            ".env.local",  # takes priority over .env
        ),
        "extra": "ignore",
    }

    # TODO: SecretStr or st.secrets
    OPENAI_API_KEY: str = Field(validation_alias="openai_api_key")

    ANTHROPIC_API_KEY: str = Field(
        default="",
        validation_alias="anthropic_api_key",
        description="Optional Anthropic API key, not used in proto ATM",
    )
    """Not used."""

    SUPABASE_URL: str = Field(validation_alias="sra_supabase_url")
    SUPABASE_KEY: str = Field(
        validation_alias="sra_supabase_service_role_key"
    )  # TODO: only for local proto, this should never be exposed client side
    SUPABASE_DIRECT_URL: PostgresDsn = Field(
        default="postgresql://postgres:postgres@127.0.0.1:54322/postgres",  # pyright: ignore
        validation_alias="sra_supabase_direct_url",
    )
    DATABASE_URL: PostgresDsn = Field(
        default="postgresql://postgres:postgres@127.0.0.1:54322/postgres",  # pyright: ignore
        description="For test and prototype envs this should be a `postgresql+psycopg:// useing pooled SupaBase DSN.",
        validation_alias="sra_database_url",
    )
    """For test and prototype envs this should be a `postgresql+psycopg:// useing pooled SupaBase DSN."""
    NCBI_EMAIL: str = Field(validation_alias="ncbi_email")
    NCBI_API_KEY: str = Field(validation_alias="ncbi_api_key")
    IN_STREAMLIT: bool = Field(
        description="Whether we're executed via streamlit run or not.",
        default_factory=ut.in_streamlit,
        validation_alias="in_streamlit",
    )

    # TODO: implement debug logging
    debug: bool = Field(
        default=False,
        description="Enable debug mode in app and asyncio, LangChain, etc. NOTE: logs variable values, do not use in production.",
    )
    """Currently doesn't do anything, just a placeholder for future use."""
    # TODO: use level enum
    log_level: str = Field(default="INFO")
    env: t.Literal["local", "test", "prototype"] = Field(
        default="prototype",
        validation_alias="environment",
        description="Environment to run in (local, test, prototype). test is for integration tests and uses test SupaBase DB.",
    )

    # TODO: remove these, not used
    llm_1: str = Field(default=ModelId.gpt_4o, min_length=2)
    llm_2: str = Field(default=ModelId.chatgpt, min_length=2)
    llm_3: str = Field(default=ModelId.sonnet_35, min_length=2)
    data_extraction_llm: str = Field(default=ModelId.gpt_4o, min_length=2)

    llm_1_temperature: float = Field(default=0.0, ge=0.0, le=1.0)
    llm_2_temperature: float = Field(default=0.0, ge=0.0, le=1.0)
    llm_3_temperature: float = Field(default=0.0, ge=0.0, le=1.0)


_cache_fn = st.cache_data if ut.in_streamlit() else cache


@_cache_fn
def get_settings() -> Settings:
    return Settings()  # pyright: ignore [reportCallIssue]
