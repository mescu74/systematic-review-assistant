from __future__ import annotations

from typing import Literal  # pyright: ignore [reportShadowedImports]

import streamlit as st
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
)
from pydantic.networks import PostgresDsn  # noqa: TC002
from pydantic_settings import BaseSettings

from core.constants import ModelId


class ScreeningTemperature(BaseModel):
    """Temperature for first, second, and third screeners.

    If using different models, keep at 0.
    """

    first: float = Field(default=0.0, ge=0.0, le=1.0)
    second: float = Field(default=0.0, ge=0.0, le=1.0)
    third: float = Field(default=0.0, ge=0.0, le=1.0)


class Settings(BaseSettings):
    """Application settings.

    These are loaded from environment variables.

    # TODO: add types and defaults
    Attributes:
        OPENAI_API_KEY (str): OpenAI API key, from OPENAI_API_KEY env var or APIKEY
        ANTHROPIC_API_KEY (str): Anthropic API key, from ANTHROPIC_API_KEY env var
        SUPABASE_URL: Supabase URL, from SRA_SUPABASE_URL env var
        SUPABASE_KEY: Supabase key (for local prototype we abuse role key)
        SUPABASE_DIRECT_URL: Supabase direct URL
        NCBI_EMAIL: NCBI email
        NCBI_API_KEY: NCBI API key
        debug: Debug mode
        env: Environment
        llm_1: First LLM model
        llm_2: Second LLM model
        llm_3: Third LLM model
        data_extraction_llm: LLM for data extraction (out of scope for prototype)
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
    ANTHROPIC_API_KEY: str = Field(validation_alias="anthropic_api_key")
    SUPABASE_URL: str = Field(validation_alias="sra_supabase_url")
    SUPABASE_KEY: str = Field(
        validation_alias="sra_supabase_service_role_key"
    )  # TODO: only for local proto, this should never be exposed client side
    SUPABASE_DIRECT_URL: PostgresDsn = Field(
        default="postgresql://postgres:postgres@127.0.0.1:54322/postgres",  # pyright: ignore
        validation_alias="sra_supabase_direct_url",
    )
    NCBI_EMAIL: str = Field(validation_alias="ncbi_email")
    NCBI_API_KEY: str = Field(validation_alias="ncbi_api_key")

    debug: bool = Field(default=False)
    env: Literal["local", "prototype"] = Field(
        default="local", validation_alias="environment"
    )

    llm_1: str = Field(default=ModelId.gpt_4o, min_length=2)
    llm_2: str = Field(default=ModelId.chatgpt, min_length=2)
    llm_3: str = Field(default=ModelId.sonnet_35, min_length=2)
    data_extraction_llm: str = Field(default=ModelId.gpt_4o, min_length=2)

    llm_1_temperature: float = Field(default=0.0, ge=0.0, le=1.0)
    llm_2_temperature: float = Field(default=0.0, ge=0.0, le=1.0)
    llm_3_temperature: float = Field(default=0.0, ge=0.0, le=1.0)


@st.cache_data
def get_settings() -> Settings:
    return Settings()  # pyright: ignore [reportCallIssue]
