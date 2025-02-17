from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelId:
    """Model API names."""

    gpt_4 = "gpt-4-turbo"
    gpt_4o = "gpt-4o"
    gpt_4o_mini = "gpt-4o-mini"
    chatgpt = "chatgpt-4o-latest"
    o1_preview = "o1-preview"
    o1_mini = "o1-mini"
    sonnet_35 = "claude-3-5-sonnet-latest"
    haiku_35 = "claude-3-5-haiku-latest"
    openai_embedding = "text-embedding-3-large"


EMBEDDING_DIMENSION = 1024
