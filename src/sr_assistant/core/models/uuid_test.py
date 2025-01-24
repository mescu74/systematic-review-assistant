from __future__ import annotations

from sqlmodel import Field, SQLModel

from sr_assistant.core.models.base import CreatedAtMixin, IdMixin


class UuidTest(IdMixin, CreatedAtMixin, SQLModel, table=True):
    """Test model for UUID handling."""

    __tablename__ = "uuid_test" # type: ignore pyright: ignore

    data: str = Field(...)
