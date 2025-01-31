# src/sr_assistant/core/models/base.py
"""Common models and mixins.

In the prototype we just use one module for all models. Moving forward,
we may want to split these into separate modules.

Examples:
    >>> ret = (
    ...     supabase.table("reviews")
    ...     .insert(
    ...         json.loads(review.model_dump_json(exclude=["created_at", "updated_at"]))
    ...     )
    ...     .execute()
    ... )
    >>> review = Review.model_validate(ret.data[0])
"""

from __future__ import annotations

import uuid
from datetime import datetime  # noqa: TC003

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sa_pg
from pydantic import ConfigDict
from sqlmodel import Field, SQLModel, String  # type: ignore

from sr_assistant.core.types import ScreeningDecisionType


class SQLModelBase(SQLModel):
    model_config = ConfigDict( # type: ignore
        from_attributes=True,
        validate_assignment=True,
    )


class Review(SQLModelBase, table=True):
    """Systematic review model.

    A systematic review project containing the research question, background,
    and inclusion/exclusion criteria.
    """

    __tablename__ = "reviews"  # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        title="Created At",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
            name="created_at",
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            onupdate=sa.func.utcnow(),
            nullable=True,
        )
    )
    background: str = Field(
        default="",
        description="Background context for the systematic review",
    )
    question: str = Field(
        description="The research question being investigated",
    )
    inclusion_criteria: str = Field(
        description="The inclusion criteria for studies",
    )
    exclusion_criteria: str = Field(
        description="The exclusion criteria for studies",
    )

    # TODO: NOT WORKING!
    # pubmed_results: list[PubMedResult] = Relationship(back_populates="review")
    # pubmed_results: Mapped[list[PubMedResult]] = Relationship(back_populates="review")
    # pubmed_results: Mapped[List[PubMedResult]] = Relationship(back_populates="review")
    # pubmed_results: list[PubMedResult] = Relationship(back_populates="review")
    # pubmed_results: list["PubMedResult"] | None = Relationship(
    #    back_populates="review",
    #    sa_relationship_kwargs={"lazy": "selectin"}
    # )


class PubMedResult(SQLModelBase, table=True):
    """PubMed search result.

    Stores both search info and result in one table for prototype simplicity.
    """

    __tablename__ = "pubmed_results"  # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        title="Created At",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
            name="created_at",
        ),
    )
    review_id: uuid.UUID = Field(
        title="Review ID", foreign_key="reviews.id", index=True
    )
    query: str = Field(
        description="PubMed search query",
        title="Query",
        schema_extra={"example": "(cancer) AND (immunotherapy) AND (clinical trial)"},
    )

    # Article info
    pmid: str = Field(index=True, title="PMID", schema_extra={"example": "39844819"})
    pmc: str | None = Field(
        default=None, title="PMC ID", schema_extra={"example": "PMC11753851"}
    )
    doi: str | None = Field(
        default=None, title="DOI", schema_extra={"example": "10.1155/jimr/5845167"}
    )
    title: str = Field(title="Title")
    abstract: str = Field(title="Abstract")
    journal: str = Field(title="Journal")
    year: str = Field(title="Publishing Year")

    # TODO: NOT WORKING!
    # review: Review = Relationship(back_populates="pubmed_results")
    # review: Mapped[Review] = Relationship(back_populates="pubmed_results")
    # review: Mapped[Review] = Relationship(back_populates="pubmed_results", sa_relationship_kwargs={"lazy": "joined"})
    # review: Mapped["Review"] = Relationship(back_populates="pubmed_results")


def enum_values(enum_class: type[StrEnum]) -> list:
    """Get values for enum."""
    return [status.value for status in enum_class]

# a separate schema sr_assistant.core.schemas.screening.ScreeningResponse is used to
# coerce the model response to schema. These responses must be mapped to this model.
# TODO: how to?
class AbstractScreeningResult(SQLModelBase, table=True):
    """Abstract screening decision model.

    Tracks the screening decisions for each paper in the review.
    """

    __tablename__ = "abstract_screening_results" # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        title="Created At",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
            name="created_at",
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            onupdate=sa.func.utcnow(),
            nullable=True,
        )
    )
    review_id: uuid.UUID = Field(
        title="Review ID", foreign_key="reviews.id", index=True
    )
    search_result_id: uuid.UUID = Field(
        description="This is what we're screening",
        title="PubMed Result ID", foreign_key="pubmed_results.id", index=True
    )
    decision: ScreeningDecisionType = Field(
        title="Screening Decision",
        description="Whether to include or exclude the paper, or mark as uncertain",
        sa_column=sa.Column(
            type_=sa_pg.ENUM(ScreeningDecisionType, values_callable=enum_values),
            nullable=False,
        ),
    )
    confidence_score: float = Field(
        ge=0.0,
        le=1.0,
        description="The confidence score for the decision. [0.0, 1.0].",
    )
    rationale: str | None = Field(
        default=None,
        title="Decision Rationale",
        description="Rationale for the screening decision",
    )
    extracted_quotes: list[str] | None = Field(
        default=None,
        description="Supporting quotes from the title/abstract. Can be omitted if uncertain.",
        sa_column=sa.Column(sa_pg.ARRAY(String())),
    )
    exclusion_reason_categories: list[str] | None = Field(
        default=None,
        description="Omit if the decision is 'include'. If the decision is 'exclude' or 'uncertain', The PRISMA exclusion reason categories for the decision. This complements the 'rationale' field.",
        sa_column=sa.Column(sa_pg.ARRAY(String())),
    )


class UuidTest(SQLModelBase, table=True):
    """Test model for UUID handling."""

    __tablename__ = "uuid_test"  # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        title="Created At",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
        ),
    )
    data: str = Field(...)

