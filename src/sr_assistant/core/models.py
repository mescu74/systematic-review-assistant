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

import enum
import uuid
from datetime import datetime  # noqa: TC003

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sa_pg
from pydantic import ConfigDict, JsonValue
from sqlmodel import Field, Relationship, SQLModel, String  # type: ignore

from sr_assistant.core.schemas import ExclusionReasons
from sr_assistant.core.types import ScreeningDecisionType, ScreeningStrategyType


def enum_values(enum_class: type[enum.StrEnum]) -> list[str]:
    """Get values for enum."""
    return [member.value for member in enum_class]


class SQLModelBase(SQLModel):
    model_config = ConfigDict(  # type: ignore
        from_attributes=True,
        validate_assignment=True,
    )


class SystematicReview(SQLModelBase, table=True):
    """Systematic review model.

    A systematic review project containing the research question, background,
    and inclusion/exclusion criteria. These are combined with search results to
    generate input to screening chain/agent.

    Notes:
        - In Python, pubmed_results field links to list of PubMedResult which again
          have conservative_result and comprehensive_result fields pointing to the
          abstracts screening results table.

    TODO:
        - PICO, etc., need to figure out UI updating by agent ...
    """

    __tablename__ = "systematic_reviews"  # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            onupdate=sa.func.utcnow(),
            nullable=True,
        ),
    )
    background: str | None = Field(
        default=None,
        description="Background context for the systematic review",
    )
    research_question: str = Field(
        description="The research question",
    )
    # TODO: sort this out, update ui
    #criteria_framework: str = Field(
    #    description="The inclusion/exclusion criteria framework",
    #)
    #criteria_framework_answers: dict | None = Field(
    #    default=None,
    #    description="The answers to the inclusion/exclusion criteria framework",
    #    sa_column=sa.Column(sa_pg.JSONB(), nullable=True),
    #)
    inclusion_criteria: str = Field(
        description="Inclusion criteria",
    )
    exclusion_criteria: str = Field(
        description="The exclusion criteria for studies",
    )
    # One review has many PubMedResults
    pubmed_results: list[PubMedResult] | None = Relationship(
        back_populates="review",
        sa_relationship_kwargs={"lazy": "selectin"},
    )


class PubMedResult(SQLModelBase, table=True):
    """PubMed search result.

    Stores both search info and result in one table for prototype simplicity.
    """

    __tablename__ = "pubmed_results"  # type: ignore pyright: ignore

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            onupdate=sa.func.utcnow(),
            nullable=True,
        ),
    )
    query: str = Field(
        description="PubMed search query",
        title="Query",
        schema_extra={"example": "(cancer) AND (immunotherapy) AND (clinical trial)"},
    )
    # Article info
    pmid: str = Field(
        index=True,
        title="PMID",
        schema_extra={"example": "39844819"},
    )
    pmc: str | None = Field(
        default=None,
        title="PMC ID",
        schema_extra={"example": "PMC11753851"},
    )
    doi: str | None = Field(
        default=None,
        title="DOI",
        schema_extra={"example": "10.1155/jimr/5845167"},
    )
    title: str = Field(title="Title")
    abstract: str = Field(title="Abstract")
    journal: str = Field(title="Journal")
    year: str = Field(title="Publishing Year")

    # One PubMedResult belongs to one SystematicReview
    review_id: uuid.UUID = Field(
        foreign_key="systematic_reviews.id", ondelete="CASCADE", index=True
    )
    review: SystematicReview = Relationship(
        back_populates="pubmed_results",
        cascade_delete=True,
        sa_relationship_kwargs={"lazy": "selectin"},
    )

    # One PubMedResult has two ScreenAbstractResults (conservative and comprehensive)
    conservative_result_id: uuid.UUID | None = Field(
        default=None, foreign_key="screen_abstract_results.id", index=True
    )
    comprehensive_result_id: uuid.UUID | None = Field(
        default=None, foreign_key="screen_abstract_results.id", index=True
    )
    conservative_result: ScreenAbstractResultModel = Relationship(
        sa_relationship_kwargs={
            "foreign_keys": "[PubMedResult.conservative_result_id]",
            "lazy": "selectin",
        }
    )
    comprehensive_result: ScreenAbstractResultModel = Relationship(
        sa_relationship_kwargs={
            "foreign_keys": "[PubMedResult.comprehensive_result_id]",
            "lazy": "selectin",
        }
    )

# ScreeningResult schema
# response fields:
# From model:
# - decision ScreeningDecisionType
# - confidence_score float
# - rationale str
# - extracted_quotes list[str]
# - exclusion_reason_categories ExclusionReasons
# Injected:
# - id uuid.UUID (run id)
# - review_id uuid.UUID
# - pubmed_result_id uuid.UUID
# - trace_id uuid.UUID
# - screening_strategy ScreeningStrategyType (comprehensive/conservative)
# - model_name str
# - start_time datetime (UTC)
# - end_time datetime (UTC)
# - response_metadata dict[str, JsonValue], JSONB in postgres
class ScreenAbstractResultModel(SQLModelBase, table=True):
    """Abstract screening decision model.

    Tracks the screening decisions for each paper in the review.
    """

    __tablename__ = "screen_abstract_results"  # type: ignore pyright: ignore

    id: uuid.UUID | None = Field(
        default_factory=uuid.uuid4, primary_key=True,
        description="Should be the run_id of the child-child `langsmith.RunTree`, populated by the chain `on_end` listener."
    )
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.func.utcnow(),
            onupdate=sa.func.utcnow(),
            nullable=True,
        ),
    )
    # Screening decision fields
    decision: ScreeningDecisionType = Field(
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
    rationale: str = Field(
        description="Rationale for the screening decision",
    )
    extracted_quotes: list[str] | None = Field(
        default=None,
        description="Supporting quotes from the title/abstract. Can be omitted if uncertain.",
        sa_column=sa.Column(sa_pg.ARRAY(String()), nullable=True),
    )
    exclusion_reason_categories: ExclusionReasons | None = Field(
        default=None,
        description="The PRISMA exclusion reason categories for the decision. This complements the 'rationale' field.",
        sa_column=sa.Column(sa_pg.JSONB(), nullable=True),
    )
    # Injected fields
    # Screening decision fields
    start_time: datetime | None = Field(
        description="Chain/graph invocation start time.",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            nullable=True,
        ),
    )
    end_time: datetime | None = Field(
        description="Chain/graph invocation end time.",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            nullable=True,
        ),
    )
    screening_strategy: ScreeningStrategyType = Field(
        description="Currently either 'compherensive' or 'conservative'. Maps to kind of prompt used.",
        sa_column=sa.Column(
            type_=sa_pg.ENUM(ScreeningStrategyType, values_callable=enum_values),
            nullable=False,
        ),
    )
    model_name: str = Field(
        description="The LangChain model name that generated this response.",
    )
    response_metadata: dict[str, JsonValue] | None = Field(
        default=None,
        description="Metadata from the chain/graph invocation.",
        sa_column=sa.Column(sa_pg.JSONB(), nullable=True),
    )

    # Screening Result is associated with one review
    review_id: uuid.UUID = Field(foreign_key="systematic_reviews.id", index=True)
    # Each ScreenAbstractResult belongs to one Review
    review: SystematicReview = Relationship(
        back_populates="screen_abstract_results",
        sa_relationship_kwargs={"lazy": "selectin"},
    )
    # Screening Result is associated with one PubMedResult
    pubmed_result_id: uuid.UUID = Field(
        description="The PubMed result being screened",
        foreign_key="pubmed_results.id",
        index=True,
    )
    # Each ScreenAbstractResult belongs to one PubMedResult. Relationship other way is
    # configured in PubMedResult.
    pubmed_result: PubMedResult = Relationship(
        sa_relationship_kwargs={"lazy": "selectin"},
    )
