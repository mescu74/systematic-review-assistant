# src/sr_assistant/core/models/base.py
"""SQLModel/Alchemy models for SRA.

These are mostly for DB R/W and abstracted away from the
application layer which uses Pydantic schemas. SQLModel
has limitations and moving forward we'll switch to
SQLAlchemy and fully separate data and schema/validation/
serialization.
"""

import enum
import typing as t
import uuid
from collections.abc import Mapping, MutableMapping, Sequence
from datetime import datetime

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sa_pg
from pydantic import ConfigDict, JsonValue
from pydantic.types import PositiveInt
from sqlalchemy.ext.asyncio import AsyncAttrs
from sqlalchemy.orm import Mapped
from sqlmodel import Field, Relationship, SQLModel  # type: ignore

from sr_assistant.core.types import (
    CriteriaFramework,
    LogLevel,
    ScreeningDecisionType,
    ScreeningStrategyType,
    SearchDatabaseSource,
    UtcDatetime,
)


def add_gin_index(table_name: str, column_name: str) -> sa.Index:
    """Add a GIN index to a JSONB column."""
    return sa.Index(
        f"ix_{table_name}_on_{column_name}_gin",
        sa.text(column_name),
        postgresql_using="gin",
    )


def enum_values(enum_class: type[enum.Enum]) -> list[str]:
    """Get values for enum."""
    return [member.value for member in enum_class]


class SQLModelBase(AsyncAttrs, SQLModel):
    """Base SQLModel with ``awaitable_attrs`` mixin attribute.

    Attributes:
        awaitable_attrs: SQLAlchemy proxy mixin that makes all attributes awaitable.
    """

    model_config = ConfigDict(  # type: ignore
        from_attributes=True,
        validate_assignment=True,
        arbitrary_types_allowed=True,
        use_attribute_docstrings=True,
    )


Base = SQLModelBase


class SystematicReview(SQLModelBase, table=True):
    """Systematic review model.

    A systematic review project containing the research question, background,
    and inclusion/exclusion criteria. These are combined with search results to
    generate input to screening chain/agent.

    Notes:
        - In Python, search_results field links to list of SearchResult which again
          have conservative_result and comprehensive_result fields pointing to the
          abstracts screening results table.

    Todo:
        - PICO, etc., need to figure out UI updating by agent ...
    """

    _tablename: t.ClassVar[t.Literal["systematic_reviews"]] = "systematic_reviews"

    __tablename__ = _tablename  # pyright: ignore # type: ignore

    __table_args__ = (
        add_gin_index(_tablename, "criteria_framework_answers"),
        add_gin_index(_tablename, "review_metadata"),
    )

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            nullable=True,
        ),
    )
    """Database generated UTC timestamp when `SystematicReview` was created."""

    updated_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            onupdate=sa.func.now(),
            nullable=True,
        ),
    )
    """Database generated UTC timestamp when this `SystematicReview` was updated."""

    background: str = Field(
        default="",
        sa_column=sa.Column(sa.Text()),
    )
    """Background context for the systematic review."""

    research_question: str = Field(
        sa_column=sa.Column(sa.Text(), nullable=False),
    )
    """The research question."""

    # TODO: sort this out, update ui, these replace inclusion_criteria
    # This enum should have properties returning UI values and descriptions such that
    # a tabbed/selectbox form can be created dynamically, and preferably drafted by AI
    criteria_framework: CriteriaFramework | None = Field(
        default=None,
        description="The criteria framework.",
        sa_column=sa.Column(
            type_=sa_pg.ENUM(
                CriteriaFramework,
                name="criteriaframework_enum",
                values_callable=enum_values,
            ),
            nullable=True,
        ),
    )
    """The criteria framework."""

    criteria_framework_answers: MutableMapping[str, JsonValue] = Field(
        default_factory=dict,
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    """Answers to the criteria framework."""

    inclusion_criteria: str | None = Field(
        description="Inclusion criteria. To be replaced by above framework fields.",
        sa_column=sa.Column(sa.Text(), nullable=True),
    )
    """Inclusion criteria."""

    exclusion_criteria: str = Field(
        description="The exclusion criteria for studies",
        sa_column=sa.Column(sa.Text(), nullable=False),
    )
    """The exclusion criteria for studies."""

    # screen_abstract_results: Mapped[list["ScreenAbstractResult"]] = Relationship(
    #    back_populates="review",
    #    sa_relationship_kwargs={"lazy": "selectin"},
    # )
    """Screening results for the review."""

    log_records: Mapped[list["LogRecord"]] = Relationship(
        sa_relationship_kwargs={"lazy": "selectin"},
    )
    """Log records for the review."""

    # screening_resolutions: Mapped[list["ScreeningResolution"]] = Relationship(
    #    back_populates="review"
    # )
    """Screening resolutions associated with this review."""

    review_metadata: MutableMapping[str, JsonValue] = Field(
        default_factory=dict,
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    """Any metadata associated with the review."""

    # Relationship to the new SearchResult model


class SearchResult(SQLModelBase, table=True):
    """Unified search result model for various databases (PubMed, Scopus, etc.)."""

    _tablename: t.ClassVar[t.Literal["search_results"]] = "search_results"
    __tablename__ = _tablename  # pyright: ignore # type: ignore

    __table_args__ = (
        sa.UniqueConstraint(
            "review_id", "source_db", "source_id", name="uq_review_source_id"
        ),
        add_gin_index(_tablename, "keywords"),
        add_gin_index(_tablename, "raw_data"),
        add_gin_index(_tablename, "source_metadata"),
    )

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)

    created_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            onupdate=sa.func.now(),
            nullable=True,
        ),
    )
    review_id: uuid.UUID = Field(foreign_key="systematic_reviews.id", index=True)

    source_db: SearchDatabaseSource = Field(
        sa_column=sa.Column(
            sa_pg.ENUM(
                SearchDatabaseSource,
                name="searchdatabasesource_enum",
                values_callable=enum_values,
            ),
            index=True,
            nullable=False,
        ),
        description="Source database (e.g., 'PubMed', 'Scopus')",
    )
    source_id: str = Field(
        sa_column=sa.Column(sa.Text(), index=True, nullable=False),
        description="Unique identifier within the source database (e.g., PMID, Scopus ID)",
    )
    doi: str | None = Field(
        default=None,
        index=True,
        description="Digital Object Identifier",
    )
    title: str = Field(description="Article title")
    abstract: str | None = Field(default=None, sa_column=sa.Column(sa.Text()))
    journal: str | None = Field(default=None, description="Journal name")
    year: str | None = Field(default=None, index=True, description="Publication year")
    authors: list[str] | None = Field(
        default=None,
        description="List of authors",
        sa_column=sa.Column(sa_pg.ARRAY(sa.Text()), nullable=True),
    )  # Consider JSON list later
    keywords: list[str] | None = Field(
        default=None,
        sa_column=sa.Column(sa_pg.ARRAY(sa.Text()), index=True, nullable=True),
        description="List of keywords",
    )
    raw_data: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
        description="Original raw record from the source API",
    )
    source_metadata: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
        description="Source-specific metadata field",
    )

    # --- Foreign Keys for Screening Results ---
    resolution_id: uuid.UUID | None = Field(
        default=None, foreign_key="screening_resolutions.id", index=True, nullable=True
    )
    conservative_result_id: uuid.UUID | None = Field(
        default=None,
        foreign_key="screen_abstract_results.id",
        index=True,
        nullable=True,
    )
    comprehensive_result_id: uuid.UUID | None = Field(
        default=None,
        foreign_key="screen_abstract_results.id",
        index=True,
        nullable=True,
    )

    # --- Relationships ---
    # review: Mapped["SystematicReview"] = Relationship(
    #    back_populates="search_results",
    #    sa_relationship_kwargs={
    #        "lazy": "selectin",
    #        "foreign_keys": "[SearchResult.review_id]",
    #        "post_update": True,
    #    },
    # )
    # resolution: Mapped["ScreeningResolution"] | None = Relationship(
    #    sa_relationship_kwargs={
    #        "lazy": "selectin",
    #        "foreign_keys": "[SearchResult.resolution_id]",
    #        "uselist": False,
    #        "post_update": True,
    #    },
    # )
    # conservative_result: Mapped["ScreenAbstractResult"] | None = Relationship(
    #    sa_relationship_kwargs={
    #        "foreign_keys": "[SearchResult.conservative_result_id]",
    #        "lazy": "selectin",
    #        "post_update": True,
    #        "uselist": False,
    #    }
    # )
    # comprehensive_result: Mapped["ScreenAbstractResult"] | None = Relationship(
    #    sa_relationship_kwargs={
    #        "foreign_keys": "[SearchResult.comprehensive_result_id]",
    #        "lazy": "selectin",
    #        "post_update": True,
    #        "uselist": False,
    #    }
    # )


class ScreeningResolution(SQLModelBase, table=True):
    _tablename = "screening_resolutions"

    __tablename__ = _tablename  # pyright: ignore # type: ignore

    __table_args__ = (add_gin_index(_tablename, "response_metadata"),)

    # --- Fields ---
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            onupdate=sa.func.now(),
            nullable=True,
        ),
    )
    search_result_id: uuid.UUID = Field(foreign_key="search_results.id", index=True)
    review_id: uuid.UUID = Field(foreign_key="systematic_reviews.id", index=True)
    resolver_decision: ScreeningDecisionType = Field(
        sa_column=sa.Column(
            sa_pg.ENUM(
                ScreeningDecisionType,
                name="screeningdecisiontype",
                values_callable=enum_values,
                create_type=False,
            ),
            nullable=False,
            index=True,
        )
    )
    resolver_reasoning: str = Field(sa_column=sa.Column(sa.Text(), nullable=False))
    resolver_confidence_score: float = Field(ge=0.0, le=1.0, nullable=False)
    resolver_model_name: str | None = Field(default=None)
    response_metadata: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    start_time: datetime | None = Field(
        default=None, sa_column=sa.Column(sa.DateTime(timezone=True), nullable=True)
    )
    end_time: datetime | None = Field(
        default=None, sa_column=sa.Column(sa.DateTime(timezone=True), nullable=True)
    )
    trace_id: uuid.UUID | None = Field(default=None, index=True, nullable=True)

    # --- Relationships ---
    search_result: Mapped["SearchResult"] = Relationship(
        sa_relationship_kwargs={
            "lazy": "selectin",
            "foreign_keys": "[ScreeningResolution.search_result_id]",
            "post_update": True,
        },
    )
    review: Mapped["SystematicReview"] = Relationship(
        sa_relationship_kwargs={
            "lazy": "selectin",
            "foreign_keys": "[ScreeningResolution.review_id]",
            "post_update": True,
        },
    )


class ScreenAbstractResult(SQLModelBase, table=True):
    """Abstract screening decision model.

    Tracks the screening decisions for each paper in the review.

    Attributes:
        decision: ScreeningDecisionType
        confidence_score: float
        rationale: str
        extracted_quotes: list[str]
        exclusion_reason_categories: ExclusionReasons
        id: uuid.UUID (run id)
        review_id: uuid.UUID
        search_result_id: uuid.UUID
        trace_id: uuid.UUID
        model_name: str
        start_time: datetime (UTC)
        end_time: datetime (UTC)
        response_metadata: dict[str, JsonValue]
    """

    _tablename: t.ClassVar[t.Literal["screen_abstract_results"]] = (
        "screen_abstract_results"
    )

    __tablename__ = _tablename  # pyright: ignore # type: ignore

    __table_args__ = (add_gin_index(_tablename, "exclusion_reason_categories"),)

    id: uuid.UUID = Field(
        primary_key=True,
        description="Should be the run_id of the child-child `langsmith.RunTree`, populated by the chain `on_end` listener.",
    )
    created_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            nullable=True,
        ),
    )
    updated_at: datetime | None = Field(
        default=None,
        description="Database generated UTC timestamp",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            server_default=sa.text("TIMEZONE('UTC', CURRENT_TIMESTAMP)"),
            onupdate=sa.func.now(),
            nullable=True,
        ),
    )
    # Screening decision fields
    decision: ScreeningDecisionType = Field(
        description="Whether to include or exclude the paper, or mark as uncertain",
        sa_column=sa.Column(
            type_=sa_pg.ENUM(
                ScreeningDecisionType,
                name="screeningdecisiontype",
                values_callable=enum_values,
            ),
            nullable=False,
            index=True,
        ),
    )
    confidence_score: float = Field(
        ge=0.0,
        le=1.0,
        description="The confidence score for the decision. [0.0, 1.0].",
    )
    rationale: str = Field(
        description="Rationale for the screening decision",
        sa_column=sa.Column(sa.Text(), nullable=False),
    )
    extracted_quotes: Sequence[str] | None = Field(
        default=None,
        description="Supporting quotes from the title/abstract. Can be omitted if uncertain.",
        sa_column=sa.Column(sa_pg.ARRAY(sa.Text()), nullable=True),
    )
    exclusion_reason_categories: Mapping[str, list[str]] = Field(
        default_factory=dict,
        description="PRISMA exclusion reason categories.",
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    # Injected fields
    trace_id: uuid.UUID | None = Field(
        default=None,
        description="trace_id of the `RunTree` with which the listener was invoked. Populated by the chain `on_end` listener.",
        index=True,
    )
    # Screening decision fields
    start_time: datetime | None = Field(
        default=None,
        description="Chain/graph invocation start time.",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            nullable=True,
        ),
    )
    end_time: datetime | None = Field(
        default=None,
        description="Chain/graph invocation end time.",
        sa_column=sa.Column(
            sa.DateTime(timezone=True),
            nullable=True,
        ),
    )
    screening_strategy: ScreeningStrategyType = Field(
        description="Currently either 'compherensive' or 'conservative'. Maps to kind of prompt used.",
        sa_column=sa.Column(
            type_=sa_pg.ENUM(
                ScreeningStrategyType,
                name="screeningstrategytype",
                values_callable=enum_values,
            ),
            nullable=False,
            index=True,
        ),
    )
    model_name: str = Field(
        description="The LangChain model name that generated this response.",
    )
    response_metadata: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        description="Metadata from the LLM invocation.",
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )

    # Relationships
    review_id: uuid.UUID = Field(foreign_key="systematic_reviews.id", index=True)
    # review: Mapped["SystematicReview"] | None = Relationship(
    #    back_populates="screen_abstract_results",
    #    sa_relationship_kwargs={
    #        "lazy": "selectin",
    #        "foreign_keys": "[ScreenAbstractResult.review_id]",
    #        "post_update": True,
    #    },
    # )


class LogRecord(SQLModelBase, table=True):
    """Model for storing app log records."""

    _tablename: t.ClassVar[t.Literal["log_records"]] = "log_records"

    __tablename__ = _tablename  # pyright: ignore # type: ignore

    __table_args__ = (
        add_gin_index(_tablename, "extra"),
        add_gin_index(_tablename, "exception"),
    )

    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    timestamp: UtcDatetime = Field(
        sa_column=sa.Column(sa_pg.TIMESTAMP(timezone=True), nullable=False, index=True)
    )
    level: LogLevel = Field(
        sa_column=sa.Column(
            type_=sa_pg.ENUM(
                LogLevel, name="loglevel_enum", values_callable=enum_values
            ),
            nullable=False,
            index=True,
        ),
    )
    message: str = Field(sa_column=sa.Column(sa.Text(), nullable=False))
    module: str | None = Field(
        default=None, sa_column=sa.Column(sa.String(100), default=None, nullable=True)
    )
    name: str | None = Field(
        default=None, sa_column=sa.Column(sa.String(100), default=None, nullable=True)
    )
    function: str | None = Field(
        default=None, sa_column=sa.Column(sa.String(200), default=None, nullable=True)
    )
    line: PositiveInt | None = Field(
        default=None, sa_column=sa.Column(sa.Integer(), default=None, nullable=True)
    )
    thread: str | None = Field(default=None, description="{thread.name}:{thread.id}")
    process: str | None = Field(default=None, description="{process.name}:{process.id}")
    extra: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        description="User provided extra context.",
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    exception: Mapping[str, JsonValue] | None = Field(
        default=None,
        description="Exception information.",
        sa_column=sa.Column(sa_pg.JSONB, nullable=True),
    )
    record: Mapping[str, JsonValue] = Field(
        default_factory=dict,
        description="Full log record context.",
        sa_column=sa.Column(sa_pg.JSONB, nullable=False),
    )
    review_id: uuid.UUID | None = Field(
        default=None,
        description="Read from extra by the model_validator or st state, or by sink. None if not related to a review.",
        sa_column=sa.Column(
            sa.UUID,
            sa.ForeignKey("systematic_reviews.id"),
            default=None,
            index=True,
            nullable=True,
        ),
    )
