# src/sr_assistant/core/schemas.py
# ruff: noqa: TC001 [many warnings are wrong, types are needed for SQLModel]

"""Core schemas for SR assistant.

Todo:
    - Migrate field descriptions to attr docstrings.
"""

from __future__ import annotations

import uuid
from collections.abc import MutableMapping  # noqa: TC003  # needed for tests
from typing import TypedDict

from pydantic import BaseModel, ConfigDict, Field, JsonValue, field_validator
from pydantic.types import AwareDatetime, PositiveInt  # noqa: TC002

from sr_assistant.core.types import (
    ComparisonExclusionReason,
    CriteriaFramework,
    InterventionExclusionReason,
    OutcomeExclusionReason,
    PopulationExclusionReason,
    ReportingExclusionReason,
    ScreeningDecisionType,
    ScreeningStrategyType,
    SearchDatabaseSource,
    StudyDesignExclusionReason,
)


class BaseSchema(BaseModel):  # noqa: D101
    model_config = ConfigDict(
        populate_by_name=True,  # use field names AND aliases
        arbitrary_types_allowed=True,
        validate_assignment=True,
        use_attribute_docstrings=True,  # populate '.description' from attr docstring
        revalidate_instances="always",
    )


class EventDetailBase(BaseSchema):
    """Base for all event details.

    All events are :class:`core.models.base.Event` and stored in the ``events`` table.
    It carries a ``StrEnum`` ``type`` field along with an event specific ``detail``
    schema.

    - Events must be persisted to a domain :class:`core.models.base.Event` table to
      populate ``id`` (UUIDv7) and ``created_at`` fields.
    - The database is the single source of truth for time and IDs. All timestamps are :
      class:`datetime.datetime` in UTC on Python schemas and models, and IDs are of type
      ``uuid6.UUID`..
    .- You can get a given event type's detail schema via the domain's
      ``<Domain>EventType`` ``StrEnum``'s ``.detail_schema`` member property.
      E.g.,:

      .. code-block:: python
          from core.schemas.protocol.events import ProtocolEventType

          schema = ProtocolEventType.InclusionCriteriaReviewed.detail_schema
          # -> ``core.schema.protocol.events.InclusionCriteriaReviewedDetail``
          evnet = InclusionCriteriaReviewed(detail=schema(**data))
          # repository injected in :class:`app.containers.Container`
          event = protocol_event_repository.save(event)
          # protocol_event_repository.get_by_type(event.type, limit=1)
          # protocol_event_repository.get(event.id)
          # app.utils.ts_from_id(event.id) == event.created_at
          # ... publish on channel or pass around, or let ``events`` table CDC actor
          # handle it.

    - This ``StrEnum``'s docstrings contain mermaid diagrams and domain knowledge,
       and document what the sources and commands are, and so forth. There's one
      EventType ``StrEnum`` for each domain consolidating this documentation as we
      have no event specific schemas apart from the details.

    See Also:
        - :mod:`core.schemas.protocol.events` for Protocol domain event schemas.
        - :mod:`core.schemas.search_strategy.events` for Search Strategy domain event
          schemas.
        - :mod:`core.schemas.screening.events` for Screening domain events detail
          schemas.
        - Data Extraction, Analysis, Reporting domains out of scope for prototype.

        And respective command modules.

    Attributes:
        version (int): Event version. Defaults to 1. Must be incremented in detail schema if mutated.
    """

    version: PositiveInt = Field(
        default=1,
        description="Event version. Must be updated in corresponding event detail schema if schema mutated.",
        title="Event Version",
    )


class ExclusionReasons(BaseSchema):
    """Exclusion reasons by category for PRISMA."""

    population_exclusion_reasons: list[PopulationExclusionReason] | None = Field(
        default=None
    )
    """Population exclusion reasons for PRISMA."""

    intervention_exclusion_reasons: list[InterventionExclusionReason] | None = Field(
        default=None
    )
    """Intervention exclusion reasons for PRISMA."""

    comparison_exclusion_reasons: list[ComparisonExclusionReason] | None = Field(
        default=None
    )
    """Comparison exclusion reasons for PRISMA."""

    outcome_exclusion_reasons: list[OutcomeExclusionReason] | None = Field(default=None)
    """Outcome exclusion reasons for PRISMA."""

    reporting_exclusion_reasons: list[ReportingExclusionReason] | None = Field(
        default=None
    )
    """Reporting exclusion reasons for PRISMA."""

    study_design_exclusion_reasons: list[StudyDesignExclusionReason] | None = Field(
        default=None
    )
    """Study design exclusion reasons for PRISMA."""


# The JSONSchema of this model is passed as a tool-call schema for the model.
# OpenAI has specific limitations, like we can't set a minimum or maximum value for the
# confidence_score and have to rely on the description.
class ScreeningResponse(BaseSchema):
    """Your systematic review screening response/decision for the given study.

    Consider the given context from the protocol carefully. Screening abstracts can be
    tricky, so assign confidence score and decision accordingly.
    """

    decision: ScreeningDecisionType = Field(...)
    """The systematic review abstracts screening decision. Should this study be included
    or not? Or are you uncertain? If your confidence score is below 0.8, assign 'uncertain'.
    """

    confidence_score: float = Field(...)
    """The systematic review abstracts screening confidence score [0.0, 1.0]. Should this study
    be included or not? Or are you uncertain? If your confidence score is below 0.8, assign 'uncertain'.
    """

    rationale: str = Field(...)
    """The rationale for the decision. Be specific. Explainable AI is important. Don't
    just reiterate the category/categories, but explain how and why you reached this conclusion.
    """

    extracted_quotes: list[str] | None = Field(default=None)
    """The supporting quotes from the title/abstract. Can be omitted if uncertain."""

    exclusion_reason_categories: ExclusionReasons | None = Field(default=None)
    """The PRISMA exclusion reason categories for the decision. Must be set if decision is 'exclude'.
    Omit if the decision is 'include'. If the decision is 'exclude' or 'uncertain', this complements
    the 'rationale' field.
    """


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
# - search_result_id uuid.UUID
# - screening_strategy ScreeningStrategyType
# - trace_id uuid.UUID
# - model_name str (passed to langchain wrapper)
# - start_time datetime (UTC)
# - end_time datetime (UTC)
# - response_metadata dict[str, JsonValue], JSONB in postgres
class ScreeningResult(ScreeningResponse):
    """ScreeningResponse processed and hydrated by the chain/graph listener/receiver.

    Currently screening uses LCEL chain with a listener callback for `on_end` event
    which receives the ScreeningResponse model and uses it to instantiate and return
    this model. `review_id` etc. are passed to the chain as metadata on invocation and
    the listener reads these and passes to this model.

    This model is further used to populate the abstracts screening SQLModel.

    We don't let the listener write directly to st.session_state because screening runs
    in threads and nothing in Streamlit is thread-safe. It's also simpler to let the
    invoker also handle the response directly.
    """

    id: uuid.UUID = Field(
        default_factory=uuid.uuid4,
    )
    """Screening result ID. `langsmith.RunTree`.id.

Populated by the chain/graph by instantiating this model. \
Should be the ``id`` of the invocation `langsmith.RunTree (e.g., run.chid.child.id when \
running in parallel with `langchain_core.runnables.RunnableParallel` and structured \
output)`.
"""

    review_id: uuid.UUID
    """Review ID. Read from RunnableConfig ``metadata['review_id']``.\
Can also be read from `st.session_state.review_id."""

    # TODO: should we remove this to break the circular relationship? Alembic warns
    # about it.
    search_result_id: uuid.UUID
    """SearchResult.id associated with this screening result.

Read from RunnableConfig metadata['search_result_id'] passed by the invoker."""

    trace_id: uuid.UUID
    """LangSmith `langsmith.RunTree`.trace_id, per input in batch, so for given \
screening input, the two reviewers share the `trace_id`. `id` is unique to the reviewer."""

    model_name: str
    """LangChain side model name."""

    screening_strategy: ScreeningStrategyType = Field(
        examples=["conservative", "comprehensive"],
    )
    """Screening strategy used to generate this response. Determines which prompt to use."""

    start_time: AwareDatetime
    """UTC start timestamp read from the chain/graph events."""

    end_time: AwareDatetime
    """UTC end timestamp read from the chain/graph events."""

    response_metadata: dict[str, JsonValue] = Field(
        default_factory=dict,
    )
    """Metadata read from the chain/graph events. inputs, invocation params, \
resp metadata, token usage, etc. JSONB in Postgres."""


class ScreeningResultCreate(BaseSchema):
    """Schema for creating a ScreenAbstractResult database record via the service layer.

    Used by ScreeningService.add_screening_decision to populate the
    ScreenAbstractResult SQLModel.
    """

    # Fields from ScreeningResponse
    decision: ScreeningDecisionType
    """The screening decision (INCLUDE, EXCLUDE, UNCERTAIN)."""

    confidence_score: float
    """The confidence score [0.0, 1.0] for the decision."""

    rationale: str
    """Specific rationale for the decision."""

    extracted_quotes: list[str] | None = None
    """Supporting quotes from the title/abstract."""

    exclusion_reason_categories: ExclusionReasons | None = None
    """PRISMA exclusion reason categories."""

    # Fields from ScreeningResult (contextual/metadata)
    id: uuid.UUID  # This is the LLM run_id, used as PK for ScreenAbstractResult model
    """Screening result ID, from the LangSmith Run ID of the reviewer invocation."""

    review_id: uuid.UUID
    """ID of the SystematicReview this screening pertains to."""

    trace_id: uuid.UUID | None = None
    """LangSmith trace_id, shared between reviewers for the same SearchResult input."""

    model_name: str
    """Name of the LLM used for this screening decision."""

    screening_strategy: ScreeningStrategyType
    """The strategy used (e.g., conservative, comprehensive)."""

    start_time: AwareDatetime | None = None
    """UTC timestamp when the LLM chain invocation started. Must be timezone-aware if not None."""

    end_time: AwareDatetime | None = None
    """UTC timestamp when the LLM chain invocation ended. Must be timezone-aware if not None."""

    response_metadata: dict[str, JsonValue] = Field(default_factory=dict)
    """Additional metadata from the LLM invocation (inputs, token usage, etc.)."""


class ScreeningResultUpdate(BaseSchema):
    """Schema for updating an existing ScreenAbstractResult database record.

    Used by ScreeningService.update_screening_decision for partial updates
    to a ScreenAbstractResult record. All fields are optional.
    """

    decision: ScreeningDecisionType | None = None
    """The screening decision."""

    confidence_score: float | None = None
    """The confidence score."""

    rationale: str | None = None
    """Specific rationale."""

    extracted_quotes: list[str] | None = None
    """Supporting quotes."""

    exclusion_reason_categories: ExclusionReasons | None = None
    """PRISMA exclusion categories."""

    trace_id: uuid.UUID | None = None
    """LangSmith trace_id."""

    model_name: str | None = None
    """Name of the LLM."""

    screening_strategy: ScreeningStrategyType | None = None
    """The strategy used."""

    start_time: AwareDatetime | None = None
    """Invocation start time. Must be timezone-aware if not None."""

    end_time: AwareDatetime | None = None
    """Invocation end time. Must be timezone-aware if not None."""

    response_metadata: dict[str, JsonValue] | None = None
    """Additional LLM metadata."""


class ScreeningResultRead(BaseSchema):
    """Schema for returning ScreenAbstractResult data from the service layer.

    Used for sending full ScreenAbstractResult data back to the caller.
    """

    id: uuid.UUID
    """Unique ID of this screening record (typically the LLM run_id)."""

    created_at: AwareDatetime | None
    """Timestamp of when this record was created in the database (UTC). Must be timezone-aware."""

    updated_at: AwareDatetime | None
    """Timestamp of when this record was last updated (UTC). Must be timezone-aware."""

    decision: ScreeningDecisionType
    """The screening decision."""

    confidence_score: float
    """The confidence score."""

    rationale: str
    """Rationale for the decision."""

    extracted_quotes: list[str] | None
    """Supporting quotes."""

    exclusion_reason_categories: dict[str, list[str]]  # Based on model's JSONB storage
    """PRISMA exclusion reason categories."""

    trace_id: uuid.UUID | None
    """LangSmith trace_id."""

    start_time: AwareDatetime | None
    """LLM chain invocation start time. Must be timezone-aware if not None."""

    end_time: AwareDatetime | None
    """LLM chain invocation end time. Must be timezone-aware if not None."""

    screening_strategy: ScreeningStrategyType
    """Screening strategy used."""

    model_name: str
    """LLM model name used."""

    response_metadata: dict[str, JsonValue]
    """Metadata from the LLM invocation."""

    review_id: uuid.UUID
    """ID of the parent SystematicReview."""


# --- Agent/Step Schemas --- # Added section


class PicosSuggestions(BaseSchema):
    """Structured PICO suggestions from LLM."""

    population: str
    """Suggested Population/Problem description based on context."""

    intervention: str
    """Suggested Intervention/Exposure description."""

    comparison: str
    """Suggested Comparison/Control description."""

    outcome: str
    """Suggested Outcome description."""

    general_critique: str
    """General critique and suggestions for the overall protocol."""


class SuggestionResult(TypedDict):
    """Return type for suggestion agent."""

    pico: PicosSuggestions | None
    raw_response: str


class ResolverOutputSchema(BaseSchema):
    """Schema for resolver model output that maps to ScreeningResolution database model.

    This is the structured output from the resolver model which will be used to
    create a ScreeningResolution model instance.
    """

    resolver_decision: ScreeningDecisionType = Field(...)
    """The decision made by the resolver for this search result. Should be one of: INCLUDE, EXCLUDE, or UNCERTAIN."""

    resolver_reasoning: str = Field(...)
    """Detailed reasoning for the resolver's decision. Explain your thought process clearly."""

    resolver_confidence_score: float = Field(...)
    """Confidence score for the resolver's decision, between 0.0 and 1.0. Higher values indicate more confidence."""

    @field_validator("resolver_confidence_score")
    @classmethod
    def check_confidence_score_range(cls, v: float) -> float:
        """Validate that resolver_confidence_score is between 0.0 and 1.0."""
        if not (0.0 <= v <= 1.0):
            msg = "resolver_confidence_score must be between 0.0 and 1.0"
            raise ValueError(msg)
        return v

    contributing_strategies: list[ScreeningStrategyType] = Field(default_factory=list)
    """List of original screening strategies (e.g., conservative, comprehensive) that the resolver found to have
    contributed meaningfully to its final decision, or aligns with. Can be empty.
    This field helps understand if the resolver sided more with one approach or found elements from multiple.
    """

    # These fields will be populated by the calling code
    review_id: uuid.UUID | None = Field(default=None)
    """ID of the systematic review. Will be populated by the caller."""

    search_result_id: uuid.UUID | None = Field(default=None)
    """ID of the search search result. Will be populated by the caller."""

    conservative_result_id: uuid.UUID | None = Field(default=None)
    """ID of the conservative screening result. Will be populated by the caller."""

    comprehensive_result_id: uuid.UUID | None = Field(default=None)
    """ID of the comprehensive screening result. Will be populated by the caller."""


# --- Systematic Review Schemas ---


class SystematicReviewBase(BaseSchema):
    """Base Pydantic schema for SystematicReview, encapsulating common fields."""

    research_question: str = Field(
        ..., description="The main research question of the systematic review."
    )
    background: str | None = Field(
        None, description="Background information and rationale for the review."
    )
    criteria_framework: CriteriaFramework | None = Field(
        default=CriteriaFramework.PICO,
        description="Framework used for criteria (e.g., PICO, PECO).",
    )
    criteria_framework_answers: MutableMapping[str, JsonValue] | None = Field(
        default=None,
        description="Answers related to the criteria framework (e.g., PICO elements as a dict).",
    )
    inclusion_criteria: str | None = Field(
        None, description="Detailed inclusion criteria for study selection."
    )
    exclusion_criteria: str | None = Field(
        None, description="Detailed exclusion criteria for study selection."
    )
    review_metadata: MutableMapping[str, JsonValue] | None = Field(
        default=None,
        description="Additional metadata for the review (e.g., tags, notes) as JSON.",
    )


class SystematicReviewCreate(SystematicReviewBase):
    """Schema for creating a new systematic review."""

    id: uuid.UUID | None = Field(
        default=None, description="Optional client-provided ID for the new review."
    )
    # Override to make it mandatory for creation, aligning with DB model
    exclusion_criteria: str = Field(
        ...,
        description="Detailed exclusion criteria for study selection. This is mandatory for new reviews.",
    )
    # research_question is inherited as mandatory from SystematicReviewBase


class SystematicReviewUpdate(SystematicReviewBase):
    """Schema for updating a Systematic Review (all fields optional)."""

    research_question: str | None = Field(
        default=None, description="The main research question of the systematic review."
    )
    background: str | None = Field(
        default=None, description="Background information and rationale for the review."
    )
    criteria_framework: CriteriaFramework | None = Field(
        default=None, description="Framework used for criteria (e.g., PICO, PECO)."
    )
    criteria_framework_answers: MutableMapping[str, JsonValue] | None = Field(
        default=None,
        description="Answers related to the criteria framework (e.g., PICO elements as a dict).",
    )
    inclusion_criteria: str | None = Field(
        default=None, description="Detailed inclusion criteria for study selection."
    )
    exclusion_criteria: str | None = Field(
        default=None, description="Detailed exclusion criteria for study selection."
    )
    review_metadata: MutableMapping[str, JsonValue] | None = Field(
        default=None,
        description="Additional metadata for the review (e.g., tags, notes) as JSON.",
    )


class SystematicReviewRead(SystematicReviewBase):
    """Schema for reading/returning a Systematic Review."""

    id: uuid.UUID

    created_at: AwareDatetime | None = None

    updated_at: AwareDatetime | None = None
    # Add relationships if needed for read operations, using nested schemas
    # e.g., search_results: list["SearchResultRead"] = []


# --- SearchResult Schemas ---
# (Define Read/Create/Update schemas if needed)


class SearchResultFilter(BaseSchema):
    """Schema for filtering SearchResult records in repository queries.

    All fields are optional to allow flexible querying.
    """

    review_id: uuid.UUID | None = Field(
        default=None, description="Filter by systematic review ID."
    )
    source_db: SearchDatabaseSource | None = Field(
        default=None, description="Filter by source database."
    )
    source_id: str | None = Field(
        default=None, description="Filter by source ID within the database."
    )
    doi: str | None = Field(
        default=None, description="Filter by Digital Object Identifier."
    )
    title: str | None = Field(
        default=None, description="Filter by title (case-insensitive partial match)."
    )  # Actual matching logic in repo
    year: str | None = Field(default=None, description="Filter by publication year.")
    # Consider adding more fields if needed, e.g., authors, keywords (though these might need more complex query logic)


class SearchResultRead(BaseSchema):
    """Schema for reading/returning SearchResult data from the service layer."""

    id: uuid.UUID
    """Unique identifier for the search result."""

    review_id: uuid.UUID
    """Identifier of the systematic review this search result belongs to."""

    source_db: SearchDatabaseSource
    """The source database from which this result was obtained (e.g., PubMed, Scopus)."""

    source_id: str
    """The unique identifier of this record within its source database (e.g., PMID)."""

    doi: str | None
    """Digital Object Identifier, if available."""

    title: str
    """The title of the publication."""

    abstract: str | None
    """The abstract of the publication."""

    journal: str | None
    """The name of the journal in which the publication appeared."""

    year: str | None
    """The publication year (as a string)."""

    authors: list[str] | None
    """A list of author names."""

    keywords: list[str] | None
    """A list of keywords associated with the publication."""

    raw_data: MutableMapping[str, JsonValue]
    """The original raw data record fetched from the source API."""

    source_metadata: MutableMapping[str, JsonValue]
    """Additional source-specific metadata."""

    created_at: AwareDatetime | None
    """Timestamp of when this search result was first stored (UTC). Must be timezone-aware."""

    updated_at: AwareDatetime | None
    """Timestamp of when this search result was last updated (UTC). Must be timezone-aware."""

    final_decision: ScreeningDecisionType | None = None
    """The final screening decision after any conflict resolution. Null if not yet resolved or no conflict."""

    resolution_id: uuid.UUID | None = None
    """Identifier of the ScreeningResolution record, if a conflict was resolved for this search result."""

    # Add relationships if needed
    # resolution: Optional["ScreeningResolutionRead"] = None
    # conservative_result: Optional["ScreenAbstractResultRead"] = None
    # comprehensive_result: Optional["ScreenAbstractResultRead"] = None


class SearchResultUpdate(BaseSchema):
    """Schema for partially updating an existing search result. All fields are optional."""

    doi: str | None = None
    """Digital Object Identifier."""

    title: str | None = None
    """The title of the publication."""

    abstract: str | None = None
    """The abstract of the publication."""

    journal: str | None = None
    """The name of the journal."""

    year: str | None = None
    """The publication year (as a string)."""

    authors: list[str] | None = None
    """A list of author names."""

    keywords: list[str] | None = None
    """A list of keywords."""

    raw_data: MutableMapping[str, JsonValue] | None = None
    """The original raw data record."""

    source_metadata: MutableMapping[str, JsonValue] | None = None
    """Additional source-specific metadata."""

    final_decision: ScreeningDecisionType | None = None
    """The final screening decision after any conflict resolution."""

    resolution_id: uuid.UUID | None = None
    """Identifier of the ScreeningResolution record."""


class ScreeningResolutionCreate(BaseSchema):
    """Schema for creating a ScreeningResolution database record.

    Used by ScreeningService.store_resolution_results to create a
    new ScreeningResolution record in the database.
    """

    search_result_id: uuid.UUID
    """ID of the SearchResult this resolution applies to."""

    review_id: uuid.UUID
    """ID of the SystematicReview."""

    # IDs of the original screening results that led to this resolution
    # FIXME: These two are NOT optional. There can be no resolution without at least two parties involved. Unless we use this for when human overrides a single LLM decision.
    conservative_result_id: uuid.UUID | None = None
    """ID of the conservative screening result, if applicable."""

    comprehensive_result_id: uuid.UUID | None = None
    """ID of the comprehensive screening result, if applicable."""

    # Data from the resolver LLM output
    resolver_decision: ScreeningDecisionType
    """The final decision made by the resolver."""

    resolver_reasoning: str
    """Detailed reasoning from the resolver."""

    resolver_confidence_score: float
    """Confidence score from the resolver."""

    resolver_model_name: str | None = None
    """Name of the LLM used for resolution."""

    response_metadata: dict[str, JsonValue] = Field(default_factory=dict)
    """Metadata from the resolver LLM invocation."""

    start_time: AwareDatetime | None = None
    """UTC timestamp of resolver LLM start. Must be timezone-aware if not None."""

    end_time: AwareDatetime | None = None
    """UTC timestamp of resolver LLM end. Must be timezone-aware if not None."""

    trace_id: uuid.UUID | None = None
    """LangSmith trace_id for the resolver invocation."""


class ScreeningResolutionRead(BaseSchema):
    """Schema for returning ScreeningResolution data from the service layer.

    Used for retrieving full ScreeningResolution data.
    """

    id: uuid.UUID
    """Unique ID of this screening resolution record."""

    created_at: AwareDatetime | None
    """Timestamp of when this record was created (UTC). Must be timezone-aware."""

    updated_at: AwareDatetime | None
    """Timestamp of when this record was last updated (UTC). Must be timezone-aware."""

    search_result_id: uuid.UUID
    """ID of the SearchResult this resolution applies to."""

    review_id: uuid.UUID
    """ID of the SystematicReview."""

    resolver_decision: ScreeningDecisionType
    """The final decision made by the resolver."""

    resolver_reasoning: str
    """Detailed reasoning from the resolver."""

    resolver_confidence_score: float
    """Confidence score from the resolver."""

    resolver_model_name: str | None
    """Name of the LLM used for resolution."""

    response_metadata: dict[str, JsonValue]
    """Metadata from the resolver LLM invocation."""

    start_time: AwareDatetime | None
    """UTC timestamp of resolver LLM start. Must be timezone-aware if not None."""

    end_time: AwareDatetime | None
    """UTC timestamp of resolver LLM end. Must be timezone-aware if not None."""

    trace_id: uuid.UUID | None
    """LangSmith trace_id for the resolver invocation."""
