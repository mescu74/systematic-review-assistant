# ruff: noqa: TC001 [many warnings are wrong, types are needed for SQLModel]

from __future__ import annotations

import typing as t
from datetime import datetime, timezone
import uuid

import streamlit as st
from langchain_core.runnables import RunnableConfig
from pydantic import BaseModel, ConfigDict, Field
from pydantic.types import AwareDatetime, JsonValue, PositiveInt  # ruff: noqa: TC002
#from typing_extensions import ReadOnly, Required

from sr_assistant.core import models
from sr_assistant.core.types import (
    ComparisonExclusionReason,
    InterventionExclusionReason,
    OutcomeExclusionReason,
    PopulationExclusionReason,
    ReportingExclusionReason,
    ScreeningDecisionType,
    ScreeningStrategyType,
    StudyDesignExclusionReason,
)


class EventDetailBase(BaseModel):
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

class ExclusionReasons(BaseModel):
    """Exclusion reasons by category for PRISMA."""

    model_config = ConfigDict(
        populate_by_name=True,
    )

    population_exclusion_reasons: list[PopulationExclusionReason] | None = Field(
        default=None,
        title="PopulationExclusionReasons",
        description="Population exclusion reasons for PRISMA.",
        examples=[
            ["Wrong age range for inclusion criteria", "Non-human subjects"],
            ["Condition or disease mismatch"],
        ],
    )
    intervention_exclusion_reasons: list[InterventionExclusionReason] | None = Field(
        default=None,
        title="InterventionExclusionReasons",
        description="Intervention exclusion reasons for PRISMA.",
        examples=[
            ["Intervention timing does not match criteria"],
            [
                "Dosage outside specified range",
                "Protocol deviation from inclusion criteria",
            ],
        ],
    )
    comparison_exclusion_reasons: list[ComparisonExclusionReason] | None = Field(
        default=None,
        title="ComparisonExclusionReasons",
        description="Comparison exclusion reasons for PRISMA.",
        examples=[
            [
                "Inappropriate control group",
                "Missing baseline data",
            ],
            [
                "Incorrect comparator",
            ],
        ],
    )
    outcome_exclusion_reasons: list[OutcomeExclusionReason] | None = Field(
        default=None,
        title="OutcomeExclusionReasons",
        description="Outcome exclusion reasons for PRISMA.",
        examples=[["Invalid outcome measurement method"]],
    )
    reporting_exclusion_reasons: list[ReportingExclusionReason] | None = Field(
        default=None,
        title="ReportingExclusionReasons",
        description="Reporting exclusion reasons for PRISMA.",
        examples=[
            [["Non-included language"], "Duplicate publication"],
            ["Pre-print publication"],
        ],
    )
    study_design_exclusion_reasons: list[StudyDesignExclusionReason] | None = Field(
        default=None,
        title="StudyDesignExclusionReasons",
        description="Study design exclusion reasons for PRISMA.",
        examples=[
            ["Wrong study design type", "Study not pre-registered"],
            ["Study duration too short"],
        ],
    )

# The JSONSchema of this model is passed as a tool-call schema for the model.
# OpenAI has specific limitations, like we can't set a minimum or maximum value for the
# confidence_score and have to rely on the description.
class ScreeningResponse(BaseModel):
    """Your systematic review screening response/decision for the given study.

    Consider the given context from the protocol carefully. Screening abstracts can be
    tricky, so assign confidence score and decision accordingly.
    """

    decision: ScreeningDecisionType = Field(
        ...,
        description="The systematic review abstracts screening decision. Should this study be included or not? Or are you uncertain? If your confidence score is below 0.8, assign 'uncertain'.",
    )
    confidence_score: float = Field(
        ...,
        description="The confidence score for the decision. [0.0, 1.0].",
    )
    rationale: str = Field(
        ...,
        description="The rationale for the decision. Be specific. Explainable AI is imporant. Don't just reiterate the category/categories, but explain how and why you reached this conclusion.",
    )
    extracted_quotes: list[str] | None = Field(
        default=None,
        description="Supporting quotes from the title/abstract. Can be omitted if uncertain.",
    )
    exclusion_reason_categories: ExclusionReasons | None = Field(
        default=None,
        description="The PRISMA exclusion reason categories for the decision. Must be set if decision is 'exclude'. Omit if the decision is 'include'. If the decision is 'exclude' or 'uncertain',  This complements the 'rationale' field.",
    )
# response fields:
# From model:
# - decision ScreeningDecisionType
# - confidence_score float
# - rationale str
# - extracted_quotes list[str]
# - exclusion_reason_categories ExclusionReasons
# Injected:
# - id uuid.UUID (trace id)
# - review_id uuid.UUID
# - pubmed_result_id uuid.UUID
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
        description="Screening result ID. RunTree.id. Populated by the chain/graph by instantiating this model. Should be the `id` of the invocation RunTree.",
    )
    review_id: uuid.UUUID = Field(
        description="Reviews ID. Read from RunnableConfig metadata['review_id'].Can also be read from `st.session_state.review_id",
    )
    pubmed_result_id: uuid.UUID = Field(
        description="PubMedResult.id associated with this screening result. Read from RunnableConfig metadata['pubmed_result_id'] passed by the invoker.",
    )
    trace_id: uuid.UUID = Field(
        description="LangSmith RunTree.trace_id, per input in batch, so for given screening input, the two reviewers share the trace_id. `id` is unique to the reviewer."
    )
    model_name: str = Field(
        description="LangChain side model name.",
    )
    screening_strategy: ScreeningStrategyType = Field(
        description="Screening strategy used to generate this response. Determines which prompt to use.",
        examples=["conservative", "comprehensive"],
    )
    start_time: AwareDatetime = Field(description="UTC start timestamp read from the chain/graph events.")
    end_time: AwareDatetime = Field(description="UTC end timestamp read from the chain/graph events.")
    response_metadata: dict[str, JsonValue] = Field(
        default_factory=dict,
        description="Metadata read from the chain/graph events. inputs, invocation params, resp metadata, token usage, etc. JSONB in Postgres.",
    )
