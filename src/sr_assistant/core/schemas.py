from __future__ import annotations

import typing as t

from pydantic import BaseModel, Field

from pydantic.types import PositiveInt

from sr_assistant.core.types import ExclusionReasonType


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

# The JSONSchema of this model is passed as a tool-call schema for the model.
# OpenAI has specific limitations, like we can't set a minimum or maximum value for the
# confidence_score and have to rely on the description.
class ScreeningResponse(BaseModel):
    """Your systematic review screening response/decision for the given study.

    Consider the given context from the protocol carefully. Screening abstracts can be
    tricky, so assign confidence score and decision accordingly.
    """
    decision: t.Literal["include", "exclude", "uncertain"] = Field(
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
    exclusion_reason_categories: list[ExclusionReasonType] | None = Field(
        default=None,
        description="Omit if the decision is 'include'. If the decision is 'exclude' or 'uncertain', The PRISMA exclusion reason categories for the decision. This complements the 'rationale' field.",
    )
