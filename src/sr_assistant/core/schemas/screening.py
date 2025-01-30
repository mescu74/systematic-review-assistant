from __future__ import annotations

import typing as t

from pydantic import BaseModel, Field, field_validator

from sr_assistant.core.types.screening import ExclusionReasonType


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
        description="The confidence score for the decision. [0.0, 1.0]. If your confidence score is below 0.9, assign 'uncertain' to decision. Confidence scrores are very important and are used to calibrate the system and address differences between models.",
    )
    rationale: str = Field(
        ...,
        description="The rationale for the decision. Be specific. Explainable AI is imporant. Don't just reiterate the category/categories, but explain how and why you reached this conclusion.",
    )
    exclusion_reason_categories: list[ExclusionReasonType] | None = Field(
        default=None,
        description="Omit if the decision is 'include'. If the decision is 'exclude' or 'uncertain', The PRISMA exclusion reason categories for the decision. This complements the 'rationale' field.",
    )
