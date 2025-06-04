"""Protocol definition page for the SRA.

- Allow user to define research question, PICO, inclusion/exclusion criteria.
- Provide LLM suggestions for PICO criteria.
- Persist protocol to database.
"""

from __future__ import annotations

import uuid
from datetime import UTC, datetime
from typing import TYPE_CHECKING, cast

import streamlit as st
from dotenv import find_dotenv, load_dotenv
from langchain_core.runnables import RunnableConfig
from loguru import logger
from streamlit.delta_generator import DeltaGenerator

from sr_assistant.app.services import ReviewService
from sr_assistant.core import schemas as app_schemas
from sr_assistant.core.models import CriteriaFramework, SystematicReview
from sr_assistant.core.schemas import PicosSuggestions, SuggestionResult
from sr_assistant.step1.suggestion_agent import SuggestionAgent

if TYPE_CHECKING:
    from streamlit.delta_generator import DeltaGenerator

    from sr_assistant.core.schemas import PicosSuggestions, SuggestionResult

load_dotenv(find_dotenv(), override=True)


def init_agent_config() -> RunnableConfig:
    """Initializes or retrieves the Langchain RunnableConfig for the agent.

    Ensures a consistent thread_id for configurable LLM runs within the session.
    """
    if "agent_config" in st.session_state:
        return cast("RunnableConfig", st.session_state.agent_config)

    config = RunnableConfig(
        recursion_limit=1000, configurable={"thread_id": str(uuid.uuid4())}
    )
    st.session_state.agent_config = config
    return config


def init_suggestion_agent() -> SuggestionAgent:
    """Initialize the suggestion agent."""
    return SuggestionAgent(model="gpt-4o", temperature=0.0)


@logger.catch(Exception, message="Error building review model from PICO")
def build_review_model_from_pico(response_widget: DeltaGenerator) -> SystematicReview:
    """Builds the SystematicReview model using PICO fields."""
    inclusion_criteria_parts = []

    if st.session_state.pico_population:
        inclusion_criteria_parts.append(
            f"Population: {st.session_state.pico_population}"
        )
    if st.session_state.pico_intervention:
        inclusion_criteria_parts.append(
            f"Intervention: {st.session_state.pico_intervention}"
        )
    if st.session_state.pico_comparison:
        inclusion_criteria_parts.append(
            f"Comparison: {st.session_state.pico_comparison}"
        )
    if st.session_state.pico_outcome:
        inclusion_criteria_parts.append(f"Outcome: {st.session_state.pico_outcome}")

    inclusion_criteria = (
        "; ".join(inclusion_criteria_parts)
        if inclusion_criteria_parts
        else "Not specified"
    )
    exclusion_criteria_value = st.session_state.get("exclusion_criteria")
    exclusion_criteria = (
        exclusion_criteria_value if exclusion_criteria_value else "Not specified"
    )

    logger.info(
        "Building review model with PICO derived criteria: Incl={}, Excl={}",
        inclusion_criteria,
        exclusion_criteria,
    )

    review_model_instance = SystematicReview(
        id=st.session_state.review_id,
        background=st.session_state.background,
        research_question=st.session_state.research_question,
        inclusion_criteria=inclusion_criteria,
        exclusion_criteria=exclusion_criteria,
        criteria_framework=CriteriaFramework.PICO,
        criteria_framework_answers={
            "population": st.session_state.pico_population,
            "intervention": st.session_state.pico_intervention,
            "comparison": st.session_state.pico_comparison,
            "outcome": st.session_state.pico_outcome,
        },
    )
    response_widget.success("SystematicReview model built from PICO successfully")
    return review_model_instance


def persist_review(review_model: SystematicReview) -> SystematicReview:
    """Persists the SystematicReview model using ReviewService."""
    review_service = ReviewService()

    try:
        existing_review = review_service.get_review(review_model.id)

        if existing_review:
            logger.debug(
                "Updating existing review {} via ReviewService", review_model.id
            )
            update_data = app_schemas.SystematicReviewUpdate(
                background=review_model.background,
                research_question=review_model.research_question,
                criteria_framework=review_model.criteria_framework,
                criteria_framework_answers=review_model.criteria_framework_answers,
                inclusion_criteria=review_model.inclusion_criteria,
                exclusion_criteria=review_model.exclusion_criteria,
                review_metadata=review_model.review_metadata,
            )
            persisted_review = review_service.update_review(
                review_id=review_model.id, review_update_data=update_data
            )
        else:
            logger.debug("Adding new review {} via ReviewService", review_model.id)
            create_data = app_schemas.SystematicReviewCreate(
                id=review_model.id,
                research_question=review_model.research_question
                if review_model.research_question
                else "Default RQ",
                exclusion_criteria=review_model.exclusion_criteria
                if review_model.exclusion_criteria
                else "Default Exclusion",
                background=review_model.background,
                criteria_framework=review_model.criteria_framework,
                criteria_framework_answers=review_model.criteria_framework_answers,
                inclusion_criteria=review_model.inclusion_criteria,
                review_metadata=review_model.review_metadata,
            )
            persisted_review = review_service.create_review(create_data)

        logger.info("Persisted review {} via ReviewService", persisted_review.id)
        return persisted_review
    except Exception:
        logger.exception(
            "Error persisting review {} using ReviewService", review_model.id
        )
        raise


if "review_status" not in st.session_state:
    st.session_state["review_status"] = "draft"

st.title("Define Research Question & Criteria")
st.write("Review status: ", st.session_state.review_status)

if "review_id" not in st.session_state:
    st.session_state["review_id"] = uuid.uuid4()
if "logger_extra_configured" not in st.session_state:
    logger.configure(extra={"review_id": st.session_state.review_id})
    st.session_state.logger_extra_configured = True

if "pico_population" not in st.session_state:
    st.session_state["pico_population"] = ""
if "pico_intervention" not in st.session_state:
    st.session_state["pico_intervention"] = ""
if "pico_comparison" not in st.session_state:
    st.session_state["pico_comparison"] = ""
if "pico_outcome" not in st.session_state:
    st.session_state["pico_outcome"] = ""
if "background" not in st.session_state:
    st.session_state["background"] = ""
if "research_question" not in st.session_state:
    st.session_state["research_question"] = ""
if "exclusion_criteria" not in st.session_state:
    st.session_state["exclusion_criteria"] = ""

if "llm_suggestions" not in st.session_state:
    st.session_state["llm_suggestions"] = ""
if "suggestions_agent" not in st.session_state:
    st.session_state["suggestions_agent"] = init_suggestion_agent()

col1, col2 = st.columns(2)

with col1:
    notification_widget = st.empty()
    review_model_widget = st.empty()
    st.subheader("Define Your Protocol")
    st.session_state.background = st.text_area(
        label="Brief introduction/Background",
        value=st.session_state.background,
        height=150,
        help="Provide context for the review, why it's needed, and the overall aim.",
    )
    st.session_state.research_question = st.text_area(
        label="Research Question",
        value=st.session_state.research_question,
        height=100,
        placeholder="E.g. 'Does intervention X improve outcome Y in population Z?'",
        help="The central question the systematic review aims to answer.",
    )

    st.divider()
    st.subheader("PICO Criteria")
    st.session_state.pico_population = st.text_area(
        "**P**opulation / Problem",
        value=st.session_state.pico_population,
        height=75,
        placeholder="E.g., Adults > 18 with type 2 diabetes",
        help="Describe the patient population or problem being addressed.",
    )
    st.session_state.pico_intervention = st.text_area(
        "**I**ntervention / Exposure",
        value=st.session_state.pico_intervention,
        height=75,
        placeholder="E.g., Metformin 1000mg daily",
        help="Describe the intervention, treatment, or exposure of interest.",
    )
    st.session_state.pico_comparison = st.text_area(
        "**C**omparison / Control",
        value=st.session_state.pico_comparison,
        height=75,
        placeholder="E.g., Placebo or standard care",
        help="Describe the alternative or control group being compared to the intervention.",
    )
    st.session_state.pico_outcome = st.text_area(
        "**O**utcome",
        value=st.session_state.pico_outcome,
        height=75,
        placeholder="E.g., Reduction in HbA1c levels",
        help="Describe the primary outcome(s) being measured.",
    )

    st.divider()
    st.session_state.exclusion_criteria = st.text_area(
        label="Explicit Exclusion Criteria (Optional)",
        value=st.session_state.get("exclusion_criteria", ""),
        placeholder="E.g., Non-human studies, Case reports, Studies not in English",
        height=100,
        help="Specific criteria used to exclude studies (complementary to PICO).",
    )

    validate_button = st.button("Validate & Suggest Improvements")

with col2:
    st.subheader("LLM Suggestions")
    suggestions_spinner = st.empty()
    if st.session_state["llm_suggestions"]:
        st.markdown(st.session_state["llm_suggestions"])
    else:
        st.write("No suggestions yet. Click the 'Validate' button to get feedback.")

if validate_button:
    with suggestions_spinner.container():
        st.info("Generating suggestions...")
        temp_review_data = SystematicReview(
            id=st.session_state.review_id,
            background=st.session_state.get("background", ""),
            research_question=st.session_state.get("research_question", ""),
            inclusion_criteria="",
            exclusion_criteria="",
        )
        with st.spinner("Analyzing protocol and suggesting PICO..."):
            suggestion_result: SuggestionResult = st.session_state[
                "suggestions_agent"
            ].get_suggestions(temp_review_data)

            pico_suggestions: PicosSuggestions | None = suggestion_result.get("pico")
            general_suggestions: str = suggestion_result.get(
                "raw_response", "No suggestions provided."
            )

            if pico_suggestions:
                if not st.session_state.pico_population and pico_suggestions.population:
                    st.session_state.pico_population = pico_suggestions.population
                    logger.info("Prefilled PICO Population.")
                if (
                    not st.session_state.pico_intervention
                    and pico_suggestions.intervention
                ):
                    st.session_state.pico_intervention = pico_suggestions.intervention
                    logger.info("Prefilled PICO Intervention.")
                if not st.session_state.pico_comparison and pico_suggestions.comparison:
                    st.session_state.pico_comparison = pico_suggestions.comparison
                    logger.info("Prefilled PICO Comparison.")
                if not st.session_state.pico_outcome and pico_suggestions.outcome:
                    st.session_state.pico_outcome = pico_suggestions.outcome
                    logger.info("Prefilled PICO Outcome.")

            st.session_state["llm_suggestions"] = general_suggestions
            st.rerun()

st.subheader("Save Protocol")
st.write(
    "Once you are happy with the question and criteria, click the 'Save Protocol' button."
)


if st.button("Save Protocol"):
    if not st.session_state.research_question:
        notification_widget.warning("Please define a research question")
        st.stop()
    elif not st.session_state.pico_population:
        notification_widget.warning("Please define the Population/Problem (PICO)")
        st.stop()
    else:
        review = build_review_model_from_pico(review_model_widget)
        if not review:
            st.session_state.review_status = "validation error"
            st.rerun()
        st.session_state.review = review
        st.session_state.review_status = f"saved at {datetime.now(UTC)}"

        with st.status("Saving protocol to database ...", expanded=True):
            try:
                review = persist_review(review)
                st.session_state.review = review
                st.write(
                    "Saved protocol to database, review status: ",
                    st.session_state.review_status,
                )
                st.json(review.model_dump(mode="json"))
                st.page_link(
                    "pages/search.py",
                    label="Next: Database Search",
                    icon=":material/search:",
                )
            except Exception as e:
                logger.exception("Failed to persist review to database.")
                st.error(f"Failed to save protocol: {e!r}")
