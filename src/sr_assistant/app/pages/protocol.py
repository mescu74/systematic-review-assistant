from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import cast

import streamlit as st
from dotenv import find_dotenv, load_dotenv
from langchain_core.runnables import RunnableConfig
from loguru import logger
from sqlmodel import select
from streamlit.delta_generator import DeltaGenerator

from sr_assistant.core.models import SystematicReview
from sr_assistant.step1.suggestion_agent import SuggestionAgent

load_dotenv(find_dotenv())


def init_agent_config() -> RunnableConfig:
    if "agent_config" in st.session_state:
        return cast(RunnableConfig, st.session_state.agent_config)

    config = RunnableConfig(
        recursion_limit=1000, configurable={"thread_id": str(uuid.uuid4().hex)}
    )
    st.session_state.agent_config = config
    return config


def init_suggestion_agent() -> SuggestionAgent:
    """Initialize the suggestion agent."""
    return SuggestionAgent(model="gpt-4o", temperature=0.0)


@logger.catch(Exception, message="Error building review model")
def build_review_model(response_widget: DeltaGenerator) -> SystematicReview:
    review = SystematicReview(
        id=st.session_state.review_id,
        background=st.session_state.background,
        research_question=st.session_state.research_question,
        inclusion_criteria=st.session_state.inclusion_criteria,
        exclusion_criteria=st.session_state.exclusion_criteria,
    )
    response_widget.success("SystematicReview model built successfully")
    return review


def persist_review(model: SystematicReview) -> SystematicReview:
    with st.session_state.session_factory() as session:
        # Check if record exists in DB
        stmt = select(SystematicReview).where(SystematicReview.id == model.id)
        db_model = session.exec(stmt).first()

        if db_model:
            # Update existing record
            db_model.background = model.background
            db_model.research_question = model.research_question
            db_model.inclusion_criteria = model.inclusion_criteria
            db_model.exclusion_criteria = model.exclusion_criteria
            persisted_model = db_model
        else:
            # Add as new record
            session.add(model)
            persisted_model = model

        session.commit()
        session.refresh(persisted_model)
        return persisted_model


if "review_status" not in st.session_state:
    st.session_state["review_status"] = "draft"

# Set up page
st.title("Define Research Question & Criteria")
st.write("Review status: ", st.session_state.review_status)

# Initialize session state variables
# HACK: This should come from the Review model, but because of streamlit and needing to
# keep this constant, we've to initialize the reviews_id here. It's like a session ID.
if "review_id" not in st.session_state:
    st.session_state["review_id"] = uuid.uuid4()
if "background" not in st.session_state:
    st.session_state["background"] = ""
if "research_question" not in st.session_state:
    st.session_state["research_question"] = ""
if "inclusion_criteria" not in st.session_state:
    st.session_state["inclusion_criteria"] = ""
if "exclusion_criteria" not in st.session_state:
    st.session_state["exclusion_criteria"] = ""
if "llm_suggestions" not in st.session_state:
    st.session_state["llm_suggestions"] = ""
if "suggestions_agent" not in st.session_state:
    st.session_state["suggestions_agent"] = init_suggestion_agent()

# Layout for input fields and LLM feedback
col1, col2 = st.columns(2)

with col1:
    notification_widget = st.empty()
    review_model_widget = st.empty()
    st.write("TODO: PICO(T/S), FINER, chat, ...")
    st.subheader("Define Your Protocol")
    st.session_state.background = st.text_area(
        label="Brief introduction to the subject of the review, including rationale for undertaking the review and overall aim.",
        value=st.session_state.background,
        height=200,
        # key="background",
    )
    # st.markdown(
    #    "> The most common pitfalls when developing research questions are that the questions incorporate the methods or the study’s expected outcomes ([Mayo, Asano, and Pamela Barbic 2013](https://ehsanx.github.io/Scientific-Writing-for-Health-Research/research-question.html#ref-mayo2013research)). Furthermore, the clarity of the research question can be impeded by the lack of a clear parameter to assess the relationship or association between exposure and outcome ([Mayo, Asano, and Pamela Barbic 2013](https://ehsanx.github.io/Scientific-Writing-for-Health-Research/research-question.html#ref-mayo2013research))."
    # )

    # NOTE: In order for the text fields to retain content when navigating between pages,
    # we can't set the key but have to assing the widget to state (for some reason).
    st.session_state.research_question = st.text_area(
        label="Research Question",
        value=st.session_state.research_question,
        # key="question",
        placeholder="E.g. 'Does intervention X improve outcome Y in population Z?'",
    )
    st.session_state.inclusion_criteria = st.text_area(
        label="Inclusion Criteria",
        value=st.session_state.inclusion_criteria,
        # key="inclusion_criteria",
        placeholder="E.g. Adults >18, Intervention X, Outcome Y, RCTs, English",
    )
    st.session_state.exclusion_criteria = st.text_area(
        label="Exclusion Criteria",
        value=st.session_state.exclusion_criteria,
        # key="exclusion_criteria",
        placeholder="E.g. Children <18, Non-human studies, Case reports",
    )
    # keywords_input = st.text_input("Search Keywords (separated by commas)")

    ## Split the input string into a list and store in session state
    # if keywords_input:
    #    st.session_state.keywords = [kw.strip() for kw in keywords_input.split(",")]
    # else:
    #    st.session_state.keywords = []

    validate_button = st.button("Validate & Suggest Improvements")

with col2:
    st.subheader("LLM Suggestions")
    suggestions_spinner = st.empty()
    if st.session_state["llm_suggestions"]:
        st.markdown(st.session_state["llm_suggestions"])
    else:
        st.write("No suggestions yet. Click the 'Validate' button to get feedback.")

# LLM Interaction
if validate_button:
    with suggestions_spinner.container():
        st.info("Validating ...")
        review = build_review_model(review_model_widget)
        with st.spinner():
            st.session_state["llm_suggestions"] = st.session_state[
                "suggestions_agent"
            ].get_suggestions(review)
            st.rerun()

st.subheader("Save Protocol")
st.write(
    "Once you are happy with the question and criteria, click the 'Save Protocol' button."
)


if st.button("Save Protocol"):
    if not st.session_state.research_question:
        notification_widget.warning("Please define a research question")
        st.stop()
    elif not st.session_state.inclusion_criteria:
        notification_widget.warning("Please define inclusion criteria")
        st.stop()
    elif not st.session_state.exclusion_criteria:
        notification_widget.warning("Please define exclusion criteria")
        st.stop()
    else:
        review = build_review_model(review_model_widget)
        if not review:
            st.session_state.review_status = "validation error"
            st.rerun()
        st.session_state.review = review
        st.session_state.review_status = f"saved at {datetime.now(tz=timezone.utc)}"

        with st.status("Saving protocol to database ...", expanded=True):
            review = persist_review(review)
            st.session_state.review = review
            st.write(
                "Saved protocol to database, review status: ",
                st.session_state.review_status,
            )
            st.json(review.model_dump(mode="json"))
