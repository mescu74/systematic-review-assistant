from __future__ import annotations

import json
import uuid
from typing import cast

import streamlit as st
from dotenv import find_dotenv, load_dotenv
from langchain_core.runnables import RunnableConfig
from supabase import Client

from sr_assistant.core.models.base import Review
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


# Set up page
st.title("Define Research Question & Criteria")

# Initialize session state variables
# HACK: This should come from the Review model, but because of streamlit and needing to
# keep this constant, we've to initialize the reviews_id here. It's like a session ID.
if "review_id" not in st.session_state:
    st.session_state["review_id"] = uuid.uuid4()
if "background" not in st.session_state:
    st.session_state["background"] = ""
if "question" not in st.session_state:
    st.session_state["question"] = ""
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
    st.session_state.question = st.text_area(
        label="Research Question",
        value=st.session_state.question,
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
    if st.session_state["llm_suggestions"]:
        st.markdown(st.session_state["llm_suggestions"])
    else:
        st.write("No suggestions yet. Click the 'Validate' button to get feedback.")

# LLM Interaction
if validate_button:
    review = Review(
        id=st.session_state.review_id,
        background=st.session_state.background,
        question=st.session_state.question,
        inclusion_criteria=st.session_state.inclusion_criteria,
        exclusion_criteria=st.session_state.exclusion_criteria,
    )
    st.session_state.review = review
    st.session_state["llm_suggestions"] = st.session_state[
        "suggestions_agent"
    ].get_suggestions(review)
    st.rerun()

st.subheader("Finalize & Save")
st.write(
    "Once you are happy with the question and criteria, click the 'Save Protocol' button."
)

save_button = st.button("Save Protocol")
if save_button:
    if "review" not in st.session_state:
        st.session_state.review = Review(
            id=st.session_state.review_id,
            background=st.session_state.background,
            question=st.session_state.question,
            inclusion_criteria=st.session_state.inclusion_criteria,
            exclusion_criteria=st.session_state.exclusion_criteria,
        )
    review = st.session_state.review
    supabase: Client = st.session_state.supabase
    # handle UUIDs, let Supabase be source for time
    payload = json.loads(review.model_dump_json(exclude={"created_at", "updated_at"}))
    # TODO: handle updates/upsert
    ret = supabase.table("reviews").insert(payload).execute()
    saved_review = Review.model_validate(ret.data[0])
    st.session_state.review = saved_review

    st.success("Criteria saved")
    # we get int overflows if using the model directly
    st.json(json.loads(st.session_state.review.model_dump_json()))
