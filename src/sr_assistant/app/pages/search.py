from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

import streamlit as st
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI
from loguru import logger
from pydantic import BaseModel, Field

from sr_assistant.app.database import session_factory
from sr_assistant.app.services import SearchService
from sr_assistant.core.models import CriteriaFramework, SystematicReview
from sr_assistant.core.repositories import (
    SystematicReviewRepository,
)

# Define TypeVar for structured output types
T = TypeVar("T")
StructuredLLMOutput = Callable[[dict[str, Any]], T]


class PubMedQuery(BaseModel):
    """PubMed query based on review protocol."""

    query: str = Field(
        ...,
        title="PubMed Query",
        description="PubMed query tailored to the review protocol, using PICO elements.",
    )


query_context = """\
Given the following review protocol, generate a valid PubMed query suitable for an initial search in PubMed. Use PubMed's field tags (like [tiab], [mh], [pt]) and boolean operators (AND, OR, NOT) effectively. Combine the PICO elements logically.

Background:
{background}

Research question:
{research_question}

PICO Criteria:
Population: {population}
Intervention: {intervention}
Comparison: {comparison}
Outcome: {outcome}

Explicit Exclusion Criteria (apply using NOT if appropriate):
{exclusion_criteria}"""

query_draft_prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            "You are an expert systematic review PubMed search query builder. Given the PICO-based protocol below, generate an effective PubMed query.",
        ),
        ("user", query_context),
    ]
)


def init_query_chain():
    llm = ChatOpenAI(model="gpt-4o", temperature=0.0).with_structured_output(
        PubMedQuery
    )
    chain = query_draft_prompt | llm
    st.session_state.query_chain = chain


def get_query(review: SystematicReview) -> str:
    """Generates a PubMed query string based on the review's PICO criteria."""
    if (
        review.criteria_framework != CriteriaFramework.PICO
        or not review.criteria_framework_answers
    ):
        logger.warning(
            "Review criteria framework is not PICO or answers are missing. Falling back to basic query gen."
        )
        fallback_prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    "Generate a basic PubMed query based on the research question.",
                ),
                ("user", "Research Question: {research_question}"),
            ]
        )
        fallback_chain = fallback_prompt | ChatOpenAI(
            model="gpt-4o", temperature=0.0
        ).with_structured_output(PubMedQuery)
        result_object = fallback_chain.invoke(
            {"research_question": review.research_question}
        )
        # LangChain's result type is complex, but we know it implements a 'query' attribute
        return result_object.query  # type: ignore

    pico = review.criteria_framework_answers
    logger.info(f"Generating query from PICO: {pico}")
    result_object = st.session_state.query_chain.invoke(
        {
            "background": review.background or "",
            "research_question": review.research_question,
            "population": pico.get("population", ""),
            "intervention": pico.get("intervention", ""),
            "comparison": pico.get("comparison", ""),
            "outcome": pico.get("outcome", ""),
            "exclusion_criteria": review.exclusion_criteria or "",
        }
    )
    # LangChain's result type is complex, but we know it implements a 'query' attribute
    return result_object.query  # type: ignore


def gen_query_cb() -> None:
    st.session_state.query_value = get_query(st.session_state.review)


def init_review_repository() -> SystematicReviewRepository:
    if "repo_review" not in st.session_state:
        repo = SystematicReviewRepository()
        st.session_state.repo_review = repo
    return st.session_state.repo_review


def search_page() -> None:
    """PubMed search page."""
    st.title("PubMed Search")
    st.markdown(
        "This page allows you to search PubMed for relevant articles based on your systematic review protocol."
    )

    # 1. Ensure a review_id is in session state (presumably set by protocol page or navigation)
    if "review_id" not in st.session_state or not st.session_state.review_id:
        st.error(
            "No active review selected. Please go back to the protocol page and select or create a review."
        )
        st.stop()

    active_review_id = st.session_state.review_id

    # 2. Load the review object if not already loaded or if ID changed
    if (
        "review" not in st.session_state
        or st.session_state.review.id != active_review_id
    ):
        repo = init_review_repository()  # Ensure repo is initialized
        # The get_by_id method from BaseRepository expects a session.
        # We need to create/get a session here if init_review_repository doesn't provide one directly for this call.
        # For simplicity, assuming init_review_repository handles its own session or service does.
        # However, repository methods themselves expect an explicit session.
        # Let's get a session from the factory if the repo isn't already session-scoped in its methods.
        # The init_review_repository() stores the repo in session_state, but its methods still need a session.
        # This part of the logic assumes a service layer would handle session scope.
        # Given direct repo usage here, we must provide a session.
        with session_factory() as session:  # Get a session for this operation
            review_object = repo.get_by_id(session, active_review_id)

        if not review_object:
            st.error(
                f"Could not load the active review with ID: {active_review_id}. Please return to protocol page."
            )
            st.session_state.pop(
                "review", None
            )  # Clear potentially stale review object
            st.session_state.pop("review_id", None)  # Clear problematic review_id
            st.stop()
        st.session_state.review = review_object

    current_review = st.session_state.review

    # Ensure logger is configured with the now-validated review_id from current_review
    if (
        "logger_extra_configured_search" not in st.session_state
        or st.session_state.get("logger_review_id_search") != current_review.id
    ):
        logger.configure(extra={"review_id": current_review.id})
        st.session_state.logger_extra_configured_search = True
        st.session_state.logger_review_id_search = current_review.id

    init_query_chain()
    search_service = SearchService()

    # Initial query generation using the loaded review
    if "query_value" not in st.session_state or st.session_state.query_value is None:
        with st.spinner("Generating initial query based on PICO..."):
            st.session_state.query_value = get_query(current_review)  # Pass the review

    st.session_state.query_value = st.text_area(
        "Accept or modify PubMed Query",
        value=st.session_state.query_value,
        height=100,
        key="query",  # Keep key if needed for direct access, though callbacks use query_value
    )
    col1, col2, col3 = st.columns(3)
    with col1:
        submitted = st.button("Search", key="search_button")
    with col2:
        st.session_state.max_results = st.slider(
            "Max results",
            min_value=1,
            max_value=500,
            value=50,
            step=5,
        )
    with col3:
        st.button("Generate query", on_click=gen_query_cb)

    if submitted:
        with st.status("Searching...", expanded=True) as status:
            try:
                logger.info("Searching PubMed: {!r}", st.session_state.query_value)
                st.write("Searching PubMed and storing results...")

                # Call SearchService
                search_results = search_service.search_pubmed_and_store_results(
                    review_id=current_review.id,
                    query=st.session_state.query_value,
                    max_results=st.session_state.max_results,
                )

                if not search_results:
                    st.warning("No results found")
                    status.update(label="No results found", state="complete")
                else:
                    logger.info(f"Search returned {len(search_results)} articles")
                    st.success(
                        f"Search complete. Found {len(search_results)} articles."
                    )
                    # Store results in session state for display
                    st.session_state.search_results = search_results

                status.update(label="Search complete", state="complete")

            except Exception as e:
                logger.exception("Search failed")
                status.update(label="Search failed", state="error")
                st.error(str(e))
                # Do not return here, allow results display below

    if st.button("Clear search results"):
        # NOTE: This currently only clears results from the session state display.
        # It does not delete any data from the database.
        # Future improvement: Add a SearchService.delete_results_by_review_id method
        # that properly manages the database transaction.
        st.session_state.search_results = []
        logger.info(f"Cleared search results for review_id: {active_review_id}")
        st.success("Search results cleared from display. Re-search to fetch again.")
        st.rerun()

    # Show existing results - now from st.session_state.search_results populated by the service call
    if "search_results" not in st.session_state:
        st.session_state.search_results = []  # Initialize if not present

    existing = st.session_state.search_results

    if not existing:
        st.info("No search results available. Perform a search to fetch articles.")
        # st.stop() # Allow page to render even with no results

    st.divider()
    df_data = [
        {
            "Source ID": r.source_id,
            "DOI": r.doi,
            "Title": r.title,
            "Journal": r.journal,
            "Year": r.year,
        }
        for r in existing
    ]

    st.dataframe(df_data, use_container_width=True)

    if existing:
        if selected_source_id := st.selectbox(
            "Select article by Source ID:", [r.source_id for r in existing]
        ):
            article = next(
                (r for r in existing if r.source_id == selected_source_id), None
            )
            if article:
                st.subheader(article.title)
                st.text(f"{article.journal} ({article.year})")
                st.write(article.abstract)
                logger.info("Selected article: {!r}", article)
                st.json(article.model_dump(mode="json"), expanded=True)

    st.page_link(
        "pages/screen_abstracts.py",
        label="Next: Screen Abstracts",
        icon=":material/arrow_forward:",
    )


if "review" not in st.session_state:
    st.error("Please define a systematic review protocol first")
else:
    # The main call to search_page no longer needs review_id parameter
    # It will pick up review_id from session_state internally
    # if "logger_extra_configured" not in st.session_state: # This global logger config is problematic
    #     logger.configure(extra={"review_id": st.session_state.review.id})
    #     st.session_state.logger_extra_configured = True
    search_page()
