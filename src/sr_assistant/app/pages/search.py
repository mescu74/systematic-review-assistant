from __future__ import annotations

import uuid

import streamlit as st
from Bio import Entrez
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI
from loguru import logger
from pydantic import BaseModel, Field

from sr_assistant.core.models import CriteriaFramework, SystematicReview
from sr_assistant.core.repositories import (
    PubMedResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.step2.pubmed_integration import pubmed_fetch_details, pubmed_search


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
        # Handle cases where PICO isn't set (though it should be from protocol page)
        logger.warning(
            "Review criteria framework is not PICO or answers are missing. Falling back to basic query gen."
        )
        # Fallback prompt (optional, or could raise error)
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
        result_dict = fallback_chain.invoke(
            {"research_question": review.research_question}
        )
        return result_dict.get("query", "")

    pico = review.criteria_framework_answers
    logger.info(f"Generating query from PICO: {pico}")
    result_dict = st.session_state.query_chain.invoke(
        {
            "background": review.background or "",  # Ensure defaults if None
            "research_question": review.research_question,
            "population": pico.get("population", ""),
            "intervention": pico.get("intervention", ""),
            "comparison": pico.get("comparison", ""),
            "outcome": pico.get("outcome", ""),
            "exclusion_criteria": review.exclusion_criteria or "",  # Ensure default
        }
    )
    return result_dict.get("query", "")


def gen_query_cb() -> None:
    st.session_state.query_value = get_query(st.session_state.review)


def expand_query_cb() -> None:
    Entrez.email = st.session_state.config.NCBI_EMAIL
    Entrez.api_key = st.session_state.config.NCBI_API_KEY.get_secret_value()

    with Entrez.esearch(db="pubmed", term=st.session_state.query_value) as handle:
        result = Entrez.read(handle)  # type: ignore
        res = result.get("QueryTranslation", "")  # type: ignore
    st.session_state.query_value = res


def init_review_repository() -> SystematicReviewRepository:
    if "repo_review" not in st.session_state:
        repo = SystematicReviewRepository()
        st.session_state.repo_review = repo
    return st.session_state.repo_review


def init_pubmed_repository() -> PubMedResultRepository:
    if "repo_pubmed" not in st.session_state:
        repo = PubMedResultRepository()
        st.session_state.repo_pubmed = repo
    return st.session_state.repo_pubmed


def search_page(review_id: uuid.UUID | None = None) -> None:
    """PubMed search page."""
    st.title("PubMed Search")
    st.markdown(
        "This page allows you to search PubMed for relevant articles based on your systematic review protocol."
    )
    init_review_repository()
    if (
        not review_id
        or "review_id" not in st.session_state
        or st.session_state.review_id != review_id
    ):
        # Added check if the review_id in session matches the one passed (or expected)
        st.error(
            "Invalid or missing review context. Please go back to the protocol page."
        )
        st.stop()

    if "review" not in st.session_state or st.session_state.review.id != review_id:
        review = st.session_state.repo_review.get_by_id(review_id)
        if not review:
            st.error(f"Could not load review with ID: {review_id}")
            st.stop()
        st.session_state.review = review
    # Ensure the loaded review is used
    current_review = st.session_state.review

    init_query_chain()
    repo = init_pubmed_repository()

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
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        submitted = st.button("Search")
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
    with col4:
        st.button("Expand query", on_click=expand_query_cb)

    if submitted:
        with st.status("Searching...", expanded=True) as status:
            try:
                # Search PubMed
                logger.info("Searching PubMed: {!r}", st.session_state.query_value)
                st.write("Searching PubMed...")
                pmids = pubmed_search(
                    st.session_state.query_value, st.session_state.max_results
                )
                if not pmids:
                    st.warning("No results found")
                    status.update(
                        label="No results found", state="complete"
                    )  # Update status
                    st.stop()  # Stop further processing

                st.write(f"Found {len(pmids)} results, fetching study details ...")

                # Fetch details
                records = pubmed_fetch_details(pmids)
                logger.info(f"Fetched {len(records['PubmedArticle'])} study details")
                st.write(f"Fetched {len(records['PubmedArticle'])} study details")

                # Store in Supabase
                results = repo.store_results(
                    review_id, st.session_state.query_value, records
                )
                logger.info(f"Stored {len(results)} articles")
                st.success(f"Stored {len(results)} articles")
                status.update(label="Search complete", state="complete")

            except Exception as e:
                logger.exception("Search failed")
                status.update(label="Search failed", state="error")
                st.error(str(e))
                # Do not return here, allow results display below

    if st.button("Clear search results"):
        repo.delete_by_review_id(review_id)
        st.session_state.pubmed_results = []  # Clear local state too
        st.success("Search results cleared")
        st.rerun()  # Rerun to reflect cleared results

    # Show existing results
    existing = repo.get_by_review_id(review_id)
    if not existing:
        st.info("No search results available for this review yet.")
        st.stop()  # Stop if no results

    st.session_state.pubmed_results = existing

    st.divider()
    df_data = [
        {
            "PMID": r.pmid,
            "DOI": r.doi,
            "PMC": r.pmc,
            "Title": r.title,
            "Journal": r.journal,
            "Year": r.year,
            "Query": r.query,
        }
        for r in existing
    ]
    st.dataframe(df_data, use_container_width=True)

    if existing:
        if pmid := st.selectbox("Select article:", [r.pmid for r in existing]):
            article = next((r for r in existing if r.pmid == pmid), None)
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
    if "logger_extra_configured" not in st.session_state:
        logger.configure(extra={"review_id": st.session_state.review.id})
        st.session_state.logger_extra_configured = True
    search_page(review_id=st.session_state.review.id)
