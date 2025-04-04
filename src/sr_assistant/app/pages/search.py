from __future__ import annotations

import uuid

import streamlit as st
from Bio import Entrez
from langchain_core.prompts import ChatPromptTemplate
from langchain_openai import ChatOpenAI
from loguru import logger
from pydantic import BaseModel, Field

from sr_assistant.core.models import SystematicReview
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
        description="PubMed query tailored to the review protocol",
    )

query_context = """\
Given the following review protocol, generate a valid PubMed query to search PubMed.

Background:
{background}

Research question:
{research_question}

Inclusion criteria:
{inclusion_criteria}

Exclusion criteria:
{exclusion_criteria}"""

query_draft_prompt = ChatPromptTemplate.from_messages([
    ("system", "You are an expert systematic review PubMed search query builder. Given the below protocol, generate a query to search PubMed."),
    ("user", query_context),
])

def init_query_chain():
    llm = ChatOpenAI(model="gpt-4o", temperature=0.0).with_structured_output(PubMedQuery)
    chain = query_draft_prompt | llm
    st.session_state.query_chain = chain

def get_query(review: SystematicReview) -> str:
    return st.session_state.query_chain.invoke(
        {
            "background": review.background,
            "research_question": review.research_question,
            "inclusion_criteria": review.inclusion_criteria,
            "exclusion_criteria": review.exclusion_criteria,
        }
    ).query

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
    st.markdown("This page allows you to search PubMed for relevant articles based on your systematic review protocol.")
    init_review_repository()
    if not review_id or "review_id" not in st.session_state:
        st.error("Please complete and save a review protocol first")
        st.stop()

    if "review" not in st.session_state:
        review = st.session_state.repo_review.get_by_id(review_id)
        st.session_state.review = review

    init_query_chain()
    repo = init_pubmed_repository()

    if "query_value" not in st.session_state or st.session_state.query_value is None:
        with st.spinner("Generating query..."):
            st.session_state.query_value = get_query(st.session_state.review)

    st.session_state.query_value = st.text_area(
        "Accept or modify PubMed Query",
        value=st.session_state.query_value,
        height=100,
        key="query",
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
                pmids = pubmed_search(st.session_state.query_value, st.session_state.max_results)
                if not pmids:
                    st.warning("No results found")
                    return

                st.write(f"Found {len(pmids)} results, fetching study details ...")

                # Fetch details
                records = pubmed_fetch_details(pmids)
                logger.info(f"Fetched {len(records['PubmedArticle'])} study details")
                st.write(f"Fetched {len(records['PubmedArticle'])} study details")

                # Store in Supabase
                results = repo.store_results(review_id, st.session_state.query_value, records)
                logger.info(f"Stored {len(results)} articles")
                st.success(f"Stored {len(results)} articles")
                status.update(label="Search complete", state="complete")

            except Exception as e:
                logger.exception("Search failed")
                status.update(label="Search failed", state="error")
                st.error(str(e))
                return

    if st.button("Clear search results"):
        repo.delete_by_review_id(review_id)
        st.success("Search results cleared")
    # Show existing results
    existing = repo.get_by_review_id(review_id)
    if not existing:
        return

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

    if pmid := st.selectbox("Select article:", [r.pmid for r in existing]):
        article = next(r for r in existing if r.pmid == pmid)
        st.subheader(article.title)
        st.text(f"{article.journal} ({article.year})")
        st.write(article.abstract)
        logger.info("Selected article: {!r}", article)
        st.json(article.model_dump(mode="json"), expanded=True)

    st.page_link("pages/screen_abstracts.py", label="Next: Screen Abstracts", icon=":material/arrow_forward:")


if "review" not in st.session_state:
    st.error("Please define a systematic review protocol first")
else:
    if "logger_extra_configured" not in st.session_state:
        logger.configure(extra={"review_id": st.session_state.review.id})
        st.session_state.logger_extra_configured = True
    search_page(review_id=st.session_state.review.id)
