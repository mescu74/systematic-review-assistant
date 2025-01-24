from __future__ import annotations

from uuid import UUID

import streamlit as st

from sr_assistant.core.repositories.pubmed import PubMedRepository
from sr_assistant.step2.pubmed_integration import pubmed_fetch_details, pubmed_search


def search_page(review_id: UUID | None = None) -> None:
    """PubMed search page."""
    if not review_id:
        st.error("Please select a review first")
        return

    repo = PubMedRepository(st.session_state.supabase)

    with st.form("pubmed_search"):
        query = st.text_input(
            "Enter PubMed Query:",
            value="(cancer) AND (immunotherapy) AND (clinical trial)",
        )
        max_results = st.number_input("Max Results", value=500, min_value=1)
        submitted = st.form_submit_button("Search")

    if submitted:
        with st.status("Searching...", expanded=True) as status:
            try:
                # Search PubMed
                st.write("Searching PubMed...")
                pmids = pubmed_search(query, max_results)
                if not pmids:
                    st.warning("No results found")
                    return

                st.write(f"Found {len(pmids)} results, fetching study details ...")

                # Fetch details
                records = pubmed_fetch_details(pmids)
                st.write(f"Fetched {len(records["PubmedArticle"])} study details")

                # Store in Supabase
                results = repo.store_results(review_id, query, records)

                st.success(f"Stored {len(results)} articles")
                status.update(label="Search complete", state="complete")

            except Exception as e:
                status.update(label="Search failed", state="error")
                st.error(str(e))
                return

    # Show existing results
    existing = repo.get_review_results(review_id)
    if not existing:
        return

    st.divider()
    df_data = [{
        "PMID": r.pmid,
        "Title": r.title,
        "Journal": r.journal,
        "Year": r.year,
        "Query": r.query,
    } for r in existing]
    st.dataframe(df_data, use_container_width=True)

    if pmid := st.selectbox("Select article:", [r.pmid for r in existing]):
        article = next(r for r in existing if r.pmid == pmid)
        st.subheader(article.title)
        st.text(f"{article.journal} ({article.year})")
        st.write(article.abstract)

if "review_id" not in st.session_state:
    st.error("Please create a review first")
else:
    search_page(review_id=st.session_state.review_id)
