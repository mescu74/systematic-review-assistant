from __future__ import annotations

import streamlit as st

from sr_assistant.app.agents import abstract_screening_chain
from sr_assistant.core.repositories import (
    PubMedRepository,
)

if "abstract_screening_chain" not in st.session_state:
    st.session_state.abstract_screening_chain = abstract_screening_chain
if "review_id" not in st.session_state:
    st.error("Please select a review first")
    st.stop()

review_id = st.session_state.review_id
if "pubmed_repository" not in st.session_state:
    st.session_state.pubmed_repository = PubMedRepository(st.session_state.supabase)
# TODO: repo for abstract screening results, and SQLModel
# - read search results from pubmed repo
# - batch invoke the chain, map responses to search results
#   here additional batching needed to update ui inbetween
# - store results in a repo
#
