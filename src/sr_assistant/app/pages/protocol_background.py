from __future__ import annotations

import streamlit as st

st.title("Background to Review")

st.text_area(
    label="Brief introduction to the subject of the review, including rationale for undertaking the review and overall aim.",
    height=300,
    key="protocol_bg",
)
st.page_link(
    "pages/protocol_question.py",
    label="Next: Research question and eligibility criteria",
    icon=":material/question_mark:",
)
