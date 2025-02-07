from __future__ import annotations

import sys

# HACK: streamlit hacks to make this work, otherwise the `app` pkg doesn't resolve
# fmt: off
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
# fmt: on

import streamlit as st
from sqlalchemy.orm import sessionmaker
from sqlmodel import Session as SQLModelSession
from sqlmodel import create_engine
from supabase import create_client

# If done properly, the DI container would assemble and wire everything together
from app.config import get_settings

st.set_page_config(layout="wide")


def main() -> None:
    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False

    def login() -> None:
        if st.button("Log in"):
            st.session_state.logged_in = True
            st.rerun()

    def logout() -> None:
        if st.button("Log out"):
            st.session_state.logged_in = False
            st.rerun()

    login_page = st.Page(login, title="Log in", icon=":material/login:")
    logout_page = st.Page(logout, title="Log out", icon=":material/logout:")

    # bg_page = st.Page(
    #    "pages/protocol_background.py",
    #    title="Background",
    #    icon=":material/wallpaper:",
    #    default=True,
    # )
    # question_page = st.Page(
    #    "pages/protocol_question.py",
    #    title="Question",
    #    icon=":material/question_mark:",
    # )
    protocol_page = st.Page(
        "pages/protocol.py",
        title="Review Protocol",
        icon=":material/wallpaper:",
    )
    search_page = st.Page(
        "pages/search.py",
        title="PubMed Search",
        icon=":material/search:",
    )
    abstracts_page = st.Page(
        "pages/screening_abstracts.py",
        title="Abstracts Screening",
        icon=":material/preview:",
    )

    if st.session_state.logged_in:
        pg = st.navigation(
            {
                "Account": [logout_page],
                "Protocol": [protocol_page],
                "Search Strategy": [search_page],
                "Screening": [abstracts_page],
            }
        )
        st.session_state.config = get_settings()
        if "supabase" not in st.session_state:
            st.session_state.supabase = create_client(
                supabase_url=st.session_state.config.SUPABASE_URL,
                supabase_key=st.session_state.config.SUPABASE_KEY,
            )
        if "engine" not in st.session_state:
            st.session_state.engine = create_engine(
                str(st.session_state.config.DATABASE_URL)
            )
        if "session_factory" not in st.session_state:
            st.session_state.session_factory = sessionmaker(
                class_=SQLModelSession, autoflush=False, bind=st.session_state.engine
            )
    else:
        pg = st.navigation([login_page])

    pg.run()


if __name__ == "__main__":
    main()
