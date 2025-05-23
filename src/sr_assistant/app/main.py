# sr_assistant/app/main.py
from __future__ import annotations

import streamlit as st
from dotenv import find_dotenv, load_dotenv
from loguru import logger
from supabase import create_client

# If done properly, the DI container would assemble and wire everything together
from sr_assistant.app.config import get_settings

st.set_page_config(layout="wide")

from sr_assistant.app.database import (  # noqa: E402 [page config has to come before]
    asession_factory,
    async_engine,
    engine,
    session_factory,
)
from sr_assistant.app.logging import configure_logging  # noqa: E402

load_dotenv(find_dotenv(), override=True)


def main() -> None:  # noqa: C901
    # this gets called on every rerun ...
    if "logging_initialized" not in st.session_state:
        configure_logging()
        st.session_state.logging_initialized = True

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
        "pages/screen_abstracts.py",
        title="Abstracts Screening",
        icon=":material/preview:",
    )
    fulltext_screening_page = st.Page(
        "pages/screen_fulltext.py",
        title="Fulltext Screening",
        icon=":material/preview:",
    )

    # Benchmark pages
    human_benchmark_page = st.Page(
        "pages/human_benchmark_page.py",
        title="Protocol Review",
        icon=":material/description:",
    )
    benchmark_tool_page = st.Page(
        "pages/benchmark_tool.py",
        title="Run Benchmark",
        icon=":material/analytics:",
    )

    if st.session_state.logged_in:
        pg = st.navigation(
            {
                "Account": [logout_page],
                "Protocol": [protocol_page],
                "Search Strategy": [search_page],
                "Screening": [abstracts_page, fulltext_screening_page],
                "Benchmark": [human_benchmark_page, benchmark_tool_page],
            }
        )
        if "config" not in st.session_state:
            st.session_state.config = get_settings()
        if "supabase" not in st.session_state:
            st.session_state.supabase = create_client(
                supabase_url=st.session_state.config.SUPABASE_URL,
                supabase_key=st.session_state.config.SUPABASE_KEY.get_secret_value(),
            )
        if "engine" not in st.session_state:
            st.session_state.engine = engine
        if "async_engine" not in st.session_state:
            st.session_state.async_engine = async_engine
        if "session_factory" not in st.session_state:
            st.session_state.session_factory = session_factory
        if "asession_factory" not in st.session_state:
            st.session_state.asession_factory = asession_factory

        # Remove review_id dependent logger configuration from main.py
        # Individual pages will handle their review-specific logging context.
        # if "logger_extra_configured" not in st.session_state:
        #     # This was problematic as st.session_state.review might not be valid here
        #     # logger.configure(extra={"review_id": st.session_state.review.id})
        #     st.session_state.logger_extra_configured = True
        logger.info(
            "main() initialized"
        )  # This log will no longer try to add a potentially invalid review_id
    else:
        pg = st.navigation([login_page])

    pg.run()


if __name__ == "__main__":
    main()
