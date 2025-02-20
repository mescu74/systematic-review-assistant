from __future__ import annotations

import typing as t
import uuid
from statistics import mean

import pandas as pd
import streamlit as st
from loguru import logger

import sr_assistant.app.utils as ut
from sr_assistant.app.agents.screening_agents import (
    # ScreeningResponse,
    ScreeningError,
    ScreeningResult,
    ScreeningStrategyType,
    screen_abstracts_batch,
    screen_abstracts_chain,
)
from sr_assistant.core.models import (
    PubMedResult,
    SystematicReview,
)
from sr_assistant.core.repositories_old import (
    PubMedRepository,
    ScreeningAbstractRepository,
)
from sr_assistant.core.schemas import ScreeningDecisionType

if t.TYPE_CHECKING:
    from supabase import Client


def init_pubmed_repository(supabase: Client) -> PubMedRepository:
    if "repo_pubmed" not in st.session_state:
        repo = PubMedRepository(supabase)
        st.session_state.repo_pubmed = repo
    return st.session_state.repo_pubmed


def init_screen_abstracts_repository(supabase: Client) -> ScreeningAbstractRepository:
    if "repo_screen_abstracts" not in st.session_state:
        repo = ScreeningAbstractRepository(supabase)
        st.session_state.repo_screen_abstracts = repo
    return st.session_state.repo_screen_abstracts


def init_screen_abstracts_chain():
    if "screen_abstracts_chain" not in st.session_state:
        st.session_state.screen_abstracts_chain = screen_abstracts_chain
    return st.session_state.screen_abstracts_chain


# @st.fragment(run_every=st.session_state.run_every)
# def show_latest_data():
#    last_timestamp = st.session_state.data.index[-1]
#    st.session_state.data = pd.concat(
#        [st.session_state.data, get_recent_data(last_timestamp)]
#    )
#    st.session_state.data = st.session_state.data[-100:]
#    st.line_chart(st.session_state.data)


def is_screening_result(result: t.Any) -> bool:
    return isinstance(result, ScreeningResult)


def is_screening_error(result: t.Any) -> bool:
    return isinstance(result, ScreeningError)


@st.fragment
def render_df() -> None:
    """Render results dataframe.

    Todo:
        - optimise this
        - make editable (Decision column, maybe a notes column?)
        - highlight disagreements
        - add second summary df with computed values like mean confidence, percentiles, etc.
    """
    # exclusion categories, compute metrics
    results = t.cast(
        list[tuple[PubMedResult, ScreeningResult]],
        st.session_state.screen_abstracts_results,
    )
    df_data = [
        {
            "PMID": pubmed_result.pmid,
            # "Title": pubmed_result.title,
            "Strategy": screening_result.screening_strategy,
            "Decision": screening_result.decision,
            "Confidence": screening_result.confidence_score,
            "Model": screening_result.model_name,
            "Rationale": screening_result.rationale,
            "Extracted quotes": screening_result.extracted_quotes,
            "Result ID": str(screening_result.id),
        }
        for pubmed_result, screening_result in results
    ]
    if len(df_data) > 0:
        st.subheader("Results")
        st.dataframe(df_data, use_container_width=True, hide_index=True)

        if pmid := st.selectbox("View details:", [pm.pmid for pm, _ in results]):
            pubmed_result, conservative_screening_result = next(
                (pubmed_result, screening_result)
                for pubmed_result, screening_result in results
                if pubmed_result.pmid == pmid
                and screening_result.screening_strategy
                == ScreeningStrategyType.CONSERVATIVE
            )
            pubmed_result, comprehensive_screening_result = next(
                (pubmed_result, screening_result)
                for pubmed_result, screening_result in results
                if pubmed_result.pmid == pmid
                and screening_result.screening_strategy
                == ScreeningStrategyType.COMPREHENSIVE
            )
            st.subheader(f"{pubmed_result.pmid}: {pubmed_result.title}")
            st.text(f"{pubmed_result.journal} ({pubmed_result.year})")
            st.markdown(
                f"[PubMedResult](https://pubmed.ncbi.nlm.nih.gov/{pubmed_result.pmid})"
            )
            st.json(pubmed_result.model_dump(mode="json"))
            if conservative_screening_result:
                st.markdown(
                    f"[Conservative Screening Result](https://pubmed.ncbi.nlm.nih.gov/{pubmed_result.pmid})"
                )
                st.json(conservative_screening_result.model_dump(mode="json"))
            else:
                st.warning("No conservative screening result")
            if comprehensive_screening_result:
                st.markdown(
                    f"[Comprehensive Screening Result](https://pubmed.ncbi.nlm.nih.gov/{pubmed_result.pmid})"
                )
                st.json(comprehensive_screening_result.model_dump(mode="json"))
            else:
                st.warning("No comprehensive screening result")

    df_conflicts_data = [
        {
            "PMID": pm.pmid,
            "Conservative Decision": conservative_result.decision,
            "Conservative confidence": conservative_result.confidence_score,
            "Conservative Rationale": conservative_result.rationale,
            "Comprehensive Decision": comprehensive_result.decision,
            "Comprehensive confidence": comprehensive_result.confidence_score,
            "Comprehensive Rationale": comprehensive_result.rationale,
        }
        for pm, conservative_result, comprehensive_result in st.session_state.get(
            "screen_abstracts_conflicts", []
        )
    ]
    if len(df_conflicts_data) > 0:
        df_conflicts = pd.DataFrame(df_conflicts_data)
        st.subheader("Conflicts")
        st.dataframe(df_conflicts, use_container_width=True, hide_index=True)


@st.fragment
def render_errors() -> None:
    errors = st.session_state.screen_abstracts_errors
    if len(errors) > 0:
        st.divider()
        st.subheader("Errors")
        for err in errors:
            if is_screening_error(err):
                # st.json(err.model_dump(mode="json"))
                st.write(f"PMID: {err.search_result.pmid}")
                st.error(repr(err.error))
            else:
                st.write(repr(err))


@st.fragment
def render_metrics() -> None:
    col1, col2, col3 = st.columns(3)
    batch_size = st.session_state.get("screen_abstracts_batch_size", 0)
    total = st.session_state.get("screen_abstracts_to_be_screened", 0)
    screened_count = st.session_state.get("screen_abstracts_screened", 0)
    total_tokens = st.session_state.get("screen_abstracts_total_tokens", 0)
    total_cost = st.session_state.get("screen_abstracts_total_cost", 0.0)
    successful_requests = st.session_state.get(
        "screen_abstracts_successful_requests", 0
    )
    error_count = len(st.session_state.get("screen_abstracts_errors", []))
    included = st.session_state.get("screen_abstracts_included", 0)
    excluded = st.session_state.get("screen_abstracts_excluded", 0)
    uncertain = st.session_state.get("screen_abstracts_uncertain", 0)
    conflicts = len(st.session_state.get("screen_abstracts_conflicts", []))

    latency = _calculate_latency(results=st.session_state.screen_abstracts_results)
    with col1.container():
        st.metric("Total abstracts", total)
        st.metric("Pending abstracts", total - screened_count)
        st.metric("Screening errors", error_count)
        # st.metric("Screening latency", f"{latency:.2f} ms")
    with col2.container():
        st.metric("Included abstracts", included)
        st.metric("Excluded abstracts", excluded)
        st.metric("Uncertain abstracts", uncertain)
        st.metric("Conflicts", conflicts)
    with col3.container():
        st.metric("Total tokens", total_tokens)
        st.metric("Total cost", f"${total_cost:.2f}")
        st.metric("Successful API requests", successful_requests)
        st.metric("Screening latency", f"{latency:.2f} s")


def _calculate_latency(results: list[tuple[PubMedResult, ScreeningResult]]) -> float:
    durations = 0
    durations = [(r.end_time - r.start_time).total_seconds() for _, r in results]
    if not durations:
        return 0.0
    return mean(durations)


@st.fragment
def render_progress_bar() -> None:
    batch_size = st.session_state.get("screen_abstracts_batch_size", 0)
    total = st.session_state.get("screen_abstracts_to_be_screened", 0)
    cur_batch_idx = st.session_state.get("screen_abstracts_batch_idx", 0)
    screened_count = st.session_state.get("screen_abstracts_screened", 0)

    with st.container():
        st.progress(screened_count / total)
        if screened_count == total:
            st.success(f"Screening complete! {screened_count} abstracts screened.")
            return
        st.markdown(
            f"Screening abstracts **{cur_batch_idx * batch_size} - {cur_batch_idx * batch_size + batch_size}** of **{total}** abstracts."
        )


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def screen_abstracts(  # noqa: C901
    search_results: list[PubMedResult],
    review: SystematicReview,
) -> None:
    ut.init_state_key("screen_abstracts_results", [])
    ut.init_state_key("screen_abstracts_errors", [])
    ut.init_state_key("screen_abstracts_chain", screen_abstracts_chain)
    ut.init_state_key("screen_abstracts_total_tokens", 0)
    ut.init_state_key("screen_abstracts_successful_requests", 0)
    ut.init_state_key("screen_abstracts_total_cost", 0.0)
    ut.init_state_key("screen_abstracts_included", 0)
    ut.init_state_key("screen_abstracts_excluded", 0)
    ut.init_state_key("screen_abstracts_uncertain", 0)
    ut.init_state_key("screen_abstracts_to_be_screened", len(search_results))
    ut.init_state_key("screen_abstracts_screened", 0)
    ut.init_state_key("screen_abstracts_conflicts", [])

    # status_text = st.status("Screening abstracts ...", expanded=True)

    batch_size = ut.init_state_key("screen_abstracts_batch_size", 10)
    batch_idx = ut.init_state_key("screen_abstracts_batch_idx", -1)
    with st.spinner("Screening abstracts ..."):
        screening_status = st.status("Screening abstracts status ...", expanded=True)
        for i in range(0, len(search_results), batch_size):
            batch = search_results[i : i + batch_size]
            batch_size = len(batch)
            batch_idx += 1
            st.session_state.screen_abstracts_batch_idx = batch_idx
            # runs in executor
            logger.info(f"Screening batch {batch_idx} of size {batch_size} ...")
            screening_status.text(
                f"Screening batch {batch_idx} of size {batch_size} ..."
            )
            batch_output = screen_abstracts_batch(batch, batch_idx, review)
            logger.info(f"Screening batch {batch_idx} of size {batch_size} ... done")
            if not batch_output:
                msg = f"Screening batch {batch_idx} of size {batch_size} returned no results"
                logger.error(msg)
                screening_status.text(msg)
                # st.write(msg)
                st.session_state.screen_abstracts_errors.extend([None] * batch_size)
                return
            cb = batch_output.cb
            results = batch_output.results
            st.session_state.screen_abstracts_total_tokens += cb.total_tokens
            st.session_state.screen_abstracts_total_cost += cb.total_cost
            st.session_state.screen_abstracts_successful_requests += (
                cb.successful_requests
            )
            for j, res in enumerate(results):
                pm_r = res.search_result

                if pm_r.id != batch[j].id:
                    msg = f"Search result id mismatch: {pm_r.id} != {batch[j].id}, skipping ... res: {res!r}"
                    logger.bind(output=pm_r.id, input=batch[j].id).error(msg)
                    # status_text.write(msg)
                    screening_status.text(msg)
                    st.session_state.screen_abstracts_errors.append(
                        ScreeningError(search_result=pm_r, error=res, message=msg)
                    )
                    continue

                if not isinstance(res.conservative_result, ScreeningResult):
                    msg = f"Conservative screening error: {type(res.conservative_result)}: {res.conservative_result!r}"
                    logger.error(msg)
                    # st.write(msg)
                    screening_status.text(msg)
                    st.session_state.screen_abstracts_errors.append(
                        res.conservative_result
                    )
                else:
                    st.session_state.screen_abstracts_results.append(
                        (pm_r, res.conservative_result)
                    )
                    match res.conservative_result.decision:
                        case ScreeningDecisionType.INCLUDE:
                            st.session_state.screen_abstracts_included += 1
                        case ScreeningDecisionType.EXCLUDE:
                            st.session_state.screen_abstracts_excluded += 1
                        case ScreeningDecisionType.UNCERTAIN:
                            st.session_state.screen_abstracts_uncertain += 1

                if not isinstance(res.comprehensive_result, ScreeningResult):
                    msg = f"Comprehensive screening error: {type(res.comprehensive_result)}: {res.comprehensive_result!r}"
                    logger.error(msg)
                    # st.write(msg)
                    screening_status.text(msg)
                    st.session_state.screen_abstracts_errors.append(
                        res.comprehensive_result
                    )
                else:
                    st.session_state.screen_abstracts_results.append(
                        (pm_r, res.comprehensive_result)
                    )
                    match res.comprehensive_result.decision:
                        case ScreeningDecisionType.INCLUDE:
                            st.session_state.screen_abstracts_included += 1
                        case ScreeningDecisionType.EXCLUDE:
                            st.session_state.screen_abstracts_excluded += 1
                        case ScreeningDecisionType.UNCERTAIN:
                            st.session_state.screen_abstracts_uncertain += 1

                st.session_state.screen_abstracts_screened += 1

                if (
                    isinstance(res.conservative_result, ScreeningResult)
                    and isinstance(res.comprehensive_result, ScreeningResult)
                    and res.conservative_result.decision
                    != res.comprehensive_result.decision
                ):
                    st.session_state.screen_abstracts_conflicts.append(
                        (pm_r, res.conservative_result, res.comprehensive_result)
                    )

                # FIXME: these don't update ...
                # render_metrics()
                # render_progress_bar()
                # render_df()
                # render_errors()


def screen_abstracts_page(review_id: uuid.UUID | None = None) -> None:
    """Abstract screening page."""
    # init review_id
    if not review_id and "review_id" not in st.session_state:
        st.error("Please complete and save the review protocol first.")
        st.stop()
    if not review_id:
        review_id = st.session_state.review_id

    review = st.session_state.get("review")
    if not isinstance(review, SystematicReview):
        st.error("Please complete and save the review protocol first.")
        st.stop()

    # Init session state
    # TODO: switch to new repos, use screening repo
    init_pubmed_repository(st.session_state.supabase)
    init_screen_abstracts_chain()

    st.header("Abstract Screening")
    screening_notification_widget = st.empty()
    # Get search results and existing screening decisions
    if "pubmed_results" not in st.session_state:
        st.session_state.pubmed_results = (
            st.session_state.repo_pubmed.get_search_results(review_id)
        )

    if st.button("Screen Abstracts"):
        pubmed_results = st.session_state.pubmed_results
        if not pubmed_results:
            st.error("No PubMed results found, please run a PubMed search first.")
            st.page_link("pages/search.py", label="PubMed Search")
            st.stop()

        # renders metrics, progress bar, dataframe fragments for results, and errors
        screen_abstracts(pubmed_results, review)

        result_count = len(st.session_state.screen_abstracts_results)
        error_count = len(st.session_state.screen_abstracts_errors)
        screening_notification_widget.success(
            f"Screening complete! {result_count} results, {error_count} errors."
        )
        st.session_state.screen_abstracts_done = True

    # Display results table with filtering
    if "screen_abstracts_done" in st.session_state:
        render_metrics()
        render_progress_bar()
        render_df()
        render_errors()


# For quick testing
if "review_id" not in st.session_state:
    st.session_state.review_id = uuid.UUID(hex="393f5e13-175a-4314-b47e-45b804ab3d6f")
if "logger_extra_configured" not in st.session_state:
    logger.configure(extra={"review_id": st.session_state.review_id})
    st.session_state.logger_extra_configured = True
if "review" not in st.session_state:
    with st.session_state.session_factory() as session:
        st.session_state.review = session.get(
            SystematicReview, st.session_state.review_id
        )
init_pubmed_repository(st.session_state.supabase)
if "pubmed_results" not in st.session_state:
    st.session_state.pubmed_results = st.session_state.repo_pubmed.get_search_results(
        st.session_state.review_id
    )
screen_abstracts_page()
