from __future__ import annotations

import typing as t
import uuid
from statistics import mean

import pandas as pd
import streamlit as st
from loguru import logger

import sr_assistant.app.utils as ut
from sr_assistant.app import services
from sr_assistant.app.agents.screening_agents import (
    ScreenAbstractResultTuple,
    ScreeningError,
    ScreeningResult,
    ScreeningStrategyType,
    invoke_resolver_chain,
)
from sr_assistant.app.database import session_factory
from sr_assistant.core.models import (
    SearchResult,
    SystematicReview,
)
from sr_assistant.core.repositories import (
    RecordNotFoundError,
    ScreenAbstractResultRepository,
    ScreeningResolutionRepository,
    SearchResultRepository,
    SystematicReviewRepository,
)
from sr_assistant.core.schemas import ScreeningDecisionType


def init_pubmed_repository() -> SearchResultRepository:
    if "search_repo" not in st.session_state:
        repo = SearchResultRepository()
        st.session_state.search_repo = repo
    return st.session_state.search_repo


def init_review_repository() -> SystematicReviewRepository:
    if "repo_review" not in st.session_state:
        repo = SystematicReviewRepository()
        st.session_state.repo_review = repo
    return st.session_state.repo_review


def init_screen_abstracts_repository() -> ScreenAbstractResultRepository:
    if "repo_screen_abstracts" not in st.session_state:
        repo = ScreenAbstractResultRepository()
        st.session_state.repo_screen_abstracts = repo
    return st.session_state.repo_screen_abstracts


def init_screen_resolution_repository() -> ScreeningResolutionRepository:
    """Initialize ScreeningResolutionRepository in session state."""
    if "repo_screen_resolution" not in st.session_state:
        repo = ScreeningResolutionRepository()
        st.session_state.repo_screen_resolution = repo
    return st.session_state.repo_screen_resolution


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
        "list[tuple[SearchResult, ScreeningResult]]",
        st.session_state.screen_abstracts_results,
    )
    df_data = [
        {
            "PMID": search_result.source_id,
            "Title": search_result.title,
            "Strategy": screening_result.screening_strategy,
            "Decision": screening_result.decision,
            "Confidence": screening_result.confidence_score,
            "Model": screening_result.model_name,
            "Rationale": screening_result.rationale,
            "Extracted quotes": screening_result.extracted_quotes,
            "Result ID": str(screening_result.id),
            "Final Decision": search_result.final_decision,
            "Resolved": bool(search_result.resolution_id),
            "Resolver Reasoning": None,
        }
        for search_result, screening_result in results
    ]
    if len(df_data) > 0:
        st.subheader("Results")
        st.dataframe(df_data, use_container_width=True, hide_index=True)

        if pmid_selection := st.selectbox(
            "View details:",
            [item["PMID"] for item in df_data],
            key="details_select_main_df",
        ):
            selected_sr_object = None
            cons_selected_screening_res = None
            comp_selected_screening_res = None

            for sr, scr_res in st.session_state.screen_abstracts_results:
                if sr.source_id == pmid_selection:
                    if selected_sr_object is None:
                        selected_sr_object = sr
                    if scr_res.screening_strategy == ScreeningStrategyType.CONSERVATIVE:
                        cons_selected_screening_res = scr_res
                    elif (
                        scr_res.screening_strategy
                        == ScreeningStrategyType.COMPREHENSIVE
                    ):
                        comp_selected_screening_res = scr_res

            if selected_sr_object:
                st.subheader(
                    f"{selected_sr_object.source_id}: {selected_sr_object.title}"
                )
                st.text(f"{selected_sr_object.journal} ({selected_sr_object.year})")
                st.markdown(
                    f"[SearchResult](https://pubmed.ncbi.nlm.nih.gov/{selected_sr_object.source_id})"
                )
                st.json(selected_sr_object.model_dump(mode="json"))
                if cons_selected_screening_res:
                    st.markdown(
                        f"[Conservative Screening Result](https://pubmed.ncbi.nlm.nih.gov/{selected_sr_object.source_id})"
                    )
                    st.json(cons_selected_screening_res.model_dump(mode="json"))
                else:
                    st.warning(f"No conservative screening result for {pmid_selection}")
                if comp_selected_screening_res:
                    st.markdown(
                        f"[Comprehensive Screening Result](https://pubmed.ncbi.nlm.nih.gov/{selected_sr_object.source_id})"
                    )
                    st.json(comp_selected_screening_res.model_dump(mode="json"))
                else:
                    st.warning(
                        f"No comprehensive screening result for {pmid_selection}"
                    )
            else:
                st.warning(
                    f"Details for selected Source ID {pmid_selection} not fully found."
                )

    df_conflicts_data = [
        {
            "PMID": pm.source_id,
            "Conservative Decision": conservative_result.decision,
            "Conservative Confidence": conservative_result.confidence_score,
            "Conservative Rationale": conservative_result.rationale,
            "Comprehensive Decision": comprehensive_result.decision,
            "Comprehensive Confidence": comprehensive_result.confidence_score,
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

        if pmid_conflict_selection := st.selectbox(
            "View details",
            [row["PMID"] for row in df_conflicts_data],
            key="conflict_details_select",
        ):
            selected_conflict_tuple = next(
                (
                    (sr, cons_res, comp_res)
                    for sr, cons_res, comp_res in st.session_state.get(
                        "screen_abstracts_conflicts", []
                    )
                    if sr.source_id == pmid_conflict_selection
                ),
                None,
            )

            if selected_conflict_tuple:
                (
                    search_result,
                    conservative_screening_result,
                    comprehensive_screening_result,
                ) = selected_conflict_tuple
                st.subheader(f"{search_result.source_id}: {search_result.title}")
                st.text(f"{search_result.journal} ({search_result.year})")
                st.markdown(
                    f"[SearchResult](https://pubmed.ncbi.nlm.nih.gov/{search_result.source_id})"
                )
                st.json(search_result.model_dump(mode="json"))

                if conservative_screening_result:
                    st.markdown(
                        f"[Conservative Screening Result](https://pubmed.ncbi.nlm.nih.gov/{search_result.source_id})"
                    )
                    st.json(conservative_screening_result.model_dump(mode="json"))
                else:
                    st.warning("No conservative screening result")
                if comprehensive_screening_result:
                    st.markdown(
                        f"[Comprehensive Screening Result](https://pubmed.ncbi.nlm.nih.gov/{search_result.source_id})"
                    )
                    st.json(comprehensive_screening_result.model_dump(mode="json"))
                else:
                    st.warning("No comprehensive screening result")
            else:
                st.warning(
                    f"Conflict details for selected Source ID {pmid_conflict_selection} not found."
                )


@st.fragment
def render_errors() -> None:
    errors = st.session_state.screen_abstracts_errors
    if len(errors) > 0:
        st.divider()
        st.subheader("Errors")
        for err in errors:
            if is_screening_error(err):
                # st.json(err.model_dump(mode="json"))
                st.write(f"PMID: {err.search_result.source_id}")
                st.error(repr(err.error))
            else:
                st.write(repr(err))


@st.fragment
def render_metrics() -> None:
    col1, col2, col3 = st.columns(3)
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
    with col2.container():
        st.metric("Included abstracts", included)
        st.metric("Excluded abstracts", excluded)
        st.metric("Uncertain abstracts", uncertain)
        st.metric("Conflicts", conflicts)
    with col3.container():
        st.metric("Total tokens", total_tokens)
        st.metric("Total cost", f"${total_cost:.2f}")
        st.metric("Successful API requests", successful_requests)
        st.metric("Mean Screening latency", f"{latency:.2f} s")


def _calculate_latency(results: list[tuple[SearchResult, ScreeningResult]]) -> float:
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


def init_screening_service() -> services.ScreeningService:
    if "screening_service" not in st.session_state:
        # Repositories are already initialized and available in session_state by this point usually
        # or can be initialized here if preferred.
        # For simplicity, assuming they are available or ScreeningService defaults work.
        st.session_state.screening_service = services.ScreeningService(
            # factory=session_factory, # Uses default
            # screen_repo=st.session_state.repo_screen_abstracts, # Can be passed if already init
            # search_repo=st.session_state.search_repo,
            # review_repo=st.session_state.repo_review
        )
    return st.session_state.screening_service


@logger.catch(onerror=lambda exc: st.error(exc) if ut.in_streamlit() else None)  # pyright: ignore [reportArgumentType]
def screen_abstracts(  # noqa: C901
    search_results_to_process: list[SearchResult],
    review: SystematicReview,
    screening_service: services.ScreeningService,
) -> None:
    ut.init_state_key("screen_abstracts_results", [])
    ut.init_state_key("screen_abstracts_errors", [])
    # ut.init_state_key("screen_abstracts_chain", screen_abstracts_chain) # No longer directly used here
    # Removed token/cost tracking that relied on direct agent cb
    # ut.init_state_key("screen_abstracts_total_tokens", 0)
    # ut.init_state_key("screen_abstracts_successful_requests", 0)
    # ut.init_state_key("screen_abstracts_total_cost", 0.0)
    ut.init_state_key("screen_abstracts_included", 0)
    ut.init_state_key("screen_abstracts_excluded", 0)
    ut.init_state_key("screen_abstracts_uncertain", 0)
    ut.init_state_key("screen_abstracts_screened", 0)
    ut.init_state_key("screen_abstracts_conflicts", [])
    ut.init_state_key("screen_abstracts_resolved", [])

    # batch_size = st.session_state.screen_abstracts_batch_size # Service handles batching internally or not at all for this method
    # batch_idx = st.session_state.screen_abstracts_batch_idx = -1
    # total_batches = (
    #     st.session_state.get("screen_abstracts_to_be_screened", 0) // batch_size
    # )
    total_to_screen = len(search_results_to_process)
    st.session_state.screen_abstracts_to_be_screened = (
        total_to_screen  # Update this for progress bar
    )

    screening_status = st.status(
        f"Screening {total_to_screen} abstracts...",
        expanded=True,
    )

    with screening_status:
        search_result_ids = [sr.id for sr in search_results_to_process]
        if not search_result_ids:
            st.error("No valid search result IDs to screen.")
            screening_status.update(label="Screening failed: No IDs.", state="error")
            return

        logger.info(
            f"Calling ScreeningService for {len(search_result_ids)} search results."
        )
        screening_status.text(
            f"Invoking screening service for {len(search_result_ids)} abstracts..."
        )

        try:
            # Call the service
            # The service handles fetching fresh models by ID, persistence, and transactions.
            service_results: list[ScreenAbstractResultTuple] = (
                screening_service.perform_batch_abstract_screening(
                    review_id=review.id, search_result_ids_to_screen=search_result_ids
                )
            )
            logger.info(
                f"ScreeningService returned {len(service_results)} result tuples."
            )

            # Reset counters before processing, as service might return partials or all
            st.session_state.screen_abstracts_included = 0
            st.session_state.screen_abstracts_excluded = 0
            st.session_state.screen_abstracts_uncertain = 0
            st.session_state.screen_abstracts_screened = (
                0  # Will be incremented per search result processed
            )
            st.session_state.screen_abstracts_results = []  # Clear previous results from UI perspective
            st.session_state.screen_abstracts_errors = []
            st.session_state.screen_abstracts_conflicts = []

            # Process results returned by the service
            for res_tuple in service_results:
                # The search_result in res_tuple is the one processed by the agent.
                # It should have conservative_result_id and comprehensive_result_id set if screening was successful.
                pm_r = res_tuple.search_result

                # Note: The service ensures that ScreenAbstractResult records are persisted
                # and SearchResult records are updated with linkage IDs (conservative_result_id, comprehensive_result_id).
                # The UI's main job here is to display these results and count metrics.

                if not isinstance(res_tuple.conservative_result, ScreeningResult):
                    # This is a ScreeningError for conservative path
                    msg = f"Conservative screening error for {pm_r.source_id}: {res_tuple.conservative_result!r}"
                    logger.error(msg)
                    st.session_state.screen_abstracts_errors.append(
                        res_tuple.conservative_result
                    )
                else:
                    st.session_state.screen_abstracts_results.append(
                        (pm_r, res_tuple.conservative_result)
                    )
                    match res_tuple.conservative_result.decision:
                        case ScreeningDecisionType.INCLUDE:
                            st.session_state.screen_abstracts_included += 1
                        case ScreeningDecisionType.EXCLUDE:
                            st.session_state.screen_abstracts_excluded += 1
                        case ScreeningDecisionType.UNCERTAIN:
                            st.session_state.screen_abstracts_uncertain += 1

                if not isinstance(res_tuple.comprehensive_result, ScreeningResult):
                    # This is a ScreeningError for comprehensive path
                    msg = f"Comprehensive screening error for {pm_r.source_id}: {res_tuple.comprehensive_result!r}"
                    logger.error(msg)
                    st.session_state.screen_abstracts_errors.append(
                        res_tuple.comprehensive_result
                    )
                else:
                    st.session_state.screen_abstracts_results.append(
                        (pm_r, res_tuple.comprehensive_result)
                    )
                    # Avoid double counting for overall included/excluded/uncertain if one strategy errored
                    # but we still want to capture the successful strategy's decision for metrics.
                    # This logic might need refinement based on how metrics should reflect partial success.
                    # For simplicity, let's assume if one strategy is successful, its decision contributes to metrics.
                    # A search result is considered "screened" if at least one strategy produced a result (error or success).
                    if not isinstance(
                        res_tuple.conservative_result, ScreeningResult
                    ):  # only add if conservative errored
                        match res_tuple.comprehensive_result.decision:
                            case ScreeningDecisionType.INCLUDE:
                                st.session_state.screen_abstracts_included += 1
                            case ScreeningDecisionType.EXCLUDE:
                                st.session_state.screen_abstracts_excluded += 1
                            case ScreeningDecisionType.UNCERTAIN:
                                st.session_state.screen_abstracts_uncertain += 1

                st.session_state.screen_abstracts_screened += (
                    1  # Count as screened if we got a tuple for it
                )

                # Conflict detection
                if (
                    isinstance(res_tuple.conservative_result, ScreeningResult)
                    and isinstance(res_tuple.comprehensive_result, ScreeningResult)
                    and res_tuple.conservative_result.decision
                    != res_tuple.comprehensive_result.decision
                ):
                    st.session_state.screen_abstracts_conflicts.append(
                        (
                            pm_r,
                            res_tuple.conservative_result,
                            res_tuple.comprehensive_result,
                        )
                    )

        except RecordNotFoundError as e:
            logger.error(f"Service error (RecordNotFound): {e}")
            st.error(f"Error during screening: {e}")
            screening_status.update(label=f"Screening failed: {e}", state="error")
            return  # Stop processing
        except services.ServiceError as e:
            logger.error(f"Service error: {e}")
            st.error(f"An unexpected error occurred during screening: {e}")
            screening_status.update(label=f"Screening failed: {e}", state="error")
            return
        except Exception as e:  # Catch any other unexpected error from the service call
            logger.error(
                f"Unexpected error calling screening service: {e}", exc_info=True
            )
            st.error(f"A critical error occurred: {e}")
            screening_status.update(
                label=f"Screening failed critically: {e}", state="error"
            )
            return

    screening_status.update(
        label=f"Screening of {st.session_state.screen_abstracts_screened}/{total_to_screen} abstracts processed by service.",
        state="complete",
        expanded=False,
    )

    # --- Resolve Conflicts --- # (This logic can largely remain, as it depends on st.session_state.screen_abstracts_conflicts)
    conflicts_to_resolve = st.session_state.screen_abstracts_conflicts
    if conflicts_to_resolve:
        resolution_status = st.status(
            f"Resolving {len(conflicts_to_resolve)} conflicts...", expanded=True
        )
        with resolution_status:
            resolved_count = 0
            for pm_result, cons_res, comp_res in conflicts_to_resolve:
                try:
                    st.write(f"Resolving conflict for PMID: {pm_result.source_id}")
                    resolution_output = invoke_resolver_chain(
                        search_result=pm_result,
                        review=review,
                        conservative_result=cons_res,
                        comprehensive_result=comp_res,
                    )
                    if resolution_output:
                        logger.info(
                            f"Simulated resolution for PMID {pm_result.source_id}: {resolution_output.resolver_decision.value}"
                        )
                        st.session_state.screen_abstracts_resolved.append(
                            resolution_output
                        )
                        resolved_count += 1
                        st.write(
                            f"Resolved PMID {pm_result.source_id}: {resolution_output.resolver_decision.value}"
                        )
                    else:
                        st.error(
                            f"Failed to get resolution for PMID {pm_result.source_id}. Check logs."
                        )

                except Exception as e:
                    msg = (
                        f"Error resolving conflict for PMID {pm_result.source_id}: {e}"
                    )
                    logger.exception(msg)
                    st.error(msg)

            resolution_status.update(
                label=f"Conflict resolution finished. Resolved {resolved_count}/{len(conflicts_to_resolve)}.",
                state="complete",
                expanded=False,
            )
            # Update conflicts metric after resolution attempt
            st.session_state.screen_abstracts_conflicts_remaining = (
                len(conflicts_to_resolve) - resolved_count
            )


def screen_abstracts_page(review_id: uuid.UUID | None = None) -> None:
    """Abstract screening page."""
    st.title("Abstract Screening")
    # init review_id
    if not review_id and "review_id" not in st.session_state:
        st.error("Please complete and save the review protocol first.")
        st.stop()
    if not review_id:
        review_id = st.session_state.review_id

    # Init session state
    init_review_repository()
    init_pubmed_repository()
    init_screen_abstracts_repository()
    init_screen_resolution_repository()
    screening_service = init_screening_service()  # Initialize ScreeningService

    # init review
    review = st.session_state.get("review")
    # Ensure review is fetched if not in session state or not the correct type
    if not isinstance(review, SystematicReview):
        review_repo = st.session_state.repo_review  # Get the initialized repo
        with session_factory() as db_session:  # Obtain a session
            # Call get_by_id with both session and id arguments
            review = review_repo.get_by_id(session=db_session, id=review_id)
            st.session_state.review = review  # Store the fetched review

    # Add a guard in case the review could not be fetched
    if not review:
        st.error(
            f"Systematic Review with ID {review_id} not found. Please ensure it exists."
        )
        st.stop()

    screening_notification_widget = st.empty()

    # Always fetch fresh SearchResult *model* instances for screening from the database.
    # This ensures we have the SQLModel instances that can be modified by screen_abstracts_batch,
    # overriding any SearchResultRead instances that might be in session_state from the search page.
    # The review_id used here is confirmed to be non-None by guards at the start of screen_abstracts_page.
    confirmed_review_id = (
        st.session_state.review_id
    )  # This is guaranteed to be a UUID by prior checks
    logger.info(
        f"Screening page: Fetching SearchResult models for review_id: {confirmed_review_id}"
    )
    search_repo = init_pubmed_repository()
    with session_factory() as db_session:
        # Pass the session to the repository method
        # And use the confirmed_review_id which is known to be a UUID
        retrieved_search_results = search_repo.get_by_review_id(
            session=db_session, review_id=confirmed_review_id
        )
        st.session_state.search_results = list(
            retrieved_search_results
        )  # Convert Sequence to list

    if not st.session_state.search_results:
        logger.warning(
            f"No search results found in DB for review_id: {confirmed_review_id} on screening page."
        )

    st.subheader("Screening review")
    st.json(review.model_dump(mode="json"), expanded=False)

    ut.init_state_key(
        "screen_abstracts_to_be_screened", len(st.session_state.search_results)
    )
    col1, col2, col3 = st.columns(3)
    col1.write(f"Abstracts to be screened: {len(st.session_state.search_results)}")
    st.session_state.screen_abstracts_batch_size = col2.slider(
        "Batch size", min_value=1, max_value=40, value=10
    )
    if col3.button("Start Abstract Screening"):
        st.session_state.screen_abstracts_screened = 0
        # st.session_state.screen_abstracts_batch_idx = 0 # Not needed
        st.session_state.screen_abstracts_done = False
        st.session_state.screen_abstracts_results = []
        st.session_state.screen_abstracts_errors = []
        st.session_state.screen_abstracts_conflicts = []
        st.session_state.screen_abstracts_resolved = []
        st.session_state.screen_abstracts_included = 0
        st.session_state.screen_abstracts_excluded = 0
        st.session_state.screen_abstracts_uncertain = 0

        search_results_for_screening = list(st.session_state.search_results)
        if not search_results_for_screening:
            st.error(
                "No PubMed results found in the database for this review. Please run a PubMed search first."
            )
            st.stop()

        # Call the refactored screen_abstracts function, passing the service
        screen_abstracts(
            search_results_for_screening, review, screening_service
        )  # Pass service

        result_count = len(
            st.session_state.screen_abstracts_results
        )  # This is now 2x num_articles if both strats succeed
        actual_articles_processed = st.session_state.screen_abstracts_screened
        conflict_count = len(st.session_state.screen_abstracts_conflicts)
        error_count = len(st.session_state.screen_abstracts_errors)
        screening_notification_widget.success(
            f"Screening processing complete! {actual_articles_processed} articles processed. {result_count} individual screening results. {conflict_count} conflicts. {error_count} errors."
        )
        st.session_state.screen_abstracts_done = True

    # Display results table with filtering
    if (
        "screen_abstracts_done" in st.session_state
        and st.session_state.screen_abstracts_done
    ):
        render_metrics()
        # render_progress_bar() # FIXME
        render_df()
        render_errors()


# Ensure the main page logic is called when the script is run by AppTest or Streamlit.
# The test fixture app_test_env will have set up st.session_state.review_id and st.session_state.review.
# For direct streamlit run, review_id might need to be handled (e.g. from query params or a default).
screen_abstracts_page(st.session_state.get("review_id"))
