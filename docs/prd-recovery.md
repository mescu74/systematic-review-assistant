# Product Requirements Document (PRD) for MPH SR Prototype Recovery Project

## Status: Draft

## Intro

The Systematic Review Assistant (MPH SR Prototype) application is currently in a significantly unstable and partially dysfunctional state due to incomplete, concurrent development efforts: a search functionality refactoring with a new service layer, and the implementation of an automated screening conflict resolver. These unfinished initiatives have resulted in numerous errors, data model inconsistencies, and a non-functional screening workflow. This document outlines the requirements for a recovery project to stabilize the application, complete the refactoring and resolver feature, and establish a solid foundation for future development, as detailed in `docs/SRA_Recovery_Project_Brief.md`.

## Goals and Context

- **Project Objectives:**
    - Restore the MPH SR Prototype to a stable, robust, and reliable state.
    - Correctly implement a generic search infrastructure, initially focused on PubMed.
    - Deliver a fully functional automated screening conflict resolver.
    - Establish a well-defined service layer architecture.
    - Implement comprehensive test coverage.
- **Measurable Outcomes:**
    - Elimination of all current linter errors in critical modules (`search.py`, `screen_abstracts.py`, `services.py`, `repositories.py`, and related tests).
    - Successful and reliable PubMed search functionality via `SearchService`.
    - Accurate display and processing of `SearchResult` data throughout the screening workflow.
    - Automated conflict resolver accurately processes disagreements, with results correctly stored and displayed.
    - Unit test coverage for services, repositories, and models meets defined targets (e.g., >75%).
    - All specified integration tests, including an end-to-end PubMed search-screen-resolve flow, pass reliably.
- **Success Criteria:**
    - The application is stable and key workflows (PubMed search, abstract screening, conflict resolution) are fully operational without errors.
    - The implemented service layer and data models are consistent and correctly utilized.
    - The automated conflict resolver functions as specified in `docs/ai/prd.resolver.md` (with necessary adaptations).
    - Test suites provide confidence in the stability and correctness of the application.
    - Integration testing can be performed safely and reliably against the designated test database.
- **Key Performance Indicators (KPIs):**
    - Number of linter errors in key modules (target: 0).
    - Pass rate for all unit and integration tests (target: 100%).
    - Unit test coverage percentage for `services.py`, `repositories.py`, `models.py` (target: >80%).
    - Successful completion of the end-to-end PubMed search -> screen -> resolve workflow without manual intervention or errors.

## Scope and Requirements (MVP / Current Version)

### Functional Requirements (High-Level)

- **FR1: Stabilized Search & Service Layer (PubMed Focus):**
    - `search.py` UI correctly interacts with `SearchService` for PubMed searches.
    - `SearchService` correctly manages PubMed search logic, DB sessions, and calls appropriate `SearchResultRepository` methods.
    - `SearchResultRepository` is equipped with all necessary methods (e.g., for storing, retrieving, deleting search results by review ID for PubMed).
    - PubMed search results are consistently stored and retrieved using the `SearchResult` model, with correct `source_db` ('PubMed') and `source_id` (PMID).
    - Resolve all linter errors in `search.py`, `services.py`, `repositories.py` and related tests.
    - PubMed search functionality is fully operational.
- **FR2: Completed Automated Conflict Resolver Implementation:**
    - Implement the missing `resolve_screening_conflict` function.
    - Update `SearchResult` model to include `final_decision` (with DB migration).
    - Ensure correct data flow for the resolver LLM and storage of its decisions in `ScreeningResolution` and `SearchResult`.
    - Integrate fully into `screen_abstracts.py` UI and workflow.
    - *(Note: The scope of disagreements handled by the resolver now includes cases involving UNCERTAIN decisions, as per the updated `docs/prd-resolver.md` v1.1).*
- **FR3: Comprehensive and Safe Test Suite:**
    - Develop unit tests for all services in `services.py` (e.g., `SearchService`, `ScreeningService`).
    - Repair and significantly expand unit tests for all data models (`models.py`) and repositories (`repositories.py`), reflecting service-managed session patterns.
    - Complete and ensure all integration tests for the resolver feature pass (e.g., `tests/integration/test_resolver.py`, `tests/integration/test_resolver_agent.py`).
    - Create and pass at least one end-to-end integration test covering the PubMed search -> abstract screening -> conflict resolution workflow.
    - Document and enforce safe integration testing procedures against the `sra_integration_test` Supabase database, including environment variable management (`.env.test`) and awareness of schema cleaning in `conftest.py`.

### Non-Functional Requirements (NFRs)

- **Stability:** The application must be free of critical runtime errors in the recovered workflows.
- **Reliability:** Core features (search, screen, resolve) must function consistently and predictably.
- **Data Integrity:** Data models (`SearchResult`, `ScreeningResolution`) must be consistent, and database interactions must be correct. `source_id` and `source_db` must be used correctly for PubMed data (i.e. `source_db` = 'PubMed'). The `final_decision` field must be accurately populated.
- **Maintainability:** Code should adhere to defined quality standards, be well-documented where necessary, and the service layer architecture must be clearly implemented.
    - Non-critical linter errors (e.g., the known issue with `Model.id` in `repositories.py`) may be temporarily addressed with `noqa` or `pyright: ignore` comments to focus on application-breaking issues. Critical linter errors that impede functionality must be fixed.
- **Testability:** The system must be structured to support comprehensive unit and integration testing.
- **Performance:** While not a primary focus for deep optimization in recovery, the system should perform searches and screening actions within a reasonable timeframe for user interaction. Linter errors that might impact performance should be resolved.

### User Experience (UX) Requirements (High-Level)

- The existing UI in `search.py` and `screen_abstracts.py` should correctly display data from the `SearchResult` model (using `source_id`, `source_db` where `source_db` is 'PubMed' for PubMed results).
- The UI in `screen_abstracts.py` should clearly indicate the status of screening decisions, including resolved conflicts and the resolver's reasoning.

### Integration Requirements (High-Level)

- Integration with OpenAI API for the resolver LLM agent.
- Integration with Supabase (PostgreSQL) for data persistence.
- (Implicit) Integration with PubMed API via existing mechanisms, managed by `SearchService`.

### Testing Requirements (High-Level)

- Comprehensive unit tests for services, repositories, and models.
- Robust integration tests for the conflict resolver functionality.
- At least one end-to-end integration test for the primary PubMed search-screen-resolve workflow.
- Clear documentation and procedures for safe execution of integration tests against a dedicated test database.
- *(See `docs/SRA_Recovery_Project_Brief.md` Goal 3 for detailed testing specifics and `docs/templates/testing-strategy.md` for general strategy if/when created).*

## Tech Stack

- **Languages:** Python
- **Frameworks/Libraries:** Streamlit, SQLModel, SQLAlchemy, LangChain, LangGraph, OpenAI SDK
- **Database:** PostgreSQL (via Supabase)
- **Testing:** Pytest
- *(This reflects the existing stack to be stabilized. Further details may be in `docs/tech-stack.md` if it exists or is created).*

## Epic Overview (MVP / Current Version)

- **Epic-1: Search and Service Layer Stabilization (Current)** - Goal: Fully refactor and stabilize `search.py`, `services.py`, and `repositories.py` to correctly implement and utilize the `SearchResult` model and the new service layer architecture for PubMed search operations, ensuring all related functionalities are operational and error-free.
- **Epic-2: Resolver Agent Implementation and Integration (Future)** - Goal: Complete the implementation and integration of the automated screening conflict resolver, including model updates, backend logic, and UI integration, ensuring it accurately processes conflicts and stores results.
- **Epic-3: Comprehensive Testing and Database Integrity (Future)** - Goal: Establish a robust testing framework (unit, integration, E2E for the core recovery flow) and ensure database integrity, including clear documentation for safe integration testing protocols.

## Key Reference Documents

- `docs/SRA_Recovery_Project_Brief.md`
- `docs/SRA_Recovery_PM_Prompt.md`
- `docs/ai/prd.resolver.md` (for resolver logic guidance)
- `src/sr_assistant/core/models.py`
- `src/sr_assistant/core/repositories.py`
- `src/sr_assistant/app/services.py`
- `src/sr_assistant/app/pages/search.py`
- `src/sr_assistant/app/pages/screen_abstracts.py`
- `src/sr_assistant/app/agents/screening_agents.py`
- `tests/conftest.py`
- Existing integration tests (e.g., `tests/integration/core/test_repositories.py`, `tests/integration/test_resolver.py`)

## Post-MVP / Future Enhancements

- Implementation of search functionality for Embase and Scopus.
- Development of the new Benchmarking Module.
- Enhancements to XAI and logging/tracing dashboards.
- *(These are beyond the scope of this recovery PRD and will be addressed in `docs/MPH_SR_Prototype_Main_Project_Brief.md` and subsequent PRDs).*

## Change Log

| Change          | Date       | Version | Description                                                                                                | Author               |
|-----------------|------------|---------|------------------------------------------------------------------------------------------------------------|----------------------|
| Scope Update    | 2025-05-09 | 0.2     | Clarified resolver scope in FR2 to align with updated `prd-resolver.md` (handling of UNCERTAIN cases). | Product Manager Agent |
| Initial Draft   | 2025-05-09 | 0.1     | First draft of PRD                                                                                         | Product Manager Agent |

## Initial Architect Prompt

### Technical Infrastructure

- **Starter Project/Template:** The existing "MPH SR Prototype" codebase.
- **Hosting/Cloud Provider:** Supabase (for PostgreSQL database) is the primary hosted component. The Streamlit application runs locally during this recovery phase.
- **Frontend Platform:** Streamlit (existing).
- **Backend Platform:** Python with Streamlit. FastAPI is a post-prototype consideration and not part of the current recovery effort. The application is primarily synchronous; existing asynchronous database code is not currently in active use.
- **Database Requirements:** PostgreSQL (via Supabase), using SQLModel and SQLAlchemy.

### Technical Constraints

- Adhere to the intended service layer pattern: services manage sessions and business logic; repositories handle data I/O.
- Resolver implementation should align with `docs/ai/prd.resolver.md` where feasible.
- Continue using the remote Supabase instance (`sra_integration_test`) for integration tests, with strict adherence to safety protocols defined in `tests/conftest.py` and `.env.test`.
- Ensure all LLM interactions (especially for the resolver) are robust and handle potential API errors.

### Deployment Considerations

- The immediate goal is a stable application for local development and testing. CI/CD for deployment is not in scope for this recovery project.
- Deployment to Streamlit Community Cloud is planned for a future phase, after the recovery is complete and the initial benchmarking module is developed.
- Existing GitHub Actions for running tests are in place and considered sufficient for the recovery phase.
- Ensure `.env.test` is correctly configured and used for integration tests to target `sra_integration_test` DB.

### Local Development & Testing Requirements

- Developers must be able to run the Streamlit application locally.
- Pytest should be used for running unit and integration tests.
- Clear instructions for setting up the local environment, including environment variables for development (`.env`, `.env.local`) and testing (`.env.test`).
- The schema cleaning mechanism in `tests/conftest.py` for the integration test database must be understood and maintained.
- Pre-commit hooks are part of the project but are temporarily disabled; re-enabling and refining them is a post-recovery concern.

### Other Technical Considerations

- **Data Migration:** Alembic is used for database migrations. Migrations will be generated and run by the user as needed (e.g., for adding the `final_decision` field).
- **Error Handling:** Robust error handling should be implemented, especially around service calls, database interactions, and LLM API calls.
- **Logging/Tracing:** While LangSmith is mentioned for the full vision, ensure sufficient logging within the recovery scope to debug issues.
