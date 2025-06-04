# Project Brief: MPH SR Prototype Stabilization, Refactoring Completion, and Resolver Implementation

**Version:** 1.0
**Date:** 2024-07-29
**Prepared by:** AIDE (AI Assistant) & User

## 1. Introduction / Problem Statement

The Systematic Review Assistant (MPH SR Prototype) application is currently in a significantly unstable and partially dysfunctional state. This is the result of two concurrent, large-scale development efforts undertaken previously, which were left incomplete:

* **Effort 1: Search Functionality Refactoring & Service Layer Introduction:**
    * **Objective:** To transition from a PubMed-specific search data model (`PubMedResult`) to a generic `SearchResult` model, enabling future integration of other databases like Scopus and Embase. This involved creating a new `SearchService` layer in `services.py` intended to manage database sessions and business logic, with repositories in `repositories.py` becoming simpler I/O layers.
    * **Current State:** This refactoring is incomplete.
        * The `SearchResult` model is defined, but its consistent use (e.g., `source_id` instead of `pmid`) is not enforced across the application, particularly in UI components like `search.py` and `screen_abstracts.py`.
        * The `SearchService` layer has been introduced, but its integration with `search.py` and `SearchResultRepository` is flawed, leading to missing methods in the repository (e.g., `store_results`, `delete_by_review_id` are called by `search.py` but not defined in `SearchResultRepository`).
        * Linter errors and likely runtime failures plague `search.py`, `services.py`, and `repositories.py`.

* **Effort 2: Automated Screening Conflict Resolver Implementation:**
    * **Objective:** To introduce a third LLM-driven "Resolver Agent" to automatically adjudicate disagreements between the existing "conservative" and "comprehensive" screening agents, thereby improving efficiency.
    * **Current State:** This feature is also partially implemented.
        * A PRD for the resolver exists (`docs/ai/prd.resolver.md`).
        * Backend components like `resolver_model` in `screening_agents.py` and the `ScreeningResolution` database model in `models.py` are present.
        * However, the core `resolve_screening_conflict` function, intended to orchestrate the resolution process and integrate with `screen_abstracts.py`, is missing.
        * The `SearchResult` database model is misaligned with the resolver PRD, notably missing the required `final_decision` field.
        * `screen_abstracts.py` contains UI stubs and repository calls for the resolver that are currently non-functional due to missing backend logic and incorrect model attributes.

**Consequences:**
The combination of these incomplete efforts has led to:

* Numerous linter errors across critical files.
* Likely widespread runtime failures, rendering key workflows like search and abstract screening unusable.
* Data model inconsistencies and incorrect database interactions.
* A fragile codebase where the intended architecture (service layer controlling sessions, repositories as simple I/O) is not correctly realized.
* Low and partially broken test coverage (unit and integration) which failed to prevent or highlight these issues during development.

This brief outlines a recovery project to address these critical issues, stabilize the application, complete the refactoring and the resolver feature, and establish a solid foundation for future development of the MPH SR Prototype.

## 2. Vision & Goals (for the recovery project)

* **Vision:** To restore the MPH SR Prototype to a stable, robust, and reliable application. This includes a correctly implemented generic search infrastructure (initially focused on PubMed), a fully functional automated screening conflict resolver, a well-defined service layer, and comprehensive test coverage, enabling confident future development, including the planned benchmarking module for the MPH SR Prototype.
* **Primary Goals (SMART):**
    1. **Goal 1 (Search & Service Layer Stabilization - PubMed Focus):**
        * **Action:** Fully refactor and stabilize `search.py`, `services.py`, and `repositories.py` to correctly implement and utilize the `SearchResult` model and the new service layer architecture for PubMed search operations.
        * **Specifics:** Resolve all linter errors; ensure `search.py` correctly uses `SearchService`; `SearchService` correctly manages PubMed search logic and DB sessions; `SearchResultRepository` has the necessary methods (e.g., for adding and deleting search results by review) called by the service. Ensure `source_id` and `source_db` are used consistently for PubMed data.
        * **Measureable:** All PubMed search functionalities in `search.py` are operational. Zero linter errors in related files. Unit tests for `SearchService` (PubMed path) pass. Integration tests for PubMed search pass.
        * **Achievable:** By focusing on PubMed first, the scope is manageable.
        * **Relevant:** Critical for basic application functionality.
        * **Time-bound:** Within 2 sprints.
    2. **Goal 2 (Resolver Agent Implementation & Integration):**
        * **Action:** Complete the implementation and integration of the automated screening conflict resolver.
        * **Specifics:** Develop and integrate the missing `resolve_screening_conflict` function. Align the `SearchResult` model with `prd.resolver.md` by adding the `final_decision` field (including DB migration). Ensure correct data flow to/from the resolver LLM. Store resolution outcomes in `ScreeningResolution` and link to `SearchResult`. Integrate resolver logic seamlessly into the `screen_abstracts.py` workflow and UI.
        * **Measureable:** Resolver automatically processes conflicts. `final_decision` is populated. UI displays resolver output. Integration tests for resolver pass.
        * **Achievable:** Leverages existing PRD and partial backend code.
        * **Relevant:** Addresses a key efficiency improvement identified in `prd.resolver.md`.
        * **Time-bound:** Within 3 sprints (can overlap with Goal 1).
    3. **Goal 3 (Comprehensive Testing & Database Integrity):**
        * **Action:** Establish a robust testing framework and ensure database integrity.
        * **Specifics:**
            * Write dedicated unit tests for `SearchService`, `ScreeningService` (and any other services in `services.py`).
            * Repair and significantly expand unit tests for all data models (`models.py`) and repositories (`repositories.py`), ensuring they reflect the service-managed session pattern and address existing linter errors.
            * Ensure `tests/integration/core/test_repositories.py` accurately tests repository logic under service-controlled sessions if services are used in those tests, or ensure they test repositories in isolation with manually managed sessions if that\'s the intent.
            * Complete and ensure all integration tests for the resolver feature (`tests/integration/test_resolver.py`, `tests/integration/test_resolver_agent.py`) pass.
            * Develop and pass at least one end-to-end integration test covering the PubMed search -> abstract screening -> conflict resolution workflow.
            * Thoroughly document and enforce the procedure for running integration tests safely against the designated `sra_integration_test` Supabase database, including environment variable management (`.env.test`) and awareness of the schema cleaning in `conftest.py`.
        * **Measureable:** Quantifiable increase in unit test coverage (target >75%) for critical modules. All specified integration tests pass consistently. Documentation for integration testing is clear and accessible.
        * **Achievable:** Builds on existing test structures.
        * **Relevant:** Essential for long-term stability and maintainability.
        * **Time-bound:** Ongoing throughout the project, with a final review and documentation sprint (target +1 sprint post Goal 1 & 2 completion).
* **Success Metrics (Initial Ideas):**
    * Elimination of all current linter errors in `search.py`, `screen_abstracts.py`, `services.py`, `repositories.py`, and their respective unit tests.
    * `search.py` successfully performs PubMed searches via `SearchService`, with results correctly stored as `SearchResult` instances and displayed accurately in the UI. (Focus on PubMed initially)
    * `screen_abstracts.py` correctly loads and displays `SearchResult` data (using `source_id`, `source_db`), facilitates dual-reviewer screening, and seamlessly triggers the resolver for conflicts.
    * The resolver agent accurately processes disagreements, with `final_decision` populated in `SearchResult` and `ScreeningResolution` records correctly created and linked. UI reflects resolved states and reasoning.
    * Unit test coverage for `services.py`, `repositories.py`, and `models.py` meets defined targets (e.g., >75%).
    * All integration tests, including the end-to-end PubMed search-screen-resolve flow, pass reliably.
    * No critical issues reported from integration testing related to incorrect database targeting.

## 3. Target Audience / Users

* Primary users remain **Researchers** conducting systematic reviews. The recovery project aims to restore and enhance their ability to efficiently find and screen literature.

## 4. Key Features / Scope (High-Level for Recovery MVP)

* **4.1. Stabilized Search Infrastructure (PubMed Focus):**
    * **Core Logic:** `search.py` UI correctly interacts with `SearchService` for PubMed searches. `SearchService` manages DB sessions and calls appropriate `SearchResultRepository` methods.
    * **Repository Methods:** `SearchResultRepository` equipped with necessary methods (e.g., for adding all results for a review, deleting results by review ID, getting results by review ID) used by the `SearchService`.
    * **Data Integrity:** PubMed search results are consistently stored and retrieved as `SearchResult` objects with correct `source_db` and `source_id`.
    * **UI:** `search.py` correctly displays PubMed data, including links and identifiers. Functionality for Scopus and Embase search will be deferred to future epics.
* **4.2. Functional Abstract Screening Workflow:**
    * **Core Logic:** `screen_abstracts.py` correctly fetches and processes `SearchResult` data for dual-reviewer (conservative/comprehensive) screening.
    * **Stability:** Existing screening chain logic (`screening_agents.py`) is stable and correctly interacts with updated models.
* **4.3. Operational Automated Conflict Resolver:**
    * **Core Logic:** Implementation of the `resolve_screening_conflict` function within `screen_abstracts.py` or `screening_agents.py`.
    * **Model Alignment:** `SearchResult` model in `models.py` updated to include `final_decision: ScreeningDecisionType | None` (with corresponding DB migration).
    * **Data Flow:** Correct input preparation for the resolver LLM, parsing of its output, and storage into `ScreeningResolution` and `SearchResult` tables.
    * **UI Integration:** `screen_abstracts.py` displays the `final_decision`, indicates resolved status, and provides access to resolver reasoning.
* **4.4. Robust and Reliable Test Suite:**
    * **Unit Tests:** Dedicated tests for `SearchService` and other services. Fixed and expanded tests for `models.py` and `repositories.py`.
    * **Integration Tests:** Repaired and expanded tests for resolver functionality. New end-to-end tests for the primary user workflow (PubMed search-screen-resolve).
* **4.5. Secure and Documented Integration Testing Protocol:**
    * Clear, actionable documentation on setting up and running integration tests against the `sra_integration_test` database, emphasizing environment variable management (`.env.test`) and the automated schema cleaning process in `conftest.py`. Highlighting the criticality of not running these tests against production.

## 5. Known Technical Constraints or Preferences

* **Technology Stack:** Continue with Python, Streamlit, SQLModel, SQLAlchemy, LangChain, OpenAI, and PostgreSQL (Supabase).
* **Architectural Pattern:** Adhere to the intended service layer pattern where services manage sessions and business logic, and repositories handle data I/O.
* **Resolver Design:** Align resolver implementation with the existing `prd.resolver.md` as much as feasible, adapting where necessary due to new insights.
* **Integration Testing:** Continue using the remote Supabase instance (`sra_integration_test`) for integration tests, with strict adherence to safety protocols.

## Risks

* **Refactoring Complexity:** Untangling the current state of `search.py`, `services.py`, and `repositories.py` may reveal deeper issues than currently apparent.
* **Resolver Logic Definition:** Defining the precise logic for `resolve_screening_conflict` and ensuring robust interaction with the resolver LLM may be challenging.
* **Test Development Time:** Writing comprehensive unit and integration tests will require significant developer effort.
* **Data State:** Existing data in development or test databases might be inconsistent due to previous partial runs; strategies for handling or resetting this data may be needed.
* **Integration Test Environment:** Any misconfiguration or accidental misuse of the integration testing environment could impact the shared test database. Strict adherence to protocols is key.

## 6. Relevant Research / Current State Artifacts

* **Resolver PRD:** `docs/ai/prd.resolver.md`
* **Key Code Files (illustrating current issues):**
    * Models: `src/sr_assistant/core/models.py`
    * Repositories: `src/sr_assistant/core/repositories.py` (and its unit tests `tests/unit/core/test_repositories.py`)
    * Services: `src/sr_assistant/app/services.py`
    * UI Pages: `src/sr_assistant/app/pages/search.py`, `src/sr_assistant/app/pages/screen_abstracts.py`
    * Agent Logic: `src/sr_assistant/app/agents/screening_agents.py`
* **Testing Framework & Existing Tests:**
    * Test Config: `tests/conftest.py`
    * Integration Tests: `tests/integration/core/test_repositories.py`, `tests/integration/test_resolver.py`, `tests/integration/test_resolver_agent.py` (partially implemented)

## 7. Next Steps Post-Recovery

Once this recovery project is complete and the MPH SR Prototype application is stable and functional:

* The next major planned initiative is the development of a new **Benchmarking Module**. This will likely involve 2-3 new UI pages and associated backend logic, allowing users to compare different screening strategies or models. Details for this will be outlined in a separate project brief.
* Future epics will address the implementation of search functionality for **Embase and Scopus**, building upon the generic `SearchResult` infrastructure stabilized in this recovery.

This current recovery project is foundational for these future enhancements to the MPH SR Prototype.
