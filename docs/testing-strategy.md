# Testing Strategy

## 1. Introduction & Goals

This document outlines the comprehensive testing strategy for the Systematic Review Assistant (SRA) application, particularly focusing on the functionalities being addressed during the recovery project. The primary goal of this strategy is to ensure the application's stability, reliability, data integrity, and correctness, aligning with the objectives stated in `docs/prd-recovery.md`.

Key quality goals include:
- Verification of all core user workflows (search, screening, conflict resolution).
- Validation of data persistence and relational integrity.
- Ensuring robust error handling and adherence to defined API contracts.
- Achieving high unit test coverage for critical backend components.

## 2. Scope of Testing

### In Scope:
-   **Core Functionalities:**
    -   PubMed Search: Service logic, data mapping, persistence (`SearchService`, `SearchResultRepository`, `SearchResult` model).
    -   Abstract Screening: Dual reviewer LLM chain invocation (`screen_abstracts_batch`), processing of `ScreeningResponse` and `ScreeningResult`, persistence of `ScreenAbstractResult`.
    -   Conflict Resolution: Resolver LLM chain invocation, processing of `ScreeningResolutionSchema`, persistence of `ScreeningResolution`, and updates to `SearchResult.final_decision`.
-   **Service Layer:** All public methods of `SearchService`, `ReviewService`, and `ScreeningService` as defined in `docs/api-reference.md`.
-   **Repository Layer:** All methods of `SearchResultRepository`, `SystematicReviewRepository`, `ScreenAbstractResultRepository`, and `ScreeningResolutionRepository`.
-   **Data Models & Schemas:** Validation logic inherent in Pydantic schemas (`schemas.py`) and SQLModel constraints (`models.py`).
-   **UI Interactions:** Key user flows in `search.py`, `screen_abstracts.py`, and `protocol.py` via E2E tests.

### Out of Scope (for initial recovery project testing phase):
-   Comprehensive performance testing beyond basic responsiveness checks.
-   Security penetration testing.
-   Usability testing beyond developer/QA assessment of core flows.
-   Testing of legacy code paths not part of the recovery (e.g., old Step1/Step2 modules).
-   Full test automation for all edge cases of LLM outputs (focus on primary success/failure paths and schema adherence).

## 3. Test Levels & Types

### 3.1. Unit Tests

-   **Goal:** Verify individual components (functions, classes, methods) in isolation to ensure they behave as expected.
-   **Scope:**
    -   Services (`services.py`): Business logic, mocking repository calls and external API interactions (e.g., PubMed client, LLM calls).
    -   Repositories (`repositories.py`): Data access logic, mocking database session/engine and actual DB calls. Focus on query construction and result mapping logic if any.
    -   Utility functions and helpers.
    -   Custom validation logic within Pydantic schemas (if any beyond standard Pydantic validation).
    -   Prompt templating logic.
    -   **Property-based testing of Pydantic schema validation/parsing, utility functions, and isolated algorithmic logic within services or repositories using Hypothesis.**
-   **Tools:** Pytest, `pytest-mock` (using the `mocker` fixture), **Hypothesis** (for property-based testing).
-   **Coverage Target:** Aim for >80% statement coverage for `services.py` and `repositories.py` (and `models.py` if it contains significant custom logic).

### 3.2. Integration Tests

-   **Goal:** Verify interactions between different components of the application and with external systems (database, LLM APIs).
-   **Scope:**
    -   Service <> Repository <> Database: Test the full flow from service call to data persistence/retrieval in the actual test database.
    -   Service <> LLM Agent/Chain: Test the invocation of LLM chains (e.g., `screen_abstracts_chain`, `resolver_chain`) from the service layer or UI page handlers.
        -   Some tests may mock the LLM response to ensure consistent behavior and test specific data paths.
        -   Critical path tests (e.g., resolver chain for key conflict scenarios) MUST run against live LLM APIs (using appropriate API keys and respecting quotas).
    -   Page <> Service: Verify that Streamlit page handlers correctly call service methods and handle their responses/exceptions (can be partially covered by E2E or more focused integration tests).
-   **Tools:** Pytest, live `sra_integration_test` PostgreSQL database (Supabase-hosted), live LLM APIs (OpenAI, Google AI as per `docs/tech-stack.md`).
-   **Protocol:** Adhere to `docs/RUNNING_INTEGRATION_TESTS.md` (to be created/updated by dev team as per Epic 3, Story 3.4) for setup and execution.

### 3.3. End-to-End (E2E) Tests

-   **Goal:** Verify complete user workflows from the UI perspective, simulating real user interaction.
-   **Scope:** The primary workflow: Perform PubMed search -> Select results for screening -> Initiate abstract screening (dual reviewers) -> Ensure disagreements are identified -> Trigger conflict resolution -> Verify final decisions and resolver reasoning are displayed.
-   **Tools:** Streamlit UI Testing Framework (likely driven by Pytest).
-   **Database:** Uses the `sra_integration_test` database.

### 3.4. Contract Tests (LLM Interactions)

-   **Goal:** Ensure that the inputs provided to LLM chains (prompts and their variables) and the outputs received (structured Pydantic schemas like `ScreeningResponse`, `ScreeningResolutionSchema`) adhere to their expected formats and constraints.
-   **Implementation:** These are often integrated within the integration tests for the LLM chains themselves. Assertions will be made on the structure and types of data passed to and received from the LLM components.

## 4. Testing Environment & Data

-   **Development Environment:** Local machine with Python 3.12, `uv` for package management, connection to Supabase-hosted `postgres` (the **prototype**) database.
-   **Unit Testing Environment:** Local machine or CI environment (as per `.github/workflows/pr.yml`), dependencies mocked. No live database or external API access.
-   **Integration & E2E Testing Environment:** Local machine or CI environment, connecting to the dedicated `sra_integration_test` Supabase-hosted PostgreSQL database. Requires API keys for LLMs.
-   **Test Data:**
    -   Unit tests: Use mocked data or small, self-contained fixtures.
    -   Integration/E2E tests: Use a combination of pre-defined test data fixtures (e.g., specific `SystematicReview` protocols, `SearchResult` examples known to cause conflicts) and data generated during test execution. The `sra_integration_test` database will be reset/cleaned before test runs (managed by `tests/conftest.py`).
-   **Environment Variables:** A `.env.test` file will store configurations specific to the testing environment (DB connection strings for `sra_integration_test`, API keys for test accounts if different from dev).

## 5. Test Execution & Reporting

-   **Execution Commands:**
    -   All tests: `make test.all` or `uv run pytest`
    -   Unit tests: `make test.unit` or `uv run pytest tests/unit`
    -   Integration tests: `make test.integration` or `uv run pytest -m integration tests/integration`
    -   Specific E2E tests might use custom markers/commands.
-   **Reporting:**
    -   Pytest will generate console output and potentially HTML reports.
    -   Test coverage will be reported using `pytest-cov` (e.g., `uv run pytest --cov=src/sr_assistant`).

## 6. Roles & Responsibilities (High-Level)

-   **Developers:** Responsible for writing and maintaining unit tests for their code and contributing to integration tests for features they implement.
-   **Architect:** Defines the overall testing strategy (this document), reviews test plans and significant test implementations.
-   **Scrum Master (SM):** Ensures testing tasks are included in story definitions and tracked.
-   **QA Team/Assigned Developers:** Primarily responsible for designing, implementing, and maintaining E2E tests.

## 7. Automation Strategy

-   All unit and integration tests SHOULD be fully automated and ARE executed in the CI/CD pipeline (`.github/workflows/pr.yml`) on pull requests and pushes to main.
-   E2E tests WILL be automated using the Streamlit testing framework.
-   **Future Goal:** Expand CI/CD to include E2E tests once they are stable and can run reliably in a headless environment if necessary.

## 8. Specific Focus Areas for Recovery Project

-   Thorough testing of the refactored `SearchService` (PubMed search flow, data mapping, session management).
-   Comprehensive testing of the `ScreeningService`, including the new conflict resolution methods and alignment of existing methods.
-   Verification of all defined Pydantic schemas in `schemas.py` through their usage in tests.
-   Validation of SQLModel definitions and database interactions, especially for `SearchResult`, `ScreenAbstractResult`, and `ScreeningResolution`, including the new `SearchResult.final_decision` field and FK relationships.
-   Testing all documented conflict scenarios that should trigger the resolver (as per `docs/epic3-recovery-testing-and-integrity.md`, Story 3.2).
-   **Consider applying property-based testing with Hypothesis, particularly for Pydantic schemas with complex validation logic (e.g., `ScreeningResponse`, `ScreeningResolutionSchema`) and critical utility functions, to enhance robustness against diverse inputs.**

## 9. Document Maintenance

This testing strategy document will be reviewed and updated as the project evolves, new features are added, or if significant changes to the architecture or technology stack occur.

## Change Log

| Change          | Date       | Version | Description             | Author          |
|-----------------|------------|---------|-------------------------|-----------------|
| Initial Draft   | 2025-05-12 | 0.1     | First draft of Testing Strategy. | Architect Agent |
| Tooling Update  | 2025-05-12 | 0.1.1   | Corrected Unit Test tooling to specify `pytest-mock` instead of `unittest.mock`. | Architect Agent |
| Env & CI Update | 2025-05-12 | 0.1.2   | Corrected prototype DB name and clarified CI execution for unit/integration tests. | Architect Agent |
| Added Hypothesis| 2025-05-12 | 0.2.0   | Incorporated Hypothesis for property-based testing into the strategy. | Architect Agent | 