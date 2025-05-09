# Epic 3: Comprehensive Testing and Database Integrity

**Goal:** Establish a robust testing framework (unit, integration, E2E for the core recovery flow) and ensure database integrity, including clear documentation for safe integration testing protocols. This supports the PRD objectives to "Implement comprehensive test coverage" and ensure the application is "stable, robust, and reliable."

**Deployability:** This epic is the final stage of the recovery project. Its completion signifies that the core functionalities (Search, Screen, Resolve) are not only working but are also covered by adequate tests and that data integrity has been reviewed. It ensures the application is in a maintainable and verifiable state before any new feature development (like the benchmarking module) begins.

## Epic-Specific Technical Context

This epic focuses on solidifying the testing landscape and ensuring data integrity post-refactoring. Key technical aspects include:
- **Test Coverage:** Achieving and verifying the >80% unit test coverage target for key modules (`services.py`, `repositories.py`, `models.py`).
- **Integration Test Suite:** Ensuring all integration tests, including those for `SearchService`, `SearchResultRepository`, `ScreeningResolutionRepository`, and the `resolver_chain` (with its various scenarios), are complete and pass reliably.
- **End-to-End (E2E) Testing:** Implementing at least one E2E test that simulates the primary user workflow: PubMed search -> abstract screening (with dual reviewers) -> automated conflict resolution.
- **Documentation:** Creating clear, actionable documentation for setting up and running integration tests safely against the `sra_integration_test` database, emphasizing environment variable management (`.env.test`) and schema cleaning (`conftest.py`).
- **Database Integrity:** Reviewing database schema, relationships, and potential data inconsistencies that might have arisen during the recovery process.

## Local Testability & Command-Line Access

- **Local Development:** Developers will primarily be running and debugging tests locally. The Streamlit application (`make run`) should be runnable to manually verify E2E scenarios if needed, connecting to the Supabase-hosted "prototype" database.
- **Command-Line Testing:** All tests (unit, integration, E2E) will be executed via `uv run pytest [options] <test_file_or_directory>` or the corresponding `make` targets (`make test.unit`, `make test.integration`, `make test.all`).
- **Code Formatting & Linting:** Linting errors can often be auto-fixed using `uv run ruff check --fix <path/to/file_or_directory>`. Code formatting is done using `uv run ruff format <path/to/file_or_directory>` (typically run after auto-fixing). The `make lint` (or `make ruff.fix`) and `make format` targets can also be used.
- **Environment Testing:** Unit tests mock DB interactions. Integration and E2E tests run against the Supabase-hosted `sra_integration_test` database.
- **Testing Prerequisites:**
    - Access to Supabase-hosted PostgreSQL instance (prototype and `sra_integration_test` databases).
    - Correct environment variables (`.env`, `.env.local`, `.env.test`).
    - Python 3.12 environment (`make install` or `uv sync`).
    - Test data fixtures as needed for various test scenarios.

## Story List

### Story 3.1: Finalize and Document Unit Test Suite

- **User Story / Goal:** As a developer, I want a comprehensive suite of unit tests with >80% coverage for critical backend modules, so that individual components are verified for correctness and maintainability is improved.
- **Detailed Requirements:**
  - Review and augment unit tests for all services (`services.py` - e.g., `SearchService`, and any logic extracted for resolver).
  - Review and augment unit tests for all repositories (`repositories.py` - `SearchResultRepository`, `ScreeningResolutionRepository`, etc.), ensuring database interactions are mocked appropriately.
  - Review and augment unit tests for all data models (`models.py`) if they contain custom logic or validation.
  - Ensure unit tests are fast, isolated, and do not require external dependencies (like a live database or LLM API).
  - Generate a coverage report and ensure the target of >80% is met for the specified modules.
  - Document any complex mocking strategies or setup required for unit tests.
- **Acceptance Criteria (ACs):**
  - AC1: Unit test coverage for `services.py` (specifically `SearchService` and resolver-related service logic) is >80%.
  - AC2: Unit test coverage for `repositories.py` (specifically `SearchResultRepository`, `ScreeningResolutionRepository`) is >80%.
  - AC3: Unit test coverage for `models.py` (for any custom logic) is >80%.
  - AC4: All unit tests pass reliably when run via `uv run pytest tests/unit` or `make test.unit`.
  - AC5: Coverage report confirms the >80% target for the specified modules.
  - AC6: Documentation for running unit tests and any complex setups is clear.
- **Dependencies:** Relies on the codebase stabilized in Epic 1 and Epic 2.

---

### Story 3.2: Complete and Document Integration Test Suite

- **User Story / Goal:** As a developer, I want a robust suite of integration tests that verify the interactions between different components (services, repositories, database, LLM agent), so that key workflows are confirmed to be functional in an integrated environment.
- **Detailed Requirements:**
  - Ensure all integration tests for `SearchService` and `SearchResultRepository` (from Epic 1) are complete and passing.
  - Ensure all integration tests for `ScreeningResolutionRepository`, `SearchResult` updates (resolver fields), and the `resolver_chain` (from Epic 2) are complete and passing. This includes:
    - Tests for `resolver_chain` running against the actual LLM API.
    - Tests for different conflict scenarios that **must** trigger the resolver:
        - `INCLUDE` vs. `EXCLUDE` (and vice-versa).
        - `INCLUDE` vs. `UNCERTAIN` (and vice-versa).
        - `EXCLUDE` vs. `UNCERTAIN` (and vice-versa).
        - `UNCERTAIN` vs. `UNCERTAIN`.
    - *(Note: This handling of UNCERTAIN cases supersedes the initial scope in `docs/prd-resolver.md` FR1 and implies that the resolver should attempt to make a definitive decision, or explicitly output UNCERTAIN if it cannot resolve based on the input.)*
  - Write any missing integration tests identified during Epic 1 and Epic 2 implementation.
  - Ensure integration tests correctly use the `sra_integration_test` database and manage test data/state appropriately (e.g., via fixtures, cleanup in `conftest.py`).
  - Document how to run integration tests and any specific configurations needed.
- **Acceptance Criteria (ACs):**
  - AC1: All integration tests for search functionality (service-repository-DB) pass.
  - AC2: All integration tests for resolver functionality (agent-service-repository-DB) pass, including tests against the live LLM and specific conflict scenarios.
  - AC3: Integration tests demonstrate correct data persistence and retrieval for `SearchResult` and `ScreeningResolution` including their relationships.
  - AC4: All integration tests pass reliably when run via `uv run pytest -m integration tests/integration` or `make test.integration`.
  - AC5: Documentation for running integration tests is clear and accurate.
- **Dependencies:** Relies on the codebase stabilized in Epic 1 and Epic 2.

---

### Story 3.3: Implement and Document End-to-End (E2E) Test for Core Workflow

- **User Story / Goal:** As a QA/developer, I want at least one automated end-to-end test for the primary user workflow (PubMed search -> screening -> resolution) using Streamlit's UI testing framework, so that the entire integrated system is verified from a user-like perspective.
- **Detailed Requirements:**
  - Design an E2E test scenario using Streamlit's UI testing framework that covers:
    - Performing a PubMed search with a specific query known to yield a small, predictable number of results.
    - Initiating abstract screening for these results via UI interactions.
    - Ensuring the dual-reviewer screening agents are invoked (verified by checking intermediate state or UI updates if possible).
    - Designing test data and/or agent prompt configurations that reliably lead to a disagreement scenario (e.g., `INCLUDE` vs. `EXCLUDE`, or cases involving `UNCERTAIN` as per Story 3.2) for at least one abstract.
    - Verifying that the resolver agent is triggered for the disagreement (e.g., via UI status change or final outcome).
    - Verifying that a `final_decision` is populated and a `ScreeningResolution` record is created (by querying the database post-test or checking UI).
    - Verifying that the UI in `screen_abstracts.py` reflects the resolved state for the specific abstract.
  - The E2E test must use Streamlit's UI testing framework and run headlessly if possible.
  - The test should be designed with future CICD execution in mind, although getting it running in CICD is not part of this recovery epic.
  - The E2E test should use the `sra_integration_test` database.
  - Document the E2E test scenario, setup, and execution steps.
- **Acceptance Criteria (ACs):**
  - AC1: At least one E2E test for the PubMed search -> screen -> resolve workflow is implemented.
  - AC2: The E2E test successfully executes, demonstrating the correct flow of data and operations through the system for the defined scenario.
  - AC3: The E2E test verifies correct data state in the database (e.g., `SearchResult`, `ScreeningResolution` records) at the end of the flow.
  - AC4: The E2E test passes reliably when run (e.g., via a specific pytest marker and command).
  - AC5: Documentation for the E2E test (scenario, setup, execution) is clear.
- **Dependencies:** Relies on stable functionality from Epic 1 and Epic 2, and completed integration tests (Story 3.2).

---

### Story 3.4: Formalize and Document Safe Integration Testing Protocol

- **User Story / Goal:** As a developer, I need clear, actionable documentation on setting up and running integration tests safely against the `sra_integration_test` database, so that tests are reliable and do not interfere with development or prototype data.
- **Detailed Requirements:**
  - Create/update a dedicated documentation file (e.g., `docs/RUNNING_INTEGRATION_TESTS.md` or a section within `CONTRIBUTING.md`, *not* `docs/testing-strategy.md` which is an Architect artifact) that covers:
    - Purpose of the `sra_integration_test` Supabase database.
    - Required environment variables in `.env.test` (and how to obtain necessary credentials like DB connection strings safely).
    - How the schema is managed/cleaned for integration tests (explaining the role of `tests/conftest.py`).
    - Step-by-step instructions for running all integration tests and specific integration test files/markers.
    - Emphasis on the criticality of *not* running these tests against the prototype/production database.
    - Troubleshooting common issues related to integration test setup.
- **Acceptance Criteria (ACs):**
  - AC1: A dedicated section or document for integration testing protocol is created/updated.
  - AC2: The documentation clearly explains the setup for `.env.test`.
  - AC3: The role of `conftest.py` in schema management for tests is documented.
  - AC4: Instructions for running integration tests are accurate and easy to follow.
  - AC5: Warnings about targeting the correct database are prominent.
- **Dependencies:** General project stability.

---

### Story 3.5: Database Integrity and Final Review

- **User Story / Goal:** As a project team member, I want a final review of the database schema and data integrity, so that we can ensure the recovered application is built on a solid data foundation.
- **Detailed Requirements:**
  - Review the complete database schema (all tables, columns, types, relationships, constraints) after all Epic 1 & 2 changes and Alembic migrations.
  - Verify correctness of foreign key relationships, especially those involving `SearchResult`, `ScreenAbstractResult`, and `ScreeningResolution`.
  - Check for any potential data inconsistencies that might have resulted from partial runs or bugs during the recovery development (e.g., orphaned records, inconsistent `final_decision` states if any manual testing occurred on shared prototype DBs prior to full resolver implementation).
  - Ensure all enums used in SQLModel fields are correctly defined in the database (e.g., `ScreeningDecisionType`, `SearchDatabaseSource`).
  - Perform a final check on indexing strategies for key tables to ensure reasonable query performance for common operations.
  - *(Note: Any required updates to formal database schema documentation or data model diagrams resulting from this review will be handled by the Architect, potentially in collaboration with the PM/dev team.)*
- **Acceptance Criteria (ACs):**
  - AC1: Database schema review is completed and any identified issues are documented or addressed.
  - AC2: Foreign key relationships are verified as correct and functional.
  - AC3: Enum types in the database match the definitions in `types.py` and `models.py`.
  - AC4: A check for common data inconsistencies (e.g., orphaned resolution records) is performed, and any findings are addressed or noted.
  - AC5: Indexing strategy is confirmed as reasonable for the expected query patterns of the recovered features.
  - AC6: Foreign Key constraints with cascade delete behavior (if any are intended or exist for relevant relationships like `SearchResult` to its associated `ScreenAbstractResult` or `ScreeningResolution`) are verified to be configured and functioning correctly to prevent orphaned records.
- **Dependencies:** Completion of all development and testing stories in Epic 1, 2, and 3 (up to 3.4).

---

## Change Log

| Change          | Date       | Version | Description             | Author               |
|-----------------|------------|---------|-------------------------|----------------------|
| Initial Draft   | 2025-05-09 | 0.1     | First draft of Epic 3   | Product Manager Agent | 