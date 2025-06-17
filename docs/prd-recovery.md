# SRA Recovery Project Product Requirements Document (PRD)

## Goals and Background Context

### Goals

- Restore the MPH SR Prototype to a stable, robust, and reliable state.
- Correctly implement a generic search infrastructure, initially focused on PubMed, using a well-defined service layer.
- Deliver a fully functional automated screening conflict resolver.
- Achieve comprehensive test coverage for all core components.

### Background Context

The Systematic Review Assistant (MPH SR Prototype) application is currently in a significantly unstable state due to incomplete, concurrent development efforts. A search functionality refactoring and the implementation of an automated screening conflict resolver were left unfinished, resulting in numerous errors and a non-functional screening workflow. This recovery project aims to stabilize the application, complete the features, and establish a solid architectural foundation for future development.

### Change Log

| Date       | Version | Description                   | Author |
| :--------- | :------ | :---------------------------- | :----- |
| 2024-07-30 | 1.2     | Migrated to v4 PRD Template   | AIDE   |
| 2025-05-09 | 1.1     | Clarified resolver scope      | PM Agent |
| 2025-05-09 | 1.0     | Initial Draft                 | PM Agent |

## Requirements

### Functional

- **FR1:** The `search.py` UI must correctly interact with `SearchService` for all PubMed search operations.
- **FR2:** `SearchService` must manage PubMed search logic, database sessions, and repository calls correctly.
- **FR3:** `SearchResultRepository` must contain all methods needed to manage `SearchResult` entities for PubMed.
- **FR4:** PubMed search results must be stored and retrieved consistently using the `SearchResult` model.
- **FR5:** The `resolve_screening_conflict` function must be fully implemented and integrated.
- **FR6:** The `SearchResult` model must be updated with a `final_decision` field (including DB migration).
- **FR7:** The automated conflict resolver's decisions must be correctly stored in `ScreeningResolution` and `SearchResult` tables.
- **FR8:** The resolver workflow, including handling of "UNCERTAIN" decisions, must be fully integrated into the `screen_abstracts.py` page.
- **FR9:** Unit tests must be developed for all services in `services.py`.
- **FR10:** Unit tests for `models.py` and `repositories.py` must be repaired and expanded.
- **FR11:** All integration tests for the resolver feature must be completed and pass reliably.
- **FR12:** A new end-to-end integration test for the PubMed search -> screen -> resolve workflow must be created and pass.

### Non Functional

- **NFR1:** The application must be free of critical runtime errors in the recovered search, screen, and resolve workflows.
- **NFR2:** Core features must function consistently and predictably.
- **NFR3:** Data integrity must be enforced, with correct usage of `source_id` and `source_db` for PubMed data.
- **NFR4:** Code must adhere to defined quality standards and the service layer architecture must be clearly implemented.
- **NFR5:** The system must be structured to support comprehensive and safe unit and integration testing against a dedicated test database.

## User Interface Design Goals

N/A for this recovery project. The focus is on backend stabilization and ensuring existing UI correctly reflects the fixed data flows.

## Technical Assumptions

### Repository Structure

Monorepo (existing)

### Service Architecture

Service-Oriented Monolith (existing)

### Testing requirements

- Unit tests for all services, models, and repositories.
- Integration tests for the conflict resolver.
- One end-to-end integration test for the core PubMed workflow.
- Integration tests must run safely against the `sra_integration_test` Supabase database.

### Additional Technical Assumptions and Requests

- Use existing tech stack: Python, Streamlit, SQLModel, LangChain, Pytest.
- Adhere to the intended service layer pattern where services manage DB sessions.
- Non-critical linter errors may be temporarily ignored to focus on functionality, but critical errors must be fixed.

## Epics

### Epic 1: Search and Service Layer Stabilization

**Goal:** Fully refactor and stabilize `search.py`, `services.py`, and `repositories.py` to correctly implement and utilize the `SearchResult` model and the new service layer architecture for PubMed search operations, ensuring all related functionalities are operational and error-free.

#### Story 1.1: Fix Service and Repository Layers

As a developer, I want to refactor the service and repository layers to manage database sessions correctly and align with the defined architecture, so that data operations are stable and testable.

##### Acceptance Criteria

- 1: `SearchService` correctly handles its own DB session lifecycle.
- 2: `SearchResultRepository` methods are updated to accept a session from the service layer.
- 3: All linter errors in `services.py` and `repositories.py` are resolved.

#### Story 1.2: Stabilize PubMed Search

As a researcher, I want to perform a PubMed search and see the results stored and displayed correctly, so that I can begin the screening process.

##### Acceptance Criteria

- 1: `search.py` UI successfully calls `SearchService`.
- 2: PubMed results are fetched and stored in the `search_results` table with `source_db` = 'PubMed' and `source_id` = PMID.
- 3: The UI correctly displays the fetched results from the database.

### Epic 2: Resolver Agent Implementation and Integration

**Goal:** Complete the implementation and integration of the automated screening conflict resolver, including model updates, backend logic, and UI integration, ensuring it accurately processes conflicts and stores results.

#### Story 2.1: Implement Resolver Logic

As a developer, I want to implement the `resolve_screening_conflict` function and related database changes, so that screening disagreements can be automatically resolved.

##### Acceptance Criteria

- 1: `SearchResult` model is migrated to include `final_decision`.
- 2: The resolver LLM chain is implemented and correctly processes inputs.
- 3: Resolver output is correctly stored in the `screening_resolutions` table and linked back to the `search_results` table.

#### Story 2.2: Integrate Resolver into UI

As a user, I want to see the final, resolved decisions in the screening interface, so that I know the outcome of the automated process.

##### Acceptance Criteria

- 1: The `screen_abstracts.py` workflow correctly invokes the resolver logic after the initial dual-agent screening.
- 2: The final decision and the resolver's reasoning are clearly displayed in the UI.
- 3: The workflow correctly handles cases involving 'UNCERTAIN' decisions from the primary reviewers.

### Epic 3: Comprehensive Testing and Database Integrity

**Goal:** Establish a robust testing framework (unit, integration, E2E for the core recovery flow) and ensure database integrity, including clear documentation for safe integration testing protocols.

#### Story 3.1: Build Unit Test Suite

As a developer, I want a comprehensive suite of unit tests, so that I can verify the correctness of individual components in isolation.

##### Acceptance Criteria

- 1: Unit tests for `SearchService` and `ScreeningService` are written, mocking all external dependencies.
- 2: Unit tests for all `SQLModel` models and `Pydantic` schemas are complete.
- 3: Unit tests for all `Repository` methods are complete.
- 4: Unit test coverage for core logic exceeds 80%.

#### Story 3.2: Create End-to-End Integration Test

As a developer, I want an end-to-end integration test, so that I can validate the entire PubMed search-screen-resolve workflow against a real database.

##### Acceptance Criteria

- 1: An integration test is created that performs a PubMed search, triggers screening, and verifies the conflict resolution.
- 2: The test runs successfully against the `sra_integration_test` database.
- 3: The test includes assertions to verify data integrity at each step of the process.

## Checklist Results Report

[[LLM: This section will be populated after running the `pm-checklist` task.]]

## Next Steps

[[LLM: This section will contain prompts for the Design Architect and Architect.]]
