# Automated Screening Decision Resolver PRD

## Goals and Background Context

### Goals

- Automatically resolve disagreements between the conservative and comprehensive screening decisions for abstracts.
- Improve the efficiency of the abstract screening phase by reducing manual resolution effort.
- Maintain consistency in final screening decisions by storing and displaying the resolver's decision and reasoning.
- Seamlessly integrate the resolution process into the existing screening workflow.

### Background Context

The SRA's abstract screening process uses two LLM reviewers ("conservative" and "comprehensive"). This dual-reviewer system inevitably leads to disagreements (e.g., one includes, the other excludes) that require manual resolution, creating a bottleneck. This PRD outlines an automated "resolver" LLM agent to analyze these disagreements and make a final, consistent decision, thereby streamlining the entire screening workflow.

### Change Log

| Date       | Version | Description                   | Author |
| :--------- | :------ | :---------------------------- | :----- |
| 2024-07-30 | 2.0     | Migrated to v4 PRD Template   | AIDE   |
| 2025-05-09 | 1.2     | Corrected schema reference    | PM Agent |
| 2025-05-08 | 1.1     | Clarified requirements        | PM Agent |
| 2025-04-04 | 1.0     | Initial Draft                 | AIDE   |

## Requirements

### Functional

- **FR1:** The system must identify `SearchResult` records where the initial two screening decisions disagree. A disagreement is defined as: `INCLUDE` vs. `EXCLUDE`, or any combination of `INCLUDE`/`EXCLUDE` with `UNCERTAIN`.
- **FR2:** For each disagreement, the system must prepare an input for the resolver agent containing the full `SearchResult` context, the `SystematicReview` protocol, and the outputs from both initial reviewers.
- **FR3:** The system must invoke the resolver LLM chain for each disagreement.
- **FR4:** The system must parse the structured output (`ResolverOutputSchema`) from the resolver, which contains the final decision and reasoning.
- **FR5:** A new `screening_resolutions` table must be created to store the detailed output from the resolver.
- **FR6:** The `search_results` table must be updated with `final_decision` and a foreign key (`resolution_id`) linking to the new `screening_resolutions` record.
- **FR7:** The screening UI must be updated to display the `final_decision` and provide access to the resolver's reasoning.
- **FR8:** The resolution process must be triggered automatically after the initial batch screening process completes.

### Non Functional

- **NFR1:** The resolution process must be reliable and produce consistent outcomes.
- **NFR2:** All resolution data must be stored for auditing and transparency.
- **NFR3:** The resolver workflow should not significantly degrade the user experience or performance of the screening page.

## User Interface Design Goals

N/A. Focus is on backend implementation and data integration with the existing UI.

## Technical Assumptions

### Repository Structure

Monorepo (existing)

### Service Architecture

Service-Oriented Monolith (existing)

### Testing requirements

- Unit tests for the resolver workflow logic.
- Unit tests for the new `ScreeningResolutionRepository`.
- Integration tests for the resolver LLM chain.

### Additional Technical Assumptions and Requests

- The resolver will be implemented as a LangChain Expression Language (LCEL) chain.
- Database changes will be managed via Alembic migrations.

## Epics

### Epic 1: Resolver Backend Implementation

**Goal:** Develop the complete backend functionality for the automated screening conflict resolver, including data model changes, the resolver agent itself, and the service-layer logic to orchestrate the process.

#### Story 1.1: Update Data Models

As a developer, I want to update the database schema and models to store resolution data, so that the resolver's decisions can be persisted.

##### Acceptance Criteria

- 1: `ScreeningResolution` SQLModel and table are created.
- 2: `SearchResult` model and `search_results` table are updated with `final_decision` and `resolution_id` fields.
- 3: An Alembic migration script is successfully generated and applied.
- 4: A `ScreeningResolutionRepository` is created to manage data access.

#### Story 1.2: Build Resolver Agent

As a developer, I want to build the resolver LLM chain, so that it can analyze disagreements and produce a structured final decision.

##### Acceptance Criteria

- 1: The `ResolverOutputSchema` Pydantic model is finalized.
- 2: A `resolver_prompt` is engineered to guide the LLM to provide the required output based on the provided context.
- 3: A `resolver_chain` is constructed using LCEL, combining the prompt, a structured-output model, and retries.

#### Story 1.3: Integrate Resolver into Screening Workflow

As a developer, I want to integrate the resolver chain into the screening page's backend workflow, so that it runs automatically after the initial screening is complete.

##### Acceptance Criteria

- 1: A function is created in `screen_abstracts.py` to identify disagreements from a batch.
- 2: This function prepares inputs and invokes the `resolver_chain.batch()` method.
- 3: The function processes the resolver's output, saving the data using the `ScreeningResolutionRepository` and updating the `SearchResult` records.
- 4: The function is called automatically after the initial screening batch finishes.

### Epic 2: UI and Testing

**Goal:** Update the user interface to reflect the outcomes of the automated resolution process and ensure the entire feature is covered by tests.

#### Story 2.1: Update Screening UI

As a researcher, I want to see the final resolved decision and the resolver's reasoning in the UI, so I can understand the outcome.

##### Acceptance Criteria

- 1: The screening results table displays the `final_decision` from the `SearchResult` model.
- 2: A visual indicator is present for results that have been processed by the resolver.
- 3: The resolver's reasoning is accessible, for instance, through a tooltip or an expandable section.

#### Story 2.2: Implement Tests

As a developer, I want to write unit and integration tests for the resolver feature, so I can ensure its correctness and prevent regressions.

##### Acceptance Criteria

- 1: Unit tests for the disagreement identification and workflow logic are created.
- 2: Unit tests for the `ScreeningResolutionRepository` are created.
- 3: Integration tests for the `resolver_chain` are created to validate its behavior (LLM can be mocked).

## Checklist Results Report

[[LLM: This section will be populated after running the `pm-checklist` task.]]

## Next Steps

[[LLM: This section will contain prompts for the Design Architect and Architect.]]

---
