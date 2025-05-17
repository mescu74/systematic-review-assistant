---
description: PRD for the automated screening decision resolver feature.
---

# PRD: Automated Screening Decision Resolver

**Version:** 1.1
**Status:** Draft
**Date:** 2025-05-09
**Author:** Product Manager Agent

## 1. Overview

### 1.1. Introduction

The current abstract screening process uses two LLM-driven reviewers with distinct strategies ("conservative" and "comprehensive") to assess PubMed abstracts. While this provides valuable dual perspectives, disagreements inevitably arise (one includes, the other excludes). Currently, resolving these disagreements requires manual intervention, slowing down the screening process. This feature introduces an automated third step: a "resolver" LLM agent that analyzes disagreements and makes a final inclusion/exclusion decision.

### 1.2. Goals

- Automatically resolve disagreements between the conservative and comprehensive screening decisions for abstracts.
- Improve the efficiency of the abstract screening phase by reducing manual resolution effort.
- Maintain consistency in final screening decisions.
- Store and display the resolver's decision and reasoning alongside the original reviewers' outputs.
- Integrate the resolution process seamlessly into the existing screening workflow.

### 1.3. Non-Goals

- Resolving cases where both initial reviewers are "UNCERTAIN".
- Providing a UI for manually overriding the resolver's decision (in this iteration).
- Implementing complex error handling beyond basic logging and potential retries for the resolver LLM call.
- Fine-tuning the resolver LLM model.
- Handling disagreements in full-text screening (this focuses on abstract screening).

## 2. User Stories

- **As a researcher (user),** I want disagreements between the conservative and comprehensive abstract screening results to be automatically resolved by a third AI agent, so I don't have to manually review conflicts.
- **As a researcher (user),** I want to see the final, resolved screening decision for each abstract clearly displayed on the results page.
- **As a researcher (user),** I want to be able to view the reasoning provided by the resolver agent for its final decision, especially when it overrides one of the initial reviewers.
- **As a system administrator,** I want the resolution process, including inputs and outputs, to be logged for monitoring and auditing purposes.

## 3. Functional Requirements

### FR1: Identify Screening Disagreements

- The system must identify `PubMedResult` records within a completed screening batch where the `conservative_result.decision` and `comprehensive_result.decision` differ, or where at least one decision is `UNCERTAIN` if the other is not identical.
- **Definition of Disagreement (v1.1):** A disagreement occurs and triggers the resolver if:
    - One decision is `INCLUDE` and the other is `EXCLUDE` (and vice-versa).
    - One decision is `INCLUDE` and the other is `UNCERTAIN` (and vice-versa).
    - One decision is `EXCLUDE` and the other is `UNCERTAIN` (and vice-versa).
    - Both decisions are `UNCERTAIN`.
- This identification should happen after a batch of abstracts has been processed by the initial dual-reviewer chain (`screen_abstracts_chain`).

### FR2: Batch Disagreements for Resolution

- Identified disagreements should be collected into batches suitable for processing by the resolver agent (e.g., 10-20 items per batch).

### FR3: Prepare Resolver Input

- For each disagreement, the system must gather the necessary context:
  - Full `PubMedResult` data (PMID, title, abstract, journal, year, etc.).
  - Relevant `SystematicReview` protocol data (background, research question, PICO components from `criteria_framework_answers`, exclusion criteria string).
  - The complete `ScreeningResult` object from the "conservative" reviewer.
  - The complete `ScreeningResult` object from the "comprehensive" reviewer.
- This context must be formatted according to the input schema expected by the resolver prompt (`resolver_prompt` in `screening_agents.py`).

### FR4: Invoke Resolver Chain

- The system must invoke the resolver LLM chain (using the defined `RESOLVER_MODEL` and `resolver_prompt` with structured output (targetting `ResolverOutputSchema`)) for each item in the disagreement batch.
- LLM calls should include appropriate error handling (e.g., retries as defined for the screening chains).

### FR5: Process Resolver LLM Output

- The system must parse the structured output (conforming to `ResolverOutputSchema`) from the resolver LLM for each resolved item.
- The `ResolverOutputSchema` (as defined in `src/sr_assistant/core/schemas.py`) currently includes fields like `resolver_decision`, `resolver_reasoning`, `resolver_confidence_score`, and `resolver_include`.
- **Update (from original v1.0):** The schema (and corresponding prompt) must ensure a single, final decision field: `final_decision` (or `resolver_decision` as in the current schema) of type `ScreeningDecisionType`. The prompt must instruct the LLM to provide this single decision (`INCLUDE`, `EXCLUDE`, or potentially `UNCERTAIN` if resolution isn't clear based on new FR1 scope).

### FR6: Store Resolution Data

- A new database table, `screening_resolutions`, must be created to store the resolver's output.
- **Schema:**
  - `id`: UUID (Primary Key)
  - `pubmed_result_id`: UUID (Foreign Key to `pubmed_results.id`)
  - `systematic_review_id`: UUID (Foreign Key to `systematic_reviews.id`)
  - `conservative_result_id`: UUID (Foreign Key to `screening_results.id`)
  - `comprehensive_result_id`: UUID (Foreign Key to `screening_results.id`)
  - `resolver_decision`: `ScreeningDecisionType` (Enum: INCLUDE, EXCLUDE, UNCERTAIN) - The final decision from the resolver.
  - `resolver_reasoning`: TEXT - The reasoning provided by the resolver.
  - `resolver_confidence_score`: FLOAT - The confidence score from the resolver.
  - `resolver_model_name`: VARCHAR - The name of the LLM used for resolution.
  - `created_at`: TIMESTAMPTZ (default now())
  - `trace_id`: `uuid.UUID | None` (optional, for linking to LangSmith trace)
- The `PubMedResult` table must be updated with:
  - `final_decision`: `ScreeningDecisionType | None` (Enum, nullable) - Stores the final outcome after potential resolution. Default is NULL.
  - `resolution_id`: `uuid.UUID | None` (Foreign Key to `screening_resolutions.id`, nullable) - Links to the resolution record if one exists.
- After successful resolution, the system must create a new `ScreeningResolution` record and update the corresponding `PubMedResult` record's `final_decision` and `resolution_id` fields.

### FR7: Update Screening Results UI

- The screening results display (likely on the `screen_abstracts.py` page) must be updated.
- It should prioritize displaying the `final_decision` from the `PubMedResult` if available.
- If a `resolution_id` exists, the UI should visually indicate that the decision was resolved (e.g., an icon, different styling).
- There should be a mechanism (e.g., a tooltip, an expandable section, a link) to view the `resolver_reasoning` stored in the `ScreeningResolution` record.

### FR8: Trigger Mechanism

- The resolution process should be triggered automatically after each batch screening process completes within the `screen_abstracts.py` page workflow.
- The process should identify disagreements from the just-completed batch, batch them, run the resolver, and update the database.
- This should happen asynchronously or in a way that doesn't significantly block the user from viewing the initial (unresolved) results of the batch. A status indicator for the resolution process should be shown.

## 4. Technical Design & Tasks

### TD1: Data Model Changes (`src/sr_assistant/core/models.py`)

- **Task:** Define `ScreeningResolution` SQLModel table with fields specified in FR6.
- **Task:** Add relationship from `ScreeningResolution` back to `PubMedResult` (one-to-one).
- **Task:** Add `final_decision: ScreeningDecisionType | None` field to `PubMedResult` model.
- **Task:** Add `resolution_id: uuid.UUID | None` field (FK to `screening_resolutions.id`) to `PubMedResult` model.
- **Task:** Add relationship from `PubMedResult` to `ScreeningResolution`.
- **Task:** Generate Alembic migration script for database schema changes.

### TD2: Resolver Chain Definition (`src/sr_assistant/app/agents/screening_agents.py`)

- **Task:** Refine `ResolverOutput` Pydantic schema to include `final_decision: ScreeningDecisionType` and remove/update the `include` list.
- **Task:** Update `resolver_prompt` (both system and human parts) to instruct the LLM to provide the `final_decision` and other fields required by the updated `ResolverOutput`. Ensure all necessary context fields (PICO, original results) are included as placeholders.
- **Task:** Define the complete `resolver_chain = resolver_prompt | resolver_model_struct_output`. Add retries similar to screening chains.

### TD3: Repository Layer (`src/sr_assistant/core/repositories.py`)

- **Task:** Create `ScreeningResolutionRepository` class inheriting from `BaseRepository[ScreeningResolution]`.
- **Task:** Implement `add` and potentially `add_all` methods for `ScreeningResolutionRepository`.
- **Task:** Implement `get_by_pubmed_id` method for `ScreeningResolutionRepository`.
- **Task:** Modify `PubMedResultRepository.update` method (or add a specific method) to allow updating `final_decision` and `resolution_id`.

### TD4: Workflow Integration (`src/sr_assistant/app/pages/screen_abstracts.py`)

- **Task:** Create a new function `resolve_batch_disagreements(batch_results: list[ScreenAbstractResultTuple], review: SystematicReview)`.
- **Task:** Inside `resolve_batch_disagreements`:
  - Iterate through `batch_results` to identify disagreements based on FR1 logic.
  - If disagreements exist:
    - Prepare resolver input data for each disagreement (FR3).
    - Batch the inputs.
    - Invoke `resolver_chain.batch(...)` (FR4).
    - Process results:
      - Parse `ResolverOutput` (FR5).
      - Create `ScreeningResolution` model instances.
      - Use `ScreeningResolutionRepository` to save resolution records (FR6).
      - Use `PubMedResultRepository` to update `final_decision` and `resolution_id` on the `PubMedResult` records (FR6).
    - Handle potential errors during the resolver batch call.
- **Task:** Modify the main workflow in `screen_abstracts.py`:
  - After the call to `screen_abstracts_batch` and processing its results:
    - Display a status indicator for the resolution step.
    - Call `resolve_batch_disagreements` with the results.
    - Update the status indicator upon completion or error.
    - Ensure the UI display logic (TD5) uses the potentially updated `PubMedResult` data.

### TD5: UI Update (`src/sr_assistant/app/pages/screen_abstracts.py`)

- **Task:** Modify the data loading/preparation for the results display (e.g., the dataframe) to include `final_decision` and `resolution_id` from the `PubMedResult` objects.
- **Task:** Update the display logic:
  - Show `final_decision` if it exists, otherwise show the original decisions (or indicate conflict).
  - Add a visual cue (e.g., icon, 'Resolved' tag) if `resolution_id` is present.
- **Task:** Implement a way to show `resolver_reasoning`. Options:
  - Add a column to the dataframe (might be too verbose).
  - Use a tooltip on the visual cue.
  - Add an expander or button that fetches and displays the `ScreeningResolution` details (reasoning, confidence, model) using the `resolution_id`.

### TD6: Testing

- **Task:** Write unit tests for `resolve_batch_disagreements` function logic (mocking chain invocation and repositories).
- **Task:** Write unit tests for `ScreeningResolutionRepository` methods.
- **Task:** Update unit tests for `PubMedResultRepository` if methods were changed.
- **Task:** Write integration tests for the resolver chain (`resolver_chain`) itself to ensure it parses inputs and calls the LLM correctly (can mock the LLM response).

## 5. Epic Structure

Epic-1: Resolver Implementation (Current)
Epic-2: Database Integration Testing (Future)

## 6. Story List

### Epic-1: Resolver Implementation

Story-1: Data Model and Repository Setup
Story-2: Resolver Chain Definition and Prompting
Story-3: Workflow Integration and Triggering
Story-4: UI Updates for Resolved Decisions

### Epic-2: Database Integration Testing

Story-1: Implement Resolver DB Integration Tests

## 7. Future Considerations

- Handling `UNCERTAIN` cases in disagreement logic.
- Allowing user configuration of the resolver model and prompt.
- Providing a UI for users to review and potentially override resolved decisions.
- Implementing more robust error handling and reporting for resolver failures.
- Performance monitoring and optimization for large-scale resolution.
- Adding Time (T) and Study Design (S) to PICO (PICOT/PICOS) if needed later.

## 8. Success Metrics

- Reduction in the number of abstracts requiring manual conflict resolution (measured by tracking disagreements before/after).
- Successful storage of resolution data in the `screening_resolutions` table.
- Correct display of `final_decision` and resolution indicators in the UI.
- Qualitative feedback from users on the accuracy and helpfulness of the automated resolution.
- Execution time for the resolution step remains within acceptable limits.

## 9. Change Log

| Version | Date       | Author                | Change Description                                                                                                |
|---------|------------|-----------------------|-------------------------------------------------------------------------------------------------------------------|
| 1.1     | 2025-05-08 | Product Manager Agent | Added Epic structure, clarified resolver input (FR3), output (FR5 - `ResolverOutputSchema`), storage (FR6), and UI (FR7). Removed redundant sections. |
| 1.2     | 2025-05-09 | Product Manager Agent | Corrected FR4/FR5 to reference `ResolverOutputSchema` instead of non-existent `ResolverOutput`.             |
| 1.0     | 2025-04-04 | AIDE (AI Assistant)   | Initial draft.                                                                                                    |

---
