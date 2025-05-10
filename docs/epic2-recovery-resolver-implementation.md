# Epic 2: Resolver Agent Implementation and Integration

**Goal:** Complete the implementation and integration of the automated screening conflict resolver, including model updates, backend logic, and UI integration, ensuring it accurately processes conflicts and stores results. This directly supports the PRD objectives to "Deliver a fully functional automated screening conflict resolver" and improve screening efficiency.

**Deployability:** This epic builds upon the stabilized search and service layer from Epic 1. It introduces a significant new piece of functionality (automated conflict resolution) that, once complete, will make the abstract screening process more autonomous. It is deployable in the sense that the screening workflow will now have an additional automated step for disagreements.

## Epic-Specific Technical Context

This epic implements the automated screening conflict resolver as detailed in the existing resolver PRD (`docs/prd-resolver.md` - content is relevant even if path is deprecated). Key technical aspects include:
- **Data Model Changes:** Verifying and updating the existing `ScreeningResolution` table in `models.py` and updating the `SearchResult` table to include `final_decision` and a link to the resolution record. This will require an Alembic migration.
- **Resolver Agent:** Reviewing and updating the `ScreeningResolutionSchema` Pydantic schema in `schemas.py`, updating the `resolver_prompt`, and constructing the `resolver_chain` using the `RESOLVER_MODEL`.
- **Repository Layer:** Verifying and updating the existing `ScreeningResolutionRepository` in `repositories.py` and potentially updating `SearchResultRepository` (or its service) to handle updates to `final_decision`.
- **Workflow Integration:** Implementing logic in `screen_abstracts.py` to identify disagreements, prepare input for the resolver, invoke the resolver chain, and process/store its output.
- **UI Updates:** Modifying `screen_abstracts.py` to display the `final_decision` and provide access to the resolver's reasoning.

Reference: `docs/prd-resolver.md` for detailed original specifications.

## Local Testability & Command-Line Access

- **Local Development:** Developers must be able to run the Streamlit application via `make run` (which executes `uv run streamlit run src/sr_assistant/app/main.py`) locally. The application connects to the default `postgres` database on a Supabase-hosted instance (this serves as the "prototype" database). The conflict resolution step during abstract screening should be testable by manually creating or identifying conflicts in a screening batch.
- **Command-Line Testing:** Pytest will be used for unit tests and integration tests.
    - Unit tests for new/modified services (e.g., resolver logic within `ScreeningService` or `screen_abstracts.py`), repositories (`ScreeningResolutionRepository`), and agent components (resolver prompt/output parsing).
    - Integration tests for the resolver chain itself must run against the actual LLM API (using approved quotas) to be valid. Integration tests for database interactions related to storing/retrieving resolution data will also be implemented.
    - Tests run via `uv run pytest [options] <test_file_or_directory>` or `make test.unit`/`make test.integration`.
- **Code Formatting & Linting:** Using `uv run ruff format` and `uv run ruff check --fix` or `make format`/`make lint`.
- **Environment Testing:** Local app against Supabase-hosted `postgres` (prototype) DB. Integration tests against Supabase-hosted `sra_integration_test` DB.
- **Testing Prerequisites:**
    - Access to Supabase-hosted PostgreSQL instance (containing the `postgres` prototype and `sra_integration_test` databases).
    - Correct environment variables (`.env`, `.env.local`, `.env.test`) for DB connections and OpenAI API key.
    - Python 3.12 environment via `make install` or `uv sync`.
    - Alembic for database migration (user-managed for the new `ScreeningResolution` table and `SearchResult` modifications).

## Story List

(Drawing inspiration from `docs/prd-resolver.md`)

### Story 2.1: Define and Align Core Pydantic Schemas for Reviews, Screening, Resolution, and Suggestions, and Setup Resolver Data Persistence

- **User Story / Goal:** As a developer, I need all core Pydantic schemas in `src/sr_assistant/core/schemas.py` (for Reviews, Screening, Resolution, and Suggestions) to be correctly defined, documented, and aligned with `docs/data-models.md`, and the database infrastructure for storing screening resolutions to be established, so that data interchange is robust, and resolver outputs can be reliably persisted.
- **Detailed Requirements:**
    1.  **General Pydantic Schema Standards (Apply to all schemas below):**
        *   All Pydantic models MUST inherit from `core.schemas.BaseSchema` (unless a `TypedDict`).
        *   All fields MUST use field docstrings for documentation; `description` parameter in `Field()` MUST NOT be used.
        *   Numerical range constraints for LLM output schemas (`ScreeningResponse`, `ScreeningResolutionSchema`) MUST be in docstrings, not `Field` parameters.
        *   `AwareDatetime` MUST be imported from `pydantic`.
        *   `collections.abc.Mapping` MUST be used instead of `typing.Mapping`.
        *   Relevant linter errors in `schemas.py` MUST be resolved.
    2.  **Refactor/Define `SystematicReview` Schemas in `schemas.py` (as per `docs/data-models.md`):**
        *   `SystematicReviewCreate`: Verify/Align.
        *   `SystematicReviewUpdate`: Refactor to ensure all fields are optional and linter errors are resolved.
        *   `SystematicReviewRead`: Verify/Align.
    3.  **Refactor/Define Screening Schemas (LLM I/O & Service Layer) in `schemas.py` (as per `docs/data-models.md`):**
        *   `ScreeningResponse` (LLM Output): Verify/Align, ensure `confidence_score` range is in docstring, not `Field`.
        *   `ScreeningResult` (Hydrated LLM Output): Verify/Align.
        *   `ScreeningResultCreate` (Service Input): Create new schema.
        *   `ScreeningResultUpdate` (Service Input): Create new schema.
        *   `ScreeningResultRead` (Service Output): Create new schema.
    4.  **Refactor/Define `ScreeningResolution` Schemas in `schemas.py` (as per `docs/data-models.md`):**
        *   `ScreeningResolutionSchema` (LLM Output): Verify/Align. Clarify LLM-populated vs. caller-populated fields in docstrings. Ensure `resolver_confidence_score` range is in docstring, not `Field`.
        *   `ScreeningResolutionCreate` (Service Input): Create new schema for internal service use to construct `models.ScreeningResolution`.
        *   `ScreeningResolutionRead` (Service Output): Create new schema.
    5.  **Refactor/Define Suggestion Agent Schemas in `schemas.py` (as per `docs/data-models.md`):**
        *   `PicosSuggestions` (LLM Output): Verify/Align.
        *   `SuggestionResult` (`TypedDict`): Verify/Align.
    6.  **Database Model & Repository for `ScreeningResolution` (from original Story 2.1 - Expanded for `models.py` alignment):**
        *   Verify and update the existing `ScreeningResolution` SQLModel table in `models.py` as per `docs/prd-resolver.md` (FR6) and `docs/data-models.md`.
        *   **Action:** Add `final_decision: ScreeningDecisionType | None` to the `SearchResult` SQLModel in `models.py`. This field is critical for storing the outcome of the conflict resolution. Ensure it's nullable, indexed, and uses the `ScreeningDecisionType` enum for its database column type.
        *   Verify all other fields and relationships for `SystematicReview`, `SearchResult`, `ScreenAbstractResult`, and `ScreeningResolution` in `models.py` against `docs/data-models.md` and project standards (e.g., `year` as `str | None` in `SearchResult`, FKs, timestamp handling, JSONB typing).
        *   Generate and verify an Alembic migration script for any schema changes made (especially for `SearchResult.final_decision`).
        *   Verify and update the existing `ScreeningResolutionRepository` in `repositories.py`, ensuring it inherits from `BaseRepository[ScreeningResolution]` and has necessary methods (add, get by ID, get by search_result_id).
        *   Ensure `SearchResultRepository` (or `SearchService`) can efficiently update `final_decision` and `resolution_id` on `SearchResult` records.
- **Acceptance Criteria (ACs):**
    *   AC1: All specified Pydantic schemas in `schemas.py` are created or refactored to precisely match definitions, types, optionality, and documentation standards in `docs/data-models.md`.
    *   AC2: All Pydantic schemas inherit from `core.schemas.BaseSchema` (as appropriate) and use field docstrings correctly.
    *   AC3: Constraints for LLM output schema fields (e.g., `confidence_score` ranges) are documented in docstrings and removed from `Field` parameters.
    *   AC4: All linter errors in `schemas.py` are resolved.
    *   AC5: The `SearchResult` SQLModel in `models.py` includes the `final_decision: ScreeningDecisionType | None` field, correctly typed and configured for the database.
    *   AC6: Other SQLModels (`SystematicReview`, `ScreenAbstractResult`, `ScreeningResolution`) are verified to align with `docs/data-models.md` and project standards.
    *   AC7: An Alembic migration script for the `SearchResult.final_decision` addition (and any other necessary SQLModel changes) is generated and successfully applies the schema changes to the database.
    *   AC8: The `ScreeningResolution` table schema in `models.py` is verified/updated and aligns with `docs/prd-resolver.md` (FR6).
    *   AC9: `ScreeningResolutionRepository` is verified/updated with methods to add and retrieve resolution records, and its unit/integration tests pass.
    *   AC10: `SearchResult` records can be efficiently updated with `final_decision` and `resolution_id`, verified by tests.
- **Tasks (Optional Initial Breakdown for SM):**
    *   [ ] Task 2.1.1: **SQLModel Changes (`models.py`):**
        *   Add `final_decision: ScreeningDecisionType | None` field to the `SearchResult` model, ensuring correct DB type, nullability, and indexing.
        *   Verify other SQLModels (`SystematicReview`, `ScreenAbstractResult`, `ScreeningResolution`) for alignment with `docs/data-models.md` and project standards (e.g., `SearchResult.year` as `str | None`, FKs, timestamps, JSONB types). Correct any minor misalignments.
        *   Generate a new Alembic migration script for the `SearchResult.final_decision` addition and any other SQLModel schema changes.
        *   Test the migration script thoroughly.
    *   [ ] Task 2.1.2: **Pydantic Schema Implementation/Refactoring (`schemas.py`):**
        *   Apply general Pydantic standards (BaseSchema inheritance, field docstrings, correct imports like `AwareDatetime` from `pydantic` and `collections.abc.Mapping`) across all relevant schemas.
        *   Refactor/Define `SystematicReviewCreate`, `SystematicReviewUpdate`, `SystematicReviewRead`.
        *   Refactor/Define `ScreeningResponse` (LLM output - check confidence score docstring for range).
        *   Refactor/Define `ScreeningResult` (hydrated LLM output).
        *   Create `ScreeningResultCreate`, `ScreeningResultUpdate`, `ScreeningResultRead`.
        *   Refactor/Define `ScreeningResolutionSchema` (LLM output - check confidence score docstring, clarify caller-populated fields).
        *   Create `ScreeningResolutionCreate`, `ScreeningResolutionRead`.
        *   Refactor/Define `PicosSuggestions`, `SuggestionResult`.
    *   [ ] Task 2.1.3: **Repository Verification/Updates (`repositories.py`):**
        *   Verify/Update `ScreeningResolutionRepository` methods (add, get by ID, get by search_result_id) and associated tests.
        *   Verify/Update `SearchResultRepository` (or `SearchService` interaction points) for efficient updates to `final_decision` and `resolution_id` on `SearchResult` records, and update associated tests.
    *   [ ] Task 2.1.4: **Linting:** Perform a final linter pass on `schemas.py` and `models.py` to resolve any outstanding issues.
- **Dependencies:** `docs/data-models.md` (for schema definitions), `docs/prd-resolver.md` (for `ScreeningResolution` model/table details). Epic 1 (specifically Story 1.5 for `SearchResult` Pydantic schemas) should be complete or its outputs available.

---

### Story 2.2: Define and Implement Resolver Agent (Chain)

- **User Story / Goal:** As a developer, I need a robust LLM agent (chain) that can take disagreement context as input and produce a final, reasoned screening decision, so that conflicts can be resolved automatically.
- **Detailed Requirements:**
  - Review and update the `ScreeningResolutionSchema` Pydantic schema in `src/sr_assistant/core/schemas.py`. Specifically address the `resolver_include` field: evaluate its necessity and utility in light of a definitive `final_decision` field being the primary outcome (as suggested in FR5 of `docs/prd-resolver.md`). If the `resolver_include` field (or a repurposed version) is retained, ensure its name is appropriate (e.g., `contributing_strategies` or similar if it serves a new purpose) and its type is `list[ScreeningStrategyType]` (from `src/sr_assistant/core/types.py`) instead of `list[str]`. The schema must accurately reflect the data intended for populating/deriving the `ScreeningResolution` model.
  - Update/Create the `resolver_prompt` (system and human messages) in `screening_agents.py` to clearly instruct the LLM. The prompt should take all necessary context: `SearchResult` data, `SystematicReview` protocol (PICO, exclusion criteria), and the conservative/comprehensive `ScreenAbstractResult` details. The prompt must align with the (potentially updated) `ScreeningResolutionSchema` for its output structure, especially emphasizing the `final_decision`.
  - Construct the `resolver_chain` using the `RESOLVER_MODEL`, the refined prompt, and structured output parsing (to `ScreeningResolutionSchema`).
  - Implement error handling and retries for the `resolver_chain` similar to existing screening chains.
- **Acceptance Criteria (ACs):**
  - AC1: `ScreeningResolutionSchema` in `schemas.py` is correctly defined/updated; if a field like `resolver_include` is kept, it has an appropriate name and type (`list[ScreeningStrategyType]`).
  - AC2: `resolver_prompt` is implemented and includes placeholders for all required context and aligns with the output schema.
  - AC3: `resolver_chain` is defined, uses the specified model, prompt, and structured output parsing to `ScreeningResolutionSchema`.
  - AC4: The chain correctly processes realistically structured input data (which may contain specific test values for known scenarios) and produces an output that successfully parses into a `ScreeningResolutionSchema` instance.
  - AC5: Integration test for `resolver_chain` (running against the actual LLM API, using approved quotas) confirms it can be invoked with representative real context and returns an output that parses correctly to the `ScreeningResolutionSchema`.
- **Dependencies:** Story 2.1 (for `ScreeningResolutionSchema` to align with `ScreeningResolution` storage).

---

### Story 2.3: Integrate Resolver into Screening Workflow

- **User Story / Goal:** As a developer, I want the resolver agent to be automatically triggered after the initial dual-reviewer screening completes for a batch, so that disagreements are processed without manual intervention.
- **Detailed Requirements:**
  - In `src/sr_assistant/app/pages/screen_abstracts.py`, after a batch of abstracts is screened by `screen_abstracts_chain`:
    - Identify disagreements: Iterate through results, find where conservative and comprehensive decisions differ (e.g., INCLUDE vs. EXCLUDE, as per `docs/prd-resolver.md` FR1). (Note: FR1 of `docs/prd-resolver.md` specifies that cases involving `UNCERTAIN` will not be considered disagreements for automated resolution in this iteration).
    - For each disagreement, prepare the full input context required by `resolver_prompt` (FR3 from `docs/prd-resolver.md`).
    - Invoke the `resolver_chain` (preferably in a batch call if supported and efficient, or iterate) for all identified disagreements in the current screening batch.
    - Process the output (an instance of `ScreeningResolutionSchema`) from the chain.
    - Create `ScreeningResolution` model instances.
    - Use `ScreeningResolutionRepository` to save the new resolution records.
    - Update the corresponding `SearchResult` records with the `final_decision` and `resolution_id`.
  - Implement appropriate status indicators in the UI for the resolution process (e.g., "Resolving conflicts...").
- **Acceptance Criteria (ACs):**
  - AC1: Disagreements between conservative and comprehensive screeners are correctly identified in `screen_abstracts.py` after a batch screening (specifically INCLUDE vs. EXCLUDE or EXCLUDE vs. INCLUDE).
  - AC2: For each identified disagreement, the `resolver_chain` is invoked with the correct input context.
  - AC3: `ScreeningResolution` records are successfully created and stored in the database for each resolved conflict.
  - AC4: `SearchResult` records are updated with the `final_decision` from the resolver and the `resolution_id`.
  - AC5: The process handles cases where there are no disagreements (matching FR1 criteria) in a batch gracefully (i.e., resolver is not called).
  - AC6: Unit tests for the disagreement identification and resolver invocation logic in `screen_abstracts.py` (mocking chain and repositories) pass.
  - AC7: Integration tests verify that the resolver is correctly invoked for `INCLUDE` vs. `EXCLUDE` disagreements (and vice-versa).
  - AC8: Integration tests verify that cases involving `UNCERTAIN` (e.g., `INCLUDE` vs. `UNCERTAIN`, `UNCERTAIN` vs. `UNCERTAIN`) are handled gracefully and *do not* trigger the resolver, adhering to FR1 of `docs/prd-resolver.md`.
- **Dependencies:** Story 2.1, Story 2.2.

---

### Story 2.4: Integrate Resolver Workflow and Update `screen_abstracts.py` UI & Logic

- **User Story / Goal:** As a researcher, I want the abstract screening page (`screen_abstracts.py`) to seamlessly integrate the automated conflict resolution workflow, correctly interact with the refactored `ScreeningService`, and clearly display initial screening decisions, the resolution process status, final resolved decisions, and resolver reasoning, so that the screening process is efficient and transparent.
- **Detailed Requirements:**
    1.  **Service Call Alignment (`ScreeningService`):
        *   Modify `screen_abstracts.py` to use the updated `ScreeningService.add_screening_decision` method. This involves:
            *   Calling it with `search_result_id`, `screening_strategy`, and a correctly populated `schemas.ScreeningResultCreate` object (which is derived from the `schemas.ScreeningResult` obtained from the `screen_abstracts_batch` LLM chain output).
        *   Ensure calls to methods like `ScreeningService.get_screening_decisions_for_search_result` or `get_screening_decisions_for_review` use the new signatures (no `session` parameter).
        *   Remove any direct session management from `screen_abstracts.py` related to these service calls.
    2.  **Conflict Resolution Triggering:**
        *   After a batch of abstracts is screened (i.e., after `screen_abstracts_batch` function completes and initial `ScreenAbstractResult` records are created via `ScreeningService.add_screening_decision`), `screen_abstracts.py` MUST call `ScreeningService.resolve_screening_conflicts_for_batch`.
        *   Pass the current `models.SystematicReview` object and the list of `SearchResult` IDs from the just-screened batch to this service method.
    3.  **UI Updates for Resolver Workflow (Incorporates original Story 2.4 scope):
        *   **Status Indication:** Display a status indicator in the UI while the `resolve_screening_conflicts_for_batch` method is processing (e.g., "Resolving conflicts...").
        *   **Display Final Decisions:** The results display (e.g., Streamlit DataFrame or custom layout) in `screen_abstracts.py` MUST prioritize displaying `SearchResult.final_decision` if it is populated. If `final_decision` is null, the UI should clearly show the original individual decisions from the conservative and comprehensive screeners and indicate an 'unresolved' or 'pending resolution' status.
        *   **Resolution Indicator:** If `SearchResult.resolution_id` is present (indicating a conflict was processed by the resolver), a clear visual cue (e.g., an icon, tag like "Resolved by AI") MUST be displayed next to the decision.
        *   **Resolver Reasoning Display:** Implement a mechanism (e.g., a tooltip on the resolution indicator, an expandable section, a modal popup when clicking the indicator) to allow users to view the `resolver_reasoning` (and potentially `resolver_confidence_score`, `resolver_model_name`) from the corresponding `ScreeningResolution` record.
        *   **Summary Metrics:** Display a summary for the current screening batch, such as 'X conflicts identified for resolution' and/or 'Y conflicts automatically resolved' (Story 2.4 AC5).
    4.  **Data Handling and State Management:**
        *   Ensure `screen_abstracts.py` correctly fetches and re-fetches `SearchResult` data as needed to reflect updates made by the `ScreeningService` (e.g., populated `final_decision`, `resolution_id`).
        *   Manage Streamlit session state appropriately to handle the multi-step process (initial screening, then resolution) and UI updates.
    5.  **Error Handling:** Adapt error handling in `screen_abstracts.py` to manage and display errors from the refactored `ScreeningService` calls, including the conflict resolution step.
- **Acceptance Criteria (ACs):**
    *   AC1: `screen_abstracts.py` correctly calls the refactored `ScreeningService.add_screening_decision` with appropriate `search_result_id`, `screening_strategy`, and `schemas.ScreeningResultCreate` data.
    *   AC2: `screen_abstracts.py` calls `ScreeningService.resolve_screening_conflicts_for_batch` after initial batch screening is complete.
    *   AC3: The UI correctly displays `SearchResult.final_decision` when available, otherwise shows original decisions and conflict status.
    *   AC4: A visual cue indicates decisions that have been automatically resolved by the resolver.
    *   AC5: Users can easily view the `resolver_reasoning` for resolved items.
    *   AC6: UI provides status updates during the conflict resolution phase.
    *   AC7: Summary metrics for conflict resolution in the batch are displayed.
    *   AC8: Session state and data refresh logic in `screen_abstracts.py` correctly handles the updated data post-resolution.
    *   AC9: Appropriate error handling is in place for all service interactions.
- **Tasks (Optional Initial Breakdown for SM):
    *   [ ] Task 2.4.1: Refactor `screen_abstracts.py` to use the new `ScreeningService.add_screening_decision` signature and logic. This includes ensuring `schemas.ScreeningResultCreate` is correctly prepared from `screen_abstracts_batch` output.
    *   [ ] Task 2.4.2: Integrate the call to `ScreeningService.resolve_screening_conflicts_for_batch` in `screen_abstracts.py` flow after initial screening of a batch.
    *   [ ] Task 2.4.3: Implement UI elements in `screen_abstracts.py` to display `final_decision` preferentially.
    *   [ ] Task 2.4.4: Implement UI visual indicators for resolved decisions.
    *   [ ] Task 2.4.5: Implement UI mechanism to display resolver reasoning.
    *   [ ] Task 2.4.6: Implement UI status indicators for the resolution process.
    *   [ ] Task 2.4.7: Implement UI summary metrics for conflict resolution.
    *   [ ] Task 2.4.8: Review and update Streamlit session state management and data refresh logic in `screen_abstracts.py`.
    *   [ ] Task 2.4.9: Update error handling for `ScreeningService` calls.
- **Dependencies:** Story 2.1 (for `ScreeningResolution` model/repo and Pydantic schemas), Story 2.2 (for `resolver_chain`), Story 2.5 (for refactored `ScreeningService` and new resolver methods).

---

### Story 2.5: Refactor and Implement `ScreeningService` for API Alignment and Resolver Logic

- **User Story / Goal:** As a developer, I need the `ScreeningService` in `services.py` to be refactored to align with `docs/api-reference.md` (correct method signatures, internal session management, proper Pydantic schema usage) and to implement all new methods required for the conflict resolution workflow, so that screening operations and conflict resolution are robust, maintainable, and correctly implemented.
- **Detailed Requirements:**
    1.  **Refactor Existing `ScreeningService` Methods:**
        *   **Session Management:** All public methods in `ScreeningService` (e.g., `add_screening_decision`, `update_screening_decision`, `get_screening_decision_by_id`, etc.) MUST manage their own database sessions internally. The `session: Session | None` parameter MUST be removed from all public method signatures.
        *   **`add_screening_decision`:**
            *   Signature MUST be `(self, search_result_id: uuid.UUID, screening_strategy: ScreeningStrategyType, screening_data: schemas.ScreeningResultCreate) -> models.ScreenAbstractResult`.
            *   Logic MUST correctly create `models.ScreenAbstractResult` from `screening_data` and then update the parent `models.SearchResult` (identified by `search_result_id`) by setting its `conservative_result_id` or `comprehensive_result_id` field based on `screening_strategy`.
        *   **`update_screening_decision`:**
            *   Signature MUST be `(self, screening_result_id: uuid.UUID, screening_update_data: schemas.ScreeningResultUpdate) -> models.ScreenAbstractResult`.
            *   Logic MUST fetch the `ScreenAbstractResult` by `screening_result_id` and apply updates from `screening_update_data`.
        *   **Method Consolidation/Removal:** The existing `add_screening_result` and `add_or_update_screening_result` methods in `services.py` should be refactored into or replaced by the correctly defined `add_screening_decision` and `update_screening_decision` methods.
        *   **Linter Error (Instantiation):** Fix the linter error in the current `add_or_update_screening_result` related to `models.ScreenAbstractResult(id=result_id, **update_data)`. Ensure explicit field mapping from the Pydantic input schema (`ScreeningResultCreate`) when instantiating `models.ScreenAbstractResult`.
        *   **Linter Error (Unused Variable):** Address the unused `e` variable in the `except` block of the current `get_screening_result_by_strategy`.
    2.  **Implement New Resolver-Related Methods in `ScreeningService` (as per `docs/api-reference.md`):**
        *   `identify_disagreements(self, review_id: uuid.UUID, search_result_ids: list[uuid.UUID]) -> list[models.SearchResult]`
        *   `prepare_resolver_inputs(self, review: models.SystematicReview, search_results_with_disagreements: list[models.SearchResult]) -> list[dict[str, Any]]`
        *   `invoke_resolver_agent_batch(self, resolver_prompt_variable_inputs: list[dict[str, Any]]) -> list[schemas.ScreeningResolutionSchema]` (This method calls the LLM chain and does not directly manage DB sessions for the call itself, but subsequent storage does).
        *   `store_resolution_results(self, review_id: uuid.UUID, search_result_id_to_resolution_data: dict[uuid.UUID, schemas.ScreeningResolutionSchema]) -> list[models.ScreeningResolution]` (This method handles DB writes and manages its session).
        *   `resolve_screening_conflicts_for_batch(self, review: models.SystematicReview, search_result_ids_in_batch: list[uuid.UUID]) -> None` (Orchestration method, manages sessions for its sequence of DB operations).
    3.  **General Standards:**
        *   All new and refactored methods must manage sessions internally as appropriate.
        *   Use defined Pydantic schemas (`ScreeningResultCreate`, `ScreeningResultUpdate`, input dicts and `ScreeningResolutionSchema` for resolver) for data interchange as per `docs/api-reference.md` and `docs/data-models.md`.
        *   Implement robust error handling for all operations.
        *   Resolve all linter errors for `ScreeningService` in `services.py`.
- **Acceptance Criteria (ACs):**
    *   AC1: All public methods of `ScreeningService` manage their own DB sessions and no longer accept `session` parameters (where DB interaction occurs).
    *   AC2: `add_screening_decision` method is implemented with the correct signature and logic, including linking to `SearchResult`.
    *   AC3: `update_screening_decision` method is implemented with the correct signature and logic.
    *   AC4: Existing methods `add_screening_result` and `add_or_update_screening_result` are correctly refactored/replaced.
    *   AC5: Linter errors in existing `ScreeningService` methods (e.g., `**update_data` instantiation, unused `e`) are resolved.
    *   AC6: All new resolver-related methods (`identify_disagreements`, `prepare_resolver_inputs`, `invoke_resolver_agent_batch`, `store_resolution_results`, `resolve_screening_conflicts_for_batch`) are implemented in `ScreeningService` as per `docs/api-reference.md`.
    *   AC7: Implemented methods correctly use the Pydantic schemas defined in `docs/data-models.md` and `schemas.py`.
    *   AC8: All linter errors related to `ScreeningService` in `services.py` are resolved.
    *   AC9: Unit tests for all `ScreeningService` methods (existing refactored and new) achieve >80% coverage, mocking repositories and LLM calls where appropriate.
    *   AC10: Relevant integration tests for the resolver workflow involving `ScreeningService` pass.
- **Tasks (Optional Initial Breakdown for SM):**
    *   [ ] Task 2.5.1: Refactor existing `ScreeningService` methods (`add_screening_result`, `add_or_update_screening_result`, `get_screening_result_by_strategy`, etc.) for internal session management, correct Pydantic schema inputs (using `ScreeningResultCreate`/`Update`), and fix linter errors.
    *   [ ] Task 2.5.2: Implement `identify_disagreements` method in `ScreeningService`.
    *   [ ] Task 2.5.3: Implement `prepare_resolver_inputs` method in `ScreeningService`.
    *   [ ] Task 2.5.4: Implement `invoke_resolver_agent_batch` method in `ScreeningService` (to call resolver LLM chain).
    *   [ ] Task 2.5.5: Implement `store_resolution_results` method in `ScreeningService` (to save `ScreeningResolution` and update `SearchResult`).
    *   [ ] Task 2.5.6: Implement `resolve_screening_conflicts_for_batch` orchestration method in `ScreeningService`.
    *   [ ] Task 2.5.7: Write/update unit tests for all `ScreeningService` methods.
    *   [ ] Task 2.5.8: Run linter and fix any remaining errors for `ScreeningService` in `services.py`.
- **Dependencies:** `docs/api-reference.md`, `docs/data-models.md`, relevant Pydantic schemas created/updated via Story 2.1, `ScreeningResolutionRepository` (from Story 2.1), `resolver_chain` (from Story 2.2).

---

## Change Log

| Change          | Date       | Version | Description             | Author               |
|-----------------|------------|---------|-------------------------|----------------------|
| Initial Draft   | 2025-05-09 | 0.1     | First draft of Epic 2   | Product Manager Agent |
| Story 2.1 Update| 2025-05-10 | 0.2     | Expanded Story 2.1 to include all core Pydantic schema definitions, SQLModel updates, and resolver DB setup. | Architect Agent      |
| Added Story 2.5 | 2025-05-10 | 0.3     | Added Story 2.5 for `ScreeningService` refactoring and resolver method implementation. | Architect Agent      |
| Story 2.4 Update| 2025-05-10 | 0.4     | Expanded Story 2.4 to cover `screen_abstracts.py` service call alignment, resolver integration, and all UI logic. | Architect Agent      | 