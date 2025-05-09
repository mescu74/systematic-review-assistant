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

### Story 2.1: Data Model and Repository Setup for Resolver

- **User Story / Goal:** As a developer, I need the existing database schema and repository layer to be verified and updated to correctly store and handle resolution decisions, so that resolver outputs can be reliably persisted and linked to search results.
- **Detailed Requirements:**
  - Verify and update the existing `ScreeningResolution` SQLModel table in `models.py` as per `docs/prd-resolver.md` (FR6), ensuring all fields like `search_result_id`, `resolver_decision`, `resolver_reasoning`, `resolver_confidence_score`, `resolver_model_name`, `trace_id` are correctly defined and typed. Review if any field analogous to `resolver_include` from `ScreeningResolutionSchema` exists or needs to be represented and ensure consistency or appropriate mapping if the schema is adjusted.
  - Verify and update the `SearchResult` model (`models.py`) to ensure it correctly includes `final_decision: ScreeningDecisionType | None` and `resolution_id: uuid.UUID | None` (FK to `ScreeningResolution`).
  - Generate and verify an Alembic migration script if any schema changes are made to `ScreeningResolution` or `SearchResult`.
  - Verify and update the existing `ScreeningResolutionRepository` in `repositories.py`, ensuring it inherits from `BaseRepository[ScreeningResolution]` and has necessary methods.
  - Ensure `add` (or `add_all`) and retrieval methods (e.g., `get_by_search_result_id`) in `ScreeningResolutionRepository` are correctly implemented.
  - Ensure `SearchResultRepository` (or `SearchService`) can efficiently update `final_decision` and `resolution_id` on `SearchResult` records.
- **Acceptance Criteria (ACs):**
  - AC1: The `ScreeningResolution` table schema in `models.py` is verified/updated and aligns with the requirements in `docs/prd-resolver.md` (FR6).
  - AC2: The `SearchResult` model in `models.py` correctly includes `final_decision` and `resolution_id` fields.
  - AC3: If schema changes were necessary, the Alembic migration script successfully applies them.
  - AC4: `ScreeningResolutionRepository` is verified/updated with methods to add and retrieve resolution records.
  - AC5: `SearchResult` records can be efficiently updated with `final_decision` and `resolution_id`.
  - AC6: Unit tests for `ScreeningResolutionRepository` methods pass.
  - AC7: Integration tests for `ScreeningResolutionRepository` methods (e.g., add, get by ID, get by search_result_id) and `SearchResult` updates (for resolver fields like `final_decision` and `resolution_id`) pass, ensuring correct CRUD operations and relationships.
- **Dependencies:** None internal to this epic initially. Epic 1 (Search/Service Layer Stability) should be largely complete.

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

### Story 2.4: Update UI to Display Resolved Decisions

- **User Story / Goal:** As a researcher, I want to clearly see the final resolved screening decision and be able to access the resolver's reasoning, so that I understand the outcome of the automated resolution.
- **Detailed Requirements:**
  - In `src/sr_assistant/app/pages/screen_abstracts.py`, modify the results display (Streamlit DataFrame):
    - Prioritize displaying `SearchResult.final_decision` if it's populated. If not, show the original conservative/comprehensive decisions (or indicate conflict).
    - Add a visual indicator (e.g., icon, tag like "Resolved") if `SearchResult.resolution_id` is present.
    - Provide a mechanism to view the `resolver_reasoning` from the linked `ScreeningResolution` record (e.g., tooltip, expandable section, modal popup when clicking the resolution indicator).
  - *Future Consideration: In a subsequent epic, this DataFrame will need to support manual overrides of the `final_decision`.*
- **Acceptance Criteria (ACs):**
  - AC1: The UI correctly displays the `final_decision` from `SearchResult` when available.
  - AC2: A visual cue indicates which decisions have been automatically resolved.
  - AC3: Users can easily view the `resolver_reasoning` for resolved items.
  - AC4: If an abstract does not have a `final_decision` (either because there was no initial disagreement meeting the resolution criteria, or the resolver did not process it), the UI clearly displays the original individual decisions from the conservative and comprehensive screeners and indicates the 'unresolved' status.
  - AC5: The UI displays a summary metric for the current screening batch, such as 'X conflicts sent for resolution' or 'X conflicts resolved by automation'.
- **Dependencies:** Story 2.1 (for data availability), Story 2.3 (for data to be populated).

---

## Change Log

| Change          | Date       | Version | Description             | Author               |
|-----------------|------------|---------|-------------------------|----------------------|
| Initial Draft   | 2025-05-09 | 0.1     | First draft of Epic 2   | Product Manager Agent | 