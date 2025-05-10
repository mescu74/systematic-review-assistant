# Epic 1: Search and Service Layer Stabilization

**Goal:** Fully refactor and stabilize `search.py`, `services.py`, and `repositories.py` to correctly implement and utilize the `SearchResult` model and the new service layer architecture for PubMed search operations, ensuring all related functionalities are operational and error-free. This directly supports the PRD objective to "Restore the MPH SR Prototype to a stable, robust, and reliable state" and "Correctly implement a generic search infrastructure, initially focused on PubMed."

**Deployability:** This epic establishes the foundational stability for the search and data access layers. Its successful completion is critical before other epics (Resolver Implementation, Comprehensive Testing) can be confidently built upon it. It ensures the core data retrieval and storage mechanism for PubMed searches is reliable.

## Epic-Specific Technical Context

This epic focuses on refactoring and stabilizing existing components. Key technical considerations include:
- Adherence to the intended service layer pattern: `SearchService` in `services.py` should manage business logic and database sessions, calling `SearchResultRepository` in `repositories.py` for data I/O.
- Consistent use of the `SearchResult` model (`models.py`) for all PubMed search-related data, ensuring fields like `source_db` (set to 'PubMed') and `source_id` (the PubMed ID) are correctly populated and used.
- Correct management of SQLAlchemy sessions within the `SearchService`.
- Ensuring `SearchResultRepository` has, and correctly implements, methods required by `SearchService` (e.g., methods for adding search results for a review, deleting results by review ID, retrieving results by review ID, checking for existing results by `source_id` for a given review).

## Local Testability & Command-Line Access

- **Local Development:** Developers must be able to run the Streamlit application (`streamlit run src/sr_assistant/app/Home.py` or `uv run streamlit run src/sr_assistant/app/main.py` as per `Makefile`) locally. The application connects to the default `postgres` database on a Supabase-hosted instance (this serves as the development database).
- **Command-Line Testing:** Pytest will be used for running unit tests and integration tests.
    - Unit tests (e.g., for `services.py`, `repositories.py`) must mock database interactions as there is no local PostgreSQL instance for unit testing. 
    - Integration tests run against a dedicated Supabase-hosted integration test database (`sra_integration_test`).
    - Developers will run tests via `uv run pytest [options] <test_file_or_directory>` command in the terminal (e.g., `uv run pytest tests/unit`, `uv run pytest -m integration tests/integration`). Alternatively, `make test.unit` and `make test.integration` can be used.
- **Code Formatting & Linting:** Code formatting is done using `uv run ruff format <path/to/file_or_directory>`. Linting errors can often be auto-fixed using `uv run ruff check --fix <path/to/file_or_directory>`. The `make format` and `make lint` (or `make ruff.fix`) targets can also be used.
- **Environment Testing:** The local Streamlit application runs against the Supabase-hosted `postgres` (development) database. Integration tests run against the Supabase-hosted `sra_integration_test` database, as configured in `.env.test`.
- **Testing Prerequisites:**
    - Access to the Supabase-hosted PostgreSQL instance containing the `postgres` (development) and `sra_integration_test` databases.
    - Correct environment variable setup (`.env`, `.env.local`, `.env.test`) for database connections and API keys (e.g., NCBI API key for PubMed searches if directly used by the service).
    - Python 3.12 environment with all project dependencies installed via `make install` (which uses `uv sync`) or directly with `uv sync`, based on `pyproject.toml`.
    - Alembic for managing any necessary schema adjustments (user-managed).

## Story List

### Story 1.1: Ensure `SearchResult` Model Consistency and Correct Service Usage in `search.py`

- **User Story / Goal:** As a developer, I want `search.py` to correctly interact with the refactored `SearchService` and consistently use the `schemas.SearchResultRead` Pydantic schema (derived from the `models.SearchResult` model) for displaying and handling PubMed search results, so that data representation is accurate, service layer encapsulation is respected, and the UI aligns with the intended generic search architecture.
- **Detailed Requirements:**
  1.  **Service Interaction:**
      *   All PubMed search operations initiated from `search.py` MUST call the `SearchService.search_pubmed_and_store_results` method.
      *   `search.py` MUST NOT pass any `session` objects to `SearchService` methods.
      *   Any logic for directly calling PubMed APIs or managing raw API records within `search.py` MUST be removed (this is now `SearchService`'s responsibility).
  2.  **Data Handling & Display:**
      *   `search.py` MUST expect a sequence of `schemas.SearchResultRead` objects (or equivalent Pydantic models that the service layer would return, based on `models.SearchResult`) from `SearchService.search_pubmed_and_store_results`.
      *   All UI components in `search.py` (e.g., Streamlit tables, detail display elements) that show search result information MUST be updated to bind to the fields of `schemas.SearchResultRead` (e.g., `source_id`, `source_db`, `title`, `abstract`, `year` as string, etc.).
      *   Ensure any actions taken on search results (e.g., selection for screening, displaying details) correctly reference and use data from these `schemas.SearchResultRead` Pydantic objects.
  3.  **Error Handling:**
      *   Update error handling in `search.py` to appropriately manage and display errors that might be raised by the refactored `SearchService`.
- **Acceptance Criteria (ACs):**
  - AC1: `search.py` calls `SearchService.search_pubmed_and_store_results` for PubMed searches and does not pass session objects.
  - AC2: PubMed search results displayed in the `search.py` UI correctly use fields from `schemas.SearchResultRead` Pydantic objects.
  - AC3: All direct PubMed API interaction logic is removed from `search.py`.
  - AC4: All references to PubMed-specific identifiers (like `pmid`) in `search.py` are handled via `SearchResultRead.source_id` (where `source_db` is 'PubMed').
  - AC5: Error handling for service calls is implemented in `search.py`.
- **Tasks (Optional Initial Breakdown for SM):**
  - [ ] Task 1.1.1: Analyze `search.py` for current PubMed search initiation logic and data handling.
  - [ ] Task 1.1.2: Refactor search initiation in `search.py` to call `SearchService.search_pubmed_and_store_results` and remove any direct PubMed API calls or session management from the page.
  - [ ] Task 1.1.3: Update data handling logic in `search.py` to work with `schemas.SearchResultRead` objects returned by the service.
  - [ ] Task 1.1.4: Update Streamlit UI components (tables, text displays) in `search.py` to bind to `schemas.SearchResultRead` attributes correctly.
  - [ ] Task 1.1.5: Implement or update error display logic in `search.py` for `SearchService` interactions.
- **Dependencies:** Relies on `SearchService` refactoring (Story 1.2) and `schemas.SearchResultRead` definition (Story 1.5).

---

### Story 1.2: Refactor `SearchService` for PubMed Logic & Session Management

- **User Story / Goal:** As a developer, I want the `SearchService` in `services.py` to correctly manage all business logic for PubMed searches, including API interaction and database session lifecycle, using architecturally sound interfaces, so that search operations are reliable, maintainable, and adhere to the service layer pattern defined in `docs/api-reference.md`.
- **Detailed Requirements:**
    1.  **Encapsulate PubMed API Interaction:** `SearchService` must implement the `search_pubmed_and_store_results` method. This method will:
        *   Accept `review_id`, `query`, and `max_results`.
        *   Internally handle all aspects of querying the PubMed API (e.g., using BioPython Entrez or a similar HTTP client mechanism).
        *   Fetch raw results from PubMed.
        *   Utilize the internal `_map_pubmed_to_search_result` method to transform raw PubMed records into `models.SearchResult` instances.
        *   Ensure `SearchResult` instances are correctly populated (e.g., `source_db='PubMed'`, `source_id=PMID`, `year` as `str | None`).
        *   Store these `SearchResult` instances using `SearchResultRepository.add_all` (or similar batch method).
    2.  **Session Management:**
        *   All public methods in `SearchService` (including the new `search_pubmed_and_store_results`, `get_search_results_by_review_id`, `get_search_result_by_source_details`, `update_search_result`, `delete_search_result`) MUST manage their own database sessions internally (e.g., using `with self.session_factory() as session:`). The `session: Session | None` parameter MUST be removed from all public method signatures.
        *   The internally managed session MUST be passed to any `SearchResultRepository` methods called.
    3.  **Update Method Signature and Logic:**
        *   The `update_search_result` method signature MUST be changed from `(self, result_update: models.SearchResult, ...)` to `(self, result_id: uuid.UUID, update_data: schemas.SearchResultUpdate, ...)` as per `docs/api-reference.md`. Internally, it will fetch the `SearchResult` by `result_id`, apply updates from `update_data`, and save.
    4.  **Method Removal:**
        *   The existing `add_api_search_results` method MUST be removed. Its functionality is superseded by the properly encapsulated `search_pubmed_and_store_results`.
    5.  **Internal Mapping Methods:**
        *   Ensure `_map_pubmed_to_search_result` and `_map_scopus_to_search_result` have correct type hints (`api_record: dict[str, Any]`).
        *   Verify mapping logic aligns with `models.SearchResult` (e.g., `year` as `str | None`). The linter error regarding `year` being `int` in the mapper vs `str` in the model needs to be resolved by ensuring the mapper provides a string.
    6.  **Error Handling:** Implement robust error handling for PubMed API interactions and database operations within the service methods.
    7.  **Linter Errors:** Resolve all linter errors within the `SearchService` code block in `services.py` after refactoring.
    8.  **`ReviewService` Session Management:**
        *   All public methods in `ReviewService` (`create_review`, `get_review`, `get_all_reviews`, `update_review`, `delete_review`) MUST manage their own database sessions internally. The `session: Session | None` parameter MUST be removed from their public signatures.
        *   The internally managed session MUST be passed to any `SystematicReviewRepository` methods called.
- **Acceptance Criteria (ACs):**
    *   AC1: `SearchService.search_pubmed_and_store_results` is implemented, encapsulates PubMed search, mapping, and storage logic, and returns a sequence of `models.SearchResult`.
    *   AC2: All public methods in `SearchService` manage sessions internally and no longer accept a `session` parameter.
    *   AC3: `SearchService.update_search_result` signature and implementation are updated to use `result_id` and `schemas.SearchResultUpdate`.
    *   AC4: `SearchService.add_api_search_results` method is removed.
    *   AC5: `SearchService` uses `SearchResultRepository` for all database interactions related to `SearchResult` objects, passing an internally managed session.
    *   AC6: `_map_pubmed_to_search_result` correctly maps PubMed data to `models.SearchResult`, including `year` as `str | None`.
    *   AC7: Linter errors in `services.py` related to `SearchService` are resolved.
    *   AC8: Unit test coverage for `SearchService` (mocking repository and PubMed API calls) achieves >80%.
    *   AC9: Relevant integration tests for PubMed search functionality via `SearchService` pass.
    *   AC10: All public methods in `ReviewService` manage sessions internally and no longer accept a `session` parameter.
    *   AC11: `ReviewService` uses `SystematicReviewRepository` for all database interactions related to `SystematicReview` objects, passing an internally managed session.
- **Tasks (Optional Initial Breakdown for SM):**
    *   [ ] Task 1.2.1: Remove `session` parameter from all public methods of `SearchService` and implement internal session management logic for each.
    *   [ ] Task 1.2.2: Implement the `search_pubmed_and_store_results` method in `SearchService`:
        *   Integrate PubMed API call logic.
        *   Use `_map_pubmed_to_search_result` for mapping.
        *   Use `SearchResultRepository` for storing results.
    *   [ ] Task 1.2.3: Refactor `update_search_result` method in `SearchService`:
        *   Change signature to accept `result_id` and `schemas.SearchResultUpdate`.
        *   Implement logic to fetch, update, and save.
    *   [ ] Task 1.2.4: Remove the `add_api_search_results` method from `SearchService`.
    *   [ ] Task 1.2.5: Review and correct internal mapping methods (`_map_pubmed_to_search_result`, `_map_scopus_to_search_result`) for type hints and `year` mapping logic.
    *   [ ] Task 1.2.6: Implement error handling for API and DB operations.
    *   [ ] Task 1.2.7: Write/update unit tests for `SearchService` methods.
    *   [ ] Task 1.2.8: Run linter and fix errors for `SearchService` in `services.py`.
    *   [ ] Task 1.2.9: Refactor `ReviewService` public methods in `services.py`:
        *   Remove `session` parameter from `create_review`, `get_review`, `get_all_reviews`, `update_review`, `delete_review`.
        *   Implement internal session management logic for each of these methods.
        *   Ensure repository calls within these methods use the internally managed session.
    *   [ ] Task 1.2.10: Write/update unit tests for `ReviewService` methods to reflect session encapsulation (mocking the session factory and repository).
    *   [ ] Task 1.2.11: Verify `src/sr_assistant/app/pages/protocol.py` interactions with `ReviewService`:
        *   Confirm `build_review_model_from_pico` correctly prepares data that aligns with the input expected by `ReviewService.create_review` (which takes `schemas.SystematicReviewCreate`).
        *   Confirm the `persist_review` function in `protocol.py` correctly uses `ReviewService.update_review` (passing `schemas.SystematicReviewUpdate`) or `ReviewService.create_review` (passing `schemas.SystematicReviewCreate`) for saving changes, instead of direct model manipulation and commit after fetching from the service.
        *   Ensure `protocol.py` does not pass session objects to `ReviewService` methods.
- **Dependencies:** `SearchResultRepository` (Story 1.3), `schemas.SearchResultUpdate` (defined via Story 1.5), `SystematicReviewRepository` (implicitly, for `ReviewService`), `schemas.SystematicReviewCreate`, `schemas.SystematicReviewUpdate` (defined via Epic 2, Story 2.1).

---

### Story 1.3: Repair and Equip `SearchResultRepository`

- **User Story / Goal:** As a developer, I want `SearchResultRepository` in `repositories.py` to have all necessary methods for `SearchResult` CRUD operations as required by `SearchService`, and for these methods to be correctly implemented and error-free, so that data persistence for search results is reliable.
- **Detailed Requirements:**
  - Ensure `SearchResultRepository` includes methods for:
    - Adding a list of `SearchResult` objects for a given review (`add_all_for_review_id`).
    - Retrieving all `SearchResult` objects for a given review (`get_by_review_id`).
    - Deleting all `SearchResult` objects for a given review (`delete_by_review_id`).
    - Checking if a `SearchResult` with a specific `source_id` and `source_db` already exists for a review (`exists_by_source_id`).
  - All repository methods must accept a `Session` object from the calling service and use it for database operations.
  - Repository methods should not manage transactions (commit/rollback) themselves; this is the responsibility of the service layer.
  - Resolve linter errors in `repositories.py` related to `SearchResultRepository` and its base class (e.g., using `noqa` for known `Model.id` issues if deemed non-critical for functionality).
- **Acceptance Criteria (ACs):**
  - AC1: `SearchResultRepository` contains the defined methods (`add_all_for_review_id`, `get_by_review_id`, `delete_by_review_id`, `exists_by_source_id`).
  - AC2: These methods are correctly implemented, perform the intended database operations using the provided session, and handle potential `SQLAlchemyError` exceptions by raising `RepositoryError` or its sub-classes.
  - AC3: Linter errors specifically within `SearchResultRepository` and its interactions with the `BaseRepository` concerning `SearchResult` are resolved or appropriately ignored.
  - AC4: Unit tests for `SearchResultRepository` methods pass.
  - AC5: Each implemented repository method in `SearchResultRepository` has corresponding integration tests that pass.
- **Tasks (Optional Initial Breakdown):**
  - [ ] Define the signatures for the required repository methods.
  - [ ] Implement each method using SQLAlchemy operations via the passed session.
  - [ ] Add error handling and logging.
  - [ ] Write/repair unit tests for these methods.
- **Dependencies:** None for method definition; `SearchService` (Story 1.2) will be a consumer.

---

### Story 1.4: End-to-End PubMed Search Workflow Stabilization

- **User Story / Goal:** As a user, I want to perform a PubMed search using the `search.py` interface and have the results correctly fetched, stored, and displayed without errors, so that I can reliably gather articles for my systematic review.
- **Detailed Requirements:**
  - The entire flow from entering a search query in `search.py`, triggering `SearchService`, fetching from PubMed, transforming to `SearchResult` objects, storing via `SearchResultRepository`, and displaying in `search.py` must be functional.
  - All linter errors in `search.py`, `services.py` (related to search), and `repositories.py` (related to `SearchResultRepository`) that impact this workflow must be resolved.
  - The system should prevent duplicate storage of the same PubMed article for the same review.
- **Acceptance Criteria (ACs):**
  - AC1: User can successfully execute a PubMed search from the `search.py` UI.
  - AC2: Search results are displayed correctly in the UI using `SearchResult` data.
  - AC3: `SearchResult` instances are correctly persisted in the database via the service and repository layers.
  - AC4: No Python errors or exceptions occur during the PubMed search and result display workflow.
  - AC5: Attempting to store the same PubMed article for the same review multiple times does not create duplicate entries (verified via `exists_by_source_id` or similar logic in service).
- **Tasks (Optional Initial Breakdown):**
  - [ ] Perform manual end-to-end testing of the PubMed search workflow.
  - [ ] Debug and fix any errors encountered at any stage of the workflow.
  - [ ] Ensure linter errors in the relevant files are addressed.
  - [ ] Implement or verify duplicate prevention logic in `SearchService`.
- **Dependencies:** Story 1.1, Story 1.2, Story 1.3.

---

### Story 1.5: Define and Align `SearchResult` Pydantic Schemas

- **User Story / Goal:** As a developer, I need `SearchResultRead` to be corrected and `SearchResultUpdate` to be defined in `src/sr_assistant/core/schemas.py` according to `docs/data-models.md`, so that data interchange for `SearchResult` entities is robust, well-documented, and aligns with architectural standards.
- **Detailed Requirements:**
    1.  **Refactor `SearchResultRead` in `src/sr_assistant/core/schemas.py`:**
        *   Ensure it inherits from `core.schemas.BaseSchema`.
        *   Verify all fields match the definition in `docs/data-models.md#SearchResultRead`, including:
            *   `year: str | None`
            *   `raw_data: collections.abc.Mapping[str, JsonValue]`
            *   `source_metadata: collections.abc.Mapping[str, JsonValue]`
            *   `final_decision: ScreeningDecisionType | None = None`
            *   `resolution_id: uuid.UUID | None = None`
            *   `created_at: AwareDatetime | None` (imported from `pydantic`)
            *   `updated_at: AwareDatetime | None` (imported from `pydantic`)
        *   Ensure all fields use field docstrings for documentation (no `description` in `Field()`).
    2.  **Create `SearchResultUpdate` in `src/sr_assistant/core/schemas.py**:
        *   Define this new schema as specified in `docs/data-models.md#SearchResultUpdate`.
        *   It must inherit from `core.schemas.BaseSchema`.
        *   All its fields must be optional (e.g., `doi: str | None = None`).
        *   Ensure all fields use field docstrings.
    3.  **General Standards:**
        *   Apply all "General Standards" and "Pydantic Field Documentation Standard" from `docs/data-models.md#Guiding-Principles` to these schemas.
        *   Resolve any linter errors specifically related to these two schemas after changes.
- **Acceptance Criteria (ACs):**
    *   AC1: `SearchResultRead` in `schemas.py` is updated to precisely match the field definitions, types, optionality, and field docstring standards specified in `docs/data-models.md`.
    *   AC2: `SearchResultUpdate` is created in `schemas.py` and precisely matches the field definitions, types, optionality, and field docstring standards specified in `docs/data-models.md`.
    *   AC3: Both `SearchResultRead` and `SearchResultUpdate` inherit from `core.schemas.BaseSchema`.
    *   AC4: `AwareDatetime` is correctly imported from `pydantic` for these schemas.
    *   AC5: `collections.abc.Mapping` is used where appropriate.
    *   AC6: All linter errors related to `SearchResultRead` and `SearchResultUpdate` in `schemas.py` are resolved.
- **Tasks (Optional Initial Breakdown for SM):**
    *   [ ] Task 1.5.1: Modify `schemas.py` to update the `SearchResultRead` schema definition.
        *   Verify/correct field types (year, raw_data, source_metadata, datetimes).
        *   Add missing fields (`final_decision`, `resolution_id`).
        *   Implement field docstrings for all fields.
        *   Ensure `BaseSchema` inheritance.
        *   Correct imports.
    *   [ ] Task 1.5.2: Modify `schemas.py` to add the new `SearchResultUpdate` schema definition.
        *   Implement all fields as optional.
        *   Implement field docstrings for all fields.
        *   Ensure `BaseSchema` inheritance.
        *   Correct imports.
    *   [ ] Task 1.5.3: Run linter on `schemas.py` and fix any errors specifically related to `SearchResultRead` and `SearchResultUpdate`.
- **Dependencies:** Relies on `docs/data-models.md` being stable and accurate for these schema definitions.

---

## Change Log

| Change          | Date       | Version | Description             | Author               |
|-----------------|------------|---------|-------------------------|----------------------|
| Initial Draft   | 2025-05-09 | 0.1     | First draft of Epic 1   | Product Manager Agent |
| Added Story 1.5 | 2025-05-10 | 0.2     | Added Story 1.5 for SearchResult Pydantic schema definition and alignment. | Architect Agent      |
| Story 1.2 Update| 2025-05-10 | 0.3     | Expanded Story 1.2 with detailed plan for SearchService refactoring (API alignment, session mgt, new methods) and ReviewService session mgt. | Architect Agent      |
| Story 1.1 Update| 2025-05-10 | 0.4     | Expanded Story 1.1 with detailed plan for `search.py` UI refactoring (service calls, data handling). | Architect Agent      | 