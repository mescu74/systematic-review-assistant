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

### Story 1.1: Ensure `SearchResult` Model Consistency in `search.py`

- **User Story / Goal:** As a developer, I want `search.py` to consistently use the `SearchResult` model for displaying and handling PubMed search results, so that data representation is accurate and aligns with the intended generic search architecture.
- **Detailed Requirements:**
  - Review all instances in `src/sr_assistant/app/pages/search.py` where PubMed search results are fetched, processed, or displayed.
  - Ensure that data previously mapped from `PubMedResult` (if applicable) is now correctly handled using `SearchResult` attributes (e.g., `source_id` instead of `pmid`, `source_db` correctly set to 'PubMed').
  - Update UI components in `search.py` (e.g., Streamlit tables, display elements) to reflect `SearchResult` fields.
  - Ensure any actions taken on search results (e.g., selection for screening) correctly reference `SearchResult` instances.
- **Acceptance Criteria (ACs):**
  - AC1: PubMed search results displayed in the `search.py` UI correctly use fields from the `SearchResult` model.
  - AC2: When fetching and processing PubMed results, `search.py` relies on `SearchResult` objects provided by the `SearchService`.
  - AC3: All references to PubMed-specific identifiers (like `pmid`) in `search.py` are replaced with or correctly mapped to `SearchResult.source_id` where `SearchResult.source_db` is 'PubMed'.
- **Tasks (Optional Initial Breakdown):**
  - [ ] Analyze `search.py` for `PubMedResult` usage or direct PubMed API parsing.
  - [ ] Identify UI elements in `search.py` displaying search result data.
  - [ ] Refactor data handling logic to use `SearchResult` objects.
  - [ ] Update Streamlit UI components to bind to `SearchResult` attributes.
- **Dependencies:** Relies on `SearchService` providing `SearchResult` objects.

---

### Story 1.2: Refactor `SearchService` for PubMed Logic & Session Management

- **User Story / Goal:** As a developer, I want the `SearchService` in `services.py` to correctly manage all business logic for PubMed searches and handle database session lifecycle, so that search operations are reliable and adhere to the service layer pattern.
- **Detailed Requirements:**
  - `SearchService` must contain the primary logic for executing PubMed searches (e.g., calling `Entrez` or other PubMed API clients).
  - `SearchService` must correctly manage SQLAlchemy session creation, usage (passing session to repository methods), and closing/committing for PubMed search operations.
  - It should transform raw results from PubMed into `SearchResult` model instances.
  - It should interact with `SearchResultRepository` to store/retrieve `SearchResult` data, passing the managed session.
  - Implement error handling for PubMed API interactions and database operations.
- **Acceptance Criteria (ACs):**
  - AC1: `SearchService` encapsulates PubMed search execution logic.
  - AC2: `SearchService` correctly creates `SearchResult` instances from PubMed API responses, populating `source_db` as 'PubMed' and `source_id` with the PMID.
  - AC3: `SearchService` uses `SearchResultRepository` for all database interactions related to `SearchResult` objects, passing a valid session.
  - AC4: Database sessions are managed appropriately by the `SearchService` (e.g., using `with session_factory()` or `with session_factory.begin()`).
  - AC5: Linter errors in `services.py` related to `SearchService` are resolved.
  - AC6: Unit test coverage for `SearchService` achieves >80%.
  - AC7: All necessary integration tests for `SearchService` functionality are implemented and pass.
- **Tasks (Optional Initial Breakdown):**
  - [ ] Review existing PubMed search logic in `search.py` or `services.py`.
  - [ ] Consolidate PubMed search execution logic into `SearchService`.
  - [ ] Implement robust session management within `SearchService` methods.
  - [ ] Ensure `SearchResult` instances are correctly created and populated.
  - [ ] Integrate calls to `SearchResultRepository` methods.
- **Dependencies:** `SearchResultRepository` (Story 1.3).

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

## Change Log

| Change          | Date       | Version | Description             | Author               |
|-----------------|------------|---------|-------------------------|----------------------|
| Initial Draft   | 2025-05-09 | 0.1     | First draft of Epic 1   | Product Manager Agent | 