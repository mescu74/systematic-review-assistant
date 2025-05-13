# API Reference

This document provides a reference for the key APIs within the Systematic Review Assistant application, focusing on the service layer and other critical interfaces.

## Table of Contents

- [Service Layer APIs](#service-layer-apis)
  - [SearchService](#searchservice)
  - [ReviewService](#reviewservice)
  - [ScreeningService](#screeningservice)
- [LLM Agent APIs (Conceptual)](#llm-agent-apis-conceptual)
- [Data Models (Reference to `docs/data-models.md`)](#data-models)

## Service Layer APIs

The service layer acts as an intermediary between the UI/application logic and the data access layer (repositories). Services are responsible for business logic, orchestrating calls to repositories, and managing database sessions. All service methods are synchronous.

Refer to `src/sr_assistant/app/services.py` for the source code.
Refer to `docs/data-models.md` for details on data structures used in requests and responses.

### SearchService

Manages search operations against external academic databases (e.g., PubMed, Scopus). It encapsulates the logic for performing searches, fetching raw data, mapping it to the internal `SearchResult` domain model, and storing these results using `SearchResultRepository`.

#### Key Principles for this Service API
- The service abstracts away the specifics of interacting with different external search APIs.
- Callers provide search parameters and receive standardized `SearchResult` objects or relevant status information.
- Raw API records from external sources are *not* exposed by the service's public interface.
- **Session Management:** Each public service method call is responsible for managing its own database session and transaction. Session objects are not exposed to or required from the calling layer for standard operations.

#### Dependencies
- `repositories.SearchResultRepository`
- (Internal) HTTP client or specific API client libraries (e.g., BioPython for Entrez)

#### Revised Key Methods

- **`__init__(self, factory: sessionmaker[Session] = session_factory, search_repo: SearchResultRepository | None = None, http_client: httpx.AsyncClient | None = None)`**
  - Initializes the service. The `factory` parameter allows injection of a SQLAlchemy `sessionmaker` for database sessions, typically defaulting to the application's main `session_factory` from `sr_assistant.app.database`.

- **`search_pubmed_and_store_results(self, review_id: uuid.UUID, query: str, max_results: int = 100) -> Sequence[schemas.SearchResultRead]`**
  - **Purpose:** Performs a PubMed search, maps results to `schemas.SearchResultRead`, and stores them (as `models.SearchResult`).
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** Sequence of `schemas.SearchResultRead` objects.
  - **Session Handling:** Manages its own session and transaction internally.

- **`search_scopus_and_store_results(self, review_id: uuid.UUID, query: str, max_results: int = 100) -> Sequence[schemas.SearchResultRead]`** (Conceptual)
  - **Purpose:** Performs a Scopus search, maps, and stores results, returning `schemas.SearchResultRead`.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** Sequence of `schemas.SearchResultRead` objects from Scopus.
  - **Session Handling:** Manages its own session and transaction internally.

- **`_map_pubmed_to_search_result(self, review_id: uuid.UUID, api_record: dict[str, Any]) -> models.SearchResult | None`** (Internal)
  - *Internal method.* Maps a raw PubMed API record to a `SearchResult` model.
  - **Note:** `api_record` type hint should be `dict[str, Any]`.
  - **Note:** Reconcile `SearchResult.year` type (see previous notes: recommend `int | None`).

- **`_map_scopus_to_search_result(self, review_id: uuid.UUID, api_record: dict[str, Any]) -> models.SearchResult | None`** (Internal)
  - *Internal method.* Maps a raw Scopus API record to a `SearchResult` model.
  - **Note:** `api_record` type hint should be `dict[str, Any]`.

- **`get_search_results_by_review_id(self, review_id: uuid.UUID) -> Sequence[schemas.SearchResultRead]`**
  - **Purpose:** Retrieves all search results for a review as `schemas.SearchResultRead` objects.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** Sequence of `schemas.SearchResultRead` objects.
  - **Session Handling:** Manages its own session internally.

- **`get_search_result_by_source_details(self, review_id: uuid.UUID, source_db: SearchDatabaseSource, source_id: str) -> schemas.SearchResultRead | None`**
  - **Purpose:** Retrieves a specific search result as a `schemas.SearchResultRead` object.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** `schemas.SearchResultRead` object or `None`.
  - **Session Handling:** Manages its own session internally.

- **`update_search_result(self, result_id: uuid.UUID, update_data: schemas.SearchResultUpdate) -> schemas.SearchResultRead`**
  - **Purpose:** Updates an existing search result and returns the updated data as `schemas.SearchResultRead`.
  - **Parameters:** (as previously defined, `session` parameter removed)
    - `update_data`: A `schemas.SearchResultUpdate` Pydantic object.
  - **Returns:** The updated `schemas.SearchResultRead` object.
  - **Session Handling:** Manages its own session and transaction internally.

- **`delete_search_result(self, result_id: uuid.UUID) -> None`**
  - **Purpose:** Deletes a search result.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Session Handling:** Manages its own session and transaction internally.

#### Key Removals / Changes & Refinements
- **REMOVED:** `add_api_search_results`. This method is fundamentally flawed as it breaks service encapsulation. The new `search_pubmed_and_store_results` (and similar for other sources) replaces its intended functionality correctly.
- `prd-recovery.md` (FR1) emphasizes stabilizing PubMed first. The `search_pubmed_and_store_results` method directly addresses this. Ensure all PubMed-related paths are robust internally.
- A new Pydantic schema `schemas.SearchResultUpdate` should be defined in `src/sr_assistant/core/schemas.py` for the `update_search_result` method.

### ReviewService

Manages systematic review data, interacting with `SystematicReviewRepository`.

#### Key Principles for this Service API
- The service provides a clear interface for CRUD operations on systematic reviews.
- Input for creation and updates uses specific Pydantic schemas for validation and clarity.
- **Session Management:** Each public service method call is responsible for managing its own database session and transaction. Session objects are not exposed to or required from the calling layer.

#### Dependencies
- `repositories.SystematicReviewRepository`

#### Revised Key Methods

- **`__init__(self, factory: sessionmaker[Session] = session_factory, review_repo: SystematicReviewRepository | None = None)`**
  - Initializes the service. The `factory` parameter allows injection of a SQLAlchemy `sessionmaker`, typically defaulting to the application's main `session_factory` from `sr_assistant.app.database`.

- **`create_review(self, review_data: schemas.SystematicReviewCreate) -> schemas.SystematicReviewRead`**
  - **Purpose:** Creates a new systematic review.
  - **Parameters:**
    - `review_data`: Pydantic schema `schemas.SystematicReviewCreate` with review details.
  - **Returns:** The created `schemas.SystematicReviewRead` object.
  - **Session Handling:** Manages its own session and transaction internally.

- **`get_review(self, review_id: uuid.UUID) -> schemas.SystematicReviewRead | None`**
  - **Purpose:** Retrieves a review by its ID.
  - **Parameters:**
    - `review_id`: ID of the review.
  - **Returns:** The `schemas.SystematicReviewRead` object or `None`.
  - **Session Handling:** Manages its own session internally.

- **`get_all_reviews(self) -> Sequence[schemas.SystematicReviewRead]`**
  - **Purpose:** Retrieves all reviews.
  - **Returns:** A sequence of `schemas.SystematicReviewRead` objects.
  - **Session Handling:** Manages its own session internally.

- **`update_review(self, review_id: uuid.UUID, review_update_data: schemas.SystematicReviewUpdate) -> schemas.SystematicReviewRead`**
  - **Purpose:** Updates an existing review.
  - **Parameters:**
    - `review_id`: ID of the review to update.
    - `review_update_data`: Pydantic schema `schemas.SystematicReviewUpdate` with fields to update.
  - **Returns:** The updated `schemas.SystematicReviewRead` object.
  - **Session Handling:** Manages its own session and transaction internally.

- **`delete_review(self, review_id: uuid.UUID) -> None`**
  - **Purpose:** Deletes a review by its ID.
  - **Parameters:**
    - `review_id`: ID of the review to delete.
  - **Session Handling:** Manages its own session and transaction internally.

#### Refinements (based on PRDs)
- This service appears largely unaffected by the immediate recovery tasks in `prd-recovery.md` and `prd-resolver.md`, but its API is now aligned with the stricter encapsulation principles.

### ScreeningService

Manages screening decisions (from conservative and comprehensive LLM reviewers) and the conflict resolution process. Interacts with `ScreenAbstractResultRepository`, `ScreeningResolutionRepository`, and `SearchResultRepository`.

#### Key Principles for this Service API
- Methods for adding or updating screening information will expect structured Pydantic schemas (e.g., `schemas.ScreeningResultCreate` or `schemas.ScreeningResult`) as input, not raw LLM response objects like `schemas.ScreeningResponse`.
- The service is responsible for validating these inputs and mapping them correctly to the `models.ScreenAbstractResult` SQLModel for database interaction.
- **Session Management:** Each public service method call is responsible for managing its own database session and transaction. Session objects are not exposed to or required from the calling layer for standard operations.

#### Dependencies
- `repositories.ScreenAbstractResultRepository`
- `repositories.ScreeningResolutionRepository`
- `repositories.SearchResultRepository`

#### Revised Key Methods

- **`__init__(self, factory: sessionmaker[Session] = session_factory, screen_repo: ScreenAbstractResultRepository | None = None, resolution_repo: ScreeningResolutionRepository | None = None, search_repo: SearchResultRepository | None = None)`**
  - Initializes the service with necessary repositories and a session factory. The `factory` parameter allows injection of a SQLAlchemy `sessionmaker`, typically defaulting to the application's main `session_factory` from `sr_assistant.app.database`.

- **`add_screening_decision(self, search_result_id: uuid.UUID, screening_strategy: ScreeningStrategyType, screening_data: schemas.ScreeningResultCreate) -> schemas.ScreeningResultRead`**
  - **Purpose:** Adds a new screening decision to the system. It creates a `ScreenAbstractResult` record and then updates the corresponding `SearchResult` (identified by `search_result_id`) to link to this new record via its `conservative_result_id` or `comprehensive_result_id` field, based on the `screening_strategy`.
  - **Parameters:**
    - `search_result_id`: The ID of the `SearchResult` to which this screening decision pertains.
    - `screening_strategy`: The strategy (CONSERVATIVE or COMPREHENSIVE) used for this decision. This determines which field (`conservative_result_id` or `comprehensive_result_id`) on the `SearchResult` is updated.
    - `screening_data`: A `schemas.ScreeningResultCreate` Pydantic object containing the data for the new screening result (e.g., decision, rationale, LLM run ID used as `ScreenAbstractResult.id`).
  - **Returns:** The created `schemas.ScreeningResultRead` object (representing the `ScreenAbstractResult`).
  - **Session Handling:** Manages its own session and transaction internally.
  - **Logic:**
    1. Validates `screening_data` (which includes `review_id` and the LLM `run_id` to be used as `ScreenAbstractResult.id`).
    2. Creates a `models.ScreenAbstractResult` instance using data from `screening_data`.
    3. Uses `ScreenAbstractResultRepository.add` to persist the new screening record.
    4. Fetches the `models.SearchResult` using `search_result_id` via `SearchResultRepository.get_by_id`.
    5. Based on `screening_strategy`, sets either `SearchResult.conservative_result_id` or `SearchResult.comprehensive_result_id` to the `id` of the newly created `ScreenAbstractResult`.
    6. Uses `SearchResultRepository.update` (or adds to session) to save changes to the `SearchResult`.
    7. Commits the session and refreshes the created `ScreenAbstractResult` instance.

- **`update_screening_decision(self, screening_result_id: uuid.UUID, screening_update_data: schemas.ScreeningResultUpdate) -> schemas.ScreeningResultRead`**
  - **Purpose:** Updates an existing screening decision.
  - **Parameters:**
    - `screening_result_id`: The ID of the `ScreenAbstractResult` to update.
    - `screening_update_data`: A `schemas.ScreeningResultUpdate` Pydantic object.
  - **Returns:** The updated `schemas.ScreeningResultRead` object.
  - **Session Handling:** Manages its own session and transaction internally.

- **`get_screening_decision_by_id(self, screening_result_id: uuid.UUID) -> schemas.ScreeningResultRead | None`**
  - **Purpose:** Retrieves a specific screening decision by its ID.
  - **Parameters:**
    - `screening_result_id`: The ID of the screening result.
  - **Returns:** The `schemas.ScreeningResultRead` object or `None`.
  - **Session Handling:** Manages its own session internally.

- **`get_screening_decisions_for_search_result(self, search_result_id: uuid.UUID) -> Sequence[schemas.ScreeningResultRead]`**
  - **Purpose:** Retrieves all screening decisions associated with a specific `SearchResult`.
  - **Parameters:**
    - `search_result_id`: The ID of the `SearchResult`.
  - **Returns:** A sequence of `schemas.ScreeningResultRead` objects.
  - **Session Handling:** Manages its own session internally.

- **`get_screening_decisions_for_review(self, review_id: uuid.UUID) -> Sequence[schemas.ScreeningResultRead]`**
  - **Purpose:** Retrieves all screening results for a specific review.
  - **Parameters:**
    - `review_id`: ID of the review.
  - **Returns:** Sequence of `schemas.ScreeningResultRead` objects.
  - **Session Handling:** Manages its own session internally.

#### Needed Methods for Conflict Resolution (based on `prd-resolver.md` FR1-FR8 & TD4)

- **`identify_disagreements(self, review_id: uuid.UUID, search_result_ids: list[uuid.UUID]) -> list[schemas.SearchResultRead]`** (New)
  - **Purpose:** Identifies search results that require conflict resolution.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** List of `schemas.SearchResultRead` objects with disagreements.
  - **Session Handling:** Manages its own session internally.

- **`prepare_resolver_inputs(self, review: models.SystematicReview, search_results_with_disagreements: list[models.SearchResult]) -> list[dict[str, Any]]`** (New)
  - **Purpose:** Prepares the input data (dictionary of prompt variables) for the resolver LLM chain for a list of search results that have disagreements (`prd-resolver.md` FR3).
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Returns:** A list of dictionaries, where each dictionary contains the prompt variables required by the `resolver_prompt` (e.g., `{'search_result_title': '...', 'search_result_abstract': '...', 'conservative_decision': '...', ...}`). The exact structure is determined by `resolver_chain.get_input_schema()`.
  - **Logic:**
    1. For each `SearchResult`:
       a. Gather its full data.
       b. Get the associated CONSERVATIVE and COMPREHENSIVE `ScreenAbstractResult` objects.
       c. Get relevant `SystematicReview` protocol data.
       d. Construct the dictionary of prompt variables.

- **`invoke_resolver_agent_batch(self, resolver_prompt_variable_inputs: list[dict[str, Any]]) -> list[schemas.ScreeningResolutionSchema]`** (New)
  - **Purpose:** Invokes the resolver LLM agent/chain in batches with the prepared prompt variable inputs (`prd-resolver.md` FR4).
  - **Parameters:**
    - `resolver_prompt_variable_inputs`: A list of dictionaries, where each dictionary represents the input variables for a single invocation of the resolver chain.
  - **Returns:** A list of `schemas.ScreeningResolutionSchema` Pydantic objects (as defined in `src/sr_assistant/core/schemas.py` and `prd-resolver.md` FR5).
  - **Logic:**
    1. Get the resolver LLM chain (e.g., from `screening_agents.py`).
    2. Invoke the chain with the batch of inputs (`resolver_chain.batch(resolver_prompt_variable_inputs, config=...)`).
    3. Handle LLM call errors (retries, logging).

- **`store_resolution_results(self, review_id: uuid.UUID, search_result_id_to_resolution_data: dict[uuid.UUID, schemas.ScreeningResolutionSchema]) -> list[schemas.ScreeningResolutionRead]`** (New)
  - **Purpose:** Stores outcomes from the resolver agent.
  - **Parameters:**
    - `review_id`: ID of the current systematic review.
    - `search_result_id_to_resolution_data`: A dictionary mapping `SearchResult.id` to its corresponding `schemas.ScreeningResolutionSchema` data from the resolver agent.
  - **Returns:** List of created `schemas.ScreeningResolutionRead` objects.
  - **Session Handling:** Manages its own session and transaction internally.
  - **Logic:**
    1. Iterate through `search_result_id_to_resolution_data`.
    2. For each item, create `models.ScreeningResolution` using data from `schemas.ScreeningResolutionSchema`. Ensure all necessary IDs (e.g., `review_id`, `search_result_id`, `conservative_result_id`, `comprehensive_result_id`) are correctly populated. The `conservative_result_id` and `comprehensive_result_id` might need to be fetched or passed through if not part of `ScreeningResolutionSchema` directly (though the current schema has them as optional and populated by caller).
    3. Add to session via `ScreeningResolutionRepository.add`.
    4. Update the corresponding `models.SearchResult.final_decision` and `resolution_id`.
    5. Commit session.

- **`resolve_screening_conflicts_for_batch(self, review: models.SystematicReview, search_result_ids_in_batch: list[uuid.UUID]) -> None`** (New - Orchestration Method)
  - **Purpose:** Orchestrates the entire conflict resolution process for a batch of search results.
  - **Parameters:** (as previously defined, `session` parameter removed)
  - **Session Handling:** Manages its own session(s) and transaction(s) internally for the sequence of operations.
  - **Logic:** (as previously defined, using the revised methods above)

#### General Notes for `ScreeningService` Implementation
- Ensure Pydantic schemas for service inputs (e.g., `ScreeningResultCreate`, `ScreeningResultUpdate`) are clearly defined in `src/sr_assistant/core/schemas.py`.
- The output of the resolver chain is `schemas.ScreeningResolutionSchema`.
- The input to the resolver chain (`invoke_resolver_agent_batch`) is a list of dictionaries, each matching the input variables of the `resolver_prompt`.
- **Pydantic Model Instantiation:** When creating SQLModel instances (e.g., `models.ScreenAbstractResult`) from Pydantic schema data (e.g., `schemas.ScreeningResultCreate`), use the `ModelName.model_validate(data_dict)` method (e.g., `models.ScreenAbstractResult.model_validate(screening_create_data.model_dump())`) as per `docs/coding-standards.md`. Avoid direct instantiation with `**data_dict` if it leads to type errors or bypasses Pydantic v2 validation strictness.
- Address the unused `e` variable linter error in the current `services.py` (`get_screening_result_by_strategy`).
- Type hints for generic `dict` should be `dict[str, Any]` or a more specific `TypedDict` if applicable.

#### Further Considerations for `ScreeningService`
- The definition of Pydantic schemas used for service method inputs (e.g., `schemas.ScreeningResultCreate`, `schemas.ScreeningResultUpdate`) needs to be solidified and placed in `src/sr_assistant/core/schemas.py`. The `schemas.ScreeningResolutionSchema` is the output from the resolver LLM chain.
- Error handling within these new methods needs to be robust.

## LLM Agent APIs (Conceptual)

This section outlines the conceptual interfaces for the key LLM agent chains and their primary invocation functions used within the application. These describe the data contracts for invoking these functionalities.

### Abstract Screening (`screen_abstracts_batch` function in `sr_assistant.app.agents.screening_agents`)

This function orchestrates the screening of a batch of search results using a parallel LLM chain (`screen_abstracts_chain`) that internally runs two review strategies: "conservative" and "comprehensive".

- **Primary Invocation:** `screen_abstracts_batch(batch: list[models.SearchResult], batch_idx: int, review: models.SystematicReview) -> ScreenAbstractsBatchOutput | None`
- **Core LLM Chain Input (for each item in batch, processed by `make_screen_abstracts_chain_input`):**
  - `background`: str | None (Background of the systematic review)
  - `research_question`: str (Research question)
  - `inclusion_criteria`: str (Inclusion criteria)
    - **IMPORTANT NOTE:** The `screen_abstracts_chain` and its underlying prompts currently expect the inclusion criteria as a single string via this `inclusion_criteria` field. It does *not* directly use the `SystematicReview.criteria_framework_answers` dictionary (which stores structured PICO etc. data). The `src/sr_assistant/app/pages/protocol.py` page currently populates `SystematicReview.inclusion_criteria` by joining the PICO fields into a string, and this is the string used by the screening chain. This is a known area for future refactoring to make criteria handling more robust and directly use structured framework answers.
  - `exclusion_criteria`: str (Exclusion criteria)
  - `title`: str (Title of the search result)
  - `journal`: str (Journal of the search result)
  - `year`: str (Year of the search result)
  - `abstract`: str (Abstract of the search result)
  (These correspond to `ScreenAbstractsChainInputDict`)

- **Core LLM Chain Output (from `screen_abstracts_chain`):**
  - The `screen_abstracts_chain` (a `RunnableParallel` instance) returns a Python `dict`. 
  - Due to an interaction where `.with_types()` can override `.with_listeners()`, the chain does not formally return an instance of the `ScreenAbstractsChainOutputDict` TypedDict. However, the output `dict` structurally conforms to `ScreenAbstractsChainOutputDict` (defined in `sr_assistant.app.agents.screening_agents`).
  - The keys of this output dictionary are `"conservative"` and `"comprehensive"`.
  - The value for each key is determined by the success of the corresponding sub-chain and the `screen_abstracts_chain_on_end_cb` listener:
    - Ideally, `schemas.ScreeningResult`: If the respective sub-chain (e.g., conservative reviewer) successfully produced a `schemas.ScreeningResponse` and the `screen_abstracts_chain_on_end_cb` listener successfully processed this into a `schemas.ScreeningResult`.
    - Potentially, `schemas.ScreeningResponse`: If the sub-chain produced a `schemas.ScreeningResponse`, but the `screen_abstracts_chain_on_end_cb` listener encountered an error during its processing. In this case, the original `schemas.ScreeningResponse` from the LLM might be present.
    - Potentially, `t.Any`: If the sub-chain itself failed to produce a parsable `schemas.ScreeningResponse` (e.g., due to an LLM or tool error within the sub-chain before the listener could act).
  - *The `screen_abstracts_batch` function is responsible for handling these varied outcomes and packaging them into `ScreenAbstractResultTuple`s, using `ScreeningError` where appropriate.*

- **`screen_abstracts_batch` Function Output:** `schemas.ScreenAbstractsBatchOutput` (a NamedTuple containing `results: list[ScreenAbstractResultTuple]` and `cb: OpenAICallbackHandler`).
  - Each `ScreenAbstractResultTuple` contains:
    - `search_result`: `models.SearchResult`
    - `conservative_result`: `schemas.ScreeningResult` or `schemas.ScreeningError`.
    - `comprehensive_result`: `schemas.ScreeningResult` or `schemas.ScreeningError`.

### Resolver Chain (`resolver_chain` - conceptual, to be built based on `prd-resolver.md`)

This chain is responsible for resolving disagreements between the conservative and comprehensive screening decisions.

- **Input (for each disagreement, prepared by `ScreeningService.prepare_resolver_inputs`):**
  - A dictionary of prompt variables, conceptually including:
    - Full `models.SearchResult` data (PMID, title, abstract, etc.)
    - Relevant `models.SystematicReview` protocol data (research question, PICO, exclusion criteria)
    - The `schemas.ScreeningResult` from the "conservative" reviewer.
    - The `schemas.ScreeningResult` from the "comprehensive" reviewer.
  - The exact input variable names and structure will be determined by the `resolver_prompt` and can be inspected via `resolver_chain.get_input_schema()` once the chain is defined.
- **Output (from the chain):** `schemas.ScreeningResolutionSchema` (Pydantic model containing `resolver_decision`, `resolver_reasoning`, etc. as defined in `src/sr_assistant/core/schemas.py`)

## Data Models

Refer to [Data Models Document](data-models.md) for detailed definitions of SQLModel database models (e.g., `SearchResult`, `SystematicReview`, `ScreenAbstractResult`, `ScreeningResolution`) and core Pydantic schemas (e.g., `ScreeningResultCreate`, `ScreeningResultUpdate`, `SearchResultUpdate`, `ScreeningResponse`, `ScreeningResolutionSchema`) used throughout the application.

## Change Log

| Date       | Version | Description                                                                                                                               | Author          |
|------------|---------|-------------------------------------------------------------------------------------------------------------------------------------------|-----------------|
| 2025-05-10 | 0.1     | Initial draft created (Placeholder or basic outline from project recovery commit).                                                           | System/User     |
| 2025-05-13 | 0.2     | Populated initial service API outlines (Search, Review, Screening) based on `services.py`.                                                | Architect Agent |
| 2025-05-13 | 0.3     | Major Refactor: Corrected SearchService (encapsulation), ScreeningService (inputs, resolver methods). Defined all key Pydantic Read/Update schemas. | Architect Agent |
| 2025-05-13 | 0.4     | Session Mgt Update: Removed session parameters from all public service methods; enforced internal session management & factory injection.     | Architect Agent |
| 2025-05-13 | 0.5     | LLM I/O & Formatting: Clarified LLM chain I/O (`ScreenAbstractsChainOutputDict`, `ScreeningResolutionSchema`). Corrected Markdown H4/bolding. Restored sections. Added this changelog. | Architect Agent |
| 2025-05-13 | 0.6     | Service Return Types: Updated all service methods to return Pydantic `Read` schemas instead of SQLModels. Added `model_validate` note.          | Architect Agent |