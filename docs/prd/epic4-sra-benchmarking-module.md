# Epic 4: SRA Benchmarking Module Implementation

**Status:** Proposed
**Date Created:** 2025-05-16
**Last Updated:** 2025-05-20
**Related PRD(s):** [`docs/prd-benchmark-may.md`](/docs/prd-benchmark-may.md)
**Architect:** Architect Agent

## 1. Goal

To design, develop, and integrate a benchmarking module into the SRA application that allows for the evaluation of the AI-driven T&A screening pipeline (direct calls to `screen_abstracts_batch` and `invoke_resolver_chain`) against a human-annotated dataset. This includes data preparation and seeding, persistence of benchmark runs and results, calculation and display of performance metrics, and mechanisms for results export, all as detailed in [`docs/prd-benchmark-may.md`](/docs/prd-benchmark-may.md) (Version 2.3).

## 2. Scope

The functional and non-functional requirements outlined in [`docs/prd-benchmark-may.md`](/docs/prd-benchmark-may.md) (sections 3.1 Functional Requirements and 3.2 Non-Functional Requirements) constitute the primary scope for this epic. The module will reside in `src/sr_assistant/benchmark/` and operate independently of the main application's service layer for screening execution, directly calling agent functions.

**Key Architectural Documents:**

- System Architecture: [`docs/architecture.md`](/docs/architecture.md)
- Data Models: [`docs/data-models.md`](/docs/data-models.md)
- Benchmark Workflows: [`docs/benchmark/benchmark-workflows.md`](/docs/benchmark/benchmark-workflows.md)

## 3. Key Features (Summary)

- Benchmark Data Seeding (Protocol from `docs/benchmark/bechmark-protocol.md` & Dataset from `docs/benchmark/human-reviewer-results-to-bench-against.csv`).
- New Database Models (`BenchmarkRun`, `BenchmarkResultItem`) and corresponding Alembic Migration script.
- Benchmark UI (Streamlit page(s) in `src/sr_assistant/benchmark/pages/`):
    - Load and Display Seeded Benchmark Protocol.
    - Trigger New Benchmark Run (utilizing direct agent calls: `screen_abstracts_batch`, `invoke_resolver_chain`).
    - Display Progress of the benchmark run.
- Automated Metrics Calculation, Persistence (in `BenchmarkRun` table), and UI Display.
- Results Export (Detailed `BenchmarkResultItem` data and Summary `BenchmarkRun` metrics as CSV files).
- Adherence to defined Logging Standards (`py/python-logging-rules.mdc` and `docs/coding-standards.md`).

## 4. User Stories (Initial Set)

---
**US4.1: Seed Benchmark Protocol Data**

- **As a:** Developer
- **I want:** A script (`tools/seed_benchmark_data.py`)
- **So that:** I can parse the benchmark protocol definition (from `docs/benchmark/bechmark-protocol.md` with user-guided XML PICO tagging and newline-delimited criteria) and create/update a corresponding `SystematicReview` record in the database.
- **Acceptance Criteria:**
    1. Script correctly parses PICO elements and stores them with XML tags in the `criteria_framework_answers` field of the `SystematicReview` model.
    2. Script correctly parses and stores `inclusion_criteria` and `exclusion_criteria` strings using newline delimiters.
    3. Script is idempotent for protocol seeding (e.g., checks for existing benchmark protocol before creating a new one).
- **Dependencies:** PRD FR1.1
- **Notes:** The script should clearly indicate success or failure and log its actions.

---
**US4.2: Seed Benchmark Dataset**

- **As a:** Developer
- **I want:** The script (`tools/seed_benchmark_data.py`) to parse the benchmark dataset CSV (`docs/benchmark/human-reviewer-results-to-bench-against.csv`)
- **So that:** I can create/update corresponding `SearchResult` records linked to the benchmark `SystematicReview`.
- **Acceptance Criteria:**
    1. Script correctly infers `source_db` for each `SearchResult` (defaulting to `OTHER` if not determinable from CSV).
    2. Script correctly stores the human ground truth decision (e.g., `true`/`false`/`null`) in `SearchResult.source_metadata.benchmark_human_decision`.
    3. Script is idempotent for dataset seeding (e.g., based on a unique identifier from the CSV if available, or clears existing benchmark search results before adding).
- **Dependencies:** PRD FR1.2, US4.1 (benchmark `SystematicReview` must exist)
- **Notes:** Consider how `source_id` for these benchmark `SearchResult` records will be handled if not present in the CSV.

---
**US4.3: Define `BenchmarkRun` Model and Schemas**

- **As a:** Developer
- **I want:** To define the `BenchmarkRun` SQLModel and corresponding Pydantic schemas (`Base`, `Create`, `Update`, `Read`)
- **So that:** Benchmark execution details and all summary result metrics (as per `docs/sr_metrics.md` and PRD FR3.1) can be robustly persisted and accessed.
- **Acceptance Criteria:**
    1. `BenchmarkRun` SQLModel includes `id`, `created_at`, `updated_at`, `benchmark_review_id` (FK), `config_details` (JSONB), `run_notes` (nullable text).
    2. `BenchmarkRun` SQLModel includes individual typed columns for all required metrics (TP, FP, FN, TN, Sensitivity, Specificity, Accuracy, PPV, NPV, F1, MCC, Cohen's Kappa, PABAK, LR+, LR-), all nullable.
    3. Pydantic schemas (`BenchmarkRunCreate`, `BenchmarkRunRead`, `BenchmarkRunUpdate`) are defined, inheriting from a common `BenchmarkRunBase`.
- **Dependencies:** PRD FR3.1, FR3.3; `docs/data-models.md` update.

---
**US4.4: Define `BenchmarkResultItem` Model and Schemas**

- **As a:** Developer
- **I want:** To define the `BenchmarkResultItem` SQLModel and corresponding Pydantic schemas (`Base`, `Create`, `Read`)
- **So that:** Individual item screening outcomes from a benchmark run (human decision vs. all AI decisions, classification) can be persisted and accessed.
- **Acceptance Criteria:**
    1. `BenchmarkResultItem` SQLModel includes `id`, `created_at`, `updated_at`, `benchmark_run_id` (FK), `search_result_id` (FK), `human_decision` (nullable boolean).
    2. Model includes fields for `conservative_decision`, `conservative_confidence`, `conservative_rationale`.
    3. Model includes fields for `comprehensive_decision`, `comprehensive_confidence`, `comprehensive_rationale`.
    4. Model includes fields for `resolver_decision`, `resolver_confidence`, `resolver_reasoning`.
    5. Model includes `final_decision` (representing SRA's output for the run) and `classification` (where classification can be values like "TP", "FP", "TN", "FN").
    6. Pydantic schemas (`BenchmarkResultItemCreate`, `BenchmarkResultItemRead`) are defined, inheriting from `BenchmarkResultItemBase`.
- **Dependencies:** PRD FR3.2, FR3.3; `docs/data-models.md` update.

---
**US4.5: Create Database Migration for Benchmark Tables**

- **As a:** Developer
- **I want:** To create and verify an Alembic migration script
- **So that:** The new `benchmark_runs` and `benchmark_result_items` tables, along with their foreign key constraints, are correctly added to the database schema.
- **Acceptance Criteria:**
    1. Alembic migration script is generated.
    2. `uv run alembic upgrade head` applies the migration successfully.
    3. `uv run alembic downgrade -1` successfully reverts the migration.
    4. Database schema inspection confirms tables and columns as defined in US4.3 & US4.4.
- **Dependencies:** US4.3, US4.4

---
**US4.6: Benchmark UI - Load and Display Benchmark Protocol**

- **As a:** Benchmark User
- **I want:** A UI page (e.g., `src/sr_assistant/benchmark/pages/human_benchmark_page.py`) that loads and displays the details of the seeded benchmark protocol
- **So that:** I can understand the context (research question, criteria) of the benchmark before running it.
- **Acceptance Criteria:**
    1. UI fetches the benchmark `SystematicReview` record from the database.
    2. UI displays research question, background (if any), inclusion criteria, and exclusion criteria.
    3. PICO elements from `criteria_framework_answers` are displayed in a readable format (XML tags stripped).
- **Dependencies:** PRD FR4.1, FR4.2; US4.1 (seeded protocol)

---
**US4.7: Benchmark UI - Trigger Benchmark Run & Process Items**

- **As a:** Benchmark User
- **I want:** A button "Run New Benchmark Screening" on the UI page
- **So that:** I can initiate the AI screening pipeline (direct agent calls) for the entire benchmark dataset.
- **Acceptance Criteria:**
    1. Clicking the button creates a new `BenchmarkRun` record (timestamped, `config_details` populated if applicable).
    2. For each benchmark `SearchResult` linked to the benchmark protocol:
        a.  `screen_abstracts_batch()` is called directly with item data and protocol.
        b.  If decisions conflict/uncertain (per resolver logic), `invoke_resolver_chain()` is called directly.
        c.  The SRA's `final_decision_for_benchmark_item` is determined.
        d.  A new `BenchmarkResultItem` is created and stored with `benchmark_run_id`, `search_result_id`, `human_decision`, all AI decisions (conservative, comprehensive, resolver (if any), `final_decision` representing SRA's output for the run), and `classification` (TP/FP/etc.).
    3. UI displays meaningful progress updates during the run (e.g., "Processing item X of Y", "Calculating metrics...").
- **Dependencies:** PRD FR4.3; US4.1, US4.2, US4.3, US4.4 (seeded data & models)

---
**US4.8: Automated Metrics Calculation and Persistence**

- **As a:** System
- **I want:** After all items in a benchmark run are processed and `BenchmarkResultItem` records are created
- **So that:** Performance metrics (as per `docs/sr_metrics.md`) are automatically calculated by comparing `BenchmarkResultItem.final_decision` (SRA's output for the run) against `BenchmarkResultItem.human_decision` for all items in that run, and the corresponding `BenchmarkRun` record is updated by populating its individual typed metric columns.
- **Acceptance Criteria:**
    1. All metrics (TP, FP, FN, TN, Sensitivity, Specificity, Accuracy, PPV, NPV, F1, MCC, Cohen's Kappa, PABAK, LR+, LR-) are calculated correctly based on the definitions in `docs/sr_metrics.md`.
    2. The `BenchmarkRun` table row for the completed run is updated with all calculated metric values in their respective columns.
- **Dependencies:** PRD FR5.1, FR5.2; US4.7

---
**US4.9: Benchmark UI - Display Summary Performance Metrics**

- **As a:** Benchmark User
- **I want:** The UI page to display all summary performance metrics from a completed `BenchmarkRun` record
- **So that:** I can evaluate the SRA's overall screening performance for that run.
- **Acceptance Criteria:**
    1. UI fetches the `BenchmarkRun` record.
    2. UI clearly displays all populated metric values (TP, FP, ..., LR-) from the record.
    3. If multiple runs exist, UI allows selection of a run to view its metrics (MVP can default to latest).
- **Dependencies:** PRD FR5.3; US4.8

---
**US4.10: Benchmark UI - Display AI Confidence Score Statistics**

- **As a:** Benchmark User
- **I want:** The UI page to display statistics for AI confidence scores
- **So that:** I can understand the distribution and general confidence levels of the AI reviewers for a selected benchmark run.
- **Acceptance Criteria:**
    1. UI fetches relevant `BenchmarkResultItem` data for the selected run.
    2. UI displays summary statistics (e.g., mean, median, min, max, quartiles, or a simple histogram/distribution plot) for:
        a.  `conservative_confidence`
        b.  `comprehensive_confidence`
        c.  `resolver_confidence` (if applicable and available for enough items)
- **Dependencies:** PRD FR5.4; US4.7

---
**US4.11: Benchmark Results Export - Detailed Items**

- **As a:** Benchmark User
- **I want:** A download button on the UI for a selected `BenchmarkRun`
- **So that:** I can download a CSV file of the detailed `BenchmarkResultItem` data for that run.
- **Acceptance Criteria:**
    1. UI provides a button to download "Detailed Results CSV".
    2. The CSV contains columns: `search_result_source_id` (from `SearchResult`), `title` (from `SearchResult`), `human_decision`, `conservative_decision`, `comprehensive_decision`, `resolver_decision` (if applicable), `final_decision` (SRA's output for the run), `classification`.
    3. One row per `BenchmarkResultItem` in the selected run.
- **Dependencies:** PRD FR6.1, FR6.2; US4.7

---
**US4.12: Benchmark Results Export - Summary Metrics**

- **As a:** Benchmark User
- **I want:** A download button on the UI for a selected `BenchmarkRun`
- **So that:** I can download a CSV file of the summary performance metrics for that run.
- **Acceptance Criteria:**
    1. UI provides a button to download "Summary Metrics CSV".
    2. The CSV contains one data row.
    3. Columns correspond to all the metric fields in the `BenchmarkRun` table (e.g., `tp, fp, fn, tn, sensitivity, specificity, ... , lr_minus, run_notes`).
- **Dependencies:** PRD FR6.3; US4.8

---
**US4.13: Adherence to Logging Standards**

- **As a:** Developer
- **I want:** All new and modified Python code related to the benchmark module (including `tools/seed_benchmark_data.py` and code in `src/sr_assistant/benchmark/`)
- **So that:** It adheres to logging standards outlined in `py/python-logging-rules.mdc` and general `docs/coding-standards.md` (e.g., use `logger.exception()` correctly).
- **Acceptance Criteria:**
    1. Code reviews confirm adherence to logging standards.
    2. Appropriate log levels are used.
    3. Exceptions are logged using `logger.exception()` where applicable.
    4. Relevant contextual information is included in log messages.
- **Dependencies:** PRD FR7.1

---
**US4.14: Improve Benchmark Screening Prompts for Conservative Decision-Making**

- **As a:** Developer
- **I want:** To fine-tune the screening agent prompts to be more conservative in exclusion decisions
- **So that:** The AI screening system minimizes false negatives (incorrectly excluding relevant papers) in the benchmark evaluations.
- **Acceptance Criteria:**
    1. System prompts are updated to emphasize that "the absence of evidence is not evidence of absence" (e.g., papers where the country of study is unclear should be included for full-text review rather than excluded).
    2. Conservative and comprehensive agent prompts focus on applying the inclusion/exclusion criteria directly rather than trying to integrate PICO framework concepts.
    3. Prompts explicitly state that the cost of false negatives (excluding relevant papers) is higher than false positives (including irrelevant papers).
    4. Resolver agent prompt is updated to be consistent with the more conservative approach.
    5. Benchmark re-run shows improved sensitivity (reduced false negatives) compared to baseline results.
- **Dependencies:** US4.7 (benchmark execution must be functional), baseline benchmark results for comparison
- **Notes:** This story addresses specific feedback from human reviewers who observed the AI being too aggressive in exclusions during benchmark testing.

## 5. Out of Scope (for this Epic - MVP)

- User uploads for benchmark datasets (uses pre-defined CSV and protocol files).
- UI for comparing multiple `BenchmarkRun` records side-by-side.
- User-configurable AI models/prompts for benchmark runs (uses SRA defaults).
- Fixing the multi-abstract screening bug in the *core application's* `ScreeningService` (FR2) - this Epic assumes the *direct agent calls* `screen_abstracts_batch` and `invoke_resolver_chain` are functional for batch processing as per NFR1 and the PRD's allowance for direct calls for the benchmark module. The fix to `ScreeningService` itself is a prerequisite for other parts of the main app but not a direct blocker for this Epic's benchmark execution path.

## 6. Dependencies

- Availability of functional `screen_abstracts_batch()` and `invoke_resolver_chain()` capable of batch processing (from `src/sr_assistant/app/agents/screening_agents.py`).
- Defined performance metrics in `docs/sr_metrics.md`.
- Access to `docs/benchmark/bechmark-protocol.md` and `docs/benchmark/human-reviewer-results-to-bench-against.csv`.
- Existing core models (`SystematicReview`, `SearchResult`) and their repositories.
- Project coding and logging standards.
- **LangChain Versioning:** The successful execution of benchmark runs via direct calls to `screen_abstracts_batch` depends on the behavior of this function as observed with the specific LangChain versions pinned in `pyproject.toml` (currently `langchain-core==0.3.37` and `langchain==0.3.19`, as listed in `pyproject.toml`). Specifically, the `on_end` callback mechanism is expected to provide data that can be correctly processed into `ScreeningResult` or equivalent detailed structures.

## 7. Risks and Mitigation

- **Risk 1 (Medium):** Direct agent calls (`screen_abstracts_batch`, `invoke_resolver_chain`) might still have underlying issues with batch processing not discovered previously.
    - **Mitigation:** Early, focused testing of these direct agent calls with a small batch from the benchmark dataset.
- **Risk 2 (Low):** Complexity in accurately parsing and mapping the human-annotated CSV to `SearchResult` models, especially inferring `source_db`.
    - **Mitigation:** Clear definition of CSV structure; default `source_db` to `OTHER` if inference is too complex for MVP.
- **Risk 3 (Low):** Metric calculation logic errors.
    - **Mitigation:** Thoroughly test metric calculations against known small datasets with expected outcomes; leverage existing libraries like scikit-learn where appropriate if they simplify implementation without adding undue complexity.
- **Risk 4 (Medium):** Future LangChain upgrades may alter the behavior of the `screen_abstracts_batch` (and its internal callbacks) in a way that breaks the benchmark module's data processing logic if the SRA's core screening chain is not refactored to accommodate such changes.
    - **Mitigation:** Monitor LangChain updates closely. If breaking changes are confirmed and not bugs, the benchmark module may need to adapt its data handling or the SRA's screening chain will require refactoring, which is outside this Epic's direct scope but impacts its long-term viability with newer LangChain versions.

## 8. Release Plan / Milestones (Conceptual)

- **M1:** Data models (`BenchmarkRun`, `BenchmarkResultItem`) defined, migration script created and tested. Data seeding script (`tools/seed_benchmark_data.py`) implemented and functional.
- **M2:** Benchmark UI page created, capable of loading and displaying protocol. "Run New Benchmark" button initiates processing loop (direct agent calls), creates `BenchmarkRun` and `BenchmarkResultItem` records.
- **M3:** Metrics calculation logic implemented and correctly populates `BenchmarkRun` table. Metrics displayed in UI.
- **M4:** CSV export functionality for detailed and summary results implemented. Logging adherence verified.
- **M5:** Final testing and documentation review.

---
