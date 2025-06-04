# Product Requirements Document (PRD): SRA Benchmarking Module (MVP)

**Version:** 2.3
**Status:** Draft
**Date:** 2025-05-16
**Author:** Product Manager (Gemini)

### 1. Intro

This document outlines the requirements for a new Benchmarking Module within the Systematic Review Assistant (SRA) application. The primary purpose of this module is to evaluate the performance of the SRA's AI-driven Title & Abstract (T&A) screening pipeline (including initial dual review and subsequent conflict resolution) against a predefined, human-annotated dataset ([`docs/benchmark/human-reviewer-results-to-bench-against.csv`](/docs/benchmark/human-reviewer-results-to-bench-against.csv)). This will allow us to measure key accuracy metrics (as defined in [`/docs/sr_metrics.md`](/docs/sr_metrics.md)), identify areas for improvement in the screening agents/resolver, and provide quantitative data for academic reporting. Results of each benchmark run **must** be persisted to the database.

**Note:** This benchmark focuses exclusively on T&A screening. Full-text screening is not yet implemented in the SRA and is out of scope for this module.

### 2. Goals and Context

- **Project Objectives:**
    1. **Data Foundation:** Establish a reliable method to load the benchmark dataset and its corresponding systematic review protocol within the SRA's existing data models (`SystematicReview`, `SearchResult`), accurately reflecting the PICO and exclusion criteria outlined in [`docs/benchmark/bechmark-protocol.md`](/docs/benchmark/bechmark-protocol.md) and refined by user input (with XML tagging for PICO elements in the protocol stored in the database).
    2. **Bug Resolution (Critical Prerequisite):** Ensure the SRA's core abstract screening mechanism can correctly process multiple abstracts in a single batch within the main application (`ScreeningService`).
    3. **Benchmarking Execution & Persistence:** Enable users to trigger the SRA's AI screening pipeline (initial dual review via `screen_abstracts_batch()` and conflict resolution logic via `invoke_resolver_chain()` for relevant cases) on the benchmark dataset. Persist the overall run details (including individual metric columns) and individual item screening outcomes (AI vs. Human) to new database tables.
    4. **Performance Evaluation:** Automatically calculate and display key screening performance metrics by comparing the SRA's final AI decisions against the human-annotated ground truth from the persisted benchmark results.
    5. **Results Export:** Provide a mechanism for users to export both the detailed comparison of AI vs. human decisions and the summary performance metrics from a selected benchmark run.
- **Measurable Outcomes:**
    1. The benchmark protocol (derived from [`docs/benchmark/bechmark-protocol.md`](/docs/benchmark/bechmark-protocol.md) with user-specified PICO/exclusions using XML tags) is successfully converted and stored as a `SystematicReview` record.
    2. The benchmark dataset ([`docs/benchmark/human-reviewer-results-to-bench-against.csv`](/docs/benchmark/human-reviewer-results-to-bench-against.csv)) is successfully converted and stored as `SearchResult` records, linked to the benchmark review, with human ground truth decisions and inferred original `source_db` correctly captured.
    3. The multi-abstract screening bug in the main application is resolved, allowing the benchmark module to process its dataset efficiently.
    4. The benchmark module UI allows users to initiate AI screening of the dataset and view all specified performance metrics.
    5. Users can download detailed comparison data and summary metrics in CSV format.
- **Success Criteria:**
    1. A one-time data seeding process populates the database with the benchmark review protocol (using user-specified PICO/exclusions with XML tags based on details from [`docs/benchmark/bechmark-protocol.md`](/docs/benchmark/bechmark-protocol.md)) and associated search results from [`docs/benchmark/human-reviewer-results-to-bench-against.csv`](/docs/benchmark/human-reviewer-results-to-bench-against.csv) (including human ground truth and inferred original `source_db`).
    2. The multi-result screening bug is fixed in the main application (`screen_abstracts.py` and/or `services.py`), verified by testing.
    3. The benchmark module UI correctly loads and presents the benchmark protocol information.
    4. The AI screening process (using `screen_abstracts_batch` directly, followed by resolver logic) can be initiated from the UI and runs to completion on the full benchmark dataset, with results stored in `BenchmarkRun` (with individual metric columns) and `BenchmarkResultItem` tables.
    5. All defined performance metrics are accurately calculated from stored benchmark results and displayed.
    6. Detailed and summary results are exportable in the specified CSV format.
- **Key Performance Indicators (KPIs):**
    1. 100% of benchmark dataset items processed by the SRA's AI screening pipeline (dual review + resolver logic) and results stored in the database.
    2. Successful calculation and display of all required performance metrics from [`/docs/sr_metrics.md`](/docs/sr_metrics.md).
    3. Ability to export results in the specified formats.

### 3. Scope and Requirements (MVP / Current Version)

#### Functional Requirements (FRs)

- **FR1: Benchmark Data Preparation & Seeding (One-Time)**
    - FR1.1: Develop/update a script (`tools/seed_benchmark_data.py`) to parse the refined PICO elements and exclusion criteria (as provided by user with XML tag instructions, based on [`docs/benchmark/bechmark-protocol.md`](/docs/benchmark/bechmark-protocol.md)) and create a corresponding `SystematicReview` database record. PICO elements **MUST** be stored with XML tags in the `criteria_framework_answers` field of the `SystematicReview` model. `inclusion_criteria` and `exclusion_criteria` strings **MUST** use newline delimiters.
    - FR1.2: The script must parse [`docs/benchmark/human-reviewer-results-to-bench-against.csv`](/docs/benchmark/human-reviewer-results-to-bench-against.csv) and create `SearchResult` database records, linked to the benchmark `SystematicReview`. `source_db` **MUST** be inferred (default `OTHER`). `source_metadata` must store `{"benchmark_human_decision": true/false/null}`.
    - FR1.3: The seeding script must be idempotent.
- **FR2: Fix Multi-Abstract Screening Capability (Core Application - Prerequisite)**
    - FR2.1: The SRA development team must identify and resolve the bug in the existing SRA application that causes errors when screening multiple search results in a batch. This primarily involves ensuring correct transaction management (e.g., single commit per batch) within the `ScreeningService.perform_batch_abstract_screening` method.
- **FR3: Benchmark Database Models & Schemas**
    - FR3.1: Create new SQLModel `BenchmarkRun` (id, timestamp, config_details JSONB). The `BenchmarkRun` table **MUST** include individual typed columns for all key summary metrics defined in [`/docs/sr_metrics.md`](/docs/sr_metrics.md) (e.g., `tp INT`, `fp INT`, `fn INT`, `tn INT`, `sensitivity FLOAT`, `specificity FLOAT`, `accuracy FLOAT`, `f1_score FLOAT`, `mcc_score FLOAT`, `cohen_kappa FLOAT`, `pabak FLOAT`, `lr_plus FLOAT`, `lr_minus FLOAT`).
    - FR3.2: Create new SQLModel `BenchmarkResultItem` (id, benchmark_run_id FK, search_result_id FK, human_decision BOOLEAN, conservative_decision TEXT, comprehensive_decision TEXT, final_decision TEXT, classification TEXT (TP/TN/FP/FN)). Note: `final_decision` in `BenchmarkResultItem` refers to the SRA's overall decision for that item in that run.
    - FR3.3: Define corresponding Pydantic schemas for these models in `schemas.py`.
    - FR3.4: Update [`/docs/data-models.md`](/docs/data-models.md) with these new models and schemas.
- **FR4: Benchmark Module UI & Workflow (`benchmark_tool.py`)**
    - FR4.1: UI button "Load Benchmark Data" (fetches seeded `SystematicReview` and `SearchResult`s).
    - FR4.2: Display benchmark protocol details (stripping XML for UI readability).
    - FR4.3: UI button "Run New Benchmark Screening" that:
        - Creates a new `BenchmarkRun` record (initially without summary metrics).
        - For each benchmark `SearchResult`:
            a.  Invokes `screen_abstracts_batch()` (directly) to get conservative and comprehensive decisions.
            b.  If decisions conflict or involve uncertainty (as per SRA's resolver logic, e.g., one INCLUDE, one EXCLUDE; or one definitive, one UNCERTAIN; or both UNCERTAIN), invoke `invoke_resolver_chain()` (directly) to get a `resolver_decision`.
            c.  Determine the `final_decision_for_benchmark_item` (this is the SRA's overall output: resolver's if available, else agreed upon if no conflict, else derived based on specific handling of UNCERTAIN cases if resolver not triggered).
            d.  Store a new `BenchmarkResultItem` record with all relevant decisions (human_decision from `SearchResult.source_metadata`, conservative_decision, comprehensive_decision, `final_decision_for_benchmark_item` as `final_decision`) and link to the `BenchmarkRun`.
        - Display progress.
- **FR5: Metrics Calculation, Display & Persistence (Benchmark Module)**
    - FR5.1: After all items in a benchmark run are processed and `BenchmarkResultItem` records are created, calculate performance metrics (as per [`/docs/sr_metrics.md`](/docs/sr_metrics.md)) by comparing `BenchmarkResultItem.final_decision` against the human ground truth (`SearchResult.source_metadata.benchmark_human_decision`) for all items associated with that `BenchmarkRun`.
    - FR5.2: Update the corresponding `BenchmarkRun` record by populating its individual metric columns (e.g., `tp`, `sensitivity`, etc.) with the calculated summary metrics.
    - FR5.3: Display these metrics in the UI from the `BenchmarkRun` record.
    - FR5.4: Display statistics for AI confidence scores (conservative, comprehensive, resolver if available) from the `BenchmarkResultItem` records for the selected run.
- **FR6: Results Export (Benchmark Module)**
    - FR6.1: UI to select a `BenchmarkRun` (if multiple runs are supported in future, for MVP can default to latest, or list existing runs by timestamp).
    - FR6.2: Download button for a CSV of detailed `BenchmarkResultItem` data for the selected run: `search_result_source_id`, `title`, `human_decision`, `conservative_decision`, `comprehensive_decision`, `final_decision` (from BenchmarkResultItem, representing SRA's final output for the item), `classification`.
    - FR6.3: Download button for a CSV of summary metrics from the selected `BenchmarkRun` (reading from its individual metric columns).
- **FR7: Adherence to Logging Standards**
    - FR7.1: All new and modified Python code related to the benchmark module and data seeding **MUST** adhere to logging standards outlined in `py/python-logging-rules.mdc` and general [`/docs/coding-standards.md`](/docs/coding-standards.md) (e.g., use `logger.exception()` correctly).

#### Non-Functional Requirements (NFRs)

- **NFR1: Simplicity & Speed of Development:** Implementation must prioritize simplicity and reuse. Direct agent calls for benchmark execution are preferred.
- **NFR2: Correctness of Metrics:** Metric calculations must align with [`/docs/sr_metrics.md`](/docs/sr_metrics.md).
- **NFR3: Data Integrity:** Accurate data seeding of the benchmark protocol (with XML-tagged PICO) and dataset (with correctly inferred original source DB). Metrics stored in `BenchmarkRun` must be typed correctly.
- **NFR4: User Feedback:** Clear UI feedback for all operations.

#### Out of Scope for this MVP

- User uploads for benchmark datasets.
- UI for comparing multiple `BenchmarkRun` records side-by-side.
- User-configurable AI models/prompts for benchmark runs (uses SRA defaults).

### 4. Tech Stack (Leveraging Existing)

- Python, Streamlit, SQLModel, SQLAlchemy, LangChain, Pandas, scikit-learn, Supabase (PostgreSQL).

### 5. Initial Architect & Developer Guidance

- **Data Seeding:** Robust script (`tools/seed_benchmark_data.py`) for protocol (XML PICO, newline criteria) and CSV data (infer `source_db`).
- **Multi-Result Bug:** The core development team needs to address the batch screening bug in `ScreeningService` by ensuring proper session and transaction management (single commit per batch).
- **New DB Models:** `BenchmarkRun` (with individual metric columns), `BenchmarkResultItem` in `models.py`, `schemas.py`, update [`/docs/data-models.md`](/docs/data-models.md). Alembic migration needed.
- **Benchmark Tool (`benchmark_tool.py`):** Orchestrate agent calls (`screen_abstracts_batch`, `invoke_resolver_chain`), persist to new tables, calculate metrics from DB, display, export.
- **Logging:** Strictly follow project logging rules.

### 6. Key Reference Documents

- This Product Requirements Document (PRD) for the SRA Benchmarking Module.
- Systematic Review Metrics: [`/docs/sr_metrics.md`](/docs/sr_metrics.md)
- Benchmark Protocol Definition: [`/docs/benchmark/bechmark-protocol.md`](/docs/benchmark/bechmark-protocol.md)
- Benchmark Ground Truth Dataset: [`/docs/benchmark/human-reviewer-results-to-bench-against.csv`](/docs/benchmark/human-reviewer-results-to-bench-against.csv)
- SRA Architecture Document: [`/docs/architecture.md`](/docs/architecture.md)
- SRA Data Models (to be updated): [`/docs/data-models.md`](/docs/data-models.md)
- SRA Coding Standards: [`/docs/coding-standards.md`](/docs/coding-standards.md)
- SRA Testing Standards: [`docs/testing-strategy.md`](/docs/testing-strategy.md)
- Streamlit UI testing guide: [`docs/streamlit-testing-framework.md`](/docs/streamlit-testing-framework.md)
- Python Logging Rules: [`py/python-logging-rules.mdc`](/.cursor/rules/py/python-logging-rules.mdc) (as referenced in FR7.1)

### 7. Epics

This PRD primarily defines the scope for the following new epic:

- **Epic 4: SRA Benchmarking Module Implementation**
    - **Goal:** To design, develop, and integrate a benchmarking module into the SRA application that allows for the evaluation of the AI-driven T&A screening pipeline against a human-annotated dataset. This includes data preparation and seeding, persistence of benchmark runs and results, calculation and display of performance metrics, and mechanisms for results export, all as detailed in this PRD.
    - **Scope:** The functional and non-functional requirements outlined in sections 3.1 and 3.2 of this PRD constitute the primary scope for this epic. Detailed user stories will be derived from these requirements during backlog grooming and sprint planning.

This PRD aims to guide the development of a functional and valuable benchmarking module within the specified constraints.
