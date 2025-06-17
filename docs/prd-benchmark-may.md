# SRA Benchmarking Module PRD

## Goals and Background Context

### Goals

- Establish a reliable method to load a benchmark dataset and its corresponding systematic review protocol.
- Enable users to trigger the SRA's AI screening pipeline on the benchmark dataset.
- Persist benchmark run details and individual screening outcomes to the database.
- Automatically calculate and display key screening performance metrics against a human-annotated ground truth.
- Provide a mechanism to export detailed and summary benchmark results.

### Background Context

To quantitatively measure the performance of the SRA's AI-driven Title & Abstract screening pipeline, a new Benchmarking Module is required. This module will evaluate the entire AI workflow (dual review + conflict resolution) against a predefined, human-annotated dataset. This will allow for the calculation of key accuracy metrics, help identify areas for improvement in the AI agents, and provide verifiable data for reporting.

### Change Log

| Date       | Version | Description                   | Author |
| :--------- | :------ | :---------------------------- | :----- |
| 2024-07-30 | 3.0     | Migrated to v4 PRD Template   | AIDE   |
| 2025-05-16 | 2.3     | Initial Draft                 | PM Agent |

## Requirements

### Functional

- **FR1:** A data seeding script (`tools/seed_benchmark_data.py`) must be created to parse a benchmark protocol and dataset CSV, creating `SystematicReview` and `SearchResult` records.
- **FR2:** The `SystematicReview` record must store PICO elements with XML tags and criteria with newline delimiters.
- **FR3:** The `SearchResult` records must store the human ground truth decision in their `source_metadata` field.
- **FR4:** The seeding script must be idempotent.
- **FR5:** A critical bug causing errors when screening multiple abstracts in a batch must be fixed in the core application's `ScreeningService`.
- **FR6:** A new `BenchmarkRun` database model must be created to store run details and summary metrics in dedicated, typed columns.
- **FR7:** A new `BenchmarkResultItem` database model must be created to store the outcome of each individual screening comparison (AI vs. Human).
- **FR8:** A new UI page for the benchmark module must be created.
- **FR9:** The UI must allow a user to load the seeded benchmark data and initiate a new benchmark run.
- **FR10:** Initiating a run will create a `BenchmarkRun` record and then process each item by calling the SRA's `screen_abstracts_batch()` and `invoke_resolver_chain()` functions directly.
- **FR11:** A `BenchmarkResultItem` record will be stored for each processed item, containing all AI decisions and the human ground truth.
- **FR12:** After a run completes, all summary performance metrics must be calculated and persisted into the appropriate columns of the `BenchmarkRun` record.
- **FR13:** The UI must display all summary metrics from a completed `BenchmarkRun`.
- **FR14:** The UI must provide buttons to download the detailed `BenchmarkResultItem` data and the summary `BenchmarkRun` metrics as separate CSV files.

### Non Functional

- **NFR1:** The implementation should prioritize simplicity and reuse of existing SRA components (e.g., direct agent calls).
- **NFR2:** Metric calculations must be correct and align with the definitions in `docs/sr_metrics.md`.
- **NFR3:** Data seeding and storage must ensure high data integrity.
- **NFR4:** All new and modified code must adhere to project logging standards.

## User Interface Design Goals

A simple, functional UI is required for the MVP.
- A single page for the benchmark module.
- Buttons to "Load Benchmark Data" and "Run New Benchmark Screening".
- A clear display area for the benchmark protocol details.
- A table or list to display summary metrics from a completed run.
- Download buttons for exporting results.

## Technical Assumptions

### Repository Structure

Monorepo (existing)

### Service Architecture

The benchmark module will be a new sub-system within the existing service-oriented monolith, with its own UI pages and logic that directly call the core SRA agent functions.

### Testing requirements

- Unit tests for the data seeding script.
- Unit tests for the metric calculation logic.
- Integration tests for the benchmark workflow, mocking the LLM calls.

## Epics

### Epic 1: Data Foundation and Core Logic

**Goal:** Establish the database models, schemas, and data seeding script required for the benchmarking module.

#### Story 1.1: Create Benchmark DB Models

As a developer, I want to create the `BenchmarkRun` and `BenchmarkResultItem` SQLModels, so that I can persist all data from a benchmark run.

##### Acceptance Criteria

- 1: `BenchmarkRun` model is created with typed columns for all summary metrics.
- 2: `BenchmarkResultItem` model is created to store individual comparisons.
- 3: Pydantic schemas for the new models are defined.
- 4: An Alembic migration is created and applied successfully.

#### Story 1.2: Implement Data Seeding Script

As a developer, I want a script to populate the database with the benchmark protocol and dataset, so that the module has data to operate on.

##### Acceptance Criteria

- 1: `tools/seed_benchmark_data.py` is created.
- 2: The script correctly parses the protocol markdown and dataset CSV.
- 3: The script creates `SystematicReview` and `SearchResult` records with the correct data and metadata.
- 4: The script is idempotent and can be run multiple times without creating duplicate data.

### Epic 2: Benchmark Execution and UI

**Goal:** Build the user interface and orchestration logic to execute a benchmark run, calculate metrics, and display/export the results.

#### Story 2.1: Fix Core Screening Bug

As a developer, I must fix the multi-abstract screening bug in the main application, so that the benchmark module can process its dataset reliably.

##### Acceptance Criteria

- 1: The root cause of the batch screening error in `ScreeningService` is identified.
- 2: The bug is fixed, ensuring correct transaction management.
- 3: A test is created to verify that batch screening now works as expected.

#### Story 2.2: Build Benchmark UI

As a user, I want a simple interface to load data, run a benchmark, and see the results, so that I can evaluate the SRA's performance.

##### Acceptance Criteria

- 1: A new Streamlit page for the benchmark module is created.
- 2: UI buttons exist to load data and start a new run.
- 3: The benchmark protocol is displayed correctly on the page.
- 4: A status indicator provides feedback during a run.

#### Story 2.3: Implement Benchmark Orchestration

As a developer, I want to write the logic that orchestrates the benchmark run, so that it correctly processes the data and calculates metrics.

##### Acceptance Criteria

- 1: The "Run New Benchmark" button triggers the orchestration logic.
- 2: The logic correctly iterates through the dataset, calls the SRA screening agents, and stores `BenchmarkResultItem` records.
- 3: After the run, summary metrics are calculated correctly.
- 4: The summary metrics are persisted to the `BenchmarkRun` record in the database.

#### Story 2.4: Display and Export Results

As a user, I want to view the summary metrics and download the full results, so that I can analyze the SRA's performance in detail.

##### Acceptance Criteria

- 1: The summary metrics from the latest `BenchmarkRun` are displayed in the UI.
- 2: A download button provides a CSV of the detailed `BenchmarkResultItem` data.
- 3: A second download button provides a CSV of the summary metrics.

## Checklist Results Report

[[LLM: This section will be populated after running the `pm-checklist` task.]]

## Next Steps

[[LLM: This section will contain prompts for the Design Architect and Architect.]]
