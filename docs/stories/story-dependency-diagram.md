# Story Dependency Diagram

This diagram outlines the dependencies between user stories across different epics.

```mermaid
graph TD
    %% Epic 1
    subgraph "Epic 1: Search & Service Layer Stabilization"
        E1S1_5["1.5 Define SearchResult Schemas (Done)"]
        E1S1_3["1.3 Impl. SearchResultRepo Methods (Done)"]
        E1S1_2["1.2 Refactor SearchService (Done)"]
        E1S1_1["1.1 search.py Usage of Service (Done)"]
        E1S1_4["1.4 E2E Search Workflow Stablzn. (Done)"]

        E1S1_5 --> E1S1_2
        E1S1_5 --> E1S1_1
        E1S1_3 --> E1S1_2
        E1S1_2 --> E1S1_1
        E1S1_1 --> E1S1_4
        E1S1_2 --> E1S1_4
        E1S1_3 --> E1S1_4
    end

    %% Epic 2
    subgraph "Epic 2: Resolver Agent & Integration"
        E2S2_1["2.1 Define Core Schemas & Resolver Persistence (Done)"]
        E2S2_2["2.2 Impl. Resolver Agent (Done)"]
        E2S2_5("2.5 Refactor ScreeningService & Resolver Logic (Planned)")
        E2S2_3("2.3 Integrate Resolver (Original Logic) (Planned)")
        E2S2_4("2.4 Integrate Resolver Workflow & UI (Planned)")

        E1S1_5 --> E2S2_1
        E2S2_1 --> E2S2_2
        E2S2_1 --> E2S2_5
        E2S2_2 --> E2S2_5
        E2S2_5 --> E2S2_3
        E2S2_2 --> E2S2_3
        E2S2_1 --> E2S2_3
        E2S2_3 --> E2S2_4
        E2S2_5 --> E2S2_4
    end

    %% Epic 3
    subgraph "Epic 3: Comprehensive Testing & DB Integrity"
        E3S3_1("3.1 Finalize Unit Test Suite (Planned)")
        E3S3_2("3.2 Complete Integration Test Suite (Planned)")
        E3S3_3("3.3 Impl. E2E Test for Core Workflow (Planned)")
        E3S3_4("3.4 Document Safe Integration Test Protocol (Planned)")
        E3S3_5("3.5 Database Integrity & Final Review (Planned)")

        E1S1_4 --> E3S3_1
        E2S2_4 --> E3S3_1
        E1S1_4 --> E3S3_2
        E2S2_4 --> E3S3_2
        E3S3_2 --> E3S3_3
        E3S3_4
        E3S3_3 --> E3S3_5
        E3S3_4 --> E3S3_5
    end

    %% Epic 4
    subgraph "Epic 4: SRA Benchmarking Module Implementation"
        E4S4_1["4.1 Seed Benchmark Protocol Data (Done)"]
        E4S4_2["4.2 Seed Benchmark Dataset (Done)"]
        E4S4_3("4.3 Define BenchmarkRun Model & Schemas (Planned)")
        E4S4_4("4.4 Define BenchmarkResultItem Model & Schemas (Planned)")
        E4S4_5("4.5 Create DB Migration for Benchmark Tables (Planned)")
        E4S4_6("4.6 Benchmark UI - Load & Display Protocol (Planned)")
        E4S4_7("4.7 Benchmark UI - Trigger Run & Process Items (Planned)")
        E4S4_8("4.8 Automated Metrics Calculation & Persistence (Planned)")
        E4S4_9("4.9 Benchmark UI - Display Summary Metrics (Planned)")
        E4S4_10("4.10 Benchmark UI - Display AI Confidence Stats (Planned)")
        E4S4_11("4.11 Benchmark Results Export - Detailed (Planned)")
        E4S4_12("4.12 Benchmark Results Export - Summary (Planned)")
        E4S4_13("4.13 Adherence to Logging Standards (Planned)")

        %% Dependencies for Epic 4
        E2S2_2 --> E4S4_1 
        E2S2_1 --> E4S4_1
        E4S4_1 --> E4S4_2 
        E4S4_1 --> E4S4_6 
        E4S4_2 --> E4S4_7 
        E4S4_6 --> E4S4_7 %% Added: Trigger run UI depends on display protocol UI
        
        E2S2_1 --> E4S4_3 
        E4S4_2 --> E4S4_3 
        E4S4_3 --> E4S4_5 
        
        E2S2_1 --> E4S4_4 
        E4S4_2 --> E4S4_4 
        E4S4_3 --> E4S4_4
        E4S4_4 --> E4S4_5

        E4S4_7 --> E4S4_8 
        E4S4_8 --> E4S4_9 
        E4S4_7 --> E4S4_10 
        E4S4_9 --> E4S4_10 %% Added: Confidence stats display builds upon summary metrics display page
        E4S4_7 --> E4S4_11 
        E4S4_9 --> E4S4_11 %% Added: Detailed export UI is part of the results display page
        E4S4_8 --> E4S4_12 
        E4S4_9 --> E4S4_12 %% Added: Summary export UI is part of the results display page
        
        E1S1_4 --> E4S4_7 %% Benchmark run needs stable search
        E2S2_2 --> E4S4_7 %% Benchmark run needs resolver agent
        E2S2_1 --> E4S4_7 %% Benchmark run needs core schemas

        E2S2_2 --> E4S4_13
        E4S4_12 --> E4S4_13 %% Logging standards apply to all preceding Epic 4 stories implementations
    end

    %% Styling notes from mermaid-synxtax-always.mdc:
    %% - Node labels in quotes: A["Label"]
    %% - GitHub handles styling, so no classDef or style.
    %% Using default shape for "Done" stories: ["ID Text"]
    %% Using rounded rectangle shape for "Planned" stories: ("ID Text")
```

## Diagram Key

- `[Story ID "Description (Done)"]`: Represents a completed story (rectangle shape).
- `(Story ID "Description (Planned)")`: Represents a planned story (rounded rectangle shape).
- Arrows (`-->`) indicate a direct dependency where the preceding story is a prerequisite.
