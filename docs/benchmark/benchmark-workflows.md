# Benchmark Sub-System Workflows

This document contains the Mermaid sequence diagrams for the SRA Benchmarking Sub-System.

## Workflow 1: Benchmark Data Seeding

```mermaid
sequenceDiagram
    actor DevUser as "Developer/User"
    participant SeedingScript as "tools/seed_benchmark_data.py"
    participant ProtocolFile as "docs/benchmark/bechmark-protocol.md"
    participant DatasetFile as "docs/benchmark/human-reviewer-results-to-bench-against.csv"
    participant CoreLogic as "Core Logic (Models, Repositories)"
    participant Database

    DevUser->>SeedingScript: Execute script
    activate SeedingScript

    SeedingScript->>ProtocolFile: Read protocol definition
    ProtocolFile-->>SeedingScript: Protocol data
    note right of SeedingScript: User provides PICO XML tagging guidance (conceptually)
    SeedingScript->>SeedingScript: Parse protocol, prepare SystematicReview data

    SeedingScript->>CoreLogic: Create/Get SystematicReview (for benchmark)
    activate CoreLogic
    CoreLogic->>Database: INSERT/SELECT SystematicReview
    Database-->>CoreLogic: "Benchmark Review ID"
    deactivate CoreLogic
    SeedingScript->>SeedingScript: Store "Benchmark Review ID"

    SeedingScript->>DatasetFile: Read human-annotated dataset CSV
    DatasetFile-->>SeedingScript: CSV data rows

    loop For each CSV row
        SeedingScript->>SeedingScript: Map row to SearchResult data (incl. human_decision to source_metadata, infer source_db)
        SeedingScript->>CoreLogic: Create SearchResult linked to "Benchmark Review ID"
        activate CoreLogic
        CoreLogic->>Database: INSERT SearchResult
        Database-->>CoreLogic: Success
        deactivate CoreLogic
    end

    SeedingScript-->>DevUser: "Seeding complete"
    deactivate SeedingScript
```

## Workflow 2: Benchmark Execution & Metrics Calculation

```mermaid
sequenceDiagram
    actor User
    participant BenchmarkUI as "Benchmark UI (human_benchmark_page.py)"
    participant BenchmarkLogic as "Benchmark Logic (orchestration.py)"
    participant LLMAgents as "LLM Agent Layer (screening_agents.py)"
    participant CoreLogicModels as "Core Logic (Models, Repos for BenchmarkRun/Item)"
    participant Database

    User->>BenchmarkUI: Clicks "Run New Benchmark"
    BenchmarkUI->>BenchmarkLogic: initiate_benchmark_run()
    activate BenchmarkLogic

    BenchmarkLogic->>CoreLogicModels: Create BenchmarkRun (initial)
    activate CoreLogicModels
    CoreLogicModels->>Database: INSERT BenchmarkRun
    Database-->>CoreLogicModels: "New BenchmarkRun ID"
    deactivate CoreLogicModels
    BenchmarkLogic->>BenchmarkLogic: Store "BenchmarkRun ID"

    BenchmarkLogic->>CoreLogicModels: Get Benchmark SearchResults
    activate CoreLogicModels
    CoreLogicModels->>Database: SELECT SearchResults for benchmark review
    Database-->>CoreLogicModels: List[SearchResult]
    deactivate CoreLogicModels

    loop For each SearchResult
        BenchmarkLogic->>LLMAgents: screen_abstracts_batch(item, protocol)
        activate LLMAgents
        LLMAgents-->>BenchmarkLogic: "Conservative & Comprehensive Decisions"
        deactivate LLMAgents

        alt Decisions Conflict/Uncertain
            BenchmarkLogic->>LLMAgents: invoke_resolver_chain(item, decisions)
            activate LLMAgents
            LLMAgents-->>BenchmarkLogic: "Resolver Decision"
            deactivate LLMAgents
        end

        BenchmarkLogic->>BenchmarkLogic: Determine "SRA's Final Decision"
        BenchmarkLogic->>BenchmarkLogic: Determine "Classification (TP/FP/TN/FN)"

        BenchmarkLogic->>CoreLogicModels: Create BenchmarkResultItem (with all decisions, classification)
        activate CoreLogicModels
        CoreLogicModels->>Database: INSERT BenchmarkResultItem
        Database-->>CoreLogicModels: Success
        deactivate CoreLogicModels
    end

    BenchmarkLogic->>BenchmarkLogic: Calculate All Summary Metrics (from BenchmarkResultItems)
    BenchmarkLogic->>CoreLogicModels: Update BenchmarkRun with all metrics
    activate CoreLogicModels
    CoreLogicModels->>Database: UPDATE BenchmarkRun
    Database-->>CoreLogicModels: Success
    deactivate CoreLogicModels

    BenchmarkLogic-->>BenchmarkUI: "Benchmark Run Complete (with BenchmarkRun data)"
    deactivate BenchmarkLogic

    BenchmarkUI->>User: Display "Summary Metrics"
    User->>BenchmarkUI: Requests CSV export
    BenchmarkUI->>BenchmarkLogic: get_detailed_csv_data(run_id) or get_summary_csv_data(run_id)
    activate BenchmarkLogic
    BenchmarkLogic-->>BenchmarkUI: CSV Data
    deactivate BenchmarkLogic
    BenchmarkUI-->>User: Provides CSV file download
```
