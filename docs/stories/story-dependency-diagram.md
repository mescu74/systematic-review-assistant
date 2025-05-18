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
        E2S2_2 --> E2S2_3
        E2S2_1 --> E2S2_3
        E2S2_5 --> E2S2_4
        E2S2_2 --> E2S2_4
        E2S2_1 --> E2S2_4
    end

    %% Epic 3
    subgraph "Epic 3: Comprehensive Testing & DB Integrity"
        E3S3_1("3.1 Finalize Unit Test Suite (Planned)")
        E3S3_2("3.2 Complete Integration Test Suite (Planned)")
        E3S3_3("3.3 Impl. E2E Test for Core Workflow (Planned)")
        E3S3_4("3.4 Document Safe Integration Test Protocol (Planned)")
        E3S3_5("3.5 Database Integrity & Final Review (Planned)")

        %% Dependencies on Epics for Epic 3 stories
        %% Conceptual links to show Epic 3 stories depend on prior epics being largely done
        E1S1_4 --> E3S3_1
        E2S2_4 --> E3S3_1
        E1S1_4 --> E3S3_2
        E2S2_4 --> E3S3_2
        E3S3_2 --> E3S3_3
        E3S3_4
        E3S3_3 --> E3S3_5
        E3S3_4 --> E3S3_5
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
