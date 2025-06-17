# Documentation Index
This document serves as the central hub for all project documentation. It is organized to help you quickly find information related to project management, architecture, development practices, and historical context.

---

## 🏛️ System Architecture

### Architecture Components
- **[API Reference](./architecture/api-internal.md)**: This document provides a reference for the key APIs within the Systematic Review Assistant application, focusing on the service layer and other critical interfaces..
- **[Backend Architecture](./architecture/backend-architecture.md)**: ** Placeholder for backend service patterns.
- **[Components](./architecture/components.md)**: ** Placeholder for component documentation.
- **[Core Workflows](./architecture/core-workflows.md)**: This document provides details on core workflows.
- **[Data Models & Schemas](./architecture/data-models.md)**: ** Details the Pydantic schemas for API operations and LLM interactions.
- **[Database Schema](./architecture/database-schema.md)**: ** An overview of the database schema, including an ERD and table definitions.
- **[Environment Variables](./architecture/environment-vars.md)**: This document lists the environment variables used by the Systematic Review Assistant (SRA) application. These variables are crucial for configuring database connections, API keys for...
- **[External APIs](./architecture/external-apis.md)**: ** Placeholder for third-party API integration details.
- **[Frontend Architecture](./architecture/frontend-architecture.md)**: ** Placeholder for frontend patterns.
- **[Naming Conventions](./architecture/naming-conventions.md)**: This document establishes clear and consistent naming conventions for key entities, Pydantic schemas, database models, classes, functions, and other significant components within the Systematic Review...
- **[Project Structure](./architecture/unified-project-structure.md)**: ** Defines the organization of code to promote separation of concerns.
- **[REST API Specification](./architecture/rest-api-spec.md)**: ** Placeholder for the formal REST API specification.
- **[Technology Stack](./architecture/tech-stack.md)**: ** Specification of the technologies and libraries used in the project.
- **[UI/UX Specification](./architecture/ui-ux-spec.md)**: ** Placeholder for UI/UX specifications.

### Overviews
- **[Architecture Overview](./architecture.md)**: This document outlines the overall project architecture for the Systematic Review Assistant (SRA), including backend systems, shared services, and non-UI specific concerns.

---

## 💻 Development & Operations

### Benchmarking & Metrics
- **[Benchmark Sub-System Workflows](./benchmark/benchmark-workflows.md)**: This document contains the Mermaid sequence diagrams for the SRA Benchmarking Sub-System..
- **[Calculating Accuracy Metrics for Systematic Review Screening: Formulas and Applications](./sr-metrics.md)**: Systematic reviews require rigorous screening of literature to identify relevant studies. Evaluating the performance of this screening process—whether conducted manually or with automated tools—relies on...
- **[MPH AI Screening Prototype - Test Protocol](./benchmark/benchmark-protocol.md)**: homeless OR “insecure housing” OR “unstable housing” OR “precarious housing”.
- **[O3 Xml Criteria Homelessness](./benchmark/o3-xml-criteria-homelessness.md)**: <ScreeningCriteria version="2025-06-04"> <Domain name="Population"> <Include> <Rule id="POP1">Participants are <emphasis>currently homeless</emphasis> (roofless, houseless, emergency or other temporary accommodation, sofa-surfing).</Rule> </Include> <Exclude> <Rule id="POPX1">Participants are solely <emphasis>formerly</emphasis>...

### Core Development Guides
- **[A Guide to Advanced Prompt Engineering for AI Agents](./guides/prompt-engineering-guide.md)**: Prompt engineering is the systematic process of designing, refining, and optimizing input queries—known as prompts—to elicit desired, accurate, and relevant outputs from large language models...
- **[Coding Standards](./architecture/coding-standards.md)**: ** The Python style guide and project conventions.
- **[In app.py](./guides/streamlit-apptest-framework-guide.md)**: ﻿ Guide to UI Testing in Streamlit with its Built-in Framework.
- **[Testing Strategy](./architecture/testing-strategy.md)**: ** The comprehensive testing strategy for the application.

---

## 📋 Project & Product Management

### Product Requirements Documents (PRDs)
- **[Benchmarking Module PRD](./prd-benchmark-may.md)**: ** Requirements for the SRA benchmarking module.
- **[Epic 1: Search and Service Layer Stabilization](./prd/epic1-recovery-search-stabilization.md)**: Goal: Fully refactor and stabilize search.py, services.py, and repositories.py to correctly implement and utilize the SearchResult model and the new service layer architecture for PubMed...
- **[Epic 2: Resolver Agent Implementation and Integration](./prd/epic2-recovery-resolver-implementation.md)**: Goal: Complete the implementation and integration of the automated screening conflict resolver, including model updates, backend logic, and UI integration, ensuring it accurately processes conflicts...
- **[Epic 3: Comprehensive Testing and Database Integrity](./prd/epic3-recovery-testing-and-integrity.md)**: Urgent Directive from CTO: Postponed and deprioritized if app is functional after epics 1, 2. Benchmarking module is a P1 project (see project briefs in...
- **[Epic 4: SRA Benchmarking Module Implementation](./prd/epic4-sra-benchmarking-module.md)**: Status: Proposed Date Created: 2025-05-16 Last Updated: 2025-05-20 Related PRD(s): [docs/prd-benchmark-may.md](/docs/prd-benchmark-may.md) Architect: Architect Agent.
- **[Overall Product Requirements (PRD)](./prd.md)**: ** The main PRD for the MPH SR Prototype.
- **[Recovery Project PRD](./prd-recovery.md)**: ** Requirements for restoring the prototype to a stable state.
- **[Resolver Agent PRD](./prd-resolver.md)**: ** Requirements for the automated screening decision resolver.

### Project Briefs & Overviews
- **[Project Brief: MPH SR Prototype](./project-brief.md)**: Version: 1.0 Date: 2024-07-29 Prepared by: AIDE (AI Assistant) & User.
- **[Project Brief: MPH SR Prototype Stabilization, Refactoring Completion, and Resolver Implementation](./project-brief.recovery.md)**: Version: 1.0 Date: 2024-07-29 Prepared by: AIDE (AI Assistant) & User.

### User Stories
- **[Story 1.1: Ensure `SearchResult` Model Consistency and Correct Service Usage in `search.py`](./stories/1.1.story.md)**: Status: Done.
- **[Story 1.2: Refactor `SearchService` for PubMed Logic & Session Management](./stories/1.2.story.md)**: Status: Done.
- **[Story 1.3: Implement Specific Methods for `SearchResultRepository`, Define `SearchResultFilter` for Queries, and Verify Inherited CRUD](./stories/1.3.story.md)**: Status: Done.
- **[Story 1.4: End-to-End PubMed Search Workflow Stabilization](./stories/1.4.story.md)**: Status: Done.
- **[Story 1.5: Define and Align `SearchResult` Pydantic Schemas](./stories/1.5.story.md)**: Status: Done.
- **[Story 1.6: Enhance Protocol Page with Manual Inclusion Criteria and PICO Drafting](./stories/1.6.story.md)**: - As a systematic review researcher - I want to have the option to manually define inclusion criteria and then use AI to draft PICO...
- **[Story 2.1: Define and Align Core Pydantic Schemas and Setup Resolver Data Persistence](./stories/2.1.story.md)**: Status: Done.
- **[Story 2.2: Define and Implement Resolver Agent (Chain)](./stories/2.2.story.md)**: Status: Done.
- **[Story 2.3: Integrate Resolver into Screening Workflow](./stories/2.3.story.md)**: Status: Ready for Development.
- **[Story 2.4: Integrate Resolver Workflow and Update screen_abstracts.py UI & Logic](./stories/2.4.story.md)**: Status: SUPERSEDED by Story 2.4A.
- **[Story 2.4A: Update screen_abstracts.py UI Integration](./stories/2.4A.story.md)**: Status: Ready for Development Priority: MEDIUM - UI integration layer.
- **[Story 2.5: Refactor and Implement ScreeningService for API Alignment and Resolver Logic](./stories/2.5.story.md)**: Status: SUPERSEDED by Stories 2.5A and 2.5B.
- **[Story 2.5A: Refactor ScreeningService Core Methods (Foundation Task)](./stories/2.5A.story.md)**: Status: Ready for Development Priority: IMMEDIATE - Foundation for all other Epic 2 work.
- **[Story 2.5B: Implement ScreeningService Resolver Methods](./stories/2.5B.story.md)**: Status: Ready for Development Priority: HIGH - Core resolver functionality.
- **[Story 2.6: Refactor Screening Workflow for Service Layer Integration & UI Error Resolution](./stories/2.6.story.md)**: Status: Complete ✅.
- **[Story 4.10: Benchmark UI - Display AI Confidence Stats](./stories/4.10.story.md)**: Status: Draft (NOTE: this story is skipped and deprioritized).
- **[Story 4.11: Benchmark Results Export - Detailed Items](./stories/4.11.story.md)**: Status: Done.
- **[Story 4.12: Benchmark Results Export - Summary Metrics](./stories/4.12.story.md)**: Status: Done.
- **[Story 4.13: Adherence to Logging Standards for Benchmark Module](./stories/4.13.story.md)**: Status: Done.
- **[Story 4.14: Improve Benchmark Screening Prompts for Conservative Decision-Making](./stories/4.14.story.md)**: Status: Done.
- **[Story 4.1: Seed Benchmark Protocol Data](./stories/4.1.story.md)**: Status: Done.
- **[Story 4.2: Seed Benchmark Dataset](./stories/4.2.story.md)**: Status: Done.
- **[Story 4.3: Define `BenchmarkRun` Model and Schemas](./stories/4.3.story.md)**: Status: Done.
- **[Story 4.4: Define `BenchmarkResultItem` Model and Schemas](./stories/4.4.story.md)**: Status: Done.
- **[Story 4.5: Create Database Migration for Benchmark Tables](./stories/4.5.story.md)**: Status: Done.
- **[Story 4.6: Benchmark UI - Load and Display Protocol](./stories/4.6.story.md)**: Status: Done.
- **[Story 4.7: Benchmark UI - Trigger Benchmark Run & Process Items](./stories/4.7.story.md)**: Status: Done.
- **[Story 4.8: Automated Metrics Calculation and Persistence](./stories/4.8.story.md)**: Status: Done.
- **[Story 4.9: Benchmark UI - Display Summary Performance Metrics](./stories/4.9.story.md)**: Status: Review.
- **[Story Dependency Diagram](./stories/story-dependency-diagram.md)**: ** Visualizes the dependencies between user stories.

---

## 📜 Archives & Miscellaneous

### Archived Chat Logs
- **[1 Po Validation](./chats/1-po-validation.md)**: This completes my validation. Please let me know if you have any questions or require further clarification on any part of this review..
- **[Code Reusability and Implementation Reference](./chats/2.2.dev-session.md)**: Exported on 17/05/2025 at 12:41:17 BST from Cursor (0.50.4).
- **[Code Review for Story @2.2](./chats/2.2.sm-session.md)**: Exported on 17/05/2025 at 12:40:26 BST from Cursor (0.50.4).
- **[Development Workflow for Cross-Cutting Story](./chats/1.2.story.impl-session-full.md)**: Exported on 13/05/2025 at 8:40:58 BST from Cursor (0.50.1).
- **[Ready to Go!](./chats/1.1.story.impl-session.md)**: Exported on 13/05/2025 at 14:18:56 BST from Cursor (0.50.3).
- **[Real Code Review for Story 2.1](./chats/1.2.story.sm-session-full.md)**: Exported on 13/05/2025 at 9:25:42 BST from Cursor (0.50.1).
- **[Updates on Story 1.2 and SearchService Changes](./chats/1.2.story.refactoring-session.md)**: Exported on 13/05/2025 at 12:18:52 BST from Cursor (0.50.3).

### Other Documents
- **[PO Recovery Approval](./po-recovery-approval.md)**: ** The validation summary from the Product Owner regarding the recovery plan.