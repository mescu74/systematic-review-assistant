# MPH SR Prototype Product Requirements Document (PRD)

## Goals and Background Context

### Goals

- **Vision:** To create a best-in-class, AI-powered Systematic Review Assistant (the "MPH SR Prototype") that significantly reduces the time and effort required for conducting high-quality systematic reviews.
- Develop and validate an AI-driven T&A screening module that achieves performance comparable to or exceeding human expert reviewers.
- Implement features supporting the entire SR lifecycle, from protocol drafting to PRISMA flow chart generation.
- Integrate an always-present, context-aware AI chat assistant for intelligent guidance.
- Build a robust Benchmarking Module for performance evaluation against standard datasets like SYNERGY.
- Ensure a scalable, maintainable, and extensible architecture for future productization.
- Guarantee full transparency and auditability through comprehensive logging and LLM tracing.

### Background Context

Systematic reviews (SRs) are foundational for evidence-based practice but are incredibly time-consuming and labor-intensive, especially the literature screening phase. Current manual processes are prone to human error and can take months or even years. The "MPH SR Prototype" aims to address this by developing an AI-assisted platform to streamline the systematic review workflow. The initial focus is on revolutionizing the efficiency and accuracy of the title and abstract (T&A) screening process, while providing robust tools for performance evaluation and ensuring full auditability.

### Change Log

| Date       | Version | Description              | Author |
| :--------- | :------ | :----------------------- | :----- |
| 2024-07-30 | 1.1     | Migrated to v4 PRD Template | AIDE   |
| 2024-07-29 | 1.0     | Initial draft            | AIDE   |

## Requirements

### Functional

- **FR1:** The system shall provide AI-assisted guidance for developing a PICO(T/S) framework-based protocol.
- **FR2:** The system shall support searching multiple databases, with PubMed as the initial focus.
- **FR3:** The system shall include AI-assisted query building tools (keyword suggestion, Boolean logic).
- **FR4:** The system shall automatically retrieve search results and perform embedding-based deduplication.
- **FR5:** A dual-agent (Conservative and Comprehensive) AI system shall screen titles and abstracts.
- **FR6:** An automated Conflict Resolver agent shall handle disagreements between the primary screening agents.
- **FR7:** The UI shall allow users to review all AI decisions, including rationales.
- **FR8:** The Benchmarking Module shall calculate and display a full suite of performance metrics (Sensitivity, Specificity, F1, MCC, etc.).
- **FR9:** The Benchmarking Module shall support loading local datasets and the SYNERGY dataset for evaluation.
- **FR10:** The system shall automatically generate a PRISMA flow diagram based on screening numbers.
- **FR11:** A contextual AI chat assistant shall be available throughout the application.
- **FR12:** The system shall integrate with LangSmith for complete LLM call tracing.
- **FR13:** All user and AI system events must be logged to a Supabase database for a complete audit trail.
- **FR14:** The system must support user management and project organization.

### Non Functional

- **NFR1:** The system must be designed with a scalable and extensible architecture to support future growth.
- **NFR2:** All AI and human actions must be 100% auditable and reproducible.
- **NFR3:** The system must prioritize explainable AI (XAI) principles in its design.
- **NFR4:** The system must be capable of handling large datasets for search results and benchmarking.
- **NFR5:** The system must adhere to ethical AI considerations and ensure data privacy.

## User Interface Design Goals

### Overall UX Vision

A clean, intuitive, and professional interface that guides researchers through the systematic review process with clarity and efficiency. The UI should build trust by making AI actions transparent and easily verifiable.

### Key Interaction Paradigms

- Wizard-like guided workflows for complex tasks like protocol creation and search strategy.
- Dashboard-centric views for monitoring project progress and screening analytics.
- Integrated chat panel for seamless access to the AI assistant.

### Core Screens and Views

- Project Dashboard
- Protocol Definition Screen
- Search & Retrieval Module
- Abstract Screening Interface (with AI rationales)
- Benchmarking & Analytics Dashboard
- User/Project Management Pages

### Accessibility

WCAG 2.1 AA

### Branding

Modern, clean, and academic. Focus on readability and data visualization clarity.

### Target Device and Platforms

Web Responsive (Desktop-first)

## Technical Assumptions

### Repository Structure

Monorepo

### Service Architecture

- **Initial Prototype Stack:** Python, Streamlit, Supabase, LangChain, LangGraph.
- **Target Production Stack:**
    - Frontend: Next.js 15 (App Router), React 19, Tailwind 4, Shadcn/UI, Radix UI.
    - Backend: FastAPI, Python, Supabase, retaining LangChain and LangGraph for LLM orchestration.

### Testing requirements

- Unit, Integration, and E2E tests are required.
- Manual testing convenience methods should be available.

### Additional Technical Assumptions and Requests

- Utilize LangSmith for full LLM call tracing.
- Log all significant application activity to a dedicated Supabase table for a complete audit trail.

## Epics

- **Epic 1: Foundational Setup & Protocol Definition:** Establish the project infrastructure, user management, and the AI-assisted protocol development module.
- **Epic 2: Search, Retrieval & Deduplication:** Implement the multi-database search module, result retrieval, and automated deduplication.
- **Epic 3: AI-Powered Screening & Conflict Resolution:** Develop the dual-agent screening system, the conflict resolver, and the user-facing screening interface.
- **Epic 4: Benchmarking & Performance Analytics:** Build the comprehensive benchmarking module for performance evaluation against various datasets.
- **Epic 5: Auditability, Logging & PRISMA Generation:** Implement the end-to-end audit trail, logging infrastructure, and automated PRISMA diagram generation.

## Epic 1 Foundational Setup & Protocol Definition

Establish project infrastructure and the core module for users to define their systematic review protocol with AI assistance.

### Story 1.1 Project Initialization

As a developer,
I want to set up the Monorepo with CI/CD,
so that we have a stable foundation for development.

#### Acceptance Criteria

- 1: Monorepo structure is created.
- 2: Initial FastAPI backend and Next.js frontend apps are scaffolded.
- 3: Basic CI/CD pipeline for testing and linting is configured.

### Story 1.2 User & Project Management

As a researcher,
I want to create an account and organize my work into projects,
so that I can manage multiple systematic reviews securely.

#### Acceptance Criteria

- 1: Users can sign up, log in, and log out.
- 2: Users can create, view, and manage distinct projects.
- 3: Project data is stored securely in Supabase.

### Story 1.3 AI-Assisted Protocol Development

As a researcher,
I want the system to guide me through creating a PICO(T/S) protocol,
so that I can define the scope of my review accurately.

#### Acceptance Criteria

- 1: A form is available to input PICO(T/S) elements.
- 2: An AI assistant provides suggestions for refining each protocol element.
- 3: The completed protocol is saved and associated with the current project.

---
*(Further epics and stories would be detailed in this manner)*

## Checklist Results Report

[[LLM: This section will be populated after running the `pm-checklist` task.]]

## Next Steps

### Design Architect Prompt

[[LLM: Generate a prompt for the Design Architect based on this PRD.]]

### Architect Prompt

[[LLM: Generate a prompt for the Architect based on this PRD.]]
