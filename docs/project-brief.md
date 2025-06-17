# Project Brief: MPH SR Prototype

**Version:** 1.0
**Date:** 2024-07-29
**Prepared by:** AIDE (AI Assistant) & User

## Introduction Problem Statement

Systematic reviews (SRs) are foundational for evidence-based practice but are incredibly time-consuming and labor-intensive, especially the literature screening phase. Current manual processes are prone to human error and can take months or even years. There's a pressing need for tools that can accelerate SRs while maintaining or improving rigor and transparency.

The "MPH SR Prototype" aims to address this by developing an AI-assisted platform to streamline the systematic review workflow, with an initial primary focus on revolutionizing the efficiency and accuracy of the title and abstract (T&A) screening process. The prototype will not only assist with screening but also provide tools for performance evaluation using robust statistical metrics, enabling comparison with human performance and other AI models, all while ensuring full auditability and transparency.

## Vision Goals

- **Vision:** To create a best-in-class, AI-powered Systematic Review Assistant (the "MPH SR Prototype") that significantly reduces the time and effort required for conducting high-quality systematic reviews. The prototype will initially focus on demonstrating superior performance in T&A screening, eventually evolving into a holistic SR platform. A core component will be a sophisticated **Benchmarking Module capable of evaluating the prototype's performance against large-scale, public datasets like the SYNERGY dataset ([https://github.com/asreview/synergy-dataset](https://github.com/asreview/synergy-dataset))**, thereby contributing to the broader research community on SR automation and fostering transparent, reproducible research in AI-assisted review methodologies, **backed by comprehensive audit trails and explainable AI principles.**
- **Primary Goals (for the full prototype, post-current recovery):**
  1. **Goal 1 (Advanced Screening Automation):** Develop and validate an AI-driven T&A screening module that achieves performance (measured by Sensitivity/Recall, Specificity, F1, MCC, etc.) comparable to or exceeding human expert reviewers, initially targeting a recall of >0.99.
  2. **Goal 2 (Comprehensive SR Workflow Support):** Implement features supporting the SR lifecycle, including AI-assisted protocol drafting (PICO/T/S elements), search strategy generation (block building, query expansion), multi-database searching (PubMed, Scopus, Embase), embedding-based deduplication, and PRISMA flow chart generation.
  3. **Goal 3 (Intelligent Assistance & Validation):** Integrate an always-present, context-aware AI chat assistant. Provide AI-driven validation for protocols and search queries, offering suggestions for improvement.
  4. **Goal 4 (Robust Benchmarking Capability):** Implement a Benchmarking Module to:
      - Calculate and display a comprehensive suite of screening performance metrics (Sensitivity, Specificity, PPV, NPV, Accuracy, F1 Score, MCC, Cohen\'s Kappa, PABAK, Likelihood Ratios).
      - Support evaluation against user-provided, human-annotated datasets.
      - Enable benchmarking against the SYNERGY dataset, including ingestion of its OpenAlex Work object format.
      - Allow comparison of different MPH SR Prototype configurations (LLMs, prompts, resolver strategies).
  5. **Goal 5 (Scalable & Extensible Architecture):** Design and build the prototype with a technology stack that can scale and be maintained for future development into a full-fledged product.
  6. **Goal 6 (Transparency & Auditability):** Implement comprehensive logging, LLM tracing (LangSmith), and data capture mechanisms to ensure a full audit trail of all AI and human actions, supporting explainable AI and reproducible research outputs.
- **Success Metrics (Initial Ideas for Full Prototype):**
    - Achieve target screening performance metrics (e.g., Recall > 0.99, high MCC) on standard datasets like SYNERGY.
    - Demonstrable reduction in time/effort for key SR tasks (e.g., screening, query generation) compared to manual methods.
    - Positive user feedback on usability, AI assistance quality, and trustworthiness of results.
    - Successful generation of PRISMA flow charts.
    - Benchmarking module provides insightful and actionable performance data.
    - Complete audit trails are exportable and verifiable.

## Target Audience Users

- Researchers, academics, and clinicians across various disciplines (medicine, social sciences, etc.) conducting systematic reviews.
- Organizations performing evidence synthesis (e.g., healthcare guideline developers, policy institutes).

## Key Features Scope High-Level Ideas for Full Prototype MVP

- **AI-Assisted Protocol Development:** Guided PICO(T/S) framework input, AI suggestions for refining protocol elements.
- **Advanced Search Module:**
    - Multi-database search (PubMed, Scopus, Embase - with PubMed as initial focus).
    - AI-assisted query building (keyword suggestion, Boolean logic, query expansion).
    - Automated query validation.
- **Efficient T&A Retrieval & Deduplication:** Fetching results, embedding-based duplicate removal.
- **High-Performance Dual-Agent Screening:**
    - Conservative and Comprehensive AI reviewers.
    - Automated Conflict Resolver agent.
    - User interface for reviewing decisions and rationales.
- **Screening Performance Analytics & Benchmarking Module:**
    - Calculation and display of confusion matrix, Sensitivity, Specificity, PPV, NPV, Accuracy, F1, MCC, Cohen\'s Kappa, PABAK, LR+/LR-.
    - Ability to load and process local datasets and the SYNERGY dataset.
    - Comparison of different screening configurations.
- **PRISMA Flow Diagram Generation:** Automated creation based on screening numbers.
- **Contextual AI Chat Assistant:** Integrated throughout the application.
- **Robust Auditing & Explainability Infrastructure:**
    - Integration with LangSmith for LLM call tracing.
    - Comprehensive logging of all user and AI system events to Supabase.
    - Functionality to export a complete audit trail for a review.
    - (Future) Dashboards for log analytics and system monitoring.
- **User Management & Project Organization.**

## Known Technical Constraints or Preferences

- **Initial Prototype Stack:** Python, Streamlit, Supabase, LangChain, LangGraph.
- **Target Production Stack:**
    - Frontend: Next.js 15 (App Router), React 19, Tailwind 4, Shadcn/UI, Radix UI.
    - Backend: FastAPI, Python, Supabase, **retaining LangChain and LangGraph for LLM orchestration.**
- **Comprehensive Tracing & Logging:** Utilize LangSmith for full LLM call tracing. Log all significant application activity and user/AI events to a dedicated Supabase table (or tables) to ensure a complete audit trail.
- **Auditability & Reproducibility:** Design the system to allow for the export of a 100% auditable trail of all AI and human activities undertaken during the systematic review process. This, combined with detailed logging and LangSmith traces, aims to ensure full transparency and reproducibility of results.
- **Explainable AI (XAI) & Analytics:** Strive for explainable AI decision-making processes wherever feasible. Plan for future capabilities to provide dashboards and analytics based on the logged activity and LLM traces.
- Must be capable of handling large datasets for search results and benchmarking.
- Prioritize ethical AI considerations, transparency in AI decision-making, and data privacy.

## Relevant Research Optional

- SYNERGY Dataset: [https://github.com/asreview/synergy-dataset](https://github.com/asreview/synergy-dataset)
- Literature on AI for systematic reviews, screening performance metrics, and LLM evaluation.
- (Internal: Link to the "MPH SR Prototype Stabilization, Refactoring Completion, and Resolver Implementation" brief for foundational work).

## PM Prompt for this full project vision

- (Content to be placed in `docs/MPH_SR_Prototype_Full_Vision_PM_Prompt.md`)
