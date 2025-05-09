# PM Prompt: Initiate Planning for Full "MPH SR Prototype" Vision

**Subject: Initiate Planning for Full "MPH SR Prototype" Vision**

Hi [PM Name],

Following the ongoing recovery and stabilization efforts for the "MPH SR Prototype" (detailed in the "MPH SR Prototype Stabilization, Refactoring Completion, and Resolver Implementation" brief), we need to begin comprehensive planning for the full vision of the prototype.

The "MPH SR Prototype" is envisioned as an AI-assisted platform to significantly accelerate systematic reviews, with an initial primary focus on achieving and validating high-performance T&A screening. Key long-term capabilities include:

*   **End-to-End SR Workflow Support:** Protocol drafting (AI-assisted PICO/T/S), advanced search (multi-database, AI query building/validation), T&A retrieval, embedding-based deduplication, and PRISMA chart generation.
*   **High-Performance Screening:** Leveraging dual AI agents (conservative/comprehensive) and an automated conflict resolver, with a strong emphasis on achieving measurable performance metrics (Sensitivity >0.95, high MCC, F1, etc.).
*   **Advanced Benchmarking Module:** A critical component for calculating a comprehensive suite of screening metrics (Sensitivity, Specificity, PPV, NPV, Accuracy, F1, MCC, Cohen\'s Kappa, PABAK, Likelihood Ratios). This module must support evaluation against user-provided datasets and, importantly, the **SYNERGY dataset ([https://github.com/asreview/synergy-dataset](https://github.com/asreview/synergy-dataset))**, including its OpenAlex Work object format. The goal is to enable robust comparison of MPH SR Prototype configurations and contribute to SR automation research.
*   **Intelligent Assistance:** An always-present, context-aware AI chat assistant.
*   **Transparency & Auditability:** Commitment to full transparency, auditability, and reproducibility through comprehensive logging to Supabase, LangSmith tracing for all LLM calls, explainable AI practices, and exportable audit trails.
*   **Technology Evolution:** While the current prototype uses Streamlit/LangChain/LangGraph, the long-term vision includes a transition to a Next.js 15/React 19 frontend and a FastAPI Python backend (retaining LangChain and LangGraph for core LLM workflow orchestration) for scalability and production readiness.

Please initiate the creation of a comprehensive PRD for the full "MPH SR Prototype". This PRD should:
1.  Detail all features outlined above.
2.  Define user stories and functional requirements for each capability.
3.  Consider the phased approach: building upon the stabilized foundation post-recovery, then tackling new modules like advanced search for Scopus/Embase, the full benchmarking suite with SYNERGY integration, etc.
4.  Lay the groundwork for the eventual technology stack transition.
5.  Detail requirements for comprehensive system logging, LangSmith integration for LLM tracing, and the generation of a complete audit trail for all review activities to support explainability and reproducibility.

The main "Project Brief: MPH SR Prototype" (located at `docs/MPH_SR_Prototype_Main_Project_Brief.md`) provides more detailed context on the overall vision. This effort will guide subsequent architectural design and epic/story creation. 