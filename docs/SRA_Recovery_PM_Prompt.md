# PM Prompt: MPH SR Prototype Recovery Project

The MPH SR Prototype is currently in a critical state, requiring a focused recovery project to address an incomplete architectural refactoring (search functionality and service layer) and an unfinished feature (automated screening conflict resolver). This has resulted in a non-functional screening workflow, numerous errors, and data model inconsistencies.

We need a detailed PRD for a recovery project with three primary pillars:

1. **Stabilize Search & Service Layer (PubMed Focus):**
    * Tasks: Correct `SearchResult` model usage in `search.py`. Ensure `SearchService` correctly manages PubMed search logic and DB sessions. Repair `SearchResultRepository` with methods required by the service (e.g., for storing/deleting results). Resolve all related linter errors and make PubMed search fully functional.
2. **Complete Resolver Agent Implementation:**
    * Tasks: Implement the missing `resolve_screening_conflict` function. Update `SearchResult` model to include `final_decision` (with DB migration). Ensure correct data flow for the resolver LLM and storage of its decisions in `ScreeningResolution` and `SearchResult`. Integrate fully into `screen_abstracts.py` UI and workflow.
3. **Build Comprehensive & Safe Test Suite:**
    * Tasks: Develop unit tests for all services. Fix and expand unit tests for models and repositories. Complete integration tests for the resolver. Create an end-to-end integration test for PubMed search -> screening -> resolution. Document and enforce safe integration testing against the `sra_integration_test` DB.

This recovery is foundational. Once stable, the next planned major initiative is the development of a new benchmarking module (2-3 new UI pages, details to follow in a separate brief). For now, the priority is to make the MPH SR Prototype robust and complete the currently broken/unfinished work.
