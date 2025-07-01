# Propuesta de mejora: Fine-tuning, RAG e interfaz para revisores médicos

Esta propuesta describe los cambios necesarios para alinear el proyecto con los requisitos de tener un sistema con capacidades de fine-tuning, recuperación aumentada (RAG) y una interfaz para que personal médico valide la información.

## Objetivos principales

1. **Integrar un pipeline de fine-tuning** de modelos de lenguaje para personalizar los agentes de revisión.
2. **Incorporar un subsistema de Recuperación Aumentada con Generación (RAG)** para mejorar la calidad de las respuestas y la deduplicación.
3. **Diseñar una interfaz para personal médico** que facilite la revisión y validación de resultados.

## Cambios propuestos

### 1. Arquitectura y backend

- Migrar la lógica actual de Streamlit a un backend FastAPI que exponga endpoints REST. Esto permitirá reutilizar la lógica de servicios con diferentes frontends (por ejemplo, Next.js o la misma app de Streamlit).
- Añadir módulos para la gestión de embeddings y almacenamiento vectorial. Se recomienda usar `pgvector` en PostgreSQL y `LangChain` o `LlamaIndex` para la capa de RAG.
- Definir tareas asíncronas (por ejemplo con Celery o RQ) para manejar el proceso de fine-tuning y la indexación de documentos.

### 2. Pipeline de fine-tuning

- Incorporar un módulo que permita
  - preparar datasets de entrenamiento a partir de los resultados de revisiones previas,
  - lanzar experimentos de fine-tuning (por ejemplo en Hugging Face) y
  - almacenar los modelos resultantes.
- Documentar los pasos para ejecutar el fine-tuning y registrar métricas de evaluación.

### 3. Subsistema RAG

- Implementar un flujo de recuperación de documentos relevantes mediante embeddings.
- Integrar la generación de respuestas apoyadas en los documentos recuperados, usando LangChain para la orquestación de prompts.
- Mantener un proceso de actualización periódica del índice para habilitar la "living review".

### 4. Interfaz para revisores médicos

- Diseñar una interfaz sencilla (puede partir de Streamlit o un frontend en Next.js) que permita:
  - Visualizar los artículos y resultados sugeridos por los agentes.
  - Marcar relevancia, aprobar o rechazar información.
  - Descargar informes generados automáticamente.
- Conectar esta interfaz a la API del backend para registrar el feedback y alimentar el proceso de fine-tuning.

### 5. Documentación y QA

- Completar la documentación de arquitectura en `docs/architecture/backend-architecture.md` y `docs/architecture/frontend-architecture.md` describiendo estos cambios.
- Añadir pruebas de integración y E2E que cubran el pipeline de fine-tuning, la consulta RAG y la interfaz de revisión.

## Siguientes pasos

1. Crear un roadmap detallado con tareas para cada módulo propuesto.
2. Priorizar la migración a FastAPI y la creación del subsistema de embeddings.
3. Desarrollar prototipos pequeños de fine-tuning y RAG antes de integrarlos completamente en la aplicación.

Estas acciones permitirán avanzar hacia un sistema de revisión sistemática automatizada que combine la potencia de modelos ajustados con un flujo de verificación amigable para expertos en salud.
