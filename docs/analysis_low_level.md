# Análisis de código a bajo nivel


## Estructura general

- El proyecto se organiza bajo `src/sr_assistant` con submódulos para la aplicación principal (`app`), un módulo de benchmarking y ejemplos antiguos (`step1`, `step2`).
- La aplicación principal se ejecuta con Streamlit (`app/main.py`) y gestiona la navegación y autenticación mediante `st.session_state`.
- No se observa un backend basado en FastAPI; toda la lógica está embebida en Streamlit.

## Componentes principales

### Servicios y acceso a datos

- Los modelos y repositorios se implementan con `SQLModel` y `SQLAlchemy`. La configuración de la base de datos se realiza en `database.py`, donde se define un `engine` síncrono y asíncrono.
- `services.py` contiene la capa de servicios. Destaca `SearchService.search_pubmed_and_store_results` que consulta PubMed usando `Bio.Entrez`. Requiere las variables de entorno `NCBI_EMAIL` y `NCBI_API_KEY`.
- Varias funciones de servicio tienen comentarios `TODO` (por ejemplo, resolución de conflictos y servicios de revisión), lo que indica que hay partes incompletas.

### Agentes y orquestación

- En `agents/` se definen cadenas de LangChain para realizar el cribado de abstracts. Se emplean modelos de OpenAI y Gemini.
- El repositorio cuenta con un módulo de benchmarking (página `benchmark_run_page.py` y documentación en `docs/benchmark`).

### Documentación

- El README describe los comandos `make` para instalar dependencias, formatear código, ejecutar tests y otros pasos de desarrollo.
- `docs/architecture/tech-stack.md` lista las tecnologías utilizadas (Streamlit, SQLModel, LangChain, etc.) y menciona extensiones como `pgvector` para deduplicación futura.
- Archivos como `docs/architecture/backend-architecture.md` solo contienen el título y no detallan la arquitectura.

## Hallazgos y problemas potenciales

1. **Dependencia de variables de entorno**: en la búsqueda de PubMed se leen directamente `NCBI_EMAIL` y `NCBI_API_KEY` con `os.getenv`. Si faltan, se lanza una excepción. Sería preferible centralizar la configuración con pydantic `BaseSettings` (ya se usa en `config.py`).
2. **Falta de separación entre backend y frontend**: actualmente todo se ejecuta en Streamlit. No hay un backend independiente ni API REST para otras interfaces.
3. **Código heredado sin usar**: los directorios `step1` y `step2` contienen scripts de versiones anteriores, lo que puede generar confusión.
4. **Documentación parcial**: algunos archivos de arquitectura están vacíos, dificultando la comprensión global del sistema.
5. **Pendientes en servicios**: existen funciones marcadas con `TODO` en `services.py` y otros módulos, indicando que no están terminadas.
6. **Flujo de pruebas**: el proyecto incluye `pytest` y configuraciones de `pre-commit`, pero las instrucciones para ejecutar las pruebas pueden requerir dependencias adicionales.

## Recomendaciones generales

- Completar la documentación de arquitectura y clarificar la evolución planificada del proyecto.
- Evaluar la creación de un backend separado (por ejemplo, FastAPI) para exponer APIs y mantener la lógica desacoplada de la interfaz de Streamlit.
- Eliminar o archivar los directorios `step1` y `step2` si ya no son necesarios.
- Unificar la gestión de configuración mediante el módulo de settings y reducir el uso directo de `os.getenv`.
- Implementar las funciones pendientes en `services.py` para que los flujos de revisión estén completos.

