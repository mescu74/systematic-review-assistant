# SR Assistant Prototype

Systematic Review Assistant Prototype.

## Getting Started

1. [Install uv](https://docs.astral.sh/uv/getting-started/installation/),
   [uv](https://github.com/astral-sh/uv) is a modern Python package manager written in Rust. It has a
   `pip` compatible interface under `uv pip` command and resolves dependencies up to 10x faster than
   [Poetry](https://python-poetry.org/). `uv` is to Python package management what
   [ruff](https://docs.astral.sh/ruff/) is to formatting and linting.

   - `brew install uv`
   - `pipx install uv`
   - `pip install uv`

2. Install project local Python 3.12 if not already present: `uv python install`.
   This doesn't install a global `python`, and is only used with `uv` commands and
   in the virtual environment. Takes a few seconds.

3. Create a venv and install the project: `uv sync`

4. Run the app stages:

   ```sh
   uv run streamlit run src/step1/app_step1.py
   uv run streamlit run src/step2/app_step2.py
   ```

   This will likely change if we move to a [multi-page app](https://docs.streamlit.io/develop/concepts/multipage-apps/page-and-navigation) with one entrypoint.

**Other useful commands:**

- `uv run ruff format .` - Formats the codebase using `ruff`.
- `uv run ruff check . --fix` - Checks the codebase for style issues using `ruff` and fixes them.
- `uv run pytest` - Runs tests using `pytest`.
