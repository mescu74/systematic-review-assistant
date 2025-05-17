# Coding Standards

The coding standards for this project are defined by the "Python Style Guide for AI Agents (Based on Ruff)". This document outlines the conventions for writing clean, modern, efficient, and maintainable Python code, specifically targeting Python 3.13. Adherence to these guidelines, which are based on Ruff's comprehensive rule set, will ensure consistency and high quality in the codebase, making it more accessible and manageable for both human developers and AI agents.

---

# Python Style Guide for AI Agents (Based on Ruff)

## Introduction

This document provides a concise style guide for AI agents writing Python code, specifically targeting Python 3.13. Its primary goal is to ensure the generation of clean, modern, efficient, and maintainable code that adheres to current best practices. The guidelines presented here are distilled from the comprehensive set of rules enforced by Ruff, an extremely fast Python linter and formatter. By following these conventions, AI agents can produce code that is less prone to errors, easier to understand, and aligns with the standards expected in contemporary Python development. This guide aims for clarity and conciseness, presenting the essence of Ruff's rules in a format suitable for automated code generation and review, following principles of effective technical communication.

## Code Style

This section outlines fundamental code style conventions based primarily on PEP 8, enforced by Ruff rules from `pycodestyle` (E, W), `flake8-quotes` (Q), `flake8-commas` (COM), and `pep8-naming` (N). Adhering to these styles improves code readability and consistency. Many of these rules are automatically fixable using Ruff's `--fix` option.

### Formatting and Linting

This project uses `ruff` with most rules enabled for linting and formatting. Refer to `pyproject.toml` for Ruff configuration. If you think a rule should be ignored, always confirm with the user before making changes to ruff configuration. Follow existing commenting conventions.

- After you've applied all edits to a file, you MUST run the following workflow:
  1. `uv run ruff check --fix edited/file.py`, then fix any remaining linter issues if related to your work. Let ruff do the fixing before attempting yourself.
  2. `uvx pyupgrade  --py312-plus --keep-runtime-typing edited/file.py`
  3. `uv run ruff format edited/file.py`

### Naming Conventions

Clear and consistent naming is vital for AI agents to understand and modify code. Follow PEP 8 naming conventions (N-series rules).

- **General:** Avoid single-letter names unless they are common iterators (e.g., `i`, `j`, `k` in loops) or in mathematical contexts where they are well-understood. Be descriptive.
- **Packages and Modules:** Use short, all-lowercase names. Underscores can be used if it improves readability (e.g., `my_package`).
- **Classes:** Use CapWords (CamelCase) convention (N801). Example: `MyClass`.
- **Type Variables:** Use CapWords, preferably short (N818). Example: `_T = TypeVar("_T")`, `_UserT = TypeVar("_UserT", bound=User)`.
- **Exceptions:** Use CapWords and end the name with `Error` (N806, N818 if it's a type variable). Example: `class MyCustomError(Exception): ...`.
- **Functions and Methods:** Use lowercase with words separated by underscores (snake_case) (N802, N807). Example: `def my_function(): ...`.
- **Method Arguments:**
  - `self` for the first argument to instance methods (N804).
  - `cls` for the first argument to class methods (N805).
- **Constants:** Use all uppercase with words separated by underscores (UPPER_SNAKE_CASE) (N818). Example: `MAX_OVERFLOW = 100`.
- **Variables:** Use lowercase with words separated by underscores (snake_case) (N803, N806). Example: `my_variable = 5`.
- **Function and Method Argument Names:** Use snake_case (N803).
- **Private/Internal:** Use a leading underscore (`_my_internal_var`) for internal variables/functions/methods. For name mangling (rarely needed), use two leading underscores (`__my_private_var`). Avoid using `__name__` style for non-magic methods (N807).

### Comments

Comments should be clear, concise, and primarily explain _why_ something is done, not _what_ is being done (the code itself should explain the what).

- Use TODO comments (TODO:, HACK:, PERF:, BUG:, WARN:, etc.). Indent following lines so they align with the first line.
- Do NOT add comments like "# add variable foobar", it is useless noise.
- **Docstrings:** Write docstrings for all public modules, functions, classes, and methods (D100-D107).
  - Follow Google docstring style. """Use imperative voice end first line in a period."""
  - For multi-line docstrings, the summary line should be on the first line, followed by a blank line, then the more detailed description (D210). The closing quotes should be on a line by themselves (D200).
  - Follow PEP 257 conventions. For example, use imperative mood for the first line: `"""Calculate the sum of two numbers."""` not `"""Calculates the sum..."""` (D400, D401).
  - Document parameters, return values, and any exceptions raised, especially for complex functions or APIs. Use a consistent format (e.g., Google style, NumPy style - D107 recommends a style, and specific rules like D402 ensure sections like `Args:`, `Returns:` are present if params are documented).
- **TODO Comments:** Use `# TODO:` for tasks that need to be done. Consider adding a reference to an issue tracker or your name (TD002, TD003, TD004).

## Modern Python Practices (Mostly `pyupgrade` - UP)

Utilize modern Python features for cleaner and more efficient code. Ruff's `pyupgrade` (UP) rules help enforce this.

- **Type Hinting:**
  - Use modern type hints (PEP 585, PEP 604). For example, use `list[int]` instead of `typing.List[int]` (UP006), and `int | str` instead of `typing.Union[int, str]` (UP007).
  - Use `collections.abc` for abstract base classes in type hints, e.g., `collections.abc.Sequence` instead of `typing.Sequence` (UP035).
  - Prefer `typing.Annotated` for an argument's type hint if it is present on a line by itself (UP040).
- **F-strings:** Prefer f-strings for string formatting over `%-formatting`, `str.format()`, or `string.Template` (UP032, FBT003).
  - Avoid unnecessary f-strings if they don't contain any expressions (UP030, RUF010).
  - Use `str()` calls inside f-strings for non-string variables where clarity is needed, though often implicit conversion is fine (RUF010).
- **Super Calls:** Use `super()` without arguments in Python 3 (UP008). Example: `super().__init__()`.
- **Class Definitions:** Do not inherit from `object` explicitly in Python 3 (UP004). Example: `class MyClass: ...` not `class MyClass(object): ...`.
- **Octal Literals:** Use `0o` prefix for octal literals (UP003). Example: `0o755`.
- **File I/O:** Use `with open(...)` for file handling to ensure files are closed properly. Avoid manual `f.close()` calls that aren't guaranteed.
- **Type Checking:** Use `isinstance()` and `issubclass()` for type checking instead of direct type comparisons like `type(obj) is MyClass` (PT001-PT026 cover various Pytest specific antipatterns, but E721 covers `type(obj) == MyClass` which is similar).
- **Unpacking:** Use iterable unpacking for assignments where appropriate. Example: `first, *rest = my_list`.
- **Dictionary Merging (Python 3.9+):** Use the `|` operator for merging dictionaries: `merged_dict = dict1 | dict2`.
- **Walrus Operator (Python 3.8+):** Use assignment expressions `:=` where they simplify logic, e.g., in `while` loops or comprehensions, but avoid overusing them if they reduce readability.
- **`datetime.timezone.utc`:** Use `datetime.timezone.utc` instead of `pytz.UTC` (UP017).
- **Yield/Return:** Avoid `yield` or `yield from` in `__init__` and `__new__` methods (ASYNC100).
- **Boolean Traps:** Avoid using `bool` as a type hint for `Optional` arguments with a default value of `True` or `False`. (PYI034)

## Error Prevention and Best Practices

These rules help prevent common bugs and improve code robustness.

- **Exception Handling:**
  - Be specific when catching exceptions. Avoid bare `except:` clauses (E722). Catch `Exception` if you need to catch most things, but preferably catch specific errors.
  - Use `raise ... from ...` to chain exceptions and preserve context (B904).
  - Avoid using `assert` for runtime checks that should be proper error handling (S101). `assert` is for debugging and can be disabled.
  - Do not use `logging.exception()` in an `except` block that is silencing or re-raising a different exception (LOG003).
  - Avoid `try...except...pass` (TRY002) and `try...finally` with `return`/`break`/`continue` in the `finally` block (TRY004, PLC0301).
  - Avoid raising `NotImplementedError` or `NotImplemented` (PIE794, PIE807). Prefer `raise TypeError` or `return NotImplemented` from binary magic methods.
- **Comparisons:**
  - Use `is` and `is not` for comparing to singletons like `None`, `True`, and `False` (E711, E712). Example: `if x is None: ...`.
  - Do not compare types directly, e.g. `type(a) == type(b)` (E721). Use `isinstance()` or `issubclass()`.
  - Avoid "yoda conditions" (SIM300): `if "value" == variable:` should be `if variable == "value":`.
  - Do not use `len(sequence) == 0` or `len(sequence) != 0`. Use `if not sequence:` or `if sequence:` (PLC1701, RUF004).
- **Mutable Default Arguments:** Do not use mutable default arguments like `[]` or `{}` in function definitions (B006, B008). Initialize them to `None` and create the mutable object inside the function.
  ```python
  def my_func(param=None):
      if param is None:
          param = []
      # ...
  ```
- **Comprehensions:**
  - Use comprehensions (list, dict, set) where they are more readable than `map()` or `filter()` or explicit loops.
  - Avoid unnecessary list comprehensions if a generator would suffice, especially for large datasets.
  - Be cautious with side effects in comprehensions.
- **`isort` / Import Order:** Ensure imports are sorted correctly (I001). Ruff can do this.
- **Unused Code:** Remove unused variables, imports, and functions (F841, F401, etc.). Ruff will report these.
- **Complexity:**
  - Keep functions and methods short and focused on a single task.
  - Avoid overly nested code blocks. Refactor using helper functions or by inverting conditions.
  - Ruff includes McCabe complexity checking (C901). Aim to keep complexity low.
- **Boolean Zen (SIM101, SIM103, SIM108, SIM109, SIM110, SIM111, SIM117, SIM118):**
  - Simplify boolean expressions. For example, `if cond: return True else: return False` should be `return cond`.
  - `if not a:` is preferred over `if a == False:`.
  - `if a:` is preferred over `if a == True:`.
  - Use `a and b` instead of `if a: if b: ...`.
  - Use `a or b` instead of ternary for default values if simpler: `result = value or default_value`.
- **Context Managers:** Use `with` statement for resources that need to be managed (files, locks, database connections). (PERF401, PERF402 recommend `try...finally` over `with` for `threading.Lock` in specific performance-critical contexts if the lock is already acquired).
- **flake8-bugbear (B series):** Many `B` rules catch common bugs. Examples:
  - Do not call `getattr` with a constant string (B009). Use direct attribute access.
  - Do not call `setattr` with a constant string (B010). Use direct attribute assignment.
  - Do not use `zip()` without `strict=True` if iterables may have different lengths (B905).
  - Avoid `try`-`except`-`else` that can be simplified by moving `else` logic outside (B012).
- **flake8-simplify (SIM series):** Many `SIM` rules help simplify code. Examples:
  - Use `key in dict` instead of `key in dict.keys()` (SIM118).
  - Use `dict.get(key, default)` instead of `if key in dict: val = dict[key] else: val = default` (SIM106).
  - Combine `isinstance` calls: `isinstance(obj, (A, B))` instead of `isinstance(obj, A) or isinstance(obj, B)` (SIM102).
- **Performance (PERF series):**
  - Avoid unnecessary list comprehensions for `sum`, `min`, `max` if a generator expression would work. (PERF401)
  - When checking for list emptiness, prefer `if not my_list:` over `if len(my_list) == 0:`. (RUF004)
  - In class `__slots__`, prefer `tuple` over `list` for `__slots__` definition. (PERF101)
  - Avoid `isinstance` checks in `try-except-else` blocks. (PERF201)
- **Security (S series - via `flake8-bandit`):**
  - Avoid `eval()` (S307), `exec()` (S102), `pickle` (S301, S302), `shelve` (S308), `subprocess.run` with `shell=True` without careful sanitization (S602, S603).
  - Be careful with XML parsing (S310-S320, S400-S413), ensure secure parsers are used if processing untrusted XML.
  - Do not use weak cryptographic functions (e.g., MD5, SHA1 for hashing if security is critical - S324).
  - Ensure temporary files are created securely (S306).

## Pydantic Best Practices

- Pydantic models are referred to as "schemas" in this project. Database models are referred to as "models".
- **Use Model.model_validate(obj)**:
  - **DO NOT** use `Model(**dict)`, it will cause Pyright type validation errors and bypasses proper Pydantic validation.
  - **ALWAYS** use `Model.model_validate(obj)` instead.
- **Schemas should inherit from `sra_assistant.core.schemas.BaseSchema`:**
  - _Rationale:_ `BaseSchema` configures `model_config` with settings we want to use for all schemas unless there is a special reason not to. It's also fine to define a `FooBase(BaseSchema)` and `Bar(FooBase)`.

## AI Agent Specific Best Practices

- **Clarity for LLMs:**
  - Be explicit. Avoid overly complex "one-liners" if a multi-line version is clearer.
  - Ensure type hints are comprehensive and correct. This is crucial for LLMs to understand data structures.
  - Use descriptive variable and function names, even if they are slightly longer.
- **Modularity:** Break down complex tasks into smaller, well-defined functions or classes. This makes it easier for an AI agent to understand, modify, or generate specific pieces of functionality.
- **Idempotency:** Where applicable, design functions to be idempotent (calling them multiple times with the same input produces the same result without unintended side effects).
- **State Management:** Be explicit about state. Avoid hidden state or reliance on global variables where possible. Pass state explicitly as arguments or manage it within class instances.
- **Error Reporting:** Provide clear error messages that can help an AI (or human) understand what went wrong.
- **Configuration:** Externalize configuration (e.g., API keys, model parameters) rather than hardcoding it. Use Pydantic for settings management if appropriate.

## Tooling: Ruff Configuration (`pyproject.toml`)

A typical Ruff configuration in `pyproject.toml` might look like this. This selects a broad set of useful rules.

```toml
[tool.ruff]
# Python version to target for features and syntax
target-version = "py313"

# Maximum line length
line-length = 88 # Or 79

# Select rule codes to enable.
# See https://docs.astral.sh/ruff/rules/ for all rules.
# Using "ALL" and then ignoring is often too noisy. A curated set is better.
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # Pyflakes
    "C90", # McCabe complexity
    "I",   # isort
    "N",   # pep8-naming
    "D",   # pydocstyle
    "UP",  # pyupgrade
    "ANN", # flake8-annotations (ANN101, ANN102, ANN401 are often good to have)
    "S",   # flake8-bandit (security)
    "BLE", # flake8-blind-except
    "B",   # flake8-bugbear
    "A",   # flake8-builtins
    "COM", # flake8-commas
    "DTZ", # flake8-datetimez
    "T10", # flake8-debugger (breakpoint)
    "EM",  # flake8-errmsg
    "EXE", # flake8-executable
    "ISC", # flake8-implicit-str-concat
    "ICN", # flake8-import-conventions
    "G",   # flake8-logging-format
    "INP", # flake8-no-pep420
    "PIE", # flake8-pie
    "T20", # flake8-print
    "PYI", # flake8-pyi
    "PT",  # flake8-pytest-style
    "Q",   # flake8-quotes
    "RET", # flake8-return
    "SLF", # flake8-self
    "SLOT",# flake8-slots
    "SIM", # flake8-simplify
    "TID", # flake8-tidy-imports
    "TCH", # flake8-type-checking
    "ARG", # flake8-unused-arguments
    "PTH", # flake8-use-pathlib
    "PERF",# Perflint
    "RUF", # Ruff-specific rules
    # "PL",  # Pylint (use specific PLxxx codes, can be very opinionated)
    # "TRY", # tryceratops
    # "LOG"  # flake8-logging
]

ignore = [
    "D100", "D101", "D102", "D103", "D104", "D105", "D106", "D107", # Ignore missing docstrings for now
    "ANN101", # Missing type annotation for `self` in method
    "ANN102", # Missing type annotation for `cls` in classmethod
    "ANN401", # Dynamically typed expressions (typing.Any) are disallowed
    "S101",   # Use of `assert` detected
    "S603",   # `subprocess` call: check for execution of untrusted input
    "S607",   # Starting a process with a partial path
    "BLE001", # Do not catch blind exception: `Exception`
    "EM101",  # Exception must not use a string literal, assign to variable first
    "EM102",  # Exception must not use an f-string
    # Add other specific rules to ignore as needed based on project decisions
]

[tool.ruff.pydocstyle]
convention = "google" # Or "numpy", "pep257"

[tool.ruff.isort]
known-first-party = ["sr_assistant"] # Replace with your project's src directory name
force-sort-within-sections = true
combine-as-imports = true

# If using per-file-ignores for specific directories (e.g., tests)
# [tool.ruff.per-file-ignores]
# "tests/*" = ["S101", "ANN401", "D103"]
# "**/__init__.py" = ["F401"]
```

## General Python Best Practices

- **Builtins:** Avoid shadowing built-in names like `list`, `dict`, `id` (A001, A002, A003).
- **Datetime:** Use timezone-aware datetimes where appropriate. Avoid naive `datetime.now()` or `datetime.utcnow()` (DTZ001-DTZ007).
- **Debugging:** Remove `print` statements (T201) and `breakpoint()` calls (T100) before committing/deploying.
- **TODOs:** Format `# TODO` comments consistently, potentially including author or issue links (TD001-TD007).
- **Unused Arguments:** Remove unused function/method arguments unless they are required by an interface (e.g., overriding a method) or explicitly marked as unused (e.g., `_unused_arg`) (ARG001-ARG005).
- **Implicit String Concat:** Use `+` or `str.join` instead of implicit string literal concatenation across lines if it harms readability (ISC001, ISC002).
- **Logging:** Avoid common logging mistakes like using string formatting instead of lazy logging arguments (LOG001-LOG009), or incorrect logging levels (G001-G010).
- **`__slots__`:** Consider using `__slots__` on classes where memory usage is critical and dynamic attribute assignment is not needed, but avoid on subclasses of `str` or `tuple` (SLOT000-SLOT002).
- **Private Member Access:** Avoid accessing private members (e.g., `_private_var` or `__mangled_var`) from outside their class (B018 for `except BaseException`, related to visibility). This is more of a convention but good practice.
- **Path Handling:** Prefer `pathlib.Path` over `os.path` for path manipulations (PTH series).

## Change Log

| Change        | Date       | Version | Description   | Author          |
| ------------- | ---------- | ------- | ------------- | --------------- |
| Initial draft | 2025-05-09 | 0.1     | Initial draft | Architect Agent |
| ...           | ...        | ...     | ...           | ...             |
