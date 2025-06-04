#!/usr/bin/env bash

# Copyright 2025 Gareth Morgan <garethm@lifework.health>.
# SPDX-License-Identifier: MIT

uv run ruff check --fix src/ tests/ tools/ >/dev/null 2>&1
uv run pyupgrade --py312-plus --keep-runtime-typing src/ tests/ tools/ >/dev/null 2>&1
uv run ruff format src/ tests/ tools/ >/dev/null 2>&1
