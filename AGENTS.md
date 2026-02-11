# AGENTS.md

This file provides context and instructions for AI agents working on the `gocam-py` project.

## Project Context

- **Name**: `gocam` (GO-CAM Data Model Python)
- **Purpose**: A Python implementation of the Gene Ontology Causal Activity Model (GO-CAM) schema.
- **Tech Stack**: Python (>=3.10), Pydantic v2, LinkML, NetworkX.
- **Key Features**:
    - Pydantic datamodels for GO-CAM.
    - Translation between GO-CAM and other formats (Minerva models, CX2, NetworkX).
    - CLI for model conversion.

## Setup Commands

- **Install only base dependencies**: `uv sync --no-dev`
- **Install with extras**: `uv sync --no-dev --all-extras`
- **Full development environment**: `uv sync --all-extras`

## Testing Instructions

- **Run Python tests**: `make test-python`
- **Run all tests (including Model examples)**: `make test`
- **Note**: Pytest is configured to use `importlib` mode.

## Code Style & Conventions

- **Linting**: Use `make lint-python` to check code style and formatting. Autofix issues with
  `make lint-fix-python`.
- **Type Checking**: Use `make type-check`
- **Exclusions**:
    - `src/gocam/datamodel/gocam.py` is a **generated file**. Do NOT include it in code reviews or
      attempt to manually edit it.
    - The `project/` directory contains LinkML artifacts and should generally be treated as
      read-only/generated.

## Pydantic Generation

- The `src/gocam/datamodel/gocam.py` file is generated from the LinkML schema in
  `src/gocam/schema/gocam.yaml`.
- To regenerate the Pydantic models after changes to the schema, run:
  `make src/gocam/datamodel/gocam.py`

## PR Instructions

- Ensure all tests pass (`make test-python`) before submitting.
- Ensure linting passes (`make lint-python`).
- Ensure type checking passes (`make type-check`).
- Do not commit changes to generated files unless you are intentionally re-running the generator.
