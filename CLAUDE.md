# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

This project explores LLPS (Liquid-Liquid Phase Separation) propensity in human membrane proteins. The primary workflow is Jupyter notebooks for analysis and discovery; the end goal is a simple browser-based dashboard (`dashboard/`) that non-coders can use to explore the data without running Python.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

Python ≥ 3.11 required.

## Key Commands

```bash
# Run analysis notebooks in sequence
jupyter notebook

# Run tests
pytest tests/

# Preview the Shinylive dashboard locally
shinylive export dashboard/ /tmp/llps_dashboard
python3 -m http.server --directory /tmp/llps_dashboard 8008

# Run the server-based Shiny dashboard
shiny run scripts/shiny_app.py --reload --port 8000

# Pre-cache STRING interactions (for offline/slow-network use)
python scripts/analysis/generate_string_cache.py --threshold 0.7 --score 700
```

## Architecture

### Primary workflow: notebooks → `llps/` package → `results/`

The numbered notebooks (`01`–`07`) are the main analysis entry point. They import from `llps_functions` (a thin backward-compat shim over the `llps/` package) and write CSVs/JSONs to `results/`.

The `llps/` package modules:

| Module | Responsibility |
|---|---|
| `data.py` | Load XLSX, classify proteins into High/Medium/Low pLLPS tiers |
| `location.py` | Parse UniProt subcellular location strings |
| `functional.py` | Classify proteins into functional groups via `data/functional_classification_terms.yaml` |
| `string_api.py` | Query STRING REST API in batches; cache to `data/string_cache_{threshold}.json` |
| `network.py` | Build NetworkX graphs; compute topology/enrichment metrics |
| `enrichment.py` | Chi-squared enrichment tests (binary or 3×3 High/Medium/Low matrix) |
| `visualization.py` | Seaborn heatmaps, matplotlib network reports |
| `io.py` | JSON/CSV/pickle save-load for results; STRING cache management |
| `constants.py` | STRING API URLs; `StringQueryConfig` dataclass |

### Dashboard: `dashboard/app.py`

The Shinylive app runs entirely in the browser (Pyodide/WASM). It **cannot import from the `llps/` package**, so it duplicates the location/function parsing logic inline. It reads `dashboard/full_dataset.csv` (or `sample_data.csv` as fallback) and provides filters, histograms, scatter plots, and CSV export.

The dashboard is auto-deployed to GitHub Pages on push to `main` (`.github/workflows/deploy-dashboard.yml`), but only when files under `dashboard/` change.

### Dashboard data file

`dashboard/full_dataset.csv` is the dataset the dashboard reads. When analysis results change, regenerate it from notebook 01 or copy `results/full_dataset.csv`. The dashboard and the `llps/` package must stay in sync on column names (`Entry`, `Entry name`, `Protein names`, `p(LLPS)`, `Length`, `Function [CC]`, `Subcellular location [CC]`).

## Key Conventions

**`llps_functions.py` is a shim.** All notebooks use `from llps_functions import X`. Don't add logic there; add it to the relevant `llps/` module and re-export from `llps/__init__.py`.

**Functional classification is YAML-driven.** Regex rules for all protein groups live in `data/functional_classification_terms.yaml`. Edit that file to change how proteins are categorised; the dashboard has its own hardcoded copy of the same rules (`FUNCTION_CATEGORIES` dict in `dashboard/app.py`) that must be kept in sync.

**STRING identity mapping.** STRING uses gene names; the dataset uses UniProt accessions. Gene names are extracted from `Entry name` by splitting on `_` (format: `GENENAME_HUMAN`).

**STRING caching.** `fetch_string_interactions()` checks for `data/string_cache_{score_threshold}.json` before hitting the network. Rate-limit hits (HTTP 429) trigger a 30 s retry sleep.

**`networkx` is a soft dependency.** `network.py` sets `_HAS_NETWORKX = False` if the import fails; `analyze_network()` raises `ImportError` only when called. Tests use `pytest.importorskip("networkx")`.

**`results/` is committed.** Pre-computed outputs are checked in so collaborators can explore without running the pipeline.

**`deprecated/` is archived.** Do not reference notebooks there from active code.
