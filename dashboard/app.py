"""
LLPS Protein Data Explorer — Shinylive Dashboard

A static, browser-based Shiny dashboard for exploring Liquid-Liquid Phase
Separation (LLPS) protein data.  No server required: everything runs via
Pyodide (WebAssembly Python).

Usage (local preview):
    shinylive export dashboard/ /tmp/dashboard_export
    python3 -m http.server --directory /tmp/dashboard_export 8008
    open http://localhost:8008

Data files (place in the dashboard/ directory so they are bundled on export):
    full_dataset.csv  — full protein dataset (preferred default)
    sample_data.csv   — small fallback sample
"""

from __future__ import annotations

import base64
import json
import re
from pathlib import Path
from typing import Any

import altair as alt
import numpy as np
import pandas as pd
from shiny import App, reactive, render, ui, req

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Ordered list of (keywords, canonical_name) for subcellular location mapping.
# The first matching keyword wins for each canonical category; a location string
# can map to multiple canonical categories.
_LOCATION_CANONICAL: list[tuple[list[str], str]] = [
    (["plasma membrane", "cell membrane"], "Plasma membrane"),
    (["nucleus", "nucleoplasm", "nucleolus", "nuclear pore", "nuclear envelope"], "Nucleus"),
    (["cytoplasm", "cytosol"], "Cytoplasm"),
    (["endoplasmic reticulum"], "ER"),
    (["golgi apparatus", "golgi"], "Golgi apparatus"),
    (["mitochondri"], "Mitochondrion"),
    (["lysosom"], "Lysosome"),
    (["early endosome", "late endosome", "endosome"], "Endosome"),
    (["peroxisom"], "Peroxisome"),
    (["synapse", "presynaptic", "postsynaptic", "synaptic vesicle"], "Synapse"),
    (["cell junction", "tight junction", "adherens junction", "gap junction"], "Cell junction"),
    (["secreted", "extracellular space", "extracellular region"], "Extracellular"),
    (["cell surface"], "Cell surface"),
    (["stress granule"], "Stress granule"),
    (["p-body", "processing body"], "P-body"),
    (["lipid droplet"], "Lipid droplet"),
    (["autophagosom", "phagosom"], "Autophagosome"),
    (["vesicle"], "Vesicle"),
    (["centrosom", "centriole", "spindle pole"], "Centrosome"),
    (["axon", "dendrite", "neurite", "growth cone"], "Neuronal process"),
    (["cilium", "flagellum", "basal body"], "Cilium"),
    (["chromatin", "chromosome", "kinetochore"], "Chromatin/Chromosome"),
    (["focal adhesion", "focal contact"], "Focal adhesion"),
    (["sarcomere", "z disc", "myofibril"], "Sarcomere"),
    (["exosom"], "Exosome"),
]

# Functional categories matching the YAML-driven classification used by the
# llps/ package. Patterns are applied to the combined Function [CC] + Protein
# names text (lowercased, evidence codes stripped).
FUNCTION_CATEGORIES: dict[str, list[str]] = {
    "Ion channel": [
        r"\bion\s*channel\b", r"\bsodium\s+channel\b", r"\bpotassium\s+channel\b",
        r"\bcalcium\s+channel\b", r"\bchloride\s+channel\b", r"\bchannel\s+protein\b",
        r"\bvoltage.gated\b", r"\bligand.gated\b", r"\bcation\s+channel\b",
        r"\banion\s+channel\b", r"\baquaporin\b", r"\bgap\s+junction\b",
        r"\bconnexin\b", r"\bpannexin\b",
    ],
    "GPCR": [
        r"\bg\s+protein.coupled\s+receptor\b", r"\bgpcr\b",
        r"\b7.transmembrane\b", r"\b7tm\b", r"\bmetabotropic\s+receptor\b",
        r"\bserpentine\s+receptor\b",
    ],
    "Receptor tyrosine kinase": [
        r"\breceptor\s+tyrosine\s+kinase\b", r"\brtk\b",
        r"\btyrosine.protein\s+kinase.*receptor\b",
        r"\begf\s+receptor\b", r"\bfgf\s+receptor\b", r"\bpdgf\s+receptor\b",
        r"\bvegf\s+receptor\b", r"\binsulin\s+receptor\b",
    ],
    "Transporter": [
        r"\btransporter\b", r"\babc\s+transporter\b", r"\bsolute\s+carrier\b",
        r"\bslc\b", r"\bsymporter\b", r"\bantiporter\b", r"\buniporter\b",
        r"\bpermease\b", r"\bcarrier\s+protein\b", r"\bimporter\b", r"\bexporter\b",
    ],
    "Kinase": [
        r"\bkinase\b", r"protein\s+kinase", r"\bphosphorylates\b",
        r"\bcdk\b", r"\bcyclin.dependent\s+kinase\b", r"\bmapk\b",
    ],
    "Phosphatase": [
        r"\bphosphatase\b", r"\bprotein\s+phosphatase\b", r"\bdephosphorylates\b",
        r"\bpten\b",
    ],
    "Protease": [
        r"\bprotease\b", r"\bpeptidase\b", r"\bcaspase\b",
        r"\bproteolytic\b", r"\bmetalloprotease\b", r"\bserine\s+protease\b",
        r"\bcysteine\s+protease\b", r"\bprotein\s+cleavage\b",
    ],
    "Ligase": [
        r"\bligase\b", r"\bubiquitin\s+ligase\b", r"\be3\s+ligase\b",
        r"\bubiquitination\b", r"\bsumo\s+ligase\b",
    ],
    "Synthetase": [
        r"\bsynthetase\b", r"\bsynthase\b", r"\baminoacyl.trna\s+synthetase\b",
    ],
    "Transferase": [
        r"\btransferase\b", r"\bmethyltransferase\b", r"\bacetyltransferase\b",
        r"\bglycosyltransferase\b", r"\bphosphotransferase\b",
    ],
    "Oxidoreductase": [
        r"\boxidoreductase\b", r"\bdehydrogenase\b", r"\boxidase\b",
        r"\breductase\b", r"\bperoxidase\b", r"\bcytochrome\b",
    ],
    "Hydrolase": [
        r"\bhydrolase\b", r"\besterase\b", r"\blipase\b",
        r"\baTPase\b", r"\bhydrolytic\b",
    ],
    "Receptor": [
        r"\breceptor\b", r"\breceptor\s+activity\b",
    ],
    "Nuclear receptor": [
        r"\bnuclear\s+receptor\b", r"\bsteroid.*receptor\b",
        r"\bestrogen\s+receptor\b", r"\bandrogen\s+receptor\b",
        r"\bretinoic\s+acid\s+receptor\b", r"\bvitamin\s+d\s+receptor\b",
    ],
    "Transcription factor": [
        r"\btranscription\s+factor\b", r"\btranscriptional.*regulator\b",
        r"\btranscriptional.*activator\b", r"\btranscriptional.*repressor\b",
        r"\bdna.binding.*transcription\b",
    ],
    "GTPase": [
        r"\bgtpase\b", r"\bras\s+protein\b", r"\brho\s+protein\b",
        r"\brab\s+protein\b", r"\bguanine\s+nucleotide.binding\b",
    ],
    "Structural": [
        r"\bcytoskeleton\b", r"\bactin\b", r"\btubulin\b",
        r"\bintermediate\s+filament\b", r"\bcollagen\b", r"\blaminin\b",
        r"\bfibronectin\b", r"\bspectrin\b", r"\bscaffold\b",
    ],
    "Adhesion": [
        r"\badhesion\b", r"\bcadherin\b", r"\bintegrin\b",
        r"\bselectin\b", r"\bcell.cell.*junction\b", r"\bcell.*adhesion\b",
    ],
    "Chaperone": [
        r"\bchaperone\b", r"\bheat\s+shock\s+protein\b", r"\bhsp\b",
        r"\bprotein\s+folding\b", r"\bchaperonin\b", r"\bco.chaperone\b",
    ],
    "RNA processing": [
        r"\brna.binding\b", r"\brna\s+binding\b", r"\brna\s+processing\b",
        r"\brna\s+splicing\b", r"\bpre.mrna.*splicing\b", r"\bsplicing\s+factor\b",
        r"\bribosom\w+\b", r"\brna\s+helicase\b", r"\btranslation.*factor\b",
        r"\btranslational.*regulator\b",
    ],
    "DNA repair": [
        r"\bdna\s+repair\b", r"\bdna\s+damage\b", r"\bdna\s+replication\b",
        r"\bhomologous\s+recombination\b", r"\bgenome.*stability\b",
        r"\bdna\s+damage\s+response\b", r"\bdouble.strand\s+break\b",
    ],
    "Chromatin remodeling": [
        r"\bchromatin\s+remodel\w*\b", r"\bhistone.*methyltransferase\b",
        r"\bhistone.*acetyltransferase\b", r"\bhistone.*deacetylase\b",
        r"\bhistone.*demethylase\b", r"\bepigenetic\b", r"\bdna\s+methyltransferase\b",
        r"\bnucleosome.*remodeling\b",
    ],
}

_PLLPS_COLOR_SCALE = alt.Scale(
    domain=["High", "Medium", "Low"],
    range=["#e74c3c", "#f39c12", "#3498db"],
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _parse_location(raw: str | None) -> list[str]:
    """Map a UniProt subcellular location string to canonical location categories.

    Uses keyword matching against _LOCATION_CANONICAL rather than raw string
    splitting, producing clean deduplicated categories suitable for filtering.
    """
    if not raw or not isinstance(raw, str) or raw.strip() == "":
        return []

    text = re.sub(r"SUBCELLULAR LOCATION:", "", raw, flags=re.IGNORECASE)
    text = re.sub(r"\{[^}]*\}", "", text)   # evidence codes e.g. {ECO:...}
    text = re.sub(r"\[[^\]]*\]", "", text)  # square bracket content
    text = re.sub(r"\([^)]*\)", "", text)   # parenthetical e.g. (By similarity)
    text = re.sub(r"Note=.*", "", text, flags=re.IGNORECASE)
    text = text.lower()

    found: list[str] = []
    for keywords, canonical in _LOCATION_CANONICAL:
        for kw in keywords:
            if kw in text:
                if canonical not in found:
                    found.append(canonical)
                break
    return found


def _parse_functions(func_str: str | None, name_str: str | None = None) -> list[str]:
    """Return matched FUNCTION_CATEGORIES for a protein."""
    parts: list[str] = []
    if pd.notna(func_str) and func_str:
        parts.append(str(func_str))
    if pd.notna(name_str) and name_str:
        parts.append(str(name_str))
    if not parts:
        return []

    combined = " ".join(parts)
    combined = re.sub(r"\{[^}]*\}", "", combined).lower()

    found: list[str] = []
    for category, patterns in FUNCTION_CATEGORIES.items():
        for pat in patterns:
            if re.search(pat, combined, re.IGNORECASE):
                found.append(category)
                break
    return found


def _count_tm_domains(domain_str: str | None) -> int:
    """Count domain spans (e.g. '37..60') in a UniProt Transmembrane/Intramembrane string."""
    if not domain_str or not isinstance(domain_str, str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', domain_str))


def _enrich_df(df: pd.DataFrame) -> pd.DataFrame:
    """Add Location Categories, Function Categories, pLLPS_class, and TMD_count columns."""
    df = df.copy()

    if "Subcellular location [CC]" in df.columns:
        df["Location Categories"] = df["Subcellular location [CC]"].apply(_parse_location)
    else:
        df["Location Categories"] = [[] for _ in range(len(df))]

    has_func = "Function [CC]" in df.columns
    has_name = "Protein names" in df.columns
    df["Function Categories"] = df.apply(
        lambda r: _parse_functions(
            r.get("Function [CC]") if has_func else None,
            r.get("Protein names") if has_name else None,
        ),
        axis=1,
    )

    if "p(LLPS)" in df.columns:
        df["pLLPS_class"] = pd.cut(
            df["p(LLPS)"],
            bins=[-np.inf, 0.4, 0.7, np.inf],
            labels=["Low", "Medium", "High"],
        )

    if "TMD_count" not in df.columns:
        tmd = (
            df["Transmembrane"].apply(_count_tm_domains)
            if "Transmembrane" in df.columns
            else pd.Series(0, index=df.index)
        )
        imd = (
            df["Intramembrane"].apply(_count_tm_domains)
            if "Intramembrane" in df.columns
            else pd.Series(0, index=df.index)
        )
        total = tmd + imd
        if total.sum() > 0:
            df["TMD_count"] = total

    return df


def _all_values(df: pd.DataFrame, col: str) -> list[str]:
    """Return sorted unique scalar values from a list-valued column."""
    vals: set[str] = set()
    for cell in df[col]:
        if isinstance(cell, list):
            vals.update(v for v in cell if v)
    return sorted(vals)


def _load_default() -> pd.DataFrame:
    base = Path(__file__).parent
    full = base / "full_dataset.csv"
    csv_path = full if full.exists() else base / "sample_data.csv"
    return _enrich_df(pd.read_csv(csv_path))


def _df_to_csv_url(df: pd.DataFrame) -> alt.Data:
    """Encode a DataFrame as a base64 CSV data URL (Pyodide-safe, bypasses MaxRowsError)."""
    b64 = base64.b64encode(df.to_csv(index=False).encode()).decode("ascii")
    return alt.Data(url=f"data:text/csv;base64,{b64}", format=alt.DataFormat(type="csv"))


def _chart_to_div(chart: alt.Chart, div_id: str) -> str:
    """Render an Altair chart as a div+script snippet.

    Assumes vega/vega-lite/vega-embed are loaded in the page head.
    The IIFE wrapper prevents vegaEmbed promise races on rapid filter changes.
    """
    spec = json.dumps(chart.to_dict(), separators=(",", ":"))
    return (
        f'<div id="{div_id}" style="width:100%;"></div>\n'
        f"<script>(function(){{"
        f'vegaEmbed("#{div_id}",{spec},{{mode:"vega-lite",renderer:"svg",'
        f"actions:{{export:true,source:false,compiled:false,editor:false}}}})"
        f".catch(console.error);}})();</script>"
    )


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

app_ui = ui.page_fluid(
    ui.head_content(
        ui.tags.title("LLPS Protein Data Explorer"),
        ui.tags.script(src="https://cdn.jsdelivr.net/npm/vega@6"),
        ui.tags.script(src="https://cdn.jsdelivr.net/npm/vega-lite@6.4.1"),
        ui.tags.script(src="https://cdn.jsdelivr.net/npm/vega-embed@7"),
        ui.tags.style("""
            body { font-family: 'Segoe UI', sans-serif; }
            .sidebar-panel { background: #f8f9fa; padding: 1rem; border-radius: 6px; }
            .metric-box { background: #fff; border: 1px solid #dee2e6; border-radius: 6px;
                          padding: 14px 18px; text-align: center; }
            .metric-label { font-size: 12px; color: #6c757d; font-weight: 600;
                            text-transform: uppercase; }
            .metric-val   { font-size: 28px; font-weight: 700; color: #212529; }
        """),
    ),
    ui.div(
        {"class": "container-fluid", "style": "max-width:1400px"},
        ui.h2("LLPS Protein Data Explorer", style="margin-top:1rem"),
        ui.p(
            "Explore protein phase-separation (pLLPS) data. "
            "Use the full dataset or upload your own XLSX/CSV file."
        ),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Data Source"),
                ui.input_radio_buttons(
                    "data_source",
                    None,
                    {"sample": "Full dataset", "upload": "Upload file"},
                    selected="sample",
                ),
                ui.panel_conditional(
                    "input.data_source === 'upload'",
                    ui.input_file(
                        "file_upload",
                        "XLSX or CSV file",
                        accept=[".xlsx", ".csv"],
                        multiple=False,
                    ),
                ),
                ui.output_ui("data_status"),
                ui.hr(),
                ui.h5("Filters"),
                ui.p(
                    {"style": "font-size:12px; color:#6c757d; margin-bottom:8px"},
                    "All plots and the table update together when filters change.",
                ),
                ui.output_ui("filter_controls"),
                width=300,
            ),
            ui.div(
                {"class": "mt-3"},
                # --- Metrics ---
                ui.output_ui("metrics_row"),
                ui.hr(),
                # --- Table ---
                ui.h4("Protein Table"),
                ui.input_text(
                    "search_text",
                    "Search by name / entry ID",
                    placeholder="e.g. PCLO, Q9Y6V0, kinase …",
                ),
                ui.output_data_frame("data_table"),
                ui.hr(),
                # --- Distribution plots ---
                ui.h4("Score & Feature Distributions"),
                ui.row(
                    ui.column(4, ui.output_ui("plot_pllps_dist")),
                    ui.column(4, ui.output_ui("plot_length_dist")),
                    ui.column(4, ui.output_ui("plot_tmd_dist")),
                ),
                ui.hr(),
                # --- pLLPS by Category box plots ---
                ui.h4("p(LLPS) by Category"),
                ui.p(
                    {"style": "color:#6c757d; font-size:13px"},
                    "Box plots (median line, IQR box, 1.5×IQR whiskers) of p(LLPS) "
                    "score across biological categories. Sorted by median p(LLPS) "
                    "descending. Hover over outlier dots for protein details.",
                ),
                ui.row(
                    ui.column(6, ui.output_ui("plot_pllps_by_location")),
                    ui.column(6, ui.output_ui("plot_pllps_by_function")),
                ),
                ui.output_ui("plot_pllps_by_tmd"),
                ui.hr(),
                # --- Scatter plot ---
                ui.h4("Scatter"),
                ui.row(
                    ui.column(
                        3,
                        ui.input_select(
                            "scatter_x",
                            "X axis",
                            choices=["p(LLPS)", "Length", "TMD_count", "n(DPR=> 25)"],
                        ),
                    ),
                    ui.column(
                        3,
                        ui.input_select(
                            "scatter_y",
                            "Y axis",
                            choices=["Length", "p(LLPS)", "TMD_count", "n(DPR=> 25)"],
                        ),
                    ),
                ),
                ui.output_ui("plot_scatter"),
                ui.hr(),
                # --- Category count bars ---
                ui.h4("Category Overview"),
                ui.row(
                    ui.column(6, ui.output_ui("plot_locations")),
                    ui.column(6, ui.output_ui("plot_functions")),
                ),
                ui.hr(),
                # --- Export ---
                ui.h4("Download filtered data"),
                ui.output_ui("export_info"),
                ui.download_button(
                    "download_csv",
                    "Download as CSV",
                    class_="btn btn-primary",
                ),
            ),
        ),
    ),
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------


def server(input: Any, output: Any, session: Any) -> None:
    raw_data: reactive.Value[pd.DataFrame | None] = reactive.Value(None)
    filtered: reactive.Value[pd.DataFrame | None] = reactive.Value(None)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    @reactive.Effect
    @reactive.event(input.data_source, input.file_upload)
    def _load() -> None:
        if input.data_source() == "sample":
            df = _load_default()
        else:
            info = input.file_upload()
            if not info:
                return
            path = info[0]["datapath"]
            df = _enrich_df(
                pd.read_excel(path, engine="openpyxl")
                if path.endswith(".xlsx")
                else pd.read_csv(path)
            )
        raw_data.set(df)
        filtered.set(df)

    @reactive.Effect
    def _initial_load() -> None:
        if raw_data() is None:
            raw_data.set(_load_default())
            filtered.set(_load_default())

    # ------------------------------------------------------------------
    # Filtering — single source of truth for all plots and the table
    # ------------------------------------------------------------------

    @reactive.Effect
    @reactive.event(
        input.filter_pllps,
        input.filter_length,
        input.filter_tmd,
        input.filter_locations,
        input.filter_functions,
        input.filter_class,
        input.search_text,
    )
    def _apply_filters() -> None:
        df = raw_data()
        if df is None:
            return

        if "p(LLPS)" in df.columns:
            pllps_range = input.filter_pllps()
            if pllps_range is not None and len(pllps_range) == 2:
                lo, hi = pllps_range
                df = df[df["p(LLPS)"].between(lo, hi)]

        if "Length" in df.columns:
            length_range = input.filter_length()
            if length_range is not None and len(length_range) == 2:
                lo, hi = length_range
                df = df[df["Length"].between(lo, hi)]

        if "TMD_count" in df.columns:
            try:
                tmd_range = input.filter_tmd()
                if tmd_range is not None and len(tmd_range) == 2:
                    lo, hi = int(tmd_range[0]), int(tmd_range[1])
                    df = df[df["TMD_count"].between(lo, hi)]
            except Exception:
                pass

        if "Location Categories" in df.columns:
            locs = input.filter_locations()
            if locs:
                df = df[df["Location Categories"].apply(lambda x: any(l in x for l in locs))]

        if "Function Categories" in df.columns:
            funcs = input.filter_functions()
            if funcs:
                df = df[df["Function Categories"].apply(lambda x: any(f in x for f in funcs))]

        if "pLLPS_class" in df.columns:
            classes = input.filter_class()
            if classes:
                df = df[df["pLLPS_class"].astype(str).isin(classes)]

        q = input.search_text().strip().lower()
        if q:
            mask = pd.Series(False, index=df.index)
            for col in ("Entry", "Entry name", "Protein names"):
                if col in df.columns:
                    mask |= df[col].astype(str).str.lower().str.contains(q, na=False)
            df = df[mask]

        filtered.set(df)

    # ------------------------------------------------------------------
    # Sidebar outputs
    # ------------------------------------------------------------------

    @output
    @render.ui
    def data_status() -> Any:
        df = raw_data()
        if df is None:
            return ui.div(ui.tags.div("No data loaded", class_="alert alert-warning mt-2"))
        label = "Full dataset" if input.data_source() == "sample" else "Uploaded data"
        return ui.div(
            ui.tags.div(f"{label}: {len(df)} proteins", class_="alert alert-success mt-2")
        )

    @output
    @render.ui
    def filter_controls() -> Any:
        df = raw_data()
        if df is None:
            return ui.div()

        controls: list[Any] = []

        if "p(LLPS)" in df.columns:
            lo, hi = float(df["p(LLPS)"].min()), float(df["p(LLPS)"].max())
            controls.append(
                ui.input_slider("filter_pllps", "p(LLPS) range", min=lo, max=hi,
                                value=[lo, hi], step=0.01)
            )

        if "Length" in df.columns:
            lo, hi = int(df["Length"].min()), int(df["Length"].max())
            controls.append(
                ui.input_slider("filter_length", "Protein length", min=lo, max=hi,
                                value=[lo, hi], step=10)
            )

        if "TMD_count" in df.columns:
            lo, hi = int(df["TMD_count"].min()), int(df["TMD_count"].max())
            controls.append(
                ui.input_slider("filter_tmd", "Transmembrane domains", min=lo, max=hi,
                                value=[lo, hi], step=1)
            )

        if "pLLPS_class" in df.columns:
            controls.append(
                ui.input_checkbox_group(
                    "filter_class",
                    "pLLPS class",
                    choices=["High", "Medium", "Low"],
                    selected=["High", "Medium", "Low"],
                )
            )

        if "Location Categories" in df.columns:
            locs = _all_values(df, "Location Categories")
            if locs:
                controls.append(
                    ui.input_selectize(
                        "filter_locations", "Subcellular location",
                        choices=locs, multiple=True
                    )
                )

        if "Function Categories" in df.columns:
            funcs = _all_values(df, "Function Categories")
            if funcs:
                controls.append(
                    ui.input_selectize(
                        "filter_functions", "Functional category",
                        choices=funcs, multiple=True
                    )
                )

        return ui.TagList(*controls)

    # ------------------------------------------------------------------
    # Metrics row
    # ------------------------------------------------------------------

    @output
    @render.ui
    def metrics_row() -> Any:
        df = filtered()
        if df is None:
            return ui.div()

        total = len(df)
        high = int((df["pLLPS_class"].astype(str) == "High").sum()) if "pLLPS_class" in df.columns else "N/A"
        avg_p = f"{df['p(LLPS)'].mean():.3f}" if "p(LLPS)" in df.columns else "N/A"
        avg_l = f"{df['Length'].mean():.0f}" if "Length" in df.columns else "N/A"

        def box(label: str, value: Any) -> Any:
            return ui.column(
                3,
                ui.div(
                    ui.div(label, class_="metric-label"),
                    ui.div(str(value), class_="metric-val"),
                    class_="metric-box",
                ),
            )

        return ui.row(
            box("Proteins shown", total),
            box("High pLLPS", high),
            box("Mean p(LLPS)", avg_p),
            box("Mean length (aa)", avg_l),
        )

    # ------------------------------------------------------------------
    # Protein table — sidebar filters are the single filter source
    # ------------------------------------------------------------------

    @output
    @render.data_frame
    def data_table() -> render.DataGrid:
        df = filtered()
        req(df is not None)
        display_cols = [
            c for c in ["Entry", "Entry name", "Protein names", "p(LLPS)",
                        "pLLPS_class", "Length", "TMD_count", "Organism"]
            if c in df.columns
        ]
        return render.DataGrid(df[display_cols], height="420px")

    # ------------------------------------------------------------------
    # Distribution charts
    # ------------------------------------------------------------------

    @output
    @render.ui
    def plot_pllps_dist() -> Any:
        df = filtered()
        if df is None or "p(LLPS)" not in df.columns:
            return ui.div("p(LLPS) column not available")

        bins = pd.cut(df["p(LLPS)"], bins=30)
        agg = (
            df.assign(_bin=bins)
            .groupby(["_bin", "pLLPS_class"], observed=True)
            .size()
            .reset_index(name="Count")
        )
        agg["bin_mid"] = agg["_bin"].apply(lambda x: round(x.mid, 4))
        agg["pLLPS_class"] = agg["pLLPS_class"].astype(str)

        chart = (
            alt.Chart(agg[["bin_mid", "pLLPS_class", "Count"]])
            .mark_bar(stroke="white", strokeWidth=0.5)
            .encode(
                x=alt.X("bin_mid:Q", title="p(LLPS) score", axis=alt.Axis(format=".2f")),
                y=alt.Y("Count:Q"),
                color=alt.Color("pLLPS_class:N", scale=_PLLPS_COLOR_SCALE, title="pLLPS class"),
                tooltip=[
                    alt.Tooltip("bin_mid:Q", title="p(LLPS)", format=".3f"),
                    alt.Tooltip("pLLPS_class:N", title="Class"),
                    alt.Tooltip("Count:Q"),
                ],
            )
            .properties(title="Distribution of p(LLPS) scores", width="container", height=280)
        )
        return ui.HTML(_chart_to_div(chart, "pllps_hist"))

    @output
    @render.ui
    def plot_length_dist() -> Any:
        df = filtered()
        if df is None or "Length" not in df.columns:
            return ui.div("Length column not available")

        bins = pd.cut(df["Length"], bins=30)
        agg = (
            df.assign(_bin=bins)
            .groupby(["_bin", "pLLPS_class"], observed=True)
            .size()
            .reset_index(name="Count")
        )
        agg["bin_mid"] = agg["_bin"].apply(lambda x: round(x.mid, 0))
        agg["pLLPS_class"] = agg["pLLPS_class"].astype(str)

        chart = (
            alt.Chart(agg[["bin_mid", "pLLPS_class", "Count"]])
            .mark_bar(stroke="white", strokeWidth=0.5)
            .encode(
                x=alt.X("bin_mid:Q", title="Protein length (aa)"),
                y=alt.Y("Count:Q"),
                color=alt.Color("pLLPS_class:N", scale=_PLLPS_COLOR_SCALE, title="pLLPS class"),
                tooltip=[
                    alt.Tooltip("bin_mid:Q", title="Length (aa)", format=".0f"),
                    alt.Tooltip("pLLPS_class:N", title="Class"),
                    alt.Tooltip("Count:Q"),
                ],
            )
            .properties(title="Distribution of protein lengths", width="container", height=280)
        )
        return ui.HTML(_chart_to_div(chart, "length_hist"))

    @output
    @render.ui
    def plot_tmd_dist() -> Any:
        df = filtered()
        if df is None or "TMD_count" not in df.columns:
            return ui.div("TMD_count column not available")

        agg = (
            df.groupby(["TMD_count", "pLLPS_class"], observed=True)
            .size()
            .reset_index(name="Count")
        )
        agg["pLLPS_class"] = agg["pLLPS_class"].astype(str)

        chart = (
            alt.Chart(agg)
            .mark_bar()
            .encode(
                x=alt.X("TMD_count:O", title="Number of TM domains"),
                y=alt.Y("Count:Q", title="Protein count"),
                color=alt.Color("pLLPS_class:N", scale=_PLLPS_COLOR_SCALE, title="pLLPS class"),
                tooltip=[
                    alt.Tooltip("TMD_count:O", title="TM domains"),
                    alt.Tooltip("pLLPS_class:N", title="Class"),
                    alt.Tooltip("Count:Q"),
                ],
            )
            .properties(title="Transmembrane domain count", width="container", height=280)
        )
        return ui.HTML(_chart_to_div(chart, "tmd_hist"))

    # ------------------------------------------------------------------
    # p(LLPS) by Category — box plots
    # ------------------------------------------------------------------

    @output
    @render.ui
    def plot_pllps_by_location() -> Any:
        df = filtered()
        if df is None or "Location Categories" not in df.columns or "p(LLPS)" not in df.columns:
            return ui.div("Location/pLLPS data not available")

        exploded = df.explode("Location Categories").dropna(subset=["Location Categories"])
        exploded = exploded[exploded["Location Categories"] != ""]
        if len(exploded) < 3:
            return ui.div("Insufficient location data in filtered set")

        top_locs = exploded["Location Categories"].value_counts().head(15).index.tolist()
        plot_df = (
            exploded[exploded["Location Categories"].isin(top_locs)]
            [["Location Categories", "p(LLPS)"]]
            .copy()
        )

        # Sort categories by median pLLPS descending
        order = (
            plot_df.groupby("Location Categories")["p(LLPS)"]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )

        chart = (
            alt.Chart(_df_to_csv_url(plot_df))
            .mark_boxplot(extent=1.5)
            .encode(
                x=alt.X(
                    "Location Categories:N",
                    sort=order,
                    axis=alt.Axis(labelAngle=-40, labelLimit=120),
                    title=None,
                ),
                y=alt.Y(
                    "p(LLPS):Q",
                    scale=alt.Scale(domain=[0, 1]),
                    title="p(LLPS)",
                ),
                color=alt.Color("Location Categories:N", legend=None),
            )
            .properties(
                title="p(LLPS) by subcellular location (top 15)",
                width="container",
                height=320,
            )
        )
        return ui.HTML(_chart_to_div(chart, "pllps_by_loc"))

    @output
    @render.ui
    def plot_pllps_by_function() -> Any:
        df = filtered()
        if df is None or "Function Categories" not in df.columns or "p(LLPS)" not in df.columns:
            return ui.div("Function/pLLPS data not available")

        exploded = df.explode("Function Categories").dropna(subset=["Function Categories"])
        exploded = exploded[exploded["Function Categories"] != ""]
        if len(exploded) < 3:
            return ui.div("Insufficient function data in filtered set")

        plot_df = exploded[["Function Categories", "p(LLPS)"]].copy()

        order = (
            plot_df.groupby("Function Categories")["p(LLPS)"]
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )

        chart = (
            alt.Chart(_df_to_csv_url(plot_df))
            .mark_boxplot(extent=1.5)
            .encode(
                x=alt.X(
                    "Function Categories:N",
                    sort=order,
                    axis=alt.Axis(labelAngle=-40, labelLimit=140),
                    title=None,
                ),
                y=alt.Y(
                    "p(LLPS):Q",
                    scale=alt.Scale(domain=[0, 1]),
                    title="p(LLPS)",
                ),
                color=alt.Color("Function Categories:N", legend=None),
            )
            .properties(
                title="p(LLPS) by functional category",
                width="container",
                height=320,
            )
        )
        return ui.HTML(_chart_to_div(chart, "pllps_by_func"))

    @output
    @render.ui
    def plot_pllps_by_tmd() -> Any:
        df = filtered()
        if df is None or "TMD_count" not in df.columns or "p(LLPS)" not in df.columns:
            return ui.div("TMD/pLLPS data not available")

        # Cap at 20 TM domains for readability; show 0 separately as non-membrane
        plot_df = df[df["TMD_count"] <= 20][["TMD_count", "p(LLPS)"]].copy()
        if len(plot_df) < 3:
            return ui.div("Insufficient data")

        n_proteins = len(plot_df)
        chart = (
            alt.Chart(_df_to_csv_url(plot_df))
            .mark_boxplot(extent=1.5)
            .encode(
                x=alt.X("TMD_count:O", title="Number of TM domains"),
                y=alt.Y(
                    "p(LLPS):Q",
                    scale=alt.Scale(domain=[0, 1]),
                    title="p(LLPS)",
                ),
                color=alt.Color(
                    "TMD_count:O",
                    legend=None,
                    scale=alt.Scale(scheme="viridis"),
                ),
            )
            .properties(
                title=f"p(LLPS) by TM domain count (n={n_proteins:,})",
                width="container",
                height=300,
            )
        )
        return ui.HTML(_chart_to_div(chart, "pllps_by_tmd"))

    # ------------------------------------------------------------------
    # Scatter chart
    # ------------------------------------------------------------------

    @output
    @render.ui
    def plot_scatter() -> Any:
        df = filtered()
        req(df is not None)
        x_col = input.scatter_x()
        y_col = input.scatter_y()
        if x_col not in df.columns or y_col not in df.columns:
            return ui.div(f"Column(s) not available: {x_col}, {y_col}")

        cols = [c for c in [x_col, y_col, "pLLPS_class", "Entry", "Protein names"] if c in df.columns]
        scatter_df = df[cols].copy()
        if "pLLPS_class" in scatter_df.columns:
            scatter_df["pLLPS_class"] = scatter_df["pLLPS_class"].astype(str)

        tooltip = [alt.Tooltip(f"{x_col}:Q"), alt.Tooltip(f"{y_col}:Q")]
        if "Entry" in scatter_df.columns:
            tooltip.append(alt.Tooltip("Entry:N"))
        if "Protein names" in scatter_df.columns:
            tooltip.append(alt.Tooltip("Protein names:N"))

        chart = (
            alt.Chart(_df_to_csv_url(scatter_df))
            .mark_circle(opacity=0.5, size=15)
            .encode(
                x=alt.X(f"{x_col}:Q"),
                y=alt.Y(f"{y_col}:Q"),
                color=alt.Color("pLLPS_class:N", scale=_PLLPS_COLOR_SCALE, title="pLLPS class"),
                tooltip=tooltip,
            )
            .properties(title=f"{y_col} vs {x_col}", width="container", height=350)
        )
        return ui.HTML(_chart_to_div(chart, "scatter_plot"))

    # ------------------------------------------------------------------
    # Category count bar charts (overview)
    # ------------------------------------------------------------------

    @output
    @render.ui
    def plot_locations() -> Any:
        df = filtered()
        if df is None or "Location Categories" not in df.columns:
            return ui.div("Location data not available")
        exploded = df.explode("Location Categories").dropna(subset=["Location Categories"])
        exploded = exploded[exploded["Location Categories"] != ""]
        if len(exploded) == 0:
            return ui.div("No location data in filtered set")

        counts = exploded["Location Categories"].value_counts().reset_index()
        counts.columns = ["Location", "Count"]

        chart = (
            alt.Chart(counts)
            .mark_bar()
            .encode(
                x=alt.X("Location:N", sort="-y", axis=alt.Axis(labelAngle=-40, labelLimit=120), title=None),
                y=alt.Y("Count:Q"),
                color=alt.Color("Count:Q", scale=alt.Scale(scheme="blues"), legend=None),
                tooltip=["Location:N", "Count:Q"],
            )
            .properties(title="Proteins per subcellular location", width="container", height=300)
        )
        return ui.HTML(_chart_to_div(chart, "loc_plot"))

    @output
    @render.ui
    def plot_functions() -> Any:
        df = filtered()
        if df is None or "Function Categories" not in df.columns:
            return ui.div("Function data not available")
        exploded = df.explode("Function Categories").dropna(subset=["Function Categories"])
        exploded = exploded[exploded["Function Categories"] != ""]
        if len(exploded) == 0:
            return ui.div("No function data in filtered set")

        counts = exploded["Function Categories"].value_counts().reset_index()
        counts.columns = ["Function", "Count"]

        chart = (
            alt.Chart(counts)
            .mark_bar()
            .encode(
                x=alt.X("Function:N", sort="-y", axis=alt.Axis(labelAngle=-40, labelLimit=140), title=None),
                y=alt.Y("Count:Q"),
                color=alt.Color("Count:Q", scale=alt.Scale(scheme="greens"), legend=None),
                tooltip=["Function:N", "Count:Q"],
            )
            .properties(title="Proteins per functional category", width="container", height=300)
        )
        return ui.HTML(_chart_to_div(chart, "func_plot"))

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    @output
    @render.ui
    def export_info() -> Any:
        df = filtered()
        if df is None:
            return ui.p("No data available")
        return ui.p(f"{len(df)} proteins · {len(df.columns)} columns (after filters)")

    @render.download(filename="llps_proteins_filtered.csv")
    def download_csv() -> str:
        df = filtered()
        if df is not None:
            export_df = df.drop(
                columns=[c for c in ("Location Categories", "Function Categories") if c in df.columns]
            )
            return export_df.to_csv(index=False)
        return ""


app = App(app_ui, server)
