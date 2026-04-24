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

FUNCTION_CATEGORIES: dict[str, list[str]] = {
    "Hydrolase": [r"\bhydrolase\b", r"\besterase\b", r"\blipase\b", r"\bnuclease\b"],
    "Kinase": [r"\bkinase\b", r"phosphotransferase", r"protein\s+kinase"],
    "Ligase": [r"\bligase\b", r"\bsynthetase\b", r"ubiquitin[- ]protein\s+ligase"],
    "Oxidoreductase": [r"oxidoreductase", r"\bdehydrogenase\b", r"\boxidase\b"],
    "Transferase": [r"\btransferase\b", r"methyltransferase", r"acetyltransferase"],
    "Protease": [r"\bprotease\b", r"\bpeptidase\b", r"proteolytic"],
    "Phosphatase": [r"\bphosphatase\b", r"protein\s+phosphatase"],
    "DNA-binding": [r"dna[- ]binding", r"binds\s+(to\s+)?dna"],
    "RNA-binding": [r"rna[- ]binding", r"binds\s+(to\s+)?rna"],
    "Receptor": [r"\breceptor\b", r"receptor\s+activity"],
    "Ion channel": [r"ion\s*channel", r"cation\s*channel"],
    "Transporter": [r"\btransporter\b", r"transport\s+protein"],
    "Chaperone": [r"\bchaperone\b", r"heat\s+shock\s+protein"],
    "Transcription regulation": [r"transcription\s*(factor|regulator|regulation)"],
}

_PLLPS_COLOR_SCALE = alt.Scale(
    domain=["High", "Medium", "Low"],
    range=["#e74c3c", "#f39c12", "#3498db"],
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _parse_location(raw: str | None) -> list[str]:
    """Extract normalised location tokens from a Subcellular location [CC] string."""
    if not raw or not isinstance(raw, str) or raw.strip() == "":
        return []

    text = re.sub(r"SUBCELLULAR LOCATION:", "", raw, flags=re.IGNORECASE)
    text = re.sub(r"\{[^}]*\}", "", text)  # remove evidence codes

    tokens: list[str] = []
    for part in re.split(r"[.;,]", text):
        part = part.strip()
        if part and len(part) > 1:
            part = re.sub(r"\s+", " ", part)
            tokens.append(part.title())
    return list(dict.fromkeys(tokens))  # deduplicate while preserving order


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
            if re.search(pat, combined):
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

    # Compute TMD_count from raw annotation columns if not already present
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
                ui.output_ui("filter_controls"),
                width=280,
            ),
            ui.div(
                {"class": "mt-3"},
                # --- Metrics ---
                ui.output_ui("metrics_row"),
                ui.hr(),
                # --- Table ---
                ui.input_text(
                    "search_text",
                    "Search",
                    placeholder="Protein name, entry ID, keyword …",
                ),
                ui.output_data_frame("data_table"),
                ui.hr(),
                # --- Distribution plots ---
                ui.h4("Distributions"),
                ui.row(
                    ui.column(4, ui.output_ui("plot_pllps_dist")),
                    ui.column(4, ui.output_ui("plot_length_dist")),
                    ui.column(4, ui.output_ui("plot_tmd_dist")),
                ),
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
                # --- Location & Function plots ---
                ui.h4("Locations and Functions"),
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

    # Bootstrap sample data on startup
    @reactive.Effect
    def _initial_load() -> None:
        if raw_data() is None:
            raw_data.set(_load_default())
            filtered.set(_load_default())

    # ------------------------------------------------------------------
    # Filtering
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
        return render.DataGrid(df[display_cols], filters=True, height="420px")

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
    # Location & function charts
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

        counts = exploded["Location Categories"].value_counts().head(20).reset_index()
        counts.columns = ["Location", "Count"]

        chart = (
            alt.Chart(counts)
            .mark_bar()
            .encode(
                x=alt.X("Location:N", sort="-y", axis=alt.Axis(labelAngle=-40), title=None),
                y=alt.Y("Count:Q"),
                color=alt.Color("Count:Q", scale=alt.Scale(scheme="blues"), legend=None),
                tooltip=["Location:N", "Count:Q"],
            )
            .properties(title="Top 20 subcellular locations", width="container", height=300)
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
                x=alt.X("Function:N", sort="-y", axis=alt.Axis(labelAngle=-30), title=None),
                y=alt.Y("Count:Q"),
                color=alt.Color("Count:Q", scale=alt.Scale(scheme="greens"), legend=None),
                tooltip=["Function:N", "Count:Q"],
            )
            .properties(title="Functional categories", width="container", height=300)
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
            # Drop list-valued columns before CSV export
            export_df = df.drop(
                columns=[c for c in ("Location Categories", "Function Categories") if c in df.columns]
            )
            return export_df.to_csv(index=False)
        return ""


app = App(app_ui, server)
