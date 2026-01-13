"""
LLPS Protein Data Explorer Dashboard - Shiny for Python

An interactive Shiny dashboard for exploring Liquid-Liquid Phase Separation
(LLPS) protein data. Supports XLSX file upload and provides interactive 
visualizations for data exploration.

Usage:
    shiny run shiny_app.py --reload --port 8000
"""

from shiny import App, reactive, render, ui, req
from shiny.types import FileInfo
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import re
import numpy as np
from scipy import stats
from htmltools import HTML, css
import io
import sys

# Add parent directory to path to import llps_functions
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import STRING interaction functions from separate module
from llps_functions import (
    fetch_string_interactions,
    match_interactions_to_pllps,
    analyze_interaction_enrichment
)

# ============================================================================
# HELPER FUNCTIONS (from app.py)
# ============================================================================

def parse_location(location_str):
    """Parse a subcellular location string from UniProt"""
    if pd.isna(location_str) or location_str == '':
        return []
    
    location_str = str(location_str)
    location_str = re.sub(r'\{[^}]*\}', '', location_str)
    location_str = re.sub(r'\[[^\]]*\]:?', '', location_str)
    location_str = re.sub(r'Isoform\s+[^:]+:\s*', '', location_str, flags=re.IGNORECASE)
    location_str = re.sub(r'\([^)]*\)', '', location_str)
    location_str = re.sub(r'^SUBCELLULAR LOCATION:\s*', '', location_str, flags=re.IGNORECASE)
    location_str = re.sub(r'\s*Note=.*', '', location_str, flags=re.IGNORECASE)
    
    parts = re.split(r'[;,.]', location_str)
    locations = []
    for part in parts:
        cleaned = part.strip()
        if len(cleaned) >= 2 and cleaned not in locations:
            locations.append(cleaned)
    
    return locations


def add_location_columns(df):
    """Add parsed location columns to the dataframe"""
    if 'Subcellular location [CC]' not in df.columns:
        return df
    df = df.copy()
    df['Location Categories'] = df['Subcellular location [CC]'].apply(parse_location)
    return df


def get_all_locations(df):
    """Get all unique location categories from the dataframe"""
    if 'Location Categories' not in df.columns:
        return []
    all_locations = set()
    for locations in df['Location Categories']:
        if isinstance(locations, list):
            all_locations.update(locations)
    return sorted(all_locations)


# Predefined function categories
FUNCTION_CATEGORIES = {
    'Hydrolase': [r'\bhydrolase\b', r'\besterase\b', r'\blipase\b', r'\bnuclease\b'],
    'Kinase': [r'\bkinase\b', r'phosphotransferase', r'protein\s+kinase'],
    'Ligase': [r'\bligase\b', r'\bsynthetase\b', r'ubiquitin[- ]protein\s+ligase'],
    'Oxidoreductase': [r'oxidoreductase', r'\bdehydrogenase\b', r'\boxidase\b'],
    'Transferase': [r'\btransferase\b', r'methyltransferase', r'acetyltransferase'],
    'Protease': [r'\bprotease\b', r'\bpeptidase\b', r'proteolytic'],
    'Phosphatase': [r'\bphosphatase\b', r'protein\s+phosphatase'],
    'DNA-binding': [r'dna[- ]binding', r'binds\s+(to\s+)?dna'],
    'RNA-binding': [r'rna[- ]binding', r'binds\s+(to\s+)?rna'],
    'Receptor': [r'\breceptor\b', r'receptor\s+activity'],
    'Ion channel': [r'ion\s*channel', r'cation\s*channel'],
    'Transporter': [r'\btransporter\b', r'transport\s+protein'],
    'Chaperone': [r'\bchaperone\b', r'heat\s+shock\s+protein'],
    'Transcription regulation': [r'transcription\s*(factor|regulator|regulation)'],
}


def parse_function(function_str, protein_name_str=None):
    """Parse function annotations"""
    categories = []
    text_parts = []
    
    if pd.notna(function_str) and function_str != '':
        text_parts.append(str(function_str))
    if pd.notna(protein_name_str) and protein_name_str != '':
        text_parts.append(str(protein_name_str))
    
    if not text_parts:
        return categories
    
    combined_text = ' '.join(text_parts)
    combined_text = re.sub(r'\{[^}]*\}', '', combined_text)
    combined_text_lower = combined_text.lower()
    
    for category, patterns in FUNCTION_CATEGORIES.items():
        for pattern in patterns:
            if re.search(pattern, combined_text_lower):
                if category not in categories:
                    categories.append(category)
                break
    
    return categories


def add_function_columns(df):
    """Add parsed function category columns to the dataframe"""
    df = df.copy()
    has_function = 'Function [CC]' in df.columns
    has_protein_name = 'Protein names' in df.columns
    
    if not has_function and not has_protein_name:
        return df
    
    def extract_functions(row):
        func_str = row.get('Function [CC]') if has_function else None
        name_str = row.get('Protein names') if has_protein_name else None
        return parse_function(func_str, name_str)
    
    df['Function Categories'] = df.apply(extract_functions, axis=1)
    return df


def get_all_functions(df):
    """Get all unique function categories from the dataframe"""
    if 'Function Categories' not in df.columns:
        return []
    all_functions = set()
    for functions in df['Function Categories']:
        if isinstance(functions, list):
            all_functions.update(functions)
    return sorted(all_functions)


# ============================================================================
# SHINY UI
# ============================================================================

app_ui = ui.page_fluid(
    ui.head_content(
        ui.tags.title("LLPS Protein Data Explorer"),
        ui.tags.style("""
            .navbar { background-color: #2c3e50 !important; }
            .value-box { background: #ecf0f1; padding: 20px; border-radius: 5px; margin: 10px 0; }
            .metric-title { font-size: 14px; color: #7f8c8d; font-weight: 600; }
            .metric-value { font-size: 32px; color: #2c3e50; font-weight: bold; }
        """)
    ),
    
    ui.div(
        {"class": "container-fluid"},
        ui.h1("🧬 LLPS Protein Data Explorer", style="margin-top: 20px;"),
        ui.p("This dashboard allows you to explore protein Liquid-Liquid Phase Separation (LLPS) data. Upload your XLSX file or use the sample data to get started."),
        
        ui.layout_sidebar(
            ui.sidebar(
                ui.card(
                    ui.h4("📁 Data Source"),
                    ui.input_radio_buttons(
                        "data_source",
                        "Choose data source:",
                        {"full": "Use full dataset", "sample": "Use sample data", "upload": "Upload XLSX file"},
                        selected="full"
                    ),
                    ui.panel_conditional(
                        "input.data_source === 'upload'",
                        ui.input_file("file_upload", "Upload XLSX file", accept=[".xlsx"], multiple=False)
                    ),
                    ui.output_ui("data_status"),
                ),
                ui.card(
                    ui.h4("🔧 Filters"),
                    ui.output_ui("filter_controls"),
                ),
                width=300,
            ),
            
            ui.navset_tab(
                ui.nav_panel(
                    "📊 Data Explorer",
                    ui.div(
                        {"class": "mt-3"},
                        ui.h3("📊 Data Overview"),
                        ui.output_ui("data_overview"),
                        ui.hr(),
                        ui.h3("📋 Data Table"),
                        ui.input_text("search_text", "🔍 Search proteins", placeholder="Enter protein name, entry, or keyword..."),
                        ui.output_data_frame("data_table"),
                        ui.hr(),
                        ui.h3("📈 Visualizations"),
                        ui.navset_tab(
                            ui.nav_panel("Distribution", ui.output_ui("plot_distribution")),
                            ui.nav_panel("Scatter Plot", ui.output_ui("plot_scatter")),
                            ui.nav_panel("Location Analysis", ui.output_ui("plot_locations")),
                            ui.nav_panel("Function Analysis", ui.output_ui("plot_functions")),
                            ui.nav_panel("Length Analysis", ui.output_ui("plot_length")),
                        ),
                    )
                ),
                ui.nav_panel(
                    "🔗 Protein Interactions",
                    ui.div(
                        {"class": "mt-3"},
                        ui.h3("🔗 Protein Interaction Analysis"),
                        ui.markdown("""
                        Analyze protein-protein interactions from the STRING database to understand:
                        - Whether high pLLPS proteins preferentially interact with each other
                        - Network properties of high pLLPS proteins
                        - Interaction patterns across pLLPS classes
                        """),
                        ui.h4("⚙️ Analysis Settings"),
                        ui.row(
                            ui.column(4, ui.input_slider("pllps_threshold", "High pLLPS Threshold", min=0.5, max=0.95, value=0.7, step=0.05)),
                            ui.column(4, ui.input_slider("string_score", "STRING Confidence Score", min=400, max=900, value=700, step=100)),
                            ui.column(4, ui.input_numeric("max_proteins", "Max Proteins to Analyze", value=100, min=10, max=500)),
                        ),
                        ui.output_ui("interaction_status"),
                        ui.input_action_button("fetch_interactions", "🚀 Fetch Interactions from STRING", class_="btn btn-primary btn-lg mt-3"),
                        ui.output_ui("interaction_results"),
                    )
                ),
                ui.nav_panel(
                    "💾 Export",
                    ui.div(
                        {"class": "mt-3"},
                        ui.h3("💾 Export Data"),
                        ui.output_ui("export_controls"),
                    )
                ),
            ),
        ),
    ),
)


# ============================================================================
# SHINY SERVER
# ============================================================================

def server(input, output, session):
    # Reactive values
    data = reactive.Value(None)
    filtered_data = reactive.Value(None)
    interaction_data = reactive.Value(None)
    
    @reactive.Effect
    @reactive.event(input.data_source, input.file_upload)
    def load_data():
        """Load data from file upload, full dataset, or sample data"""
        df = None
        
        if input.data_source() == "full":
            # Load full dataset
            full_path = Path(__file__).parent / "Human Phase separation data.xlsx"
            if full_path.exists():
                df = pd.read_excel(full_path, engine='openpyxl')
        elif input.data_source() == "sample":
            sample_path = Path(__file__).parent / "data" / "sample_data.xlsx"
            if sample_path.exists():
                df = pd.read_excel(sample_path, engine='openpyxl')
        elif input.data_source() == "upload":
            file_info: list[FileInfo] | None = input.file_upload()
            if file_info:
                df = pd.read_excel(file_info[0]["datapath"], engine='openpyxl')
        
        if df is not None:
            df = add_location_columns(df)
            df = add_function_columns(df)
            data.set(df)
            filtered_data.set(df)
    
    @output
    @render.ui
    def data_status():
        df = data()
        if df is not None:
            dataset_name = "Full dataset" if input.data_source() == "full" else ("Sample data" if input.data_source() == "sample" else "Uploaded data")
            return ui.div(ui.tags.div(f"✅ {dataset_name} loaded ({len(df)} proteins)", class_="alert alert-success"))
        return ui.div(ui.tags.div("⏳ No data loaded", class_="alert alert-warning"))
    
    @output
    @render.ui
    def filter_controls():
        df = data()
        if df is None:
            return ui.div()
        
        controls = []
        
        if 'p(LLPS)' in df.columns:
            pllps_min = float(df['p(LLPS)'].min())
            pllps_max = float(df['p(LLPS)'].max())
            controls.append(
                ui.input_slider("filter_pllps", "p(LLPS) Range", min=pllps_min, max=pllps_max, value=[pllps_min, pllps_max], step=0.01)
            )
        
        if 'Length' in df.columns:
            len_min = int(df['Length'].min())
            len_max = int(df['Length'].max())
            controls.append(
                ui.input_slider("filter_length", "Protein Length Range", min=len_min, max=len_max, value=[len_min, len_max], step=10)
            )
        
        if 'Location Categories' in df.columns:
            locations = get_all_locations(df)
            if locations:
                controls.append(
                    ui.input_selectize("filter_locations", "Subcellular Location", choices=locations, multiple=True)
                )
        
        if 'Function Categories' in df.columns:
            functions = get_all_functions(df)
            if functions:
                controls.append(
                    ui.input_selectize("filter_functions", "Functional Category", choices=functions, multiple=True)
                )
        
        return ui.TagList(*controls)
    
    @reactive.Effect
    @reactive.event(input.filter_pllps, input.filter_length, input.filter_locations, input.filter_functions)
    def apply_filters():
        df = data()
        if df is None:
            return
        
        filtered = df.copy()
        
        if 'p(LLPS)' in df.columns and input.filter_pllps():
            pllps_range = input.filter_pllps()
            filtered = filtered[(filtered['p(LLPS)'] >= pllps_range[0]) & (filtered['p(LLPS)'] <= pllps_range[1])]
        
        if 'Length' in df.columns and input.filter_length():
            len_range = input.filter_length()
            filtered = filtered[(filtered['Length'] >= len_range[0]) & (filtered['Length'] <= len_range[1])]
        
        if input.filter_locations():
            selected_locs = input.filter_locations()
            def location_matches(row):
                if isinstance(row.get('Location Categories'), list):
                    return any(loc in selected_locs for loc in row['Location Categories'])
                return False
            filtered = filtered[filtered.apply(location_matches, axis=1)]
        
        if input.filter_functions():
            selected_funcs = input.filter_functions()
            def function_matches(row):
                if isinstance(row.get('Function Categories'), list):
                    return any(func in selected_funcs for func in row['Function Categories'])
                return False
            filtered = filtered[filtered.apply(function_matches, axis=1)]
        
        filtered_data.set(filtered)
    
    @output
    @render.ui
    def data_overview():
        df = filtered_data()
        if df is None or len(df) == 0:
            return ui.div("No data available")
        
        metrics = []
        metrics.append(ui.column(4, ui.div(ui.div("Total Proteins", class_="metric-title"), ui.div(str(len(df)), class_="metric-value"), class_="value-box")))
        
        if 'Length' in df.columns:
            metrics.append(ui.column(4, ui.div(ui.div("Avg. Protein Length", class_="metric-title"), ui.div(f"{df['Length'].mean():.0f}", class_="metric-value"), class_="value-box")))
        
        if 'p(LLPS)' in df.columns:
            metrics.append(ui.column(4, ui.div(ui.div("Avg. p(LLPS)", class_="metric-title"), ui.div(f"{df['p(LLPS)'].mean():.3f}", class_="metric-value"), class_="value-box")))
        
        return ui.row(*metrics)
    
    @output
    @render.data_frame
    def data_table():
        df = filtered_data()
        if df is None:
            return pd.DataFrame()
        
        if input.search_text():
            search_term = input.search_text().lower()
            search_cols = ['Entry', 'Entry name', 'Protein names']
            mask = pd.Series([False] * len(df))
            for col in search_cols:
                if col in df.columns:
                    mask = mask | df[col].astype(str).str.lower().str.contains(search_term, na=False)
            df = df[mask]
        
        return render.DataGrid(df, height="400px", width="100%", filters=True)
    
    @output
    @render.ui
    def plot_distribution():
        df = filtered_data()
        if df is None or 'p(LLPS)' not in df.columns:
            return ui.div("p(LLPS) column not found")
        
        fig = px.histogram(df, x='p(LLPS)', nbins=50, title="Distribution of LLPS Probability", labels={'p(LLPS)': 'Probability of LLPS'})
        return ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="dist_plot"))
    
    @output
    @render.ui
    def plot_scatter():
        df = filtered_data()
        if df is None:
            return ui.div("No data available")
        
        numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
        if len(numeric_cols) < 2:
            return ui.div("Need at least 2 numeric columns")
        
        x_col = 'Length' if 'Length' in numeric_cols else numeric_cols[0]
        y_col = 'p(LLPS)' if 'p(LLPS)' in numeric_cols else numeric_cols[1]
        
        fig = px.scatter(df, x=x_col, y=y_col, title=f"{y_col} vs {x_col}",
                        hover_data=['Entry', 'Entry name'] if all(c in df.columns for c in ['Entry', 'Entry name']) else None)
        return ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="scatter_plot"))
    
    @output
    @render.ui
    def plot_locations():
        df = filtered_data()
        if df is None or 'Location Categories' not in df.columns:
            return ui.div("Location data not available")
        
        df_locations = df.explode('Location Categories')
        df_locations = df_locations[df_locations['Location Categories'].notna()]
        
        if len(df_locations) == 0:
            return ui.div("No location data")
        
        location_counts = df_locations['Location Categories'].value_counts().head(20)
        fig = px.bar(x=location_counts.index, y=location_counts.values, title="Top 20 Protein Locations", labels={'x': 'Location', 'y': 'Count'})
        fig.update_layout(xaxis_tickangle=-45)
        return ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="loc_plot"))
    
    @output
    @render.ui
    def plot_functions():
        df = filtered_data()
        if df is None or 'Function Categories' not in df.columns:
            return ui.div("Function data not available")
        
        df_functions = df.explode('Function Categories')
        df_functions = df_functions[df_functions['Function Categories'].notna()]
        df_functions = df_functions[df_functions['Function Categories'] != '']
        
        if len(df_functions) == 0:
            return ui.div("No function data")
        
        function_counts = df_functions['Function Categories'].value_counts().head(20)
        fig = px.bar(x=function_counts.index, y=function_counts.values, title="Top 20 Functional Categories", labels={'x': 'Function', 'y': 'Count'})
        fig.update_layout(xaxis_tickangle=-45)
        return ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="func_plot"))
    
    @output
    @render.ui
    def plot_length():
        df = filtered_data()
        if df is None or 'Length' not in df.columns:
            return ui.div("Length data not available")
        
        fig = px.histogram(df, x='Length', nbins=50, title="Distribution of Protein Lengths", labels={'Length': 'Protein Length (amino acids)'})
        return ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="length_plot"))
    
    @output
    @render.ui
    def interaction_status():
        df = data()
        if df is None:
            return ui.div(ui.tags.div("Load data first", class_="alert alert-warning"))
        
        threshold = input.pllps_threshold()
        high_pllps = df[df['p(LLPS)'] >= threshold]
        low_pllps = df[df['p(LLPS)'] < threshold]
        
        return ui.TagList(
            ui.p(f"High pLLPS proteins in dataset: {len(high_pllps)} (analyzing top {min(len(high_pllps), input.max_proteins())})"),
            ui.p(f"Low pLLPS proteins in dataset: {len(low_pllps)}")
        )
    
    @reactive.Effect
    @reactive.event(input.fetch_interactions)
    def fetch_interactions():
        df = data()
        req(df is not None)
        
        threshold = input.pllps_threshold()
        high_pllps_df = df[df['p(LLPS)'] >= threshold].sort_values('p(LLPS)', ascending=False)
        high_pllps_ids = high_pllps_df['Entry'].head(input.max_proteins()).tolist()
        
        req(len(high_pllps_ids) > 0)
        
        # Create progress message
        progress_msg = reactive.Value("Starting...")
        
        def update_progress(msg):
            progress_msg.set(msg)
        
        interactions, errors = fetch_string_interactions(
            high_pllps_ids, 
            score_threshold=input.string_score(),
            progress_callback=update_progress
        )
        
        if len(interactions) > 0:
            matched = match_interactions_to_pllps(interactions, df)
            interaction_data.set({
                'matched': matched, 
                'threshold': threshold,
                'errors': errors,
                'num_interactions': len(interactions)
            })
        else:
            interaction_data.set({
                'matched': pd.DataFrame(), 
                'threshold': threshold,
                'errors': errors + ["No interactions retrieved"],
                'num_interactions': 0
            })
    
    @output
    @render.ui
    def interaction_results():
        results = interaction_data()
        if results is None:
            return ui.div()
        
        matched_df = results['matched']
        threshold = results['threshold']
        errors = results.get('errors', [])
        num_interactions = results.get('num_interactions', 0)
        
        # Show errors if any
        error_display = []
        if errors:
            error_display.append(ui.div(
                ui.tags.div(
                    ui.h5("⚠️ Warnings/Errors:"),
                    ui.tags.ul(*[ui.tags.li(err) for err in errors[:5]]),  # Show first 5 errors
                    class_="alert alert-warning"
                )
            ))
        
        # Show info about fetched interactions
        if num_interactions > 0:
            error_display.append(ui.div(
                ui.tags.div(f"✅ Retrieved {num_interactions} interactions from STRING", class_="alert alert-info")
            ))
        
        if len(matched_df) == 0:
            return ui.TagList(
                *error_display,
                ui.div(ui.tags.div("No interaction data available. Please try fetching again or adjust parameters.", class_="alert alert-warning"))
            )
        
        enrichment = analyze_interaction_enrichment(matched_df, threshold)
        
        if enrichment is None:
            return ui.div(ui.tags.div("No complete interaction pairs for analysis", class_="alert alert-warning"))
        
        metrics = ui.row(
            ui.column(3, ui.div(ui.div("Total Interactions", class_="metric-title"), ui.div(str(enrichment['total']), class_="metric-value"), class_="value-box")),
            ui.column(3, ui.div(ui.div("High-High", class_="metric-title"), ui.div(f"{enrichment['high_high']} ({100*enrichment['high_high']/enrichment['total']:.1f}%)", class_="metric-value"), class_="value-box")),
            ui.column(3, ui.div(ui.div("High-Low", class_="metric-title"), ui.div(f"{enrichment['high_low']} ({100*enrichment['high_low']/enrichment['total']:.1f}%)", class_="metric-value"), class_="value-box")),
            ui.column(3, ui.div(ui.div("Low-Low", class_="metric-title"), ui.div(f"{enrichment['low_low']} ({100*enrichment['low_low']/enrichment['total']:.1f}%)", class_="metric-value"), class_="value-box")),
        )
        
        comparison_data = pd.DataFrame({
            'Type': ['High-High', 'High-Low', 'Low-Low'] * 2,
            'Count': [
                enrichment['high_high'], enrichment['high_low'], enrichment['low_low'],
                enrichment['expected_hh'] * enrichment['total'] / 100,
                enrichment['expected_hl'] * enrichment['total'] / 100,
                enrichment['expected_ll'] * enrichment['total'] / 100
            ],
            'Category': ['Observed']*3 + ['Expected']*3
        })
        
        fig = px.bar(comparison_data, x='Type', y='Count', color='Category', barmode='group', title="Observed vs Expected Interaction Counts")
        
        result_text = f"**Enrichment:** {enrichment['enrichment']:.2f}x"
        if enrichment['p_value'] is not None:
            result_text += f"\n\n**Chi-squared:** χ² = {enrichment['chi2']:.2f}, p = {enrichment['p_value']:.2e}"
            if enrichment['p_value'] < 0.05:
                if enrichment['enrichment'] > 1:
                    result_text += "\n\n✅ **Significant:** High pLLPS proteins preferentially interact!"
                else:
                    result_text += "\n\n✅ **Significant:** High pLLPS proteins avoid each other."
            else:
                result_text += "\n\n⚠️ **Not significant** (p >= 0.05)"
        
        return ui.TagList(
            *error_display,
            ui.hr(),
            ui.h4("📊 Enrichment Analysis Results"),
            metrics,
            ui.br(),
            ui.markdown(result_text),
            ui.HTML(fig.to_html(include_plotlyjs="cdn", div_id="enrichment_plot"))
        )
    
    @output
    @render.ui
    def export_controls():
        df = filtered_data()
        if df is None:
            return ui.div("No data to export")
        
        return ui.TagList(
            ui.download_button("download_csv", "📥 Download Filtered Data as CSV", class_="btn btn-primary"),
            ui.p(f"Rows: {len(df)}, Columns: {len(df.columns)}", class_="mt-3")
        )
    
    @session.download(filename="llps_proteins_filtered.csv")
    def download_csv():
        df = filtered_data()
        if df is not None:
            return df.to_csv(index=False)
        return ""


# Create the app
app = App(app_ui, server)
