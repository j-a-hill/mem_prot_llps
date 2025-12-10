"""
LLPS Protein Data Explorer Dashboard

An interactive Streamlit dashboard for exploring Liquid-Liquid Phase Separation
(LLPS) protein data. Supports XLSX file upload and provides interactive 
visualizations for data exploration.

Usage:
    streamlit run app.py
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import re
import numpy as np
import requests
import time
from scipy import stats

# Page configuration
st.set_page_config(
    page_title="LLPS Protein Data Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Title and description
st.title("🧬 LLPS Protein Data Explorer")
st.markdown("""
This dashboard allows you to explore protein Liquid-Liquid Phase Separation (LLPS) data.
Upload your XLSX file or use the sample data to get started.
""")

def parse_location(location_str):
    """
    Parse a subcellular location string from UniProt and return a list of location terms.
    
    Process:
    1. Remove curly bracket annotations like {ECO:xxx} for cleaner output
    2. Remove normal bracket annotations like (By similarity) for cleaner output
    3. Remove isoform tags like [Isoform 1]: or Isoform 1:
    4. Remove everything after 'Note=' for cleaner output
    5. Split by common separators (comma, semicolon, period)
    6. Extract unique location tags for analysis and plotting
    
    The original column can be referenced for full annotation details if needed.
    """
    if pd.isna(location_str) or location_str == '':
        return []
    
    location_str = str(location_str)
    
    # Remove curly bracket content like {ECO:0000269|PubMed:12345}
    location_str = re.sub(r'\{[^}]*\}', '', location_str)
    
    # Remove all square bracket content, optionally followed by colon
    # This handles [Isoform 1]: and [Note: ...]
    location_str = re.sub(r'\[[^\]]*\]:?', '', location_str)
    
    # Remove Isoform ...: tags (without brackets)
    location_str = re.sub(r'Isoform\s+[^:]+:\s*', '', location_str, flags=re.IGNORECASE)
    
    # Remove normal bracket content like (By similarity) or (Potential)
    location_str = re.sub(r'\([^)]*\)', '', location_str)
    
    # Remove "SUBCELLULAR LOCATION:" prefix if present
    location_str = re.sub(r'^SUBCELLULAR LOCATION:\s*', '', location_str, flags=re.IGNORECASE)
    
    # Remove everything after 'Note=' (case-insensitive)
    location_str = re.sub(r'\s*Note=.*', '', location_str, flags=re.IGNORECASE)
    
    # Split by common separators and clean up
    # UniProt uses semicolons, commas, and periods as separators
    parts = re.split(r'[;,.]', location_str)
    
    locations = []
    for part in parts:
        # Clean up whitespace and filter empty/very short strings
        cleaned = part.strip()
        if len(cleaned) < 2:  # Skip empty or single-char remnants
            continue
        if cleaned not in locations:
            locations.append(cleaned)
    
    return locations


def add_location_columns(df):
    """Add parsed location columns to the dataframe."""
    if 'Subcellular location [CC]' not in df.columns:
        return df
    
    df = df.copy()
    # Add column with list of all parsed locations
    df['Location Categories'] = df['Subcellular location [CC]'].apply(parse_location)
    return df


def get_all_locations(df):
    """Get all unique location categories from the dataframe."""
    if 'Location Categories' not in df.columns:
        return []
    
    all_locations = set()
    for locations in df['Location Categories']:
        if isinstance(locations, list):
            all_locations.update(locations)
    
    # Return locations sorted alphabetically
    return sorted(all_locations)


# Predefined function categories for protein classification
# Based on UniProt's hierarchical keyword system for Molecular function and Biological process
# Reference: https://www.uniprot.org/keywords/KW-9992 (Molecular function)
# Reference: https://www.uniprot.org/keywords/KW-9999 (Biological process)
FUNCTION_CATEGORIES = {
    # =========================================================================
    # MOLECULAR FUNCTION (UniProt KW-9992)
    # =========================================================================
    
    # --- Enzyme activity (UniProt enzyme classification) ---
    'Hydrolase': [
        r'\bhydrolase\b',
        r'\besterase\b',
        r'\blipase\b',
        r'\bnuclease\b',
        r'\bglycosidase\b',
    ],
    'Kinase': [
        r'\bkinase\b',
        r'phosphotransferase',
        r'protein\s+kinase',
        r'serine/threonine[- ]protein\s+kinase',
        r'tyrosine[- ]protein\s+kinase',
    ],
    'Ligase': [
        r'\bligase\b',
        r'\bsynthetase\b',
        r'ubiquitin[- ]protein\s+ligase',
        r'e3\s+ligase',
    ],
    'Lyase': [
        r'\blyase\b',
        r'\bdecarboxylase\b',
        r'\bdehydratase\b',
        r'\baldolase\b',
    ],
    'Isomerase': [
        r'\bisomerase\b',
        r'\bracemase\b',
        r'\bmutase\b',
        r'\bepimerase\b',
    ],
    'Oxidoreductase': [
        r'oxidoreductase',
        r'\bdehydrogenase\b',
        r'\boxidase\b',
        r'\breductase\b',
        r'\bperoxidase\b',
        r'\bcatalase\b',
        r'\boxygenase\b',
    ],
    'Transferase': [
        r'\btransferase\b',
        r'methyltransferase',
        r'acetyltransferase',
        r'glycosyltransferase',
        r'aminotransferase',
    ],
    'Protease': [
        r'\bprotease\b',
        r'\bpeptidase\b',
        r'\bendopeptidase\b',
        r'\bexopeptidase\b',
        r'proteolytic',
        r'\bcaspase\b',
    ],
    'Phosphatase': [
        r'\bphosphatase\b',
        r'protein\s+phosphatase',
    ],
    
    # --- Nucleotide-binding proteins ---
    'ATP-binding': [
        r'atp[- ]binding',
        r'\batpase\b',
        r'atp\s+hydrolysis',
        r'abc\s+transporter',
    ],
    'GTP-binding': [
        r'gtp[- ]binding',
        r'\bgtpase\b',
        r'gtp\s+hydrolysis',
        r'small\s+gtpase',
        r'\bras\b',
        r'\brho\b',
        r'\brab\b',
    ],
    
    # --- Nucleic acid binding (UniProt KW-0238, KW-0694) ---
    'DNA-binding': [
        r'dna[- ]binding',
        r'binds\s+(to\s+)?dna',
        r'dna\s+binding',
        r'sequence[- ]specific\s+dna',
    ],
    'RNA-binding': [
        r'rna[- ]binding',
        r'binds\s+(to\s+)?rna',
        r'rna\s+binding',
        r'mrna[- ]binding',
    ],
    
    # --- Receptor activity (UniProt KW-0675) ---
    'Receptor': [
        r'\breceptor\b',
        r'receptor\s+activity',
        r'signal\s+receptor',
    ],
    'G-protein coupled receptor': [
        r'g[- ]protein[- ]coupled\s+receptor',
        r'\bgpcr\b',
        r'seven[- ]transmembrane',
        r'7[- ]?tm\s+receptor',
        r'rhodopsin[- ]like',
    ],
    
    # --- Ion channel activity (UniProt KW-0407) ---
    'Ion channel': [
        r'ion\s*channel',
        r'cation\s*channel',
        r'anion\s*channel',
        r'chloride\s*channel',
        r'sodium\s*channel',
        r'potassium\s*channel',
        r'calcium\s*channel',
    ],
    'Voltage-gated channel': [
        r'voltage[- ]gated',
        r'voltage[- ]dependent',
        r'voltage[- ]sensitive',
    ],
    'Ligand-gated channel': [
        r'ligand[- ]gated',
        r'receptor[- ]operated\s+channel',
        r'ionotropic\s+receptor',
    ],
    
    # --- Transporter activity (UniProt KW-0813) ---
    'Transporter': [
        r'\btransporter\b',
        r'transport\s+protein',
        r'carrier\s+protein',
        r'solute\s+carrier',
        r'\bslc\d',
    ],
    'Symporter': [
        r'\bsymporter\b',
        r'cotransporter',
        r'co[- ]transporter',
    ],
    'Antiporter': [
        r'\bantiporter\b',
        r'exchanger',
    ],
    
    # --- Chaperone activity (UniProt KW-0143) ---
    'Chaperone': [
        r'\bchaperone\b',
        r'chaperonin',
        r'heat\s+shock\s+protein',
        r'\bhsp\d',
        r'protein\s+folding',
    ],
    
    # --- Structural molecule activity ---
    'Structural protein': [
        r'structural\s+protein',
        r'structural\s+molecule',
        r'\bscaffold\b',
        r'cytoskelet',
    ],
    
    # --- Transcription regulation (UniProt KW-0805) ---
    'Transcription regulation': [
        r'transcription\s*(factor|regulator|regulation)',
        r'transcriptional\s*(activator|repressor|regulator)',
        r'dna[- ]binding\s+transcription',
    ],
    
    # =========================================================================
    # BIOLOGICAL PROCESS (UniProt KW-9999)
    # =========================================================================
    
    # --- Cell cycle (UniProt KW-0131) ---
    'Cell cycle': [
        r'cell\s+cycle',
        r'\bmitosis\b',
        r'\bmeiosis\b',
        r'cell\s+division',
        r'\bcytokinesis\b',
        r'chromosome\s+segregation',
    ],
    
    # --- Cell death (UniProt KW-0053) ---
    'Apoptosis': [
        r'\bapoptosis\b',
        r'\bapoptotic\b',
        r'programmed\s+cell\s+death',
        r'\bcaspase\b',
    ],
    
    # --- DNA damage/repair (UniProt KW-0227) ---
    'DNA damage': [
        r'dna\s+damage',
        r'dna\s+repair',
        r'double[- ]strand\s+break',
        r'mismatch\s+repair',
        r'base\s+excision\s+repair',
        r'nucleotide\s+excision',
    ],
    
    # --- DNA replication (UniProt KW-0235) ---
    'DNA replication': [
        r'dna\s+replication',
        r'dna\s+synthesis',
        r'\breplicon\b',
        r'replication\s+fork',
    ],
    
    # --- Transcription (UniProt KW-0804) ---
    'Transcription': [
        r'\btranscription\b',
        r'rna\s+polymerase',
        r'gene\s+expression',
        r'mrna\s+synthesis',
    ],
    
    # --- Translation (UniProt KW-0810) ---
    'Translation': [
        r'\btranslation\b',
        r'protein\s+biosynthesis',
        r'\bribosom',
        r'trna',
        r'mrna\s+translation',
    ],
    
    # --- Signal transduction (UniProt KW-0597) ---
    'Signal transduction': [
        r'signal\s+transduction',
        r'signaling\s+pathway',
        r'signalling\s+pathway',
        r'signal\s+cascade',
    ],
    
    # --- Immunity (UniProt KW-0391) ---
    'Immunity': [
        r'\bimmune\b',
        r'\bimmunity\b',
        r'innate\s+immun',
        r'adaptive\s+immun',
        r'\binflammation\b',
        r'\binterferon\b',
        r'\bcytokine\b',
        r'antigen\s+presentation',
    ],
    
    # --- Stress response (UniProt KW-0346) ---
    'Stress response': [
        r'stress\s+response',
        r'oxidative\s+stress',
        r'heat\s+shock',
        r'unfolded\s+protein\s+response',
    ],
    
    # --- Lipid metabolism (UniProt KW-0443) ---
    'Lipid metabolism': [
        r'lipid\s+metabol',
        r'fatty\s+acid',
        r'phospholipid',
        r'cholesterol',
        r'sphingolipid',
    ],
    
    # --- Carbohydrate metabolism (UniProt KW-0119) ---
    'Carbohydrate metabolism': [
        r'carbohydrate\s+metabol',
        r'\bglycolysis\b',
        r'gluconeogenesis',
        r'pentose\s+phosphate',
    ],
    
    # --- Protein modification (UniProt) ---
    'Ubl conjugation': [
        r'\bubiquitin',
        r'\bsumo',
        r'ubl\s+conjugation',
        r'protein\s+ubiquitination',
        r'sumoylation',
    ],
    'Phosphorylation': [
        r'\bphosphorylation\b',
        r'protein\s+phosphorylation',
    ],
    
    # --- Chromatin regulation (UniProt KW-0156) ---
    'Chromatin regulator': [
        r'\bchromatin\b',
        r'histone\s+modif',
        r'\bnucleosome\b',
        r'chromatin\s+remodel',
        r'histone\s+acetyl',
        r'histone\s+methyl',
    ],
    
    # --- mRNA processing (UniProt KW-0507) ---
    'mRNA processing': [
        r'mrna\s+processing',
        r'mrna\s+splicing',
        r'\bspliceosome\b',
        r'pre[- ]mrna',
        r'mrna\s+decay',
        r'deadenylation',
    ],
    
    # =========================================================================
    # CELLULAR COMPONENT / MEMBRANE (UniProt KW-9998)
    # =========================================================================
    'Transmembrane': [
        r'transmembrane',
        r'integral\s+membrane',
        r'multi[- ]?pass\s+membrane',
        r'single[- ]?pass\s+membrane',
    ],
}


def parse_function(function_str, protein_name_str=None):
    """
    Parse function annotations and protein names to extract functional categories.
    
    This function analyzes the 'Function [CC]' column and optionally the 'Protein names'
    column to identify and categorize proteins based on their function.
    
    Process:
    1. Combine function text and protein name for comprehensive matching
    2. Remove curly bracket annotations like {ECO:xxx} for cleaner matching
    3. Match against predefined functional category patterns
    4. Return a list of matching functional categories
    
    Args:
        function_str: The function annotation string from UniProt 'Function [CC]' column
        protein_name_str: Optional protein name string for additional context
    
    Returns:
        List of matched functional categories
    """
    categories = []
    
    # Combine function and protein name for matching
    text_parts = []
    if pd.notna(function_str) and function_str != '':
        text_parts.append(str(function_str))
    if pd.notna(protein_name_str) and protein_name_str != '':
        text_parts.append(str(protein_name_str))
    
    if not text_parts:
        return categories
    
    combined_text = ' '.join(text_parts)
    
    # Remove curly bracket annotations like {ECO:0000269|PubMed:12345}
    combined_text = re.sub(r'\{[^}]*\}', '', combined_text)
    
    # Convert to lowercase for case-insensitive matching
    combined_text_lower = combined_text.lower()
    
    # Match against each category's patterns
    for category, patterns in FUNCTION_CATEGORIES.items():
        for pattern in patterns:
            if re.search(pattern, combined_text_lower):
                if category not in categories:
                    categories.append(category)
                break  # Move to next category once matched
    
    return categories


def add_function_columns(df):
    """Add parsed function category columns to the dataframe."""
    df = df.copy()
    
    # Check for required columns
    has_function = 'Function [CC]' in df.columns
    has_protein_name = 'Protein names' in df.columns
    
    if not has_function and not has_protein_name:
        return df
    
    # Apply function parsing with both columns if available
    def extract_functions(row):
        func_str = row.get('Function [CC]') if has_function else None
        name_str = row.get('Protein names') if has_protein_name else None
        return parse_function(func_str, name_str)
    
    df['Function Categories'] = df.apply(extract_functions, axis=1)
    return df


def get_all_functions(df):
    """Get all unique function categories from the dataframe."""
    if 'Function Categories' not in df.columns:
        return []
    
    all_functions = set()
    for functions in df['Function Categories']:
        if isinstance(functions, list):
            all_functions.update(functions)
    
    # Return functions sorted alphabetically
    return sorted(all_functions)


@st.cache_data
def load_data(uploaded_file):
    """Load data from uploaded XLSX file."""
    return pd.read_excel(uploaded_file, engine='openpyxl')


@st.cache_data
def load_sample_data():
    """Load sample data if available."""
    sample_path = Path(__file__).parent / "data" / "sample_data.xlsx"
    if sample_path.exists():
        return pd.read_excel(sample_path, engine='openpyxl')
    return None


def display_data_overview(df):
    """Display overview statistics of the dataset."""
    st.header("📊 Data Overview")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Total Proteins", len(df))
    
    with col2:
        if 'Length' in df.columns and df['Length'].notna().any():
            st.metric("Avg. Protein Length", f"{df['Length'].mean():.0f}")
        else:
            st.metric("Rows", len(df))
    
    with col3:
        if 'p(LLPS)' in df.columns and df['p(LLPS)'].notna().any():
            st.metric("Avg. p(LLPS)", f"{df['p(LLPS)'].mean():.3f}")
        else:
            st.metric("Data Points", df.size)


def display_data_table(df):
    """Display interactive data table with filtering."""
    st.header("📋 Data Table")
    
    # Search functionality
    search_col1, search_col2 = st.columns([3, 1])
    
    with search_col1:
        search_term = st.text_input(
            "🔍 Search proteins", 
            placeholder="Enter protein name, entry, or keyword..."
        )
    
    with search_col2:
        search_columns = st.multiselect(
            "Search in columns",
            options=df.columns.tolist(),
            default=[col for col in ['Entry', 'Entry name', 'Protein names'] 
                    if col in df.columns]
        )
    
    # Filter data based on search
    filtered_df = df.copy()
    if search_term and search_columns:
        mask = pd.Series([False] * len(df))
        for col in search_columns:
            if col in df.columns:
                mask = mask | df[col].astype(str).str.contains(
                    search_term, case=False, na=False
                )
        filtered_df = df[mask]
    
    st.write(f"Showing {len(filtered_df)} of {len(df)} proteins")
    st.dataframe(filtered_df, use_container_width=True, height=400)
    
    return filtered_df


def display_visualizations(df):
    """Display interactive visualizations."""
    st.header("📈 Visualizations")
    
    # Create tabs for different visualization types
    viz_tabs = st.tabs([
        "Distribution", 
        "Scatter Plot", 
        "Location Analysis",
        "Function Analysis",
        "Combined Analysis",
        "Length Analysis",
        "Statistical Analysis"
    ])
    
    with viz_tabs[0]:
        st.subheader("Distribution of p(LLPS)")
        if 'p(LLPS)' in df.columns:
            fig = px.histogram(
                df, 
                x='p(LLPS)',
                nbins=50,
                title="Distribution of LLPS Probability",
                labels={'p(LLPS)': 'Probability of LLPS'},
                color_discrete_sequence=['#1f77b4']
            )
            fig.update_layout(
                xaxis_title="p(LLPS)",
                yaxis_title="Count",
                showlegend=False
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Column 'p(LLPS)' not found in data.")
    
    with viz_tabs[1]:
        st.subheader("Scatter Plot")
        
        col1, col2 = st.columns(2)
        numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
        
        if len(numeric_cols) >= 2:
            with col1:
                x_axis = st.selectbox(
                    "X-axis", 
                    options=numeric_cols,
                    index=numeric_cols.index('Length') if 'Length' in numeric_cols else 0
                )
            with col2:
                y_axis = st.selectbox(
                    "Y-axis",
                    options=numeric_cols,
                    index=numeric_cols.index('p(LLPS)') if 'p(LLPS)' in numeric_cols else min(1, len(numeric_cols)-1)
                )
            
            # Optional color coding
            color_options = ['None'] + df.columns.tolist()
            default_color_idx = 0
            
            color_col = st.selectbox(
                "Color by (optional)",
                options=color_options,
                index=default_color_idx
            )
            
            fig = px.scatter(
                df,
                x=x_axis,
                y=y_axis,
                color=None if color_col == 'None' else color_col,
                hover_data=['Entry', 'Entry name'] if all(c in df.columns for c in ['Entry', 'Entry name']) else None,
                title=f"{y_axis} vs {x_axis}"
            )
            fig.update_traces(marker=dict(size=8, opacity=0.7))
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Need at least 2 numeric columns for scatter plot.")
    
    with viz_tabs[2]:
        st.subheader("Subcellular Location Analysis")
        if 'Location Categories' in df.columns:
            # Explode Location Categories to count each location
            df_locations = df.explode('Location Categories')
            df_locations = df_locations[df_locations['Location Categories'].notna()]
            
            # Exclude membrane topology categories (Single-pass, Multi-pass)
            exclude_pattern = 'Single-pass|Multi-pass'
            df_locations = df_locations[~df_locations['Location Categories'].astype(str).str.contains(exclude_pattern, case=False, regex=True)]
            
            # Get top 20 locations by count (reused for both charts)
            location_counts = df_locations['Location Categories'].value_counts()
            top_20_locations = location_counts.head(20).index.tolist()
            location_counts_top20 = location_counts.head(20)
            
            # Location distribution bar chart - Top 20 locations by count
            fig = px.bar(
                x=location_counts_top20.index,
                y=location_counts_top20.values,
                title="Protein Count by Subcellular Location - Top 20 (proteins may appear in multiple)",
                labels={'x': 'Subcellular Location', 'y': 'Number of Proteins'},
                color=location_counts_top20.index,
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            fig.update_layout(showlegend=False, xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)
            
            # p(LLPS) by location box plot - Top 20 locations by count
            if 'p(LLPS)' in df.columns:
                df_locations_top20 = df_locations[df_locations['Location Categories'].isin(top_20_locations)]
                
                fig2 = px.box(
                    df_locations_top20,
                    x='Location Categories',
                    y='p(LLPS)',
                    title="p(LLPS) Distribution by Subcellular Location - Top 20",
                    color='Location Categories',
                    color_discrete_sequence=px.colors.qualitative.Set2
                )
                fig2.update_layout(showlegend=False, xaxis_tickangle=-45)
                st.plotly_chart(fig2, use_container_width=True)
            
            # Show location breakdown table
            with st.expander("📊 Location Statistics"):
                if 'p(LLPS)' in df.columns:
                    location_stats = df_locations.groupby('Location Categories').agg({
                        'Entry': 'count',
                        'p(LLPS)': ['mean', 'std', 'min', 'max']
                    }).round(3)
                    location_stats.columns = ['Count', 'Mean p(LLPS)', 'Std p(LLPS)', 'Min p(LLPS)', 'Max p(LLPS)']
                    location_stats = location_stats.sort_values('Count', ascending=False)
                    st.dataframe(location_stats, use_container_width=True)
                else:
                    location_stats = df_locations['Location Categories'].value_counts().to_frame('Count')
                    st.dataframe(location_stats, use_container_width=True)
        else:
            st.info("Column 'Subcellular location [CC]' not found in data.")
    
    with viz_tabs[3]:
        st.subheader("Protein Function Analysis")
        if 'Function Categories' in df.columns:
            # Explode Function Categories to count each function
            df_functions = df.explode('Function Categories')
            df_functions = df_functions[df_functions['Function Categories'].notna()]
            df_functions = df_functions[df_functions['Function Categories'] != '']
            
            if len(df_functions) > 0:
                # Get function counts
                function_counts = df_functions['Function Categories'].value_counts()
                top_20_functions = function_counts.head(20).index.tolist()
                function_counts_top20 = function_counts.head(20)
                
                # Function distribution bar chart - Top 20 functions by count
                fig = px.bar(
                    x=function_counts_top20.index,
                    y=function_counts_top20.values,
                    title="Top 20 Functional Categories by Protein Count",
                    labels={'x': 'Functional Category', 'y': 'Number of Proteins'},
                    color=function_counts_top20.index,
                    color_discrete_sequence=px.colors.qualitative.Dark24
                )
                fig.update_layout(showlegend=False, xaxis_tickangle=-45)
                st.plotly_chart(fig, use_container_width=True)
                
                # p(LLPS) by function box plot - Top 20 functions by count
                if 'p(LLPS)' in df.columns and len(top_20_functions) > 0:
                    df_functions_top20 = df_functions[df_functions['Function Categories'].isin(top_20_functions)]
                    
                    fig2 = px.box(
                        df_functions_top20,
                        x='Function Categories',
                        y='p(LLPS)',
                        title="p(LLPS) Distribution by Functional Category - Top 20",
                        color='Function Categories',
                        color_discrete_sequence=px.colors.qualitative.Dark24
                    )
                    fig2.update_layout(showlegend=False, xaxis_tickangle=-45)
                    st.plotly_chart(fig2, use_container_width=True)
                
                # Show function breakdown table
                with st.expander("📊 Function Statistics"):
                    if 'p(LLPS)' in df.columns:
                        function_stats = df_functions.groupby('Function Categories').agg({
                            'Entry': 'count',
                            'p(LLPS)': ['mean', 'std', 'min', 'max']
                        }).round(3)
                        function_stats.columns = ['Count', 'Mean p(LLPS)', 'Std p(LLPS)', 'Min p(LLPS)', 'Max p(LLPS)']
                        function_stats = function_stats.sort_values('Count', ascending=False)
                        st.dataframe(function_stats, use_container_width=True)
                    else:
                        function_stats = df_functions['Function Categories'].value_counts().to_frame('Count')
                        st.dataframe(function_stats, use_container_width=True)
                
                # Show available functional categories
                with st.expander("ℹ️ Available Function Categories"):
                    st.markdown("""
                    Functional categories are automatically detected from the 'Function [CC]' and 
                    'Protein names' columns using pattern matching. The following categories are recognized:
                    """)
                    
                    # Display categories in columns
                    category_list = sorted(FUNCTION_CATEGORIES.keys())
                    cols = st.columns(3)
                    for i, cat in enumerate(category_list):
                        cols[i % 3].write(f"• {cat}")
            else:
                st.info("No proteins with recognized functional categories found in the current selection.")
        else:
            st.info("Function categories not available. Ensure 'Function [CC]' or 'Protein names' columns exist in data.")
    
    with viz_tabs[4]:
        st.subheader("Combined Location & Function Analysis")
        st.markdown("""
        Explore p(LLPS) differences across combinations of subcellular location and protein function.
        For example, compare membrane-bound kinases vs cytosolic kinases.
        """)
        
        if 'Location Categories' in df.columns and 'Function Categories' in df.columns and 'p(LLPS)' in df.columns:
            # Get available locations and functions
            available_locations = get_all_locations(df)
            available_functions = get_all_functions(df)
            
            if available_locations and available_functions:
                # Selection UI
                st.markdown("### Select Categories to Compare")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    selected_locations_combined = st.multiselect(
                        "Subcellular Locations",
                        options=available_locations,
                        default=available_locations[:3] if len(available_locations) >= 3 else available_locations,
                        help="Select locations to include in the comparison",
                        key="combined_locations"
                    )
                
                with col2:
                    selected_functions_combined = st.multiselect(
                        "Functional Categories", 
                        options=available_functions,
                        default=available_functions[:3] if len(available_functions) >= 3 else available_functions,
                        help="Select functions to include in the comparison",
                        key="combined_functions"
                    )
                
                if selected_locations_combined and selected_functions_combined:
                    # Create combined location-function dataset
                    # Explode both location and function categories
                    df_combined = df.copy()
                    df_combined = df_combined.explode('Location Categories')
                    df_combined = df_combined.explode('Function Categories')
                    df_combined = df_combined[df_combined['Location Categories'].notna()]
                    df_combined = df_combined[df_combined['Function Categories'].notna()]
                    df_combined = df_combined[df_combined['Function Categories'] != '']
                    
                    # Filter to selected categories
                    df_combined = df_combined[
                        df_combined['Location Categories'].isin(selected_locations_combined) &
                        df_combined['Function Categories'].isin(selected_functions_combined)
                    ]
                    
                    # Warn if dataset is very large after explosion
                    if len(df_combined) > 10000:
                        st.warning(f"Large dataset ({len(df_combined)} data points). Visualizations may be slow.")
                    
                    if len(df_combined) > 0:
                        # Create combined label for grouping (using arrow for cleaner separation)
                        df_combined['Location-Function'] = df_combined['Location Categories'] + ' → ' + df_combined['Function Categories']
                        
                        # Grouped Box Plot
                        st.markdown("### p(LLPS) Distribution by Location and Function")
                        
                        fig_grouped = px.box(
                            df_combined,
                            x='Function Categories',
                            y='p(LLPS)',
                            color='Location Categories',
                            title="p(LLPS) by Function Category, Grouped by Location",
                            labels={
                                'Function Categories': 'Functional Category',
                                'Location Categories': 'Location',
                                'p(LLPS)': 'Probability of LLPS'
                            },
                            color_discrete_sequence=px.colors.qualitative.Set2
                        )
                        fig_grouped.update_layout(
                            xaxis_tickangle=-45,
                            legend=dict(
                                orientation="h",
                                yanchor="bottom",
                                y=1.02,
                                xanchor="right",
                                x=1
                            )
                        )
                        st.plotly_chart(fig_grouped, use_container_width=True)
                        
                        # Heatmap of mean p(LLPS)
                        st.markdown("### Mean p(LLPS) Heatmap: Location vs Function")
                        
                        # Create pivot table for heatmap
                        pivot_data = df_combined.pivot_table(
                            values='p(LLPS)',
                            index='Location Categories',
                            columns='Function Categories',
                            aggfunc='mean'
                        )
                        
                        if not pivot_data.empty:
                            fig_heatmap = px.imshow(
                                pivot_data,
                                title="Mean p(LLPS) by Location and Function",
                                labels=dict(x="Function", y="Location", color="Mean p(LLPS)"),
                                color_continuous_scale="Viridis",
                                aspect="auto"
                            )
                            fig_heatmap.update_layout(
                                xaxis_tickangle=-45
                            )
                            st.plotly_chart(fig_heatmap, use_container_width=True)
                        
                        # Count heatmap
                        st.markdown("### Protein Count Heatmap: Location vs Function")
                        
                        # Use groupby + size for robust counting regardless of null values
                        count_series = df_combined.groupby(
                            ['Location Categories', 'Function Categories']
                        ).size()
                        count_data = count_series.unstack(fill_value=0)
                        
                        if not count_data.empty:
                            fig_count = px.imshow(
                                count_data,
                                title="Protein Count by Location and Function",
                                labels=dict(x="Function", y="Location", color="Count"),
                                color_continuous_scale="Blues",
                                aspect="auto"
                            )
                            fig_count.update_layout(
                                xaxis_tickangle=-45
                            )
                            st.plotly_chart(fig_count, use_container_width=True)
                        
                        # Statistical summary table
                        with st.expander("📊 Detailed Statistics by Location-Function Combination"):
                            combo_stats = df_combined.groupby(['Location Categories', 'Function Categories']).agg(
                                Count=('p(LLPS)', 'size'),
                                Mean_pLLPS=('p(LLPS)', 'mean'),
                                Std_pLLPS=('p(LLPS)', 'std'),
                                Median_pLLPS=('p(LLPS)', 'median'),
                                Min_pLLPS=('p(LLPS)', 'min'),
                                Max_pLLPS=('p(LLPS)', 'max')
                            ).round(3)
                            combo_stats.columns = ['Count', 'Mean p(LLPS)', 'Std p(LLPS)', 'Median p(LLPS)', 'Min p(LLPS)', 'Max p(LLPS)']
                            combo_stats = combo_stats.sort_values('Mean p(LLPS)', ascending=False)
                            st.dataframe(combo_stats, use_container_width=True)
                        
                        # Violin plot for detailed distribution comparison
                        st.markdown("### Detailed Distribution Comparison")
                        
                        fig_violin = px.violin(
                            df_combined,
                            x='Location-Function',
                            y='p(LLPS)',
                            title="p(LLPS) Distribution by Location-Function Combination",
                            color='Location Categories',
                            box=True,
                            points="outliers",
                            color_discrete_sequence=px.colors.qualitative.Set2
                        )
                        fig_violin.update_layout(
                            xaxis_tickangle=-45,
                            showlegend=True,
                            legend=dict(
                                orientation="h",
                                yanchor="bottom",
                                y=1.02,
                                xanchor="right",
                                x=1
                            )
                        )
                        st.plotly_chart(fig_violin, use_container_width=True)
                        
                    else:
                        st.warning("No proteins match the selected location and function combination. Try selecting different categories.")
                else:
                    st.info("Select at least one location and one function category to see the combined analysis.")
            else:
                st.info("Both location and function categories are required for combined analysis.")
        else:
            st.info("Combined analysis requires 'Location Categories', 'Function Categories', and 'p(LLPS)' columns.")
    
    with viz_tabs[5]:
        st.subheader("Protein Length Analysis")
        if 'Length' in df.columns:
            fig = px.histogram(
                df,
                x='Length',
                nbins=50,
                title="Distribution of Protein Lengths",
                labels={'Length': 'Protein Length (amino acids)'},
                color_discrete_sequence=['#2ecc71']
            )
            st.plotly_chart(fig, use_container_width=True)
            
            # Length vs p(LLPS) correlation
            if 'p(LLPS)' in df.columns:
                # Check for sufficient valid data pairs
                valid_pairs = df[['Length', 'p(LLPS)']].dropna()
                if len(valid_pairs) >= 2:
                    st.write("**Correlation between Length and p(LLPS):**")
                    correlation = valid_pairs['Length'].corr(valid_pairs['p(LLPS)'])
                    if pd.notna(correlation):
                        st.write(f"Pearson correlation coefficient: {correlation:.4f}")
                    else:
                        st.write("Insufficient data variance for correlation calculation.")
        else:
            st.info("Column 'Length' not found in data.")
    
    with viz_tabs[6]:
        st.subheader("Statistical Analysis")
        
        # Correlation Analysis
        st.write("**📈 Correlation Analysis:**")
        numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
        
        if len(numeric_cols) >= 2:
            corr_matrix = df[numeric_cols].corr()
            st.write("Correlation matrix:")
            st.dataframe(corr_matrix.round(3), use_container_width=True)
            
            # Highlight specific correlations
            if 'Length' in numeric_cols and 'p(LLPS)' in numeric_cols:
                corr_length_pllps = df['Length'].corr(df['p(LLPS)'])
                if pd.notna(corr_length_pllps):
                    st.write(f"🔗 Length vs p(LLPS) correlation: {corr_length_pllps:.4f}")
        else:
            st.info("Not enough numeric columns for correlation analysis.")
        
        st.markdown("---")
        
        # Distribution Percentiles
        if 'p(LLPS)' in df.columns:
            st.write("**📊 p(LLPS) Distribution Percentiles:**")
            
            percentiles = [10, 25, 50, 75, 90, 95, 99]
            percentile_data = []
            for p in percentiles:
                value = df['p(LLPS)'].quantile(p/100)
                count_above = (df['p(LLPS)'] >= value).sum()
                percentile_data.append({
                    'Percentile': f"{p}th",
                    'Value': f"{value:.3f}",
                    'Proteins at or above': count_above
                })
            
            percentile_df = pd.DataFrame(percentile_data)
            st.dataframe(percentile_df, use_container_width=True, hide_index=True)
        else:
            st.info("Column 'p(LLPS)' not found for percentile analysis.")


def display_filtering_sidebar(df):
    """Display filtering options in sidebar."""
    st.sidebar.header("🔧 Filters")
    
    filtered_df = df.copy()
    
    # Numeric filters
    if 'p(LLPS)' in df.columns and df['p(LLPS)'].dropna().size > 0:
        pllps_min = float(df['p(LLPS)'].dropna().min())
        pllps_max = float(df['p(LLPS)'].dropna().max())
        pllps_range = st.sidebar.slider(
            "p(LLPS) Range",
            min_value=pllps_min,
            max_value=pllps_max,
            value=(pllps_min, pllps_max)
        )
        filtered_df = filtered_df[
            (filtered_df['p(LLPS)'] >= pllps_range[0]) & 
            (filtered_df['p(LLPS)'] <= pllps_range[1])
        ]
    
    if 'Length' in df.columns and df['Length'].dropna().size > 0:
        length_min = int(df['Length'].dropna().min())
        length_max = int(df['Length'].dropna().max())
        length_range = st.sidebar.slider(
            "Protein Length Range",
            min_value=length_min,
            max_value=length_max,
            value=(length_min, length_max)
        )
        filtered_df = filtered_df[
            (filtered_df['Length'] >= length_range[0]) & 
            (filtered_df['Length'] <= length_range[1])
        ]
    
    # Subcellular location filter
    if 'Location Categories' in df.columns:
        available_locations = get_all_locations(df)
        if available_locations:
            selected_locations = st.sidebar.multiselect(
                "Subcellular Location",
                options=available_locations,
                default=[],
                help="Filter proteins by subcellular location. Select one or more locations."
            )
            if selected_locations:
                # Filter by checking if any of the location categories match
                def location_matches(row):
                    if isinstance(row.get('Location Categories'), list):
                        return any(loc in selected_locations for loc in row['Location Categories'])
                    return False
                
                mask = filtered_df.apply(location_matches, axis=1)
                filtered_df = filtered_df[mask]
    
    # Function category filter
    if 'Function Categories' in df.columns:
        available_functions = get_all_functions(df)
        if available_functions:
            selected_functions = st.sidebar.multiselect(
                "Functional Category",
                options=available_functions,
                default=[],
                help="Filter proteins by functional category. Select one or more categories."
            )
            if selected_functions:
                # Filter by checking if any of the function categories match
                def function_matches(row):
                    if isinstance(row.get('Function Categories'), list):
                        return any(func in selected_functions for func in row['Function Categories'])
                    return False
                
                mask = filtered_df.apply(function_matches, axis=1)
                filtered_df = filtered_df[mask]
    
    # Disease involvement filter
    if 'Involvement in disease' in df.columns:
        disease_filter = st.sidebar.checkbox("Has disease involvement", value=False)
        if disease_filter:
            filtered_df = filtered_df[filtered_df['Involvement in disease'].notna() & 
                                      (filtered_df['Involvement in disease'] != '')]
    
    # PDB structure filter
    if 'Cross-reference (PDB)' in df.columns:
        pdb_filter = st.sidebar.checkbox("Has PDB structure", value=False)
        if pdb_filter:
            filtered_df = filtered_df[filtered_df['Cross-reference (PDB)'].notna() & 
                                      (filtered_df['Cross-reference (PDB)'] != '')]
    
    st.sidebar.markdown("---")
    st.sidebar.write(f"**Filtered proteins:** {len(filtered_df)}/{len(df)}")
    
    return filtered_df


def display_download_section(df):
    """Display download options for filtered data."""
    st.header("💾 Export Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # CSV download
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download as CSV",
            data=csv,
            file_name="llps_proteins_filtered.csv",
            mime="text/csv"
        )
    
    with col2:
        # Show column info
        with st.expander("Column Information"):
            for col in df.columns:
                dtype = df[col].dtype
                non_null = df[col].notna().sum()
                st.write(f"**{col}**: {dtype} ({non_null} non-null values)")


def fetch_string_interactions(protein_ids, species=9606, score_threshold=700, batch_size=100):
    """
    Fetch protein-protein interactions from STRING database.
    
    Args:
        protein_ids: List of UniProt IDs
        species: NCBI taxonomy ID (9606 = human)
        score_threshold: Minimum confidence (0-1000)
        batch_size: Proteins per request
    
    Returns:
        DataFrame with interactions
    """
    string_api_url = "https://string-db.org/api/json/network"
    all_interactions = []
    
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        status_text.text(f"Fetching batch {batch_num}/{total_batches}...")
        progress_bar.progress(batch_num / total_batches)
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "pllps_streamlit_app"
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
            elif response.status_code == 429:
                status_text.text("Rate limited. Waiting 30s...")
                time.sleep(30)
                response = requests.post(string_api_url, data=params, timeout=60)
                if response.status_code == 200:
                    interactions = response.json()
                    all_interactions.extend(interactions)
        except Exception as e:
            status_text.text(f"Error in batch {batch_num}: {str(e)}")
        
        time.sleep(1)  # Rate limiting
    
    progress_bar.empty()
    status_text.empty()
    
    return pd.DataFrame(all_interactions) if all_interactions else pd.DataFrame()


def match_interactions_to_pllps(interactions_df, pllps_df):
    """
    Match interaction partners to pLLPS dataset.
    
    Returns:
        DataFrame with both proteins' pLLPS values
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Create lookup dictionaries
    pllps_by_entry = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    if 'Entry name' in pllps_df.columns:
        pllps_by_name = dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)']))
    else:
        pllps_by_name = {}
    
    results = []
    for _, row in interactions_df.iterrows():
        protein_a = row.get('preferredName_A', row.get('stringId_A', ''))
        protein_b = row.get('preferredName_B', row.get('stringId_B', ''))
        score = row.get('score', 0)
        
        # Try multiple matching strategies
        pllps_a = pllps_by_entry.get(protein_a) or pllps_by_name.get(protein_a)
        pllps_b = pllps_by_entry.get(protein_b) or pllps_by_name.get(protein_b)
        
        results.append({
            'protein_a': protein_a,
            'protein_b': protein_b,
            'score': score,
            'pllps_a': pllps_a,
            'pllps_b': pllps_b
        })
    
    return pd.DataFrame(results)


def analyze_interaction_enrichment(matched_df, threshold=0.7):
    """
    Analyze high-high vs high-low interaction enrichment.
    
    Returns:
        dict with enrichment results
    """
    # Keep only interactions where both proteins have pLLPS scores
    complete = matched_df.dropna(subset=['pllps_a', 'pllps_b']).copy()
    
    if len(complete) == 0:
        return None
    
    # Classify interactions
    complete['class_a'] = np.where(complete['pllps_a'] >= threshold, 'High', 'Low')
    complete['class_b'] = np.where(complete['pllps_b'] >= threshold, 'High', 'Low')
    
    def interaction_type(row):
        if row['class_a'] == 'High' and row['class_b'] == 'High':
            return 'High-High'
        elif row['class_a'] == 'Low' and row['class_b'] == 'Low':
            return 'Low-Low'
        return 'High-Low'
    
    complete['interaction_type'] = complete.apply(interaction_type, axis=1)
    
    # Count
    counts = complete['interaction_type'].value_counts()
    total = len(complete)
    high_high = counts.get('High-High', 0)
    high_low = counts.get('High-Low', 0)
    low_low = counts.get('Low-Low', 0)
    
    # Calculate expected (null hypothesis)
    all_proteins = pd.concat([
        complete[['protein_a', 'pllps_a']].rename(columns={'protein_a': 'p', 'pllps_a': 'v'}),
        complete[['protein_b', 'pllps_b']].rename(columns={'protein_b': 'p', 'pllps_b': 'v'})
    ]).drop_duplicates(subset='p')
    
    n_high = (all_proteins['v'] >= threshold).sum()
    p_high = n_high / len(all_proteins) if len(all_proteins) > 0 else 0
    p_low = 1 - p_high
    
    expected_hh = p_high ** 2
    expected_hl = 2 * p_high * p_low
    expected_ll = p_low ** 2
    
    observed_hh = high_high / total if total > 0 else 0
    enrichment = observed_hh / expected_hh if expected_hh > 0 else 0
    
    # Chi-squared test
    observed = [high_high, high_low, low_low]
    expected = [expected_hh * total, expected_hl * total, expected_ll * total]
    
    p_value = None
    chi2 = None
    if all(e > 5 for e in expected):
        try:
            chi2, p_value = stats.chisquare(observed, expected)
        except:
            pass
    
    return {
        'high_high': high_high,
        'high_low': high_low,
        'low_low': low_low,
        'total': total,
        'enrichment': enrichment,
        'p_value': p_value,
        'chi2': chi2,
        'expected_hh': expected_hh * 100,
        'expected_hl': expected_hl * 100,
        'expected_ll': expected_ll * 100,
        'complete_df': complete
    }


def display_interaction_analysis(df):
    """Display protein interaction analysis section."""
    st.header("🔗 Protein Interaction Analysis")
    
    st.markdown("""
    Analyze protein-protein interactions from the STRING database to understand:
    - Whether high pLLPS proteins preferentially interact with each other
    - Network properties of high pLLPS proteins
    - Interaction patterns across pLLPS classes
    """)
    
    # Check if required columns exist
    if 'Entry' not in df.columns or 'p(LLPS)' not in df.columns:
        st.warning("Required columns 'Entry' and 'p(LLPS)' not found in dataset.")
        return
    
    # Settings
    st.subheader("⚙️ Analysis Settings")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        pllps_threshold = st.slider(
            "High pLLPS Threshold",
            min_value=0.5,
            max_value=0.95,
            value=0.7,
            step=0.05,
            help="Proteins with p(LLPS) >= this value are classified as 'High pLLPS'"
        )
    
    with col2:
        string_score = st.slider(
            "STRING Confidence Score",
            min_value=400,
            max_value=900,
            value=700,
            step=100,
            help="Minimum confidence for protein interactions (400=medium, 700=high, 900=highest)"
        )
    
    with col3:
        sample_size = st.number_input(
            "Max Proteins to Analyze",
            min_value=10,
            max_value=500,
            value=100,
            step=10,
            help="Number of top high pLLPS proteins to analyze (reduces API calls)"
        )
    
    # Get high pLLPS proteins
    high_pllps_df = df[df['p(LLPS)'] >= pllps_threshold].sort_values('p(LLPS)', ascending=False)
    high_pllps_ids = high_pllps_df['Entry'].head(sample_size).tolist()
    
    st.write(f"**High pLLPS proteins in dataset:** {len(high_pllps_df)} (analyzing top {len(high_pllps_ids)})")
    st.write(f"**Low pLLPS proteins in dataset:** {len(df) - len(high_pllps_df)}")
    
    # Button to fetch interactions
    if st.button("🚀 Fetch Interactions from STRING", type="primary"):
        if len(high_pllps_ids) == 0:
            st.error("No high pLLPS proteins found with the current threshold.")
            return
        
        with st.spinner("Fetching interactions from STRING database..."):
            interactions_df = fetch_string_interactions(
                high_pllps_ids,
                score_threshold=string_score,
                batch_size=100
            )
        
        if len(interactions_df) == 0:
            st.error("No interactions found. Check network connectivity or try different settings.")
            return
        
        st.success(f"✅ Retrieved {len(interactions_df)} interactions")
        
        # Match to pLLPS data
        with st.spinner("Matching interactions to pLLPS dataset..."):
            matched_df = match_interactions_to_pllps(interactions_df, df)
        
        both_matched = (matched_df['pllps_a'].notna() & matched_df['pllps_b'].notna()).sum()
        st.write(f"**Interactions with both proteins in dataset:** {both_matched}/{len(matched_df)} ({100*both_matched/len(matched_df):.1f}%)")
        
        # Store in session state
        st.session_state['interaction_results'] = {
            'matched_df': matched_df,
            'threshold': pllps_threshold
        }
    
    # Display results if available
    if 'interaction_results' in st.session_state:
        st.markdown("---")
        st.subheader("📊 Enrichment Analysis Results")
        
        matched_df = st.session_state['interaction_results']['matched_df']
        threshold = st.session_state['interaction_results']['threshold']
        
        results = analyze_interaction_enrichment(matched_df, threshold)
        
        if results is None:
            st.warning("No complete interaction pairs for analysis.")
            return
        
        # Display metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Interactions", results['total'])
        
        with col2:
            st.metric("High-High", f"{results['high_high']} ({100*results['high_high']/results['total']:.1f}%)")
        
        with col3:
            st.metric("High-Low", f"{results['high_low']} ({100*results['high_low']/results['total']:.1f}%)")
        
        with col4:
            st.metric("Low-Low", f"{results['low_low']} ({100*results['low_low']/results['total']:.1f}%)")
        
        # Enrichment visualization
        st.markdown("### Interaction Type Distribution")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Observed vs Expected
            comparison_data = pd.DataFrame({
                'Type': ['High-High', 'High-Low', 'Low-Low'] * 2,
                'Count': [
                    results['high_high'], results['high_low'], results['low_low'],
                    results['expected_hh'] * results['total'] / 100,
                    results['expected_hl'] * results['total'] / 100,
                    results['expected_ll'] * results['total'] / 100
                ],
                'Category': ['Observed']*3 + ['Expected']*3
            })
            
            fig = px.bar(
                comparison_data,
                x='Type',
                y='Count',
                color='Category',
                barmode='group',
                title="Observed vs Expected Interaction Counts",
                color_discrete_map={'Observed': '#2ecc71', 'Expected': '#95a5a6'}
            )
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            # Enrichment factor
            st.markdown("#### Enrichment Analysis")
            st.metric(
                "High-High Enrichment",
                f"{results['enrichment']:.2f}x",
                help="Ratio of observed to expected High-High interactions"
            )
            
            if results['p_value'] is not None:
                st.write(f"**Chi-squared test:**")
                st.write(f"χ² = {results['chi2']:.2f}")
                st.write(f"p-value = {results['p_value']:.2e}")
                
                if results['p_value'] < 0.05:
                    if results['enrichment'] > 1:
                        st.success("✅ **Significant:** High pLLPS proteins preferentially interact with each other!")
                    else:
                        st.info("✅ **Significant:** High pLLPS proteins tend to avoid each other.")
                else:
                    st.warning("⚠️ **Not significant:** No clear interaction preference detected (p >= 0.05)")
            else:
                st.info("⚠️ Chi-squared test not applicable (sample size too small)")
        
        # Network visualization
        st.markdown("### Interaction Network Visualization")
        
        # Create network data for visualization
        complete_df = results['complete_df']
        
        # Create a simple edge list visualization
        if len(complete_df) > 0:
            # Sample for visualization if too many
            if len(complete_df) > 200:
                st.info(f"Showing 200 of {len(complete_df)} interactions for visualization")
                viz_df = complete_df.sample(200, random_state=42)
            else:
                viz_df = complete_df
            
            # Create scatter plot colored by interaction type
            fig = px.scatter(
                viz_df,
                x='pllps_a',
                y='pllps_b',
                color='interaction_type',
                title="Interaction Pairs by pLLPS Values",
                labels={'pllps_a': 'Protein A p(LLPS)', 'pllps_b': 'Protein B p(LLPS)'},
                color_discrete_map={
                    'High-High': '#e74c3c',
                    'High-Low': '#f39c12',
                    'Low-Low': '#3498db'
                },
                hover_data=['protein_a', 'protein_b', 'score']
            )
            
            # Add threshold lines
            fig.add_hline(y=threshold, line_dash="dash", line_color="gray", opacity=0.5)
            fig.add_vline(x=threshold, line_dash="dash", line_color="gray", opacity=0.5)
            
            st.plotly_chart(fig, use_container_width=True)
        
        # Data export
        st.markdown("### 💾 Export Interaction Data")
        
        csv = complete_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Interaction Data (CSV)",
            data=csv,
            file_name="protein_interactions.csv",
            mime="text/csv"
        )


def main():
    """Main application function."""
    
    # Sidebar - File upload
    st.sidebar.header("📁 Data Source")
    
    data_source = st.sidebar.radio(
        "Choose data source:",
        ["Upload XLSX file", "Use sample data"]
    )
    
    df = None
    
    if data_source == "Upload XLSX file":
        uploaded_file = st.sidebar.file_uploader(
            "Upload your XLSX file",
            type=['xlsx'],
            help="Upload an Excel file containing protein LLPS data"
        )
        if uploaded_file is not None:
            try:
                df = load_data(uploaded_file)
                st.sidebar.success(f"✅ Loaded {len(df)} proteins")
            except Exception as e:
                st.sidebar.error(f"Error loading file: {str(e)}")
    else:
        df = load_sample_data()
        if df is not None:
            st.sidebar.success(f"✅ Sample data loaded ({len(df)} proteins)")
        else:
            st.sidebar.warning("No sample data available. Please upload a file.")
    
    if df is not None:
        # Add parsed location columns
        df = add_location_columns(df)
        # Add parsed function category columns
        df = add_function_columns(df)
        
        # Apply filters from sidebar
        filtered_df = display_filtering_sidebar(df)
        
        # Create main tabs
        main_tabs = st.tabs([
            "📊 Data Explorer",
            "🔗 Protein Interactions",
            "💾 Export"
        ])
        
        with main_tabs[0]:
            # Main content
            display_data_overview(filtered_df)
            
            st.markdown("---")
            
            # Data table with search
            filtered_df = display_data_table(filtered_df)
            
            st.markdown("---")
            
            # Visualizations
            display_visualizations(filtered_df)
        
        with main_tabs[1]:
            # Protein Interaction Analysis
            # Use full dataset (not filtered) for interaction analysis
            display_interaction_analysis(df)
        
        with main_tabs[2]:
            # Download section
            display_download_section(filtered_df)
    else:
        # Show instructions when no data is loaded
        st.info("""
        ### Getting Started
        
        1. **Upload your data**: Use the sidebar to upload your XLSX file containing protein LLPS data
        2. **Expected columns**: The dashboard works best with these columns:
           - `Entry` - UniProt entry ID
           - `Entry name` - UniProt entry name
           - `Protein names` - Full protein names
           - `p(LLPS)` - Probability of LLPS
           - `n(DPR=> 25)` - Number of dipeptide repeats
           - `Length` - Protein sequence length
           - `Function [CC]` - Function annotation
           - `Subcellular location [CC]` - Subcellular location
           - `Involvement in disease` - Disease associations
           - `Cross-reference (PDB)` - PDB structure references
        3. **Explore**: Use filters, search, and visualizations to explore your data
        """)
        
        # Show expected file format
        with st.expander("📝 Example Data Format"):
            example_data = pd.DataFrame({
                'Entry': ['P12345', 'Q67890'],
                'Entry name': ['PROT1_HUMAN', 'PROT2_HUMAN'],
                'Protein names': ['Example protein 1', 'Example protein 2'],
                'p(LLPS)': [0.85, 0.32],
                'n(DPR=> 25)': [5, 2],
                'Length': [450, 320]
            })
            st.dataframe(example_data)
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #888;'>
        <small>LLPS Protein Data Explorer | Built with Streamlit</small>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
