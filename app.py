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
    
    with viz_tabs[5]:
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
        
        # Main content
        display_data_overview(filtered_df)
        
        st.markdown("---")
        
        # Data table with search
        filtered_df = display_data_table(filtered_df)
        
        st.markdown("---")
        
        # Visualizations
        display_visualizations(filtered_df)
        
        st.markdown("---")
        
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
