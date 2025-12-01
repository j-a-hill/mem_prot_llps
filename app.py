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

# Common subcellular location keywords for matching
# These help match UniProt location annotations to standard categories
COMMON_LOCATION_KEYWORDS = [
    'Nucleus',
    'Cytoplasm', 
    'Membrane',
    'Mitochondrion',
    'Endoplasmic reticulum',
    'Golgi apparatus',
    'Secreted',
    'Extracellular',
    'Cell membrane',
    'Lysosome',
    'Peroxisome',
    'Cytoskeleton'
]


def parse_location(location_str):
    """
    Parse a subcellular location string from UniProt and return a list of location terms.
    
    Process:
    1. Split by comma/semicolon separators
    2. Match against common keywords OR preserve original term
    3. Return all unique locations found
    
    No 'Other' category - all terms are preserved for full traceability.
    """
    if pd.isna(location_str) or location_str == '':
        return []
    
    location_str = str(location_str)
    # Split by common separators and clean up
    parts = [p.strip() for p in location_str.replace(';', ',').split(',')]
    
    categories = []
    for part in parts:
        if not part:  # Skip empty strings (already stripped in list comprehension)
            continue
        part_lower = part.lower()
        matched = False
        # Try to match against common keywords
        for keyword in COMMON_LOCATION_KEYWORDS:
            if keyword.lower() in part_lower or part_lower in keyword.lower():
                if keyword not in categories:
                    categories.append(keyword)
                matched = True
                break
        # If no keyword match, preserve the original term
        if not matched and part not in categories:
            categories.append(part)
    
    return categories


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
    
    # Sort by predefined order, then alphabetically for any additional terms
    sorted_locations = []
    for loc in COMMON_LOCATION_KEYWORDS:
        if loc in all_locations:
            sorted_locations.append(loc)
    # Add any additional locations not in predefined list
    for loc in sorted(all_locations):
        if loc not in sorted_locations:
            sorted_locations.append(loc)
    
    return sorted_locations


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
        "Length Analysis"
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
            
            # Location distribution bar chart
            location_counts = df_locations['Location Categories'].value_counts()
            
            fig = px.bar(
                x=location_counts.index,
                y=location_counts.values,
                title="Protein Count by Subcellular Location (proteins may appear in multiple locations)",
                labels={'x': 'Subcellular Location', 'y': 'Number of Proteins'},
                color=location_counts.index,
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            fig.update_layout(showlegend=False, xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)
            
            # p(LLPS) by location box plot
            if 'p(LLPS)' in df.columns:
                fig2 = px.box(
                    df_locations,
                    x='Location Categories',
                    y='p(LLPS)',
                    title="p(LLPS) Distribution by Subcellular Location",
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
