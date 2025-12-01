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
            color_col = st.selectbox(
                "Color by (optional)",
                options=['None'] + df.columns.tolist(),
                index=0
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
