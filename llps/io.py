"""
Sections 7 + 8 - Export, Caching, and Result Saving/Loading.

Functions for caching STRING interactions, exporting protein lists, and
saving/loading analysis results to/from disk.
"""

import pickle  # nosec B403 - only used for project-internal result files
import json
import pandas as pd
from pathlib import Path
from typing import List, Any

from llps.string_api import fetch_string_interactions  # noqa: F401 – re-exported alias


def save_interactions_to_cache(
    interactions_df: pd.DataFrame,
    score_threshold: int = 700,
    output_dir: str = "data"
) -> str:
    """
    Save STRING interactions to cache file for offline use.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame with interaction data from STRING
    score_threshold : int, optional
        Score threshold used (for filename)
    output_dir : str, optional
        Directory to save cache file (default: 'data')
    
    Returns
    -------
    str
        Path to saved cache file
    
    Examples
    --------
    >>> cache_path = save_interactions_to_cache(interactions_df, score_threshold=700)
    >>> print(f"Cache saved to: {cache_path}")
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    cache_file = output_path / f"string_cache_{score_threshold}.json"
    
    # Convert DataFrame to list of dicts for JSON
    interactions_list = interactions_df.to_dict('records')
    
    with open(cache_file, 'w') as f:
        json.dump(interactions_list, f, indent=2)
    
    return str(cache_file)


def export_protein_list(
    protein_ids: List[str],
    output_file: str = "high_pllps_proteins.txt"
) -> str:
    """
    Export protein IDs to a file for use with STRING web interface.
    
    Parameters
    ----------
    protein_ids : List[str]
        List of UniProt protein IDs.
    output_file : str
        Output filename.
    
    Returns
    -------
    str
        Path to the output file.
    """
    with open(output_file, 'w') as f:
        f.write('\n'.join(protein_ids))
    
    print(f"\n📄 Exported {len(protein_ids)} protein IDs to {output_file}")
    print(f"   Use this file with STRING web interface:")
    print(f"   1. Go to https://string-db.org/cgi/input")
    print(f"   2. Select 'Multiple proteins' tab")
    print(f"   3. Paste or upload the protein list")
    print(f"   4. Select 'Homo sapiens' as organism")
    print(f"   5. Download network as TSV for analysis")
    
    return output_file


def get_string_interactions(*args, **kwargs):
    """Alias for fetch_string_interactions (backward compatibility)"""
    return fetch_string_interactions(*args, **kwargs)


def save_analysis_result(data: Any, filename: str, results_dir: str = "results", 
                         format: str = "csv") -> Path:
    """
    Save analysis results to file.
    
    Parameters
    ----------
    data : Any
        Data to save (DataFrame, dict, list, etc.)
    filename : str
        Base filename (without extension)
    results_dir : str
        Directory to save results in (default: "results")
    format : str
        Format to save in: "csv", "json", "pickle" (default: "csv")
    
    Returns
    -------
    Path
        Path to saved file
    """
    results_path = Path(results_dir)
    results_path.mkdir(exist_ok=True)
    
    if format == "csv":
        filepath = results_path / f"{filename}.csv"
        if isinstance(data, pd.DataFrame):
            data.to_csv(filepath, index=False)
        else:
            raise ValueError("CSV format requires DataFrame input")
    elif format == "json":
        filepath = results_path / f"{filename}.json"
        with open(filepath, 'w') as f:
            if isinstance(data, (dict, list)):
                json.dump(data, f, indent=2, default=str)
            elif isinstance(data, pd.DataFrame):
                json.dump(data.to_dict(orient='records'), f, indent=2, default=str)
            else:
                raise ValueError("JSON format requires dict, list, or DataFrame")
    elif format == "pickle":
        filepath = results_path / f"{filename}.pkl"
        if isinstance(data, pd.DataFrame):
            data.to_pickle(filepath)
        else:
            with open(filepath, 'wb') as f:
                pickle.dump(data, f)  # nosec B301 - only stores project-internal result files
    else:
        raise ValueError(f"Unknown format: {format}")
    
    print(f"✅ Saved {format.upper()} to: {filepath}")
    return filepath


def load_analysis_result(filename: str, results_dir: str = "results", 
                         format: str = "csv") -> Any:
    """
    Load analysis results from file.
    
    Parameters
    ----------
    filename : str
        Base filename (with or without extension)
    results_dir : str
        Directory to load results from (default: "results")
    format : str
        Format to load: "csv", "json", "pickle" (default: "csv")
    
    Returns
    -------
    Any
        Loaded data (DataFrame, dict, list, etc.)
    """
    results_path = Path(results_dir)
    
    # Add extension if not present
    if not any(filename.endswith(ext) for ext in ['.csv', '.json', '.pkl']):
        if format == "csv":
            filename = f"{filename}.csv"
        elif format == "json":
            filename = f"{filename}.json"
        elif format == "pickle":
            filename = f"{filename}.pkl"
    
    filepath = results_path / filename
    
    if not filepath.exists():
        raise FileNotFoundError(f"Result file not found: {filepath}")
    
    if format == "csv" or filename.endswith('.csv'):
        data = pd.read_csv(filepath)
        print(f"✅ Loaded CSV from: {filepath} ({len(data)} rows)")
    elif format == "json" or filename.endswith('.json'):
        with open(filepath, 'r') as f:
            data = json.load(f)
        print(f"✅ Loaded JSON from: {filepath}")
    elif format == "pickle" or filename.endswith('.pkl'):
        with open(filepath, 'rb') as f:
            data = pickle.load(f)  # nosec B301 - only loads project-internal result files
        print(f"✅ Loaded pickle from: {filepath}")
    else:
        raise ValueError(f"Cannot determine format for: {filename}")
    
    return data


def list_saved_results(results_dir: str = "results") -> List[str]:
    """
    List all saved result files.
    
    Parameters
    ----------
    results_dir : str
        Directory containing results (default: "results")
    
    Returns
    -------
    List[str]
        List of result filenames
    """
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"⚠️  Results directory not found: {results_path}")
        return []
    
    files = []
    for ext in ['*.csv', '*.json', '*.pkl']:
        files.extend([f.name for f in results_path.glob(ext)])
    
    if files:
        print(f"📁 Found {len(files)} result files in {results_path}:")
        for f in sorted(files):
            file_path = results_path / f
            size = file_path.stat().st_size
            size_str = f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/(1024*1024):.1f} MB"
            print(f"   - {f} ({size_str})")
    else:
        print(f"📁 No result files found in {results_path}")
    
    return sorted(files)
