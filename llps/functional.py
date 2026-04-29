"""
Section 9 - Functional Classification Utilities.

Functions for classifying proteins by functional category using GO IDs
and for identifying membrane proteins using the UniProt subcell ontology.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import List, Dict, Iterable

import pandas as pd
import requests
import yaml


# Cache for loaded functional terms
_FUNCTIONAL_TERMS_CACHE = None
_GO_DAG_CACHE: dict[str, object] = {}
_GO_SLIM_DAG_CACHE: dict[str, object] = {}
_GO_DESC_CACHE: dict[str, set[str]] = {}
_GO_SLIM_CACHE: dict[str, tuple[set[str], set[str]]] = {}
_GROUP_DESC_CACHE: dict[str, dict[str, set[str]]] = {}

_GO_BASIC_URL = "https://current.geneontology.org/ontology/go-basic.obo"
_GO_SLIM_GENERIC_URL = "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"


def _load_functional_terms(yaml_path: str | None = None) -> Dict:
    """
    Load functional classification terms from YAML file.
    
    Parameters
    ----------
    yaml_path : str, optional
        Path to YAML file. If None, uses default location.
    
    Returns
    -------
    Dict
        Dictionary with functional_groups and GO term IDs
    """
    global _FUNCTIONAL_TERMS_CACHE
    
    if _FUNCTIONAL_TERMS_CACHE is not None:
        return _FUNCTIONAL_TERMS_CACHE
    
    if yaml_path is None:
        # Try to find the default YAML file
        possible_paths = [
            Path('data/functional_classification_terms.yaml'),
            Path(__file__).parent.parent / 'data' / 'functional_classification_terms.yaml',
        ]
        
        for path in possible_paths:
            if path.exists():
                yaml_path = path
                break
        
        if yaml_path is None:
            raise FileNotFoundError(
                "Could not find functional_classification_terms.yaml. "
                "Please provide the path explicitly."
            )
    
    with open(yaml_path, 'r') as f:
        terms = yaml.safe_load(f)
    
    _FUNCTIONAL_TERMS_CACHE = terms
    return terms


def _download_obo(url: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def _load_go_dag(go_obo_path: str | Path | None = None, allow_download: bool = True):
    if go_obo_path is None:
        go_obo_path = Path(__file__).parent.parent / "data" / "go" / "go-basic.obo"
    path = Path(go_obo_path)
    if not path.exists():
        if not allow_download:
            raise FileNotFoundError(f"GO OBO file not found: {path}")
        _download_obo(_GO_BASIC_URL, path)

    cache_key = str(path.resolve())
    if cache_key not in _GO_DAG_CACHE:
        try:
            from goatools.obo_parser import GODag
        except ImportError as exc:
            raise ImportError("goatools is required for GO-based classification") from exc
        _GO_DAG_CACHE[cache_key] = GODag(str(path))
    return _GO_DAG_CACHE[cache_key]


def _load_go_slim_dag(
    slim_obo_path: str | Path | None = None,
    allow_download: bool = True,
    slim_url: str = _GO_SLIM_GENERIC_URL,
):
    if slim_obo_path is None:
        slim_obo_path = Path(__file__).parent.parent / "data" / "go" / "goslim_generic.obo"
    path = Path(slim_obo_path)
    if not path.exists():
        if not allow_download:
            raise FileNotFoundError(f"GO Slim OBO file not found: {path}")
        _download_obo(slim_url, path)

    cache_key = str(path.resolve())
    if cache_key not in _GO_SLIM_DAG_CACHE:
        try:
            from goatools.obo_parser import GODag
        except ImportError as exc:
            raise ImportError("goatools is required for GO Slim mapping") from exc
        _GO_SLIM_DAG_CACHE[cache_key] = GODag(str(path))
    return _GO_SLIM_DAG_CACHE[cache_key]


def parse_go_ids(go_ids: Iterable[str] | str | None) -> List[str]:
    if go_ids is None or (isinstance(go_ids, float) and pd.isna(go_ids)):
        return []
    if isinstance(go_ids, (list, tuple, set)):
        return sorted({str(go_id) for go_id in go_ids if go_id})
    text = str(go_ids)
    return sorted(set(re.findall(r"GO:\d{7}", text)))


def _expand_go_descendants(go_ids: Iterable[str], go_dag) -> set[str]:
    expanded: set[str] = set()
    for go_id in go_ids:
        if go_id in _GO_DESC_CACHE:
            expanded |= _GO_DESC_CACHE[go_id]
            continue
        if go_id not in go_dag:
            _GO_DESC_CACHE[go_id] = {go_id}
            expanded.add(go_id)
            continue
        desc = set(go_dag[go_id].get_all_children())
        desc.add(go_id)
        _GO_DESC_CACHE[go_id] = desc
        expanded |= desc
    return expanded


def map_go_ids_to_slim(
    go_ids: Iterable[str] | str | None,
    go_dag=None,
    slim_dag=None,
    use_direct: bool = True,
) -> List[str]:
    go_ids_parsed = parse_go_ids(go_ids)
    if not go_ids_parsed:
        return []

    go_dag = go_dag or _load_go_dag()
    slim_dag = slim_dag or _load_go_slim_dag()

    try:
        import goatools.mapslim as mapslim
    except ImportError as exc:
        raise ImportError("goatools is required for GO Slim mapping") from exc

    slim_ids: set[str] = set()
    for go_id in go_ids_parsed:
        if go_id not in go_dag:
            continue
        if go_id in _GO_SLIM_CACHE:
            direct_ids, all_ids = _GO_SLIM_CACHE[go_id]
        else:
            direct_ids, all_ids = mapslim.mapslim(go_id, go_dag, slim_dag)
            _GO_SLIM_CACHE[go_id] = (direct_ids, all_ids)
        slim_ids |= direct_ids if use_direct else all_ids

    return sorted(slim_ids)


def _get_group_descendants(yaml_path: str | None, go_dag) -> dict[str, set[str]]:
    key = f"{yaml_path or 'default'}::{id(go_dag)}"
    if key in _GROUP_DESC_CACHE:
        return _GROUP_DESC_CACHE[key]

    terms = _load_functional_terms(yaml_path)
    functional_groups = terms.get("functional_groups", {})
    group_desc: dict[str, set[str]] = {}
    for category, group_data in functional_groups.items():
        parent_ids = group_data.get("go_ids", [])
        group_desc[category] = _expand_go_descendants(parent_ids, go_dag)

    _GROUP_DESC_CACHE[key] = group_desc
    return group_desc


def parse_function_categories(
    function_str: str | None = None,
    protein_name_str: str | None = None,
    yaml_path: str | None = None,
    go_ids: Iterable[str] | str | None = None,
    go_dag=None,
) -> List[str]:
    """
    Classify proteins into functional categories using GO IDs.

    Parameters
    ----------
    function_str : str, optional
        Ignored unless it contains GO IDs (kept for backward compatibility).
    protein_name_str : str, optional
        Deprecated, retained for signature compatibility.
    yaml_path : str, optional
        Path to YAML file with GO-driven functional groups.
    go_ids : iterable or str, optional
        GO IDs for the protein (list or a string containing GO:XXXXXXX terms).
    go_dag : GODag, optional
        Pre-loaded GO DAG.

    Returns
    -------
    List[str]
        List of functional categories that match the protein.
    """
    if go_ids is None and function_str and "GO:" in str(function_str):
        go_ids = function_str

    go_ids_parsed = parse_go_ids(go_ids)
    if not go_ids_parsed:
        warnings.warn("GO IDs not provided; returning empty functional categories.")
        return []

    go_dag = go_dag or _load_go_dag()
    group_desc = _get_group_descendants(yaml_path, go_dag)

    categories: list[str] = []
    for category, descendants in group_desc.items():
        if any(go_id in descendants for go_id in go_ids_parsed):
            categories.append(category)

    return categories


def is_membrane_protein(
    function_str: str,
    protein_name_str: str = None,
    location_str: str = None,
    yaml_path: str = None,
    tmd_count: int | None = None,
    use_function: bool = False,
    location_ids: Iterable[str] | None = None,
    ontology=None,
) -> bool:
    """
    Check if a protein is a membrane protein based on annotations.
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation
    protein_name_str : str, optional
        Protein names
    location_str : str, optional
        Subcellular location [CC] annotation
    yaml_path : str, optional
        Deprecated; retained for signature compatibility
    location_ids : iterable, optional
        UniProt SL accessions for location terms
    ontology : SubcellOntology, optional
        Parsed UniProt subcellular location ontology
    
    Returns
    -------
    bool
        True if protein is classified as a membrane protein
    """
    location_ids_final: list[str] = []
    if location_ids is not None:
        location_ids_final = list(location_ids)
    elif location_str:
        try:
            from llps.location import parse_location as _parse_location
        except Exception:
            _parse_location = None
        if _parse_location is not None:
            location_ids_final = _parse_location(
                location_str,
                ontology=ontology,
                return_accessions=True,
            )

    if location_ids_final:
        try:
            from llps.location import is_membrane_localized
        except Exception:
            is_membrane_localized = None
        if is_membrane_localized is not None:
            if is_membrane_localized(location_ids_final, ontology=ontology):
                return True

    if tmd_count is not None and isinstance(tmd_count, int) and tmd_count > 0:
        return True

    return False


def count_tm_domains(domain_str: "str | float") -> int:
    """Count domain spans (e.g. '37..60') in a UniProt Transmembrane/Intramembrane string."""
    if pd.isna(domain_str) or not domain_str:
        return 0
    return len(re.findall(r'\d+\.\.\d+', str(domain_str)))


def classify_protein_function(
    function_str: str | None = None,
    protein_name_str: str | None = None,
    go_ids: Iterable[str] | str | None = None,
) -> List[str]:
    """
    Classify protein into functional categories based on annotations.
    
    DEPRECATED: Use parse_function_categories() instead.
    Kept for backward compatibility.
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation
    protein_name_str : str, optional
        Protein names
    
    Returns
    -------
    List[str]
        List of functional categories the protein belongs to
    """
    return parse_function_categories(
        function_str=function_str,
        protein_name_str=protein_name_str,
        go_ids=go_ids,
    )


def filter_membrane_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter DataFrame to only membrane proteins.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with protein data including 'Function [CC]', 'Protein names',
        and 'Subcellular location [CC]' columns
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only membrane proteins
    """
    def check_membrane(row):
        # Prefer explicit transmembrane annotation / TMD_count when available
        tmd = None
        if 'TMD_count' in row.index and pd.notna(row.get('TMD_count')):
            try:
                tmd = int(row.get('TMD_count') or 0)
            except Exception:
                tmd = None
        else:
            # Try to infer from Transmembrane / Intramembrane string columns
            tmd_str = row.get('Transmembrane', '') or row.get('Intramembrane', '')
            if isinstance(tmd_str, str) and tmd_str:
                tmd = count_tm_domains(tmd_str)

        return is_membrane_protein(
            row.get('Function [CC]', ''),
            row.get('Protein names', ''),
            row.get('Subcellular location [CC]', ''),
            yaml_path=None,
            tmd_count=tmd,
            use_function=False,
        )
    
    membrane_mask = df.apply(check_membrane, axis=1)
    membrane_df = df[membrane_mask].copy()
    
    print(f"🔍 Filtered to {len(membrane_df)} membrane proteins ({len(membrane_df)/len(df)*100:.1f}%)")
    
    return membrane_df


def add_membrane_flag(df: pd.DataFrame, col_name: str = 'Is_Membrane') -> pd.DataFrame:
    """
    Add a boolean membrane flag column to the DataFrame using available
    TMD annotations and functional/location text.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    col_name : str
        Column name to add (default: 'Is_Membrane')

    Returns
    -------
    pd.DataFrame
        Copy of df with the membrane flag column added
    """
    if df is None or not isinstance(df, pd.DataFrame):
        return df

    df = df.copy()

    def _flag_row(row):
        tmd = None
        if 'TMD_count' in row.index and pd.notna(row.get('TMD_count')):
            try:
                tmd = int(row.get('TMD_count') or 0)
            except Exception:
                tmd = None
        else:
            tmd_str = row.get('Transmembrane', '') or row.get('Intramembrane', '')
            if isinstance(tmd_str, str) and tmd_str:
                tmd = count_tm_domains(tmd_str)

        return is_membrane_protein(
            row.get('Function [CC]', ''),
            row.get('Protein names', ''),
            row.get('Subcellular location [CC]', ''),
            yaml_path=None,
            tmd_count=tmd,
            use_function=False,
        )

    df[col_name] = df.apply(_flag_row, axis=1)
    return df


def add_functional_categories(
    df: pd.DataFrame,
    go_col: str | None = None,
    yaml_path: str | None = None,
    go_dag=None,
) -> pd.DataFrame:
    """
    Add functional category columns to DataFrame.
    
    Similar to add_location_columns(), this function:
    1. Parses GO ID fields
    2. Adds a 'Functional_Categories' column with list of categories
    3. Optionally adds binary columns for each category
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a GO ID column (e.g., 'GO_IDs' or 'Gene Ontology IDs')
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added 'Functional_Categories' column
    
    Examples
    --------
    >>> df_with_categories = add_functional_categories(df)
    >>> print(df_with_categories['Functional_Categories'].head())
    """
    if go_col is None:
        for candidate in ("GO_IDs", "Gene Ontology IDs"):
            if candidate in df.columns:
                go_col = candidate
                break

    if not go_col:
        print("Warning: no GO ID column found; Functional_Categories not added")
        return df

    df = df.copy()
    go_dag = go_dag or _load_go_dag()
    group_desc = _get_group_descendants(yaml_path, go_dag)

    def _classify(go_value):
        go_list = parse_go_ids(go_value)
        if not go_list:
            return []
        categories: list[str] = []
        for category, descendants in group_desc.items():
            if any(go_id in descendants for go_id in go_list):
                categories.append(category)
        return categories

    df["Functional_Categories"] = df[go_col].apply(_classify)
    
    # Add binary columns for each category
    all_categories = set()
    for cats in df['Functional_Categories']:
        all_categories.update(cats)
    
    for category in sorted(all_categories):
        col_name = f'Is_{category.replace(" ", "_")}'
        df[col_name] = df['Functional_Categories'].apply(lambda x: category in x)
    
    print(f"Added functional categories to {len(df)} proteins")
    print(f"Categories found: {sorted(all_categories)}")
    print(f"Total category assignments: {sum(len(cats) for cats in df['Functional_Categories'])}")

    return df


def add_go_slim_categories(
    df: pd.DataFrame,
    go_col: str | None = None,
    output_col: str = "GO_Slim_Categories",
    go_dag=None,
    slim_dag=None,
    use_direct: bool = True,
) -> pd.DataFrame:
    if go_col is None:
        for candidate in ("GO_IDs", "Gene Ontology IDs"):
            if candidate in df.columns:
                go_col = candidate
                break

    if not go_col:
        print("Warning: no GO ID column found; GO Slim categories not added")
        return df

    df = df.copy()
    go_dag = go_dag or _load_go_dag()
    slim_dag = slim_dag or _load_go_slim_dag()

    def _slim(go_value):
        slim_ids = map_go_ids_to_slim(go_value, go_dag=go_dag, slim_dag=slim_dag, use_direct=use_direct)
        return [slim_dag[go_id].name for go_id in slim_ids if go_id in slim_dag]

    df[output_col] = df[go_col].apply(_slim)
    return df
