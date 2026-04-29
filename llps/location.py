"""
Section 2 - Location Parsing and Analysis.

Functions for parsing UniProt subcellular location strings and analysing
interaction preferences broken down by subcellular location.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Union, Iterable
import re

import pandas as pd
import requests


_SUBCELL_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/docs/subcell.txt"
)
_SUBCELL_CACHE: dict[str, "SubcellOntology"] = {}


def _clean_term_name(text: str) -> str:
    text = text.strip()
    if text.endswith("."):
        text = text[:-1]
    return text.strip()


def _normalize_key(text: str) -> str:
    text = _clean_term_name(text).lower()
    text = re.sub(r"[-_]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _split_synonyms(text: str) -> list[str]:
    parts = [p.strip() for p in text.split(";")]
    return [_clean_term_name(p) for p in parts if p]


@dataclass
class SubcellOntology:
    terms: dict[str, dict[str, Any]]
    name_to_accession: dict[str, str]
    children: dict[str, set[str]]
    _desc_cache: dict[str, set[str]]

    def lookup(self, name: str) -> str | None:
        return self.name_to_accession.get(_normalize_key(name))

    def get_descendants(self, accession: str, include_self: bool = True) -> set[str]:
        if accession in self._desc_cache:
            cached = self._desc_cache[accession]
            return set(cached) if include_self else set(cached) - {accession}

        descendants: set[str] = {accession}
        stack = list(self.children.get(accession, []))
        while stack:
            child = stack.pop()
            if child in descendants:
                continue
            descendants.add(child)
            stack.extend(self.children.get(child, []))
        self._desc_cache[accession] = descendants
        return set(descendants) if include_self else set(descendants) - {accession}


def download_subcell_terms(
    cache_path: str | Path | None = None,
    url: str = _SUBCELL_URL,
    force: bool = False,
) -> Path:
    if cache_path is None:
        cache_path = Path(__file__).parent.parent / "data" / "subcell.txt"
    path = Path(cache_path)
    if path.exists() and not force:
        return path

    path.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def _parse_subcell_txt(path: Path) -> SubcellOntology:
    terms: dict[str, dict[str, Any]] = {}
    parent_name_map: dict[str, set[str]] = {}
    name_to_accession: dict[str, str] = {}
    current: dict[str, Any] | None = None

    def _finalize(record: dict[str, Any] | None) -> None:
        if not record:
            return
        accession = record.get("accession")
        name = record.get("name")
        if not accession or not name:
            return
        terms[accession] = {
            "name": name,
            "synonyms": sorted(record.get("synonyms", set())),
            "parents": set(),
            "kind": record.get("kind", "location"),
        }
        parent_name_map[accession] = set(record.get("parent_names", set()))

        for key in [name, *record.get("synonyms", [])]:
            norm = _normalize_key(key)
            if norm and norm not in name_to_accession:
                name_to_accession[norm] = accession

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("//"):
                _finalize(current)
                current = None
                continue
            if line.startswith("ID") or line.startswith("IT") or line.startswith("IO"):
                _finalize(current)
                current = {
                    "name": _clean_term_name(line[5:]),
                    "kind": "topology" if line.startswith("IT") else "orientation" if line.startswith("IO") else "location",
                    "synonyms": set(),
                    "parent_names": set(),
                }
                continue
            if current is None:
                continue
            if line.startswith("AC"):
                current["accession"] = line[5:].strip()
            elif line.startswith("SY"):
                current["synonyms"].update(_split_synonyms(line[5:]))
            elif line.startswith("HI") or line.startswith("HP"):
                current["parent_names"].add(_clean_term_name(line[5:]))

    _finalize(current)

    for accession, parent_names in parent_name_map.items():
        parents = set()
        for parent_name in parent_names:
            parent_acc = name_to_accession.get(_normalize_key(parent_name))
            if parent_acc:
                parents.add(parent_acc)
        terms[accession]["parents"] = parents

    children: dict[str, set[str]] = {acc: set() for acc in terms}
    for acc, meta in terms.items():
        for parent_acc in meta["parents"]:
            children.setdefault(parent_acc, set()).add(acc)

    return SubcellOntology(
        terms=terms,
        name_to_accession=name_to_accession,
        children=children,
        _desc_cache={},
    )


def load_subcell_ontology(
    subcell_path: str | Path | None = None,
    allow_download: bool = True,
) -> SubcellOntology:
    if subcell_path is None:
        subcell_path = Path(__file__).parent.parent / "data" / "subcell.txt"
    path = Path(subcell_path)
    if not path.exists():
        if not allow_download:
            raise FileNotFoundError(f"Subcell vocabulary not found: {path}")
        path = download_subcell_terms(path)

    cache_key = str(path.resolve())
    if cache_key not in _SUBCELL_CACHE:
        _SUBCELL_CACHE[cache_key] = _parse_subcell_txt(path)
    return _SUBCELL_CACHE[cache_key]


def _clean_location_text(location_str: str) -> list[str]:
    location_str = str(location_str)
    location_str = re.sub(r"\{[^}]*\}", "", location_str)
    location_str = re.sub(r"\[[^\]]*\]:?", "", location_str)
    location_str = re.sub(r"Isoform\s+[^:]+:\s*", "", location_str, flags=re.IGNORECASE)
    location_str = re.sub(r"\([^)]*\)", "", location_str)
    location_str = re.sub(r"^SUBCELLULAR LOCATION:\s*", "", location_str, flags=re.IGNORECASE)
    location_str = re.sub(r"\s*Note=.*", "", location_str, flags=re.IGNORECASE)
    parts = re.split(r"[;,.]", location_str)
    tokens: list[str] = []
    for part in parts:
        cleaned = part.strip()
        if len(cleaned) < 2:
            continue
        tokens.append(_clean_term_name(cleaned))
    return tokens


def parse_location(
    location_str: Union[str, float, None],
    ontology: SubcellOntology | None = None,
    return_accessions: bool = False,
) -> List[str]:
    """
    Parse a subcellular location string from UniProt and return a list of location terms.

    Parameters
    ----------
    location_str : str
        Subcellular location string from UniProt
    ontology : SubcellOntology, optional
        Parsed UniProt subcellular location vocabulary.
    return_accessions : bool, optional
        If True, return SL accessions instead of display names.

    Returns
    -------
    list
        List of parsed location terms
    """
    if pd.isna(location_str) or location_str == "":
        return []

    ontology = ontology or load_subcell_ontology()
    tokens = _clean_location_text(location_str)

    accessions: list[str] = []
    for token in tokens:
        acc = ontology.lookup(token)
        if acc and acc not in accessions:
            accessions.append(acc)

    if return_accessions:
        return accessions

    names = [ontology.terms[acc]["name"] for acc in accessions if acc in ontology.terms]
    return names


def is_membrane_localized(
    location_ids: Iterable[str] | None,
    ontology: SubcellOntology | None = None,
    membrane_root: str = "Membrane",
) -> bool:
    if not location_ids:
        return False

    ontology = ontology or load_subcell_ontology()
    root_acc = ontology.lookup(membrane_root)
    if not root_acc:
        return False

    membrane_desc = ontology.get_descendants(root_acc, include_self=True)
    return any(loc in membrane_desc for loc in location_ids)


def categorize_location_to_compartment(
    location_str: Union[str, float, None],
    is_membrane: bool = False
) -> str:
    """
    Categorize a UniProt subcellular location string into major cellular compartments.

    Uses the is_membrane flag to disambiguate locations like "Endoplasmic reticulum"
    which could mean the membrane or the lumen. This avoids overcounting membrane proteins.

    Parameters
    ----------
    location_str : str
        Subcellular location string from UniProt
    is_membrane : bool, optional
        Whether the protein is classified as a membrane protein. Default is False.
        If True and location contains membrane-associated keywords, categorizes accordingly.

    Returns
    -------
    str
        Major compartment category (e.g., 'Plasma Membrane', 'Cytosol', 'ER Lumen', 'Mitochondrion')

    Notes
    -----
    - Compartments ending in "membrane" are only classified as membrane if is_membrane=True
    - Proteins in ER/Golgi/etc without the membrane flag are categorized as "_Lumen"
    - This prevents overcounting of membrane proteins
    """
    if pd.isna(location_str) or location_str == '':
        return 'Unknown'

    location_str = str(location_str).lower()

    # Plasma/cell membrane - check first for membrane proteins (takes precedence over cytosolic keywords)
    if any(term in location_str for term in ['plasma membrane', 'cell membrane']):
        return 'Plasma Membrane'

    # Check for cytosolic/nuclear proteins (typically non-membrane, but skip if membrane protein with other membrane locations)
    if any(term in location_str for term in ['cytoplasm', 'cytosol', 'nucleus', 'nucleoplasm']):
        # If it's a membrane protein, check if it also has other membrane locations to use instead
        if is_membrane and any(term in location_str for term in ['membrane', 'endoplasmic reticulum', 'er ', 'golgi', 'mitochondri', 'peroxisom', 'lysosom', 'vacuole']):
            # Has both cytoplasm/nucleus AND a specific membrane location, continue to check for more specific compartment
            pass
        else:
            return 'Cytosol'

    # Mitochondrion - check if membrane or matrix/lumen
    if 'mitochondri' in location_str:
        if is_membrane:
            return 'Mitochondrial Membrane'
        else:
            return 'Mitochondrial Matrix'

    # Endoplasmic reticulum - disambiguate with is_membrane flag
    if 'endoplasmic reticulum' in location_str or 'er ' in location_str or location_str.endswith('er'):
        if is_membrane:
            return 'ER Membrane'
        else:
            return 'ER Lumen'

    # Golgi apparatus - disambiguate with is_membrane flag
    if 'golgi' in location_str:
        if is_membrane:
            return 'Golgi Membrane'
        else:
            return 'Golgi Lumen'

    # Peroxisome
    if 'peroxisom' in location_str:
        if is_membrane:
            return 'Peroxisomal Membrane'
        else:
            return 'Peroxisomal Matrix'

    # Lysosome/Vacuole
    if any(term in location_str for term in ['lysosom', 'vacuole']):
        if is_membrane:
            return 'Lysosomal Membrane'
        else:
            return 'Lysosomal Lumen'

    # General membrane (but not plasma) - only if is_membrane=True
    if 'membrane' in location_str:
        if is_membrane:
            return 'Other Membrane'
        else:
            return 'Membrane-Associated Lumen'

    # If nothing matched but is_membrane is True, probably a membrane protein
    if is_membrane:
        return 'Other Membrane'

    # Final fallback for non-membrane proteins that had cytosolic keywords but no specific compartment found
    if any(term in location_str for term in ['cytoplasm', 'cytosol', 'nucleus', 'nucleoplasm']):
        return 'Cytosol'

    return 'Other'


def add_location_columns(
    df: pd.DataFrame,
    ontology: SubcellOntology | None = None,
    add_accessions: bool = False,
    accessions_col: str = "Location_IDs",
) -> pd.DataFrame:
    """
    Add parsed location columns to the dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'Subcellular location [CC]' column

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'Location Categories' column
    """
    if 'Subcellular location [CC]' not in df.columns:
        return df

    df = df.copy()
    ontology = ontology or load_subcell_ontology()

    df["Location Categories"] = df["Subcellular location [CC]"].apply(
        lambda value: parse_location(value, ontology=ontology, return_accessions=False)
    )

    if add_accessions:
        df[accessions_col] = df["Subcellular location [CC]"].apply(
            lambda value: parse_location(value, ontology=ontology, return_accessions=True)
        )

    return df
