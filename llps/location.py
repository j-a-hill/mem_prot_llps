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
    is_membrane: bool = False,
    ontology: "SubcellOntology | None" = None,
) -> str:
    """
    Categorize a UniProt subcellular location string into major cellular compartments.

    Uses the SubcellOntology hierarchy to resolve location terms, then walks ancestry to
    assign the protein to a major compartment.  The is_membrane flag disambiguates organelles
    that have both membrane and lumenal sub-compartments (e.g. Mitochondrion → Mitochondrial
    Membrane vs Mitochondrial Matrix).

    Parameters
    ----------
    location_str : str
        Subcellular location string from UniProt
    is_membrane : bool, optional
        Whether the protein has transmembrane domain annotation. Default False.
    ontology : SubcellOntology, optional
        Pre-loaded ontology (loaded from data/subcell.txt if omitted).

    Returns
    -------
    str
        Major compartment label, e.g. 'Plasma Membrane', 'Cytosol', 'ER Lumen'.
    """
    if pd.isna(location_str) or location_str == "":
        return "Unknown"

    ontology = ontology or load_subcell_ontology()
    accessions = set(parse_location(location_str, ontology=ontology, return_accessions=True))

    if not accessions:
        return "Unknown"

    def _has(root_name: str) -> bool:
        acc = ontology.lookup(root_name)
        if not acc:
            return False
        return bool(accessions & ontology.get_descendants(acc, include_self=True))

    # Plasma/cell membrane always takes priority regardless of is_membrane flag
    if _has("cell membrane"):
        return "Plasma Membrane"

    # Cytosol/Nucleus: only return early when there is no more-specific organelle membrane
    has_cytosol_or_nucleus = _has("cytoplasm") or _has("nucleus")
    if has_cytosol_or_nucleus:
        has_organelle = is_membrane and (
            _has("mitochondrion")
            or _has("endoplasmic reticulum")
            or _has("golgi apparatus")
            or _has("peroxisome")
            or _has("lysosome")
        )
        if not has_organelle:
            return "Cytosol"

    if _has("mitochondrion"):
        return "Mitochondrial Membrane" if is_membrane else "Mitochondrial Matrix"

    if _has("endoplasmic reticulum"):
        return "ER Membrane" if is_membrane else "ER Lumen"

    if _has("golgi apparatus"):
        return "Golgi Membrane" if is_membrane else "Golgi Lumen"

    if _has("peroxisome"):
        return "Peroxisomal Membrane" if is_membrane else "Peroxisomal Matrix"

    if _has("lysosome"):
        return "Lysosomal Membrane" if is_membrane else "Lysosomal Lumen"

    # Any remaining membrane term
    if _has("membrane"):
        return "Other Membrane" if is_membrane else "Membrane-Associated Lumen"

    if is_membrane:
        return "Other Membrane"

    if has_cytosol_or_nucleus:
        return "Cytosol"

    return "Other"


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


def parse_sl_ids(sl_ids: "str | list | None") -> List[str]:
    """Parse a semicolon-joined SL accession string (as returned by the UniProt JSON API)
    into a deduplicated list of accession strings."""
    if sl_ids is None or (isinstance(sl_ids, float) and pd.isna(sl_ids)):
        return []
    if isinstance(sl_ids, list):
        return [s for s in sl_ids if s]
    return [s.strip() for s in str(sl_ids).split(";") if s.strip().startswith("SL-")]


def compartment_from_sl_ids(
    sl_ids: "str | list | None",
    is_membrane: bool = False,
    ontology: "SubcellOntology | None" = None,
) -> str:
    """
    Assign a major cellular compartment label directly from UniProt SL accession IDs.

    No text parsing — uses the SL IDs fetched from the UniProt JSON API and walks
    the SubcellOntology hierarchy to assign a compartment.

    Parameters
    ----------
    sl_ids : str or list
        Semicolon-joined SL accession string or list of accessions.
    is_membrane : bool
        Whether the protein has transmembrane domain annotation.
    ontology : SubcellOntology, optional
        Pre-loaded ontology (loaded from data/subcell.txt if omitted).

    Returns
    -------
    str
        Compartment label, e.g. 'Plasma Membrane', 'Mitochondrion', 'Cytoplasm'.
    """
    accessions = set(parse_sl_ids(sl_ids))
    if not accessions:
        return "Unknown"

    ontology = ontology or load_subcell_ontology()

    def _has(root_name: str) -> bool:
        acc = ontology.lookup(root_name)
        if not acc:
            return False
        return bool(accessions & ontology.get_descendants(acc, include_self=True))

    if _has("cell membrane"):
        return "Plasma Membrane"

    has_cytosol_or_nucleus = _has("cytoplasm") or _has("nucleus")
    if has_cytosol_or_nucleus:
        has_organelle = is_membrane and (
            _has("mitochondrion")
            or _has("endoplasmic reticulum")
            or _has("golgi apparatus")
            or _has("peroxisome")
            or _has("lysosome")
        )
        if not has_organelle:
            return "Cytoplasm" if _has("cytoplasm") else "Nucleus"

    if _has("mitochondrion"):
        return "Mitochondrion"

    if _has("endoplasmic reticulum"):
        return "Endoplasmic Reticulum"

    if _has("golgi apparatus"):
        return "Golgi Apparatus"

    if _has("peroxisome"):
        return "Peroxisome"

    if _has("lysosome"):
        return "Lysosome"

    if _has("secreted") or _has("extracellular region"):
        return "Extracellular"

    if _has("membrane"):
        return "Other Membrane"

    if is_membrane:
        return "Other Membrane"

    if has_cytosol_or_nucleus:
        return "Cytoplasm" if _has("cytoplasm") else "Nucleus"

    return "Other"
