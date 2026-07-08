# %% Imports and paths
import re
import time
import warnings
from ast import literal_eval
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

ROOT     = Path(__file__).parent
RAW      = ROOT / "Human Phase separation data.xlsx"
TM_CACHE = ROOT / "data/uniprot_tm_cache.csv"
GO_CACHE = ROOT / "data/uniprot_go_cache.csv"
OUT      = ROOT / "output/full_dataset.csv"

_UNIPROT_URL   = "https://rest.uniprot.org/uniprotkb/search"
_GO_BASIC_URL  = "https://current.geneontology.org/ontology/go-basic.obo"
_GO_SLIM_URL   = "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"


# %% Load raw data
df = pd.read_excel("Human Phase separation data.xlsx", engine="openpyxl")
print(f"{len(df):,} proteins  |  p(LLPS) {df['p(LLPS)'].min():.2f}–{df['p(LLPS)'].max():.2f}")
df.head(2)


# %% Fetch GO annotations from UniProt (batched, cached)
def _fetch_uniprot(entry_ids, fields, batch_size=100, cache_path=None):
    if cache_path and Path(cache_path).exists():
        print(f"Loading cache: {cache_path}")
        return pd.read_csv(cache_path)
    n = (len(entry_ids) + batch_size - 1) // batch_size
    print(f"Fetching {len(entry_ids)} proteins in {n} batches...")
    records = []
    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i : i + batch_size]
        params = {
            "query": "accession:(" + " OR ".join(batch) + ")",
            "fields": fields,
            "format": "tsv",
            "size": batch_size,
        }
        try:
            resp = requests.get(_UNIPROT_URL, params=params, timeout=30)
            if resp.status_code == 429:
                warnings.warn("Rate limited, waiting 10s...")
                time.sleep(10)
                resp = requests.get(_UNIPROT_URL, params=params, timeout=30)
            if resp.status_code == 200:
                records.append(pd.read_csv(StringIO(resp.text), sep="\t"))
                print(f"  batch {i // batch_size + 1}/{n}: {len(records[-1])} proteins")
            else:
                warnings.warn(f"  batch {i // batch_size + 1}/{n}: HTTP {resp.status_code}")
        except requests.RequestException as e:
            warnings.warn(str(e))
        time.sleep(0.5)
    result = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)
    return result


go_raw = _fetch_uniprot(
    df["Entry"].tolist(),
    fields="accession,go_id,cc_subcellular_location",
    cache_path=GO_CACHE,
)

# Normalise column names
col_map = {}
for col in go_raw.columns:
    low = col.lower()
    if low == "entry":
        col_map[col] = "Entry"
    elif "gene ontology ids" in low or low == "go ids":
        col_map[col] = "GO_IDs"
    elif "subcellular location" in low:
        col_map[col] = "Subcellular location [CC]"
go_raw = go_raw.rename(columns=col_map)
for c in ("GO_IDs", "Subcellular location [CC]"):
    if c not in go_raw.columns:
        go_raw[c] = np.nan

df = df.merge(go_raw[["Entry", "GO_IDs", "Subcellular location [CC]"]], on="Entry", how="left")

# Collapse duplicate location columns if they were already in the raw xlsx
if "Subcellular location [CC]_x" in df.columns:
    df["Subcellular location [CC]"] = df["Subcellular location [CC]_x"].fillna(
        df["Subcellular location [CC]_y"]
    )
    df = df.drop(columns=["Subcellular location [CC]_x", "Subcellular location [CC]_y"])

df[["Entry", "GO_IDs"]].head(2)


# %% pLLPS classification (Low / Medium / High)
df["pLLPS_Class"] = pd.cut(
    df["p(LLPS)"],
    bins=[-float("inf"), 0.4, 0.7, float("inf")],
    labels=["Low", "Medium", "High"],
).astype(str)

df["pLLPS_Class"].value_counts()


# %% GO helpers: parse IDs, load DAGs, map to GO slim
def parse_go_ids(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        return sorted({str(g) for g in value if g})
    return sorted(set(re.findall(r"GO:\d{7}", str(value))))


_GO_DAG_CACHE      = {}
_GO_SLIM_DAG_CACHE = {}
_GO_SLIM_MAP_CACHE = {}
_GO_DESC_CACHE     = {}


def _download_obo(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(requests.get(url, timeout=60).content)
    return path


def load_go_dag(path=None):
    path = Path(path or ROOT / "data/go/go-basic.obo")
    if not path.exists():
        _download_obo(_GO_BASIC_URL, path)
    key = str(path.resolve())
    if key not in _GO_DAG_CACHE:
        from goatools.obo_parser import GODag
        _GO_DAG_CACHE[key] = GODag(str(path))
    return _GO_DAG_CACHE[key]


def load_go_slim_dag(path=None):
    path = Path(path or ROOT / "data/go/goslim_generic.obo")
    if not path.exists():
        _download_obo(_GO_SLIM_URL, path)
    key = str(path.resolve())
    if key not in _GO_SLIM_DAG_CACHE:
        from goatools.obo_parser import GODag
        _GO_SLIM_DAG_CACHE[key] = GODag(str(path))
    return _GO_SLIM_DAG_CACHE[key]


def map_to_slim(go_ids_parsed, go_dag, slim_dag):
    import goatools.mapslim as mapslim
    slim_ids = set()
    for gid in go_ids_parsed:
        if gid not in go_dag:
            continue
        if gid not in _GO_SLIM_MAP_CACHE:
            direct, _ = mapslim.mapslim(gid, go_dag, slim_dag)
            _GO_SLIM_MAP_CACHE[gid] = direct
        slim_ids |= _GO_SLIM_MAP_CACHE[gid]
    return sorted(slim_ids)


go_dag   = load_go_dag()
slim_dag = load_go_slim_dag()


# %% Split GO IDs by namespace and compute GO slim columns
def split_go_namespaces(value):
    bp, mf, cc = [], [], []
    for gid in parse_go_ids(value):
        if gid not in go_dag:
            continue
        ns = getattr(go_dag[gid], "namespace", None)
        if ns == "biological_process":
            bp.append(gid)
        elif ns == "molecular_function":
            mf.append(gid)
        elif ns == "cellular_component":
            cc.append(gid)
    return bp, mf, cc


df[["GO_BP", "GO_MF", "GO_CC"]] = df["GO_IDs"].apply(
    lambda v: pd.Series(split_go_namespaces(v))
)


def slim_names_for(go_list):
    ids = map_to_slim(go_list, go_dag, slim_dag)
    return [slim_dag[sid].name for sid in ids if sid in slim_dag]


df["Function_Slim"] = df["GO_BP"].apply(slim_names_for) + df["GO_MF"].apply(slim_names_for)
df["Location_Slim"] = df["GO_CC"].apply(slim_names_for)
df["Function_Top"]  = df["Function_Slim"].apply(lambda l: l[0] if l else "Other")
df["Location_Top"]  = df["Location_Slim"].apply(lambda l: l[0] if l else "Other")

df[["Entry", "Function_Top", "Location_Top"]].head(5)


# %% Functional categories (YAML-driven, GO-based)
def _expand_descendants(go_ids, go_dag):
    out = set()
    for gid in go_ids:
        if gid in _GO_DESC_CACHE:
            out |= _GO_DESC_CACHE[gid]
        elif gid in go_dag:
            desc = set(go_dag[gid].get_all_children()) | {gid}
            _GO_DESC_CACHE[gid] = desc
            out |= desc
    return out


with open(ROOT / "data/functional_classification_terms.yaml") as f:
    _func_terms = yaml.safe_load(f)

_group_descendants = {
    cat: _expand_descendants(data.get("go_ids", []), go_dag)
    for cat, data in _func_terms.get("functional_groups", {}).items()
}


def classify_function(go_value):
    ids = set(parse_go_ids(go_value))
    return [cat for cat, desc in _group_descendants.items() if ids & desc]


df["Functional_Categories"] = df["GO_IDs"].apply(classify_function)
df["Functional_Categories"].apply(len).describe()


# %% GO slim categories (flat list, all namespaces)
df["GO_Slim_Categories"] = df["GO_IDs"].apply(
    lambda v: slim_names_for(parse_go_ids(v))
)


# %% TMD count from local UniProt TM cache
def count_tm_spans(domain_str):
    if pd.isna(domain_str) or not domain_str:
        return 0
    return len(re.findall(r"\d+\.\.\d+", str(domain_str)))


if TM_CACHE.exists():
    tm_raw = _fetch_uniprot(
        df["Entry"].tolist(),
        fields="accession,ft_transmem,ft_intramem",
        cache_path=TM_CACHE,
    )
    col_map = {}
    for col in tm_raw.columns:
        low = col.lower()
        if low == "entry":
            col_map[col] = "Entry"
        elif "intramembrane" in low:
            col_map[col] = "Intramembrane"
        elif "transmembrane" in low:
            col_map[col] = "Transmembrane"
    tm_raw = tm_raw.rename(columns=col_map)
    df = df.merge(tm_raw[["Entry", "Transmembrane", "Intramembrane"]], on="Entry", how="left")
    df["TMD_count"] = (
        df["Transmembrane"].apply(count_tm_spans) +
        df["Intramembrane"].apply(count_tm_spans)
    )
else:
    df["TMD_count"] = 0
    print("TM cache not found — TMD_count set to 0. Run scripts/analysis/enrich_dataset_with_tmd.py to populate.")

df["TMD_count"].value_counts().sort_index().head(10)


# %% Derive Is_Membrane and Compartment; drop legacy columns
df["Is_Membrane"] = df["TMD_count"].fillna(0).astype(int) > 0
df["Compartment"]  = df["Location_Top"]

_drop = [
    "All_Functional_Groups", "Function Categories", "Functional Group", "Functional Slim",
    "Transmembrane", "Intramembrane",
    "Subcellular location [CC]_x", "Subcellular location [CC]_y",
    "Location Categories", "Location_IDs", "Location Primary",
    "Location_SL_IDs", "Location_Names",
    *[c for c in df.columns if c.startswith("Is_") and c != "Is_Membrane"],
]
df = df.drop(columns=[c for c in _drop if c in df.columns])

print(f"Is_Membrane: {df['Is_Membrane'].sum():,} proteins")
df["Compartment"].value_counts()


# %% Save
OUT.parent.mkdir(exist_ok=True)
df.to_csv(OUT, index=False)
print(f"Saved {len(df):,} proteins  |  {df.shape[1]} columns → {OUT}")
print(df.columns.tolist())


# %% Export membrane protein UniProt ID lists
_mem = df[df["Is_Membrane"]]

_out_all = OUT.parent / "membrane_proteins_all.csv"
_mem[["Entry"]].to_csv(_out_all, index=False)
print(f"Membrane proteins (all):       {len(_mem):,} → {_out_all}")

_out_high = OUT.parent / "membrane_proteins_pLLPS_0.6plus.csv"
_mem_high = _mem[_mem["p(LLPS)"] >= 0.6]
_mem_high[["Entry"]].to_csv(_out_high, index=False)
print(f"Membrane proteins (p≥0.6):    {len(_mem_high):,} → {_out_high}")
