"""
bin_fasta_for_submission.py

Splits the topology region FASTA files (output/topology_fasta/) into fixed-size
batches for web servers with a per-submission sequence limit — e.g. PSPHunter
accepts 5 sequences at a time.

Usage: python bin_fasta_for_submission.py --tool PSPHunter --batch-size 5

Outputs
-------
  output/topology_fasta/{tool}_batches/{region}_batch{NN}.fasta
  output/topology_fasta/{tool}_batches_manifest.csv
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT      = Path(__file__).parent
FASTA_DIR = ROOT / "output" / "topology_fasta"
REGIONS   = ["Cytoplasmic", "Transmembrane", "Extracellular_Lumenal", "Whole"]


def read_fasta(path):
    records = []
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", required=True, help="Web server name, used for the output folder")
    parser.add_argument("--batch-size", type=int, required=True)
    args = parser.parse_args()

    out_dir = FASTA_DIR / f"{args.tool}_batches"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for region in REGIONS:
        fasta_path = FASTA_DIR / f"{region}.fasta"
        if not fasta_path.exists():
            print(f"[skip] {fasta_path} not found")
            continue
        records = read_fasta(fasta_path)
        n_batches = -(-len(records) // args.batch_size)  # ceil division
        for i in range(n_batches):
            batch = records[i * args.batch_size:(i + 1) * args.batch_size]
            batch_path = out_dir / f"{region}_batch{i + 1:02d}.fasta"
            with open(batch_path, "w") as f:
                for header, seq in batch:
                    f.write(f">{header}\n{seq}\n")
            manifest_rows.append({
                "Region": region,
                "Batch file": str(batch_path.relative_to(ROOT)),
                "UniProt_IDs": ";".join(h for h, _ in batch),
                "n_seqs": len(batch),
                "Scores (paste here)": "",
            })
        print(f"{region}: {len(records)} sequences -> {n_batches} batches of <= {args.batch_size}")

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = FASTA_DIR / f"{args.tool}_batches_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"\nSaved {manifest_path.relative_to(ROOT)}  ({len(manifest)} batches total)")


if __name__ == "__main__":
    main()
