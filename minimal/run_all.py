"""
Run the whole pipeline in order.

    python run_all.py

Stops at the first step that fails, and tells you which one. Takes about a minute.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

STEPS = [
    "01_download.py",
    "02_build_master.py",
    "03_explore.py",
    "04_stats.py",              # topology / hydropathy
    "05_figures.py",            # figures for 04
    "06_background.py",         # the offset and the background flip
    "07_thresholds.py",         # published cutoffs and the shortlist
    "08_figures_background.py"  # figures for 06 and 07
]

for step in STEPS:
    print("\n" + "=" * 70)
    print(f"  {step}")
    print("=" * 70)
    r = subprocess.run([sys.executable, str(HERE / step)], cwd=HERE)
    if r.returncode != 0:
        sys.exit(f"\n{step} failed -- fix it before running the later steps.")

print("\n" + "=" * 70)
print("  done")
print("=" * 70)
