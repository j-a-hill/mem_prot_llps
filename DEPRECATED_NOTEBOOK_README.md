# Notebook Deprecation Notice

## ⚠️ DEPRECATED NOTEBOOK

**This notebook (`pllps_interaction_analysis_DEPRECATED.ipynb`) has been deprecated due to its large size (18MB) which causes crashes and makes it difficult to use.**

## Why was it deprecated?

- **Too large:** 18MB file size causes Jupyter to crash
- **Embedded outputs:** All plots and data were embedded in the notebook
- **Hard to navigate:** Over 50 cells in one file
- **Slow loading:** Takes minutes to load
- **Memory intensive:** Often crashes when running

## What should you use instead?

The analysis has been **split into 6 focused notebooks** that are fast, reliable, and easy to use:

1. **`01_data_loading_and_classification.ipynb`** (8.5 KB)
2. **`02_string_interactions.ipynb`** (12 KB)
3. **`03_enrichment_analysis.ipynb`** (11 KB)
4. **`04_network_analysis.ipynb`** (13 KB)
5. **`05_functional_groups.ipynb`** (12 KB)
6. **`06_visualization_summary.ipynb`** (16 KB)

**Total size: 72 KB** (250x smaller!)

## Benefits of the new notebooks

✅ **Small file sizes** - No more crashes!  
✅ **Results saved to disk** - Reusable outputs in `results/` directory  
✅ **Faster loading** - Loads instantly  
✅ **Modular** - Run only what you need  
✅ **Clear workflow** - Each notebook has a focused purpose  

## How to migrate

See **[`NOTEBOOK_GUIDE.md`](NOTEBOOK_GUIDE.md)** for:
- Detailed guide to each notebook
- How to run the analysis step-by-step
- How to load and reuse saved results
- Comparison with old notebook structure

## Can I still use this old notebook?

You can, but we **strongly recommend** using the new modular notebooks instead. This old notebook:
- May crash due to large size
- Takes a long time to load
- Re-computes everything (slow)
- Does not save results for reuse

The new notebooks are faster, more reliable, and easier to use.

---

**Questions?** See [`NOTEBOOK_GUIDE.md`](NOTEBOOK_GUIDE.md) or open an issue on GitHub.
