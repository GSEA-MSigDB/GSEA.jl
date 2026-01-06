# Gene set enrichment analysis (GSEA) 🏔️

This is the new GSEA.

We reimplemented the `S0` and `S0a` algorithms.  
They run 1,000 times faster, reproduce all results, and create prettier plots.

We also implemented a new `D2` algorithm.  
It uses information theory to deliver the most accurate, interpretable, and robust gene-set scores.

## Install

Download the latest [release](https://github.com/GSEA-MSigDB/GSEA.jl/releases/latest), extract it, and find the executable at `gsea/bin/gsea`.

Or clone this repository and build

```bash
git clone https://github.com/GSEA-MSigDB/GSEA.jl

cd GSEA.jl

julia --project deps/build.jl app tarball

# Find the executable at `build/gsea/bin/gsea`.
```

## Get started

```bash
gsea --help
```

![GSEA command-line interface screenshot](screenshot.png)

Run the sarcopenia example

```bash
gsea metric-rank \
    ~/Downloads \
    in/ex.target.tsv \
    in/ex.data.tsv \
    in/ex.set.json \
    --number-of-permutations 10 \
    --more-plots "WP_DNA_MISMATCH_REPAIR;WP_CELL_CYCLE"
```

## Contact us

If you have any questions, issues, or concerns, feel free to [open a GitHub issue](https://github.com/GSEA-MSigDB/GSEA.jl/issues/new/choose).

Every report helps make GSEA better.

---

Made by [Kata](https://github.com/KwatMDPhD/Kata.jl) ✅
