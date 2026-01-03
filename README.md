# Gene set enrichment analysis (GSEA) 🏔️

This is the new GSEA.
Rebuilt from scratch.
It runs 1,000 times faster, delivers 99% more accurate results, and creates beautiful plots.
It’s amazing!

## Get started

```bash
gsea --help
```

Run the sarcopenia example

```bash
gsea metric-rank \
  ~/Downloads \
  ex/target.tsv \
  ex/data.tsv \
  ex/set.json \
  --standard-deviation 3 \
  --number-of-permutations 10 \
  --more-plots "WP_DNA_MISMATCH_REPAIR;WP_CELL_CYCLE"
```

Convert older file formats (.cls, .gct, and .gmt) to .tsv and .json

```bash
gsea cls ~/Downloads/1.tsv data/1.cls

gsea gct ~/Downloads/1.tsv data/1.gct

gsea gmt ~/Downloads/12.json data/h.all.v7.1.symbols.gmt data/c2.all.v7.1.symbols.gmt
```

## Install

We plan to sign this app for macOS soon.
In the meantime, [enable third-party apps](https://support.apple.com/en-us/102445#openanyway).

Download and extract the latest [release](https://github.com/GSEA-MSigDB/GSEA.jl/releases/latest).

```bash
PATH=$(pwd)/gsea/bin:$PATH
```

## Build

```bash
cd test

julia --project ../deps/build.jl app tarball
```

## Contact us

If you have any questions, issues, or concerns, feel free to [open a GitHub issue](https://github.com/GSEA-MSigDB/GSEA.jl/issues/new/choose).

---

Made by [Kata](https://github.com/KwatMDPhD/Kata.jl) ✅
