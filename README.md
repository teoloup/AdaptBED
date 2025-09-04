# AdaptBED — adaptive-sampling BED formatter

Adaptive sampling BED helper built on pybedtools (0.12.0) with chrom sizes from bioframe.

## Features

- Accepts BED3/6/12; normalizes to BED6.
- Sorts and merges ROIs (not strand-aware) to produce a clean target set.
- Computes stats for ROI and final targets (count, total bp, mean/median bp, % genome).
- Expands targets by a buffer on each side.
- Global cap to keep total target size ≤ 3× total ROI size (S + 2·N·buffer ≤ 3·S).
- Optional strand-aware buffering (+ left, - right) with strand-aware merge.
- Clamps to chromosome bounds using bioframe chromsizes.
- Optional buffer auto-calculation from mean and CV using a log-normal model (Nxx); prints and records N50/N25/N15/N10/N01 and picks N15 by default.
- Final adaptive bed keeps distinct names in column 4.
- Based on guidelines from ONT Adaptive sampling guide (https://nanoporetech.com/document/adaptive-sampling#where-to-find-create-and-modify-fasta-and-bed-files)

## Install

Requires Python 3.8+.

```powershell
# From the project root
pip install -U pip
pip install -e .
```

This will install dependencies:
- pybedtools==0.12.0
- bioframe>=0.6

## Usage

```
# Basic
adaptbed -i input.bed -o outdir --prefix sample

# With explicit buffer-size (bp)
adaptbed -i input.bed -o outdir --prefix sample -b 8000

# Strand-aware buffering (+ left, - right; merged per strand)
adaptbed -i input.bed -o outdir --prefix sample -s

# Auto-calculate buffer from read mean and CV (percent)
# Prints and records N50/N25/N15/N10/N01; uses N15 as the buffer
adaptbed -i input.bed -o outdir --prefix sample --calculate-buffer --mean-length 12000 --cv 30

# Ignore the global 3× cap (may greatly increase targets)
adaptbed -i input.bed -o outdir --prefix sample --ignore-cap
```

Outputs in `outdir`:
- `<prefix>.roi.bed`: merged ROI BED6
- `<prefix>.adaptive_sampling.bed`: final targets BED6 (distinct names in col 4; strand-aware if `-s`)
- `<prefix>.stats.tsv`: summary + meta (buffer used; when `--calculate-buffer`, mean, cv, Nxx values and choice)

## How buffer and cap work

Given merged ROI size S and count N, requested buffer B:
- Target size T(B) = S + 2·N·B (no double-count of ROI bases)
- Cap rule (unless `--ignore-cap`): S + 2·N·B ≤ 3·S ⇒ B ≤ floor(S/N)
- We use `B_used = min(B, floor(S/N))` and then clamp to chromosome ends.

## Notes

- Genome default is `hg38`. Use `--genome` to select others supported by bioframe.
- Stats for strand-aware targets are computed on the unstranded union to avoid double-counting.
- If any intervals touch chromosome ends, clamping warnings are printed.

## Troubleshooting

- Missing bioframe: `pip install bioframe`
- pybedtools import errors: ensure Python 3.8+ and reinstall with `pip install -e .`

## License

MIT

