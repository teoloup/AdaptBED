"""
Pipeline steps:
- Accept BED3/6/12 input and normalize to BED6 for processing (name if exists; score=0; strand as specified).
- Parse, sort, and merge (not strand-aware) to produce the ROI BED (bed6) with strand set to '.'; keep names if present (distinct names collapsed).
- Compute stats on ROI: number of regions, total size (bp), mean and median target size, and percent genome covered
	using bioframe.fetch_chromsizes(genome) with default hg38 if not provided.
	- Compute a GLOBAL buffer cap unless --ignore-cap so that TOTAL target size ≤ 3× TOTAL ROI size:
	  Let N be ROI interval count and S the total ROI size (bp). The uncapped target size is S + 2*N*buffer.
	  Enforce S + 2*N*buffer ≤ 3*S ⇒ buffer ≤ floor(S/N). Apply this once globally.
- Build adaptive bed by, for each ROI interval, creating two intervals (no original):
  - chrom  (start - buffer)  end  name  0  +
  - chrom  start  (end + buffer)  name  0  -
  Clamp starts to >= 0. Strand is as shown.
- Sort and merge the adaptive bed strand-aware (-s) and keep strand (-c 6 -o distinct). Output as BED6 with '.' name and 0 score.
- Compute stats again on adaptive bed.

Notes:
- Starts are clamped to 0 and ends are clamped to chromosome sizes from bioframe.fetch_chromsizes.
"""

from __future__ import annotations
from dataclasses import dataclass
import math
from pathlib import Path
from statistics import mean, median
from typing import Iterable, List, Tuple, Dict, Optional
import bioframe as bf
import pybedtools  
from pybedtools import BedTool  
from statistics import NormalDist


# ---------- Small helpers ----------

def _interval_len(fields: List[str]) -> int:
	return int(fields[2]) - int(fields[1])


def _to_bed6_lines(bt: BedTool) -> List[str]:
	"""Convert any BED3/6/12 BedTool to bed6 lines.

	- name: use column 4 if present, else '.'
	- score: '0'
	- strand: '.' for ROI stage (we'll set explicit strands later for adaptive stage)
	"""

	lines: List[str] = []
	for f in bt:
		fields = list(f.fields)
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		name = fields[3] if len(fields) >= 4 and fields[3] else "."
		score = "0"
		strand = "."
		lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")
	return lines


def _bedtool_from_lines(lines: List[str]) -> BedTool:
	# Construct a BedTool from in-memory lines
	data = "\n".join(lines) + ("\n" if lines and not lines[-1].endswith("\n") else "")
	return BedTool(data, from_string=True)


def _merge_roi_bed6(bt6: BedTool) -> BedTool:
	"""Sort and merge (not strand-aware) returning bed6 with collapsed distinct names.

	Output columns: chrom, start, end, name, score=0, strand='.'.
	"""

	# Always sort before merging for correctness
	bt_sorted = bt6.sort()

	# Merge overlapping intervals; collapse distinct names if present.
	# Input is bed6, so column 4 holds names. After merge with c=4, o=distinct, the output
	# has bed3 + a 4th column containing the distinct names.
	merged = bt_sorted.merge(c=[4], o=["distinct"])  # type: ignore[arg-type]

	# Reconstruct bed6 lines with name from merged 4th column, score=0, strand='.'
	lines: List[str] = []
	for f in merged:
		mfields = list(f.fields)
		chrom, start, end = mfields[0], int(mfields[1]), int(mfields[2])
		name = mfields[3] if len(mfields) >= 4 and mfields[3] else "."
		lines.append(f"{chrom}\t{start}\t{end}\t{name}\t0\t.")
	return _bedtool_from_lines(lines)


def _compute_stats(bt: BedTool, genome: str) -> Dict[str, float]:
	"""Compute stats for a BedTool (assumes at least bed3 columns present).

	Returns keys: count, total_bp, mean_bp, median_bp, genome_size_bp, percent_genome
	"""

	lens = [int(f.end) - int(f.start) for f in bt]
	count = len(lens)
	total_bp = int(sum(lens)) if lens else 0
	mean_bp = float(mean(lens)) if lens else 0.0
	median_bp = float(median(lens)) if lens else 0.0

	chromsizes = bf.fetch_chromsizes(genome)
	# chromsizes is a pandas Series; use sum() for total genome size
	try:
		genome_size_bp = int(chromsizes.sum())
	except Exception:
		# Fallback: iterate values if sum() fails for any reason
		vals = getattr(chromsizes, "values", [])
		genome_size_bp = int(sum(int(v) for v in vals))
	percent_genome = (total_bp / genome_size_bp * 100.0) if genome_size_bp else 0.0

	return {
		"count": float(count),
		"total_bp": float(total_bp),
		"mean_bp": mean_bp,
		"median_bp": median_bp,
		"genome_size_bp": float(genome_size_bp),
		"percent_genome": percent_genome,
	}


# ---------- Buffer calculation (log-normal) ----------

def nxx_from_mean_cv(mean_bp: float, cv_percent: float, xx: float) -> float:
	"""
	Base-weighted Nxx from mean & CV assuming log-normal read lengths.
	Nxx = minimum length L such that reads >= L contain xx% of bases.
	"""
	c = cv_percent / 100.0
	a = math.log(1.0 + c * c)  # a = ln(1 + CV^2)
	sigma = math.sqrt(a)
	p = 1.0 - (xx / 100.0)  # p-quantile of length-biased log-normal
	z = NormalDist().inv_cdf(p)
	return mean_bp * math.sqrt(1.0 + c * c) * math.exp(sigma * z)


def compute_buffer_from_mean_cv(mean_bp: float, cv_percent: float) -> int:
	"""Return the suggested buffer size using N15 (rounded to nearest int)."""
	n15 = nxx_from_mean_cv(mean_bp, cv_percent, 15.0)
	# round to nearest integer for bp
	return int(round(n15))


def _cap_buffer_globally(roi_bt6: BedTool, requested_buffer: int, ignore_cap: bool) -> Tuple[int, bool]:
	"""Cap buffer size globally so that TOTAL target size ≤ 3× TOTAL ROI size.

	Given merged ROI bed6, let N = number of intervals and S = sum of their lengths.
	Target size (without double-counting ROI sequence) is S + 2 * N * buffer.
	Cap rule: S + 2*N*buffer ≤ 3*S ⇒ buffer ≤ floor(S / N) when N > 0.

	Returns (buffer_to_use, capped_applied?). If ignore_cap or N == 0, return requested as-is.
	"""
	if ignore_cap:
		return requested_buffer, False

	lens = [int(f.end) - int(f.start) for f in roi_bt6]
	n = len(lens)
	if n == 0:
		return requested_buffer, False
	s = int(sum(lens))
	if s <= 0:
		return requested_buffer, False
	# Max buffer that keeps total target size within 3x ROI size
	max_buffer = s // n
	if requested_buffer > max_buffer:
		return max_buffer, True
	return requested_buffer, False


def _build_adaptive_from_roi(
	bt_roi6: BedTool,
	buffer_size: int,
	chromsizes: Dict[str, int],
	strand_aware: bool,
) -> Tuple[BedTool, List[str]]:
	"""From ROI bed6, build adaptive bed6 by expanding regions.

	Returns the final bed6 and a list of clamp warnings.
	"""

	out_lines: List[str] = []
	clamp_warnings: List[str] = []

	for f in bt_roi6:
		fields = list(f.fields)
		chrom, a, b = fields[0], int(fields[1]), int(fields[2])
		name = fields[3] if len(fields) >= 4 and fields[3] else "."
		chrom_size = chromsizes.get(chrom)

		if strand_aware:
			# Two entries: left (+) and right (-)
			left_start = max(0, a - buffer_size)
			left_end = b
			right_start = a
			right_end = b + buffer_size
			if chrom_size is not None:
				if left_end > int(chrom_size):
					clamp_warnings.append(
						f"Clamped to chrom end: {chrom}:{left_start}-{left_end} -> {chrom}:{left_start}-{int(chrom_size)}"
					)
				if right_end > int(chrom_size):
					clamp_warnings.append(
						f"Clamped to chrom end: {chrom}:{right_start}-{right_end} -> {chrom}:{right_start}-{int(chrom_size)}"
					)
				left_end = min(left_end, int(chrom_size))
				right_end = min(right_end, int(chrom_size))
			if left_start < left_end:
				out_lines.append(f"{chrom}\t{left_start}\t{left_end}\t{name}\t0\t+")
			if right_start < right_end:
				out_lines.append(f"{chrom}\t{right_start}\t{right_end}\t{name}\t0\t-")
		else:
			# Single entry expanded both sides; strand '.'; keep name
			start_exp = max(0, a - buffer_size)
			end_exp = b + buffer_size
			if chrom_size is not None:
				if end_exp > int(chrom_size):
					clamp_warnings.append(
						f"Clamped to chrom end: {chrom}:{start_exp}-{end_exp} -> {chrom}:{start_exp}-{int(chrom_size)}"
					)
				end_exp = min(end_exp, int(chrom_size))
			if start_exp < end_exp:
				out_lines.append(f"{chrom}\t{start_exp}\t{end_exp}\t{name}\t0\t.")

	bt_adaptive = _bedtool_from_lines(out_lines)
	bt_adaptive_sorted = bt_adaptive.sort()

	if strand_aware:
		# Merge strand-aware and retain distinct names (col 4) and strand (col 6)
		merged = bt_adaptive_sorted.merge(s=True, c=[4, 6], o=["distinct", "distinct"])  # type: ignore[arg-type]
		# Reconstruct proper bed6 with distinct names and retained strand
		lines: List[str] = []
		for f in merged:
			mfields = list(f.fields)
			chrom, start, end = mfields[0], int(mfields[1]), int(mfields[2])
			name = mfields[3] if len(mfields) >= 4 and mfields[3] else "."
			strand = mfields[4] if len(mfields) >= 5 and mfields[4] else "."
			lines.append(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}")
		return _bedtool_from_lines(lines), clamp_warnings
	else:
		# Not strand-aware: simple merge; collapse distinct names
		merged = bt_adaptive_sorted.merge(c=[4], o=["distinct"])  # type: ignore[arg-type]
		lines: List[str] = []
		for f in merged:
			mfields = list(f.fields)
			chrom, start, end = mfields[0], int(mfields[1]), int(mfields[2])
			name = mfields[3] if len(mfields) >= 4 and mfields[3] else "."
			lines.append(f"{chrom}\t{start}\t{end}\t{name}\t0\t.")
		return _bedtool_from_lines(lines), clamp_warnings


# ---------- Public API ----------

@dataclass
class PipelineOutputs:
	roi_bed6: BedTool
	roi_stats: Dict[str, float]
	adaptive_bed6: BedTool
	adaptive_stats: Dict[str, float]
	buffer_size_used: int
	global_cap_applied: bool
	clamp_warnings: List[str]


def run_pipeline(
	input_bed_path: Path,
	genome: str,
	buffer_size: int,
	ignore_cap: bool,
	strand_aware_buffer: bool,
) -> PipelineOutputs:
	"""Execute the full pipeline and return BedTools and stats.

	- Normalize to bed6 with name if present, score=0, strand='.'
	- Sort and merge (not strand-aware) -> ROI bed6
	- Stats on ROI
	- Build adaptive bed6 with buffer cap per ROI unless ignore_cap
	- Sort/merge strand-aware and stats
	"""

	bt_raw = BedTool(str(input_bed_path))
	bt6 = _bedtool_from_lines(_to_bed6_lines(bt_raw))
	roi_bed6 = _merge_roi_bed6(bt6)
	roi_stats = _compute_stats(roi_bed6, genome=genome)

	# Apply global cap based on ROI stats
	buffer_to_use, cap_applied = _cap_buffer_globally(roi_bed6, buffer_size, ignore_cap)

	# Fetch chromsizes once for clamping and pass to builder
	chromsizes_series = bf.fetch_chromsizes(genome)
	# Normalize values to int for safety (Series -> dict)
	chromsizes_int: Dict[str, int] = {str(k): int(v) for k, v in chromsizes_series.items()}
	adaptive_bed6, clamp_warnings = _build_adaptive_from_roi(
		roi_bed6, buffer_to_use, chromsizes_int, strand_aware_buffer
	)
	# Compute stats without double-counting ROI bases: use unstranded union if strand-aware
	if strand_aware_buffer:
		# Convert to bed3 lines and merge to union coverage
		bed3_lines: List[str] = [f"{f.chrom}\t{int(f.start)}\t{int(f.end)}" for f in adaptive_bed6]
		union_bt = _bedtool_from_lines(bed3_lines).sort().merge()
		adaptive_stats = _compute_stats(union_bt, genome=genome)
	else:
		adaptive_stats = _compute_stats(adaptive_bed6, genome=genome)

	return PipelineOutputs(
		roi_bed6=roi_bed6,
		roi_stats=roi_stats,
		adaptive_bed6=adaptive_bed6,
		adaptive_stats=adaptive_stats,
		buffer_size_used=buffer_to_use,
		global_cap_applied=cap_applied,
		clamp_warnings=clamp_warnings,
	)

