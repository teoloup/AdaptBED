"""
BED processing utilities implemented exactly as specified.

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
- Uses bioframe for BED operations (merge, overlap, subtract) for reliability and cross-platform compatibility.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from statistics import mean, median
from typing import Iterable, List, Tuple, Dict, Optional
import bioframe as bf
import pandas as pd
from statistics import NormalDist
from pybedtools import BedTool  # type: ignore

import os
import gzip
import tempfile
import requests


# ---------- Bioframe-based BED operations (reliable and cross-platform) ----------

def _bedframe_from_bedtool(bt):
    """Convert BedTool to bioframe-compatible DataFrame."""
    rows = []
    for feature in bt:
        rows.append({
            'chrom': feature.chrom,
            'start': int(feature.start),
            'end': int(feature.end),
            'name': feature.name if hasattr(feature, 'name') and feature.name != '.' else f'interval_{len(rows)}',
            'score': getattr(feature, 'score', 0),
            'strand': getattr(feature, 'strand', '.')
        })
    return pd.DataFrame(rows)


def _bedtool_from_bedframe(df: pd.DataFrame):
    """Convert bioframe DataFrame back to BedTool format."""
    lines = []
    for _, row in df.iterrows():
        # Handle missing columns with defaults
        name = row.get('name', '.')
        score = row.get('score', 0)
        strand = row.get('strand', '.')
        
        line = f"{row['chrom']}\t{int(row['start'])}\t{int(row['end'])}\t{name}\t{score}\t{strand}"
        lines.append(line)
    
    data = "\n".join(lines)
    return BedTool(data, from_string=True)


def _get_cached_chromsizes(genome: str, custom_cache_dir: Optional[Path] = None):
    """Get chromosome sizes, using cache when available."""
    cache_dir = _get_cache_dir(custom_cache_dir)
    chromsizes_file = cache_dir / f"{genome}_chromsizes.tsv"

    # Try to load from cache first
    if chromsizes_file.exists():
        try:
            chromsizes_df = pd.read_csv(chromsizes_file, sep='\t', index_col=0, header=None, names=['size'])
            chromsizes_series = chromsizes_df['size']
            print(f"Using cached chromosome sizes for {genome}")
            return chromsizes_series
        except Exception as e:
            print(f"Warning: Could not load cached chromsizes ({e}), fetching fresh data...")

    # Fetch from bioframe and cache
    print(f"Fetching chromosome sizes for {genome} from bioframe...")
    chromsizes_series = bf.fetch_chromsizes(genome)

    # Cache the result
    try:
        chromsizes_df = chromsizes_series.to_frame()
        chromsizes_df.to_csv(chromsizes_file, sep='\t', header=False)
        print(f"Cached chromosome sizes to {chromsizes_file}")
    except Exception as e:
        print(f"Warning: Could not cache chromsizes ({e})")

    return chromsizes_series

def _get_cache_dir(custom_cache_dir: Optional[Path] = None) -> Path:
    """Get the cache directory for storing genomic data."""
    if custom_cache_dir is not None:
        cache_dir = Path(custom_cache_dir)
    else:
        # Use cache directory within the project folder
        # More robust approach: find the script's directory and go up to project root
        script_dir = Path(__file__).parent  # src/bed_buffer_tool
        project_root = script_dir.parent.parent  # Go up to project root
        cache_dir = project_root / ".adaptbed_cache"

    cache_dir.mkdir(exist_ok=True, parents=True)
    return cache_dir


def _download_repeatmasker(genome: str, custom_cache_dir: Optional[Path] = None) -> Path:
    """Download RepeatMasker data from UCSC and cache it locally."""
    cache_dir = _get_cache_dir(custom_cache_dir)
    genome_dir = cache_dir / genome
    genome_dir.mkdir(exist_ok=True)
    rmsk_file = genome_dir / "rmsk.bed.gz"

    if rmsk_file.exists():
        print(f"Using cached RepeatMasker data: {rmsk_file}")
        return rmsk_file

    # Download from UCSC with progress reporting
    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz"
    print(f"Downloading RepeatMasker data for {genome} from UCSC...")

    # Use streaming download to show progress
    response = requests.get(url, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))
    downloaded_size = 0

    # Save compressed with progress reporting
    with open(rmsk_file, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded_size += len(chunk)

                # Progress reporting
                if total_size > 0:
                    progress = (downloaded_size / total_size) * 100
                    print(f"\rDownload progress: {progress:.1f}%", end='', flush=True)
                else:
                    # If no content-length header, show bytes downloaded
                    if downloaded_size % (1024 * 1024) == 0:  # Every MB
                        mb_downloaded = downloaded_size / (1024 * 1024)
                        print(f"\rDownloaded: {mb_downloaded:.1f} MB", end='', flush=True)

    print("\nDownload complete!")
    return rmsk_file


def _load_repeatmasker_bed(genome: str, min_size: int = 400, custom_cache_dir: Optional[Path] = None) -> BedTool:
    """Load RepeatMasker data as BedTool, filtering by minimum size."""
    rmsk_file = _download_repeatmasker(genome, custom_cache_dir)

    print(f"Loading repeat data from {rmsk_file}...")

    # Read the gzipped TSV file with progress reporting
    with gzip.open(rmsk_file, 'rt') as f:
        # Skip header if present, parse TSV
        lines = []
        processed_count = 0
        filtered_count = 0

        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            processed_count += 1

            # Progress reporting every 100k lines
            if processed_count % 100000 == 0:
                print(f"Processed {processed_count:,} repeat annotations, filtered {filtered_count:,}...")

            chrom, start, end = fields[5], fields[6], fields[7]  # Correct field indices
            size = int(end) - int(start)

            if size >= min_size:
                filtered_count += 1
                # Convert to BED6: chrom, start, end, name, score, strand
                name = fields[10] if len(fields) > 10 else "repeat"  # repName
                score = fields[1] if len(fields) > 1 else "0"  # swScore
                strand = fields[8] if len(fields) > 8 else "."  # strand
                lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")

        print(f"Completed: processed {processed_count:,} annotations, kept {filtered_count:,} repeats >= {min_size}bp")

    # Create BedTool from lines
    if lines:
        print(f"Creating BedTool with {len(lines):,} repeat regions...")
        return _bedtool_from_lines(lines)
    else:
        print("No repeat regions found matching criteria")
        # Return empty BedTool
        return _bedtool_from_lines([])


def _mask_repeats_in_bed(adaptive_bed: BedTool, repeats_bed: BedTool) -> Tuple[BedTool, int]:
    """Mask repeats in adaptive BED by subtracting overlaps. Returns masked BED and bases removed."""
    if len(repeats_bed) == 0:
        return adaptive_bed, 0

    print(f"Starting repeat masking: {len(adaptive_bed)} adaptive regions, {len(repeats_bed)} repeat regions")

    # Convert to bioframe DataFrames for operations
    print("Converting to bioframe DataFrames...")
    adaptive_df = _bedframe_from_bedtool(adaptive_bed)
    repeats_df = _bedframe_from_bedtool(repeats_bed)

    # Get original size
    original_size = int((adaptive_df['end'] - adaptive_df['start']).sum())
    print(f"Original adaptive BED size: {original_size:,} bp")

    # Use bioframe subtract
    print("Performing bioframe subtract operation...")
    masked_df = bf.subtract(adaptive_df, repeats_df)

    # Calculate bases removed
    masked_size = int((masked_df['end'] - masked_df['start']).sum())
    bases_removed = original_size - masked_size

    print(f"Repeat masking completed: removed {bases_removed:,} bp, final size: {masked_size:,} bp")

    # Convert back to BedTool
    masked_bed = _bedtool_from_bedframe(masked_df)

    return masked_bed, bases_removed
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
		fields = f.fields
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

	# Convert to bioframe DataFrame
	df = _bedframe_from_bedtool(bt6)
	
	# Sort the DataFrame
	df_sorted = bf.sort_bedframe(df)
	
	# First, merge overlapping intervals (only chrom, start, end)
	merged_df = bf.merge(df_sorted)
	
	# For each merged interval, collect names from overlapping original intervals
	final_rows = []
	for _, merged_row in merged_df.iterrows():
		chrom, start, end = merged_row['chrom'], merged_row['start'], merged_row['end']
		
		# Find original intervals that overlap with this merged interval
		overlaps = bf.overlap(df_sorted, pd.DataFrame({
			'chrom': [chrom], 'start': [start], 'end': [end]
		}))
		
		# Collect unique names from overlapping intervals
		names = overlaps['name'].unique() if len(overlaps) > 0 else []
		name_str = ','.join(names) if len(names) > 0 else '.'
		
		final_rows.append({
			'chrom': chrom,
			'start': int(start),
			'end': int(end),
			'name': name_str,
			'score': 0,
			'strand': '.'
		})
	
	final_df = pd.DataFrame(final_rows)
	return _bedtool_from_bedframe(final_df)


def _compute_stats(bt: BedTool, genome: str, custom_cache_dir: Optional[Path] = None) -> Dict[str, float]:
	"""Compute stats for a BedTool (assumes at least bed3 columns present).

	Returns keys: count, total_bp, mean_bp, median_bp, genome_size_bp, percent_genome
	"""

	lens = [int(f.end) - int(f.start) for f in bt]
	count = len(lens)
	total_bp = int(sum(lens)) if lens else 0
	mean_bp = float(mean(lens)) if lens else 0.0
	median_bp = float(median(lens)) if lens else 0.0

	chromsizes = _get_cached_chromsizes(genome, custom_cache_dir)
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
		fields = f.fields
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
	repeat_regions_count: int = 0
	repeat_bases_masked: int = 0
	repeat_warnings: List[str] = None
	
	def __post_init__(self):
		if self.repeat_warnings is None:
			self.repeat_warnings = []


def run_pipeline(
	input_bed_path: Path,
	genome: str,
	buffer_size: int,
	ignore_cap: bool,
	strand_aware_buffer: bool,
	mask_repeats: bool = False,
	chunk_size: int = 400,
	custom_cache_dir: Optional[Path] = None,
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
	roi_stats = _compute_stats(roi_bed6, genome=genome, custom_cache_dir=custom_cache_dir)

	# Apply global cap based on ROI stats
	buffer_to_use, cap_applied = _cap_buffer_globally(roi_bed6, buffer_size, ignore_cap)

	# Fetch chromsizes once for clamping and pass to builder
	chromsizes_series = _get_cached_chromsizes(genome, custom_cache_dir)
	# Normalize values to int for safety (Series -> dict)
	chromsizes_int: Dict[str, int] = {str(k): int(v) for k, v in chromsizes_series.items()}
	adaptive_bed6, clamp_warnings = _build_adaptive_from_roi(
		roi_bed6, buffer_to_use, chromsizes_int, strand_aware_buffer
	)
	# Compute stats without double-counting ROI bases: use unstranded union if strand-aware
	if strand_aware_buffer:
		# Convert to bed3 lines and merge to union coverage using bioframe
		bed3_lines: List[str] = [f"{f.chrom}\t{int(f.start)}\t{int(f.end)}" for f in adaptive_bed6]
		bed3_df = _bedframe_from_bedtool(_bedtool_from_lines(bed3_lines))
		merged_df = bf.merge(bed3_df)
		merged_df = bf.sort_bedframe(merged_df)
		union_bt = _bedtool_from_bedframe(merged_df)
		adaptive_stats = _compute_stats(union_bt, genome=genome, custom_cache_dir=custom_cache_dir)
	else:
		adaptive_stats = _compute_stats(adaptive_bed6, genome=genome, custom_cache_dir=custom_cache_dir)

	# Handle repeat masking
	repeat_regions_count = 0
	repeat_bases_masked = 0
	repeat_warnings = []
	
	if mask_repeats:
		try:
			repeats_bed = _load_repeatmasker_bed(genome, min_size=chunk_size, custom_cache_dir=custom_cache_dir)
			repeat_regions_count = len(repeats_bed)
			
			# Check for overlaps using bioframe
			print("Checking for overlaps between adaptive regions and repeats...")
			adaptive_df = _bedframe_from_bedtool(adaptive_bed6)
			repeats_df = _bedframe_from_bedtool(repeats_bed)
			overlaps_df = bf.overlap(adaptive_df, repeats_df)
			overlapping_regions = len(overlaps_df)
			print(f"Found {overlapping_regions:,} overlapping regions")
			
			if overlapping_regions > 0:
				# Mask the repeats
				adaptive_bed6, repeat_bases_masked = _mask_repeats_in_bed(adaptive_bed6, repeats_bed)
				repeat_warnings.append(
					f"Masked {overlapping_regions} regions overlapping repeats, "
					f"{repeat_bases_masked} bases removed from adaptive BED."
				)
				# Recompute stats after masking
				if strand_aware_buffer:
					bed3_lines = [f"{f.chrom}\t{int(f.start)}\t{int(f.end)}" for f in adaptive_bed6]
					bed3_df = _bedframe_from_bedtool(_bedtool_from_lines(bed3_lines))
					merged_df = bf.merge(bed3_df)
					merged_df = bf.sort_bedframe(merged_df)
					union_bt = _bedtool_from_bedframe(merged_df)
					adaptive_stats = _compute_stats(union_bt, genome=genome)
				else:
					adaptive_stats = _compute_stats(adaptive_bed6, genome=genome)
			else:
				repeat_warnings.append("No overlaps with repeat regions found.")
		except Exception as e:
			repeat_warnings.append(f"Failed to process repeat data: {e}")
	else:
		# Still check for potential issues without masking
		try:
			repeats_bed = _load_repeatmasker_bed(genome, min_size=chunk_size, custom_cache_dir=custom_cache_dir)
			repeat_regions_count = len(repeats_bed)
			
			# Check for overlaps using bioframe
			print("Checking for overlaps between adaptive regions and repeats (warning only)...")
			adaptive_df = _bedframe_from_bedtool(adaptive_bed6)
			repeats_df = _bedframe_from_bedtool(repeats_bed)
			overlaps_df = bf.overlap(adaptive_df, repeats_df)
			overlapping_regions = len(overlaps_df)
			print(f"Found {overlapping_regions:,} overlapping regions")
			
			if overlapping_regions > 0:
				repeat_warnings.append(
					f"Warning: {overlapping_regions} regions overlap with repeats. "
					"Use --mask-repeat-regions to remove them."
				)
		except Exception as e:
			repeat_warnings.append(f"Failed to check repeat data: {e}")

	return PipelineOutputs(
		roi_bed6=roi_bed6,
		roi_stats=roi_stats,
		adaptive_bed6=adaptive_bed6,
		adaptive_stats=adaptive_stats,
		buffer_size_used=buffer_to_use,
		global_cap_applied=cap_applied,
		clamp_warnings=clamp_warnings,
		repeat_regions_count=repeat_regions_count,
		repeat_bases_masked=repeat_bases_masked,
		repeat_warnings=repeat_warnings,
	)

