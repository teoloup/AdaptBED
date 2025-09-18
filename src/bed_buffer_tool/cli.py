import argparse
from pathlib import Path
from typing import Optional, List

from .bed_ops import run_pipeline, compute_buffer_from_mean_cv, nxx_from_mean_cv  # import processing


# Parse a flexible boolean from string input for Python 3.8 compatibility.
def str2bool(s: str) -> bool:
    """
    Accept common true/false spellings for CLI flags that take explicit values.
    Examples: true/false, yes/no, y/n, 1/0.
    """
    v = s.strip().lower()
    if v in {"true", "t", "yes", "y", "1"}:
        return True
    if v in {"false", "f", "no", "n", "0"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {s}. Use true/false.")


def positive_int(value: str) -> int:
    try:
        iv = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer")
    if iv <= 0:
        raise argparse.ArgumentTypeError("Must be > 0")
    return iv


def non_negative_int(value: str) -> int:
    try:
        iv = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer")
    if iv < 0:
        raise argparse.ArgumentTypeError("Must be >= 0")
    return iv


def percent_0_100(value: str) -> float:
    try:
        fv = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a number")
    if fv < 0.0 or fv > 100.0:
        raise argparse.ArgumentTypeError("Must be between 0 and 100 (percent)")
    return fv


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="adaptbed",
        description=(
            "AdaptBED — adaptive-sampling BED formatter: parse/merge ROIs, compute stats, and build targets."
        ),
    )

    # Required
    p.add_argument("-i", "--input", required=True, help="Path to input BED file")
    p.add_argument("-o", "--output", required=True, help="Output directory (created if missing)")
    p.add_argument("--prefix", required=True, help="Prefix for naming output files")

    # Optional
    p.add_argument(
        "-s",
        "--strand-aware-buffer",
        action="store_true",
        help=(
            "If set, create two entries per ROI (+ and -) then merge strand-aware; "
            "otherwise expand symmetrically with strand '.'"
        ),
    )
    p.add_argument("--genome", default="hg38", help="Genome build (default: hg38)")
    p.add_argument(
        "-b",
        "--buffer-size",
        type=non_negative_int,
        default=5000,
        help="Padding to add to regions in bp (default: 5000)",
    )
    p.add_argument(
        "--calculate-buffer",
        action="store_true",
        help=(
            "Calculate buffer size from --mean-length and --cv "
            "(overrides --buffer-size)"
        ),
    )
    p.add_argument(
        "--mean-length",
        type=positive_int,
        required=False,
        help="Mean library length (bp). Required if --calculate-buffer is set",
    )
    p.add_argument(
        "--cv",
        type=percent_0_100,
        required=False,
        help=(
            "Coefficient of variation in percent [0-100]. "
            "Required if --calculate-buffer is set"
        ),
    )
    p.add_argument(
        "--ignore-cap",
        action="store_true",
        help="Ignore buffer 3x cap (may greatly increase target size)",
    )
    p.add_argument(
        "--mask-repeat-regions",
        action="store_true",
        help="Mask repeat regions in adaptive BED that overlap with regions >= chunk-size",
    )
    p.add_argument(
        "--chunk-size",
        type=positive_int,
        default=400,
        help="Minimum repeat region size to consider for masking (default: 400 bp)",
    )
    p.add_argument(
        "--cache-dir",
        type=str,
        help="Custom directory for caching genomic data (RepeatMasker, chromsizes). "
             "If not specified, uses .adaptbed_cache in project directory",
    )

    return p


def validate_args(args: argparse.Namespace) -> None:
    # Validate file paths
    in_path = Path(args.input)
    if not in_path.exists() or not in_path.is_file():
        raise SystemExit(f"Error: --input file not found: {in_path}")

    # Ensure output directory exists (mkdir -p behavior)
    out_dir = Path(args.output)
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Folder created successfully: {out_dir}")
    except Exception as e:
        raise SystemExit(f"Error: cannot create --output directory '{out_dir}': {e}")

    # calculate-buffer dependency checks
    if args.calculate_buffer:
        if args.mean_length is None or args.cv is None:
            raise SystemExit(
                "Error: --calculate-buffer requires both --mean-length and --cv."
            )
    # If not calculating buffer, buffer-size is enforced by arg type (>=0)


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Validation and side-effect: ensure output directory exists
    validate_args(args)

    # Determine requested buffer size
    nxx_meta = {}
    if args.calculate_buffer:
        mean_len = float(args.mean_length)
        cv_pct = float(args.cv)
        # Compute Nxx values
        N50 = int(round(nxx_from_mean_cv(mean_len, cv_pct, 50.0)))
        N25 = int(round(nxx_from_mean_cv(mean_len, cv_pct, 25.0)))
        N15 = int(round(nxx_from_mean_cv(mean_len, cv_pct, 15.0)))
        N10 = int(round(nxx_from_mean_cv(mean_len, cv_pct, 10.0)))
        N01 = int(round(nxx_from_mean_cv(mean_len, cv_pct, 1.0)))
        requested_buffer = N15
        # Print a concise summary
        print(f"Calculated buffer-size (N15) from mean={int(mean_len)}bp, cv={cv_pct}%: {requested_buffer} bp")
        print(f"Nxx estimates (bp): N50={N50}, N25={N25}, N15={N15}, N10={N10}, N01={N01}")
        print("Nanopore suggests choosing between N25–N10; we use N15 by default.")
        # Prepare meta for stats file
        nxx_meta = {
            "mean_length": int(mean_len),
            "cv_percent": cv_pct,
            "N50": N50,
            "N25": N25,
            "N15": N15,
            "N10": N10,
            "N01": N01,
            "suggested_range": "N25-N10",
            "chosen": "N15",
        }
    else:
        requested_buffer = args.buffer_size
        print(f"Requested buffer-size: {requested_buffer} bp")

    # Run processing pipeline (it may apply a global cap)
    custom_cache_dir = Path(args.cache_dir) if args.cache_dir else None
    outputs = run_pipeline(
        input_bed_path=Path(args.input),
        genome=args.genome or "hg38",
        buffer_size=requested_buffer,
        ignore_cap=bool(args.ignore_cap),
        strand_aware_buffer=bool(args.strand_aware_buffer),
        mask_repeats=bool(args.mask_repeat_regions),
        chunk_size=args.chunk_size,
        custom_cache_dir=custom_cache_dir,
    )

    # Report buffer actually used
    if outputs.global_cap_applied and not args.ignore_cap:
        print(f"Buffer-size capped globally to {outputs.buffer_size_used} bp to keep total targets ≤ 3x ROI size")
    else:
        print(f"Buffer-size used: {outputs.buffer_size_used} bp")

    # Write outputs
    out_dir = Path(args.output)
    roi_path = out_dir / f"{args.prefix}.roi.bed"
    adaptive_path = out_dir / f"{args.prefix}.adaptive_sampling.bed"
    stats_path = out_dir / f"{args.prefix}.stats.tsv"

    outputs.roi_bed6.saveas(str(roi_path))
    outputs.adaptive_bed6.saveas(str(adaptive_path))

    # Print stats to stdout
    def print_stats(label: str, stats: dict) -> None:
        print(f"[{label}] count={int(stats['count'])} total_bp={int(stats['total_bp'])} "
              f"mean_bp={stats['mean_bp']:.2f} median_bp={stats['median_bp']:.2f} "
              f"genome_size_bp={int(stats['genome_size_bp'])} percent_genome={stats['percent_genome']:.6f}")

    print_stats("roi", outputs.roi_stats)
    print_stats("adaptive", outputs.adaptive_stats)

    # Assemble stats output: ROI then adaptive
    def fmt_stats(prefix: str, stats: dict) -> str:
        return "\n".join(
            [
                f"{prefix}\tcount\t{int(stats['count'])}",
                f"{prefix}\ttotal_bp\t{int(stats['total_bp'])}",
                f"{prefix}\tmean_bp\t{stats['mean_bp']:.2f}",
                f"{prefix}\tmedian_bp\t{stats['median_bp']:.2f}",
                f"{prefix}\tgenome_size_bp\t{int(stats['genome_size_bp'])}",
                f"{prefix}\tpercent_genome\t{stats['percent_genome']:.6f}",
            ]
        )

    stats_lines = [
        "metric\tkey\tvalue",
        fmt_stats("roi", outputs.roi_stats),
        fmt_stats("adaptive", outputs.adaptive_stats),
        f"meta\tbuffer_size_used\t{outputs.buffer_size_used}",
        f"meta\trepeat_regions_count\t{outputs.repeat_regions_count}",
        f"meta\trepeat_bases_masked\t{outputs.repeat_bases_masked}",
    ]

    # Include Nxx meta if calculated
    if args.calculate_buffer and nxx_meta:
        stats_lines.extend(
            [
                f"meta\tmean_length\t{nxx_meta['mean_length']}",
                f"meta\tcv_percent\t{nxx_meta['cv_percent']}",
                f"meta\tN50\t{nxx_meta['N50']}",
                f"meta\tN25\t{nxx_meta['N25']}",
                f"meta\tN15\t{nxx_meta['N15']}",
                f"meta\tN10\t{nxx_meta['N10']}",
                f"meta\tN01\t{nxx_meta['N01']}",
                f"meta\tsuggested_range\t{nxx_meta['suggested_range']}",
                f"meta\tchosen\t{nxx_meta['chosen']}",
            ]
        )

    stats_path.write_text("\n".join(stats_lines) + "\n")

    # Warnings
    if outputs.global_cap_applied and not args.ignore_cap:
        print(
            "Warning: global buffer cap applied so total targets ≤ 3x total ROI. "
            "Use --ignore-cap to disable this behavior.")

    for w in outputs.clamp_warnings:
        print(f"Warning: {w}")

    for w in outputs.repeat_warnings:
        print(f"Repeat: {w}")

    print(f"Written: {roi_path}")
    print(f"Written: {adaptive_path}")
    print(f"Written: {stats_path}")
    return 0