# AdaptBED — Adaptive-Sampling BED Formatter

Adaptive sampling BED helper built on bioframe with chrom sizes and pybedtools for BED operations. Creates optimized target regions for adaptive sampling in ONT sequencing by expanding regions of interest (ROIs) with intelligent buffering while respecting chromosome boundaries and optional repeat masking.

## Features

- **Flexible Input**: Accepts BED3/6/12 formats and normalizes to BED6
- **Intelligent Buffering**: Expands targets with configurable buffer sizes
- **Smart Capping**: Global cap prevents excessive target expansion (3× ROI size limit)
- **Strand-Aware Processing**: Optional strand-specific buffering and merging
- **Auto-Buffer Calculation**: Calculates optimal buffer from read length statistics using log-normal model
- **Repeat Masking**: Optional masking of repetitive regions from UCSC RepeatMasker
- **Chromosome Boundary Clamping**: Prevents intervals from extending beyond chromosome ends
- **Comprehensive Statistics**: Detailed stats for ROI and final targets (count, total bp, mean/median bp, % genome)

## Installation

### Requirements
- Python 3.8+
- Internet connection (for downloading repeat data when using `--mask-repeat-regions`)

### Install from Source

```bash
# Clone the repository
git clone <repository-url>
cd adaptive_sampling_tool

# Install in development mode
pip install -e .
```

### Dependencies
The following packages will be automatically installed:
- `pybedtools==0.12.0` - BED file operations
- `bioframe>=0.6` - Genome coordinate operations and chromsizes
- `pandas` - Data manipulation
- `requests` - HTTP downloads for repeat data

### Windows Notes
On Windows, pybedtools uses pure Python operations in this project. No external bedtools binary is required, avoiding common compilation issues.

## Basic Usage

### Simple Example
```bash
# Basic usage with default 5kb buffer
adaptbed -i input.bed -o results --prefix sample
```

### With Custom Cache Directory
```bash
# Use custom cache location for genomic data
adaptbed -i input.bed -o results --prefix sample --cache-dir /path/to/cache
```

### Auto-Calculate Buffer
```bash
# Calculate buffer from read statistics (uses N15 by default)
adaptbed -i input.bed -o results --prefix sample --calculate-buffer --mean-length 12000 --cv 30
```

### With Repeat Masking
```bash
# Mask repetitive regions to avoid enriching repeats
adaptbed -i input.bed -o results --prefix sample --mask-repeat-regions
```

## Parameters Explanation

### Required Parameters

| Parameter | Short | Description |
|-----------|-------|-------------|
| `--input` | `-i` | Input BED file path (BED3/6/12 format) |
| `--output` | `-o` | Output directory for results |
| `--prefix` | `-p` | Prefix for output files |

### Optional Parameters

#### Buffer Configuration
| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--buffer-size` | `-b` | Buffer size in base pairs to add to each side of ROIs | 5000 |
| `--calculate-buffer` | | Auto-calculate buffer from read statistics using log-normal model | False |
| `--mean-length` | | Mean read length for buffer calculation (required with `--calculate-buffer`) | None |
| `--cv` | | Coefficient of variation (%) for read lengths (required with `--calculate-buffer`) | None |
| `--ignore-cap` | | Ignore the global 3× cap on total target size | False |

#### Processing Options
| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--strand-aware-buffer` | `-s` | Use strand-aware buffering (+ left, - right) with strand-aware merging | False |
| `--mask-repeat-regions` | | Mask repeat regions from UCSC RepeatMasker to avoid repetitive sequences | False |
| `--chunk-size` | | Minimum size threshold for repeat regions (bp) | 400 |

#### Genome and Output
| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--genome` | `-g` | Genome build for chromosome sizes (supported by bioframe) | hg38 |
| `--cache-dir` | | Custom directory for caching genomic data (RepeatMasker, chromsizes) | Project's .adaptbed_cache |
| `--quiet` | `-q` | Suppress progress messages | False |

## Output Files

AdaptBED generates the following files in the output directory:

- **`<prefix>.roi.bed`**: Merged ROI BED6 file (cleaned and sorted regions of interest)
- **`<prefix>.adaptive_sampling.bed`**: Final adaptive sampling targets BED6 with distinct names
- **`<prefix>.stats.tsv`**: Comprehensive statistics and metadata

### Statistics Output
The stats file includes:
- ROI statistics (count, total bp, mean/median bp, % genome)
- Final target statistics
- Buffer size used
- Auto-calculation details (when `--calculate-buffer` is used)
- Repeat masking information (when `--mask-repeat-regions` is used)
- Clamping warnings (if any intervals were adjusted for chromosome boundaries)

## Advanced Examples

### Strand-Aware Processing
```bash
# Strand-aware buffering for directional sequencing
adaptbed -i input.bed -o results --prefix sample -s -b 3000
```

### Auto-Buffer with Custom N-Value
```bash
# Use N10 instead of default N15 for buffer calculation
adaptbed -i input.bed -o results --prefix sample --calculate-buffer --mean-length 15000 --cv 25
```

### Complete Pipeline with All Features
```bash
# Full processing pipeline
adaptbed -i input.bed -o results --prefix sample \
    --calculate-buffer --mean-length 12000 --cv 30 \
    --mask-repeat-regions --chunk-size 500 \
    --genome hg38
```

### High-Throughput Processing
```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
    adaptbed -i ${sample}.bed -o results --prefix ${sample} --mask-repeat-regions
done
```

## How Buffer and Cap Work

### Buffer Calculation
- **Manual**: Each ROI is expanded by `buffer-size` on both sides
- **Auto**: Buffer is calculated from read length distribution using log-normal model
  - N50/N25/N15/N10/N01 values are computed and recorded
  - N15 (15% of reads longer than this value) is used by default
  - Formula: `target_size = roi_size + 2 × roi_count × buffer`

### Global Cap
- Prevents excessive target expansion: `total_targets ≤ 3 × roi_size`
- When cap is reached: `buffer_used = min(requested_buffer, floor(roi_size/roi_count))`
- Can be disabled with `--ignore-cap` for specialized use cases

### Strand-Aware Processing
- **Forward strand (+)**: Buffer added to left side of intervals
- **Reverse strand (-)**: Buffer added to right side of intervals
- Merging preserves strand information and distinct names
- Statistics computed on unstranded union to avoid double-counting

## Caching Behavior

AdaptBED caches genomic data to improve performance and enable offline usage:

### Cache Location
- **Default**: `.adaptbed_cache/` in the project directory
- **Custom**: Specify with `--cache-dir /path/to/cache`
- **Contents**: RepeatMasker data and chromosome sizes

### Cached Data
- **RepeatMasker files**: Downloaded from UCSC (~150-200MB per genome)
- **Chromosome sizes**: Fetched from bioframe and cached locally
- **Automatic**: Data is cached on first use, reused on subsequent runs

### Cache Management
```bash
# Use custom cache location
adaptbed -i input.bed -o results --prefix sample --cache-dir /custom/cache/path

# Cache is automatically created and populated
# Delete .adaptbed_cache/ to clear all cached data
```

### Offline Usage
Once data is cached, AdaptBED can run without internet:
- Cached RepeatMasker data enables `--mask-repeat-regions`
- Cached chromosome sizes enable all genome operations
- Only initial run requires internet connection

## Troubleshooting

### Common Issues

**Import Errors**
```bash
# Missing dependencies
pip install bioframe pybedtools

# Reinstall package
pip install -e .
```

**Genome Not Found**
```bash
# Check available genomes
python -c "import bioframe as bf; print(bf.fetch_chromsizes())"

# Use supported genome
adaptbed -i input.bed -o results --prefix sample --genome hg19
```

**Repeat Data Download Issues**
- Ensure internet connection for initial repeat data download
- Repeat data is cached locally after first use
- Check UCSC genome browser for supported genomes

**Memory Issues with Large Files**
- Use smaller chunk sizes for repeat masking: `--chunk-size 200`
- Process large files in batches if needed

### Performance Tips
- Use `--quiet` for batch processing
- Auto-buffer calculation adds computation time
- Repeat masking requires initial download but is fast thereafter
- Strand-aware processing is more memory-intensive

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Citation

If you use AdaptBED in your research, please cite:
[Add citation information here]