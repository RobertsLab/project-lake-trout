#!/usr/bin/env python3
"""
12-pav.py - Present-Absent Variation (PAV) Analysis

This script identifies genomic regions that are present or absent in each sample
relative to the reference genome (GCF_016432855.1_SaNama_1.0).

PAV analysis uses coverage depth from aligned PacBio HiFi reads to identify:
1. Absent regions: Genomic regions with zero/very low coverage (potential deletions)
2. Present regions: Genomic regions with sufficient coverage
3. Insertions: Novel sequences in samples not in reference (from CIGAR analysis)
4. Deletions: Sequences in reference absent from samples (from CIGAR analysis)

Author: Generated for project-lake-trout
Date: 2025-12-07
"""

import os
import subprocess
import gzip
import re
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# =============================================================================
# Configuration
# =============================================================================

# Define samples by subspecies
# LEAN: bc2041, bc2069, bc2070, bc2068
# SISCOWET: bc2071, bc2073, bc2072, bc2096
SAMPLES_LEAN = ["bc2041", "bc2069", "bc2070", "bc2068"]
SAMPLES_SISCOWET = ["bc2071", "bc2073", "bc2072", "bc2096"]
SAMPLES = SAMPLES_LEAN + SAMPLES_SISCOWET  # All 8 samples

# Paths (relative to code directory)
BASE_DIR = Path(__file__).parent.parent
BAM_DIR = BASE_DIR / "analyses" / "05-pacbio-align"
OUTPUT_DIR = BASE_DIR / "analyses" / "11-pav"
REFERENCE = BASE_DIR / "data" / "GCF_016432855.1_SaNama_1.0_genomic.fa"
GENES_BED = BASE_DIR / "data" / "20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed"

# Analysis parameters
WINDOW_SIZE = 1000  # Window size for coverage calculation
MIN_COVERAGE = 1    # Minimum coverage to consider region "present"
MIN_INS_SIZE = 50   # Minimum insertion size to report
MIN_DEL_SIZE = 50   # Minimum deletion size to report
THREADS = 16        # Number of threads for tools

# =============================================================================
# Helper Functions
# =============================================================================

def run_command(cmd, description=""):
    """Run a shell command and handle errors."""
    print(f"  Running: {description or cmd[:80]}...")
    try:
        result = subprocess.run(
            cmd, shell=True, check=True,
            capture_output=True, text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: {e.stderr}")
        raise

def ensure_directory(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)

def file_exists_and_nonempty(filepath):
    """Check if file exists and is not empty."""
    path = Path(filepath)
    return path.exists() and path.stat().st_size > 0

# =============================================================================
# Step 1: Generate Reference Genome Index and Chromosome Sizes
# =============================================================================

def setup_reference():
    """Index reference genome and create chromosome sizes file."""
    print("\n" + "="*60)
    print("Step 1: Setting up reference genome files")
    print("="*60)
    
    # Index reference if not already done
    fai_file = f"{REFERENCE}.fai"
    if not Path(fai_file).exists():
        print("  Indexing reference genome...")
        run_command(f"samtools faidx {REFERENCE}", "samtools faidx")
    else:
        print("  Reference index already exists.")
    
    # Create chromosome sizes file
    chrom_sizes = OUTPUT_DIR / "genome.chrom.sizes"
    run_command(
        f"cut -f1,2 {fai_file} > {chrom_sizes}",
        "Creating chromosome sizes file"
    )
    print(f"  Chromosome sizes saved to: {chrom_sizes}")
    
    return chrom_sizes

# =============================================================================
# Step 2: Calculate Coverage with samtools/bedtools
# =============================================================================

def create_genome_windows(chrom_sizes_file, window_size):
    """Create genome windows BED file."""
    windows_file = OUTPUT_DIR / "genome_windows.bed"
    cmd = f"bedtools makewindows -g {chrom_sizes_file} -w {window_size} > {windows_file}"
    run_command(cmd, "Creating genome windows")
    return windows_file

def calculate_sample_coverage(sample, windows_file):
    """Calculate coverage for a single sample using samtools and bedtools."""
    bam_file = BAM_DIR / f"{sample}.bam"
    regions_file = OUTPUT_DIR / f"{sample}.regions.bed.gz"
    
    if not bam_file.exists():
        print(f"  WARNING: BAM file not found for {sample}: {bam_file}")
        return None
    
    # Use bedtools coverage to get mean coverage per window
    # This is more efficient than per-base depth for large files
    temp_file = OUTPUT_DIR / f"{sample}.regions.bed"
    
    cmd = f"""bedtools coverage \
        -a {windows_file} \
        -b {bam_file} \
        -mean \
        > {temp_file}"""
    
    run_command(cmd, f"bedtools coverage for {sample}")
    
    # Gzip the output to match mosdepth format
    cmd = f"gzip -f {temp_file}"
    run_command(cmd, f"Compressing coverage file for {sample}")
    
    return sample

def calculate_coverage():
    """Calculate coverage for all samples using samtools/bedtools."""
    print("\n" + "="*60)
    print("Step 2: Calculating coverage with bedtools")
    print("="*60)
    
    # First create genome windows
    chrom_sizes = OUTPUT_DIR / "genome.chrom.sizes"
    windows_file = create_genome_windows(chrom_sizes, WINDOW_SIZE)
    print(f"  Created genome windows: {windows_file}")
    
    for sample in SAMPLES:
        print(f"\n  Processing {sample}...")
        calculate_sample_coverage(sample, windows_file)
    
    print("\n  Coverage calculation complete for all samples.")

# =============================================================================
# Step 3: Identify Absent Regions (Zero Coverage)
# =============================================================================

def extract_absent_regions(sample):
    """Extract regions with zero coverage for a sample."""
    regions_file = OUTPUT_DIR / f"{sample}.regions.bed.gz"
    absent_file = OUTPUT_DIR / f"{sample}.absent_regions.bed"
    absent_merged = OUTPUT_DIR / f"{sample}.absent_regions_merged.bed"
    
    if not regions_file.exists():
        print(f"  WARNING: Regions file not found for {sample}")
        return 0
    
    # Read and filter for zero coverage
    # bedtools coverage -mean output: chrom, start, end, mean_coverage
    absent_regions = []
    with gzip.open(regions_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                try:
                    cov = float(parts[3])
                except ValueError:
                    continue
                if cov == 0:
                    absent_regions.append((chrom, start, end, "absent", 0, "."))
    
    # Write absent regions
    with open(absent_file, 'w') as f:
        for region in absent_regions:
            f.write('\t'.join(map(str, region)) + '\n')
    
    # Merge adjacent regions using bedtools
    if absent_regions:
        cmd = f"""bedtools merge -i {absent_file} | \
            awk -v OFS='\\t' '{{print $1, $2, $3, "absent_region", 0, "."}}' > {absent_merged}"""
        run_command(cmd, f"Merging absent regions for {sample}")
    else:
        # Create empty file
        open(absent_merged, 'w').close()
    
    # Count merged regions
    with open(absent_merged, 'r') as f:
        count = sum(1 for _ in f)
    
    return count

def identify_absent_regions():
    """Identify absent regions for all samples."""
    print("\n" + "="*60)
    print("Step 3: Identifying absent regions (zero coverage)")
    print("="*60)
    
    for sample in SAMPLES:
        count = extract_absent_regions(sample)
        print(f"  {sample}: {count} absent regions")

# =============================================================================
# Step 4: Identify Present Regions (With Coverage)
# =============================================================================

def extract_present_regions(sample):
    """Extract regions with coverage for a sample."""
    regions_file = OUTPUT_DIR / f"{sample}.regions.bed.gz"
    present_file = OUTPUT_DIR / f"{sample}.present_regions.bed"
    present_merged = OUTPUT_DIR / f"{sample}.present_regions_merged.bed"
    
    if not regions_file.exists():
        print(f"  WARNING: Regions file not found for {sample}")
        return 0
    
    # Read and filter for coverage >= MIN_COVERAGE
    # bedtools coverage -mean output: chrom, start, end, mean_coverage
    present_regions = []
    with gzip.open(regions_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                try:
                    cov = float(parts[3])
                except ValueError:
                    continue
                if cov >= MIN_COVERAGE:
                    present_regions.append((chrom, start, end, "present", cov, "."))
    
    # Write present regions
    with open(present_file, 'w') as f:
        for region in present_regions:
            f.write('\t'.join(map(str, region)) + '\n')
    
    # Merge adjacent regions
    if present_regions:
        cmd = f"""bedtools merge -i {present_file} | \
            awk -v OFS='\\t' '{{print $1, $2, $3, "present_region", 0, "."}}' > {present_merged}"""
        run_command(cmd, f"Merging present regions for {sample}")
    else:
        open(present_merged, 'w').close()
    
    # Count merged regions
    with open(present_merged, 'r') as f:
        count = sum(1 for _ in f)
    
    return count

def identify_present_regions():
    """Identify present regions for all samples."""
    print("\n" + "="*60)
    print("Step 4: Identifying present regions (with coverage)")
    print("="*60)
    
    for sample in SAMPLES:
        count = extract_present_regions(sample)
        print(f"  {sample}: {count} present regions")

# =============================================================================
# Step 5: Extract Insertions from BAM Alignments
# =============================================================================

def parse_cigar_for_indels(cigar_string, ref_pos, min_ins_size, min_del_size):
    """Parse CIGAR string and extract insertions and deletions."""
    insertions = []
    deletions = []
    
    # Parse CIGAR: match numbers followed by operation letters
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    current_pos = ref_pos
    for match in pattern.finditer(cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        
        if op == 'I' and length >= min_ins_size:
            insertions.append((current_pos, current_pos + 1, length))
        elif op == 'D' and length >= min_del_size:
            deletions.append((current_pos, current_pos + length, length))
        
        # Advance reference position for consuming operations
        if op in ['M', 'D', 'N', '=', 'X']:
            current_pos += length
    
    return insertions, deletions

def extract_indels_from_bam(sample):
    """Extract insertions and deletions from BAM file."""
    bam_file = BAM_DIR / f"{sample}.bam"
    ins_raw = OUTPUT_DIR / f"{sample}.insertions_raw.bed"
    del_raw = OUTPUT_DIR / f"{sample}.deletions_raw.bed"
    ins_file = OUTPUT_DIR / f"{sample}.insertions.bed"
    del_file = OUTPUT_DIR / f"{sample}.deletions.bed"
    
    if not bam_file.exists():
        print(f"  WARNING: BAM file not found for {sample}")
        return 0, 0
    
    print(f"  Extracting indels from {sample}...")
    
    insertions = []
    deletions = []
    
    # Use samtools view to read BAM
    cmd = f"samtools view {bam_file}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    
    line_count = 0
    for line in proc.stdout:
        line_count += 1
        if line_count % 100000 == 0:
            print(f"    Processed {line_count} reads...")
        
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        
        chrom = parts[2]
        if chrom == '*':
            continue
        
        try:
            pos = int(parts[3])
        except ValueError:
            continue
        
        cigar = parts[5]
        
        ins, dels = parse_cigar_for_indels(cigar, pos, MIN_INS_SIZE, MIN_DEL_SIZE)
        
        for start, end, length in ins:
            insertions.append((chrom, start, end, f"insertion_{length}bp", length, "+"))
        
        for start, end, length in dels:
            deletions.append((chrom, start, end, f"deletion_{length}bp", length, "+"))
    
    proc.wait()
    
    # Sort and write insertions
    insertions.sort(key=lambda x: (x[0], x[1]))
    with open(ins_raw, 'w') as f:
        for ins in insertions:
            f.write('\t'.join(map(str, ins)) + '\n')
    
    # Sort and write deletions
    deletions.sort(key=lambda x: (x[0], x[1]))
    with open(del_raw, 'w') as f:
        for d in deletions:
            f.write('\t'.join(map(str, d)) + '\n')
    
    # Merge insertions
    if insertions:
        cmd = f"""bedtools merge -i {ins_raw} -c 4,5 -o count,sum | \
            awk -v OFS='\\t' '{{print $1, $2, $3, "insertion_cluster", $4, "+", $5}}' > {ins_file}"""
        run_command(cmd, f"Merging insertions for {sample}")
    else:
        open(ins_file, 'w').close()
    
    # Merge deletions
    if deletions:
        cmd = f"""bedtools merge -i {del_raw} -c 4,5 -o count,sum | \
            awk -v OFS='\\t' '{{print $1, $2, $3, "deletion_cluster", $4, "+", $5}}' > {del_file}"""
        run_command(cmd, f"Merging deletions for {sample}")
    else:
        open(del_file, 'w').close()
    
    # Count results
    ins_count = sum(1 for _ in open(ins_file)) if file_exists_and_nonempty(ins_file) else 0
    del_count = sum(1 for _ in open(del_file)) if file_exists_and_nonempty(del_file) else 0
    
    return ins_count, del_count

def extract_all_indels():
    """Extract insertions and deletions for all samples."""
    print("\n" + "="*60)
    print("Step 5: Extracting insertions and deletions from BAM files")
    print("="*60)
    
    for sample in SAMPLES:
        ins_count, del_count = extract_indels_from_bam(sample)
        print(f"  {sample}: {ins_count} insertions, {del_count} deletions")

# =============================================================================
# Step 6: Create Combined PAV Feature Files
# =============================================================================

def create_combined_pav(sample):
    """Create combined PAV feature file for a sample."""
    absent_file = OUTPUT_DIR / f"{sample}.absent_regions_merged.bed"
    ins_file = OUTPUT_DIR / f"{sample}.insertions.bed"
    del_file = OUTPUT_DIR / f"{sample}.deletions.bed"
    pav_file = OUTPUT_DIR / f"{sample}.pav_features.bed"
    
    features = []
    
    # Add absent regions
    if file_exists_and_nonempty(absent_file):
        with open(absent_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    features.append((chrom, start, end, "ABSENT", end - start, "."))
    
    # Add insertions
    if file_exists_and_nonempty(ins_file):
        with open(ins_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    score = parts[4] if len(parts) > 4 else "0"
                    strand = parts[5] if len(parts) > 5 else "."
                    features.append((chrom, start, end, "INSERTION", score, strand))
    
    # Add deletions
    if file_exists_and_nonempty(del_file):
        with open(del_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    score = parts[4] if len(parts) > 4 else "0"
                    strand = parts[5] if len(parts) > 5 else "."
                    features.append((chrom, start, end, "DELETION", score, strand))
    
    # Sort features
    features.sort(key=lambda x: (x[0], x[1]))
    
    # Write combined file with header
    with open(pav_file, 'w') as f:
        f.write("#chrom\tstart\tend\tfeature_type\tscore\tstrand\n")
        for feat in features:
            f.write('\t'.join(map(str, feat)) + '\n')
    
    return len(features)

def create_all_combined_pav():
    """Create combined PAV files for all samples."""
    print("\n" + "="*60)
    print("Step 6: Creating combined PAV feature files")
    print("="*60)
    
    for sample in SAMPLES:
        count = create_combined_pav(sample)
        print(f"  {sample}: {count} total PAV features")

# =============================================================================
# Step 7: Intersect PAV Features with Gene Annotations
# =============================================================================

def intersect_with_genes(sample):
    """Intersect PAV features with gene annotations."""
    absent_file = OUTPUT_DIR / f"{sample}.absent_regions_merged.bed"
    del_file = OUTPUT_DIR / f"{sample}.deletions.bed"
    ins_file = OUTPUT_DIR / f"{sample}.insertions.bed"
    
    results = {}
    
    # Intersect absent regions with genes
    if file_exists_and_nonempty(absent_file):
        out_file = OUTPUT_DIR / f"{sample}.genes_in_absent_regions.bed"
        cmd = f"bedtools intersect -a {GENES_BED} -b {absent_file} -wa -wb > {out_file}"
        run_command(cmd, f"Intersecting absent regions with genes for {sample}")
        results['absent'] = sum(1 for _ in open(out_file)) if file_exists_and_nonempty(out_file) else 0
    
    # Intersect deletions with genes
    if file_exists_and_nonempty(del_file):
        out_file = OUTPUT_DIR / f"{sample}.genes_with_deletions.bed"
        cmd = f"bedtools intersect -a {GENES_BED} -b {del_file} -wa -wb > {out_file}"
        run_command(cmd, f"Intersecting deletions with genes for {sample}")
        results['deletions'] = sum(1 for _ in open(out_file)) if file_exists_and_nonempty(out_file) else 0
    
    # Intersect insertions with genes
    if file_exists_and_nonempty(ins_file):
        out_file = OUTPUT_DIR / f"{sample}.genes_with_insertions.bed"
        cmd = f"bedtools intersect -a {GENES_BED} -b {ins_file} -wa -wb > {out_file}"
        run_command(cmd, f"Intersecting insertions with genes for {sample}")
        results['insertions'] = sum(1 for _ in open(out_file)) if file_exists_and_nonempty(out_file) else 0
    
    return results

def intersect_all_with_genes():
    """Intersect PAV features with genes for all samples."""
    print("\n" + "="*60)
    print("Step 7: Intersecting PAV features with gene annotations")
    print("="*60)
    
    for sample in SAMPLES:
        results = intersect_with_genes(sample)
        print(f"  {sample} genes affected:")
        for feat_type, count in results.items():
            print(f"    - {feat_type}: {count}")

# =============================================================================
# Step 8: Generate Summary Statistics
# =============================================================================

def generate_summary():
    """Generate summary statistics for all samples."""
    print("\n" + "="*60)
    print("Step 8: Generating summary statistics")
    print("="*60)
    
    summary_data = []
    
    for sample in SAMPLES:
        absent_file = OUTPUT_DIR / f"{sample}.absent_regions_merged.bed"
        del_file = OUTPUT_DIR / f"{sample}.deletions.bed"
        ins_file = OUTPUT_DIR / f"{sample}.insertions.bed"
        
        # Count absent regions
        n_absent = 0
        absent_bp = 0
        if file_exists_and_nonempty(absent_file):
            with open(absent_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        n_absent += 1
                        absent_bp += int(parts[2]) - int(parts[1])
        
        # Count deletions
        n_del = 0
        del_bp = 0
        if file_exists_and_nonempty(del_file):
            with open(del_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        n_del += 1
                        del_bp += int(parts[2]) - int(parts[1])
        
        # Count insertions
        n_ins = 0
        if file_exists_and_nonempty(ins_file):
            with open(ins_file, 'r') as f:
                n_ins = sum(1 for _ in f)
        
        summary_data.append({
            'sample': sample,
            'n_absent_regions': n_absent,
            'absent_total_bp': absent_bp,
            'n_deletions': n_del,
            'deletion_total_bp': del_bp,
            'n_insertions': n_ins
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(summary_data)
    summary_file = OUTPUT_DIR / "pav_summary_stats.csv"
    df.to_csv(summary_file, index=False)
    
    print("\nSummary Statistics:")
    print(df.to_string(index=False))
    print(f"\n  Summary saved to: {summary_file}")
    
    return df

# =============================================================================
# Step 9: Compare PAV Across Samples
# =============================================================================

def compare_samples():
    """Compare PAV regions across all samples."""
    print("\n" + "="*60)
    print("Step 9: Comparing PAV across samples")
    print("="*60)
    
    # Create list of absent region files
    absent_files = [str(OUTPUT_DIR / f"{s}.absent_regions_merged.bed") for s in SAMPLES]
    
    # Check all files exist
    existing_files = [f for f in absent_files if file_exists_and_nonempty(f)]
    
    if len(existing_files) < 2:
        print("  Not enough samples with absent regions for comparison.")
        return
    
    # Run multiinter
    multiinter_out = OUTPUT_DIR / "absent_regions_multiinter.bed"
    cmd = f"""bedtools multiinter \
        -i {' '.join(existing_files)} \
        -header \
        -names {' '.join(SAMPLES)} \
        > {multiinter_out}"""
    run_command(cmd, "Running bedtools multiinter")
    
    # Filter for regions absent in all samples
    all_samples_out = OUTPUT_DIR / "absent_in_all_samples.bed"
    n_samples = len(SAMPLES)
    cmd = f"awk '$4 == {n_samples}' {multiinter_out} > {all_samples_out}"
    run_command(cmd, "Filtering regions absent in all samples")
    
    # Filter for regions absent in only one sample
    one_sample_out = OUTPUT_DIR / "absent_in_one_sample.bed"
    cmd = f"awk '$4 == 1' {multiinter_out} > {one_sample_out}"
    run_command(cmd, "Filtering sample-specific absent regions")
    
    # Count results
    all_count = sum(1 for _ in open(all_samples_out)) if file_exists_and_nonempty(all_samples_out) else 0
    one_count = sum(1 for _ in open(one_sample_out)) if file_exists_and_nonempty(one_sample_out) else 0
    
    print(f"\n  Regions absent in all samples: {all_count}")
    print(f"  Regions absent in only one sample: {one_count}")

# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main function to run the complete PAV analysis."""
    print("="*60)
    print("Present-Absent Variation (PAV) Analysis")
    print("="*60)
    print(f"\nSamples: {', '.join(SAMPLES)}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Reference genome: {REFERENCE}")
    
    # Ensure output directory exists
    ensure_directory(OUTPUT_DIR)
    
    # Run analysis steps
    setup_reference()
    calculate_coverage()
    identify_absent_regions()
    identify_present_regions()
    extract_all_indels()
    create_all_combined_pav()
    intersect_all_with_genes()
    summary_df = generate_summary()
    compare_samples()
    
    print("\n" + "="*60)
    print("PAV Analysis Complete!")
    print("="*60)
    print(f"\nOutput files are in: {OUTPUT_DIR}")
    print("\nKey output files per sample:")
    print("  - {sample}.pav_features.bed      : Combined PAV features")
    print("  - {sample}.absent_regions_merged.bed : Zero-coverage regions")
    print("  - {sample}.insertions.bed        : Insertions ≥50bp")
    print("  - {sample}.deletions.bed         : Deletions ≥50bp")
    print("  - {sample}.genes_in_absent_regions.bed : Genes in absent regions")
    print("\nCross-sample files:")
    print("  - pav_summary_stats.csv          : Summary statistics")
    print("  - absent_regions_multiinter.bed  : Multi-sample comparison")
    print("  - absent_in_all_samples.bed      : Regions absent in all samples")
    print("  - absent_in_one_sample.bed       : Sample-specific absent regions")

if __name__ == "__main__":
    main()

