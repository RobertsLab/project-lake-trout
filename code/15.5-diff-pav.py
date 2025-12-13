#!/usr/bin/env python3
"""
15.5-diff-pav.py - Stringent Differential PAV Analysis: Siscowet vs Lean Trout

This script identifies insertions and deletions that are specific to each 
ecotype (Siscowet vs Lean) using STRINGENT criteria:
  1. A PAV must be present in ALL 4 samples of one ecotype
  2. A PAV must be absent in ALL 4 samples of the other ecotype
  3. Only PAVs > 100bp are considered

Sample Classification:
    LEAN: bc2041, bc2069, bc2070, bc2068
    SISCOWET: bc2071, bc2073, bc2072, bc2096

Output files are formatted for genome browser visualization with:
- Color coding (blue for Lean-specific, red for Siscowet-specific)
- Descriptive names
- Track headers for UCSC/IGV compatibility

Author: Generated for project-lake-trout
Date: 2025-12-09
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
import tempfile

# =============================================================================
# Configuration
# =============================================================================

# Define samples by ecotype
SAMPLES_LEAN = ["bc2041", "bc2069", "bc2070", "bc2068"]
SAMPLES_SISCOWET = ["bc2071", "bc2073", "bc2072", "bc2096"]

# Paths
BASE_DIR = Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "analyses" / "11-pav"
OUTPUT_DIR = BASE_DIR / "analyses" / "15-diff-pav"

# BED file color coding (RGB format for genome browsers)
COLOR_LEAN = "0,0,255"       # Blue for Lean-specific
COLOR_SISCOWET = "255,0,0"   # Red for Siscowet-specific
COLOR_SHARED = "128,128,128" # Gray for shared variants

# Stringent parameters
MIN_OVERLAP_FRACTION = 0.5  # 50% reciprocal overlap to consider same variant
MIN_PAV_SIZE = 100          # Minimum PAV size in bp (>100bp)
REQUIRE_ALL_SAMPLES = True  # Require presence in ALL samples of ecotype

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

def get_available_samples(samples, variant_type):
    """Get list of samples that have data for the variant type."""
    available = []
    for sample in samples:
        filepath = INPUT_DIR / f"{sample}.{variant_type}.bed"
        if file_exists_and_nonempty(filepath):
            available.append(sample)
    return available

def count_bed_lines(filepath):
    """Count non-header lines in a BED file."""
    if not file_exists_and_nonempty(filepath):
        return 0
    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('track'):
                count += 1
    return count

def load_bed_file(filepath, min_size=0):
    """
    Load a BED file and return list of variants.
    Optionally filter by minimum size.
    """
    variants = []
    if not file_exists_and_nonempty(filepath):
        return variants
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                size = end - start
                
                # Filter by minimum size
                if size > min_size:
                    variants.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'size': size,
                        'line': line.strip()
                    })
    
    return variants

def write_bed_file(variants, filepath):
    """Write variants to a BED file."""
    with open(filepath, 'w') as f:
        for v in variants:
            f.write(f"{v['chrom']}\t{v['start']}\t{v['end']}\n")

def variants_overlap(v1, v2, min_overlap_frac=0.5):
    """
    Check if two variants overlap with at least min_overlap_frac reciprocal overlap.
    """
    if v1['chrom'] != v2['chrom']:
        return False
    
    # Calculate overlap
    overlap_start = max(v1['start'], v2['start'])
    overlap_end = min(v1['end'], v2['end'])
    
    if overlap_start >= overlap_end:
        return False
    
    overlap_size = overlap_end - overlap_start
    
    # Calculate overlap fraction for both variants
    frac1 = overlap_size / v1['size']
    frac2 = overlap_size / v2['size']
    
    # Reciprocal overlap: both must meet threshold
    return frac1 >= min_overlap_frac and frac2 >= min_overlap_frac

# =============================================================================
# Step 1: Find variants present in ALL samples of an ecotype (stringent)
# =============================================================================

def find_variants_in_all_samples(samples, variant_type):
    """
    Find variants that are present in ALL samples of an ecotype.
    Uses iterative intersection approach.
    Returns list of variants present in all samples.
    """
    available = get_available_samples(samples, variant_type)
    
    if len(available) < len(samples):
        print(f"  WARNING: Only {len(available)}/{len(samples)} samples available for {variant_type}")
        print(f"    Available: {available}")
        print(f"    Missing: {set(samples) - set(available)}")
        if len(available) == 0:
            return None
    
    print(f"  Processing {len(available)} samples for {variant_type}")
    
    # Load first sample as starting point (filtered by size)
    first_sample = available[0]
    first_file = INPUT_DIR / f"{first_sample}.{variant_type}.bed"
    
    # Create temporary directory for intermediate files
    temp_dir = OUTPUT_DIR / "temp"
    ensure_directory(temp_dir)
    
    # Filter first file by size and sort
    current_file = temp_dir / f"current.{variant_type}.bed"
    
    # Load and filter by size
    variants = load_bed_file(first_file, min_size=MIN_PAV_SIZE)
    print(f"    {first_sample}: {len(variants)} variants > {MIN_PAV_SIZE}bp")
    
    if not variants:
        return None
    
    # Write filtered variants
    write_bed_file(variants, current_file)
    
    # Sort
    sorted_file = temp_dir / f"current.sorted.{variant_type}.bed"
    cmd = f"sort -k1,1 -k2,2n {current_file} > {sorted_file}"
    run_command(cmd, f"Sorting {first_sample}")
    
    current_file = sorted_file
    
    # Iteratively intersect with remaining samples
    for sample in available[1:]:
        sample_file = INPUT_DIR / f"{sample}.{variant_type}.bed"
        
        # Filter sample by size
        sample_filtered = temp_dir / f"{sample}.filtered.{variant_type}.bed"
        variants = load_bed_file(sample_file, min_size=MIN_PAV_SIZE)
        print(f"    {sample}: {len(variants)} variants > {MIN_PAV_SIZE}bp")
        
        if not variants:
            print(f"    WARNING: No variants > {MIN_PAV_SIZE}bp in {sample}")
            # Clean up and return None - can't have variants in ALL samples
            return None
        
        write_bed_file(variants, sample_filtered)
        
        # Sort sample file
        sample_sorted = temp_dir / f"{sample}.sorted.{variant_type}.bed"
        cmd = f"sort -k1,1 -k2,2n {sample_filtered} > {sample_sorted}"
        run_command(cmd, f"Sorting {sample}")
        
        # Intersect with current set (keep only variants in BOTH)
        # -f and -r for reciprocal overlap
        # -wa writes the original entry from file A
        next_file = temp_dir / f"intersect.{sample}.{variant_type}.bed"
        cmd = f"""bedtools intersect \
            -a {current_file} \
            -b {sample_sorted} \
            -wa \
            -f {MIN_OVERLAP_FRACTION} \
            -r \
            -u \
            > {next_file}"""
        run_command(cmd, f"Intersecting with {sample}")
        
        # Count remaining variants
        remaining = count_bed_lines(next_file)
        print(f"    After intersecting with {sample}: {remaining} variants remain")
        
        if remaining == 0:
            print(f"    No variants present in all samples up to {sample}")
            return None
        
        current_file = next_file
    
    # Final result
    final_count = count_bed_lines(current_file)
    print(f"  Found {final_count} variants present in ALL {len(available)} samples")
    
    # Copy to output
    output_file = OUTPUT_DIR / f"stringent.all_{len(available)}_samples.{variant_type}.bed"
    cmd = f"cp {current_file} {output_file}"
    run_command(cmd, "Copying final result")
    
    return output_file

# =============================================================================
# Step 2: Find ecotype-specific variants (stringent)
# =============================================================================

def find_stringent_ecotype_specific(ecotype_all_file, other_ecotype_files, 
                                    output_file, ecotype_name, variant_type):
    """
    Find variants present in ALL samples of one ecotype and ABSENT from ALL samples 
    of the other ecotype.
    
    ecotype_all_file: BED file with variants present in ALL samples of this ecotype
    other_ecotype_files: list of BED files for each sample of the other ecotype
    """
    if ecotype_all_file is None or not file_exists_and_nonempty(ecotype_all_file):
        print(f"  No variants present in all {ecotype_name} samples")
        return None
    
    # Create temp directory
    temp_dir = OUTPUT_DIR / "temp"
    ensure_directory(temp_dir)
    
    current_file = ecotype_all_file
    
    # Subtract each sample from the other ecotype
    for other_file in other_ecotype_files:
        if not file_exists_and_nonempty(other_file):
            continue
        
        sample_name = Path(other_file).stem.split('.')[0]
        
        # Filter other sample by size first
        other_filtered = temp_dir / f"other.{sample_name}.filtered.bed"
        variants = load_bed_file(other_file, min_size=MIN_PAV_SIZE)
        
        if not variants:
            print(f"    {sample_name}: no variants > {MIN_PAV_SIZE}bp (skipping)")
            continue
        
        write_bed_file(variants, other_filtered)
        
        # Sort
        other_sorted = temp_dir / f"other.{sample_name}.sorted.bed"
        cmd = f"sort -k1,1 -k2,2n {other_filtered} > {other_sorted}"
        run_command(cmd, f"Sorting {sample_name}")
        
        # Remove any overlapping variants
        next_file = temp_dir / f"subtract.{sample_name}.bed"
        cmd = f"""bedtools intersect \
            -a {current_file} \
            -b {other_sorted} \
            -v \
            -f {MIN_OVERLAP_FRACTION} \
            -r \
            > {next_file}"""
        run_command(cmd, f"Subtracting {sample_name}")
        
        remaining = count_bed_lines(next_file)
        print(f"    After removing overlaps with {sample_name}: {remaining} variants remain")
        
        current_file = next_file
    
    # Copy final result
    final_count = count_bed_lines(current_file)
    if final_count > 0:
        cmd = f"cp {current_file} {output_file}"
        run_command(cmd, f"Saving {ecotype_name}-specific variants")
        print(f"  Found {final_count} {ecotype_name}-specific {variant_type}")
        return output_file
    else:
        print(f"  No {ecotype_name}-specific {variant_type} found")
        return None

# =============================================================================
# Step 3: Create genome browser tracks
# =============================================================================

def create_genome_browser_bed(input_file, output_file, track_name, 
                               track_description, color, variant_type):
    """
    Format BED file for genome browser visualization with track header.
    Creates BED9 format with color coding.
    """
    if not file_exists_and_nonempty(input_file):
        print(f"  Skipping empty file: {input_file}")
        return None
    
    # Read input
    variants = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                
                # Calculate variant size
                size = end - start
                
                # Create descriptive name
                if variant_type == "insertions":
                    name = f"INS_{size}bp"
                elif variant_type == "deletions":
                    name = f"DEL_{size}bp"
                else:
                    name = f"VAR_{size}bp"
                
                # Score (use count if available, else 0)
                score = parts[3] if len(parts) > 3 else "0"
                try:
                    score = min(int(float(score)), 1000)  # Cap at 1000 for BED format
                except:
                    score = 0
                
                variants.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': '.',
                    'thickStart': start,
                    'thickEnd': end,
                    'itemRgb': color
                })
    
    # Sort by position
    variants.sort(key=lambda x: (x['chrom'], x['start']))
    
    # Write output with track header
    with open(output_file, 'w') as f:
        # Write track header for UCSC genome browser
        f.write(f'track name="{track_name}" description="{track_description}" ')
        f.write(f'visibility=pack itemRgb="On"\n')
        
        for v in variants:
            f.write(f"{v['chrom']}\t{v['start']}\t{v['end']}\t{v['name']}\t")
            f.write(f"{v['score']}\t{v['strand']}\t{v['thickStart']}\t{v['thickEnd']}\t")
            f.write(f"{v['itemRgb']}\n")
    
    return output_file

# =============================================================================
# Step 4: Generate summary statistics
# =============================================================================

def generate_summary():
    """Generate summary statistics for stringent differential PAV analysis."""
    summary_data = []
    
    variant_types = ['insertions', 'deletions']
    categories = ['lean_specific', 'siscowet_specific']
    
    for vtype in variant_types:
        for category in categories:
            bed_file = OUTPUT_DIR / f"stringent.{category}.{vtype}.browser.bed"
            count = count_bed_lines(bed_file)
            
            # Calculate total bp and size stats
            total_bp = 0
            sizes = []
            if file_exists_and_nonempty(bed_file):
                with open(bed_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or line.startswith('track'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            size = int(parts[2]) - int(parts[1])
                            total_bp += size
                            sizes.append(size)
            
            min_size = min(sizes) if sizes else 0
            max_size = max(sizes) if sizes else 0
            mean_size = sum(sizes) / len(sizes) if sizes else 0
            
            summary_data.append({
                'variant_type': vtype,
                'category': category,
                'count': count,
                'total_bp': total_bp,
                'min_size': min_size,
                'max_size': max_size,
                'mean_size': round(mean_size, 1)
            })
    
    df = pd.DataFrame(summary_data)
    summary_file = OUTPUT_DIR / "stringent_diff_pav_summary.csv"
    df.to_csv(summary_file, index=False)
    
    return df

# =============================================================================
# Step 5: Create combined visualization tracks
# =============================================================================

def create_combined_tracks():
    """
    Create combined track files with all variants colored by ecotype.
    """
    print("\n" + "="*60)
    print("Creating combined visualization tracks")
    print("="*60)
    
    for vtype in ['insertions', 'deletions']:
        combined_file = OUTPUT_DIR / f"stringent.all_{vtype}_by_ecotype.browser.bed"
        
        all_variants = []
        
        # Add Lean-specific (blue)
        lean_file = OUTPUT_DIR / f"stringent.lean_specific.{vtype}.browser.bed"
        if file_exists_and_nonempty(lean_file):
            with open(lean_file, 'r') as f:
                for line in f:
                    if not line.startswith('track') and line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            parts[3] = f"LEAN_{parts[3]}"  # Prefix name
                            all_variants.append(parts)
        
        # Add Siscowet-specific (red)
        sisc_file = OUTPUT_DIR / f"stringent.siscowet_specific.{vtype}.browser.bed"
        if file_exists_and_nonempty(sisc_file):
            with open(sisc_file, 'r') as f:
                for line in f:
                    if not line.startswith('track') and line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            parts[3] = f"SISC_{parts[3]}"  # Prefix name
                            all_variants.append(parts)
        
        # Sort all variants
        all_variants.sort(key=lambda x: (x[0], int(x[1])))
        
        # Write combined file
        with open(combined_file, 'w') as f:
            vtype_cap = vtype.capitalize()
            f.write(f'track name="Stringent_All_{vtype_cap}_by_Ecotype" ')
            f.write(f'description="Stringent {vtype} (>100bp, in all 4 samples): Blue=Lean, Red=Siscowet" ')
            f.write('visibility=pack itemRgb="On"\n')
            
            for v in all_variants:
                f.write('\t'.join(v) + '\n')
        
        print(f"  Created: {combined_file}")
        print(f"    Total variants: {len(all_variants)}")

# =============================================================================
# Cleanup
# =============================================================================

def cleanup_temp_files():
    """Remove temporary files."""
    temp_dir = OUTPUT_DIR / "temp"
    if temp_dir.exists():
        import shutil
        shutil.rmtree(temp_dir)
        print("  Cleaned up temporary files")

# =============================================================================
# Main Analysis
# =============================================================================

def main():
    """Run the stringent differential PAV analysis."""
    print("="*60)
    print("STRINGENT Differential PAV Analysis: Siscowet vs Lean Trout")
    print("="*60)
    
    print(f"\nStringent Criteria:")
    print(f"  1. PAV must be present in ALL 4 samples of one ecotype")
    print(f"  2. PAV must be absent from ALL 4 samples of other ecotype")
    print(f"  3. PAV size must be > {MIN_PAV_SIZE}bp")
    print(f"  4. Reciprocal overlap threshold: {MIN_OVERLAP_FRACTION*100}%")
    
    print(f"\nLean samples: {', '.join(SAMPLES_LEAN)}")
    print(f"Siscowet samples: {', '.join(SAMPLES_SISCOWET)}")
    print(f"Input directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Create output directory
    ensure_directory(OUTPUT_DIR)
    
    # Check available data
    print("\n" + "="*60)
    print("Checking available data")
    print("="*60)
    
    lean_ins = get_available_samples(SAMPLES_LEAN, "insertions")
    lean_del = get_available_samples(SAMPLES_LEAN, "deletions")
    sisc_ins = get_available_samples(SAMPLES_SISCOWET, "insertions")
    sisc_del = get_available_samples(SAMPLES_SISCOWET, "deletions")
    
    print(f"  Lean insertions: {len(lean_ins)}/4 samples available")
    print(f"  Lean deletions: {len(lean_del)}/4 samples available")
    print(f"  Siscowet insertions: {len(sisc_ins)}/4 samples available")
    print(f"  Siscowet deletions: {len(sisc_del)}/4 samples available")
    
    # Get file paths for other ecotype samples
    lean_ins_files = [INPUT_DIR / f"{s}.insertions.bed" for s in SAMPLES_LEAN]
    lean_del_files = [INPUT_DIR / f"{s}.deletions.bed" for s in SAMPLES_LEAN]
    sisc_ins_files = [INPUT_DIR / f"{s}.insertions.bed" for s in SAMPLES_SISCOWET]
    sisc_del_files = [INPUT_DIR / f"{s}.deletions.bed" for s in SAMPLES_SISCOWET]
    
    # ==========================================================================
    # Process Insertions
    # ==========================================================================
    print("\n" + "="*60)
    print("Step 1: Processing INSERTIONS (stringent)")
    print("="*60)
    
    # Find insertions in ALL Lean samples
    print("\n  Finding insertions present in ALL Lean samples...")
    lean_ins_all = find_variants_in_all_samples(SAMPLES_LEAN, "insertions")
    
    # Find insertions in ALL Siscowet samples
    print("\n  Finding insertions present in ALL Siscowet samples...")
    sisc_ins_all = find_variants_in_all_samples(SAMPLES_SISCOWET, "insertions")
    
    # Find Lean-specific insertions (in all Lean, absent from all Siscowet)
    print("\n  Finding LEAN-specific insertions...")
    lean_specific_ins = OUTPUT_DIR / "stringent.lean_specific.insertions.bed"
    find_stringent_ecotype_specific(
        lean_ins_all, sisc_ins_files,
        lean_specific_ins, "Lean", "insertions"
    )
    
    # Find Siscowet-specific insertions (in all Siscowet, absent from all Lean)
    print("\n  Finding SISCOWET-specific insertions...")
    sisc_specific_ins = OUTPUT_DIR / "stringent.siscowet_specific.insertions.bed"
    find_stringent_ecotype_specific(
        sisc_ins_all, lean_ins_files,
        sisc_specific_ins, "Siscowet", "insertions"
    )
    
    # Create genome browser tracks
    print("\n  Creating genome browser tracks for insertions...")
    
    create_genome_browser_bed(
        lean_specific_ins,
        OUTPUT_DIR / "stringent.lean_specific.insertions.browser.bed",
        "Stringent_Lean_Insertions",
        f"Insertions in ALL 4 Lean samples, absent from ALL Siscowet (>{MIN_PAV_SIZE}bp, blue)",
        COLOR_LEAN,
        "insertions"
    )
    
    create_genome_browser_bed(
        sisc_specific_ins,
        OUTPUT_DIR / "stringent.siscowet_specific.insertions.browser.bed",
        "Stringent_Siscowet_Insertions",
        f"Insertions in ALL 4 Siscowet samples, absent from ALL Lean (>{MIN_PAV_SIZE}bp, red)",
        COLOR_SISCOWET,
        "insertions"
    )
    
    # ==========================================================================
    # Process Deletions
    # ==========================================================================
    print("\n" + "="*60)
    print("Step 2: Processing DELETIONS (stringent)")
    print("="*60)
    
    # Find deletions in ALL Lean samples
    print("\n  Finding deletions present in ALL Lean samples...")
    lean_del_all = find_variants_in_all_samples(SAMPLES_LEAN, "deletions")
    
    # Find deletions in ALL Siscowet samples
    print("\n  Finding deletions present in ALL Siscowet samples...")
    sisc_del_all = find_variants_in_all_samples(SAMPLES_SISCOWET, "deletions")
    
    # Find Lean-specific deletions (in all Lean, absent from all Siscowet)
    print("\n  Finding LEAN-specific deletions...")
    lean_specific_del = OUTPUT_DIR / "stringent.lean_specific.deletions.bed"
    find_stringent_ecotype_specific(
        lean_del_all, sisc_del_files,
        lean_specific_del, "Lean", "deletions"
    )
    
    # Find Siscowet-specific deletions (in all Siscowet, absent from all Lean)
    print("\n  Finding SISCOWET-specific deletions...")
    sisc_specific_del = OUTPUT_DIR / "stringent.siscowet_specific.deletions.bed"
    find_stringent_ecotype_specific(
        sisc_del_all, lean_del_files,
        sisc_specific_del, "Siscowet", "deletions"
    )
    
    # Create genome browser tracks
    print("\n  Creating genome browser tracks for deletions...")
    
    create_genome_browser_bed(
        lean_specific_del,
        OUTPUT_DIR / "stringent.lean_specific.deletions.browser.bed",
        "Stringent_Lean_Deletions",
        f"Deletions in ALL 4 Lean samples, absent from ALL Siscowet (>{MIN_PAV_SIZE}bp, blue)",
        COLOR_LEAN,
        "deletions"
    )
    
    create_genome_browser_bed(
        sisc_specific_del,
        OUTPUT_DIR / "stringent.siscowet_specific.deletions.browser.bed",
        "Stringent_Siscowet_Deletions",
        f"Deletions in ALL 4 Siscowet samples, absent from ALL Lean (>{MIN_PAV_SIZE}bp, red)",
        COLOR_SISCOWET,
        "deletions"
    )
    
    # ==========================================================================
    # Create combined visualization tracks
    # ==========================================================================
    create_combined_tracks()
    
    # ==========================================================================
    # Generate summary statistics
    # ==========================================================================
    print("\n" + "="*60)
    print("Generating Summary Statistics")
    print("="*60)
    
    summary_df = generate_summary()
    print("\nStringent Differential PAV Summary:")
    print(summary_df.to_string(index=False))
    
    # ==========================================================================
    # Cleanup
    # ==========================================================================
    print("\n" + "="*60)
    print("Cleaning up")
    print("="*60)
    cleanup_temp_files()
    
    # ==========================================================================
    # Final output
    # ==========================================================================
    print("\n" + "="*60)
    print("Stringent Differential PAV Analysis Complete!")
    print("="*60)
    
    print(f"\nStringent Criteria Applied:")
    print(f"  - PAV present in ALL 4 samples of one ecotype")
    print(f"  - PAV absent from ALL 4 samples of other ecotype")
    print(f"  - PAV size > {MIN_PAV_SIZE}bp")
    
    print(f"\nOutput files in: {OUTPUT_DIR}")
    print("\nGenome Browser Track Files (stringent):")
    print("  Lean-specific (Blue):")
    print("    - stringent.lean_specific.insertions.browser.bed")
    print("    - stringent.lean_specific.deletions.browser.bed")
    print("  Siscowet-specific (Red):")
    print("    - stringent.siscowet_specific.insertions.browser.bed")
    print("    - stringent.siscowet_specific.deletions.browser.bed")
    print("  Combined tracks:")
    print("    - stringent.all_insertions_by_ecotype.browser.bed")
    print("    - stringent.all_deletions_by_ecotype.browser.bed")
    print("\nSummary file:")
    print("    - stringent_diff_pav_summary.csv")
    print("\nTo view in UCSC Genome Browser:")
    print("  1. Go to https://genome.ucsc.edu/cgi-bin/hgCustom")
    print("  2. Choose the appropriate genome assembly")
    print("  3. Upload the .browser.bed files as custom tracks")


if __name__ == "__main__":
    main()
