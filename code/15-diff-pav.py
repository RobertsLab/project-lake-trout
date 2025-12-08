#!/usr/bin/env python3
"""
15-diff-pav.py - Differential PAV Analysis: Siscowet vs Lean Trout

This script identifies insertions and deletions that are specific to each 
ecotype (Siscowet vs Lean) and creates genome browser-compatible BED files
for visualization.

Sample Classification:
    LEAN: bc2041, bc2069, bc2070, bc2068
    SISCOWET: bc2071, bc2073, bc2072, bc2096

Output files are formatted for genome browser visualization with:
- Color coding (blue for Lean-specific, red for Siscowet-specific)
- Descriptive names
- Track headers for UCSC/IGV compatibility

Author: Generated for project-lake-trout
Date: 2025-12-08
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

# Overlap parameters
MIN_OVERLAP_FRACTION = 0.5  # 50% reciprocal overlap to consider same variant

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

# =============================================================================
# Step 1: Merge variants within each ecotype
# =============================================================================

def merge_ecotype_variants(samples, variant_type, output_prefix):
    """
    Merge variants from all samples of an ecotype into a single file.
    Returns path to merged file.
    """
    available = get_available_samples(samples, variant_type)
    
    if not available:
        print(f"  WARNING: No {variant_type} data available for samples: {samples}")
        return None
    
    print(f"  Available samples for {variant_type}: {', '.join(available)}")
    
    # Concatenate all files
    input_files = [str(INPUT_DIR / f"{s}.{variant_type}.bed") for s in available]
    concat_file = OUTPUT_DIR / f"{output_prefix}.{variant_type}.concat.bed"
    sorted_file = OUTPUT_DIR / f"{output_prefix}.{variant_type}.sorted.bed"
    merged_file = OUTPUT_DIR / f"{output_prefix}.{variant_type}.merged.bed"
    
    # Concatenate (skip headers)
    with open(concat_file, 'w') as outf:
        for infile in input_files:
            with open(infile, 'r') as inf:
                for line in inf:
                    if not line.startswith('#'):
                        outf.write(line)
    
    # Sort
    cmd = f"sort -k1,1 -k2,2n {concat_file} > {sorted_file}"
    run_command(cmd, f"Sorting {variant_type} for {output_prefix}")
    
    # Merge overlapping variants
    # Keep track of sample counts with -c and -o
    cmd = f"""bedtools merge \
        -i {sorted_file} \
        -c 4,5 \
        -o count,sum \
        > {merged_file}"""
    run_command(cmd, f"Merging {variant_type} for {output_prefix}")
    
    # Clean up intermediate files
    concat_file.unlink()
    sorted_file.unlink()
    
    return merged_file

# =============================================================================
# Step 2: Find ecotype-specific variants
# =============================================================================

def find_ecotype_specific(ecotype_file, other_ecotype_file, output_file, 
                          ecotype_name, variant_type, color):
    """
    Find variants present in one ecotype but absent in another.
    Uses bedtools intersect with -v flag to find non-overlapping regions.
    """
    if not file_exists_and_nonempty(ecotype_file):
        print(f"  No {variant_type} data for {ecotype_name}")
        return None
    
    if not file_exists_and_nonempty(other_ecotype_file):
        # If other ecotype has no data, all variants are specific
        cmd = f"cp {ecotype_file} {output_file}"
        run_command(cmd, f"Copying all {variant_type} as {ecotype_name}-specific")
    else:
        # Find variants that don't overlap with other ecotype
        # -f 0.5 -r requires 50% reciprocal overlap to be considered "same"
        cmd = f"""bedtools intersect \
            -a {ecotype_file} \
            -b {other_ecotype_file} \
            -v \
            -f {MIN_OVERLAP_FRACTION} \
            -r \
            > {output_file}"""
        run_command(cmd, f"Finding {ecotype_name}-specific {variant_type}")
    
    return output_file

def find_shared_variants(lean_file, siscowet_file, output_file, variant_type):
    """
    Find variants present in both ecotypes (shared).
    """
    if not file_exists_and_nonempty(lean_file) or not file_exists_and_nonempty(siscowet_file):
        print(f"  Cannot find shared {variant_type} - missing data")
        return None
    
    cmd = f"""bedtools intersect \
        -a {lean_file} \
        -b {siscowet_file} \
        -wa \
        -f {MIN_OVERLAP_FRACTION} \
        -r \
        > {output_file}"""
    run_command(cmd, f"Finding shared {variant_type}")
    
    return output_file

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

def create_bigbed_index(bed_file, chrom_sizes):
    """
    Create bigBed index for large files (optional, requires bedToBigBed).
    """
    try:
        bigbed_file = str(bed_file).replace('.bed', '.bb')
        
        # First need to remove track header for bedToBigBed
        temp_file = str(bed_file) + '.tmp'
        cmd = f"grep -v '^track' {bed_file} > {temp_file}"
        run_command(cmd, "Removing track header")
        
        cmd = f"bedToBigBed {temp_file} {chrom_sizes} {bigbed_file}"
        run_command(cmd, f"Creating bigBed for {bed_file}")
        
        Path(temp_file).unlink()
        return bigbed_file
    except:
        print(f"  Note: bedToBigBed not available, skipping bigBed creation")
        return None

# =============================================================================
# Step 4: Generate summary statistics
# =============================================================================

def generate_summary():
    """Generate summary statistics for differential PAV analysis."""
    summary_data = []
    
    variant_types = ['insertions', 'deletions']
    categories = ['lean_specific', 'siscowet_specific', 'shared']
    
    for vtype in variant_types:
        for category in categories:
            bed_file = OUTPUT_DIR / f"{category}.{vtype}.browser.bed"
            count = count_bed_lines(bed_file)
            
            # Calculate total bp
            total_bp = 0
            if file_exists_and_nonempty(bed_file):
                with open(bed_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or line.startswith('track'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            total_bp += int(parts[2]) - int(parts[1])
            
            summary_data.append({
                'variant_type': vtype,
                'category': category,
                'count': count,
                'total_bp': total_bp
            })
    
    df = pd.DataFrame(summary_data)
    summary_file = OUTPUT_DIR / "diff_pav_summary.csv"
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
        combined_file = OUTPUT_DIR / f"all_{vtype}_by_ecotype.browser.bed"
        
        all_variants = []
        
        # Add Lean-specific (blue)
        lean_file = OUTPUT_DIR / f"lean_specific.{vtype}.browser.bed"
        if file_exists_and_nonempty(lean_file):
            with open(lean_file, 'r') as f:
                for line in f:
                    if not line.startswith('track') and line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            parts[3] = f"LEAN_{parts[3]}"  # Prefix name
                            all_variants.append(parts)
        
        # Add Siscowet-specific (red)
        sisc_file = OUTPUT_DIR / f"siscowet_specific.{vtype}.browser.bed"
        if file_exists_and_nonempty(sisc_file):
            with open(sisc_file, 'r') as f:
                for line in f:
                    if not line.startswith('track') and line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            parts[3] = f"SISC_{parts[3]}"  # Prefix name
                            all_variants.append(parts)
        
        # Add shared (gray)
        shared_file = OUTPUT_DIR / f"shared.{vtype}.browser.bed"
        if file_exists_and_nonempty(shared_file):
            with open(shared_file, 'r') as f:
                for line in f:
                    if not line.startswith('track') and line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            parts[3] = f"BOTH_{parts[3]}"  # Prefix name
                            all_variants.append(parts)
        
        # Sort all variants
        all_variants.sort(key=lambda x: (x[0], int(x[1])))
        
        # Write combined file
        with open(combined_file, 'w') as f:
            vtype_cap = vtype.capitalize()
            f.write(f'track name="All_{vtype_cap}_by_Ecotype" ')
            f.write(f'description="All {vtype} colored by ecotype: Blue=Lean, Red=Siscowet, Gray=Shared" ')
            f.write('visibility=pack itemRgb="On"\n')
            
            for v in all_variants:
                f.write('\t'.join(v) + '\n')
        
        print(f"  Created: {combined_file}")
        print(f"    Total variants: {len(all_variants)}")

# =============================================================================
# Main Analysis
# =============================================================================

def main():
    """Run the differential PAV analysis."""
    print("="*60)
    print("Differential PAV Analysis: Siscowet vs Lean Trout")
    print("="*60)
    
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
    
    print(f"  Lean insertions: {len(lean_ins)} samples")
    print(f"  Lean deletions: {len(lean_del)} samples")
    print(f"  Siscowet insertions: {len(sisc_ins)} samples")
    print(f"  Siscowet deletions: {len(sisc_del)} samples")
    
    # ==========================================================================
    # Process Insertions
    # ==========================================================================
    print("\n" + "="*60)
    print("Step 1: Processing INSERTIONS")
    print("="*60)
    
    # Merge insertions within each ecotype
    lean_ins_merged = merge_ecotype_variants(SAMPLES_LEAN, "insertions", "lean")
    sisc_ins_merged = merge_ecotype_variants(SAMPLES_SISCOWET, "insertions", "siscowet")
    
    # Find ecotype-specific insertions
    if lean_ins_merged and sisc_ins_merged:
        print("\n  Finding ecotype-specific insertions...")
        
        lean_specific_ins = OUTPUT_DIR / "lean_specific.insertions.bed"
        sisc_specific_ins = OUTPUT_DIR / "siscowet_specific.insertions.bed"
        shared_ins = OUTPUT_DIR / "shared.insertions.bed"
        
        find_ecotype_specific(lean_ins_merged, sisc_ins_merged, 
                             lean_specific_ins, "Lean", "insertions", COLOR_LEAN)
        find_ecotype_specific(sisc_ins_merged, lean_ins_merged,
                             sisc_specific_ins, "Siscowet", "insertions", COLOR_SISCOWET)
        find_shared_variants(lean_ins_merged, sisc_ins_merged, 
                            shared_ins, "insertions")
        
        # Create genome browser tracks
        print("\n  Creating genome browser tracks for insertions...")
        
        create_genome_browser_bed(
            lean_specific_ins,
            OUTPUT_DIR / "lean_specific.insertions.browser.bed",
            "Lean_Specific_Insertions",
            "Insertions found in Lean trout but not Siscowet (blue)",
            COLOR_LEAN,
            "insertions"
        )
        
        create_genome_browser_bed(
            sisc_specific_ins,
            OUTPUT_DIR / "siscowet_specific.insertions.browser.bed",
            "Siscowet_Specific_Insertions",
            "Insertions found in Siscowet trout but not Lean (red)",
            COLOR_SISCOWET,
            "insertions"
        )
        
        create_genome_browser_bed(
            shared_ins,
            OUTPUT_DIR / "shared.insertions.browser.bed",
            "Shared_Insertions",
            "Insertions found in both Lean and Siscowet trout (gray)",
            COLOR_SHARED,
            "insertions"
        )
    
    # ==========================================================================
    # Process Deletions
    # ==========================================================================
    print("\n" + "="*60)
    print("Step 2: Processing DELETIONS")
    print("="*60)
    
    # Merge deletions within each ecotype
    lean_del_merged = merge_ecotype_variants(SAMPLES_LEAN, "deletions", "lean")
    sisc_del_merged = merge_ecotype_variants(SAMPLES_SISCOWET, "deletions", "siscowet")
    
    # Find ecotype-specific deletions
    if lean_del_merged and sisc_del_merged:
        print("\n  Finding ecotype-specific deletions...")
        
        lean_specific_del = OUTPUT_DIR / "lean_specific.deletions.bed"
        sisc_specific_del = OUTPUT_DIR / "siscowet_specific.deletions.bed"
        shared_del = OUTPUT_DIR / "shared.deletions.bed"
        
        find_ecotype_specific(lean_del_merged, sisc_del_merged,
                             lean_specific_del, "Lean", "deletions", COLOR_LEAN)
        find_ecotype_specific(sisc_del_merged, lean_del_merged,
                             sisc_specific_del, "Siscowet", "deletions", COLOR_SISCOWET)
        find_shared_variants(lean_del_merged, sisc_del_merged,
                            shared_del, "deletions")
        
        # Create genome browser tracks
        print("\n  Creating genome browser tracks for deletions...")
        
        create_genome_browser_bed(
            lean_specific_del,
            OUTPUT_DIR / "lean_specific.deletions.browser.bed",
            "Lean_Specific_Deletions",
            "Deletions found in Lean trout but not Siscowet (blue)",
            COLOR_LEAN,
            "deletions"
        )
        
        create_genome_browser_bed(
            sisc_specific_del,
            OUTPUT_DIR / "siscowet_specific.deletions.browser.bed",
            "Siscowet_Specific_Deletions",
            "Deletions found in Siscowet trout but not Lean (red)",
            COLOR_SISCOWET,
            "deletions"
        )
        
        create_genome_browser_bed(
            shared_del,
            OUTPUT_DIR / "shared.deletions.browser.bed",
            "Shared_Deletions",
            "Deletions found in both Lean and Siscowet trout (gray)",
            COLOR_SHARED,
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
    print("\nDifferential PAV Summary:")
    print(summary_df.to_string(index=False))
    
    # ==========================================================================
    # Final output
    # ==========================================================================
    print("\n" + "="*60)
    print("Differential PAV Analysis Complete!")
    print("="*60)
    
    print(f"\nOutput files in: {OUTPUT_DIR}")
    print("\nGenome Browser Track Files:")
    print("  Lean-specific (Blue):")
    print("    - lean_specific.insertions.browser.bed")
    print("    - lean_specific.deletions.browser.bed")
    print("  Siscowet-specific (Red):")
    print("    - siscowet_specific.insertions.browser.bed")
    print("    - siscowet_specific.deletions.browser.bed")
    print("  Shared (Gray):")
    print("    - shared.insertions.browser.bed")
    print("    - shared.deletions.browser.bed")
    print("  Combined tracks:")
    print("    - all_insertions_by_ecotype.browser.bed")
    print("    - all_deletions_by_ecotype.browser.bed")
    print("\nTo view in UCSC Genome Browser:")
    print("  1. Go to https://genome.ucsc.edu/cgi-bin/hgCustom")
    print("  2. Choose the appropriate genome assembly")
    print("  3. Upload the .browser.bed files as custom tracks")


if __name__ == "__main__":
    main()

