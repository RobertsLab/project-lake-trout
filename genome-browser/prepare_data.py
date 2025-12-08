#!/usr/bin/env python3
"""
prepare_data.py - Data Preparation for Lake Trout Genome Browser

This script prepares data files for the genome browser by:
1. Copying and processing PAV (Presence-Absence Variation) BED files
2. Processing methylation data into browser-compatible bedGraph format
3. Creating summary statistics
4. Generating chromosome size files

Run this script before deploying the genome browser to GitHub Pages.

Usage:
    python prepare_data.py

Output Structure:
    genome-browser/data/
    ├── genome/
    │   ├── GCF_016432855.1_SaNama_1.0_genomic.fa (symlink)
    │   ├── GCF_016432855.1_SaNama_1.0_genomic.fa.fai (symlink)
    │   └── chrom.sizes
    ├── annotations/
    │   └── genes.bed
    ├── pav/
    │   ├── lean_specific.insertions.bed
    │   ├── lean_specific.deletions.bed
    │   ├── siscowet_specific.insertions.bed
    │   ├── siscowet_specific.deletions.bed
    │   ├── shared.insertions.bed
    │   └── shared.deletions.bed
    └── methylation/
        ├── methylation_diff.bedGraph
        ├── lean_mean.bedGraph
        ├── siscowet_mean.bedGraph
        ├── dmrs_hyper_siscowet.bed
        └── dmrs_hypo_siscowet.bed

Author: Generated for project-lake-trout
Date: 2025-12-08
"""

import os
import sys
import shutil
import gzip
from pathlib import Path
import subprocess

# =============================================================================
# Configuration
# =============================================================================

# Paths relative to this script
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_OUTPUT_DIR = SCRIPT_DIR / "data"

# Source directories
PAV_DIR = BASE_DIR / "analyses" / "15-diff-pav"
PAV_RAW_DIR = BASE_DIR / "analyses" / "11-pav"
METH_DIR = BASE_DIR / "analyses" / "10-mCG-call"
DIFF_METH_DIR = BASE_DIR / "analyses" / "14-diff-meth"
GENOME_DIR = BASE_DIR / "data"

# Sample definitions
SAMPLES_LEAN = ["bc2041", "bc2069", "bc2070", "bc2068"]
SAMPLES_SISCOWET = ["bc2071", "bc2073", "bc2072", "bc2096"]
ALL_SAMPLES = SAMPLES_LEAN + SAMPLES_SISCOWET

# =============================================================================
# Helper Functions
# =============================================================================

def ensure_directory(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)

def copy_file(src, dst, strip_header=False):
    """Copy a file, optionally stripping track headers."""
    src = Path(src)
    dst = Path(dst)
    
    if not src.exists():
        print(f"  WARNING: Source file not found: {src}")
        return False
    
    if strip_header:
        with open(src, 'r') as f_in, open(dst, 'w') as f_out:
            for line in f_in:
                if not line.startswith('track'):
                    f_out.write(line)
    else:
        shutil.copy2(src, dst)
    
    print(f"  Copied: {src.name}")
    return True

def create_symlink(src, dst):
    """Create a symbolic link."""
    src = Path(src)
    dst = Path(dst)
    
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    
    if src.exists():
        dst.symlink_to(src.absolute())
        print(f"  Symlinked: {dst.name} -> {src}")
        return True
    else:
        print(f"  WARNING: Source not found for symlink: {src}")
        return False

# =============================================================================
# Step 1: Prepare Genome Reference Files
# =============================================================================

def prepare_genome_files():
    """Prepare genome reference files."""
    print("\n" + "="*60)
    print("Step 1: Preparing genome reference files")
    print("="*60)
    
    genome_output = DATA_OUTPUT_DIR / "genome"
    ensure_directory(genome_output)
    
    # Create symlinks to reference genome (files are too large to copy)
    ref_fa = GENOME_DIR / "GCF_016432855.1_SaNama_1.0_genomic.fa"
    ref_fai = GENOME_DIR / "GCF_016432855.1_SaNama_1.0_genomic.fa.fai"
    
    create_symlink(ref_fa, genome_output / "GCF_016432855.1_SaNama_1.0_genomic.fa")
    create_symlink(ref_fai, genome_output / "GCF_016432855.1_SaNama_1.0_genomic.fa.fai")
    
    # Copy chromosome sizes
    chrom_sizes_src = PAV_RAW_DIR / "genome.chrom.sizes"
    if chrom_sizes_src.exists():
        shutil.copy2(chrom_sizes_src, genome_output / "chrom.sizes")
        print(f"  Copied: chrom.sizes")
    else:
        # Generate from fasta index
        if ref_fai.exists():
            with open(ref_fai, 'r') as f_in, open(genome_output / "chrom.sizes", 'w') as f_out:
                for line in f_in:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        f_out.write(f"{parts[0]}\t{parts[1]}\n")
            print("  Generated: chrom.sizes from FAI")

# =============================================================================
# Step 2: Prepare Gene Annotations
# =============================================================================

def prepare_annotations():
    """Prepare gene annotation files."""
    print("\n" + "="*60)
    print("Step 2: Preparing gene annotations")
    print("="*60)
    
    anno_output = DATA_OUTPUT_DIR / "annotations"
    ensure_directory(anno_output)
    
    genes_src = GENOME_DIR / "20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed"
    
    if genes_src.exists():
        shutil.copy2(genes_src, anno_output / "genes.bed")
        print(f"  Copied: genes.bed")
        
        # Count genes
        with open(genes_src, 'r') as f:
            gene_count = sum(1 for line in f if not line.startswith('#'))
        print(f"  Total genes: {gene_count:,}")
    else:
        print(f"  WARNING: Gene annotation file not found: {genes_src}")

# =============================================================================
# Step 3: Prepare PAV Files
# =============================================================================

def prepare_pav_files():
    """Prepare PAV (Presence-Absence Variation) files."""
    print("\n" + "="*60)
    print("Step 3: Preparing PAV files")
    print("="*60)
    
    pav_output = DATA_OUTPUT_DIR / "pav"
    ensure_directory(pav_output)
    
    # PAV files to copy (from differential PAV analysis)
    pav_files = [
        ("lean_specific.insertions.browser.bed", "lean_specific.insertions.bed"),
        ("lean_specific.deletions.browser.bed", "lean_specific.deletions.bed"),
        ("siscowet_specific.insertions.browser.bed", "siscowet_specific.insertions.bed"),
        ("siscowet_specific.deletions.browser.bed", "siscowet_specific.deletions.bed"),
        ("shared.insertions.browser.bed", "shared.insertions.bed"),
        ("shared.deletions.browser.bed", "shared.deletions.bed"),
    ]
    
    for src_name, dst_name in pav_files:
        src = PAV_DIR / src_name
        dst = pav_output / dst_name
        
        if src.exists():
            # Copy and strip track header (IGV.js doesn't need it)
            copy_file(src, dst, strip_header=True)
            
            # Count variants
            with open(dst, 'r') as f:
                count = sum(1 for line in f if line.strip() and not line.startswith('#'))
            print(f"    Variants: {count:,}")
        else:
            print(f"  WARNING: PAV file not found: {src}")
            print(f"    Run 'python code/15-diff-pav.py' first to generate differential PAV files.")

# =============================================================================
# Step 4: Prepare Methylation Files
# =============================================================================

def prepare_methylation_files():
    """Prepare methylation data files for browser visualization."""
    print("\n" + "="*60)
    print("Step 4: Preparing methylation files")
    print("="*60)
    
    meth_output = DATA_OUTPUT_DIR / "methylation"
    ensure_directory(meth_output)
    
    # Check if differential methylation analysis has been run
    dmr_file = DIFF_METH_DIR / "dmrs.bed"
    sig_dmcs = DIFF_METH_DIR / "significant_dmcs.bed"
    diff_bedgraph = DIFF_METH_DIR / "methylation_diff_siscowet_vs_lean.bedGraph"
    
    # Copy differential methylation bedGraph
    if diff_bedgraph.exists():
        copy_file(diff_bedgraph, meth_output / "methylation_diff.bedGraph", strip_header=True)
    else:
        print("  Differential methylation bedGraph not found.")
        print("    Run 'python code/14-diff-meth.py' first.")
        print("    Creating placeholder methylation files...")
        create_placeholder_methylation(meth_output)
        return
    
    # Copy DMR files if they exist
    hyper_src = DIFF_METH_DIR / "dmrs_hyper_siscowet.bed"
    hypo_src = DIFF_METH_DIR / "dmrs_hypo_siscowet.bed"
    
    if hyper_src.exists():
        copy_file(hyper_src, meth_output / "dmrs_hyper_siscowet.bed")
    if hypo_src.exists():
        copy_file(hypo_src, meth_output / "dmrs_hypo_siscowet.bed")
    
    # Generate mean methylation bedGraphs for each ecotype
    print("\n  Generating ecotype mean methylation files...")
    generate_ecotype_methylation(meth_output)

def create_placeholder_methylation(output_dir):
    """Create placeholder methylation files when analysis hasn't been run."""
    # Create empty placeholder files
    (output_dir / "methylation_diff.bedGraph").touch()
    (output_dir / "dmrs_hyper_siscowet.bed").touch()
    (output_dir / "dmrs_hypo_siscowet.bed").touch()
    (output_dir / "lean_mean.bedGraph").touch()
    (output_dir / "siscowet_mean.bedGraph").touch()
    print("    Created placeholder files (empty)")

def generate_ecotype_methylation(output_dir):
    """Generate mean methylation bedGraph files for each ecotype."""
    
    # Check for bedGraph files
    lean_files = [METH_DIR / f"{s}.cov10.bedGraph" for s in SAMPLES_LEAN]
    sisc_files = [METH_DIR / f"{s}.cov10.bedGraph" for s in SAMPLES_SISCOWET]
    
    lean_exists = [f for f in lean_files if f.exists()]
    sisc_exists = [f for f in sisc_files if f.exists()]
    
    if not lean_exists or not sisc_exists:
        print("    Methylation source files not found - skipping ecotype means")
        (output_dir / "lean_mean.bedGraph").touch()
        (output_dir / "siscowet_mean.bedGraph").touch()
        return
    
    # Process Lean samples - create downsampled mean
    print(f"    Processing Lean methylation ({len(lean_exists)} samples)...")
    create_mean_methylation(lean_exists, output_dir / "lean_mean.bedGraph", sample_rate=100)
    
    # Process Siscowet samples
    print(f"    Processing Siscowet methylation ({len(sisc_exists)} samples)...")
    create_mean_methylation(sisc_exists, output_dir / "siscowet_mean.bedGraph", sample_rate=100)

def create_mean_methylation(input_files, output_file, sample_rate=100):
    """
    Create mean methylation bedGraph from multiple samples.
    Uses downsampling for file size reduction.
    """
    # Read and aggregate methylation values
    meth_data = {}
    
    for filepath in input_files:
        with open(filepath, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                # Downsample to reduce file size
                if line_count % sample_rate != 0:
                    continue
                
                if line.startswith('#') or line.startswith('track'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    try:
                        chrom = parts[0]
                        pos = int(parts[1])
                        meth = float(parts[3])
                        
                        key = (chrom, pos)
                        if key not in meth_data:
                            meth_data[key] = []
                        meth_data[key].append(meth)
                    except (ValueError, IndexError):
                        continue
    
    # Write mean values
    with open(output_file, 'w') as f:
        for (chrom, pos), values in sorted(meth_data.items()):
            mean_meth = sum(values) / len(values)
            f.write(f"{chrom}\t{pos}\t{pos+1}\t{mean_meth:.2f}\n")
    
    print(f"      Wrote {len(meth_data):,} positions to {output_file.name}")

# =============================================================================
# Step 5: Create Summary Statistics File
# =============================================================================

def create_summary():
    """Create a summary statistics JSON file for the browser."""
    print("\n" + "="*60)
    print("Step 5: Creating summary statistics")
    print("="*60)
    
    summary_file = DATA_OUTPUT_DIR / "summary.json"
    
    # Count features in each file
    stats = {
        "genome": {
            "assembly": "GCF_016432855.1",
            "name": "SaNama_1.0",
            "species": "Salvelinus namaycush"
        },
        "samples": {
            "lean": SAMPLES_LEAN,
            "siscowet": SAMPLES_SISCOWET
        },
        "tracks": {}
    }
    
    # Count PAV features
    pav_dir = DATA_OUTPUT_DIR / "pav"
    if pav_dir.exists():
        for bed_file in pav_dir.glob("*.bed"):
            with open(bed_file, 'r') as f:
                count = sum(1 for line in f if line.strip() and not line.startswith('#'))
            stats["tracks"][bed_file.stem] = {"count": count}
    
    # Write summary
    import json
    with open(summary_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"  Created: {summary_file.name}")
    print(f"  Tracks prepared: {len(stats['tracks'])}")

# =============================================================================
# Step 6: Verify Data Integrity
# =============================================================================

def verify_data():
    """Verify that all required data files exist."""
    print("\n" + "="*60)
    print("Step 6: Verifying data integrity")
    print("="*60)
    
    required_files = [
        "annotations/genes.bed",
        "pav/lean_specific.insertions.bed",
        "pav/lean_specific.deletions.bed",
        "pav/siscowet_specific.insertions.bed",
        "pav/siscowet_specific.deletions.bed",
        "methylation/methylation_diff.bedGraph",
    ]
    
    optional_files = [
        "genome/GCF_016432855.1_SaNama_1.0_genomic.fa",
        "genome/GCF_016432855.1_SaNama_1.0_genomic.fa.fai",
        "genome/chrom.sizes",
        "pav/shared.insertions.bed",
        "pav/shared.deletions.bed",
        "methylation/lean_mean.bedGraph",
        "methylation/siscowet_mean.bedGraph",
        "methylation/dmrs_hyper_siscowet.bed",
        "methylation/dmrs_hypo_siscowet.bed",
    ]
    
    missing_required = []
    missing_optional = []
    
    for f in required_files:
        path = DATA_OUTPUT_DIR / f
        if not path.exists() or path.stat().st_size == 0:
            missing_required.append(f)
            print(f"  MISSING (required): {f}")
        else:
            print(f"  OK: {f}")
    
    for f in optional_files:
        path = DATA_OUTPUT_DIR / f
        if not path.exists() or (path.stat().st_size == 0 and not path.is_symlink()):
            missing_optional.append(f)
            print(f"  MISSING (optional): {f}")
        else:
            print(f"  OK: {f}")
    
    print()
    if missing_required:
        print(f"  WARNING: {len(missing_required)} required files missing!")
        print("    The browser may not work correctly.")
    else:
        print("  All required files present.")
    
    if missing_optional:
        print(f"  Note: {len(missing_optional)} optional files missing.")

# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main function to prepare all browser data."""
    print("="*60)
    print("Lake Trout Genome Browser - Data Preparation")
    print("="*60)
    
    # Ensure base output directory exists
    ensure_directory(DATA_OUTPUT_DIR)
    
    # Run preparation steps
    prepare_genome_files()
    prepare_annotations()
    prepare_pav_files()
    prepare_methylation_files()
    create_summary()
    verify_data()
    
    print("\n" + "="*60)
    print("Data Preparation Complete!")
    print("="*60)
    print(f"\nOutput directory: {DATA_OUTPUT_DIR}")
    print("\nNext steps:")
    print("  1. If running locally:")
    print("     cd genome-browser && python -m http.server 8000")
    print("     Open http://localhost:8000 in your browser")
    print()
    print("  2. For GitHub Pages deployment:")
    print("     - Commit and push the genome-browser directory")
    print("     - Enable GitHub Pages in repository settings")
    print("     - Note: Large files (genome FASTA) should be hosted externally")
    print("       or you can use CORS-enabled genome hosting services")
    print()
    print("  3. For external genome hosting (recommended for large files):")
    print("     - Host genome files on a web server with CORS enabled")
    print("     - Update DATA_BASE_URL in js/config.js")

if __name__ == "__main__":
    main()

