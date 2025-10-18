#!/usr/bin/env python3
"""
Create a comprehensive sequencing data summary report for Lake Trout PacBio data.

This script analyzes local BAM files or generates a template summary report
based on the expected data structure from owl.fish.washington.edu.
"""

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Check for optional dependencies
try:
    import pysam
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False
    print("Note: pysam not available. BAM analysis will be limited.")


def get_file_size(path: Path) -> float:
    """Get file size in GB."""
    return path.stat().st_size / (1024**3)


def analyze_bam_with_samtools(bam_path: Path) -> Dict:
    """
    Analyze BAM file using samtools if available.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        Dictionary with statistics
    """
    stats = {
        "filename": bam_path.name,
        "file_size_gb": get_file_size(bam_path),
    }
    
    try:
        # Check if samtools is available
        result = subprocess.run(
            ["which", "samtools"],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            stats["note"] = "samtools not available for detailed analysis"
            return stats
        
        # Get basic stats with samtools flagstat
        print(f"  Running samtools flagstat on {bam_path.name}...")
        result = subprocess.run(
            ["samtools", "flagstat", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=300
        )
        
        if result.returncode == 0:
            flagstat_output = result.stdout
            # Parse flagstat output
            for line in flagstat_output.split("\n"):
                if "in total" in line:
                    stats["num_reads"] = int(line.split()[0])
                elif "mapped (" in line:
                    stats["mapped_reads"] = int(line.split()[0])
        
        # Get read length statistics
        print(f"  Analyzing read lengths...")
        result = subprocess.run(
            ["samtools", "view", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            lengths = []
            for i, line in enumerate(result.stdout.split("\n")[:1000]):  # Sample first 1000 reads
                if line and not line.startswith("@"):
                    fields = line.split("\t")
                    if len(fields) > 9:
                        lengths.append(len(fields[9]))
            
            if lengths:
                lengths.sort()
                stats["mean_read_length"] = sum(lengths) / len(lengths)
                stats["median_read_length"] = lengths[len(lengths) // 2]
                stats["min_read_length"] = min(lengths)
                stats["max_read_length"] = max(lengths)
                
                # Calculate N50
                total = sum(lengths)
                cumsum = 0
                for length in sorted(lengths, reverse=True):
                    cumsum += length
                    if cumsum >= total / 2:
                        stats["n50"] = length
                        break
    
    except subprocess.TimeoutExpired:
        stats["note"] = "Analysis timed out"
    except Exception as e:
        stats["note"] = f"Error during analysis: {str(e)}"
    
    return stats


def analyze_bam_with_pysam(bam_path: Path, sample_size: int = 10000) -> Dict:
    """
    Analyze BAM file using pysam.
    
    Args:
        bam_path: Path to BAM file
        sample_size: Number of reads to sample for statistics
        
    Returns:
        Dictionary with statistics
    """
    stats = {
        "filename": bam_path.name,
        "file_size_gb": get_file_size(bam_path),
    }
    
    if not HAVE_PYSAM:
        return stats
    
    try:
        print(f"  Opening BAM file with pysam...")
        samfile = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
        
        lengths = []
        qualities = []
        mapped_count = 0
        total_count = 0
        
        for i, read in enumerate(samfile.fetch(until_eof=True)):
            total_count += 1
            
            if i < sample_size:
                if read.query_length:
                    lengths.append(read.query_length)
                
                if read.query_qualities:
                    qualities.append(sum(read.query_qualities) / len(read.query_qualities))
                
                if not read.is_unmapped:
                    mapped_count += 1
            
            if (i + 1) % 100000 == 0:
                print(f"  Processed {i + 1:,} reads...", end="\r")
        
        print(f"  Processed {total_count:,} reads total")
        
        stats["num_reads"] = total_count
        stats["mapped_reads"] = int(mapped_count * total_count / min(sample_size, total_count))
        
        if lengths:
            lengths.sort()
            stats["mean_read_length"] = sum(lengths) / len(lengths)
            stats["median_read_length"] = lengths[len(lengths) // 2]
            stats["min_read_length"] = min(lengths)
            stats["max_read_length"] = max(lengths)
            
            # Calculate N50
            total = sum(lengths)
            cumsum = 0
            for length in sorted(lengths, reverse=True):
                cumsum += length
                if cumsum >= total / 2:
                    stats["n50"] = length
                    break
        
        if qualities:
            stats["mean_quality"] = sum(qualities) / len(qualities)
        
        samfile.close()
        
    except Exception as e:
        stats["note"] = f"Error: {str(e)}"
    
    return stats


def create_summary_report(
    bam_stats: List[Dict],
    files_found: List[str],
    output_path: Path,
    data_url: str
):
    """
    Create a comprehensive markdown summary report.
    
    Args:
        bam_stats: List of BAM file statistics
        files_found: List of data files found
        output_path: Path to save report
        data_url: URL of data source
    """
    print(f"\nCreating summary report: {output_path}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w") as f:
        # Header
        f.write("# PacBio Revio Sequencing Data Summary\n\n")
        f.write("## Lake Trout (_Salvelinus namaycush_) Genomics Project\n\n")
        f.write(f"**Report Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n")
        
        # Overview
        f.write("## Overview\n\n")
        f.write("This report summarizes PacBio Revio HiFi sequencing data generated for ")
        f.write("comparative genomics analysis of lean and siscowet lake trout subspecies.\n\n")
        
        # Data source
        f.write("## Data Source\n\n")
        f.write(f"- **Repository**: {data_url}\n")
        f.write("- **Sequencing Platform**: PacBio Revio\n")
        f.write("- **Technology**: HiFi (High-Fidelity) Circular Consensus Sequencing (CCS)\n")
        f.write("- **Species**: _Salvelinus namaycush_ (Lake Trout)\n")
        f.write("- **Subspecies**: Lean and Siscowet\n\n")
        
        # File inventory
        if files_found:
            f.write("## File Inventory\n\n")
            f.write(f"Total files found: **{len(files_found)}**\n\n")
            
            # Group files by type
            file_types = defaultdict(list)
            for file in files_found:
                filename = file if isinstance(file, str) else file.name
                if ".bam" in filename and not filename.endswith(".pbi"):
                    file_types["HiFi BAM files"].append(filename)
                elif filename.endswith(".bam.pbi"):
                    file_types["PacBio index files (.pbi)"].append(filename)
                elif filename.endswith(".xml"):
                    file_types["Metadata XML files"].append(filename)
                elif filename.endswith((".fastq.gz", ".fasta.gz")):
                    file_types["Sequence files"].append(filename)
                elif filename.endswith((".json", ".csv", ".tsv", ".txt")):
                    file_types["Data/Metadata files"].append(filename)
                else:
                    file_types["Other files"].append(filename)
            
            for file_type, file_list in sorted(file_types.items()):
                f.write(f"### {file_type}\n\n")
                f.write(f"**Count**: {len(file_list)}\n\n")
                for file in sorted(file_list)[:20]:  # Limit to first 20
                    f.write(f"- `{file}`\n")
                if len(file_list) > 20:
                    f.write(f"- ... and {len(file_list) - 20} more files\n")
                f.write("\n")
        
        # Sequencing statistics
        if bam_stats:
            f.write("## Sequencing Statistics\n\n")
            
            # Summary table
            f.write("### Summary Table\n\n")
            f.write("| Sample | File Size (GB) | Total Reads | Mean Length (bp) | N50 (bp) | Mean Quality |\n")
            f.write("|--------|----------------|-------------|------------------|----------|---------------|\n")
            
            total_reads = 0
            total_size = 0
            total_bases = 0
            
            for stats in bam_stats:
                filename = stats.get("filename", "Unknown")
                size = stats.get("file_size_gb", 0)
                reads = stats.get("num_reads", 0)
                mean_len = stats.get("mean_read_length", 0)
                n50 = stats.get("n50", 0)
                quality = stats.get("mean_quality", 0)
                
                total_reads += reads
                total_size += size
                if reads and mean_len:
                    total_bases += reads * mean_len
                
                f.write(f"| {filename} | {size:.2f} | {reads:,} | {mean_len:.0f} | {n50:,} | {quality:.1f} |\n")
            
            f.write(f"| **Total** | **{total_size:.2f}** | **{total_reads:,}** | - | - | - |\n\n")
            
            if total_bases > 0:
                f.write(f"**Total Sequencing Output**: {total_bases / 1e9:.2f} Gb\n\n")
            
            # Detailed statistics per sample
            f.write("### Detailed Read Statistics\n\n")
            for stats in bam_stats:
                f.write(f"#### {stats.get('filename', 'Unknown')}\n\n")
                
                for key, value in stats.items():
                    if key == "filename":
                        continue
                    
                    label = key.replace("_", " ").title()
                    
                    if isinstance(value, float):
                        if "gb" in key.lower():
                            f.write(f"- **{label}**: {value:.2f} GB\n")
                        elif "quality" in key.lower():
                            f.write(f"- **{label}**: {value:.1f}\n")
                        else:
                            f.write(f"- **{label}**: {value:.0f}\n")
                    elif isinstance(value, int):
                        if value > 1000:
                            f.write(f"- **{label}**: {value:,}\n")
                        else:
                            f.write(f"- **{label}**: {value}\n")
                    else:
                        f.write(f"- **{label}**: {value}\n")
                
                f.write("\n")
        
        # Sample information
        f.write("## Sample Information\n\n")
        f.write("### Subspecies\n\n")
        f.write("The sequencing data includes samples from two lake trout subspecies:\n\n")
        f.write("1. **Lean Lake Trout**\n")
        f.write("   - Morphotype: Pelagic/limnetic\n")
        f.write("   - Habitat: Open water\n")
        f.write("   - Characteristics: Streamlined body, smaller fat content\n\n")
        f.write("2. **Siscowet Lake Trout**\n")
        f.write("   - Morphotype: Benthic/profundal\n")
        f.write("   - Habitat: Deep water\n")
        f.write("   - Characteristics: Higher fat content, adapted to deep waters\n\n")
        
        # Technology information
        f.write("## Sequencing Technology\n\n")
        f.write("### PacBio Revio Platform\n\n")
        f.write("The PacBio Revio system is the latest generation of HiFi sequencing technology:\n\n")
        f.write("- **Read Type**: HiFi (High-Fidelity) reads\n")
        f.write("- **Accuracy**: >99.9% (Q30+)\n")
        f.write("- **Read Length**: Typically 10-25 kb, can exceed 30 kb\n")
        f.write("- **Chemistry**: Circular Consensus Sequencing (CCS)\n")
        f.write("- **Applications**: \n")
        f.write("  - De novo genome assembly\n")
        f.write("  - Structural variant detection\n")
        f.write("  - Full-length isoform sequencing\n")
        f.write("  - Epigenetic analysis (5mC, 6mA methylation)\n")
        f.write("  - Haplotype phasing\n\n")
        
        # Potential analyses
        f.write("## Potential Analyses\n\n")
        f.write("This dataset enables multiple types of genomic analyses:\n\n")
        f.write("### 1. Genome Assembly\n")
        f.write("- De novo assembly for each subspecies\n")
        f.write("- Comparative genomics between lean and siscowet\n")
        f.write("- Identification of subspecies-specific genomic features\n\n")
        
        f.write("### 2. Structural Variation Analysis\n")
        f.write("- Detection of large insertions/deletions\n")
        f.write("- Identification of inversions and translocations\n")
        f.write("- Copy number variation analysis\n\n")
        
        f.write("### 3. Isoform Analysis\n")
        f.write("- Full-length transcript sequencing\n")
        f.write("- Alternative splicing patterns\n")
        f.write("- Gene expression differences between subspecies\n\n")
        
        f.write("### 4. Epigenetic Analysis\n")
        f.write("- DNA methylation patterns (5mC)\n")
        f.write("- Comparison of methylation between subspecies\n")
        f.write("- Gene regulation insights\n\n")
        
        # Analysis workflows
        f.write("## Recommended Analysis Workflows\n\n")
        f.write("### Alignment\n")
        f.write("```bash\n")
        f.write("# Align HiFi reads to reference genome using pbmm2\n")
        f.write("pbmm2 align --preset CCS --sort \\\n")
        f.write("  reference.fa \\\n")
        f.write("  input.hifi_reads.bam \\\n")
        f.write("  output.aligned.bam\n")
        f.write("```\n\n")
        
        f.write("### Variant Calling\n")
        f.write("```bash\n")
        f.write("# Call variants using pbsv or DeepVariant\n")
        f.write("pbsv discover aligned.bam variants.svsig.gz\n")
        f.write("pbsv call reference.fa variants.svsig.gz variants.vcf\n")
        f.write("```\n\n")
        
        f.write("### Methylation Analysis\n")
        f.write("```bash\n")
        f.write("# Extract methylation tags using primrose\n")
        f.write("primrose aligned.bam output.bam\n")
        f.write("# Analyze with pb-CpG-tools or custom scripts\n")
        f.write("```\n\n")
        
        # References
        f.write("## References\n\n")
        f.write("### Related Data\n")
        f.write("- **NCBI BioProject**: [PRJNA674328](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA674328)\n")
        f.write("- **Reference Genome**: GCF_016432855.1 (SaNama_1.0)\n\n")
        
        f.write("### Tools and Documentation\n")
        f.write("- [PacBio SMRT Tools](https://www.pacb.com/support/software-downloads/)\n")
        f.write("- [pbmm2 Aligner](https://github.com/PacificBiosciences/pbmm2)\n")
        f.write("- [pbsv Structural Variant Caller](https://github.com/PacificBiosciences/pbsv)\n")
        f.write("- [Primrose Methylation Caller](https://github.com/PacificBiosciences/primrose)\n\n")
        
        # Footer
        f.write("---\n\n")
        f.write("*This report was generated using `create_sequencing_summary.py` from the ")
        f.write("project-lake-trout repository.*\n")
        f.write(f"\n*For questions or issues, please contact the RobertsLab team.*\n")


def scan_directory(directory: Path) -> List[Path]:
    """
    Scan a directory for sequencing data files.
    
    Args:
        directory: Directory to scan
        
    Returns:
        List of file paths found
    """
    files = []
    
    if not directory.exists():
        print(f"Warning: Directory does not exist: {directory}")
        return files
    
    # Common extensions for sequencing data
    extensions = [
        "*.bam", "*.bam.pbi", "*.fastq.gz", "*.fasta.gz",
        "*.xml", "*.json", "*.csv", "*.tsv", "*.txt"
    ]
    
    for ext in extensions:
        files.extend(directory.glob(ext))
        files.extend(directory.glob(f"**/{ext}"))
    
    return files


def main():
    parser = argparse.ArgumentParser(
        description="Create comprehensive sequencing data summary for Lake Trout PacBio data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze local BAM files
  python create_sequencing_summary.py --data-dir /path/to/data --analyze-bams
  
  # Create summary report without analyzing BAMs
  python create_sequencing_summary.py --data-dir /path/to/data
  
  # Generate template report
  python create_sequencing_summary.py --generate-template
        """
    )
    
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Directory containing sequencing data files"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("../../analyses/04-pacbio/sequencing_data_summary.md"),
        help="Output path for summary report (default: ../../analyses/04-pacbio/sequencing_data_summary.md)"
    )
    parser.add_argument(
        "--analyze-bams",
        action="store_true",
        help="Analyze BAM files in detail (requires samtools or pysam)"
    )
    parser.add_argument(
        "--max-bams",
        type=int,
        default=10,
        help="Maximum number of BAM files to analyze (default: 10)"
    )
    parser.add_argument(
        "--data-url",
        default="https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/",
        help="URL of data repository"
    )
    parser.add_argument(
        "--generate-template",
        action="store_true",
        help="Generate a template report without analyzing files"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("PacBio Revio Sequencing Data Summary Generator")
    print("Lake Trout Genomics Project")
    print("=" * 70)
    print()
    
    files_found = []
    bam_stats = []
    
    if args.generate_template:
        print("Generating template report...")
        create_summary_report([], [], args.output, args.data_url)
        print(f"\n✓ Template report created: {args.output}")
        print("\nTo populate with actual data, run with --data-dir option.")
        return
    
    if args.data_dir:
        print(f"Scanning directory: {args.data_dir}")
        files_found = scan_directory(args.data_dir)
        print(f"Found {len(files_found)} files")
        
        if args.analyze_bams and files_found:
            bam_files = [f for f in files_found if f.suffix == ".bam"]
            print(f"\nFound {len(bam_files)} BAM files")
            
            for i, bam_path in enumerate(bam_files[:args.max_bams]):
                print(f"\nAnalyzing BAM {i+1}/{min(len(bam_files), args.max_bams)}: {bam_path.name}")
                
                # Try pysam first, fall back to samtools
                if HAVE_PYSAM:
                    stats = analyze_bam_with_pysam(bam_path)
                else:
                    stats = analyze_bam_with_samtools(bam_path)
                
                bam_stats.append(stats)
    
    # Create the report
    create_summary_report(bam_stats, files_found, args.output, args.data_url)
    
    print("\n" + "=" * 70)
    print(f"✓ Summary report created: {args.output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
