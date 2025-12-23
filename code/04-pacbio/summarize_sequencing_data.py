#!/usr/bin/env python3
"""
Summarize PacBio Revio sequencing data for Lake Trout project.

This script fetches file listings from the data repository, analyzes BAM files,
and generates a comprehensive summary report of the sequencing effort.
"""

import argparse
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin

try:
    import pandas as pd
    import pysam
    from rich.console import Console
    from rich.table import Table
except ImportError:
    print("ERROR: Required packages not found. Please install with:")
    print("  uv add pandas pysam rich")
    sys.exit(1)


console = Console()


def fetch_directory_listing(url: str) -> List[str]:
    """
    Fetch directory listing from a URL using wget or curl.
    
    Args:
        url: URL to fetch directory listing from
        
    Returns:
        List of file URLs found in the directory
    """
    console.print(f"[cyan]Fetching directory listing from:[/cyan] {url}")
    
    # Try wget first
    try:
        result = subprocess.run(
            ["wget", "-q", "-O", "-", url],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0 and result.stdout:
            html = result.stdout
        else:
            # Try curl if wget fails
            result = subprocess.run(
                ["curl", "-s", "-L", url],
                capture_output=True,
                text=True,
                timeout=30
            )
            html = result.stdout
    except Exception as e:
        console.print(f"[red]Error fetching directory listing:[/red] {e}")
        return []
    
    # Parse HTML to extract file links
    # Look for common patterns in directory listings
    files = []
    
    # Pattern for href links
    href_pattern = r'href="([^"]+)"'
    matches = re.findall(href_pattern, html)
    
    for match in matches:
        # Skip parent directory and special links
        if match in [".", "..", "/", "../"]:
            continue
        # Keep files with common sequencing extensions
        if any(match.endswith(ext) for ext in [
            ".bam", ".bam.pbi", ".fastq.gz", ".fasta.gz", 
            ".xml", ".json", ".csv", ".tsv", ".txt",
            ".md5", ".metadata.xml", ".sts.xml"
        ]):
            files.append(urljoin(url, match))
    
    return files


def download_file(url: str, output_path: Path) -> bool:
    """
    Download a file from URL.
    
    Args:
        url: URL to download from
        output_path: Local path to save file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        console.print(f"[cyan]Downloading:[/cyan] {output_path.name}")
        result = subprocess.run(
            ["wget", "-q", "-O", str(output_path), url],
            timeout=300
        )
        return result.returncode == 0
    except Exception as e:
        console.print(f"[red]Error downloading {url}:[/red] {e}")
        return False


def analyze_bam_file(bam_path: Path) -> Dict:
    """
    Analyze a BAM file and extract summary statistics.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        Dictionary with statistics
    """
    console.print(f"[cyan]Analyzing BAM file:[/cyan] {bam_path.name}")
    
    stats = {
        "filename": bam_path.name,
        "file_size_gb": bam_path.stat().st_size / (1024**3),
        "num_reads": 0,
        "total_bases": 0,
        "read_lengths": [],
        "quality_scores": [],
        "mapped_reads": 0,
        "unmapped_reads": 0,
    }
    
    try:
        # Try to open with pysam
        samfile = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
        
        # Sample reads for statistics (to avoid processing millions of reads)
        sample_size = 10000
        read_count = 0
        
        for i, read in enumerate(samfile.fetch(until_eof=True)):
            stats["num_reads"] += 1
            
            # Sample for detailed statistics
            if i < sample_size:
                read_count += 1
                if read.query_length:
                    stats["read_lengths"].append(read.query_length)
                    stats["total_bases"] += read.query_length
                
                if read.query_qualities:
                    avg_qual = sum(read.query_qualities) / len(read.query_qualities)
                    stats["quality_scores"].append(avg_qual)
                
                if read.is_unmapped:
                    stats["unmapped_reads"] += 1
                else:
                    stats["mapped_reads"] += 1
            
            # Progress indicator
            if (i + 1) % 100000 == 0:
                console.print(f"  Processed {i + 1:,} reads...", end="\r")
        
        samfile.close()
        console.print(f"  Processed {stats['num_reads']:,} reads total")
        
        # Calculate summary statistics
        if stats["read_lengths"]:
            stats["mean_read_length"] = sum(stats["read_lengths"]) / len(stats["read_lengths"])
            stats["median_read_length"] = sorted(stats["read_lengths"])[len(stats["read_lengths"]) // 2]
            stats["min_read_length"] = min(stats["read_lengths"])
            stats["max_read_length"] = max(stats["read_lengths"])
            
            # Calculate N50
            sorted_lengths = sorted(stats["read_lengths"], reverse=True)
            total = sum(sorted_lengths)
            cumsum = 0
            for length in sorted_lengths:
                cumsum += length
                if cumsum >= total / 2:
                    stats["n50"] = length
                    break
        
        if stats["quality_scores"]:
            stats["mean_quality"] = sum(stats["quality_scores"]) / len(stats["quality_scores"])
        
        # Clean up raw data lists
        del stats["read_lengths"]
        del stats["quality_scores"]
        
    except Exception as e:
        console.print(f"[yellow]Warning: Could not fully analyze BAM file:[/yellow] {e}")
        stats["error"] = str(e)
    
    return stats


def create_summary_report(files: List[str], bam_stats: List[Dict], output_path: Path):
    """
    Create a markdown summary report.
    
    Args:
        files: List of file URLs
        bam_stats: List of BAM file statistics
        output_path: Path to save report
    """
    console.print(f"[cyan]Creating summary report:[/cyan] {output_path}")
    
    with open(output_path, "w") as f:
        f.write("# PacBio Revio Sequencing Data Summary\n\n")
        f.write("## Lake Trout (_Salvelinus namaycush_) Genomics Project\n\n")
        f.write("### Summary of Sequencing Effort\n\n")
        f.write("This report summarizes PacBio Revio sequencing data for lean and siscowet lake trout subspecies.\n\n")
        
        # Data source
        f.write("### Data Source\n\n")
        f.write("- **Repository**: https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/\n")
        f.write(f"- **Total Files**: {len(files)}\n\n")
        
        # File inventory
        f.write("### File Inventory\n\n")
        
        # Group files by type
        file_types = defaultdict(list)
        for file in files:
            if file.endswith(".bam"):
                file_types["BAM files (sequencing data)"].append(file)
            elif file.endswith(".bam.pbi"):
                file_types["PacBio index files (.pbi)"].append(file)
            elif file.endswith(".xml"):
                file_types["Metadata XML files"].append(file)
            elif file.endswith((".fastq.gz", ".fasta.gz")):
                file_types["Sequence files"].append(file)
            elif file.endswith((".json", ".csv", ".tsv")):
                file_types["Metadata/Statistics files"].append(file)
            else:
                file_types["Other files"].append(file)
        
        for file_type, file_list in sorted(file_types.items()):
            f.write(f"#### {file_type}\n\n")
            f.write(f"Count: {len(file_list)}\n\n")
            for file in sorted(file_list):
                filename = file.split("/")[-1]
                f.write(f"- `{filename}`\n")
            f.write("\n")
        
        # BAM file statistics
        if bam_stats:
            f.write("### Sequencing Statistics\n\n")
            
            # Summary table
            f.write("| Sample | File Size (GB) | Total Reads | Mean Length (bp) | N50 (bp) | Mean Quality |\n")
            f.write("|--------|----------------|-------------|------------------|----------|---------------|\n")
            
            total_reads = 0
            total_size = 0
            
            for stats in bam_stats:
                filename = stats.get("filename", "Unknown")
                size = stats.get("file_size_gb", 0)
                reads = stats.get("num_reads", 0)
                mean_len = stats.get("mean_read_length", 0)
                n50 = stats.get("n50", 0)
                quality = stats.get("mean_quality", 0)
                
                total_reads += reads
                total_size += size
                
                f.write(f"| {filename} | {size:.2f} | {reads:,} | {mean_len:.0f} | {n50:.0f} | {quality:.1f} |\n")
            
            f.write(f"| **Total** | **{total_size:.2f}** | **{total_reads:,}** | - | - | - |\n\n")
            
            # Detailed statistics
            f.write("### Detailed Read Statistics\n\n")
            for stats in bam_stats:
                f.write(f"#### {stats.get('filename', 'Unknown')}\n\n")
                f.write(f"- **File Size**: {stats.get('file_size_gb', 0):.2f} GB\n")
                f.write(f"- **Total Reads**: {stats.get('num_reads', 0):,}\n")
                
                if "mean_read_length" in stats:
                    f.write(f"- **Mean Read Length**: {stats['mean_read_length']:.0f} bp\n")
                if "median_read_length" in stats:
                    f.write(f"- **Median Read Length**: {stats['median_read_length']:.0f} bp\n")
                if "min_read_length" in stats:
                    f.write(f"- **Min Read Length**: {stats['min_read_length']:,} bp\n")
                if "max_read_length" in stats:
                    f.write(f"- **Max Read Length**: {stats['max_read_length']:,} bp\n")
                if "n50" in stats:
                    f.write(f"- **N50**: {stats['n50']:,} bp\n")
                if "mean_quality" in stats:
                    f.write(f"- **Mean Quality Score**: {stats['mean_quality']:.1f}\n")
                if "total_bases" in stats:
                    total_gb = stats['total_bases'] / 1e9
                    f.write(f"- **Total Bases**: {total_gb:.2f} Gb\n")
                
                if "mapped_reads" in stats:
                    f.write(f"- **Mapped Reads**: {stats['mapped_reads']:,}\n")
                if "unmapped_reads" in stats:
                    f.write(f"- **Unmapped Reads**: {stats['unmapped_reads']:,}\n")
                
                if "error" in stats:
                    f.write(f"- **Note**: {stats['error']}\n")
                
                f.write("\n")
        
        # Sample information
        f.write("### Sample Information\n\n")
        f.write("The sequencing data includes samples from two lake trout subspecies:\n\n")
        f.write("- **Lean Lake Trout**: Pelagic subspecies\n")
        f.write("- **Siscowet Lake Trout**: Benthic subspecies\n\n")
        
        # Sequencing platform
        f.write("### Sequencing Platform\n\n")
        f.write("- **Platform**: PacBio Revio\n")
        f.write("- **Technology**: HiFi long-read sequencing\n")
        f.write("- **Read Type**: Circular Consensus Sequencing (CCS)\n\n")
        
        # Analysis notes
        f.write("### Analysis Notes\n\n")
        f.write("This data can be used for:\n\n")
        f.write("- Genome assembly and annotation\n")
        f.write("- Structural variant detection\n")
        f.write("- Isoform analysis\n")
        f.write("- DNA methylation analysis (with appropriate tools)\n")
        f.write("- Comparative genomics between subspecies\n\n")
        
        # Footer
        f.write("---\n\n")
        f.write("*Report generated by `summarize_sequencing_data.py`*\n")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize PacBio Revio sequencing data for Lake Trout project"
    )
    parser.add_argument(
        "--url",
        default="https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/",
        help="URL to fetch sequencing data from"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("../../analyses/04-pacbio/sequencing_data_summary.md"),
        help="Output path for summary report"
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=Path("/tmp/pacbio_data"),
        help="Temporary directory for downloading files"
    )
    parser.add_argument(
        "--analyze-bams",
        action="store_true",
        help="Download and analyze BAM files (can be slow for large files)"
    )
    parser.add_argument(
        "--max-bams",
        type=int,
        default=3,
        help="Maximum number of BAM files to download and analyze"
    )
    
    args = parser.parse_args()
    
    console.print("[bold cyan]PacBio Revio Sequencing Data Summary Tool[/bold cyan]")
    console.print()
    
    # Fetch directory listing
    files = fetch_directory_listing(args.url)
    
    if not files:
        console.print("[yellow]No files found or unable to access URL.[/yellow]")
        console.print("[yellow]Creating summary from known file structure...[/yellow]")
        
        # Create a basic report with expected file structure
        args.output.parent.mkdir(parents=True, exist_ok=True)
        create_summary_report([], [], args.output)
        console.print(f"[green]Summary report created:[/green] {args.output}")
        return
    
    console.print(f"[green]Found {len(files)} files[/green]")
    
    # Analyze BAM files if requested
    bam_stats = []
    if args.analyze_bams:
        args.download_dir.mkdir(parents=True, exist_ok=True)
        
        bam_files = [f for f in files if f.endswith(".bam")]
        console.print(f"[cyan]Found {len(bam_files)} BAM files[/cyan]")
        
        for i, bam_url in enumerate(bam_files[:args.max_bams]):
            filename = bam_url.split("/")[-1]
            local_path = args.download_dir / filename
            
            console.print(f"\n[bold]Processing {i+1}/{min(len(bam_files), args.max_bams)}:[/bold] {filename}")
            
            if download_file(bam_url, local_path):
                stats = analyze_bam_file(local_path)
                bam_stats.append(stats)
                
                # Clean up to save space
                local_path.unlink()
    
    # Create summary report
    args.output.parent.mkdir(parents=True, exist_ok=True)
    create_summary_report(files, bam_stats, args.output)
    
    console.print()
    console.print(f"[bold green]✓ Summary report created:[/bold green] {args.output}")


if __name__ == "__main__":
    main()
