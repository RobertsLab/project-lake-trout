#!/usr/bin/env python3
"""Batch-align PacBio HiFi CCS reads to a reference genome with pbmm2.

The script discovers BAM inputs that match a configurable pattern (default:
``*.hifi_reads.bam``), invokes `pbmm2 align` for each file, and stores the
aligned, coordinate-sorted BAMs in the analyses area. Modified-base (MM/ML)
information present in the source BAMs is preserved by pbmm2.

Usage examples (run from repository root or this directory):

```
uv run python code/04-pacbio/align_hifi_pbmm2.py \
  --reads-dir data/pacbio-reads \
  --genome data/GCF_016432855.1_SaNama_1.0_genomic.fna.gz \
  --output-dir analyses/04-pacbio/alignments \
  --cpus 46
```

Set `--dry-run` to print planned commands without executing them.

Requirements:
    * `pbmm2` executable available on the PATH
    * Reference genome (FASTA/FASTA.GZ);
      pbmm2 will create and reuse a minimap2 index alongside the FASTA
    * Input BAMs must include PacBio-style tags (MM/ML) if methylation calls are
      desired downstream

The defaults assume this script lives under ``code/04-pacbio/`` inside the
project repository.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
LOGGER = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_READS_DIR = REPO_ROOT / "data"
DEFAULT_GENOME = DEFAULT_READS_DIR / "GCF_016432855.1_SaNama_1.0_genomic.fna.gz"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "analyses" / "04-pacbio" / "alignments"
DEFAULT_PATTERN = "*.hifi_reads.bam"


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align PacBio HiFi BAMs to a reference genome using pbmm2."
    )
    parser.add_argument(
        "--reads-dir",
        type=Path,
        default=DEFAULT_READS_DIR,
        help="Directory containing input HiFi BAMs (default: %(default)s)",
    )
    parser.add_argument(
        "--genome",
        type=Path,
        default=DEFAULT_GENOME,
        help="Reference genome FASTA/FASTQ to align against (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for aligned BAMs and logs (default: %(default)s)",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help="Glob pattern used to find BAMs within --reads-dir (default: %(default)s)",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of CPU threads for pbmm2 (default: %(default)s)",
    )
    parser.add_argument(
        "--preset",
        default="CCS",
        help="pbmm2 preset to use (default: %(default)s)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help="pbmm2 log level (default: %(default)s)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned pbmm2 commands without running them.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output BAMs if they are already present.",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments to append verbatim to each pbmm2 command after `--`.",
    )

    return parser.parse_args(argv)


def find_bam_files(reads_dir: Path, pattern: str) -> list[Path]:
    bam_files = sorted(reads_dir.glob(pattern))
    return [bam for bam in bam_files if bam.is_file()]


def sample_name_from_bam(bam_path: Path) -> str:
    stem = bam_path.stem
    if stem.endswith(".hifi_reads"):
        stem = stem[: -len(".hifi_reads")]
    return stem


def build_pbmm2_command(
    genome: Path,
    bam_path: Path,
    output_bam: Path,
    threads: int,
    preset: str,
    log_level: str,
    log_file: Path,
    extra_args: list[str] | None,
) -> list[str]:
    cmd = [
        "pbmm2",
        "align",
        "--sort",
        "--preset",
        preset,
        "--threads",
        str(threads),
        "--log-level",
        log_level,
        "--log-file",
        str(log_file),
        str(genome),
        str(bam_path),
        str(output_bam),
    ]
    if extra_args:
        cmd.extend(extra_args)
    return cmd


def ensure_pbmm2_available() -> None:
    if not shutil.which("pbmm2"):
        LOGGER.error("pbmm2 executable not found on PATH. Install pbmm2 and retry.")
        sys.exit(1)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)

    ensure_pbmm2_available()

    reads_dir = args.reads_dir.resolve()
    genome = args.genome.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not reads_dir.exists():
        LOGGER.error("Reads directory does not exist: %s", reads_dir)
        return 1

    if not genome.exists():
        LOGGER.error("Reference genome not found: %s", genome)
        return 1

    bam_files = find_bam_files(reads_dir, args.pattern)
    if not bam_files:
        LOGGER.warning(
            "No BAM files matched pattern '%s' in %s", args.pattern, reads_dir
        )
        return 0

    for bam_path in bam_files:
        sample_name = sample_name_from_bam(bam_path)
        output_bam = output_dir / f"{sample_name}.pbmm2.sorted.bam"
        log_file = output_dir / f"{sample_name}.pbmm2.log"

        if output_bam.exists() and not args.force:
            LOGGER.info("Skipping existing alignment: %s", output_bam)
            continue

        cmd = build_pbmm2_command(
            genome=genome,
            bam_path=bam_path,
            output_bam=output_bam,
            threads=args.cpus,
            preset=args.preset,
            log_level=args.log_level,
            log_file=log_file,
            extra_args=args.extra_args,
        )

        LOGGER.info("Running: %s", " ".join(cmd))
        if args.dry_run:
            continue

        completed = subprocess.run(cmd, check=False)
        if completed.returncode != 0:
            LOGGER.error(
                "pbmm2 failed for %s (exit code %s)", bam_path, completed.returncode
            )
            return completed.returncode

    LOGGER.info("All alignments completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
