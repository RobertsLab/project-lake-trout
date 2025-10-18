# PacBio DNA Methylation Workflow Log

## Overview
This directory captures every step of the PacBio DNA methylation exploration, including prompts, commands, scripts, and notes needed to reproduce the workflow end-to-end.

## Logging Guidelines
- Record each prompt or decision in chronological order, noting the context and intent.
- Store any scripts, notebooks, or configuration files generated during the analysis in this folder.
- Include links or relative paths to related assets (e.g., data files, figures) when referenced.

## Environment & Dependencies
- Prefer installing Python packages with [`uv`](https://github.com/astral-sh/uv) to ensure reproducible environments.
- Document any `uv` commands used (e.g., environment creation, package installation) alongside the code that depends on them.
- The local uv project is pinned to Python 3.11 and lives in this folder; a `.venv/` directory is created automatically when dependencies are installed.

## Output Destination
- Write all analysis outputs, intermediate files, and generated artifacts to `../../analyses/04-pacbio/` to keep results organized.

## Environment Setup Log

| Timestamp (UTC) | Command | Purpose |
| --- | --- | --- |
| 2025-10-17 | `uv init --python 3.11` | Scaffolded the dedicated PacBio analysis project and virtual environment managed by uv. |
| 2025-10-17 | `uv add pysam pandas numpy matplotlib seaborn rich` | Added core scientific stack for working with BAM files, tabular data, and visualization. |
| 2025-10-17 | `uv add modbampy` | Enabled parsing of modified-base tags embedded in HiFi BAM records. |
| 2025-10-17 | Created `align_hifi_pbmm2.py` | Batch-align HiFi CCS BAMs to the genome using pbmm2 with configurable paths and CPU counts. |
| 2025-10-17 | `curl https://api.github.com/.../pbmm2/releases/latest` → `wget https://github.com/PacificBiosciences/pbmm2/.../pbmm2` | Downloaded pbmm2 v1.17.0 binary, marked it executable, and symlinked it into `.venv/bin/` for uv-managed runs. |
| 2025-10-17 | Updated `align_hifi_pbmm2.py` | Added flexible pbmm2 discovery (`--pbmm2` override, local `tools/pbmm2` fallback) so the script works even when run outside the project directory. |

> **Note:** Attempted to install `pbcore` from PyPI/GitHub for direct `.pbi` parsing, but its latest release requires `numpy<=1.22.4`, which conflicts with modern `pandas` builds. For now, rely on `pysam` and `modbampy` to access CCS reads and modified-base tags directly from the BAM. If `.pbi` access becomes essential, revisit with a constrained environment or containerized PacBio SMRT Tools installation.

## Alignment Script (`align_hifi_pbmm2.py`)

- Assumes this repository layout (`code/04-pacbio/` for scripts, `data/` for raw HiFi BAMs, `analyses/04-pacbio/` for outputs). The user request referenced `code/04-bio`; the script was placed in `code/04-pacbio/` to match the established naming scheme.
- Requires the PacBio `pbmm2` executable on your `PATH`. Installed here by downloading the v1.17.0 release binary, storing it under `code/04-pacbio/tools/`, and symlinking it into `.venv/bin/` so that `uv run` picks it up automatically.
- Discovers `*.hifi_reads.bam` files (configurable) in the supplied reads directory, aligns each to the reference genome, and writes sorted BAMs plus log files to the output directory.

### Quick start

Run commands *inside* this uv project or provide `--project` so the managed `.venv` (and bundled `pbmm2`) are available. Example from the project directory:

```bash
cd code/04-pacbio
uv run python align_hifi_pbmm2.py \
	--reads-dir ../../data/pacbio-reads \
	--genome ../../data/GCF_016432855.1_SaNama_1.0_genomic.fna.gz \
	--output-dir ../../analyses/04-pacbio/alignments \
	--cpus 32
```

Or from the repository root:

```bash
uv --project code/04-pacbio run python align_hifi_pbmm2.py \
	--reads-dir ../data/pacbio-reads \
	--genome ../data/GCF_016432855.1_SaNama_1.0_genomic.fna.gz \
	--output-dir ../analyses/04-pacbio/alignments \
	--cpus 32
```

Key options:

- `--pattern` (default `*.hifi_reads.bam`) controls which BAMs are aligned.
- `--preset` (default `CCS`) lets you choose pbmm2 presets (`CCS`, `SUBREAD`, etc.).
- `--dry-run` prints planned commands without execution; `--force` overwrites existing outputs.
- `--pbmm2` lets you point to a specific pbmm2 binary; otherwise the script checks PATH and `code/04-pacbio/tools/pbmm2`.
- Pass additional pbmm2 flags by appending them after `--` (e.g., `-- --bam-index`).

Aligned BAMs and pbmm2 logs land in `analyses/04-pacbio/alignments/` by default, keeping results separate from raw data.

## Next Steps
- Add initial prompts and planned analyses.
- Create scaffolding scripts or notebooks for processing PacBio data.
- Outline quality control and validation checkpoints for DNA methylation calls.
