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

> **Note:** Attempted to install `pbcore` from PyPI/GitHub for direct `.pbi` parsing, but its latest release requires `numpy<=1.22.4`, which conflicts with modern `pandas` builds. For now, rely on `pysam` and `modbampy` to access CCS reads and modified-base tags directly from the BAM. If `.pbi` access becomes essential, revisit with a constrained environment or containerized PacBio SMRT Tools installation.

## Next Steps
- Add initial prompts and planned analyses.
- Create scaffolding scripts or notebooks for processing PacBio data.
- Outline quality control and validation checkpoints for DNA methylation calls.
