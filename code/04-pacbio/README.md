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

## Output Destination
- Write all analysis outputs, intermediate files, and generated artifacts to `../../analyses/04-pacbio/` to keep results organized.

## Next Steps
- Add initial prompts and planned analyses.
- Create scaffolding scripts or notebooks for processing PacBio data.
- Outline quality control and validation checkpoints for DNA methylation calls.
