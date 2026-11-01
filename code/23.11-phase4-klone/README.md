# Phase 4 on Klone — scripts

Companion to [`../23.11-phase4-klone-plan.md`](../23.11-phase4-klone-plan.md). Read the plan first,
especially the RNAseq caveat — it determines which evidence arm you should run.

## Order

```bash
source config.sh                        # edit ACCOUNT/PARTITION/LOCAL_HOST first
sbatch 00-stage.sbatch
sbatch 01-repeatmodeler.sbatch          # resubmit until it completes (see below)
sbatch 02-repeatmasker.sbatch
sbatch 03-rnaseq.sbatch                 # ONLY if running ARM=rnaseq
sbatch 04-braker.sbatch                 # ARM=protein by default
sbatch 05-tsebra.sbatch
sbatch 06-specific.sbatch
sbatch 07-function.sbatch
```

## Before the first submission

- `hyakalloc` — confirm account, partition, and core/RAM ceilings; the numbers in these
  scripts are starting points, not verified against a current allocation.
- `ls -l ~/.gm_key` — GeneMark licence. Expired keys fail BRAKER *after* the queue wait.
- Edit `LOCAL_HOST` in `config.sh` to the real hostname of the box holding
  `output/21.1-gene-anchored-synteny/`.

## Step 1 will not finish in one window

`01-repeatmodeler.sbatch` is resubmit-safe: it detects the newest `RM_*` round directory and
resumes with `-recoverDir`. Submit the same script again after each timeout until
`01_repeats/denovo_library.fa` exists. On the preemptible `ckpt` partition this is mandatory.

## `runs.tsv`

Pre-generated from [`../../data/SraRunTable.csv`](../../data/SraRunTable.csv):
`run<TAB>ecotype<TAB>library`, 72 rows, 36 per ecotype. Copy it to `$RNA/` before running
step 3. Note the upstream typo `Sicowet` in the non-parasitized siscowet rows — already handled.

## Arms

`ARM=protein` (default) uses OrthoDB Actinopterygii + SaNama RefSeq proteins only.
`ARM=rnaseq` additionally supplies the PRJNA316738 BAMs. Outputs are kept in separate
directory trees so both can be run and compared.
