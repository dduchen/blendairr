<div align="center">

# blendAIRR

**Build hybrid IgBLAST germline reference databases for custom or non-reference species**

[![Docker](https://img.shields.io/badge/container-ghcr.io-blue?logo=docker)](https://github.com/dduchen/blendairr/pkgs/container/blendAIRR)
[![Build](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

---

## Overview

blendAIRR merges custom species germline sequences (e.g. from a non-reference mouse strain such as MRL) with the closest available IMGT reference species (e.g. C57BL/6 mouse). Novel alleles are jointly clustered with the reference set using [PIgLET](https://bitbucket.org/yaarilab/piglet), which assigns each novel sequence an IMGT-style gene-family designation based on co-clustering with known reference alleles.

The result is a ready-to-run IgBLAST database with all required auxiliary, annotation, and pipeline files — including per-chain `igblastn → MakeDb.py` wrapper scripts and an install script for adding the reference into an existing IGDATA directory.

---

## Quick start

```bash
# Pull
docker pull ghcr.io/dduchen/blendairr:latest

# Run — mount the directory containing your input FASTAs
docker run --rm \
  -v "$(pwd)/data":/data \
  ghcr.io/dduchen/blendairr:latest \
  --species mouse \
  --input_dir /data/MRL \
  --outdir /data/mrl_ref

# Annotate sequences (generated scripts are in the output directory)
bash data/mrl_ref/hybrid_run_heavy.sh sequences.fasta out_prefix
bash data/mrl_ref/hybrid_run_light.sh sequences.fasta out_prefix
```

### Singularity / HPC

```bash
singularity pull docker://ghcr.io/dduchen/blendairr:latest

singularity run blendAIRR_latest.sif \
  --species mouse \
  --input_dir ./MRL \
  --outdir ./mrl_ref
```

---

## Input requirements

blendAIRR expects IMGT-gapped germline FASTA files (dots for gap positions, 312 nt V region) organised as follows. Missing files are filled automatically from the reference species.

```
input_dir/
  heavy/
    IGHV.fasta   ← required; IMGT-gapped V segments
    IGHD.fasta   ← optional
    IGHJ.fasta   ← optional
  light/
    IGKV.fasta   ← optional
    IGKJ.fasta   ← optional
    IGLV.fasta   ← optional
    IGLJ.fasta   ← optional
```

Files may also sit directly in `input_dir/` without subdirectories.

---

## Key arguments

| Argument | Required | Description |
|---|---|---|
| `--species` | ✓ | Closest IMGT reference species (e.g. `mouse`, `human`). Run `--list_species` to see what is available in your IGDATA. |
| `--input_dir` | ✓ | Directory containing your custom germline FASTAs (see layout above). |
| `--outdir` | ✓ | Output directory — created if absent. |
| `--prefix` | | Prefix for all output filenames (default: `hybrid`). |
| `--asc` | | Use PIgLET ASC cluster names (`IGHVFx-Gy*01`) instead of IMGT-style names derived from the closest reference allele. |
| `--skip_blast` | | Skip `makeblastdb` — annotation only. |
| `--list_species` | | List available IMGT species and exit. |

Run `build_hybrid_igblast_ref --help` (or `docker run --rm ghcr.io/dduchen/blendairr --help`) for the full argument reference.

---

## Output structure

All outputs are written under `--outdir`:

```
outdir/
  germlines/
    gapped/              ← hybrid gapped FASTAs — pass to MakeDb.py -r
    ungapped/            ← ungapped FASTAs (intermediate)
  fasta/                 ← edit_imgt_file.pl-processed FASTAs (BLAST input)
  database/              ← IgBLAST BLAST databases (makeblastdb output)
  auxiliary/             ← J-gene aux file + V-gene ndm.imgt annotation
  internal_data/<sp>/    ← internal BLAST DBs for igblastn -organism
  annotations/           ← PIgLET cluster tables, header maps, provenance TSVs
  logs/                  ← per-step log files

  hybrid_run_heavy.sh          ← igblastn + MakeDb.py pipeline for IGH
  hybrid_run_light.sh          ← igblastn + MakeDb.py pipeline for IGK + IGL
  hybrid_install_to_igdata.sh  ← install reference into an existing IGDATA
  hybrid_igblast_cmd.sh        ← low-level igblastn wrapper (all chains)
  hybrid_manifest.tsv          ← inventory of all output files
```

### Generated pipeline scripts

Each generated script contains embedded absolute paths baked in at build time — no configuration needed after `blendAIRR` finishes.

```bash
# Heavy chain (IGH): igblastn fmt7 + AIRR, then MakeDb.py
bash outdir/hybrid_run_heavy.sh sequences.fasta out_prefix [makedb_outdir]

# Light chain (IGK + IGL): same pipeline, separate igblastn calls per locus
bash outdir/hybrid_run_light.sh sequences.fasta out_prefix [makedb_outdir]

# Install into an existing IGDATA (dry-run first!)
bash outdir/hybrid_install_to_igdata.sh ~/share/igblast --dry-run
bash outdir/hybrid_install_to_igdata.sh ~/share/igblast
```

---

## Building locally

```bash
git clone https://github.com/dduchen/blendairr
cd blendAIRR

# Build
make build

# Smoke test
make test

# Run with local image
make run ARGS="--species mouse --input_dir /data/MRL --outdir /data/out"

# Build a Singularity image from the local Docker image
make singularity
```

---

## GitHub Container Registry — image tags

Images are built and pushed automatically by the CI/CD pipeline on every push to `main` and every version tag.

| Git event | Docker tags |
|---|---|
| Push to `main` | `:latest`, `:main`, `:sha-abc1234` |
| Tag `v1.2.3` | `:1.2.3`, `:1.2`, `:1`, `:latest` |
| Pull request | Build + test only — not pushed |

Pin to a specific release for reproducible workflows:
```bash
docker pull ghcr.io/dduchen/blendairr:1.0.0
```

---

## Citation

If you use blendAIRR, please cite:

- **blendAIRR** — [YOUR CITATION HERE]
- **PIgLET** — Omer et al. *(cite PIgLET paper/preprint)*
- **IgBLAST** — Ye et al. 2013, *Nucleic Acids Res.*
- **Immcantation / Change-O** — Gupta et al. 2015, *J. Immunol.*

---

## License

MIT © dduchen
