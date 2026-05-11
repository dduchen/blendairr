<div align="center">

# blendAIRR

**Build hybrid IgBLAST germline reference databases for custom or non-reference species**

[![Docker](https://img.shields.io/badge/container-ghcr.io-blue?logo=docker)](https://github.com/dduchen/blendairr/pkgs/container/blendAIRR)
[![Build](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

---

## Overview

blendAIRR merges custom germline sequences (e.g. from a non-reference mouse strain such as MRL) with the closest available IMGT reference species (e.g. C57BL/6 mouse) to produce a ready-to-run IgBLAST reference database. Two annotation modes are available:

**`--as-is-ids` (recommended)** — Input allele names are used directly. Custom sequences are merged with reference sequences, exact duplicates are removed, and duplicate allele names across strains are disambiguated by appending a strain tag (e.g. `IGHV1-11*01_C57BL/6`). No gene-family clustering is performed. This is the simpler, faster, and more transparent approach — allele names in the output directly correspond to names in your input FASTAs and the IMGT reference.

**Default (PIgLET clustering)** — Novel alleles are jointly clustered with the reference set using [PIgLET](https://bitbucket.org/yaarilab/piglet), which assigns each novel sequence an IMGT-style gene-family name based on co-clustering with known reference alleles. Novel sequences receive new allele designations (e.g. `IGHV1-24*05`). **Downstream analysis must use the allele name mapping table** (`annotations/<prefix>_name_map.tsv`) to reconcile original input names with the renamed alleles assigned by blendAIRR. If PIgLET fails to load, blendAIRR automatically falls back to `--as-is-ids` mode with a warning.

---

## Quick start

```bash
# Pull
docker pull ghcr.io/dduchen/blendairr:latest

# Run (recommended: --as-is-ids)
docker run --rm \
  -v "$(pwd)/data":/data \
  ghcr.io/dduchen/blendairr:latest \
  --species mouse \
  --input_dir /data/MRL \
  --outdir /data/mrl_ref \
  --as-is-ids

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
  --outdir ./mrl_ref \
  --as-is-ids
```

---

## Annotation modes

### `--as-is-ids` — recommended for most users

Input allele names are preserved directly. blendAIRR:

1. Loads your custom FASTAs and the IMGT reference FASTAs for the specified species.
2. Merges them (custom sequences first).
3. Removes exact-duplicate sequences (keeping the custom sequence when both are identical).
4. Renames any remaining alleles that share a name but represent distinct sequences by appending `_1`, `_2` etc.
5. Parses IMGT pipe-delimited headers (e.g. `BK063713|IGHV1-11*01|Mus_musculus_BALB/cJ|F|V-R`) into clean allele+strain tags (e.g. `IGHV1-11*01_BALB/cJ`).

The allele names in igblastn output, MakeDb.py output, and your analysis are all identical to the names in your input files and the IMGT reference — **no mapping table is needed**.

```bash
build_hybrid_igblast_ref \
  --species mouse \
  --input_dir ./MRL \
  --outdir ./mrl_ref \
  --as-is-ids
```

### Default (PIgLET clustering) — for novel allele discovery

PIgLET jointly clusters your custom sequences with the reference and assigns IMGT-style gene-family names to novel alleles. A sequence with no close reference match may be assigned a name like `IGHV1-24*05` (next available allele number in that gene family).

**Important:** every novel allele receives a new name that differs from its original input identifier. All downstream analysis must join on the name mapping table to recover the original source sequence IDs.

```bash
# Build
build_hybrid_igblast_ref \
  --species mouse \
  --input_dir ./MRL \
  --outdir ./mrl_ref

# The name mapping table:
#   annotations/hybrid_name_map.tsv
#   columns: source_id (original input name) -> new_allele (assigned IMGT name)
```

### `--asc`

Uses PIgLET ASC (Allele Sequence Cluster) names (`IGHVFx-Gy*01`) instead of IMGT-style reference-derived names. The same mapping-table caveat as the default mode applies.

---

## Input requirements

blendAIRR expects IMGT-gapped germline FASTA files (dots for gap positions, 312 nt V region) organised as follows. Missing loci are filled automatically from the reference species.

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

**FASTA header formats accepted in `--as-is-ids` mode:**

| Input header | Output name |
|---|---|
| `IGHV7-1*01_S6154` | `IGHV7-1*01_S6154` (kept as-is) |
| `IGHV7-3*04` | `IGHV7-3*04` (kept as-is) |
| `BK063713\|IGHV1-18-28*01\|Mus_musculus_BALB/cJ\|F\|V-R` | `IGHV1-18-28*01_BALB/cJ` |
| `AC090843\|IGHV1-11*01\|Mus_musculus_C57BL/6\|F\|V-REGI` | `IGHV1-11*01_C57BL/6` |

---

## Key arguments

| Argument | Required | Description |
|---|---|---|
| `--species` | ✓ | Closest IMGT reference species (e.g. `mouse`, `human`). Run `--list_species` to see what is available. |
| `--input_dir` | ✓ | Directory containing your custom germline FASTAs (see layout above). |
| `--outdir` | ✓ | Output directory — created if absent. |
| `--as-is-ids` | | **Recommended.** Use input allele names directly; skip PIgLET clustering. No name mapping table required downstream. |
| `--prefix` | | Prefix for all output filenames (default: `hybrid`). |
| `--asc` | | Use PIgLET ASC cluster names (`IGHVFx-Gy*01`) instead of IMGT-style names. Requires PIgLET. |
| `--skip_blast` | | Skip `makeblastdb` — annotation only. |
| `--list_species` | | List available IMGT species and exit. |

Run `build_hybrid_igblast_ref --help` for the full argument reference including advanced options (`--family_threshold`, `--allele_threshold`, `--chain`, `--v_trim3`, `--j_trim3`).

---

## Output structure

All outputs are written under `--outdir`:

```
outdir/
  germlines/
    gapped/              ← hybrid gapped FASTAs — pass to MakeDb.py -r
    ungapped/            ← ungapped FASTAs (intermediate)
  fasta/                 ← edit_imgt_file.pl-processed FASTAs (BLAST input)
  database/              ← IgBLAST BLAST databases (per-locus + combined)
  auxiliary/             ← J-gene aux file + V-gene ndm.imgt annotation
  internal_data/<sp>/    ← internal BLAST DBs for igblastn -organism
  annotations/           ← cluster tables, header maps, name mapping TSVs *
  logs/                  ← per-step log files

  hybrid_run_heavy.sh          ← igblastn + MakeDb.py pipeline for IGH
  hybrid_run_light.sh          ← igblastn + MakeDb.py pipeline for IGK + IGL
  hybrid_install_to_igdata.sh  ← install reference into an existing IGDATA
  hybrid_igblast_cmd.sh        ← low-level igblastn wrapper (all chains)
  hybrid_manifest.tsv          ← inventory of all output files
```

\* `annotations/` is only populated in default (PIgLET) and `--asc` modes. In `--as-is-ids` mode no name mapping table is produced because allele names are unchanged.

### Name mapping table (PIgLET mode only)

`annotations/hybrid_name_map.tsv` contains two key columns:

| Column | Description |
|---|---|
| `source_id` | Original allele name from your input FASTA |
| `new_allele` | IMGT-style name assigned by blendAIRR |

Join this table onto igblastn / MakeDb.py output using the `v_call`, `d_call`, or `j_call` columns to recover provenance:

```R
library(data.table)
name_map <- fread("annotations/hybrid_name_map.tsv")
airr_data <- fread("changeo/gather_gex_heavy_db-pass.tsv")
airr_data[name_map, source_id := i.source_id, on = .(v_call = new_allele)]
```

### Generated pipeline scripts

Each script contains embedded absolute paths baked in at build time.

```bash
# Heavy chain (IGH)
bash outdir/hybrid_run_heavy.sh sequences.fasta out_prefix [makedb_outdir]

# Light chain (IGK + IGL)
bash outdir/hybrid_run_light.sh sequences.fasta out_prefix [makedb_outdir]

# Install into an existing IGDATA (dry-run first!)
bash outdir/hybrid_install_to_igdata.sh ~/share/igblast --dry-run
bash outdir/hybrid_install_to_igdata.sh ~/share/igblast
```

---

## Installing into an existing IGDATA

The install script copies all required files into an existing IgBLAST share directory and builds the combined `ig_v`/`ig_d`/`ig_j` databases via `makeblastdb`. A `--dry-run` flag previews every file operation before committing.

```bash
bash outdir/hybrid_install_to_igdata.sh $IGDATA --dry-run
bash outdir/hybrid_install_to_igdata.sh $IGDATA

# If your germlines directory is in a non-standard location:
bash outdir/hybrid_install_to_igdata.sh $IGDATA \
  --germlines-root /path/to/germlines

# Override the organism name if needed:
bash outdir/hybrid_install_to_igdata.sh $IGDATA \
  --organism-name mrl_mouse
```

After installation, run igblastn with:

```bash
export IGDATA=~/share/igblast
igblastn \
  -organism       hybrid_mouse \
  -germline_db_V  $IGDATA/database/imgt_hybrid_mouse_IGHV \
  -germline_db_D  $IGDATA/database/imgt_hybrid_mouse_IGHD \
  -germline_db_J  $IGDATA/database/imgt_hybrid_mouse_IGHJ \
  -auxiliary_data $IGDATA/optional_file/hybrid_mouse_gl.aux \
  -domain_system  imgt -ig_seqtype Ig \
  -outfmt 19 -query sequences.fasta -out output.airr.tsv
```

---

## Building locally

```bash
git clone https://github.com/dduchen/blendairr
cd blendAIRR

# Vendor the pre-built PIgLET package from your local R library
bash vendor_piglet.sh
git add vendor/piglet-built && git commit -m "vendor piglet"

# Build
make build

# Smoke test
make test

# Run with local image
make run ARGS="--species mouse --input_dir /data/MRL --outdir /data/out --as-is-ids"
```

---

## GitHub Container Registry — image tags

| Git event | Docker tags |
|---|---|
| Push to `main` | `:latest`, `:main`, `:sha-abc1234` |
| Tag `v1.2.3` | `:1.2.3`, `:1.2`, `:1`, `:latest` |
| Pull request | Build + test only — not pushed |

Pin to a specific release for reproducible workflows:
```bash
docker pull ghcr.io/dduchen/blendairr:1.2.0
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
