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

**Default: input names used directly** — Input allele names are used as-is. Custom sequences are merged with reference sequences under a priority-aware deduplication scheme (custom/OGRDB sequences take precedence over IMGT), and IMGT pipe-delimited headers are parsed into clean allele+strain tags (e.g. `IGHV1-11*01_C57BL/6`). No gene-family clustering is performed. This is the simpler, faster, and more transparent approach — allele names in the output directly correspond to names in your input FASTAs and the IMGT reference, so **no name mapping table is needed downstream**. This mode runs by default with the most minimal invocation.

**`--piglet` (opt-in clustering)** — Novel alleles are jointly clustered with the reference set using [PIgLET](https://bitbucket.org/yaarilab/piglet), which assigns each novel sequence an IMGT-style gene-family name based on co-clustering with known reference alleles. Novel sequences receive new allele designations (e.g. `IGHV1-24*05`). **Downstream analysis must use the allele name mapping table** (`annotations/<prefix>_name_map.tsv`) to reconcile original input names with the renamed alleles. If PIgLET fails to load, blendAIRR automatically falls back to the default as-is mode with a warning.

---

## Quick start

```bash
# Pull
docker pull ghcr.io/dduchen/blendairr:latest

# Run (default mode — input names used directly)
docker run --rm \
  -v "$(pwd)/data":/data \
  ghcr.io/dduchen/blendairr:latest \
  --species mouse \
  --input_dir /data/MRL \
  --outdir /data/mrl_ref

# Install into your IgBLAST share directory
bash data/mrl_ref/hybrid_install_to_igdata.sh $IGDATA --dry-run   # preview
bash data/mrl_ref/hybrid_install_to_igdata.sh $IGDATA             # install

# Annotate sequences (generated scripts are in the output directory)
export IGDATA=~/share/igblast
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

## Annotation modes

### Default mode — input names used directly (recommended for most users)

Input allele names are preserved directly. This mode runs by default. blendAIRR:

1. Loads your custom FASTAs and the IMGT reference FASTAs for the specified species.
2. Parses IMGT pipe-delimited headers (e.g. `BK063713|IGHV1-11*01|Mus_musculus_BALB/cJ|F|V-R`) into clean allele+strain tags (e.g. `IGHV1-11*01_BALB/cJ`).
3. Merges custom and reference sequences under a **priority-aware deduplication** scheme:
   - Custom/OGRDB sequences take precedence over IMGT reference sequences.
   - An IMGT sequence identical to a custom/OGRDB sequence is dropped (the custom allele is retained), *except* for J genes.
   - **J genes are exempt** from cross-gene dedup, so both OGRDB novel J alleles and standard IMGT J alleles (e.g. `IGKJ1*01`) are retained — this is required for correct J-gene assignment.
   - Among IMGT sequences, deduplication is applied only *within* the same gene family, preserving cross-gene identical sequences that are biologically meaningful.
4. Disambiguates any remaining identical names by appending `_1`, `_2` etc.
5. OGRDB novel allele names (e.g. `IGKV0-4HMC*00`) are **preserved unchanged**.

The allele names in igblastn output, MakeDb.py output, and your analysis are all identical to the names in your input files and the IMGT reference — **no mapping table is needed**.

```bash
build_hybrid_igblast_ref \
  --species mouse \
  --input_dir ./MRL \
  --outdir ./mrl_ref
```

### `--piglet` — PIgLET clustering for novel allele discovery

PIgLET jointly clusters your custom sequences with the reference and assigns IMGT-style gene-family names to novel alleles. A sequence with no close reference match may be assigned a name like `IGHV1-24*05` (next available allele number in that gene family).

**Important:** every novel allele receives a new name that differs from its original input identifier. All downstream analysis must join on the name mapping table to recover the original source sequence IDs.

```bash
# Build with PIgLET clustering (opt-in)
build_hybrid_igblast_ref \
  --species mouse \
  --input_dir ./MRL \
  --outdir ./mrl_ref \
  --piglet

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
| `--piglet` | | Enable PIgLET clustering to assign IMGT-style names to novel alleles. Requires a name mapping table downstream. Off by default. |
| `--as-is-ids` | | Explicitly request input-names-directly mode. This is the **default**, so the flag is optional (kept for backward compatibility). |
| `--prefix` | | Prefix for all output filenames (default: `hybrid`). |
| `--asc` | | Use PIgLET ASC cluster names (`IGHVFx-Gy*01`) instead of IMGT-style names. Implies `--piglet`. |
| `--skip_blast` | | Skip `makeblastdb` — annotation only. |
| `--list_species` | | List available IMGT species and exit. |

Run `build_hybrid_igblast_ref --help` for the full argument reference including advanced options (`--family_threshold`, `--allele_threshold`, `--chain`, `--v_trim3`, `--j_trim3`).

---

## Output structure

All outputs are written under `--outdir`:

```
outdir/
  germlines/
    gapped/              ← hybrid gapped FASTAs (per-locus + ALL_V/ALL_J) — pass to MakeDb.py -r
    ungapped/            ← ungapped FASTAs (intermediate)
  fasta/                 ← edit_imgt_file.pl-processed FASTAs (BLAST input)
  database/              ← IgBLAST BLAST databases:
                           per-locus (imgt_<org>_IGHV …), combined
                           (imgt_<org>_ig_v/ig_d/ig_j), and constant (ig_c)
  auxiliary/             ← J-gene aux file + V-gene ndm.imgt annotation
  internal_data/<org>/   ← internal BLAST DBs (V spans all loci) for -organism
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

The generated scripts resolve all paths **relative to their own location** (`_SCRIPT_DIR`), so they work whether run from inside the Docker container, copied to the host, or moved to a different directory — as long as the output directory structure stays intact. They also respect a pre-set `$IGDATA` environment variable, falling back to the output directory only when `$IGDATA` is unset.

```bash
# Set IGDATA to your installed share directory (recommended)
export IGDATA=~/share/igblast

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

The install script copies all required files into an existing IgBLAST share directory and builds the combined `ig_v`/`ig_d`/`ig_j`/`ig_c` databases via `makeblastdb`. By default it **never overwrites** existing files. A `--dry-run` flag previews every file operation before committing.

```bash
# Preview (no files written)
bash outdir/hybrid_install_to_igdata.sh $IGDATA --dry-run

# Install (skips files that already exist)
bash outdir/hybrid_install_to_igdata.sh $IGDATA

# Overwrite existing files AND rebuild combined databases (e.g. after a rebuild)
bash outdir/hybrid_install_to_igdata.sh $IGDATA --force

# Preview what --force would overwrite
bash outdir/hybrid_install_to_igdata.sh $IGDATA --dry-run --force

# Point at a different source directory (e.g. host path after a Docker build)
bash outdir/hybrid_install_to_igdata.sh $IGDATA --src-dir ./data/mrl_ref
```

Flags may appear in any order, before or after the IGDATA path. `--force` overwrites existing files (shown as `[OVERWRITTEN]`) and forces a rebuild of the combined `ig_*` BLAST databases from the freshly copied FASTAs.

After installation, run igblastn. **Heavy chain (IGH):**

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

**Light chain (IGK / IGL)** — suppress D-gene search with `-num_alignments_D 0` and supply the constant `ig_d` database so igblastn does not attempt an internal-data D lookup:

```bash
igblastn \
  -organism       hybrid_mouse \
  -germline_db_V  $IGDATA/database/imgt_hybrid_mouse_IGKV \
  -germline_db_J  $IGDATA/database/imgt_hybrid_mouse_IGKJ \
  -germline_db_D  $IGDATA/database/imgt_hybrid_mouse_ig_d \
  -num_alignments_D 0 \
  -auxiliary_data $IGDATA/optional_file/hybrid_mouse_gl.aux \
  -domain_system  imgt -ig_seqtype Ig \
  -outfmt 19 -query light_sequences.fasta -out output.airr.tsv
```

The generated `hybrid_run_light.sh` handles this automatically for both IGK and IGL. Chain-type detection relies on the `internal_data/<organism>/<organism>_V` database spanning **all** loci (IGH+IGK+IGL) — blendAIRR builds this correctly so kappa/lambda queries are classified as VK/VL rather than defaulting to VH.

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

## Troubleshooting

**No J-gene calls on light chains (IGK/IGL).** igblastn determines chain type by aligning the query against the `internal_data/<organism>/<organism>_V` database. If that database contains only heavy-chain V genes, kappa/lambda queries are misclassified as VH and no light J call is made. blendAIRR builds the internal V database from an all-loci `ALL_V.fasta`; if you build the internal database manually, make sure it includes IGHV + IGKV + IGLV.

**`Duplicate seq_ids are found` during `makeblastdb`.** `edit_imgt_file.pl` strips FASTA headers to the gene name, so strain-tagged names (`IGHV1-12*01_C57BL/6`) collapse onto the base name (`IGHV1-12*01`). blendAIRR deduplicates IDs after header processing and before `makeblastdb`; if you script this yourself, add a dedup pass.

**OGRDB novel allele names.** Names such as `IGKV0-4HMC*00` are preserved unchanged. Chain-type detection works correctly through the internal V database, so no renaming is required or performed.

**Install script skips everything / reports wrong paths.** The generated scripts use paths relative to their own location. If you moved only the script without its sibling `database/`, `auxiliary/`, and `internal_data/` directories, use `--src-dir` to point at the intact output directory. Use `--force` to overwrite files from a previous install.

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
