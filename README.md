<div align="center">

# blendAIRR

**Build hybrid IgBLAST germline reference databases for custom or non-reference species**

[![Docker](https://img.shields.io/badge/container-ghcr.io-blue?logo=docker)](https://github.com/dduchen/blendairr/pkgs/container/blendAIRR)
[![Build](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/dduchen/blendairr/actions/workflows/docker-publish.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

---

## Overview

blendAIRR merges custom germline sequences (e.g. from OGRDB, or from a non-reference strain such as mouse MRL) with the closest available IMGT reference species to produce a ready-to-run IgBLAST germline reference database. A ready-to-run human OGRDB example set ships with the repository under `ogrdb_human/`. Three naming modes are available:

- **`--asc --trig_nc` (recommended)** — custom and reference sequences are jointly clustered with [PIgLET](https://bitbucket.org/yaarilab/piglet) and every allele receives a species-agnostic colon identifier (`IGHV:01:001:001`) under the scheme proposed by the TR/IG Nomenclature Review Committee. An interactive circos plot of the clustering is written for each V locus by default. Opt-in; the bare command runs the default mode.
- **default (input names used directly)** — allele names are kept verbatim; custom sequences are merged with the reference under a priority-aware deduplication (custom/OGRDB take precedence over IMGT). No name mapping table is needed downstream.
- **`--piglet`** — joint PIgLET clustering assigns IMGT-style names to novel alleles (e.g. `IGHV1-24*05`); downstream analysis joins on the emitted name mapping table.

---

## Quick start

The repository ships a ready-to-run human OGRDB example under `ogrdb_human/` (per-locus FASTAs in `heavy/`, `light/`, and `constant/`), so the commands below run as-is after cloning.

```bash
docker pull ghcr.io/dduchen/blendairr:latest

# Recommended — TR/IG NC colon nomenclature on the bundled human OGRDB set
docker run --rm \
  -v "$(pwd)":/data \
  ghcr.io/dduchen/blendairr:latest \
  --species human \
  --input_dir /data/ogrdb_human \
  --outdir /data/ogrdb_trig \
  --asc --trig_nc

# Or the simpler default mode (input names used directly, no clustering)
docker run --rm \
  -v "$(pwd)":/data \
  ghcr.io/dduchen/blendairr:latest \
  --species human \
  --input_dir /data/ogrdb_human \
  --outdir /data/ogrdb_ref

# Install into an IgBLAST share directory (dry-run first)
bash ogrdb_trig/hybrid_install_to_igdata.sh $IGDATA --dry-run
bash ogrdb_trig/hybrid_install_to_igdata.sh $IGDATA

# Annotate sequences with the generated pipeline scripts
export IGDATA=~/share/igblast
bash ogrdb_trig/hybrid_run_heavy.sh sequences.fasta out_prefix
bash ogrdb_trig/hybrid_run_light.sh sequences.fasta out_prefix
```

Singularity / HPC:

```bash
singularity pull docker://ghcr.io/dduchen/blendairr:latest
singularity run blendAIRR_latest.sif \
  --species human --input_dir ./ogrdb_human --outdir ./ogrdb_trig --asc --trig_nc
```

---

## Naming modes

### `--asc --trig_nc` — TR/IG NC colon format (recommended)

`--trig_nc` (implies `--asc`) emits colon-delimited identifiers under the scheme proposed by the TR/IG Nomenclature Review Committee, a species-agnostic alternative to historical IMGT naming.

```
LOCUS:family:gene:allele        e.g.  IGHV1-2*09  ->  IGHV:01:002:008
```

- **Family** is lifted from the closest IMGT family via joint clustering; genuinely novel families are numbered above the highest existing one.
- **Gene** is the ASC allele cluster (PIgLET 95% `cluster_id`), size-ranked within the family, so co-clustering IMGT genes merge into one gene. IMGT gene numbers are not carried through.
- **Allele** is numbered sequentially within each gene cluster.

The strain tag is preserved as a suffix (`IGHV1-20*02_C57BL/6 -> IGHV:01:00N:00M_C57BL/6`). Clustering thresholds default to the NRC values and are exposed as `--family-threshold 75` and `--allele-threshold 95`. To keep colon-format and IMGT-format references from mixing in one `IGDATA`, this mode appends an `_asc` suffix to the prefix (organism `hybrid_asc_human`); override with `--prefix`.

An **interactive circos plot** of the clustering is written per V locus by default, one self-contained HTML per locus under `annotations/`, with family and gene rings, alleles as leaf ticks, and the dendrogram in the centre. It needs the R packages `plotly` and `htmlwidgets`; the step is skipped with a note if they are missing, or with `--skip_plot`.

Name mapping lives in `annotations/<prefix>_final_name_map.tsv` and `<prefix>_full_provenance.tsv`, which round-trip `original_id → normalised_name → trignc_name → final_fasta_name` for every sequence.

> The TR-IG NRC recommendations are still under active discussion; `--trig_nc` implements the colon format and size-ranked clustering as currently proposed, for exploratory / parallel-nomenclature use.

### default — input names used directly

Allele names are preserved exactly as in your input FASTAs and the IMGT reference, so **no mapping table is needed**. Custom/OGRDB sequences take precedence over identical IMGT sequences (except J genes, which are kept so both novel and standard J alleles anchor correctly), IMGT pipe-delimited headers are parsed to clean allele+strain tags (`IGHV1-11*01_C57BL/6`), and exact name+sequence duplicates are collapsed.

```bash
build_hybrid_igblast_ref --species human --input_dir ./ogrdb_human --outdir ./ogrdb_ref
```

### `--piglet` — IMGT-style names for novel alleles

Joint PIgLET clustering assigns each novel allele the next available IMGT-style name in its gene family. Every renamed allele differs from its input id, so downstream analysis must join on `annotations/<prefix>_name_map.tsv` (`source_id → new_allele`) to recover provenance.

### `--species none` — novel species with no reference

For a species with no IMGT or IgBLAST reference, `--species none` runs clustering and TR/IG NC assignment on your input alone. Original ids are kept (illegal FASTA/`makeblastdb` characters sanitized to underscores), family and gene numbers are size-ranked from the clustering, and V/J clustering uses a relaxed 3-prime trim (override with `--v_trim3` / `--j_trim3`). `--ref_dir` and `--igdata` are optional; BLAST databases and a motif-inferred aux file are still built when the IgBLAST tools are available.

```bash
build_hybrid_igblast_ref --species none --input_dir ./novel_sp --outdir ./novel_asc_ref --asc --trig_nc
```

---

## Input layout

IMGT-gapped germline FASTAs (dots for gaps, 312 nt V region), organised as below. Missing loci are filled from the reference species. The bundled `ogrdb_human/` follows this layout.

```
input_dir/
  heavy/   IGHV.fasta (required), IGHD.fasta, IGHJ.fasta
  light/   IGKV.fasta, IGKJ.fasta, IGLV.fasta, IGLJ.fasta
  constant/ IGHC.fasta            (optional custom constant segments)
```

Files may also sit directly in `input_dir/` without subdirectories. blendAIRR reads the per-segment FASTAs; the combined `IGHVDJ.fasta` / `IGKVJ.fasta` / `IGLVJ.fasta` in `ogrdb_human/` are provided for convenience with other tools.

---

## Key arguments

| Argument | Description |
|---|---|
| `--species` | Closest IMGT reference species (e.g. `human`, `mouse`), or `none`. Run `--list_species` to list. |
| `--input_dir` | Directory of custom germline FASTAs (layout above). |
| `--outdir` | Output directory. |
| `--trig_nc` | **Recommended.** Emit TR/IG NC colon identifiers (`IGHV:01:001:001`). Implies `--asc`. |
| `--asc` | Use PIgLET ASC cluster names (`IGHVFx-Gy*01`). Implies `--piglet`. |
| `--piglet` | PIgLET clustering with IMGT-style names for novel alleles (needs a name map downstream). |
| `--prefix` | Output filename prefix (default `hybrid`; `hybrid_asc` under `--asc`/`--trig_nc`). |
| `--family-threshold` / `--allele-threshold` | Clustering identity thresholds, percent (default 75 / 95). |
| `--skip_plot` | Skip the interactive circos plots (they run by default under `--asc`/`--trig_nc`). |
| `--skip_blast` | Annotation only; skip `makeblastdb`. |

Run `build_hybrid_igblast_ref --help` for the full reference, including `--chain`, `--v_trim3`, `--j_trim3`.

---

## Running igblastn

The install script copies the reference into an existing IgBLAST share directory and builds the combined databases; it never overwrites by default (`--force` to overwrite, `--dry-run` to preview). The `-organism` and database names follow the build prefix (`hybrid_<species>`, or `hybrid_asc_<species>` under `--asc`/`--trig_nc`).

```bash
# Heavy chain (IGH)
igblastn \
  -organism       hybrid_asc_human \
  -germline_db_V  $IGDATA/database/imgt_hybrid_asc_human_IGHV \
  -germline_db_D  $IGDATA/database/imgt_hybrid_asc_human_IGHD \
  -germline_db_J  $IGDATA/database/imgt_hybrid_asc_human_IGHJ \
  -auxiliary_data $IGDATA/optional_file/hybrid_asc_human_gl.aux \
  -domain_system  imgt -ig_seqtype Ig \
  -outfmt 19 -query sequences.fasta -out output.airr.tsv

# Light chain (IGK / IGL) — add -num_alignments_D 0 and the constant ig_d database
igblastn \
  -organism       hybrid_asc_human \
  -germline_db_V  $IGDATA/database/imgt_hybrid_asc_human_IGKV \
  -germline_db_J  $IGDATA/database/imgt_hybrid_asc_human_IGKJ \
  -germline_db_D  $IGDATA/database/imgt_hybrid_asc_human_ig_d \
  -num_alignments_D 0 \
  -auxiliary_data $IGDATA/optional_file/hybrid_asc_human_gl.aux \
  -domain_system  imgt -ig_seqtype Ig \
  -outfmt 19 -query light_sequences.fasta -out output.airr.tsv
```

The generated `hybrid_run_heavy.sh` / `hybrid_run_light.sh` run these (plus MakeDb.py) automatically, resolving all paths relative to the output directory.

---

## Citation

If you use blendAIRR, please cite this repository (https://github.com/dduchen/blendairr) along with the tools it builds on:

- **PIgLET** — Ayelet Peres, William D Lees, Oscar L Rodriguez, Noah Y Lee, Pazit Polak, Ronen Hope, Meirav Kedmi, Andrew M Collins, Mats Ohlin, Steven H Kleinstein, Corey T Watson, Gur Yaari. IGHV allele similarity clustering improves genotype inference from adaptive immune receptor repertoire sequencing data. *Nucleic Acids Research*, 51(16):e86, 2023. https://doi.org/10.1093/nar/gkad603
- **IgBLAST** — Ye J, Ma N, Madden TL, Ostell JM. IgBLAST: an immunoglobulin variable domain sequence analysis tool. *Nucleic Acids Research*, 41(W1):W34–W40, 2013.
- **Change-O / Immcantation** — Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics*, 31(20):3356–3358, 2015.

---

## License

MIT © dduchen
