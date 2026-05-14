#!/usr/bin/env bash
# =============================================================================
# blendAIRR — build_hybrid_igblast_ref.sh
#
# DESCRIPTION
#   Builds a hybrid IgBLAST germline reference database by merging custom
#   germline sequences (e.g. a non-reference mouse strain) with the closest
#   IMGT reference species. By default, novel alleles are jointly clustered
#   with the reference using PIgLET to assign IMGT-style gene-family names.
#   Use --as-is-ids to skip clustering and keep input names directly.
#
# QUICK START
#   # Docker (recommended):
#   docker run --rm -v $(pwd)/data:/data ghcr.io/dduchen/blendairr:latest \
#       --species mouse --input_dir /data/MRL --outdir /data/out
#
#   # Direct:
#   bash build_hybrid_igblast_ref.sh \
#       --species mouse --input_dir ./MRL --outdir ./mrl_igblast_ref
#
#   # Then annotate:
#   bash ./mrl_igblast_ref/hybrid_run_heavy.sh sequences.fasta out_prefix
#   bash ./mrl_igblast_ref/hybrid_run_light.sh sequences.fasta out_prefix
#
# REQUIRED
#   -s, --species SPECIES   Closest IMGT reference species (e.g. mouse, human).
#                           Run --list_species to see available options.
#   -i, --input_dir DIR     Directory containing custom germline FASTAs:
#                             DIR/heavy/IGHV.fasta  (required)
#                             DIR/heavy/IGHD.fasta  DIR/heavy/IGHJ.fasta
#                             DIR/light/IGK*.fasta  DIR/light/IGL*.fasta
#                           Files may also sit directly in DIR/. All V files
#                           must be IMGT-gapped (dots = gaps, 312 nt V region).
#                           Missing loci are filled from the reference species.
#   -o, --outdir DIR        Output directory (created if absent).
#
# NAMING / CLUSTERING
#   (default)               PIgLET joint clustering assigns IMGT-style names to
#                           novel alleles. Falls back to --as-is-ids if PIgLET
#                           fails to load, with a warning in the log.
#   --as-is-ids             Skip PIgLET; use input FASTA allele names directly.
#                           Exact-duplicate sequences are removed first; then
#                           duplicate names (distinct sequences) are suffixed
#                           _1, _2 ... (e.g. IGHV1-1*01, IGHV1-1*01_1).
#   --asc                   Use PIgLET ASC cluster names (e.g. IGHVFx-Gy*01)
#                           instead of IMGT-style reference-derived names.
#
# OPTIONAL
#   -g, --igdata DIR        IgBLAST share directory (default: $IGDATA or PATH).
#   -p, --prefix STR        Output file prefix (default: hybrid).
#   -r, --rscript PATH      Path to Rscript binary (default: auto-detect).
#   -e, --edit_imgt PATH    Path to edit_imgt_file.pl (default: auto-detect).
#   --skip_blast            Annotation only; skip makeblastdb.
#   --skip_constant         Skip constant region database build.
#   --family_threshold N    PIgLET family clustering threshold (default: 75).
#   --allele_threshold N    PIgLET allele cluster threshold (default: 95).
#   --v_trim3 N             3-prime V trim for clustering (default: 318).
#   --j_trim3 N             3-prime J trim for clustering (default: 40).
#   --chain heavy|light|all Chain(s) for generated igblastn scripts (default: all).
#   --list_species          List available IMGT species and exit.
#   -h, --help              Show this help and exit.
#
# OUTPUTS  (under --outdir)
#   germlines/gapped/            Hybrid gapped FASTAs (MakeDb.py -r input)
#   database/                    IgBLAST BLAST databases (per-locus + combined)
#   auxiliary/                   J-gene aux file + V-gene ndm.imgt
#   internal_data/<organism>/    Internal BLAST DBs for igblastn -organism
#   annotations/                 Cluster tables, name maps, provenance TSVs
#   logs/                        Per-step logs
#   <prefix>_run_heavy.sh        igblastn + MakeDb.py pipeline — IGH
#   <prefix>_run_light.sh        igblastn + MakeDb.py pipeline — IGK + IGL
#   <prefix>_install_to_igdata.sh  Install reference into an IGDATA directory
#   <prefix>_manifest.tsv        File inventory
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SPECIES=""
INPUT_DIR=""
OUTDIR=""
IGDATA="${IGDATA:-}"
PREFIX="hybrid"
RSCRIPT=""
EDIT_IMGT=""
SKIP_BLAST=false
SKIP_CONSTANT=false
FAM_THRESH=75
ALLELE_THRESH=95
V_TRIM3=318
J_TRIM3=40
LIST_SPECIES=false
USE_ASC=false
AS_IS_IDS=false   # skip PIgLET clustering; use input IDs directly
CHAIN="all"   # heavy | light | all

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
BLUE='\033[0;34m'; NC='\033[0m'; BOLD='\033[1m'

info()  { echo -e "${BLUE}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
err()   { echo -e "${RED}[ERROR]${NC} $*" >&2; }
die()   { err "$*"; exit 1; }
header(){ echo -e "\n${BOLD}=== $* ===${NC}"; }

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
  # Print everything from the DESCRIPTION line to the end of the header block
  sed -n '/^# DESCRIPTION/,/^# ===/p' "$0" | grep '^#' | sed 's/^# \?//'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--species)           SPECIES="$2";        shift 2 ;;
    -i|--input_dir)         INPUT_DIR="$2";      shift 2 ;;
    -o|--outdir)            OUTDIR="$2";         shift 2 ;;
    -g|--igdata)            IGDATA="$2";         shift 2 ;;
    -p|--prefix)            PREFIX="$2";         shift 2 ;;
    -r|--rscript)           RSCRIPT="$2";        shift 2 ;;
    -e|--edit_imgt)         EDIT_IMGT="$2";      shift 2 ;;
    --skip_blast)           SKIP_BLAST=true;     shift ;;
    --skip_constant)        SKIP_CONSTANT=true;  shift ;;
    --family_threshold)     FAM_THRESH="$2";     shift 2 ;;
    --allele_threshold)     ALLELE_THRESH="$2";  shift 2 ;;
    --v_trim3)              V_TRIM3="$2";        shift 2 ;;
    --j_trim3)              J_TRIM3="$2";        shift 2 ;;
    --list_species)         LIST_SPECIES=true;   shift ;;
    --asc)                  USE_ASC=true;        shift ;;
    --as-is-ids)            AS_IS_IDS=true;      shift ;;
    --chain)                CHAIN="$2";          shift 2 ;;
    -h|--help)              usage ;;
    *) die "Unknown argument: $1" ;;
  esac
done

# ---------------------------------------------------------------------------
# Auto-detect IgBLAST
# ---------------------------------------------------------------------------
header "IgBLAST detection"

detect_igblast() {
  # Priority: explicit --igdata > $IGDATA env > common paths > igblastn in PATH
  if [[ -n "$IGDATA" && -d "$IGDATA" ]]; then
    echo "$IGDATA"; return
  fi
  for candidate in \
      /usr/local/share/igblast \
      /usr/share/igblast \
      "$HOME/share/igblast" \
      "$HOME/.local/share/igblast" \
      /opt/igblast; do
    [[ -d "$candidate" ]] && { echo "$candidate"; return; }
  done
  # Fall back to inferring from igblastn binary
  if command -v igblastn &>/dev/null; then
    bin_dir="$(dirname "$(command -v igblastn)")"
    parent="$(dirname "$bin_dir")"
    for sub in share/igblast igblast; do
      [[ -d "$parent/$sub" ]] && { echo "$parent/$sub"; return; }
    done
    echo "$parent"   # best effort
    return
  fi
  echo ""
}

IGDATA="$(detect_igblast)"
if [[ -z "$IGDATA" ]]; then
  warn "IgBLAST share directory not found."
  warn "Set \$IGDATA or pass --igdata.  makeblastdb step will be skipped."
  SKIP_BLAST=true
else
  ok "IgBLAST share: $IGDATA"
fi

# Check igblastn binary
if command -v igblastn &>/dev/null; then
  IGBLASTN_VERSION=$(igblastn -version 2>&1 | head -1)
  ok "igblastn found: $IGBLASTN_VERSION"
else
  warn "igblastn not in PATH — will still build databases but cannot run alignment"
fi

# Check makeblastdb
if command -v makeblastdb &>/dev/null; then
  ok "makeblastdb found: $(makeblastdb -version 2>&1 | head -1)"
else
  warn "makeblastdb not in PATH"; SKIP_BLAST=true
fi

# ---------------------------------------------------------------------------
# List available IMGT species and exit if requested
# ---------------------------------------------------------------------------
if $LIST_SPECIES; then
  header "Available IMGT species in IGDATA"
  if [[ -z "$IGDATA" ]]; then
    die "Cannot list species: IGDATA not found"
  fi
  germline_base="$IGDATA/germlines/imgt"
  if [[ ! -d "$germline_base" ]]; then
    # Try alternate location
    germline_base="$(find "$IGDATA" -maxdepth 3 -type d -name imgt 2>/dev/null | head -1)"
  fi
  if [[ -z "$germline_base" || ! -d "$germline_base" ]]; then
    die "IMGT germline directory not found under $IGDATA"
  fi
  echo ""
  echo "Germline base: $germline_base"
  echo ""
  echo "Available species (directories with vdj/ subdirectory):"
  find "$germline_base" -mindepth 1 -maxdepth 1 -type d | while read -r d; do
    sp="$(basename "$d")"
    if [[ -d "$d/vdj" ]]; then
      loci=$(ls "$d/vdj/"*V*.fasta 2>/dev/null | xargs -I{} basename {} .fasta \
             | sed 's/imgt_[^_]*_//' | sed 's/V$//' | tr '\n' ' ' 2>/dev/null || echo "?")
      printf "  %-20s  loci: %s\n" "$sp" "$loci"
    fi
  done
  echo ""
  exit 0
fi

# ---------------------------------------------------------------------------
# Validate required args
# ---------------------------------------------------------------------------
[[ -z "$SPECIES"   ]] && die "Missing required: --species"
[[ -z "$INPUT_DIR" ]] && die "Missing required: --input_dir"
[[ -z "$OUTDIR"    ]] && die "Missing required: --outdir"
[[ -d "$INPUT_DIR" ]] || die "Input directory not found: $INPUT_DIR"

# ---------------------------------------------------------------------------
# Locate reference IMGT germline directory for chosen species
# ---------------------------------------------------------------------------
header "Locating reference germlines for species: $SPECIES"

find_ref_vdj_dir() {
  local sp="$1"
  # Docker image sets IGBLAST_GERMLINES to the bundled germlines root
  if [[ -n "${IGBLAST_GERMLINES:-}" ]]; then
    for candidate in         "${IGBLAST_GERMLINES}/imgt/${sp}/vdj"         "${IGBLAST_GERMLINES}/imgt/${sp}"         "${IGBLAST_GERMLINES}/${sp}/vdj"; do
      [[ -d "$candidate" ]] && { echo "$candidate"; return; }
    done
  fi
  local base="$IGDATA/germlines/imgt"
  # Try common layouts alongside IGDATA
  for candidate in       "$base/$sp/vdj"       "$base/$sp"       "$IGDATA/germlines/imgt_${sp}"       "$IGDATA/../germlines/imgt/$sp/vdj"; do
    [[ -d "$candidate" ]] && { echo "$candidate"; return; }
  done
  echo ""
}

REF_VDJ_DIR="$(find_ref_vdj_dir "$SPECIES")"
if [[ -z "$REF_VDJ_DIR" ]]; then
  warn "Reference VDJ directory not found for species '$SPECIES'."
  warn "Will use custom sequences only (no joint clustering possible)."
  REF_VDJ_DIR=""
else
  ok "Reference VDJ dir: $REF_VDJ_DIR"
fi

# Locate reference constant region directory
REF_CONST_DIR=""
for candidate in \
    "$IGDATA/germlines/imgt/$SPECIES/constant" \
    "$IGDATA/germlines/imgt/$SPECIES" \
    "$IGDATA/../germlines/imgt/$SPECIES/constant"; do
  [[ -d "$candidate" ]] && { REF_CONST_DIR="$candidate"; break; }
done
[[ -n "$REF_CONST_DIR" ]] && ok "Reference constant dir: $REF_CONST_DIR" \
                           || warn "No constant region directory found for $SPECIES"

# ---------------------------------------------------------------------------
# Locate helper tools
# ---------------------------------------------------------------------------
header "Locating helper tools"

# Rscript — priority: --rscript flag > active conda env > PATH
if [[ -z "$RSCRIPT" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
    RSCRIPT="${CONDA_PREFIX}/bin/Rscript"
    ok "Rscript (conda): $RSCRIPT"
  else
    RSCRIPT="$(command -v Rscript 2>/dev/null || echo "")"
  fi
fi
[[ -z "$RSCRIPT" ]] && die "Rscript not found. Activate a conda env, install R, or use --rscript."
ok "Rscript: $RSCRIPT"

# edit_imgt_file.pl
if [[ -z "$EDIT_IMGT" ]]; then
  for candidate in \
      "$(command -v edit_imgt_file.pl 2>/dev/null)" \
      "$IGDATA/../bin/edit_imgt_file.pl" \
      "$IGDATA/bin/edit_imgt_file.pl" \
      "$(dirname "$(command -v igblastn 2>/dev/null)")/edit_imgt_file.pl" \
      "/usr/local/bin/edit_imgt_file.pl"; do
    [[ -f "$candidate" ]] && { EDIT_IMGT="$candidate"; break; }
  done
fi
if [[ -z "$EDIT_IMGT" || ! -f "$EDIT_IMGT" ]]; then
  warn "edit_imgt_file.pl not found. Attempting to download from igblast GitHub..."
  EDIT_IMGT="$SCRIPT_DIR/edit_imgt_file.pl"
  curl -fsSL \
    "https://raw.githubusercontent.com/psathyrella/igblast/master/bin/edit_imgt_file.pl" \
    -o "$EDIT_IMGT" && chmod +x "$EDIT_IMGT" \
    && ok "Downloaded edit_imgt_file.pl → $EDIT_IMGT" \
    || { warn "Download failed. BLAST db construction may fail."; EDIT_IMGT=""; }
else
  ok "edit_imgt_file.pl: $EDIT_IMGT"
fi

# Check R packages
header "Checking R package dependencies"
$RSCRIPT - <<'RCHECK'
pkgs <- c("optparse","data.table","DECIPHER","piglet","Biostrings","stringr")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(missing)) {
  cat("MISSING R packages:", paste(missing, collapse=", "), "\n")
  quit(status=1)
} else {
  cat("All required R packages found.\n")
}
RCHECK
ok "R packages OK"

# ---------------------------------------------------------------------------
# Create output directories
# ---------------------------------------------------------------------------
header "Creating output directories"
for d in \
    "$OUTDIR/germlines/gapped" \
    "$OUTDIR/germlines/ungapped" \
    "$OUTDIR/fasta" \
    "$OUTDIR/database" \
    "$OUTDIR/internal_data" \
    "$OUTDIR/auxiliary" \
    "$OUTDIR/annotations" \
    "$OUTDIR/logs"; do
  mkdir -p "$d"
done
ok "Output root: $OUTDIR"

# ---------------------------------------------------------------------------
# Pre-Step: Write a debug R script to OUTDIR before invoking R at all.
# This file exists even if the R script crashes immediately.
# ---------------------------------------------------------------------------
header "Pre-Step: Writing bash-level debug R script"

DEBUG_R="$OUTDIR/${PREFIX}_debug_session.R"

# Resolve per-locus file paths the same way the R script will
_cpath() { # find custom file: $1=locus (e.g. IGHV)
  local locus="${1:0:3}" seg="${1:3:1}"  # e.g. IGH / V
  for f in \
      "$INPUT_DIR/heavy/${locus}${seg}.fasta" \
      "$INPUT_DIR/light/${locus}${seg}.fasta" \
      "$INPUT_DIR/${locus}${seg}.fasta" \
      "$INPUT_DIR/imgt_custom_${locus}${seg}.fasta"; do
    [[ -f "$f" ]] && { echo "\"$f\""; return; }
  done
  echo "NULL"
}
_rpath() { # find reference file: $1=locus
  local locus="${1:0:3}" seg="${1:3:1}"
  local sp="$SPECIES"
  for f in \
      "$REF_VDJ_DIR/imgt_${sp}_${locus}${seg}.fasta" \
      "$REF_VDJ_DIR/${locus}${seg}.fasta" \
      "$REF_VDJ_DIR/${sp}/vdj/imgt_${sp}_${locus}${seg}.fasta"; do
    [[ -f "$f" ]] && { echo "\"$f\""; return; }
  done
  echo "NULL"
}

cat > "$DEBUG_R" <<RDEBUG
# ================================================================
# Debug / interactive session script
# Written by build_blendAIRR.sh BEFORE invoking R,
# so this file exists even when the R script crashes.
# Source in R:  source("${DEBUG_R}")
# ================================================================

suppressPackageStartupMessages({
  library(data.table); library(DECIPHER)
  library(piglet);     library(Biostrings)
})

# ---- Parameters (mirrors CLI args passed to R) ----
custom_dir    <- "$INPUT_DIR"
ref_dir       <- "${REF_VDJ_DIR:-$INPUT_DIR}"
species       <- "$SPECIES"
outdir        <- "$OUTDIR"
igdata        <- "$IGDATA"
prefix        <- "$PREFIX"
v_trim3       <- $V_TRIM3
j_trim3       <- $J_TRIM3
fam_thresh    <- $FAM_THRESH
allele_thresh <- $ALLELE_THRESH

# ---- Resolved input file paths ----
custom_path_IGHV <- $(_cpath IGHV)
ref_path_IGHV    <- $(_rpath IGHV)
custom_path_IGHD <- $(_cpath IGHD)
ref_path_IGHD    <- $(_rpath IGHD)
custom_path_IGHJ <- $(_cpath IGHJ)
ref_path_IGHJ    <- $(_rpath IGHJ)
custom_path_IGKV <- $(_cpath IGKV)
ref_path_IGKV    <- $(_rpath IGKV)
custom_path_IGKJ <- $(_cpath IGKJ)
ref_path_IGKJ    <- $(_rpath IGKJ)
custom_path_IGLV <- $(_cpath IGLV)
ref_path_IGLV    <- $(_rpath IGLV)
custom_path_IGLJ <- $(_cpath IGLJ)
ref_path_IGLJ    <- $(_rpath IGLJ)

# ---- Load sequences (mirrors original working code) ----
read_seq <- function(path) {
  if (is.null(path) || identical(path, "NULL") || !file.exists(path))
    return(NULL)
  s <- readDNAStringSet(path)
  names(s) <- gsub(" ", "", names(s))
  s
}

custom_IGHV <- read_seq(custom_path_IGHV)
ref_IGHV    <- read_seq(ref_path_IGHV)
custom_IGKV <- read_seq(custom_path_IGKV)
ref_IGKV    <- read_seq(ref_path_IGKV)
custom_IGLV <- read_seq(custom_path_IGLV)
ref_IGLV    <- read_seq(ref_path_IGLV)
custom_IGHJ <- read_seq(custom_path_IGHJ)
ref_IGHJ    <- read_seq(ref_path_IGHJ)
custom_IGKJ <- read_seq(custom_path_IGKJ)
ref_IGKJ    <- read_seq(ref_path_IGKJ)
custom_IGLJ <- read_seq(custom_path_IGLJ)
ref_IGLJ    <- read_seq(ref_path_IGLJ)
custom_IGHD <- read_seq(custom_path_IGHD)
ref_IGHD    <- read_seq(ref_path_IGHD)

cat("Loaded sequences:\n")
for (nm in c("IGHV","IGKV","IGLV","IGHJ","IGKJ","IGLJ","IGHD")) {
  cust_n <- length(get(paste0("custom_", nm)))
  ref_n  <- length(get(paste0("ref_",    nm)))
  cat(sprintf("  %-6s  custom: %3d   ref: %3d\n", nm, cust_n, ref_n))
}

# ---- PIgLET: the critical named-vector fix ----
# BUG:  as.character(DNAStringSet) drops names -> germ.dist subscript error
# FIX:  setNames(as.character(seqs), names(seqs))
# NOTE: pass GAPPED sequences; PIgLET handles internal trimming by position

run_piglet_debug <- function(seqs_gapped, trim3, label = "") {
  named_vec <- setNames(as.character(seqs_gapped), names(seqs_gapped))
  cat(sprintf("Running PIgLET on %d seqs (%s)\n", length(named_vec), label))
  piglet::inferAlleleClusters(
    germline_set             = named_vec,
    trim_3prime_side         = trim3,
    mask_5prime_side         = 0L,
    family_threshold         = fam_thresh,
    allele_cluster_threshold = allele_thresh
  )
}

# ---- Minimal sanity check ----
if (!is.null(custom_IGHV)) {
  cat("\nTesting PIgLET on custom IGHV alone...\n")
  asc_ighv <- run_piglet_debug(custom_IGHV, trim3 = v_trim3, label = "IGHV")
  print(head(asc_ighv@alleleClusterTable))
}

# ---- Joint set sanity check ----
if (!is.null(custom_IGHV) && !is.null(ref_IGHV)) {
  cat("\nBuilding joint IGHV set (ref first)...\n")
  joint_ighv     <- c(ref_IGHV, custom_IGHV)
  dup_content    <- duplicated(as.character(joint_ighv))
  dup_name       <- duplicated(names(joint_ighv))
  joint_dd       <- joint_ighv[!dup_content & !dup_name]
  cat(sprintf("  ref: %d  custom: %d  joint dedup: %d\n",
              length(ref_IGHV), length(custom_IGHV), length(joint_dd)))
  cat("  Testing PIgLET on joint set...\n")
  asc_joint <- run_piglet_debug(joint_dd, trim3 = v_trim3, label = "IGHV_joint")
  print(head(asc_joint@alleleClusterTable))
}

# ---- Inspect outputs (available after a successful run) ----
ann_dir    <- file.path(outdir, "annotations")
gapped_dir <- file.path(outdir, "germlines", "gapped")

# header_map    <- fread(file.path(ann_dir, paste0(prefix, "_header_normalisation_map.tsv")))
# cluster_annot <- fread(file.path(ann_dir, paste0(prefix, "_allele_cluster_annotation.tsv")))
# provenance    <- fread(file.path(ann_dir, paste0(prefix, "_full_provenance.tsv")))
# hybrid_IGHV   <- readDNAStringSet(file.path(gapped_dir, paste0(prefix, "_IGHV.fasta")))
RDEBUG

info "Debug R script written: $DEBUG_R"


# ---------------------------------------------------------------------------
# Organism identifier — used as the consistent prefix for ALL output files.
# Format: <prefix>_<species>  e.g. "hybrid_mouse"
# This replaces the old pattern of PREFIX for databases and PREFIX_SPECIES for aux.
# ---------------------------------------------------------------------------
ORGANISM="${PREFIX}_${SPECIES}"

# ---------------------------------------------------------------------------
# Step 1: R annotation + germline construction
# ---------------------------------------------------------------------------
header "Step 1: PIgLET joint clustering and germline construction (R)"

R_SCRIPT="${HYBRID_IGBLAST_R_SCRIPT:-$SCRIPT_DIR/R/piglet_annotate_and_build.R}"

# ---------------------------------------------------------------------------
# Run mode: --as-is-ids skips PIgLET and uses input allele names directly.
# Default (PIgLET) mode falls back to as-is automatically if PIgLET fails.
# ---------------------------------------------------------------------------
_run_piglet() {
  [[ -f "$R_SCRIPT" ]] || { warn "R script not found: $R_SCRIPT"; return 1; }
  # Verify PIgLET loads before attempting the full run
  if ! $RSCRIPT -e "library(piglet)" 2>/dev/null; then
    warn "PIgLET R package could not be loaded."
    warn "  Verify install: Rscript -e \"library(piglet)\""
    return 1
  fi
  $RSCRIPT "$R_SCRIPT" \
    --custom_dir        "$INPUT_DIR" \
    --ref_dir           "${REF_VDJ_DIR:-$INPUT_DIR}" \
    --species           "$SPECIES" \
    --outdir            "$OUTDIR" \
    --igdata            "$IGDATA" \
    --prefix            "$PREFIX" \
    --family_threshold  "$FAM_THRESH" \
    --allele_cluster_threshold "$ALLELE_THRESH" \
    --v_trim3prime      "$V_TRIM3" \
    --j_trim3prime      "$J_TRIM3" \
    $($USE_ASC && echo "--use_asc") \
    --organism          "$ORGANISM" \
    2>&1 | tee "$OUTDIR/logs/piglet_annotate.log"
}

_run_as_is() {
  info "Running in --as-is-ids mode: input allele names used directly."
  info "  Sequence deduplication: ON  |  Name deduplication: appends _1, _2 ..."
  [[ -f "$R_SCRIPT" ]] || { warn "R script not found: $R_SCRIPT"; return 1; }
  $RSCRIPT "$R_SCRIPT" \
    --custom_dir        "$INPUT_DIR" \
    --ref_dir           "${REF_VDJ_DIR:-$INPUT_DIR}" \
    --species           "$SPECIES" \
    --outdir            "$OUTDIR" \
    --igdata            "$IGDATA" \
    --prefix            "$PREFIX" \
    --organism          "$ORGANISM" \
    --as_is_ids \
    2>&1 | tee "$OUTDIR/logs/as_is_annotate.log"
}

if $AS_IS_IDS; then
  info "Mode: --as-is-ids (PIgLET clustering skipped)"
  _run_as_is || die "as-is-ids annotation step failed."
else
  info "Mode: PIgLET joint clustering (default)"
  if ! _run_piglet; then
    warn "PIgLET clustering failed. Retrying in --as-is-ids fallback mode."
    warn "  Original log: $OUTDIR/logs/piglet_annotate.log"
    warn "  To skip PIgLET intentionally, use: --as-is-ids"
    _run_as_is || die "as-is-ids fallback also failed. Check input files."
    warn "Completed in as-is-ids fallback mode. Allele names taken from input FASTAs."
  fi
fi

ok "R annotation complete"

# ---------------------------------------------------------------------------
# Step 2: Apply edit_imgt_file.pl and build BLAST databases
# ---------------------------------------------------------------------------
if $SKIP_BLAST; then
  warn "Skipping makeblastdb step (--skip_blast or tools missing)"
else
  header "Step 2: edit_imgt_file.pl + makeblastdb"

  GAPPED_DIR="$OUTDIR/germlines/gapped"
  FASTA_DIR="$OUTDIR/fasta"
  DB_DIR="$OUTDIR/database"

  # Helper: edit + makeblastdb for a single segment
  build_db() {
    local label="$1"         # e.g. IGHV
    local gapped_in="$2"     # gapped FASTA path

    if [[ ! -f "$gapped_in" ]]; then
      warn "  Skipping $label: $gapped_in not found"
      return 0
    fi

    local fasta_out="$FASTA_DIR/imgt_${ORGANISM}_${label}.fasta"
    info "  Processing $label → $fasta_out"

    if [[ -n "$EDIT_IMGT" ]]; then
      perl "$EDIT_IMGT" "$gapped_in" > "$fasta_out"
    else
      awk '/^>/{sub(/ .*/,""); print} !/^>/{gsub(/\./,""); print}' \
        "$gapped_in" > "$fasta_out"
      warn "  edit_imgt_file.pl unavailable; used fallback reformatter for $label"
    fi

    # Use the final imgt_ prefixed name directly — no aliases, no renaming needed.
    # The installed database IS the built database; no path fixups on install.
    makeblastdb \
      -parse_seqids \
      -dbtype nucl \
      -in  "$fasta_out" \
      -out "$db_dir/imgt_${ORGANISM}_${label}" \
      2>&1 | tee -a "$OUTDIR/logs/makeblastdb.log"

    ok "  Built DB: $db_dir/imgt_${ORGANISM}_${label}"
  }

  db_dir="$DB_DIR"

  # V segments — build per-locus DBs (igblastn needs separate -germline_db_V per locus,
  # and also a combined all-locus DB for AssignGenes.py / MakeDb.py)
  for LOCUS in IGHV IGKV IGLV; do
    build_db "$LOCUS" "$GAPPED_DIR/imgt_${ORGANISM}_${LOCUS}.fasta" "$LOCUS"
  done

  # Combined V (for tools that accept a single V db)
  if [[ -f "$GAPPED_DIR/imgt_${ORGANISM}_ALL_V.fasta" ]]; then
    build_db "ALL_V" "$GAPPED_DIR/imgt_${ORGANISM}_ALL_V.fasta" "ALL_V"
  fi

  # D segments
  build_db "IGHD" "$GAPPED_DIR/imgt_${ORGANISM}_IGHD.fasta" "IGHD"

  # J segments — per-locus and combined
  for LOCUS in IGHJ IGKJ IGLJ; do
    build_db "$LOCUS" "$GAPPED_DIR/imgt_${ORGANISM}_${LOCUS}.fasta" "$LOCUS"
  done
  if [[ -f "$GAPPED_DIR/imgt_${ORGANISM}_ALL_J.fasta" ]]; then
    build_db "ALL_J" "$GAPPED_DIR/imgt_${ORGANISM}_ALL_J.fasta" "ALL_J"
  fi

  # Constant regions
  if ! $SKIP_CONSTANT && [[ -n "$REF_CONST_DIR" ]]; then
    info "Building constant region database..."
    CONST_FASTA="$FASTA_DIR/${ORGANISM}_C.fasta"
    : > "$CONST_FASTA"
    for const_file in \
        "$REF_CONST_DIR/imgt_${SPECIES}_IGHC.fasta" \
        "$REF_CONST_DIR/imgt_${SPECIES}_IGKC.fasta" \
        "$REF_CONST_DIR/imgt_${SPECIES}_IGLC.fasta" \
        "$REF_CONST_DIR/"*IGHC.fasta \
        "$REF_CONST_DIR/"*IGKC.fasta \
        "$REF_CONST_DIR/"*IGLC.fasta; do
      [[ -f "$const_file" ]] || continue
      if [[ -n "$EDIT_IMGT" ]]; then
        perl "$EDIT_IMGT" "$const_file" >> "$CONST_FASTA"
      else
        awk '/^>/{sub(/ .*/,""); print} !/^>/{gsub(/\./,""); print}' \
          "$const_file" >> "$CONST_FASTA"
      fi
    done
    # Deduplicate sequences by name
    if command -v seqkit &>/dev/null; then
      seqkit rmdup "$CONST_FASTA" -o "${CONST_FASTA}.dedup" && \
        mv "${CONST_FASTA}.dedup" "$CONST_FASTA"
    fi
    if [[ -s "$CONST_FASTA" ]]; then
      makeblastdb -parse_seqids -dbtype nucl \
        -in  "$CONST_FASTA" \
        -out "$db_dir/${ORGANISM}_C" \
        2>&1 | tee -a "$OUTDIR/logs/makeblastdb.log"
      ok "Built constant region DB: $db_dir/${ORGANISM}_C"
    else
      warn "Constant region FASTA empty; skipping DB build"
    fi
  fi

  ok "makeblastdb step complete"

  # ---------------------------------------------------------------------------
  # Step 3: Build internal_data/<species>/ for igblastn
  # igblastn requires: $IGDATA/internal_data/<organism>/<organism>_V (blast db)
  #                    $IGDATA/internal_data/<organism>/<organism>.ndm.imgt
  # The internal V/D/J dbs are SEPARATE from the search dbs -- they are used
  # only for alignment coordinate mapping (one allele per gene is enough).
  # We must build them fresh with makeblastdb from processed FASTA files.
  # ---------------------------------------------------------------------------
  header "Step 3: Build internal_data/ databases"
  INTERNAL="$OUTDIR/internal_data/$ORGANISM"
  mkdir -p "$INTERNAL"

  # Helper: process a gapped FASTA through edit_imgt_file.pl (or fallback),
  # then build a BLAST db named <dest_base> inside $INTERNAL.
  # Files are named <ORGANISM>_V etc. so igblastn -organism $ORGANISM finds them.
  # Write python3 translation helper to a fixed path
  TRANSLATE_PY="$OUTDIR/logs/_translate.py"
  cat > "$TRANSLATE_PY" << 'TRANSLATE_PY_EOF'
import sys, textwrap
CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
def translate(nt):
    nt = nt.upper().replace('-','').replace('.','').replace(' ','')
    aa = [CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3)]
    return ''.join(aa).rstrip('*') or 'X'
with open(sys.argv[1]) as fin, open(sys.argv[2],'w') as fout:
    name,seq = None,[]
    for line in fin:
        line = line.rstrip()
        if line.startswith('>'):
            if name and seq:
                p = translate(''.join(seq))
                if p != 'X': fout.write(f'>{name}\n'+''.join(c+'\n' for c in textwrap.wrap(p,60)))
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line)
    if name and seq:
        p = translate(''.join(seq))
        if p != 'X': fout.write(f'>{name}\n'+''.join(c+'\n' for c in textwrap.wrap(p,60)))
TRANSLATE_PY_EOF

  build_internal_db() {
    local label="$1"
    local src_fasta="$2"
    [[ -f "$src_fasta" ]] || { warn "  internal_data $label: $src_fasta not found"; return 0; }

    local tmp_fasta="${INTERNAL}/${ORGANISM}_${label}_raw.fasta"
    local dest_base="${INTERNAL}/${ORGANISM}_${label}"
    local tmp_prot="${INTERNAL}/${ORGANISM}_${label}_prot.fasta"

    # Process gapped FASTA: strip gaps and clean headers
    if [[ -n "$EDIT_IMGT" ]]; then
      perl "$EDIT_IMGT" "$src_fasta" > "$tmp_fasta"
    else
      awk '/^>/{sub(/ .*/,""); print} !/^>/{gsub(/\./, ""); print}' \
        "$src_fasta" > "$tmp_fasta"
    fi

    # Nucleotide BLAST database (n* files)
    makeblastdb -parse_seqids -dbtype nucl \
      -in "$tmp_fasta" -out "$dest_base" \
      >> "$OUTDIR/logs/makeblastdb.log" 2>&1 \
      && ok "  internal_data ${ORGANISM}_${label} (nucl): built" \
      || warn "  internal_data ${ORGANISM}_${label} (nucl): failed"

    # Protein BLAST database (p* files) — igblastn requires both alongside nucl
    python3 "$TRANSLATE_PY" "$tmp_fasta" "$tmp_prot" \
      >> "$OUTDIR/logs/makeblastdb.log" 2>&1

    if [[ -f "$tmp_prot" && -s "$tmp_prot" ]]; then
      makeblastdb -parse_seqids -dbtype prot \
        -in "$tmp_prot" -out "$dest_base" \
        >> "$OUTDIR/logs/makeblastdb.log" 2>&1 \
        && ok "  internal_data ${ORGANISM}_${label} (prot): built" \
        || warn "  internal_data ${ORGANISM}_${label} (prot): makeblastdb failed"
      rm -f "$tmp_prot"
    else
      warn "  internal_data ${ORGANISM}_${label} (prot): translation failed; check $OUTDIR/logs/makeblastdb.log"
    fi

    rm -f "$tmp_fasta"
  }

  # Build V (combined all-locus), D, J internal databases
  # Prefer combined ALL_V/ALL_J; fall back to IGH locus only
  for v_src in \
      "$GAPPED_DIR/imgt_${ORGANISM}_ALL_V.fasta" \
      "$GAPPED_DIR/imgt_${ORGANISM}_IGHV.fasta"; do
    [[ -f "$v_src" ]] && { build_internal_db "V" "$v_src"; break; }
  done

  build_internal_db "D" "$GAPPED_DIR/imgt_${ORGANISM}_IGHD.fasta"

  for j_src in \
      "$GAPPED_DIR/imgt_${ORGANISM}_ALL_J.fasta" \
      "$GAPPED_DIR/imgt_${ORGANISM}_IGHJ.fasta"; do
    [[ -f "$j_src" ]] && { build_internal_db "J" "$j_src"; break; }
  done

  # Copy ndm.imgt into internal_data (built by the R step)
  NDM_SRC="$OUTDIR/auxiliary/${ORGANISM}.ndm.imgt"
  if [[ -f "$NDM_SRC" ]]; then
    cp "$NDM_SRC" "${INTERNAL}/${ORGANISM}.ndm.imgt"
    ok "  internal_data ${ORGANISM}.ndm.imgt: copied"
  else
    warn "  ndm.imgt not found at $NDM_SRC (R step may not have generated it yet)"
  fi

  # Create optional_file/ under OUTDIR for aux lookup when IGDATA=OUTDIR
  mkdir -p "$OUTDIR/optional_file"
  AUX_SRC_STEP3="$OUTDIR/auxiliary/${ORGANISM}_gl.aux"
  [[ -f "$AUX_SRC_STEP3" ]] && \
    cp "$AUX_SRC_STEP3" "$OUTDIR/optional_file/${ORGANISM}_gl.aux" 2>/dev/null || true

  ok "internal_data populated: $INTERNAL  (organism: $ORGANISM)"
  info "  Files in $INTERNAL:"
  ls -1 "$INTERNAL/" 2>/dev/null | sed 's/^/    /' || true


fi  # end if $SKIP_BLAST (Step 2 + Step 3 block)
# ---------------------------------------------------------------------------
# Step 4: Copy auxiliary file to IGDATA optional_file (if writable)
# ---------------------------------------------------------------------------
header "Step 4: Auxiliary file"

AUX_SRC="$OUTDIR/auxiliary/${ORGANISM}_gl.aux"
if [[ -f "$AUX_SRC" ]]; then
  ok "Auxiliary file: $AUX_SRC"
  AUX_OPT_DIR="$IGDATA/optional_file"
  if [[ -n "$IGDATA" && -d "$AUX_OPT_DIR" && -w "$AUX_OPT_DIR" ]]; then
    cp "$AUX_SRC" "$AUX_OPT_DIR/${ORGANISM}_gl.aux"
    ok "Copied to IGDATA: $AUX_OPT_DIR/${ORGANISM}_gl.aux"
  else
    info "IGDATA optional_file not writable; pass aux path explicitly in igblastn call"
  fi
else
  warn "Auxiliary file not found (R step may have failed)"
fi

# ---------------------------------------------------------------------------
# Step 5: Generate ready-to-use igblastn and MakeDb.py command templates
# ---------------------------------------------------------------------------
header "Step 5: Generating command templates"

GAPPED_DIR="$OUTDIR/germlines/gapped"
DB_DIR="$OUTDIR/database"
AUX_FILE="$OUTDIR/auxiliary/${ORGANISM}_gl.aux"
IGBLAST_CMD="$OUTDIR/${PREFIX}_igblast_cmd.sh"

# Resolve OUTDIR to an absolute path for the generated script
_ABS_OUTDIR="$(cd "$OUTDIR" && pwd)"
_ABS_DB="${_ABS_OUTDIR}/database"
_ABS_GAPPED="${_ABS_OUTDIR}/germlines/gapped"
_ABS_AUX="${_ABS_OUTDIR}/auxiliary/${ORGANISM}_gl.aux"
_ABS_NDM="${_ABS_OUTDIR}/auxiliary/${ORGANISM}.ndm.imgt"

cat > "$IGBLAST_CMD" <<'IGBLAST_EOF'
#!/usr/bin/env bash
# Auto-generated igblastn command script
# Usage: bash <this_script> <query.fasta> [out_prefix] [chain: heavy|light|all]

set -euo pipefail

QUERY_FASTA="${1:?Usage: $0 <query.fasta> [out_prefix] [chain: heavy|light|all]}"
OUT_PREFIX="${2:-igblast_out}"
CHAIN="${3:-all}"    # heavy | light | all
IGBLAST_EOF

# Bake in resolved paths and species at write time
cat >> "$IGBLAST_CMD" <<IGBLAST_EOF2

# ---- Resolved paths (baked in at build time) ----
DB_PREFIX="${_ABS_DB}/imgt_${ORGANISM}"
AUX_FILE="${_ABS_AUX}"
NDM_FILE="${_ABS_NDM}"
ORGANISM="${SPECIES}"
export IGDATA="${_ABS_OUTDIR}"

# ---- Optional args (set at runtime from baked-in availability) ----
_IGHD_ARG=""
ls "${_ABS_DB}/imgt_${ORGANISM}_IGHD".n?? &>/dev/null 2>&1 && \
    _IGHD_ARG="-germline_db_D ${_ABS_DB}/imgt_${ORGANISM}_IGHD"
_C_ARG=""
ls "${_ABS_DB}/imgt_${ORGANISM}_C".n?? &>/dev/null 2>&1 && \
    _C_ARG="-c_region_db ${_ABS_DB}/imgt_${ORGANISM}_C"
_NDM_ARG=""
[[ -f "\${NDM_FILE}" ]] && _NDM_ARG="-custom_internal_data \${NDM_FILE}"

# ---- Pre-flight: verify required databases ----
_missing=0
for _db in "\${DB_PREFIX}_IGHV" "\${DB_PREFIX}_IGHJ" \
           "\${DB_PREFIX}_IGKV" "\${DB_PREFIX}_IGKJ" \
           "\${DB_PREFIX}_IGLV" "\${DB_PREFIX}_IGLJ"; do
  ls "\${_db}".n?? &>/dev/null || { echo "[ERROR] Missing: \${_db}" >&2; _missing=1; }
done
(( _missing == 0 )) || { echo "[ERROR] Re-run build script to rebuild databases." >&2; exit 1; }

# ---- Logging ----
_LOG="\${OUT_PREFIX}_igblast_run.log"
exec > >(tee -a "\${_LOG}") 2>&1

echo "=== igblast run: \$(date) ==="
echo "Query    : \${QUERY_FASTA}"
echo "Prefix   : \${OUT_PREFIX}"
echo "Chain    : \${CHAIN}"
echo "IGDATA   : \${IGDATA}"
echo "DB       : \${DB_PREFIX}"
echo "AUX      : \${AUX_FILE}"
echo "Organism : \${ORGANISM}"
echo "NDM file : \${NDM_FILE}  \$([[ -n \${_NDM_ARG} ]] && echo '(active)' || echo '(not found - using internal)')"
[[ -n "\${_IGHD_ARG}" ]] && echo "D db     : ${_ABS_DB}/imgt_${ORGANISM}_IGHD" || echo "D db     : (none)"
[[ -n "\${_C_ARG}"    ]] && echo "C db     : ${_ABS_DB}/imgt_${ORGANISM}_C"    || echo "C db     : (none)"
echo ""

# ---- Diagnostics ----
echo "--- internal_data/${ORGANISM}/ ---"
ls -1 "\${IGDATA}/internal_data/${ORGANISM}/" 2>/dev/null | sed 's/^/  /' || echo "  (missing)"
echo ""

# ---- Helper: run igblastn with argument logging ----
_run_igblast() {
  local label="\$1"; shift
  echo "--- Running \${label} ---"
  igblastn "\$@"
  local rc=\$?
  [[ \$rc -ne 0 ]] && echo "  [ERROR] igblastn exited with code \${rc}" >&2
  return \$rc
}

# ---- IGH (heavy chain, VDJ) ----
if [[ "\${CHAIN}" == "heavy" || "\${CHAIN}" == "all" ]]; then
  _run_igblast "IGH (AIRR outfmt 19)" \
      -germline_db_V  "\${DB_PREFIX}_IGHV" \
      -germline_db_J  "\${DB_PREFIX}_IGHJ" \
      \${_IGHD_ARG} \
      -auxiliary_data "\${AUX_FILE}" \
      \${_C_ARG} \
      -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \
      \${_NDM_ARG} \
      -outfmt 19 \
      -query  "\${QUERY_FASTA}" \
      -out    "\${OUT_PREFIX}.igh.airr.tsv"
  echo "  -> \${OUT_PREFIX}.igh.airr.tsv"

  _run_igblast "IGH (fmt7 for MakeDb.py)" \
      -germline_db_V  "\${DB_PREFIX}_IGHV" \
      -germline_db_J  "\${DB_PREFIX}_IGHJ" \
      \${_IGHD_ARG} \
      -auxiliary_data "\${AUX_FILE}" \
      \${_C_ARG} \
      -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \
      \${_NDM_ARG} \
      -outfmt '7 std qseq sseq btop' \
      -query  "\${QUERY_FASTA}" \
      -out    "\${OUT_PREFIX}.igh.fmt7"
  echo "  -> \${OUT_PREFIX}.igh.fmt7"
fi

# ---- IGK (kappa light chain, VJ) ----
if [[ "\${CHAIN}" == "light" || "\${CHAIN}" == "all" ]]; then
  _run_igblast "IGK (AIRR outfmt 19)" \
      -germline_db_V  "\${DB_PREFIX}_IGKV" \
      -germline_db_J  "\${DB_PREFIX}_IGKJ" \
      -auxiliary_data "\${AUX_FILE}" \
      \${_C_ARG} \
      -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \
      \${_NDM_ARG} \
      -outfmt 19 \
      -query  "\${QUERY_FASTA}" \
      -out    "\${OUT_PREFIX}.igk.airr.tsv"
  echo "  -> \${OUT_PREFIX}.igk.airr.tsv"

# ---- IGL (lambda light chain, VJ) ----
  _run_igblast "IGL (AIRR outfmt 19)" \
      -germline_db_V  "\${DB_PREFIX}_IGLV" \
      -germline_db_J  "\${DB_PREFIX}_IGLJ" \
      -auxiliary_data "\${AUX_FILE}" \
      \${_C_ARG} \
      -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \
      \${_NDM_ARG} \
      -outfmt 19 \
      -query  "\${QUERY_FASTA}" \
      -out    "\${OUT_PREFIX}.igl.airr.tsv"
  echo "  -> \${OUT_PREFIX}.igl.airr.tsv"
fi

echo ""
echo "=== Done: \$(date) ==="
echo "Log: \${_LOG}"
IGBLAST_EOF2
chmod +x "$IGBLAST_CMD"
ok "igblastn command template: $IGBLAST_CMD"
ok "igblastn command template: $IGBLAST_CMD"


# ---------------------------------------------------------------------------
# Step 6: Generate per-chain igblast+MakeDb pipeline scripts and install script
# ---------------------------------------------------------------------------
header "Step 6: Generating pipeline scripts"

_ABS_OUTDIR="$(cd "$OUTDIR" && pwd)"
_ABS_DB="${_ABS_OUTDIR}/database"
_ABS_GAPPED="${_ABS_OUTDIR}/germlines/gapped"
_ABS_AUX="${_ABS_OUTDIR}/auxiliary/${ORGANISM}_gl.aux"
_ABS_NDM="${_ABS_OUTDIR}/auxiliary/${ORGANISM}.ndm.imgt"

# ---- Helper: write a single-chain igblast+MakeDb script ----
write_chain_script() {
  local script_path="$1"
  local chain_label="$2"
  local v_db_args="$3"
  local d_db_arg="$4"
  local ref_fastas="$5"

  # Write the script in one unquoted heredoc.
  # Static help text that must NOT expand uses literal dollar signs (\$).
  # Build-time values (paths, species) ARE expanded by the unquoted heredoc.
  cat > "$script_path" <<CHAIN_SCRIPT
#!/usr/bin/env bash
# =============================================================================
# ${chain_label} chain: igblastn -> MakeDb.py pipeline
# Generated by build_blendAIRR.sh
# =============================================================================
# DESCRIPTION
#   Runs igblastn against the hybrid ${chain_label,,} chain germline databases then
#   parses the output with MakeDb.py to produce a Change-O / AIRR database.
#   All database paths are embedded at build time -- no configuration needed.
#
# USAGE
#   bash $(basename "$script_path") <query.fasta> [out_prefix] [makedb_outdir]
#
# ARGUMENTS
#   query.fasta     REQUIRED. FASTA file of BCR sequences to annotate.
#                   Should be full-length or near-full-length V(D)J contigs.
#
#   out_prefix      Output file prefix (default: igblast_out).
#                   All output files will be named <out_prefix>.*
#
#   makedb_outdir   Directory for MakeDb.py output (default: current dir).
#
# OUTPUTS
#   <out_prefix>.${chain_label,,}.fmt7         igblastn alignment (fmt7 for MakeDb.py)
#   <out_prefix>.${chain_label,,}.airr.tsv     igblastn AIRR-format annotation
#   <out_prefix>.${chain_label,,}_igblast.log  Full run log
#   <makedb_outdir>/*_db-pass.tsv              Change-O database (pass)
#   <makedb_outdir>/*_db-fail.tsv              Sequences failing annotation
#
# EMBEDDED PATHS (baked in at build time)
#   DB prefix : ${_ABS_DB}/${PREFIX}
#   Aux file  : ${_ABS_AUX}
#   IGDATA    : ${_ABS_OUTDIR}
#   Organism  : ${SPECIES}
# =============================================================================

set -euo pipefail
[[ "\${1:-}" == "-h" || "\${1:-}" == "--help" ]] && {
  sed -n '/^# DESCRIPTION/,/^# ===/p' "\$0" | sed 's/^# \?//'; exit 0; }

QUERY_FASTA="\${1:?Usage: \$0 <query.fasta> [out_prefix] [makedb_outdir]}"
OUT_PREFIX="\${2:-igblast_out}"
MAKEDB_OUTDIR="\${3:-.}"

# Baked-in paths
export IGDATA="${_ABS_OUTDIR}"
DB_PREFIX="${_ABS_DB}/imgt_${ORGANISM}"
AUX_FILE="${_ABS_AUX}"
NDM_FILE="${_ABS_NDM}"
ORGANISM="${SPECIES}"

_NDM_ARG=""
[[ -f "\${NDM_FILE}" ]] && _NDM_ARG="-custom_internal_data \${NDM_FILE}"
_C_ARG=""
ls "${_ABS_DB}/imgt_${ORGANISM}_C".n?? &>/dev/null 2>&1 && _C_ARG="-c_region_db ${_ABS_DB}/imgt_${ORGANISM}_C"

_LOG="\${OUT_PREFIX}.${chain_label,,}_igblast.log"
exec > >(tee -a "\${_LOG}") 2>&1
echo "=== ${chain_label} chain igblast+MakeDb: \$(date) ==="
echo "Query  : \${QUERY_FASTA}"
echo "Prefix : \${OUT_PREFIX}"
echo "IGDATA : \${IGDATA}"
echo "DB     : \${DB_PREFIX}"
echo ""

# ---- igblastn fmt7 (for MakeDb.py) ----
echo "--- igblastn fmt7 ---"
igblastn \\
    ${v_db_args} \\
    ${d_db_arg} \\
    -auxiliary_data "\${AUX_FILE}" \\
    \${_C_ARG} \\
    -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \\
    \${_NDM_ARG} \\
    -outfmt '7 std qseq sseq btop' \\
    -query  "\${QUERY_FASTA}" \\
    -out    "\${OUT_PREFIX}.${chain_label,,}.fmt7"
echo "  -> \${OUT_PREFIX}.${chain_label,,}.fmt7"

# ---- igblastn AIRR (outfmt 19) ----
echo "--- igblastn AIRR ---"
igblastn \\
    ${v_db_args} \\
    ${d_db_arg} \\
    -auxiliary_data "\${AUX_FILE}" \\
    \${_C_ARG} \\
    -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \\
    \${_NDM_ARG} \\
    -outfmt 19 \\
    -query  "\${QUERY_FASTA}" \\
    -out    "\${OUT_PREFIX}.${chain_label,,}.airr.tsv"
echo "  -> \${OUT_PREFIX}.${chain_label,,}.airr.tsv"

# ---- MakeDb.py ----
echo "--- MakeDb.py ---"
mkdir -p "\${MAKEDB_OUTDIR}"
MakeDb.py igblast \\
    -i "\${OUT_PREFIX}.${chain_label,,}.fmt7" \\
    -s "\${QUERY_FASTA}" \\
    -r ${ref_fastas} \\
    --extended --failed --asis-calls --partial \\
    --outdir "\${MAKEDB_OUTDIR}"
echo "  -> \${MAKEDB_OUTDIR}/"

echo ""
echo "=== Done: \$(date) === Log: \${_LOG}"
CHAIN_SCRIPT

  chmod +x "$script_path"
  ok "  $(basename "$script_path")"
}
# Heavy chain V/D/J db args and ref FASTAs
_IGH_V_ARGS="-germline_db_V ${_ABS_DB}/imgt_${ORGANISM}_IGHV -germline_db_J ${_ABS_DB}/imgt_${ORGANISM}_IGHJ"
_IGH_D_ARG=""
ls "${_ABS_DB}/imgt_${ORGANISM}_IGHD".n?? &>/dev/null 2>&1 && \
    _IGH_D_ARG="-germline_db_D ${_ABS_DB}/imgt_${ORGANISM}_IGHD"
_IGH_REFS="${_ABS_GAPPED}/${ORGANISM}_IGHV.fasta"
[[ -f "${_ABS_GAPPED}/${ORGANISM}_IGHD.fasta" ]] && _IGH_REFS+=" ${_ABS_GAPPED}/${ORGANISM}_IGHD.fasta"
_IGH_REFS+=" ${_ABS_GAPPED}/${ORGANISM}_IGHJ.fasta"

# Light chain V/J db args and ref FASTAs
_IGL_V_ARGS="-germline_db_V ${_ABS_DB}/imgt_${ORGANISM}_IGKV -germline_db_J ${_ABS_DB}/imgt_${ORGANISM}_IGKJ -germline_db_V ${_ABS_DB}/imgt_${ORGANISM}_IGLV -germline_db_J ${_ABS_DB}/imgt_${ORGANISM}_IGLJ"
_IGL_REFS="${_ABS_GAPPED}/${ORGANISM}_IGKV.fasta ${_ABS_GAPPED}/${ORGANISM}_IGKJ.fasta ${_ABS_GAPPED}/${ORGANISM}_IGLV.fasta ${_ABS_GAPPED}/${ORGANISM}_IGLJ.fasta"

HEAVY_CMD="$OUTDIR/${PREFIX}_run_heavy.sh"
LIGHT_CMD="$OUTDIR/${PREFIX}_run_light.sh"

write_chain_script "$HEAVY_CMD" "Heavy" "$_IGH_V_ARGS" "$_IGH_D_ARG" "$_IGH_REFS"
write_chain_script "$LIGHT_CMD" "Light" \
    "-germline_db_V ${_ABS_DB}/imgt_${ORGANISM}_IGKV -germline_db_J ${_ABS_DB}/imgt_${ORGANISM}_IGKJ" \
    "" \
    "${_ABS_GAPPED}/${ORGANISM}_IGKV.fasta ${_ABS_GAPPED}/${ORGANISM}_IGKJ.fasta ${_ABS_GAPPED}/${ORGANISM}_IGLV.fasta ${_ABS_GAPPED}/${ORGANISM}_IGLJ.fasta"

# Note: light chain script runs IGK only in igblastn; IGL requires a separate invocation.
# Append IGL igblastn calls to the light script
cat >> "$LIGHT_CMD" <<'IGL_PATCH_QUOTED'

# ---- IGL lambda: separate igblastn calls (different V/J database) ----
IGL_PATCH_QUOTED
cat >> "$LIGHT_CMD" <<IGL_PATCH_BAKED
echo "--- igblastn IGL fmt7 ---"
igblastn \\
    -germline_db_V  "${_ABS_DB}/imgt_${ORGANISM}_IGLV" \\
    -germline_db_J  "${_ABS_DB}/imgt_${ORGANISM}_IGLJ" \\
    -auxiliary_data "\${AUX_FILE}" \\
    \${_C_ARG} \\
    -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \\
    \${_NDM_ARG} \\
    -outfmt '7 std qseq sseq btop' \\
    -query  "\${QUERY_FASTA}" \\
    -out    "\${OUT_PREFIX}.igl.fmt7"
echo "  -> \${OUT_PREFIX}.igl.fmt7"
echo "--- igblastn IGL AIRR ---"
igblastn \\
    -germline_db_V  "${_ABS_DB}/imgt_${ORGANISM}_IGLV" \\
    -germline_db_J  "${_ABS_DB}/imgt_${ORGANISM}_IGLJ" \\
    -auxiliary_data "\${AUX_FILE}" \\
    \${_C_ARG} \\
    -domain_system imgt -ig_seqtype Ig -organism \${ORGANISM} \\
    \${_NDM_ARG} \\
    -outfmt 19 \\
    -query  "\${QUERY_FASTA}" \\
    -out    "\${OUT_PREFIX}.igl.airr.tsv"
echo "  -> \${OUT_PREFIX}.igl.airr.tsv"
IGL_PATCH_BAKED

# ---- Install script: copy files to an existing IGBLAST_DB without overwriting ----
INSTALL_CMD="$OUTDIR/${PREFIX}_install_to_igdata.sh"
cat > "$INSTALL_CMD" <<'INSTALL_EOF_QUOTED'
#!/usr/bin/env bash
# =============================================================================
# Install hybrid reference into an existing IgBLAST IGDATA directory.
# Generated by build_blendAIRR.sh
# =============================================================================
#
# DESCRIPTION
#   Copies the hybrid reference databases, germline FASTAs, auxiliary file,
#   and V-gene annotation (ndm.imgt) into an existing IgBLAST IGDATA
#   directory so they can be used with -organism <name> across any tool
#   that calls igblastn with that IGDATA (e.g. Immcantation pipelines).
#
#   SAFETY: This script NEVER overwrites existing files. Run with --dry-run
#   first to review exactly what would be copied.
#
# USAGE
#   bash <this_script> <IGDATA_DIR> [--dry-run]
#
# ARGUMENTS
#   IGDATA_DIR   Path to your IgBLAST share directory -- the one that contains
#                internal_data/, optional_file/, and database/ subdirectories.
#                This is typically ~/share/igblast or the directory pointed to
#                by your $IGDATA environment variable.
#
#   --dry-run          Preview what would be installed without copying any files.
#                      Strongly recommended before the first real install.
#   --organism-name N  Override the organism name (default: embedded at build time).
#   --germlines-root D Override the germlines root directory (default: auto-detected
#                      relative to IGDATA_DIR as <igdata>/../germlines or
#                      $HOME/share/germlines). Use this if your germlines are in
#                      a non-standard location, e.g.:
#                        --germlines-root /usr/local/share/germlines
#
# WHAT GETS INSTALLED
#   <germlines>/imgt/<organism>/vdj/        Hybrid gapped FASTAs (IGHV/D/J/IGK/IGL)
#   <germlines>/imgt/<organism>/constant/   Constant region FASTAs (from reference)
#   <germlines>/imgt/<organism>/leader/     Leader FASTAs (from reference)
#   <germlines>/imgt/<organism>/leader_vexon/ Leader+Vexon FASTAs (from reference)
#   <germlines>/imgt/<organism>/vdj_aa/     AA FASTAs (from reference)
#   database/imgt_<organism>_IGHV.*         Per-locus BLAST databases
#   database/imgt_<organism>_ig_v.*         Combined V/D/J BLAST databases
#   optional_file/<organism>_gl.aux         J-gene auxiliary file
#   internal_data/<organism>/               Internal BLAST DBs + ndm.imgt
#
# AFTER INSTALLATION
#   Use these flags in igblastn:
#     -organism      <ORGANISM>        (embedded below)
#     -germline_db_V database/<PREFIX>_IGHV
#     -germline_db_J database/<PREFIX>_IGHJ
#     -germline_db_D database/<PREFIX>_IGHD
#     -auxiliary_data optional_file/<AUX_FILENAME>
#
# EMBEDDED VALUES (set at build time)
#   Organism : (see script body below)
#   Prefix   : (see script body below)
# =============================================================================

set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '/^# DESCRIPTION/,/^# ===.*$/p' "$0" | grep '^#' | sed 's/^# \?//'
  exit 0
fi

IGDATA_TARGET="${1:?Usage: $0 <IGDATA_DIR> [--dry-run]}"
DRY_RUN=false
[[ "${2:-}" == "--dry-run" ]] && DRY_RUN=true

if [[ ! -d "$IGDATA_TARGET" ]]; then
  echo "[ERROR] IGDATA directory not found: $IGDATA_TARGET" >&2; exit 1
fi

INSTALL_EOF_QUOTED

cat >> "$INSTALL_CMD" <<INSTALL_EOF_BAKED
# Closest reference species (e.g. mouse)
REF_SPECIES="${SPECIES}"
# Hybrid organism name used for internal_data/ and igblastn -organism.
# Defaults to <prefix>_<species> (e.g. hybrid_mouse) to avoid colliding
# with the existing reference species files in IGDATA.
ORGANISM="${PREFIX}_${SPECIES}"
SRC_DB="${_ABS_DB}"
SRC_GAPPED="${_ABS_GAPPED}"
SRC_AUX="${_ABS_AUX}"
SRC_NDM="${_ABS_NDM}"
SRC_INTERNAL="${_ABS_OUTDIR}/internal_data/${ORGANISM}"
INSTALL_EOF_BAKED

cat >> "$INSTALL_CMD" <<'INSTALL_EOF_QUOTED2'

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
n_copied=0; n_skipped=0; n_would=0

while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --dry-run)        DRY_RUN=true;                shift ;;
    --organism-name)  ORGANISM="$2";               shift 2 ;;
    --germlines-root) GERMLINES_ROOT_OVERRIDE="$2"; shift 2 ;;
    *) echo "Unknown flag: $1" >&2; exit 1 ;;
  esac
done
# Allow caller to override germlines root location
[[ -n "${GERMLINES_ROOT_OVERRIDE:-}" ]] && GERMLINES_ROOT_OVERRIDE_SET=true || GERMLINES_ROOT_OVERRIDE_SET=false

echo "=== Installing ${ORGANISM} into ${IGDATA_TARGET} ==="
echo "  Reference species  : ${REF_SPECIES}"
echo "  Hybrid organism ID : ${ORGANISM}"
$DRY_RUN && echo "  [DRY RUN — no files will be written]"
echo ""

# safe_copy: copy src to dest_dir/dest_name; skip if dest already exists
safe_copy() {
  local src="$1" dest_dir="$2" dest_name="${3:-$(basename "$1")}"
  local dest="${dest_dir}/${dest_name}"
  [[ -f "$src" ]] || { echo -e "${YELLOW}[SKIP]${NC} source missing: $src"; return; }
  if [[ -e "$dest" ]]; then
    echo -e "${YELLOW}[SKIP]${NC} already exists: $dest"
    (( n_skipped++ )) || true; return
  fi
  if $DRY_RUN; then
    echo -e "${GREEN}[WOULD COPY]${NC} $src -> $dest"
    (( n_would++ )) || true
  else
    mkdir -p "$dest_dir"
    cp "$src" "$dest"
    echo -e "${GREEN}[COPIED]${NC} $dest"
    (( n_copied++ )) || true
  fi
}

# safe_copy_renamed: copy reference file, substituting species name in filename
# e.g. imgt_mouse_IGHV.fasta -> imgt_hybrid_mouse_IGHV.fasta
safe_copy_renamed() {
  local src="$1" dest_dir="$2" prefix_old="$3" prefix_new="$4"
  local base; base=$(basename "$src")
  local dest_name="${base/$prefix_old/$prefix_new}"
  safe_copy "$src" "$dest_dir" "$dest_name"
}

# ── 1. vdj/ — hybrid FASTAs with imgt_ prefix
#    Format: imgt_<organism>_IGHV.fasta  (matches Immcantation/igblastn convention)
echo "--- germlines/imgt/${ORGANISM}/vdj/ (hybrid FASTAs) ---"

# Resolve GERMLINES_ROOT: try common locations relative to IGDATA
GERMLINES_ROOT=""
for _candidate in     "${IGDATA_TARGET}/../germlines"     "${IGDATA_TARGET}/germlines"     "$(dirname "${IGDATA_TARGET}")/germlines"     "${HOME}/share/germlines"; do
  if [[ -d "${_candidate}/imgt" ]]; then
    GERMLINES_ROOT="$(cd "${_candidate}" && pwd)"
    break
  fi
done
if [[ -z "$GERMLINES_ROOT" ]]; then
  echo -e "${YELLOW}[WARN]${NC} Could not find germlines/imgt/ directory."
  echo "  Tried locations relative to: ${IGDATA_TARGET}"
  echo "  Set GERMLINES_ROOT manually and re-run, or copy FASTAs manually to:"
  echo "  <germlines_root>/imgt/${ORGANISM}/vdj/imgt_${ORGANISM}_IGHV.fasta ..."
  GERMLINES_ROOT="${IGDATA_TARGET}/../germlines"   # fallback, may not exist
fi
# Allow explicit override via --germlines-root flag
$GERMLINES_ROOT_OVERRIDE_SET && GERMLINES_ROOT="$GERMLINES_ROOT_OVERRIDE"
echo "  Germlines root: ${GERMLINES_ROOT}"

VDJ_DEST="${GERMLINES_ROOT}/imgt/${ORGANISM}/vdj"
REF_VDJ="${GERMLINES_ROOT}/imgt/${REF_SPECIES}/vdj"
for locus in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ; do
  src="${SRC_GAPPED}/imgt_${ORGANISM}_${locus}.fasta"
  # Fallback: the file may not have imgt_ prefix (older builds)
  [[ -f "$src" ]] || src="${SRC_GAPPED}/${ORGANISM}_${locus}.fasta"
  if [[ -f "$src" ]]; then
    safe_copy "$src" "$VDJ_DEST" "imgt_${ORGANISM}_${locus}.fasta"
  else
    # Copy from reference with name substitution
    ref_src="${REF_VDJ}/imgt_${REF_SPECIES}_${locus}.fasta"
    [[ -f "$ref_src" ]] && \
      safe_copy_renamed "$ref_src" "$VDJ_DEST" \
        "imgt_${REF_SPECIES}_" "imgt_${ORGANISM}_"
  fi
done

# ── 2. constant/ — copy reference files, rename species prefix
echo ""
echo "--- germlines/imgt/${ORGANISM}/constant/ (from reference) ---"
CONST_SRC="${GERMLINES_ROOT}/imgt/${REF_SPECIES}/constant"
CONST_DEST="${GERMLINES_ROOT}/imgt/${ORGANISM}/constant"
if [[ -d "$CONST_SRC" ]]; then
  for f in "$CONST_SRC"/imgt_${REF_SPECIES}_*.fasta; do
    [[ -f "$f" ]] && safe_copy_renamed "$f" "$CONST_DEST" \
      "imgt_${REF_SPECIES}_" "imgt_${ORGANISM}_"
  done
else
  echo -e "${YELLOW}[SKIP]${NC} reference constant/ not found: $CONST_SRC"
fi

# ── 3. leader/ — copy reference files, rename
echo ""
echo "--- germlines/imgt/${ORGANISM}/leader/ (from reference) ---"
LEADER_SRC="${GERMLINES_ROOT}/imgt/${REF_SPECIES}/leader"
LEADER_DEST="${GERMLINES_ROOT}/imgt/${ORGANISM}/leader"
if [[ -d "$LEADER_SRC" ]]; then
  for f in "$LEADER_SRC"/imgt_${REF_SPECIES}_*.fasta; do
    [[ -f "$f" ]] && safe_copy_renamed "$f" "$LEADER_DEST" \
      "imgt_${REF_SPECIES}_" "imgt_${ORGANISM}_"
  done
else
  echo -e "${YELLOW}[SKIP]${NC} reference leader/ not found: $LEADER_SRC"
fi

# ── 4. leader_vexon/ — copy reference files (prefix: imgt_lv_)
echo ""
echo "--- germlines/imgt/${ORGANISM}/leader_vexon/ (from reference) ---"
LV_SRC="${GERMLINES_ROOT}/imgt/${REF_SPECIES}/leader_vexon"
LV_DEST="${GERMLINES_ROOT}/imgt/${ORGANISM}/leader_vexon"
if [[ -d "$LV_SRC" ]]; then
  for f in "$LV_SRC"/imgt_lv_${REF_SPECIES}_*.fasta; do
    [[ -f "$f" ]] && safe_copy_renamed "$f" "$LV_DEST" \
      "imgt_lv_${REF_SPECIES}_" "imgt_lv_${ORGANISM}_"
  done
else
  echo -e "${YELLOW}[SKIP]${NC} reference leader_vexon/ not found: $LV_SRC"
fi

# ── 5. vdj_aa/ — translate hybrid V FASTAs to amino acid
#    Translate the actual hybrid gapped nucleotide V sequences (in-frame,
#    frame 1, stripping gaps) to produce organism-specific AA sequences.
#    Falls back to reference copies only if translation fails.
echo ""
echo "--- germlines/imgt/${ORGANISM}/vdj_aa/ (translated from hybrid V FASTAs) ---"
AA_DEST="${GERMLINES_ROOT}/imgt/${ORGANISM}/vdj_aa"
mkdir -p "$AA_DEST"

_translate_v_fasta() {
  local src_nt="$1"    # gapped nucleotide FASTA
  local dest_aa="$2"   # output amino acid FASTA
  [[ -f "$src_nt" ]] || return 1

  python3 - "$src_nt" "$dest_aa" << 'TRANSLATE_EOF'
import sys, textwrap
CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
def translate(nt):
    nt = nt.upper().replace('-','').replace('.','').replace(' ','')
    aa = [CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3)]
    return ''.join(aa).rstrip('*') or 'X'
with open(sys.argv[1]) as fin, open(sys.argv[2],'w') as fout:
    name,seq = None,[]
    for line in fin:
        line=line.rstrip()
        if line.startswith('>'):
            if name and seq:
                p=translate(''.join(seq))
                if p!='X': fout.write(f'>{name}\n'+''.join(c+'\n' for c in textwrap.wrap(p,60)))
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line)
    if name and seq:
        p=translate(''.join(seq))
        if p!='X': fout.write(f'>{name}\n'+''.join(c+'\n' for c in textwrap.wrap(p,60)))
TRANSLATE_EOF
}

# Source FASTAs: prefer already-installed vdj/, fall back to blendAIRR output dir
_translated_any=false
for locus in IGHV IGKV IGLV; do
  src_nt="${VDJ_DEST}/imgt_${ORGANISM}_${locus}.fasta"
  [[ -f "$src_nt" ]] || src_nt="${SRC_GAPPED}/imgt_${ORGANISM}_${locus}.fasta"
  dest_aa="${AA_DEST}/imgt_aa_${ORGANISM}_${locus}.fasta"
  if [[ -e "$dest_aa" ]]; then
    echo -e "${YELLOW}[SKIP]${NC} already exists: imgt_aa_${ORGANISM}_${locus}.fasta"
    (( n_skipped++ )) || true
    _translated_any=true
  elif [[ -f "$src_nt" ]]; then
    if $DRY_RUN; then
      echo -e "${GREEN}[WOULD TRANSLATE]${NC} ${locus}: $src_nt -> $dest_aa"
      (( n_would++ )) || true
    else
      if _translate_v_fasta "$src_nt" "$dest_aa"; then
        n_seqs=$(grep -c "^>" "$dest_aa" 2>/dev/null || echo 0)
        echo -e "${GREEN}[TRANSLATED]${NC} imgt_aa_${ORGANISM}_${locus}.fasta (${n_seqs} seqs)"
        (( n_copied++ )) || true
        _translated_any=true
      else
        echo -e "${YELLOW}[WARN]${NC} translation failed for ${locus}"
      fi
    fi
  else
    echo -e "${YELLOW}[SKIP]${NC} source NT FASTA not found: $src_nt"
  fi
done

# Fallback: copy reference vdj_aa for loci not translated above
AA_SRC="${GERMLINES_ROOT}/imgt/${REF_SPECIES}/vdj_aa"
if [[ -d "$AA_SRC" && "$_translated_any" == "false" ]]; then
  echo "  Falling back to reference vdj_aa copies..."
  for f in "$AA_SRC"/imgt_aa_${REF_SPECIES}_*.fasta; do
    [[ -f "$f" ]] && safe_copy_renamed "$f" "$AA_DEST" \
      "imgt_aa_${REF_SPECIES}_" "imgt_aa_${ORGANISM}_"
  done
fi

# ── 6. BLAST databases -> database/
#    Two naming schemes installed side-by-side:
#    a) Per-locus:  imgt_<organism>_IGHV.*  (for -germline_db_V per call)
#    b) Combined:   imgt_<organism>_ig_v.*  (matches reference convention;
#                   required by some tools that auto-detect by organism name)
#
#    Locus -> segment mapping:
#      IGHV/IGKV/IGLV -> ig_v   IGHD -> ig_d   IGHJ/IGKJ/IGLJ -> ig_j
echo ""
echo "--- Germline BLAST databases -> database/ ---"

# a) Per-locus databases (prefix imgt_<organism>_LOCUS)
# Exclude .nal alias files — they embed absolute build-time paths
# and will be broken after copying to a different location.
for locus in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ; do
  for ext in nhr nin nsq nsi nsd nog njs nto ntf not nos ndb; do
    src="${SRC_DB}/imgt_${ORGANISM}_${locus}.${ext}"
    [[ -f "$src" ]] && \
      safe_copy "$src" "${IGDATA_TARGET}/database"
  done
done

# b) Combined ig_v / ig_d / ig_j databases
#    Merge per-locus FASTAs and build with makeblastdb so real .nhr/.nin/.nsq
#    index files exist. No .nal alias files — they cause issues with some tools.
DB_DEST="${IGDATA_TARGET}/database"
VDJ_SRC="${GERMLINES_ROOT}/imgt/${ORGANISM}/vdj"

_build_combined_db() {
  local seg="$1"; shift         # e.g. ig_v
  local loci=("$@")             # e.g. IGHV IGKV IGLV
  local out_name="imgt_${ORGANISM}_${seg}"
  local out_db="${DB_DEST}/${out_name}"

  if [[ -e "${out_db}.nhr" ]]; then
    echo -e "${YELLOW}[SKIP]${NC} already exists: ${out_name}.nhr"
    (( n_skipped++ )) || true; return
  fi

  # Merge available per-locus FASTAs
  local tmp_merged; tmp_merged=$(mktemp /tmp/blendairr_merged_XXXXXX.fasta)
  local found=0
  for locus in "${loci[@]}"; do
    local fa="${VDJ_SRC}/imgt_${ORGANISM}_${locus}.fasta"
    # Fallback to blendAIRR output dir if not yet installed
    [[ -f "$fa" ]] || fa="${SRC_GAPPED}/imgt_${ORGANISM}_${locus}.fasta"
    if [[ -f "$fa" ]]; then
      cat "$fa" >> "$tmp_merged"
      (( found++ )) || true
    fi
  done

  if [[ $found -eq 0 ]]; then
    echo -e "${YELLOW}[SKIP]${NC} no source FASTAs found for ${out_name}"
    rm -f "$tmp_merged"; return
  fi

  if $DRY_RUN; then
    echo -e "${GREEN}[WOULD BUILD]${NC} ${out_name} (${found} loci merged)"
    (( n_would++ )) || true
    rm -f "$tmp_merged"; return
  fi

  mkdir -p "$DB_DEST"
  # Save processed FASTA permanently alongside database files
  # (mirrors reference: imgt_mouse_ig_v.fasta next to imgt_mouse_ig_v.nhr)
  local fasta_dir; fasta_dir="$(dirname "$DB_DEST")/fasta"
  mkdir -p "$fasta_dir"
  local proc_fasta="${fasta_dir}/${out_name}.fasta"

  if command -v edit_imgt_file.pl &>/dev/null; then
    perl "$(command -v edit_imgt_file.pl)" "$tmp_merged" > "$proc_fasta"
  else
    awk '/^>/{sub(/ .*/,""); print} !/^>/{gsub(/\./, ""); print}' \
      "$tmp_merged" > "$proc_fasta"
  fi

  makeblastdb -parse_seqids -dbtype nucl \
    -in  "$proc_fasta" \
    -out "$out_db" \
    2>/dev/null \
    && echo -e "${GREEN}[BUILT]${NC} ${out_name}" \
    || echo -e "${YELLOW}[WARN]${NC} makeblastdb failed for ${out_name}"
  (( n_copied++ )) || true
  rm -f "$tmp_merged"
}

echo ""
echo "--- Combined BLAST databases (ig_v / ig_d / ig_j) ---"
_build_combined_db "ig_v" IGHV IGKV IGLV
_build_combined_db "ig_d" IGHD
_build_combined_db "ig_j" IGHJ IGKJ IGLJ

# ── 7. Auxiliary file -> optional_file/
echo ""
echo "--- Auxiliary file -> optional_file/ ---"
safe_copy "$SRC_AUX" "${IGDATA_TARGET}/optional_file"

# ── 8. V gene annotation and domain models -> internal_data/<organism>/
#
#    .ndm.imgt  — NT FWR/CDR positions in the UNGAPPED sequence.
#                Generated de novo by blendAIRR from the hybrid gapped FASTAs.
#
#    .pdm.imgt  — Protein FWR/CDR positions (IMGT numbering).
#                IMGT protein positions are conserved across all Ig V gene
#                families in a species (e.g. VH aa 1-25 = FWR1 for all IGHV),
#                so the reference .pdm.imgt is valid for hybrid sequences.
#                Generating this de novo requires the undocumented IgBLAST
#                binary format; copying from reference is correct and safe.
#
#    .ndm.kabat/.pdm.kabat — Kabat-numbered equivalents. Only consulted
#                when -domain_system kabat is used. blendAIRR always uses
#                -domain_system imgt, so these are safe reference copies.
echo ""
echo "--- V gene annotation + domain models -> internal_data/${ORGANISM}/ ---"
echo "  ndm.imgt: generated from hybrid sequences"
echo "  pdm.imgt, ndm/pdm.kabat: copied from reference (IMGT positions conserved)"
INT_DEST="${IGDATA_TARGET}/internal_data/${ORGANISM}"
REF_INT="${IGDATA_TARGET}/internal_data/${REF_SPECIES}"
safe_copy "$SRC_NDM" "$INT_DEST" "${ORGANISM}.ndm.imgt"
safe_copy "${REF_INT}/${REF_SPECIES}.pdm.imgt"  "$INT_DEST" "${ORGANISM}.pdm.imgt"
safe_copy "${REF_INT}/${REF_SPECIES}.ndm.kabat" "$INT_DEST" "${ORGANISM}.ndm.kabat"
safe_copy "${REF_INT}/${REF_SPECIES}.pdm.kabat" "$INT_DEST" "${ORGANISM}.pdm.kabat"

# ── 9. Internal V-gene BLAST databases -> internal_data/<organism>/
#    Only _V is needed — reference internal_data/ has no _D or _J databases.
#    Both n* (nucleotide) and p* (protein) required by igblastn.
echo ""
echo "--- Internal BLAST databases -> internal_data/${ORGANISM}/ (V only) ---"
for ext in nhr nin nsq nsi nsd nog njs nto ntf not nos ndb \
           phr pin psq psi psd pog pjs pto ptf pos pot pdb; do
  src="${SRC_INTERNAL}/${ORGANISM}_V.${ext}"
  [[ -f "$src" ]] && safe_copy "$src" "$INT_DEST" "${ORGANISM}_V.${ext}"
done

echo ""
if $DRY_RUN; then
  echo "Dry run: ${n_would} would be copied, ${n_skipped} already exist."
  echo "Re-run without --dry-run to install."
else
  echo "Install complete: ${n_copied} copied, ${n_skipped} skipped."
  echo ""
  echo "────────────────────────────────────────────────────────"
  echo "  To use with igblastn:"
  echo ""
  echo "    export IGDATA=${IGDATA_TARGET}"
  echo "    igblastn \\"
  echo "      -organism       ${ORGANISM} \\"
  echo "      -germline_db_V  database/${ORGANISM}_IGHV \\"
  echo "      -germline_db_J  database/${ORGANISM}_IGHJ \\"
  echo "      -germline_db_D  database/${ORGANISM}_IGHD \\"
  echo "      -auxiliary_data optional_file/$(basename "$SRC_AUX") \\"
  echo "      ..."
  echo ""
  echo "  Germline FASTAs for MakeDb.py -r:"
  echo "    germlines/imgt/${ORGANISM}/vdj/imgt_${ORGANISM}_IGHV.fasta ..."
  echo "────────────────────────────────────────────────────────"
fi
INSTALL_EOF_QUOTED2

chmod +x "$INSTALL_CMD"
ok "  $(basename "$INSTALL_CMD")"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
header "Build complete"
echo ""
echo "  Species       : $SPECIES"
echo "  Input dir     : $INPUT_DIR"
echo "  Output root   : $OUTDIR"
echo "  Prefix        : $PREFIX"
echo ""
echo "  Key outputs:"
echo "    Gapped FASTAs           : $OUTDIR/germlines/gapped/"
echo "    BLAST databases         : $OUTDIR/database/"
echo "    Aux file                : $_ABS_AUX"
echo "    ndm.imgt                : $_ABS_NDM"
echo "    Annotation tables       : $OUTDIR/annotations/"
echo "    Heavy chain pipeline    : $HEAVY_CMD"
echo "    Light chain pipeline    : $LIGHT_CMD"
echo "    Install to IGDATA       : $INSTALL_CMD"
echo "    Manifest                : $OUTDIR/${PREFIX}_manifest.tsv"
echo ""
echo "  Usage:"
echo "    Heavy chain: bash $HEAVY_CMD <query.fasta> <out_prefix> [makedb_outdir]"
echo "    Light chain: bash $LIGHT_CMD <query.fasta> <out_prefix> [makedb_outdir]"
echo "    Install:     bash $INSTALL_CMD <IGDATA_DIR> [--dry-run]"
echo ""
