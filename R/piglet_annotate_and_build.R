#!/usr/bin/env Rscript
# =============================================================================
# piglet_annotate_and_build.R
#
# PURPOSE:
#   Build a hybrid IgBLAST germline reference by merging custom (novel-species)
#   V/D/J sequences with the closest IMGT reference species.  Novel sequences
#   are re-annotated using PIgLET joint clustering so they inherit correct
#   gene-family names from co-clustering with known reference alleles.
#
# CRITICAL PIGLET API NOTE:
#   inferAlleleClusters() requires germline_set to be a NAMED character vector:
#     setNames(as.character(seqs), names(seqs))
#   Passing as.character(DNAStringSet) without names causes the internal
#   germ.dist subscript error because alleleClusterTable$imgt_allele holds the
#   input names but germ.dist rows are unnamed -> subscript out of bounds.
#
# HEADER FORMATS HANDLED (custom input files):
#   A: >IGHV1-2*01
#   B: >IGHV1-14*01_S3227          (suffix junk after allele)
#   C: >V00762|IGHJ1*01|Mus_...|   (full IMGT pipe-delimited)
#   D: >IGKV0-2HY3*00              (OGRDB alphanumeric hash, kept verbatim)
#
# OUTPUTS (all under --outdir):
#   annotations/  *_header_normalisation_map.tsv  raw -> normalised -> final
#   annotations/  *_allele_cluster_annotation.tsv  PIgLET cluster table
#   germlines/gapped/   per-locus + ALL_V/ALL_J  (for MakeDb.py -r)
#   germlines/ungapped/ per-locus                 (intermediate)
#   auxiliary/    *_<species>_hybrid_gl.aux        (IgBLAST CDR3 anchors)
#   auxiliary/    *.aux.diagnostic                 (with anchor_method column)
#   *_manifest.tsv
#
# USAGE (standalone or called by build_hybrid_igblast_ref.sh):
#   Rscript piglet_annotate_and_build.R \
#     --custom_dir /path/to/custom/ \
#     --ref_dir    /path/to/igblast/germlines/imgt/mouse/vdj \
#     --species    mouse \
#     --outdir     /path/to/output \
#     --igdata     /path/to/igblast/share
# =============================================================================

# =============================================================================
# Expand .libPaths() to include previous R version libraries and common
# system library paths. This allows the script to find packages installed
# under a different R minor version (e.g. rocker base image R 4.3.x libs
# when running under R 4.4.x), and packages installed in user or site
# libraries not on the default search path.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DECIPHER)
  library(piglet)
  library(Biostrings)
  library(stringr)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--custom_dir",               type = "character"),
  make_option("--ref_dir",                  type = "character"),
  make_option("--species",                  type = "character", default = "mouse"),
  make_option("--outdir",                   type = "character"),
  make_option("--igdata",                   type = "character"),
  make_option("--prefix",                   type = "character", default = "hybrid"),
  make_option("--family_threshold",         type = "double",  default = 75),
  make_option("--allele_cluster_threshold", type = "double",  default = 95),
  make_option("--v_trim3prime",             type = "integer", default = 318),
  make_option("--j_trim3prime",             type = "integer", default = 40),
  make_option("--use_asc",                  action = "store_true", default = FALSE,
              help = "Use PIgLET ASC names (IGHVFx-Gy*01) instead of IMGT-style names"),
  make_option("--trig_nc",                  action = "store_true", default = FALSE,
              help = "Emit TR/IG Nomenclature Review Committee colon format (IGHV:01:001:001) instead of asterisk format. Requires --use_asc."),
  make_option("--organism",                 type = "character",    default = NULL,
              help = "Organism name used as output file prefix (default: <prefix>_<species>)"),
  make_option("--as_is_ids",                action = "store_true", default = FALSE,
              help = "Accept input allele names as-is; skip PIgLET clustering entirely. Duplicate names get _1, _2 suffix. Sequence-level deduplication still applied.")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("custom_dir","ref_dir","outdir","igdata"))
  if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))
if (isTRUE(opt$trig_nc) && !isTRUE(opt$use_asc))
  stop("--trig_nc requires --use_asc (TRIG-NC naming operates on ASC clusters)")

for (d in c(opt$outdir,
            file.path(opt$outdir, "germlines","gapped"),
            file.path(opt$outdir, "germlines","ungapped"),
            file.path(opt$outdir, "fasta"),
            file.path(opt$outdir, "database"),
            file.path(opt$outdir, "annotations"),
            file.path(opt$outdir, "auxiliary")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Hybrid IgBLAST Reference Builder ===\n")
cat(sprintf("Custom dir : %s\n", opt$custom_dir))
cat(sprintf("Reference  : %s  species=%s\n", opt$ref_dir, opt$species))
cat(sprintf("Naming mode: %s\n",
            if (isTRUE(opt$trig_nc)) "TRIG-NC (colon format, ASC-based)"
            else if (isTRUE(opt$use_asc)) "ASC (PIgLET cluster names)"
            else "IMGT (reference-based)"))
cat(sprintf("ASC thresholds: family=%.0f%%  allele=%.0f%%\n",
            opt$family_threshold, opt$allele_cluster_threshold))

# file_prefix: used for ALL output file names (FASTAs, databases, aux, ndm.imgt).
# Defaults to <prefix>_<species> (e.g. "hybrid_mouse") for consistency with the
# bash script's ORGANISM variable. Can be overridden with --organism.
file_prefix <- ifelse(!is.null(opt$organism) && nchar(trimws(opt$organism)) > 0,
                      opt$organism,
                      paste0(opt$prefix, "_", opt$species))
# imgt_file_prefix: used for gapped FASTA filenames to match Immcantation convention
# e.g. imgt_hybrid_mouse_IGHV.fasta  (mirrors imgt_mouse_IGHV.fasta in reference)
imgt_file_prefix <- paste0("imgt_", file_prefix)
cat(sprintf("File prefix: %s  (gapped FASTAs: %s_*.fasta)\n\n",
            file_prefix, imgt_file_prefix))

# =============================================================================
# Shared J-gene aux coordinate helpers (used by BOTH as-is and ASC paths)
# =============================================================================

find_j_anchor <- function(nt_seq, chain_type) {
  # The conserved J-REGION CDR3 anchor is the Trp/Phe of the [WF]-G-X-G motif
  # (J-TRP for IGH, J-PHE for IGK/IGL). Use a general regex rather than an
  # enumerated list so unusual but valid motifs (e.g. FGGG, WGRG) are still
  # matched. The anchor is the nt position (0-based) of the conserved W/F codon.
  #
  # Chain preference for the leading residue:
  #   IGH -> prefer W (Trp), then F ; IGK/IGL -> prefer F (Phe), then W
  # We scan all three frames; within a frame we take the FIRST [WF]G.G match.
  lead <- switch(chain_type, IGH = c("W","F"), IGK = c("F","W"),
                 IGL = c("F","W"), c("W","F"))
  seq_obj  <- DNAString(nt_seq)
  best_pos <- NA_integer_
  best_motif <- NA_character_
  for (frame in 0L:2L) {
    sublen <- nchar(nt_seq) - frame
    sublen <- sublen - (sublen %% 3L)
    if (sublen < 3L) next
    aa <- as.character(translate(subseq(seq_obj, frame+1L, frame+sublen)))
    # Try preferred lead residue first, then the alternate.
    for (L in lead) {
      m <- regexpr(paste0(L, "G.G"), aa)   # e.g. "FG.G" / "WG.G"
      if (m > 0L) {
        # (m-1)*3+frame is the 0-based nt start of the W/F codon. The IMGT
        # reference cdr3_stop convention sits 1 nt earlier (validated against
        # 25/26 mouse reference J alleles: motif was systematically +1), so
        # subtract 1 to align exactly with the reference aux.
        best_pos   <- (m - 1L) * 3L + frame - 1L
        best_motif <- substr(aa, m, m + attr(m, "match.length") - 1L)  # e.g. "FGGG"
        break
      }
    }
    if (!is.na(best_pos)) break
  }
  # Nucleotide-codon fallback: if no AA protein motif matched (truncated,
  # divergent, or pseudogene J), locate the conserved anchor codon at the
  # nucleotide level — Trp (TGG) for IGH, Phe (TTT/TTC) for IGK/IGL. To avoid
  # grabbing a spurious downstream codon, REQUIRE the anchor codon to be
  # followed by a Glycine codon (GG[ACGT]) within one codon — reconstructing the
  # conserved [WF]-G context that the protein motif encodes. Among qualifying
  # matches, prefer the LAST (3'-most) one, which corresponds to the J-motif.
  if (is.na(best_pos)) {
    anchor_codons <- switch(chain_type, IGH = "TGG", IGK = c("TTT","TTC"),
                            IGL = c("TTT","TTC"), c("TGG","TTT","TTC"))
    hits <- integer(0L)
    for (ac in anchor_codons) {
      # anchor codon immediately followed by a Gly codon (GG.)
      m <- gregexpr(paste0(ac, "GG."), toupper(nt_seq))[[1L]]
      if (m[1L] > 0L) hits <- c(hits, m)
    }
    if (length(hits) == 0L) {
      # No [WF]-G context anywhere: fall back to the plain last anchor codon.
      codon_pat <- paste(anchor_codons, collapse = "|")
      m <- gregexpr(codon_pat, toupper(nt_seq))[[1L]]
      if (m[1L] > 0L) hits <- m
    }
    if (length(hits) > 0L) {
      hit_start  <- max(hits)   # 3'-most qualifying anchor codon
      best_pos   <- as.integer(hit_start - 2L)   # 0-based, -1 IMGT correction
      codon      <- toupper(substr(nt_seq, hit_start, hit_start + 2L))
      res        <- if (codon == "TGG") "W" else if (codon %in% c("TTT","TTC")) "F" else "?"
      best_motif <- paste0(res, "(", codon, ",nt-fallback)")
    }
  }
  # Attach the matched motif so callers can report it (diagnostic only).
  if (!is.na(best_pos)) attr(best_pos, "motif") <- best_motif
  best_pos
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Build all candidate names to try when looking up a J-gene anchor.
#' Handles: species-tag suffixes, allele stripping, legacy short names.
.aux_candidates <- function(gene_name) {
  # Build all name variants to try when looking up a J gene in the reference aux.
  # Reference aux uses both IMGT allele names (IGHJ1*01) and legacy short names (JH1).
  # Novel alleles may have a species-tag suffix: IGKJ1*02_mouse
  cands     <- gene_name                            # 1. exact
  no_tag    <- sub("_[^_*]+$", "", gene_name)       # 2. strip _species
  if (no_tag != gene_name) cands <- c(cands, no_tag)
  no_allele <- sub("\\*.*$", "", no_tag)           # 3. strip *allele
  if (no_allele != no_tag) cands <- c(cands, no_allele)
  # 4. Legacy short-form: IGHJ1->JH1, IGKJ1->JK1, IGLJ1->JL1
  short <- no_allele
  short <- sub("^IGHJ(\\d+)$", "JH\\1", short)
  short <- sub("^IGKJ(\\d+)$", "JK\\1", short)
  short <- sub("^IGLJ(\\d+)$", "JL\\1", short)
  if (short != no_allele) cands <- c(cands, short)
  unique(cands)
}

lookup_ref_anchor <- function(gene_name, ref_aux_dt, hmap_dt, annot_dt) {
  if (is.null(ref_aux_dt)) return(NULL)

  # TRIG-NC: the gene_name may be a colon-format name (IGKJ:01:001:001) whose
  # anchor lives in the reference aux under the original IMGT name (IGKJ1*01).
  # Resolve colon -> IMGT via the TRIG-NC map before candidate matching.
  resolved_names <- gene_name
  if (grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ]):[0-9]", gene_name) &&
      exists(".trignc_map_acc", envir = .GlobalEnv)) {
    tm <- get(".trignc_map_acc", envir = .GlobalEnv)
    # map: before (IMGT/original) -> after (colon). We want before for this after.
    origs <- unique(tm[after == gene_name, before])
    # Keep only IMGT-style origins (skip PIgLET F-names and other colon names)
    origs <- origs[grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ])[0-9]", origs)]
    if (length(origs) > 0L) resolved_names <- c(origs, gene_name)
  }

  # Try all candidate names derived from every resolved name. Track whether the
  # match was at the EXACT ALLELE level (name still contains *NN) or only a
  # coarser gene/legacy level, so the caller can prefer motif search over a
  # gene-level anchor that may belong to a different allele of the same gene.
  for (rn in resolved_names) {
    cands <- .aux_candidates(rn)
    for (ci in seq_along(cands)) {
      cand <- cands[ci]
      hit <- ref_aux_dt[gene == cand]
      if (nrow(hit) > 0L) {
        res <- copy(hit[1L])
        # exact-allele if the matched candidate retains an allele (*NN)
        res[, .match_level := if (grepl("\\*[0-9]", cand)) "allele" else "gene"]
        return(res)
      }
    }
  }
  # Via header normalisation map (normalised_name -> raw_header)
  if (!is.null(hmap_dt) && nrow(hmap_dt) > 0L) {
    for (rh in hmap_dt[normalised_name == gene_name, raw_header]) {
      for (cand in .aux_candidates(rh)) {
        hit <- ref_aux_dt[gene == cand]
        if (nrow(hit) > 0L) return(hit[1L])
      }
    }
  }
  # Via PIgLET annotation table (new_allele -> imgt_allele -> legacy name)
  if (!is.null(annot_dt) && nrow(annot_dt) > 0L) {
    for (ia in annot_dt[new_allele == gene_name, imgt_allele]) {
      for (cand in .aux_candidates(ia)) {
        hit <- ref_aux_dt[gene == cand]
        if (nrow(hit) > 0L) return(hit[1L])
      }
    }
  }
  NULL
}

# Modal extra_bps among reference J genes of a given chain (IGH/IGK/IGL).
# Used as the default for NOVEL J genes that have no reference row to inherit
# from, so they follow the same convention as the reference (typically 1).
.chain_default_extra_bps <- function(chain_type, ref_aux_dt) {
  jt <- switch(chain_type, IGH = "JH", IGK = "JK", IGL = "JL", NA_character_)
  if (!is.na(jt) && !is.null(ref_aux_dt) && "chain_type" %in% names(ref_aux_dt) &&
      "extra_bps" %in% names(ref_aux_dt)) {
    # Prefer allele-level entries (names with *) as the comparison class.
    vals <- ref_aux_dt[chain_type == jt & grepl("[*]", gene), extra_bps]
    vals <- vals[!is.na(vals)]
    if (!length(vals)) {
      vals <- ref_aux_dt[chain_type == jt, extra_bps]; vals <- vals[!is.na(vals)]
    }
    if (length(vals) > 0L) {
      tt <- sort(table(vals), decreasing = TRUE)
      return(as.integer(names(tt)[1L]))
    }
  }
  1L  # IG J genes conventionally carry extra_bps = 1
}

# ── Reference J sequence -> curated anchor map (for sequence-identity liftover)
# Builds a lookup from ungapped reference J sequence to its curated aux anchor,
# by joining reference J FASTAs (from the --species germlines) to the reference
# aux by gene name. Enables validating a NOVEL J gene's inferred anchor against
# the curated anchor of an IDENTICAL reference sequence, even when the names
# differ (e.g. OGRDB IGKJ0-XXXX*00 identical in sequence to IMGT IGKJ2*01).
.build_ref_jseq_anchor_map <- function(ref_j_gapped, ref_aux_dt) {
  if (is.null(ref_j_gapped) || length(ref_j_gapped) == 0L ||
      is.null(ref_aux_dt) || nrow(ref_aux_dt) == 0L)
    return(list())
  ung <- as.character(DECIPHER::RemoveGaps(ref_j_gapped, removeGaps = "all"))
  nms <- names(ref_j_gapped)
  m <- list()
  for (i in seq_along(ung)) {
    # find this ref gene's curated anchor via name candidates
    hit <- NULL
    for (cand in .aux_candidates(nms[i])) {
      h <- ref_aux_dt[gene == cand]
      if (nrow(h) > 0L) { hit <- h[1L]; break }
    }
    if (!is.null(hit)) {
      key <- toupper(ung[i])
      # keep first (allele-level preferred via .aux_candidates ordering)
      if (is.null(m[[key]]))
        m[[key]] <- list(stop = as.integer(hit$cdr3_stop),
                         frame = as.integer(hit$frame),
                         gene = nms[i])
    }
  }
  m
}

# ── Shared J-gene aux coordinate derivation ─────────────────────────────────
# Single source of truth for how a J gene's aux entry (CDR3-stop anchor, coding
# frame offset, extra_bps) is derived. Used by BOTH the --as-is-ids path and the
# --asc / --trig_nc path so the coordinates are always computed identically.
#
# Priority:
#   1. Exact-allele reference match  -> authoritative (per-allele cdr3_stop/frame)
#   2. Motif search on the sequence  -> preferred over a gene-level ref match
#   3. Gene-level reference match     -> last resort
#   4. Nothing found                  -> anchor 0, frame 0
#
# Relationships used for motif-derived (novel) genes, validated against the IMGT
# reference aux (holds for all mouse J alleles):
#   frame     = (cdr3_stop - 2) mod 3
#   extra_bps = modal extra_bps among allele-level reference J genes of the chain
#
# Returns a list: anchor, frame, extra_bps, method, motif_anchor, agrees
.derive_j_aux <- function(gene, nt, chain_type, ref_aux_dt, hmap_dt = NULL,
                          annot_dt = NULL, ref_seq_map = NULL) {
  anchor <- NA_integer_; frame <- 0L; method <- "none"; extra <- NA_integer_
  motif_anchor <- find_j_anchor(nt, chain_type)
  motif_str    <- if (!is.na(motif_anchor)) attr(motif_anchor, "motif") else NA_character_

  # ── Priority 0: exact sequence-identity match to a reference J sequence ──
  # If this J sequence is byte-identical to a reference sequence, that
  # reference's CURATED anchor is authoritative by definition (same sequence).
  # This beats motif inference — especially for pseudogenes / atypical J where
  # motif inference is unreliable. Captured up-front so it takes precedence.
  seqmatch_stop <- NA_integer_; seqmatch_gene <- NA_character_
  if (!is.null(ref_seq_map)) {
    rm <- ref_seq_map[[toupper(nt)]]
    if (!is.null(rm)) { seqmatch_stop <- rm$stop; seqmatch_gene <- rm$gene }
  }
  if (!is.na(seqmatch_stop)) {
    anchor <- as.integer(seqmatch_stop)
    frame  <- ((anchor - 2L) %% 3L + 3L) %% 3L
    extra  <- .chain_default_extra_bps(chain_type, ref_aux_dt)
    # keep any curated frame/extra from the matched reference row if available
    hitrow <- if (!is.na(seqmatch_gene)) {
      h <- NULL
      for (cand in .aux_candidates(seqmatch_gene)) {
        hh <- ref_aux_dt[gene == cand]; if (nrow(hh) > 0L) { h <- hh[1L]; break }
      }; h
    } else NULL
    if (!is.null(hitrow)) {
      frame <- as.integer(hitrow$frame)
      if ("extra_bps" %in% names(hitrow)) extra <- as.integer(hitrow$extra_bps)
    }
    method <- "seq_identity"
    return(list(anchor = anchor, frame = as.integer(frame),
                extra_bps = as.integer(extra), method = method, motif = motif_str,
                motif_anchor = if (is.na(motif_anchor)) NA_integer_ else as.integer(motif_anchor),
                agrees = if (!is.na(motif_anchor)) abs(motif_anchor - anchor) <= 3L else NA,
                seqmatch_stop = seqmatch_stop, seqmatch_gene = seqmatch_gene))
  }

  rh <- lookup_ref_anchor(gene, ref_aux_dt, hmap_dt, annot_dt)
  match_level <- if (!is.null(rh) && ".match_level" %in% names(rh))
                   rh$.match_level else if (!is.null(rh)) "allele" else NA_character_

  if (!is.null(rh) && identical(match_level, "allele")) {
    anchor <- as.integer(rh$cdr3_stop)
    frame  <- as.integer(rh$frame)
    if ("extra_bps" %in% names(rh)) extra <- as.integer(rh$extra_bps)
    method <- "reference_aux(allele)"
    if (!is.na(motif_anchor) && abs(motif_anchor - anchor) > 3L)
      message(sprintf("    [CHECK] %s (%s): allele anchor=%d, motif=%d (Δ=%d)",
                      gene, chain_type, anchor, motif_anchor, motif_anchor - anchor))
  } else if (!is.na(motif_anchor)) {
    anchor <- motif_anchor
    frame  <- ((anchor - 2L) %% 3L + 3L) %% 3L
    method <- "motif_search"
    if (!is.null(rh) && "extra_bps" %in% names(rh)) extra <- as.integer(rh$extra_bps)
  } else if (!is.null(rh)) {
    anchor <- as.integer(rh$cdr3_stop)
    frame  <- as.integer(rh$frame)
    if ("extra_bps" %in% names(rh)) extra <- as.integer(rh$extra_bps)
    method <- "reference_aux(gene)"
  } else {
    message(sprintf("    [WARN] No anchor for %s (%s); defaulting 0", gene, chain_type))
    anchor <- 0L; frame <- 0L; method <- "default"
  }

  if (is.na(extra)) extra <- .chain_default_extra_bps(chain_type, ref_aux_dt)

  # Consistent per-gene logging for motif-derived (novel) J genes, used by BOTH
  # the as-is and ASC/TRIG-NC paths.
  if (method == "motif_search") {
    cat(sprintf("  [AUX-motif] %s: motif '%s' -> stop=%d frame=%d extra_bps=%d\n",
                gene, if (is.na(motif_str)) "?" else motif_str,
                as.integer(anchor), as.integer(frame), as.integer(extra)))
  }

  list(anchor = as.integer(anchor), frame = as.integer(frame),
       extra_bps = as.integer(extra), method = method,
       motif = motif_str,
       motif_anchor = if (is.na(motif_anchor)) NA_integer_ else as.integer(motif_anchor),
       agrees = if (grepl("^reference_aux", method) && !is.na(motif_anchor))
                  abs(motif_anchor - anchor) <= 3L else NA,
       seqmatch_stop = seqmatch_stop, seqmatch_gene = seqmatch_gene)
}

# If --as_is_ids, run the lightweight pipeline and exit early
if (isTRUE(opt$as_is_ids)) {
  # Helper: truncate names to BLAST's 50-char local id limit
  .truncate_id <- function(nm, max_len = 50L) {
    ifelse(nchar(nm) > max_len, substr(nm, 1L, max_len), nm)
  }

  # Helper: deduplicate names with _1, _2 suffix
  .dedup_names <- function(nms) {
    suffix_count <- list()
    for (i in seq_along(nms)) {
      nm <- nms[i]
      if (!is.null(suffix_count[[nm]])) {
        nms[i] <- paste0(nm, "_", suffix_count[[nm]])
        suffix_count[[nm]] <- suffix_count[[nm]] + 1L
      } else {
        suffix_count[[nm]] <- 1L
      }
    }
    nms
  }



  source_as_is <- function() {
    cat("\n=== as-is-ids mode: input names used directly, merged with reference ===\n")
    cat(sprintf("Input dir  : %s\n", opt$custom_dir))
    cat(sprintf("Ref dir    : %s\n", opt$ref_dir))
    cat(sprintf("Organism   : %s\n", file_prefix))

    loci <- c("IGHV","IGHD","IGHJ","IGKV","IGKJ","IGLV","IGLJ")
    # Constant loci are handled separately — they are copied from reference
    # with name deduplication but no seq-dedup (all constant alleles are kept).
    const_loci <- c("IGHC","IGKC","IGLC")

    for (locus in loci) {
      # ── Find custom input FASTA ────────────────────────────────────────
      chain_sub <- if (startsWith(locus, "IGH")) "heavy" else "light"
      cust_fa   <- Filter(file.exists, c(
        file.path(opt$custom_dir, chain_sub, paste0(locus, ".fasta")),
        file.path(opt$custom_dir, paste0(locus, ".fasta"))))[1L]

      # ── Find reference FASTA ───────────────────────────────────────────
      ref_fa <- Filter(file.exists, c(
        file.path(opt$ref_dir, paste0("imgt_", opt$species, "_", locus, ".fasta")),
        file.path(opt$ref_dir, paste0(opt$species, "_", locus, ".fasta")),
        file.path(opt$ref_dir, paste0(locus, ".fasta"))))[1L]

      # Skip locus if neither custom nor reference exists
      if (is.na(cust_fa) && is.na(ref_fa)) {
        cat(sprintf("  [SKIP] %s: no custom or reference FASTA found\n", locus)); next
      }

      # ── Load and merge: custom first, then reference ───────────────────
      cust_seqs <- if (!is.na(cust_fa))
        Biostrings::readDNAStringSet(cust_fa) else Biostrings::DNAStringSet()
      ref_seqs  <- if (!is.na(ref_fa))
        Biostrings::readDNAStringSet(ref_fa)  else Biostrings::DNAStringSet()

      n_cust <- length(cust_seqs); n_ref <- length(ref_seqs)
      seqs   <- c(cust_seqs, ref_seqs)
      n_raw  <- length(seqs)

      cat(sprintf("  %s: %d custom + %d reference = %d total\n",
                  locus, n_cust, n_ref, n_raw))

      if (n_raw == 0L) { cat(sprintf("  [SKIP] %s: empty after merge\n", locus)); next }

      # ── Step 1: parse and normalise FASTA headers ──────────────────────
      # Three header formats are handled:
      #   (a) IMGT pipe-delimited: acc|GENE*ALLELE|Species_strain|func|...
      #       -> GENE*ALLELE_strain  e.g. IGHV1-18-28*01_BALB/cJ
      #   (b) Plain allele with tag: GENE*ALLELE_tag  (kept as-is)
      #   (c) Plain allele: GENE*ALLELE               (kept as-is)
      # All names are truncated to the BLAST 50-char local id limit after parsing.
      .parse_imgt_header <- function(h) {
        h <- sub("\\s.*$", "", h)   # strip trailing description first
        if (grepl("|", h, fixed = TRUE)) {
          parts <- strsplit(h, "\\|")[[1L]]
          gene_allele <- if (length(parts) >= 2L) parts[2L] else h
          strain <- ""
          if (length(parts) >= 3L) {
            sp_parts <- strsplit(parts[3L], "_")[[1L]]
            if (length(sp_parts) >= 3L) {
              strain <- paste(sp_parts[seq(3L, length(sp_parts))], collapse="_")
            } else if (length(sp_parts) == 2L) {
              strain <- sp_parts[2L]
            }
          }
          if (nchar(strain) > 0L) paste0(gene_allele, "_", strain) else gene_allele
        } else {
          h   # plain format — keep as-is
        }
      }
      nms <- vapply(names(seqs), .parse_imgt_header, character(1L),
                    USE.NAMES = FALSE)
      too_long <- nchar(nms) > 50L
      if (any(too_long)) {
        cat(sprintf("  [WARN] %s: %d name(s) exceed 50 chars and will be truncated\n",
                    locus, sum(too_long)))
        nms <- .truncate_id(nms)

      }
      names(seqs) <- nms  # apply parsed/truncated names before dedup
      # ── Step 2: sequence deduplication with priority rules ──────────────
      # Priority: OGRDB/custom > IMGT reference
      # Rules applied in order:
      #  a) OGRDB sequences are identified (non-IMGT naming, hex-suffix after rename)
      #  b) Any IMGT sequence identical to an OGRDB sequence is dropped (cross-gene)
      #     The OGRDB sequence is the MRL-specific allele we want to retain
      #  c) Among remaining OGRDB sequences, keep unique sequences only (cross-gene)
      #     Different OGRDB IDs may represent the same sequence
      #  d) Among remaining IMGT sequences, gene-level dedup: drop sequences
      #     identical to another allele within the SAME gene family only
      #     (preserves cross-gene identical sequences as biologically meaningful)

      seqs_char <- as.character(seqs)
      n_before  <- length(seqs)

      # Classify each sequence as OGRDB (custom novel) or IMGT (reference)
      # OGRDB alleles after renaming: IGxV[A-Z]{4}*00 pattern
      # OGRDB novel alleles use family "0-" placeholder: IGKV0-XXXX*00, IGKJ0-XXXX*00
      is_ogrdb <- grepl("^IG[HKL][VDJ]0-[A-Z0-9]", names(seqs))
      is_imgt  <- !is_ogrdb
      ogrdb_seqs <- seqs_char[is_ogrdb]
      imgt_seqs  <- seqs_char[is_imgt]

      # (b) Drop IMGT sequences identical to any OGRDB sequence (cross-gene)
      # EXCEPTION: J genes — always retain IMGT J entries even when identical
      # to an OGRDB entry. igblastn reports the J gene name from the database
      # and uses it for aux CDR3 anchor lookup. Losing standard IGKJ1-5 entries
      # means sequences that match those J genes get no J call or wrong anchor.
      is_j_gene <- grepl("^IG[HKL]J", names(seqs))
      imgt_drop_cg <- is_imgt & (seqs_char %in% ogrdb_seqs) & !is_j_gene
      n_imgt_cg_dropped <- sum(imgt_drop_cg)

      # (c) Among OGRDB: drop cross-gene duplicates (keep first = custom priority)
      ogrdb_dup <- is_ogrdb & duplicated(seqs_char)
      n_ogrdb_dup <- sum(ogrdb_dup)

      # (d) Among remaining IMGT: gene-level dedup only
      # J genes are exempt from CROSS-gene dedup (retain OGRDB + IMGT J alleles),
      # but exact duplicates — same sequence AND same name — are still removed
      # so we never end up with IGKJ1*01 and IGKJ1*01_1 that are byte-identical.
      gene_base <- sub("\\*.*$", "", names(seqs))  # strip allele
      gene_base <- sub("_[^_*]+$", "", gene_base)   # strip strain tag
      imgt_gene_dup <- logical(length(seqs))
      seen_per_gene <- list()
      for (k in seq_along(seqs)) {
        if (!is_imgt[k] || imgt_drop_cg[k] || ogrdb_dup[k]) next
        if (is_j_gene[k]) next  # J genes handled below (cross-gene exempt)
        gene_k <- gene_base[k]
        seq_k  <- seqs_char[k]
        if (is.null(seen_per_gene[[gene_k]])) {
          seen_per_gene[[gene_k]] <- seq_k
        } else if (seq_k %in% seen_per_gene[[gene_k]]) {
          imgt_gene_dup[k] <- TRUE
        } else {
          seen_per_gene[[gene_k]] <- c(seen_per_gene[[gene_k]], seq_k)
        }
      }

      # (e) Exact-duplicate removal for ALL sequences (incl. J genes):
      # collapse entries where BOTH the name and the sequence are identical.
      # This is the sequence-level dedup that prevents _1-suffixed clones of
      # byte-identical J sequences. Distinct-name identical-sequence J entries
      # (e.g. OGRDB IGKJ0-X vs IMGT IGKJ1*01) are preserved by design.
      name_seq_key <- paste0(names(seqs), "\u0001", seqs_char)
      exact_dup <- duplicated(name_seq_key)
      n_exact_dup <- sum(exact_dup)

      keep_seq <- !(imgt_drop_cg | ogrdb_dup | imgt_gene_dup | exact_dup)
      n_dropped <- sum(!keep_seq)

      if (n_dropped > 0L) {
        if (n_imgt_cg_dropped > 0L)
          cat(sprintf("  [DEDUP] %s: %d IMGT sequence(s) replaced by identical OGRDB allele(s)\n",
                      locus, n_imgt_cg_dropped))
        if (n_ogrdb_dup > 0L)
          cat(sprintf("  [DEDUP] %s: %d OGRDB cross-gene duplicate(s) removed\n",
                      locus, n_ogrdb_dup))
        n_gene_dup <- sum(imgt_gene_dup)
        if (n_gene_dup > 0L)
          cat(sprintf("  [DEDUP] %s: %d within-gene IMGT duplicate(s) removed\n",
                      locus, n_gene_dup))
        if (n_exact_dup > 0L)
          cat(sprintf("  [DEDUP] %s: %d exact duplicate(s) removed (same name + sequence)\n",
                      locus, n_exact_dup))
        seqs <- seqs[keep_seq]
        nms  <- names(seqs)
      }

      # ── Step 3: name-level deduplication ──────────────────────────────
      n_name_dup <- sum(duplicated(nms))
      if (n_name_dup > 0L) {
        nms <- .dedup_names(nms)
        names(seqs) <- nms
        cat(sprintf("  [DEDUP-name] %s: %d duplicate name(s) suffixed _1, _2 ...\n",
                    locus, n_name_dup))
      }

      cat(sprintf("  [OK] %s: %d final unique sequences\n", locus, length(seqs)))

      # ── Write gapped FASTA ─────────────────────────────────────────────
      gapped_dir <- file.path(opt$outdir, "germlines", "gapped")
      dir.create(gapped_dir, recursive=TRUE, showWarnings=FALSE)
      out_fa <- file.path(gapped_dir, paste0(imgt_file_prefix, "_", locus, ".fasta"))
      Biostrings::writeXStringSet(seqs, out_fa)
    }
    # ── Constant region FASTAs: deduplicate by name, write to gapped/ ───────
    cat("\n--- Constant region sequences (from reference, name-deduped) ---\n")
    for (locus in const_loci) {
      # Constant FASTAs live in constant/ not vdj/ — search both the constant
      # subdirectory and common IGDATA germline locations
      ref_const_dir <- sub("/vdj$", "/constant", opt$ref_dir)
      ref_fa <- Filter(file.exists, c(
        file.path(ref_const_dir, paste0("imgt_", opt$species, "_", locus, ".fasta")),
        file.path(ref_const_dir, paste0(opt$species, "_", locus, ".fasta")),
        file.path(ref_const_dir, paste0(locus, ".fasta")),
        # Also try vdj/ as fallback in case constant/ isn't separate
        file.path(opt$ref_dir, paste0("imgt_", opt$species, "_", locus, ".fasta")),
        file.path(opt$ref_dir, paste0(locus, ".fasta"))))[1L]
      if (is.na(ref_fa)) {
        cat(sprintf("  [SKIP] %s: no reference FASTA found\n", locus)); next
      }
      seqs <- Biostrings::readDNAStringSet(ref_fa)
      if (length(seqs) == 0L) {
        cat(sprintf("  [SKIP] %s: empty FASTA\n", locus)); next
      }
      # Strip IMGT pipe-delimited headers to plain allele names
      nms <- vapply(names(seqs), .parse_imgt_header, character(1L),
                    USE.NAMES=FALSE)
      nms <- .truncate_id(nms)
      names(seqs) <- nms
      # Name-level deduplication only (keep all constant alleles, just remove dups)
      n_before <- length(seqs)
      keep <- !duplicated(nms)
      seqs <- seqs[keep]
      n_dup <- n_before - length(seqs)
      cat(sprintf("  %s: %d sequences", locus, length(seqs)))
      if (n_dup > 0L) cat(sprintf(" (%d duplicate names removed)", n_dup))
      cat("\n")
      gapped_dir <- file.path(opt$outdir, "germlines", "gapped")
      dir.create(gapped_dir, recursive=TRUE, showWarnings=FALSE)
      out_fa <- file.path(gapped_dir, paste0(imgt_file_prefix, "_", locus, ".fasta"))
      Biostrings::writeXStringSet(seqs, out_fa)
    }

    cat("\n=== as-is-ids annotation complete ===\n")

    # ── Build aux file from as-is J gene sequences ──────────────────────
    cat("\n--- Building auxiliary file (as-is-ids mode) ---\n")
    .jaux_val <- list()   # accumulate motif-inference validation rows
    aux_dir <- file.path(opt$outdir, "auxiliary")
    dir.create(aux_dir, recursive=TRUE, showWarnings=FALSE)
    aux_path_ai <- file.path(aux_dir, paste0(file_prefix, "_gl.aux"))

    # Load reference aux for anchor lookup
    ref_aux_dt_ai <- NULL
    for (cand in c(
        file.path(opt$igdata, "optional_file", paste0(opt$species, "_gl.aux")),
        file.path(opt$igdata, "optional_file", paste0(opt$species, ".aux")))) {
      if (file.exists(cand)) {
        aux_raw <- readLines(cand)
        aux_data <- aux_raw[!grepl("^\\s*#", aux_raw) & nchar(trimws(aux_raw)) > 0L]
        aux_split <- strsplit(trimws(aux_data), "\\s+")
        aux_rows_ok <- aux_split[sapply(aux_split, length) == 5L]
        if (length(aux_rows_ok) > 0L) {
          ref_aux_dt_ai <- as.data.table(do.call(rbind, aux_rows_ok))
          setnames(ref_aux_dt_ai, c("gene","frame","chain_type","cdr3_stop","extra_bps"))
          ref_aux_dt_ai[, c("frame","cdr3_stop","extra_bps") :=
            .(as.integer(frame), as.integer(cdr3_stop), as.integer(extra_bps))]
          cat(sprintf("  Loaded reference aux: %s (%d entries)\n",
                      cand, nrow(ref_aux_dt_ai)))
          break
        }
      }
    }

    # Build reference-J sequence -> curated anchor map for sequence-identity
    # liftover validation (now that ref_aux_dt_ai is loaded).
    .ref_jseq_map <- {
      rj <- Biostrings::DNAStringSet()
      for (jl in c("IGHJ","IGKJ","IGLJ")) {
        rjfa <- Filter(file.exists, c(
          file.path(opt$ref_dir, paste0("imgt_", opt$species, "_", jl, ".fasta")),
          file.path(opt$ref_dir, paste0(opt$species, "_", jl, ".fasta")),
          file.path(opt$ref_dir, paste0(jl, ".fasta"))))[1L]
        if (!is.na(rjfa)) {
          s <- tryCatch(Biostrings::readDNAStringSet(rjfa),
                        error = function(e) Biostrings::DNAStringSet())
          if (length(s) > 0L) {
            names(s) <- vapply(names(s), function(h) {
              if (grepl("\\|", h)) { f <- strsplit(h, "\\|")[[1L]]
                if (length(f) >= 2L && nzchar(f[2])) f[2] else h } else sub("\\s.*$","",h)
            }, character(1L))
            rj <- c(rj, s)
          }
        }
      }
      if (length(rj) > 0L && !is.null(ref_aux_dt_ai))
        .build_ref_jseq_anchor_map(rj, ref_aux_dt_ai) else list()
    }

    con_aux <- file(aux_path_ai, open="wt")
    writeLines(c(
      "#gene/allele name, first coding frame start position, chain type, CDR3 stop, extra bps beyond J coding end.",
      "#All positions are 0-based", ""), con_aux)

    # Copy legacy short names from reference
    if (!is.null(ref_aux_dt_ai)) {
      legacy_ai <- ref_aux_dt_ai[!grepl("\\*", gene) & grepl("^J[HKLA-Z]\\d", gene)]
      for (r in seq_len(nrow(legacy_ai)))
        writeLines(paste(legacy_ai$gene[r], legacy_ai$frame[r], legacy_ai$chain_type[r],
                         legacy_ai$cdr3_stop[r], legacy_ai$extra_bps[r], sep="\t"), con_aux)
      if (nrow(legacy_ai)) writeLines("", con_aux)
    }

    chain_map_ai <- c(IGH="JH", IGK="JK", IGL="JL")
    for (j_locus in c("IGHJ","IGKJ","IGLJ")) {
      j_fa <- Filter(file.exists, c(
        file.path(opt$outdir, "germlines", "gapped",
                  paste0(imgt_file_prefix, "_", j_locus, ".fasta"))))[1L]
      if (is.na(j_fa)) next
      j_seqs <- Biostrings::readDNAStringSet(j_fa)
      # Strip gap characters (. - =) inline — no external function needed
      j_ung   <- gsub("[.=-]", "", as.character(j_seqs))
      names(j_ung) <- names(j_seqs)
      ct_key  <- sub("J$","",j_locus)            # IGHJ->IGH etc.
      ct_out  <- chain_map_ai[ct_key]
      written_j <- character(0L)
      for (i in seq_along(j_ung)) {
        nm  <- names(j_ung)[i]
        nt  <- j_ung[i]

        # ── Derive anchor / frame / extra_bps via the SHARED helper ────────
        # Same logic used by the ASC/TRIG-NC path: exact-allele reference wins,
        # else motif search, else gene-level reference. Guarantees identical
        # coordinate derivation regardless of naming mode.
        .d <- .derive_j_aux(nm, nt, ct_key, ref_aux_dt_ai,
                            hmap_dt = NULL, annot_dt = NULL,
                            ref_seq_map = .ref_jseq_map)
        anchor_val <- .d$anchor
        frame_val  <- .d$frame
        extra_val  <- .d$extra_bps
        # (.derive_j_aux already logs [AUX-motif] for novel/motif-derived genes)

        # Accumulate motif-inference validation: for reference-lifted anchors,
        # record the independent motif result on the same sequence for later
        # comparison (validates the sequence-only method against ground truth).
        # Validation row: prefer the sequence-identity reference anchor as the
        # ground truth when available (works for novel names too); else use the
        # name-lifted reference anchor. Compare our final anchor + the motif
        # inference against it.
        gt_stop <- if (!is.na(.d$seqmatch_stop)) .d$seqmatch_stop
                   else if (grepl("^reference_aux", .d$method)) as.integer(anchor_val)
                   else NA_integer_
        if (!is.na(gt_stop)) {
          .jaux_val[[length(.jaux_val) + 1L]] <- data.table(
            gene = nm, chain = ct_out, anchor_method = .d$method,
            ground_truth = if (!is.na(.d$seqmatch_stop)) "seq_identity" else "name_lift",
            gt_gene = if (!is.na(.d$seqmatch_gene)) .d$seqmatch_gene else nm,
            reference_stop = gt_stop,
            our_stop   = as.integer(anchor_val),
            motif_stop = as.integer(.d$motif_anchor),
            delta_our   = as.integer(anchor_val - gt_stop),
            delta_motif = if (is.na(.d$motif_anchor)) NA_integer_
                          else as.integer(.d$motif_anchor - gt_stop),
            agrees_pm1 = abs(anchor_val - gt_stop) <= 1L,
            motif = if (is.na(.d$motif)) NA_character_ else .d$motif,
            sequence = nt)
        }

        writeLines(paste(nm, frame_val, ct_out,
                         anchor_val, extra_val, sep="\t"), con_aux)
        written_j <- c(written_j, nm)

        # ── Base-name alias (strip strain tag) ────────────────────────────
        base_nm <- sub("(\\*\\d+)_.*$", "\\1", nm)
        if (base_nm != nm && !base_nm %in% written_j) {
          writeLines(paste(base_nm, frame_val, ct_out,
                           anchor_val, extra_val, sep="\t"), con_aux)
          written_j <- c(written_j, base_nm)
        }
      }
    }
    # ── Append any reference J entries not already written ──────────────
    # Sequence deduplication in source_as_is() removes J sequences that are
    # identical to reference sequences, so those J genes never appear in the
    # hybrid gapped FASTA and would be missing from the aux. Copy them from
    # the reference aux so igblastn can find CDR3 anchors for all J alleles.
    if (!is.null(ref_aux_dt_ai)) {
      ref_j_entries <- ref_aux_dt_ai[grepl("\\*", gene)]  # IMGT allele entries only
      n_added <- 0L
      for (r in seq_len(nrow(ref_j_entries))) {
        nm_r <- ref_j_entries$gene[r]
        # Normalise: strip strain tag for comparison
        base_r <- sub("(\\*\\d+)_.*$", "\\1", nm_r)
        if (!nm_r %in% written_j && !base_r %in% written_j) {
          writeLines(paste(nm_r,
                           ref_j_entries$frame[r],
                           ref_j_entries$chain_type[r],
                           ref_j_entries$cdr3_stop[r],
                           ref_j_entries$extra_bps[r], sep="\t"), con_aux)
          written_j <- c(written_j, nm_r, base_r)
          n_added <- n_added + 1L
        }
      }
      if (n_added > 0L)
        cat(sprintf("  [AUX] Appended %d reference J entries not in hybrid FASTAs\n",
                    n_added))
    }
    close(con_aux)
    cat(sprintf("  Wrote aux: %s (%d total entries)\n",
                aux_path_ai, length(written_j)))

    # ── J-anchor validation report (sequence-identity ground truth) ──────
    if (length(.jaux_val) > 0L) {
      vdt <- rbindlist(.jaux_val, fill = TRUE)
      n_chk <- nrow(vdt)
      n_seq <- sum(vdt$ground_truth == "seq_identity", na.rm = TRUE)
      n_our_ok <- sum(vdt$agrees_pm1, na.rm = TRUE)
      n_our_ex <- sum(vdt$delta_our == 0L, na.rm = TRUE)
      # How well does the sequence-only motif method reproduce ground truth?
      mvals <- vdt$delta_motif[!is.na(vdt$delta_motif)]
      n_motif_ok <- sum(abs(mvals) <= 1L); n_motif_ex <- sum(mvals == 0L)
      cat(sprintf("  J-anchor validation vs reference (%d genes; %d via exact sequence identity):\n",
                  n_chk, n_seq))
      cat(sprintf("    final anchor  : %d/%d within ±1nt, %d exact\n",
                  n_our_ok, n_chk, n_our_ex))
      cat(sprintf("    motif-only    : %d/%d within ±1nt, %d exact\n",
                  n_motif_ok, length(mvals), n_motif_ex))
      disc <- vdt[abs(delta_our) > 0L][order(-abs(delta_our))]
      if (nrow(disc) > 0L) {
        cat(sprintf("    %d final-anchor discrepancy(ies):\n", nrow(disc)))
        for (r in seq_len(min(nrow(disc), 20L)))
          cat(sprintf("      %-22s our=%d vs %s '%s'=%d  Δ=%+d motif='%s'\n",
                      disc$gene[r], disc$our_stop[r], disc$ground_truth[r],
                      disc$gt_gene[r], disc$reference_stop[r], disc$delta_our[r],
                      if (is.na(disc$motif[r])) "?" else disc$motif[r]))
      }
      val_out <- file.path(opt$outdir, "annotations",
                           paste0(opt$prefix, "_jaux_validation.tsv"))
      dir.create(dirname(val_out), recursive = TRUE, showWarnings = FALSE)
      fwrite(vdt[order(-abs(delta_our))], val_out, sep = "\t")
      cat(sprintf("  Wrote J-anchor validation table: %s (%d rows)\n",
                  val_out, nrow(vdt)))
    }

    # ── Build ndm.imgt from as-is V gene sequences ───────────────────────
    cat("\n--- Building ndm.imgt (as-is-ids mode) ---\n")
    ndm_path_ai <- file.path(aux_dir, paste0(file_prefix, ".ndm.imgt"))

    gapped_dir_ai <- file.path(opt$outdir, "germlines", "gapped")
    chain_ndm <- c(IGHV="VH", IGKV="VK", IGLV="VL")
    fwr1_end  <- c(IGHV=75L, IGKV=78L, IGLV=78L)

    first_nt_from <- function(chars, gpos) {
      for (j in seq(gpos, length(chars)))
        if (!chars[j] %in% c(".","=","-")) return(j)
      NA_integer_
    }
    last_nt_to <- function(chars, gpos) {
      gpos <- min(gpos, length(chars))
      for (j in seq(gpos, 1L, -1L))
        if (!chars[j] %in% c(".","=","-")) return(j)
      NA_integer_
    }

    ndm_rows_ai <- rbindlist(lapply(names(chain_ndm), function(v_locus) {
      v_fa <- file.path(gapped_dir_ai,
                        paste0(imgt_file_prefix, "_", v_locus, ".fasta"))
      if (!file.exists(v_fa)) return(data.table())
      v_seqs <- Biostrings::readDNAStringSet(v_fa)
      gene_bases <- sub("[*].*$", "", names(v_seqs))
      v_seqs <- v_seqs[!duplicated(gene_bases)]
      ct <- chain_ndm[v_locus]
      fe <- fwr1_end[v_locus]
      rbindlist(lapply(seq_along(v_seqs), function(i) {
        chars <- strsplit(as.character(v_seqs[[i]]),"")[[1L]]
        n     <- length(chars)
        # Convert gapped character index -> ungapped nucleotide count.
        # igblastn ndm.imgt requires UNGAPPED positions.
        .ung_ai <- function(gp) {
          if (is.na(gp) || gp < 1L) return(gp)
          sum(!chars[seq_len(min(gp, n))] %in% c(".","-","="))
        }
        clamp <- function(sg, eg) {
          s <- if (is.na(sg) || sg < 0L) -1L else as.integer(.ung_ai(sg))
          e <- if (is.na(eg) || eg < 0L) -1L else as.integer(.ung_ai(eg))
          if (s != -1L && e != -1L && s > e) { s <- -1L; e <- -1L }
          c(s, e)
        }
        # Compute boundaries in GAPPED coords, convert to UNGAPPED via clamp()
        fwr1_g <- c(first_nt_from(chars,1L), last_nt_to(chars,min(fe,n)))
        fwr1   <- clamp(fwr1_g[1L], fwr1_g[2L])
        cdr1_g <- c(first_nt_from(chars,min(fe+1L,n)), last_nt_to(chars,min(114L,n)))
        cdr1   <- clamp(cdr1_g[1L], cdr1_g[2L])
        fwr2_ss <- if(!is.na(cdr1_g[2L])&&cdr1_g[2L]>0L) cdr1_g[2L]+1L else 115L
        fwr2_g  <- c(first_nt_from(chars,min(fwr2_ss,n)), last_nt_to(chars,min(165L,n)))
        fwr2    <- clamp(fwr2_g[1L], fwr2_g[2L])
        cdr2_ss <- if(!is.na(fwr2_g[2L])&&fwr2_g[2L]>0L) fwr2_g[2L]+1L else 166L
        cdr2_g  <- c(first_nt_from(chars,min(cdr2_ss,n)), last_nt_to(chars,min(195L,n)))
        cdr2    <- clamp(cdr2_g[1L], cdr2_g[2L])
        fwr3_ss <- if(!is.na(cdr2_g[2L])&&cdr2_g[2L]>0L) cdr2_g[2L]+1L else 196L
        fwr3_g  <- c(first_nt_from(chars,min(fwr3_ss,n)), last_nt_to(chars,min(312L,n)))
        fwr3    <- clamp(fwr3_g[1L], fwr3_g[2L])
        nm <- names(v_seqs)[i]
        base_nm <- sub("(\\*\\d+)_.*$","\\1", nm)
        row1 <- data.table(gene=nm,
          fwr1_start=fwr1[1L],fwr1_stop=fwr1[2L],
          cdr1_start=cdr1[1L],cdr1_stop=cdr1[2L],
          fwr2_start=fwr2[1L],fwr2_stop=fwr2[2L],
          cdr2_start=cdr2[1L],cdr2_stop=cdr2[2L],
          fwr3_start=fwr3[1L],fwr3_stop=fwr3[2L],
          chain_type=ct, trailing=0L)
        if (base_nm != nm) rbindlist(list(row1, {r2<-copy(row1);r2$gene<-base_nm;r2}))
        else row1
      }))
    }), fill=TRUE)

    if (nrow(ndm_rows_ai) > 0L) {
      pos_cols <- setdiff(names(ndm_rows_ai), c("gene","chain_type","trailing"))
      for (col in pos_cols) {
        set(ndm_rows_ai, NULL, col, as.integer(ndm_rows_ai[[col]]))
        set(ndm_rows_ai, which(is.na(ndm_rows_ai[[col]])), col, -1L)
      }
      con_ndm <- file(ndm_path_ai, open="wt")
      for (r in seq_len(nrow(ndm_rows_ai))) {
        pos_vals <- as.integer(unlist(ndm_rows_ai[r, pos_cols, with=FALSE]))
        pos_vals[is.na(pos_vals)] <- -1L
        if (length(pos_vals)==10L)
          writeLines(paste(c(ndm_rows_ai$gene[r], pos_vals,
                             ndm_rows_ai$chain_type[r], ndm_rows_ai$trailing[r]),
                           collapse="\t"), con_ndm)
      }
      # Copy to internal_data/
      int_dir <- file.path(opt$outdir, "internal_data", opt$species)
      if (dir.exists(int_dir)) {
        int_nm <- file.path(int_dir, paste0(opt$species, ".ndm.imgt"))
        file.copy(ndm_path_ai, int_nm, overwrite=TRUE)
      }
      # ── Append reference V entries missing from hybrid ndm ──────────────
      # V genes removed by sequence dedup won't have ndm entries. Read the
      # reference ndm and append any gene bases not already covered.
      ref_ndm_path <- NULL
      for (cand in c(
          file.path(opt$igdata, "internal_data", opt$species,
                    paste0(opt$species, ".ndm.imgt")),
          file.path(dirname(opt$igdata), "internal_data", opt$species,
                    paste0(opt$species, ".ndm.imgt")))) {
        if (file.exists(cand)) { ref_ndm_path <- cand; break }
      }
      if (!is.null(ref_ndm_path)) {
        ref_ndm_lines <- readLines(ref_ndm_path)
        # Gene bases already written
        written_genes <- sub("[*\t].*$", "", ndm_rows_ai$gene)
        n_ref_added <- 0L
        for (ln in ref_ndm_lines) {
          gene_nm <- strsplit(ln, "\t")[[1L]][1L]
          gene_base <- sub("[*].*$", "", gene_nm)
          if (!gene_base %in% written_genes && !is.na(gene_base) && nchar(gene_base) > 0L) {
            writeLines(ln, con_ndm)
            written_genes <- c(written_genes, gene_base)
            n_ref_added <- n_ref_added + 1L
          }
        }
        if (n_ref_added > 0L)
          cat(sprintf("  [NDM] Appended %d reference V entries not in hybrid FASTAs\n",
                      n_ref_added))
      }
      close(con_ndm)
      cat(sprintf("  Wrote ndm.imgt: %s (%d V genes)\n",
                  ndm_path_ai, nrow(ndm_rows_ai)))
    }

    # ── Write combined ALL_V / ALL_J FASTAs ──────────────────────────────
    # The internal_data V database (used by igblastn for chain-type detection)
    # and MakeDb.py both need an all-locus V/J FASTA. Without ALL_V, the bash
    # Step 3 falls back to IGHV-only, so kappa/lambda queries get aligned only
    # to heavy V genes internally and are misclassified as VH (no light J call).
    gapped_dir2 <- file.path(opt$outdir, "germlines", "gapped")
    combine_loci <- function(loci, out_name) {
      all_seqs <- NULL
      for (lc in loci) {
        fa <- file.path(gapped_dir2, paste0(imgt_file_prefix, "_", lc, ".fasta"))
        if (file.exists(fa)) {
          s <- Biostrings::readDNAStringSet(fa)
          all_seqs <- if (is.null(all_seqs)) s else c(all_seqs, s)
        }
      }
      if (!is.null(all_seqs) && length(all_seqs) > 0L) {
        out_fa <- file.path(gapped_dir2, paste0(imgt_file_prefix, "_", out_name, ".fasta"))
        Biostrings::writeXStringSet(all_seqs, out_fa)
        cat(sprintf("  Wrote %s: %d sequences (%s)\n",
                    out_name, length(all_seqs), paste(loci, collapse="+")))
      }
    }
    combine_loci(c("IGHV","IGKV","IGLV"), "ALL_V")
    combine_loci(c("IGHJ","IGKJ","IGLJ"), "ALL_J")
  }

  source_as_is()
  quit(save="no", status=0L)
}
cat(sprintf("Output     : %s\n\n", opt$outdir))

# =============================================================================
# SECTION 1 -- Header normalisation
# =============================================================================
# Regex matching the gene-name body of an IMGT allele string.
.IMGT_RE       <- "(IG[HKL][VDJ]\\d+(?:-\\d+)?(?:-[A-Z0-9]+)?)(\\*\\d+)?"
# OGRDB hash IDs: IGxV<digits>[-<digits>]<UPPERCASE+digit hash>*<digits>
# e.g. IGKV0-2HY3*00  IGKJ0-4JXG*00  (hash is alphanumeric, starts uppercase)
.OGRDB_HASH_RE <- "^IG[HKL][VDJ]\\d+(-\\d+)?[A-Z][A-Z0-9]*\\*\\d+$"

.clean_nm <- function(nm) {
  nm <- trimws(gsub("\\s+", "", nm))
  nm <- gsub("-G-",   "-",  nm)
  nm <- gsub("-G\\*", "\\*", nm)
  nm <- sub("(\\*\\d+)[FP]$", "\\1", nm)
  nm
}

.match_imgt <- function(tok) {
  m <- regmatches(tok, regexpr(.IMGT_RE, tok, perl = TRUE))
  if (length(m) == 1L && nchar(m) > 0L) return(.clean_nm(m))
  ""
}

#' Extract canonical IMGT gene name from one raw FASTA header.
#' Returns c(canonical=..., flag=...)
#' Flags: standard | pipe_imgt | pipe_field | suffix_stripped | ogrdb_hash | opaque
.extract_gene_name <- function(raw, locus = "") {
  tok <- trimws(gsub("\\s+", " ", raw))

  # Format C: pipe-delimited IMGT
  if (grepl("|", tok, fixed = TRUE)) {
    fields <- strsplit(tok, "\\|")[[1]]
    if (length(fields) >= 2L) {
      f2 <- trimws(fields[2L])
      nm <- .match_imgt(f2)
      if (nchar(nm) > 0L) return(c(canonical = nm, flag = "pipe_imgt"))
      if (grepl("^IG[HKL][VDJ]", f2))
        return(c(canonical = .clean_nm(f2), flag = "pipe_imgt"))
    }
    for (f in strsplit(tok, "\\|")[[1]]) {
      nm <- .match_imgt(trimws(f))
      if (nchar(nm) > 0L) return(c(canonical = nm, flag = "pipe_field"))
    }
  }

  tok1 <- strsplit(tok, "\\s")[[1]][1L]

  # Format D: OGRDB alphanumeric hash ID -- detect BEFORE standard regex
  # The hash (e.g. HY3) would otherwise be consumed as part of the gene suffix,
  # collapsing all hash variants of the same gene to the same name.
  if (grepl(.OGRDB_HASH_RE, tok1, perl = TRUE)) {
    message(sprintf("  [HEADER] OGRDB hash ID: '%s' (locus=%s) kept verbatim", tok1, locus))
    return(c(canonical = tok1, flag = "ogrdb_hash"))
  }

  # Format A: direct IMGT match
  nm <- .match_imgt(tok1)
  if (nchar(nm) > 0L) return(c(canonical = nm, flag = "standard"))

  # Format B: trailing underscore junk (e.g. IGHV1-14*01_S3227)
  parts <- strsplit(tok1, "_")[[1]]
  if (length(parts) >= 2L) {
    for (n in seq(length(parts) - 1L, 1L)) {
      nm <- .match_imgt(paste(parts[seq_len(n)], collapse = "_"))
      if (nchar(nm) > 0L) return(c(canonical = nm, flag = "suffix_stripped"))
    }
  }

  # Opaque fallback
  if (!grepl("^IG[HKL][VDJ]", tok1))
    warning(sprintf("  [HEADER] Cannot parse '%s' (locus=%s); kept verbatim",
                    substr(raw, 1L, 80L), locus))
  c(canonical = .clean_nm(tok1), flag = "opaque")
}

#' Normalise all headers in a DNAStringSet.
#' Returns list(seqs=DNAStringSet, map=data.table(raw_header, normalised_name, parse_flag))
normalise_headers <- function(seqs, locus = "") {
  if (is.null(seqs) || length(seqs) == 0L)
    return(list(seqs = seqs, map = data.table()))

  raw_names <- names(seqs)
  result    <- vapply(raw_names, .extract_gene_name,
                      FUN.VALUE = character(2L), locus = locus)
  canonical <- result["canonical", ]
  flags     <- result["flag",      ]

  # Disambiguate duplicate normalised names within this file.
  # Track occurrence count per BASE name (ignoring _dup suffixes already assigned)
  # so each duplicate gets a unique sequential suffix: _dup2, _dup3, _dup4, ...
  name_count <- list()   # base_name -> number of times seen so far
  for (i in seq_along(canonical)) {
    nm  <- canonical[i]   # this is the base name (not yet modified at position i)
    cnt <- if (is.null(name_count[[nm]])) 0L else name_count[[nm]]
    name_count[[nm]] <- cnt + 1L
    if (cnt > 0L) {
      # cnt is the number of PRIOR occurrences; this is occurrence (cnt+1)
      new_nm       <- paste0(nm, "_dup", cnt + 1L)
      message(sprintf("  [HEADER] Duplicate '%s' at pos %d -> '%s'", nm, i, new_nm))
      canonical[i] <- new_nm
      flags[i]     <- "duplicate_resolved"
    }
  }
  names(seqs) <- canonical
  list(seqs = seqs,
       map  = data.table(raw_header = raw_names, normalised_name = canonical,
                         parse_flag = flags))
}

# =============================================================================
# SECTION 2 -- Core helpers
# =============================================================================

#' Read a FASTA, normalise headers, print a parse-flag summary.
#' Returns list(seqs, map) or NULL.
read_fasta_safe <- function(path, locus = "", source_tag = "custom") {
  if (is.null(path) || !file.exists(path)) {
    message(sprintf("  [WARN] %s (%s) not found: %s", locus, source_tag,
                    ifelse(is.null(path), "(NULL)", path)))
    return(NULL)
  }
  seqs <- tryCatch(readDNAStringSet(path), error = function(e) {
    message(sprintf("  [ERROR] Cannot read %s: %s", path, conditionMessage(e)))
    NULL
  })
  if (is.null(seqs) || length(seqs) == 0L) {
    message(sprintf("  [WARN] %s (%s): empty: %s", locus, source_tag, path))
    return(NULL)
  }
  # Remove exact-duplicate sequences (same ungapped content) before naming.
  # Keep first occurrence; later duplicates are silently dropped here;
  # the header map still records them via the normalise step below.
  raw_content <- as.character(DECIPHER::RemoveGaps(seqs, removeGaps = "all"))
  dup_seq     <- duplicated(raw_content)
  if (any(dup_seq)) {
    message(sprintf("  [DEDUP] %s (%s): removing %d exact-duplicate sequences",
                    locus, source_tag, sum(dup_seq)))
    seqs <- seqs[!dup_seq]
  }
  norm <- normalise_headers(seqs, locus = locus)
  cat(sprintf("  Loaded %-6s (%s): %d seqs\n", locus, source_tag, length(norm$seqs)))
  fs <- norm$map[, .N, by = parse_flag][order(-N)]
  for (i in seq_len(nrow(fs)))
    cat(sprintf("    %-25s : %d\n", fs$parse_flag[i], fs$N[i]))
  norm
}

ungap <- function(seqs) {
  if (is.null(seqs)) return(NULL)
  # removeGaps requires XStringSet; for character vectors strip gap chars directly.
  if (is(seqs, "XStringSet"))
    return(DECIPHER::RemoveGaps(seqs, removeGaps = "all"))
  if (is.character(seqs))
    return(gsub("[.-]", "", seqs))   # [.-] in a char class: literal dot and hyphen
  DECIPHER::RemoveGaps(seqs, removeGaps = "all")
}

normalise_name <- function(nm) {
  nm <- gsub("-G-",   "-",  nm)
  nm <- gsub("-G\\*", "\\*", nm)
  nm <- sub("(\\*\\d+)[FP]$", "\\1", nm)
  nm
}

# ---------------------------------------------------------------------------
# TRIG-NC (TR/IG Nomenclature Review Committee) colon-format naming
# ---------------------------------------------------------------------------

#' Uniformly re-derive TRIG-NC names across a MERGED (custom + reference) set.
#'
#' Size-ranks ALL genes together within each family so reference IMGT genes and
#' novel custom genes share one consistent numbering. Policy:
#'   - FAMILY: lifted from the IMGT name where the sequence already has a
#'     standard IMGT family (IGHV1-.. -> family 1), else assigned above the
#'     reference max family (novel families size-ranked, most genes = lowest).
#'   - GENE: within each family, cluster sequences into genes by the IMGT gene
#'     base for IMGT-named refs, and by identity for novel/OGRDB names; then
#'     size-rank those gene clusters (most alleles = gene 1).
#'   - ALLELE: sequential within each gene.
#'
#' @param merged_seqs DNAStringSet of the merged custom+reference sequences
#' @param ref_seqs    DNAStringSet of the reference (to compute ref max family)
#' @param label       locus label for logging
#' @return DNAStringSet with TRIG-NC colon names
.trignc_rename_merged <- function(merged_seqs, ref_seqs, label = "", annot = NULL) {
  if (is.null(merged_seqs) || length(merged_seqs) == 0L) return(merged_seqs)
  nms <- names(merged_seqs)
  n   <- length(nms)

  seg_prefix <- {
    fp <- regmatches(nms[1], regexpr("^(IG[HKL][VDJ]|TR[ABGD][VDJ])", nms[1]))
    if (length(fp) && nzchar(fp)) fp else sub("^([A-Z]+[VDJ]).*$", "\\1", nms[1])
  }

  extract_fam <- function(v) suppressWarnings(as.integer(
    sub("^(?:IG[HKL][VDJ]|TR[ABGD][VDJ])([0-9]+).*$", "\\1", v)))

  # The merged FASTA may carry PIgLET cluster names (IGHVF15-G111*02) rather
  # than IMGT names. PIgLET F-numbers are internal cluster ids, NOT IMGT
  # families. Resolve each name to its underlying IMGT allele via the
  # annotation table (new_allele -> imgt_allele) before extracting the family.
  imgt_lookup <- character(0L)
  fam_cluster_to_imgtfam <- integer(0L)   # PIgLET family_cluster -> IMGT family
  name_to_famcluster     <- character(0L) # sequence name -> family_cluster
  if (!is.null(annot) && nrow(annot) > 0L &&
      all(c("imgt_allele","new_allele") %in% names(annot))) {
    ad0 <- as.data.frame(annot, stringsAsFactors = FALSE)
    # Map both the PIgLET new_allele AND the imgt_allele to the imgt_allele,
    # so lookup works whichever name variant is on the sequence.
    imgt_lookup <- setNames(c(ad0$imgt_allele, ad0$imgt_allele),
                            c(ad0$new_allele,  ad0$imgt_allele))

    # Build family_cluster -> dominant IMGT family from REFERENCE members.
    # OGRDB/novel alleles that share a PIgLET family_cluster with reference
    # alleles inherit that cluster's dominant IMGT family (rank-sorted liftover),
    # instead of being treated as brand-new families.
    if ("family_cluster" %in% names(ad0)) {
      imgt_fam_of <- suppressWarnings(as.integer(
        sub("^(?:IG[HKL][VDJ]|TR[ABGD][VDJ])([0-9]+)-.*$", "\\1", ad0$imgt_allele)))
      # only reference rows with a real IMGT family (has a dash) count as evidence
      is_ref_imgt <- (ad0$source == "reference") & !is.na(imgt_fam_of) &
                     grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ])[0-9]+-", ad0$imgt_allele)
      if (any(is_ref_imgt)) {
        fc  <- ad0$family_cluster[is_ref_imgt]
        fam <- imgt_fam_of[is_ref_imgt]
        # dominant IMGT family per cluster = most frequent
        split_fam <- split(fam, fc)
        fam_cluster_to_imgtfam <- vapply(split_fam, function(v) {
          tt <- sort(table(v), decreasing = TRUE); as.integer(names(tt)[1L])
        }, integer(1L))
      }
      # name -> family_cluster for every allele (both name variants)
      name_to_famcluster <- setNames(
        c(as.character(ad0$family_cluster), as.character(ad0$family_cluster)),
        c(ad0$new_allele, ad0$imgt_allele))
    }
  }
  # For each FASTA name, the IMGT name to lift the family from:
  resolve_imgt <- function(x) {
    r <- imgt_lookup[x]
    ifelse(is.na(r), x, r)
  }
  imgt_names <- resolve_imgt(nms)

  # Reference max family (from the reference set, all standard IMGT families)
  ref_max_fam <- {
    rf <- extract_fam(resolve_imgt(names(ref_seqs)))
    rf <- rf[!is.na(rf)]
    if (length(rf)) max(rf) else 0L
  }

  # ── Family assignment ────────────────────────────────────────────────────
  # A name has a liftable IMGT family if it starts with locus + digits, whether
  # or not it has a dash: both "IGHV1-11*01" and bare "IGHV1" yield family 1.
  # OGRDB hash names (IGKV0-4HMC*00) have family 0 which is a placeholder, so
  # treat family 0 as NON-liftable (novel).
  # Use the RESOLVED IMGT names for family detection, not the PIgLET names.
  has_family <- grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ])[0-9]+", imgt_names)
  fam_raw    <- extract_fam(imgt_names)
  is_std_imgt <- has_family & !is.na(fam_raw) & fam_raw > 0L
  fam_num <- ifelse(is_std_imgt, fam_raw, NA_integer_)

  # Rank-sorted liftover for OGRDB/novel names that lack a direct IMGT family
  # (family 0 placeholder, PIgLET F-names): if the sequence shares a PIgLET
  # family_cluster with reference IMGT alleles, inherit that cluster's dominant
  # IMGT family. This collapses the many OGRDB IGKV0-* alleles onto real
  # families (e.g. cluster 2 -> IGKV4) instead of spawning families 60+.
  if (length(fam_cluster_to_imgtfam) > 0L && length(name_to_famcluster) > 0L) {
    still_na <- which(is.na(fam_num))
    for (i in still_na) {
      fc <- name_to_famcluster[nms[i]]
      if (is.na(fc)) fc <- name_to_famcluster[imgt_names[i]]
      if (!is.na(fc) && !is.na(fam_cluster_to_imgtfam[fc])) {
        fam_num[i] <- fam_cluster_to_imgtfam[fc]
      }
    }
    n_lifted <- sum(is.na(fam_raw) & !is.na(fam_num))
    if (n_lifted > 0L)
      cat(sprintf("  [TRIG-NC] %s: %d OGRDB/novel allele(s) lifted to reference family via family_cluster\n",
                  label, n_lifted))
  }

  # PIgLET cluster lookups (by sequence name) for novel grouping.
  #   family_cluster : 75%-threshold family grouping (novel FAMILY assignment)
  #   cluster_id     : 95%-threshold allele grouping (novel GENE assignment)
  fam_cl_lookup  <- character(0L)
  gene_cl_lookup <- character(0L)
  if (!is.null(annot) && nrow(annot) > 0L &&
      all(c("imgt_allele","family_cluster","cluster_id") %in% names(annot))) {
    ad <- as.data.frame(annot, stringsAsFactors = FALSE)
    fam_cl_lookup  <- setNames(as.character(ad$family_cluster), ad$imgt_allele)
    gene_cl_lookup <- setNames(as.character(ad$cluster_id),     ad$imgt_allele)
  }

  # Gene grouping key:
  #   - IMGT-named with a dash (IGHV1-11): use the gene base "IGHV1-11"
  #   - novel/bare: use PIgLET cluster_id (95% allele cluster) when available,
  #     else fall back to sequence identity.
  seqs_char <- as.character(DECIPHER::RemoveGaps(merged_seqs, removeGaps = "all"))
  # Gene grouping keys off the RESOLVED IMGT name. Full IMGT gene names group by
  # gene base; bare-family / novel names group by PIgLET allele cluster (95%),
  # falling back to sequence identity.
  has_dash  <- grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ])[0-9]+-[0-9]", imgt_names)
  gene_key <- character(n)
  for (i in seq_len(n)) {
    if (has_dash[i]) {
      gene_key[i] <- sub("[*].*$", "", imgt_names[i])   # IMGT gene base
    } else {
      gc <- gene_cl_lookup[nms[i]]
      if (is.na(gc)) gc <- gene_cl_lookup[imgt_names[i]]
      gene_key[i] <- if (!is.na(gc) && nzchar(gc)) paste0("gclust:", gc)
                     else paste0("seqgene:", seqs_char[i])
    }
  }

  # Assign novel families to rows with no liftable IMGT family (fam_num is NA):
  # OGRDB hash (family 0), PIgLET F-cluster names, etc. Group these into
  # families by PIgLET's family_cluster (the 75% family threshold) so related
  # novel sequences SHARE a family, rather than each becoming its own family.
  need_fb <- is.na(fam_num)
  if (any(need_fb)) {
    next_fam <- max(c(fam_num[!need_fb], ref_max_fam), 0L, na.rm = TRUE) + 1L
    # Novel family key: PIgLET family_cluster if available, else the gene_key.
    fb_famkey <- vapply(which(need_fb), function(i) {
      fc <- fam_cl_lookup[nms[i]]
      if (!is.na(fc) && nzchar(fc)) paste0("fclust:", fc) else gene_key[i]
    }, character(1L))
    # Size-rank novel families by number of distinct genes (most genes = lowest
    # new family number).
    fb_gene <- gene_key[need_fb]
    fam_gene_counts <- tapply(fb_gene, fb_famkey, function(g) length(unique(g)))
    fb_order  <- names(sort(fam_gene_counts, decreasing = TRUE))
    fb_fam_rank <- setNames(next_fam + seq_along(fb_order) - 1L, fb_order)
    fam_num[need_fb] <- fb_fam_rank[fb_famkey]
    cat(sprintf("  [TRIG-NC] %s: %d truly-novel allele(s) in %d novel famil%s -> families %d..%d (ref max %d)\n",
                label, sum(need_fb), length(fb_order),
                if (length(fb_order)==1L) "y" else "ies",
                next_fam, next_fam + length(fb_order) - 1L, ref_max_fam))
  }

  # ── Within each family: size-rank genes, number alleles ──────────────────
  out <- character(n)
  for (f in sort(unique(fam_num))) {
    rows <- which(fam_num == f)
    keys <- gene_key[rows]
    gcount <- table(keys)
    gorder <- names(sort(gcount, decreasing = TRUE))   # most alleles first
    grank  <- setNames(seq_along(gorder), gorder)
    actr   <- setNames(integer(length(gorder)), gorder)
    for (ri in rows) {
      k <- gene_key[ri]
      g <- grank[[k]]
      actr[k] <- actr[k] + 1L
      out[ri] <- sprintf("%s:%02d:%03d:%03d", seg_prefix, f, g, actr[k])
    }
  }
  # Record the IMGT/original -> TRIG-NC mapping for the provenance table.
  # `nms` are the names BEFORE conversion (post-merge IMGT/OGRDB names);
  # `out` are the TRIG-NC colon names. Store globally so the final name map
  # can join original input ids through to the colon names.
  if (exists(".trignc_map_acc", envir = .GlobalEnv)) {
    prev <- get(".trignc_map_acc", envir = .GlobalEnv)
  } else {
    prev <- data.table(before = character(0L), after = character(0L), locus = character(0L))
  }
  # Record BOTH the on-FASTA name (nms) and the resolved IMGT name (imgt_names)
  # -> colon (out), so downstream joins can match on either variant.
  assign(".trignc_map_acc",
         rbind(prev,
               data.table(before = nms,        after = out, locus = label),
               data.table(before = imgt_names,  after = out, locus = label)),
         envir = .GlobalEnv)

  names(merged_seqs) <- out
  merged_seqs
}

#' Convert a DNAStringSet to a NAMED character vector for PIgLET.
#' This is the critical fix: as.character() alone drops names.
dss_to_named_vec <- function(seqs) setNames(as.character(seqs), names(seqs))

# =============================================================================
# SECTION 3 -- PIgLET wrapper with pre-flight validation
# =============================================================================

#' Validate and sanitise a named character vector before passing to PIgLET.
#'
#' PIgLET crashes (subscript out of bounds in germ.dist) when:
#'   1. germline_set is unnamed  -> names must be explicitly set
#'   2. Duplicate names exist    -> germ.dist uses names as row/col keys
#'   3. Names contain characters PIgLET's regex doesn't handle (rare)
#'   4. Sequences are shorter than trim_3prime_side after ungapping
#'      (generates an empty/degenerate distance matrix)
#'
#' Returns a sanitised named character vector and a data.table recording
#' any sequences that were dropped, with the reason.
preflight_piglet <- function(named_vec, trim3, label = "") {
  n_in   <- length(named_vec)
  dropped <- data.table(name = character(), reason = character())

  # 1. Must be named (belt-and-suspenders check)
  if (is.null(names(named_vec)) || any(nchar(names(named_vec)) == 0L)) {
    stop(sprintf("[%s] preflight_piglet: unnamed or empty-named sequences", label))
  }

  # 2. Remove duplicate names (keep first occurrence; the joint set already
  #    puts reference first so reference alleles are retained on collision)
  dup_names <- duplicated(names(named_vec))
  if (any(dup_names)) {
    dups <- names(named_vec)[dup_names]
    message(sprintf("  [PIgLET pre-flight] %s: removing %d duplicate names: %s",
                    label, sum(dup_names),
                    paste(head(unique(dups), 5), collapse = ", ")))
    dropped <- rbind(dropped,
                     data.table(name = dups, reason = "duplicate_name"))
    named_vec <- named_vec[!dup_names]
  }

  # 3. Filter sequences by length.
  #
  #    PIgLET receives GAPPED sequences and trims internally by gapped position
  #    (substr(seq, 1, trim_3prime_side)).  Standard IMGT V segments are ~312
  #    gapped nt, well above the default trim3=318 -- PIgLET takes the full
  #    sequence when it is shorter than trim3, so these do NOT need to be removed.
  #
  #    PREVIOUS BUG: used nchar(gsub("[\\.\\-]", "", seq)) to count "ungapped"
  #    length.  This stripped IMGT gap dots AND frame-shift dashes, making every
  #    ~312-position gapped V appear as ~240 nt and causing mass false removal.
  #
  #    CORRECT behaviour: only drop truly empty sequences (len == 0).
  #    Warn about sequences shorter than trim3 in gapped length, but keep them.
  if (!is.null(trim3) && trim3 > 0L) {
    gapped_len <- nchar(named_vec)          # full string length, dots included
    too_short  <- gapped_len == 0L          # only truly empty
    if (any(too_short)) {
      short_names <- names(named_vec)[too_short]
      message(sprintf(
        "  [PIgLET pre-flight] %s: removing %d empty sequences: %s",
        label, sum(too_short),
        paste(head(short_names, 5), collapse = ", ")))
      dropped   <- rbind(dropped,
                         data.table(name = short_names, reason = "empty_sequence"))
      named_vec <- named_vec[!too_short]
    }
    # Warn (keep) sequences genuinely shorter than trim3 in gapped length
    warn_short <- gapped_len > 0L & gapped_len < trim3
    if (any(warn_short))
      message(sprintf(
        "  [PIgLET pre-flight] %s: %d seqs have gapped len < trim3=%d (kept; PIgLET uses full length): %s",
        label, sum(warn_short), trim3,
        paste(head(names(named_vec)[warn_short], 5), collapse = ", ")))
  }

  # 4. Need at least 2 sequences for clustering
  if (length(named_vec) < 2L)
    stop(sprintf("[%s] Only %d sequence(s) remain after pre-flight; cannot cluster",
                 label, length(named_vec)))

  if (length(named_vec) < n_in)
    message(sprintf("  [PIgLET pre-flight] %s: %d -> %d sequences after filtering",
                    label, n_in, length(named_vec)))

  list(vec = named_vec, dropped = dropped)
}

#' Run PIgLET::inferAlleleClusters safely.
#'
#' Following the pattern from the original working code:
#'   inferAlleleClusters(germline_set = as.character(mrl_l), ...)
#' but crucially using dss_to_named_vec() instead of plain as.character().
#'
#' Returns list(tbl = alleleClusterTable data.table,
#'              dropped = pre-flight dropped sequences data.table)
run_piglet <- function(seqs_gapped, trim3 = 318L, mask5 = 0L,
                       fam_thresh = 75, allele_thresh = 95, label = "") {
  # Convert to named character vector (THE key fix)
  named_vec <- dss_to_named_vec(seqs_gapped)

  pf <- preflight_piglet(named_vec, trim3 = trim3, label = label)
  named_vec <- pf$vec

  cat(sprintf("  Running PIgLET on %d sequences (%s)...\n",
              length(named_vec), label))

  res <- tryCatch(
    inferAlleleClusters(
      germline_set             = named_vec,   # named character vector
      trim_3prime_side         = trim3,
      mask_5prime_side         = mask5,
      family_threshold         = fam_thresh,
      allele_cluster_threshold = allele_thresh
    ),
    error = function(e) {
      message(sprintf("\n  [ERROR] PIgLET failed for %s: %s", label, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(res)) {
    # Return a minimal table so downstream code doesn't crash
    tbl <- data.table(imgt_allele = names(named_vec),
                      new_allele  = names(named_vec),
                      cluster_id  = seq_along(named_vec),
                      family_cluster = seq_along(named_vec))
    message(sprintf("  [WARN] %s: PIgLET failed; sequences keep their current names", label))
    return(list(tbl = tbl, dropped = pf$dropped, dist_mat = NULL))
  }

  # Access slots via $ (confirmed working for this PIgLET version)
  act      <- res$alleleClusterTable
  tbl      <- as.data.table(act)

  # Normalise column names across PIgLET versions. The allele-name column has
  # been variously called: imgt_allele, allele, sample_allele, germline_call,
  # or the sequence name column. The new-cluster-name column: new_allele,
  # new_call, or func_group. Detect and rename to canonical names.
  nms_lower <- tolower(names(tbl))
  # --- source (input) allele name column -> imgt_allele ---
  if (!"imgt_allele" %in% names(tbl)) {
    cand <- c("imgt_allele","allele","sample_allele","germline_call",
              "gene","seq_name","sequence_id","name")
    hit  <- cand[cand %in% names(tbl)]
    if (length(hit) == 0L) hit <- names(tbl)[nms_lower %in% cand]
    if (length(hit) > 0L) {
      setnames(tbl, hit[1L], "imgt_allele")
    } else {
      # Last resort: use the row names or first character column
      char_cols <- names(tbl)[vapply(tbl, is.character, logical(1L))]
      if (length(char_cols) > 0L) setnames(tbl, char_cols[1L], "imgt_allele")
      else stop(sprintf("run_piglet(%s): cannot find allele-name column in PIgLET output; columns are: %s",
                        label, paste(names(tbl), collapse=", ")))
    }
  }
  # --- new cluster name column -> new_allele ---
  if (!"new_allele" %in% names(tbl)) {
    cand2 <- c("new_allele","new_call","func_group","new_tag","threshold")
    hit2  <- cand2[cand2 %in% names(tbl)]
    if (length(hit2) == 0L) hit2 <- names(tbl)[tolower(names(tbl)) %in% cand2]
    if (length(hit2) > 0L) setnames(tbl, hit2[1L], "new_allele")
    else tbl[, new_allele := imgt_allele]  # fall back to identity naming
  }
  # --- ensure cluster_id / family_cluster exist (used by TRIG-NC sizing) ---
  if (!"cluster_id" %in% names(tbl)) {
    cand3 <- c("cluster_id","allele_cluster","cluster","func_group_id")
    hit3  <- cand3[cand3 %in% names(tbl)]
    if (length(hit3) > 0L) setnames(tbl, hit3[1L], "cluster_id")
    else tbl[, cluster_id := .I]
  }
  if (!"family_cluster" %in% names(tbl)) {
    cand4 <- c("family_cluster","family","fam_cluster","subgroup")
    hit4  <- cand4[cand4 %in% names(tbl)]
    if (length(hit4) > 0L) setnames(tbl, hit4[1L], "family_cluster")
    else tbl[, family_cluster := 1L]
  }
  if (!"removed_duplicated" %in% names(tbl)) {
    cand5 <- c("removed_duplicated","is_duplicate","duplicated")
    hit5  <- cand5[cand5 %in% names(tbl)]
    if (length(hit5) > 0L) setnames(tbl, hit5[1L], "removed_duplicated")
    else tbl[, removed_duplicated := FALSE]
  }

  cat(sprintf("  [PIgLET] %s: table columns = %s\n",
              label, paste(names(tbl), collapse=", ")))
  tbl[, new_allele := normalise_name(new_allele)]

  # Distance matrix for nearest-reference lookup in step 5
  dist_mat <- tryCatch(as.matrix(res$distanceMatrix), error = function(e) NULL)
  list(tbl = tbl, dropped = pf$dropped, dist_mat = dist_mat)
}

#' Apply the PIgLET alleleClusterTable to rename a DNAStringSet.
#' tbl may be a data.frame or data.table; uses base-R subsetting throughout.
rename_by_table <- function(seqs, tbl, label = "") {
  # Coerce to data.frame so [ works identically regardless of input class
  tbl_df  <- as.data.frame(tbl, stringsAsFactors = FALSE)
  renamed <- seqs
  for (i in seq_along(names(seqs))) {
    orig <- names(seqs)[i]
    rows <- tbl_df[tbl_df$imgt_allele == orig, , drop = FALSE]
    if (nrow(rows) == 0L) {
      message(sprintf("  [WARN] %s: no PIgLET mapping for '%s'; keeping name", label, orig))
    } else {
      names(renamed)[i] <- rows$new_allele[1L]
    }
  }
  renamed
}

#' Merge custom + reference.  Custom takes priority on name AND content collision.
#' Sequences with the same name: custom wins.
#' Sequences with the same content but different names: custom name wins.
#' Reference-only sequences (unique name and content): appended.
merge_with_priority <- function(custom_seqs, ref_seqs) {
  if (is.null(custom_seqs)) return(ref_seqs)
  if (is.null(ref_seqs))    return(custom_seqs)
  # Name dedup: keep custom on collision
  ref_name_only <- ref_seqs[!names(ref_seqs) %in% names(custom_seqs)]
  # Content dedup: exclude ref seqs whose ungapped content already appears in custom
  custom_content <- as.character(DECIPHER::RemoveGaps(custom_seqs, removeGaps = "all"))
  ref_content    <- as.character(DECIPHER::RemoveGaps(ref_name_only, removeGaps = "all"))
  ref_novel      <- ref_name_only[!ref_content %in% custom_content]
  merged <- c(custom_seqs, ref_novel)
  merged[!duplicated(names(merged))]  # final name safety check
}

#' Return an empty annotation data.table with all expected columns.
#' Used as a safe default when a locus has no annotation to report.
.empty_annot <- function() {
  data.table(imgt_allele    = character(0L),
             new_allele     = character(0L),
             cluster_id     = character(0L),
             family_cluster = character(0L),
             source         = character(0L),
             locus          = character(0L))
}

# =============================================================================
# SECTION 4 -- Joint PIgLET clustering with reference-family inheritance
# =============================================================================
#
# Strategy (mirrors the original working code but adds joint clustering):
#
#   Original code did:
#     asc <- inferAlleleClusters(germline_set = as.character(mrl_l), ...)
#     for(seqname in seq_along(names(mrl_l))) {
#       names(mrl_l_asc)[seqname] <- annot[annot$imgt_allele==seqname,]$new_allele
#     }
#
#   We do the same but:
#     a) cluster custom + reference JOINTLY so novel alleles inherit family IDs
#     b) use dss_to_named_vec() to preserve names
#     c) post-process: for each custom allele that co-clusters with a reference
#        allele, replace the PIgLET-assigned family with the reference family
# =============================================================================

annotate_custom_with_ref <- function(custom_gapped, ref_gapped,
                                     trim3, fam_thresh, allele_thresh,
                                     label = "") {
  if (is.null(custom_gapped))
    return(list(renamed_custom = NULL, annot_table = .empty_annot(),
                dropped = data.table()))

  ref_names    <- names(ref_gapped)
  custom_names <- names(custom_gapped)

  # --- 0. Early exit: are all custom sequences already in the reference? ----
  #   Compare by sequence content (ungapped) so gap-only differences don't
  #   trigger unnecessary clustering.  If every custom sequence is content-
  #   identical to a reference sequence (by name OR by content), skip PIgLET.
  custom_seqs_chr <- as.character(DECIPHER::RemoveGaps(custom_gapped, removeGaps = "all"))
  ref_seqs_chr    <- as.character(DECIPHER::RemoveGaps(ref_gapped,    removeGaps = "all"))

  # Compare by UNGAPPED SEQUENCE CONTENT ONLY.
  # Do NOT use name matching: OGRDB hash IDs (e.g. IGLV0-CXWW*00) normalise
  # to the same canonical name as their IMGT equivalent (IGLV1*01) even when
  # they represent distinct sequences, causing false "already in reference" hits.
  novel_mask <- !(custom_seqs_chr %in% ref_seqs_chr)
  n_novel <- sum(novel_mask)

  if (n_novel == 0L) {
    cat(sprintf("  %s: all %d custom seqs in reference; skipping PIgLET\n", label, length(custom_gapped)))
    # Build a minimal annotation table mapping each custom name to its ref name
    annot_rows <- lapply(seq_along(custom_gapped), function(i) {
      cn  <- names(custom_gapped)[i]
      seq <- custom_seqs_chr[i]
      # Find matching ref entry: prefer name match, then content match
      if (cn %in% ref_names) {
        rn <- cn
      } else {
        rn <- ref_names[ref_seqs_chr == seq][1L]
        if (is.na(rn)) rn <- cn
      }
      data.table(imgt_allele = cn, new_allele = rn,
                 cluster_id = NA_character_, family_cluster = NA_character_,
                 source = "custom")
    })
    annot_tbl <- rbindlist(annot_rows)
    # Keep custom sequences with their (possibly corrected) names
    renamed_custom <- custom_gapped
    for (i in seq_along(names(renamed_custom)))
      names(renamed_custom)[i] <- annot_tbl$new_allele[i]
    return(list(renamed_custom = renamed_custom,
                annot_table    = annot_tbl,
                dropped        = data.table()))
  }

  cat(sprintf("  %s: %d novel custom sequences; proceeding with joint PIgLET\n", label, n_novel))

  # --- 1. Joint set: reference first (reference names are stable anchors) ---
  # Pre-clean reference names: strip _dup suffixes and ensure *NN allele numbers
  # so distance matrix row/col names are proper IMGT gene names.
  ref_for_joint <- ref_gapped
  ref_nms <- names(ref_for_joint)
  # Strip _dup<N> suffixes from reference names
  ref_nms_clean <- sub("_dup\\d+$", "", ref_nms)
  # Add *01 to any name missing an allele number
  no_star <- !grepl("[*]", ref_nms_clean)
  if (any(no_star))
    ref_nms_clean[no_star] <- paste0(ref_nms_clean[no_star], "*01")
  # Re-number per gene base to ensure uniqueness
  gene_seen <- list()
  for (k in seq_along(ref_nms_clean)) {
    gb  <- sub("[*].*$", "", ref_nms_clean[k])
    cur <- suppressWarnings(as.integer(sub(".*[*]", "", ref_nms_clean[k])))
    if (is.na(cur)) cur <- 1L
    prev <- gene_seen[[gb]]
    if (is.null(prev)) {
      gene_seen[[gb]] <- cur
    } else {
      nxt <- as.integer(prev) + 1L
      ref_nms_clean[k] <- paste0(gb, "*", sprintf("%02d", nxt))
      gene_seen[[gb]]   <- nxt
    }
  }
  names(ref_for_joint) <- ref_nms_clean

  joint <- c(ref_for_joint, custom_gapped)

  # Deduplicate by SEQUENCE CONTENT -- keep reference copy on tie
  content_dup <- duplicated(as.character(joint))
  joint_dd    <- joint[!content_dup]

  # Also remove any NAME duplicates that survive content-dedup
  name_dup  <- duplicated(names(joint_dd))
  joint_dd  <- joint_dd[!name_dup]

  # Update ref_names to match the cleaned names for distance matrix lookup
  # (a named vector mapping original name -> clean name for later reconciliation)
  ref_name_map <- setNames(ref_nms_clean, ref_nms)
  ref_names_clean <- ref_nms_clean
  # Valid reference names: exclude any _dup sequences from being used as a
  # reference anchor. _dup names are internal artefacts from input deduplication
  # and should never propagate into gene-base assignments.
  ref_names_valid <- ref_names_clean[!grepl("_dup", ref_names_clean, fixed=TRUE)]
  # Reverse map: cleaned name -> original IMGT name (preserves full subtype)
  # e.g. "IGHV1-2*01" -> "IGHV1-2*01" (unchanged if clean)
  # Exclude entries where the original name is a _dup artefact — in those cases
  # the cleaned name is preferable (it has already stripped the _dup suffix).
  valid_orig <- ifelse(grepl("_dup", ref_nms, fixed=TRUE), ref_nms_clean, ref_nms)
  rev_ref_name_map <- setNames(valid_orig, ref_nms_clean)

  # --- 2. Run PIgLET on the joint deduplicated set --------------------------
  pf_res <- run_piglet(joint_dd, trim3 = trim3, fam_thresh = fam_thresh,
                       allele_thresh = allele_thresh, label = label)
  tbl    <- pf_res$tbl
  dist_mat <- if (!is.null(pf_res$dist_mat)) pf_res$dist_mat else NULL

  # --- 3. Split annotation into reference rows and custom rows --------------
  # tbl uses names from joint_dd which has cleaned ref names
  tbl_ref    <- tbl[tbl$imgt_allele %in% ref_names_clean, , drop = FALSE]
  tbl_custom <- tbl[tbl$imgt_allele %in% custom_names,    , drop = FALSE]

  # --- 3b. Handle sequences missing from tbl -----------------------------------
  # Two cases:
  #   (i)  Content-identical to a reference allele: removed by our content-dedup
  #         -> adopt the reference IMGT name directly.
  #   (ii) PIgLET internally merged it (removed_duplicated=TRUE in alleleClusterTable)
  #        or dropped it entirely (not in tbl at all): use distance matrix to find
  #        the closest sequence that IS retained, then adopt its cluster assignment.
  missing_custom <- custom_names[!custom_names %in% as.character(tbl$imgt_allele)]

  # PIgLET-removed sequences (removed_duplicated=TRUE):
  # Only re-route EXACT duplicates (distance == 0 to another retained sequence).
  # Near-duplicates (distance > 0) keep their tbl_custom entry and cluster
  # assignment — they are genuinely distinct alleles worth keeping.
  if ("removed_duplicated" %in% names(tbl) && !is.null(dist_mat)) {
    flagged <- as.character(
      tbl[tbl$imgt_allele %in% custom_names & tbl$removed_duplicated == TRUE,
          imgt_allele])
    exact_dups <- character(0L)
    for (fd in flagged) {
      if (fd %in% rownames(dist_mat)) {
        row_fd  <- dist_mat[fd, , drop = TRUE]
        # Retained sequences (not removed_duplicated) with distance == 0
        retained_names <- as.character(
          tbl[tbl$removed_duplicated == FALSE | is.na(tbl$removed_duplicated),
              imgt_allele])
        cands <- names(row_fd)[names(row_fd) %in% retained_names]
        if (any(row_fd[cands] == 0, na.rm = TRUE))
          exact_dups <- c(exact_dups, fd)
      }
    }
    if (length(exact_dups) > 0L) {
      tbl_custom  <- tbl_custom[!tbl_custom$imgt_allele %in% exact_dups, ]
      missing_custom <- union(missing_custom, exact_dups)
    }
  }

  ref_ungapped_chr <- as.character(DECIPHER::RemoveGaps(ref_gapped, removeGaps = "all"))

  for (dc in missing_custom) {
    dc_seq    <- as.character(DECIPHER::RemoveGaps(custom_gapped[dc], removeGaps = "all"))

    # Case (i): content-identical to a reference allele
    ref_match <- ref_names[ref_ungapped_chr == dc_seq]
    if (length(ref_match) > 0L) {
      imgt_ref_name <- ref_match[1L]
      ref_row <- tbl_ref[tbl_ref$imgt_allele == ref_match[1L], , drop = FALSE]
      if (nrow(ref_row) == 0L) {
        # Ref name may have been cleaned; try via ref_name_map
        cleaned <- ref_name_map[ref_match[1L]]
        if (!is.na(cleaned))
          ref_row <- tbl_ref[tbl_ref$imgt_allele == cleaned, , drop = FALSE]
      }
      if (nrow(ref_row) > 0L) {
        extra <- copy(as.data.table(ref_row))
        # Strip any _dup<N> suffix from the adopted reference name — these are
        # internal disambiguation tags for identical bare-family ref seqs and
        # must never appear in output. The TRIG-NC merged rename will re-derive
        # the final name regardless, but keep the adopted name clean here too.
        imgt_ref_clean <- sub("_dup[0-9]+$", "", imgt_ref_name)
        extra[, imgt_allele   := dc]
        extra[, new_allele    := imgt_ref_clean]
        extra[, piglet_cluster := if ("new_allele" %in% names(ref_row))
                                    ref_row$new_allele[1L] else NA_character_]
        tbl_custom <- rbind(tbl_custom, extra, fill = TRUE)
        message(sprintf("  [INFO] %s: '%s' identical to ref '%s' -> adopts IMGT name",
                        label, dc, imgt_ref_clean))
        next
      }
    }

    # Case (ii): PIgLET removed/merged this sequence — use distance matrix
    # to find the closest RETAINED sequence and inherit its cluster assignment.
    resolved <- FALSE
    if (!is.null(dist_mat) && dc %in% rownames(dist_mat)) {
      row      <- dist_mat[dc, , drop = TRUE]
      # Candidates: sequences retained in tbl (not removed)
      retained <- as.character(tbl$imgt_allele)
      cand_cols <- names(row)[names(row) %in% retained]
      if (length(cand_cols) > 0L) {
        closest   <- cand_cols[which.min(row[cand_cols])]
        # Find the cluster entry for the closest sequence
        proxy_row <- tbl[tbl$imgt_allele == closest, , drop = FALSE]
        if (nrow(proxy_row) > 0L) {
          extra <- copy(as.data.table(proxy_row))
          extra[, imgt_allele   := dc]
          extra[, piglet_cluster := proxy_row$new_allele[1L]]
          # new_allele will be determined in step 5 using the cluster/family rep
          tbl_custom <- rbind(tbl_custom, extra, fill = TRUE)
          message(sprintf(
            "  [INFO] %s: '%s' PIgLET-merged -> proxied via '%s' (dist=%.4f)",
            label, dc, closest, row[closest]))
          resolved <- TRUE
        }
      }
    }
    if (!resolved) {
      # Final fallback: adist to all reference sequences
      # Exclude _dup sequences from the adist reference pool
      ref_joint_valid <- ref_for_joint[!grepl("_dup", names(ref_for_joint), fixed=TRUE)]
      ref_ug_clean <- as.character(DECIPHER::RemoveGaps(ref_joint_valid, removeGaps = "all"))
      dists        <- as.integer(adist(dc_seq, ref_ug_clean))
      best_cleaned <- names(ref_joint_valid)[which.min(dists)]
      # Recover original IMGT name (with full subtype) from the reverse map
      best_ref <- if (exists("rev_ref_name_map") && best_cleaned %in% names(rev_ref_name_map))
                    rev_ref_name_map[[best_cleaned]] else best_cleaned
      ref_row   <- tbl_ref[tbl_ref$imgt_allele == best_cleaned, , drop = FALSE]
      if (nrow(ref_row) == 0L) {
        # Try with original name if cleaned lookup failed
        ref_row <- tbl_ref[tbl_ref$imgt_allele == best_ref, , drop = FALSE]
      }
      if (nrow(ref_row) > 0L) {
        extra <- copy(as.data.table(ref_row))
        extra[, imgt_allele   := dc]
        extra[, piglet_cluster := ref_row$new_allele[1L]]
        tbl_custom <- rbind(tbl_custom, extra, fill = TRUE)
        message(sprintf("  [INFO] %s: '%s' -> fallback adist closest ref '%s'",
                        label, dc, best_ref))
      } else {
        message(sprintf("  [WARN] %s: '%s' could not be resolved; excluded from naming",
                        label, dc))
      }
    }
  }

  # --- 4. Build cluster -> best reference representative -------------------
  #   Use the ORIGINAL IMGT allele name (imgt_allele), NOT the PIgLET-renamed
  #   new_allele, as the representative.  This ensures novel allele names are
  #   built from real IMGT gene names (IGKV8-28) not PIgLET cluster names
  #   (IGKVF25-G102).
  #   "Best" = reference allele with lowest original allele number in the cluster.
  if (nrow(tbl_ref) > 0L && "cluster_id" %in% names(tbl_ref)) {
    # Prefer the most specific IMGT name in each cluster (most dashes in gene base).
    # More dashes = more specific subtype: IGHV1-2 (1 dash) > IGHV1 (0 dashes).
    # Within equal specificity, prefer the lowest allele number.
    tbl_ref[, n_dashes := nchar(gsub("[^-]", "", sub("[*].*$", "", imgt_allele)))]
    tbl_ref[, allele_num_sort := suppressWarnings(
      as.integer(sub(".*[*]", "", imgt_allele)))]
    # Exclude _dup-originated alleles from cluster representatives
    tbl_ref_valid <- tbl_ref[!grepl("_dup", imgt_allele, fixed=TRUE)]
    cluster_ref_rep <- if (nrow(tbl_ref_valid) > 0L) {
      tbl_ref_valid[
        order(-n_dashes, allele_num_sort),
        .(ref_rep = imgt_allele[1L]),
        by = cluster_id
      ]
    } else {
      data.table(cluster_id=character(0L), ref_rep=character(0L))
    }
    tbl_ref[, c("n_dashes", "allele_num_sort") := NULL]
  } else {
    cluster_ref_rep <- data.table(cluster_id = character(0L),
                                  ref_rep    = character(0L))
  }

  # --- 5. Assign final IMGT-style names to custom alleles --------------------
  #
  # Rules:
  #   a) Co-clusters with a reference allele (ref_rep set):
  #      Name = <IMGT_ref_gene_base>*<next_after_ref_max>_<species>
  #      e.g. IGKV8-28*02 in ref, next = *03 -> IGKV8-28*03_mouse
  #
  #   b) Custom-only cluster (no ref co-member), but same family_cluster
  #      as a reference allele (fam_rep set):
  #      Use that reference gene as the base -> same scheme as (a).
  #
  #   c) Completely novel (no ref representative at any level):
  #      Extract locus+family number from PIgLET name F<N> group.
  #      Name = <locus><family_num>-N*<next>_<species>  e.g. IGHV1-N*01_mouse
  #
  # Never produce _dup suffixes, bare family names (no *), or IGKVS/IGHVS names.

  tbl_cust <- copy(tbl_custom)
  tbl_cust[, piglet_cluster := new_allele]

  # Attach cluster-level reference representative (from step 4)
  if ("cluster_id" %in% names(tbl_cust) && nrow(cluster_ref_rep) > 0L) {
    tbl_cust <- merge(tbl_cust, cluster_ref_rep, by = "cluster_id", all.x = TRUE)
  } else {
    tbl_cust[, ref_rep := NA_character_]
  }

  # Attach family-level reference representative (fallback for custom-only clusters)
  if ("family_cluster" %in% names(tbl_cust) && nrow(tbl_ref) > 0L &&
      "family_cluster" %in% names(tbl_ref)) {
    fam_ref_rep <- as.data.table(tbl_ref)[
      order(-(nchar(gsub("[^-]","",sub("[*].*$","",imgt_allele)))),
            suppressWarnings(as.integer(sub(".*[*]", "", imgt_allele)))),
      .(fam_rep = imgt_allele[1L]),
      by = family_cluster]
    tbl_cust <- merge(tbl_cust, fam_ref_rep, by = "family_cluster", all.x = TRUE)
  } else {
    tbl_cust[, fam_rep := NA_character_]
  }

  # Build lookup: IMGT gene base -> max allele number in reference
  ref_allele_max <- list()
  if (nrow(tbl_ref) > 0L) {
    tbl_ref_tmp <- as.data.table(tbl_ref)
    tbl_ref_tmp[, gene_base := sub("[*].*$", "", imgt_allele)]
    tbl_ref_tmp[, allele_n  := suppressWarnings(
      as.integer(sub(".*[*]", "", imgt_allele)))]
    for (gb in unique(tbl_ref_tmp$gene_base))
      ref_allele_max[[gb]] <- as.integer(max(tbl_ref_tmp[gene_base == gb, allele_n],
                                             na.rm = TRUE))
  }

  next_allele_n <- list()
  .next_allele <- function(gene_base) {
    if (is.null(next_allele_n[[gene_base]])) {
      base_max <- as.integer(if (!is.null(ref_allele_max[[gene_base]]))
        ref_allele_max[[gene_base]] else 0L)
      next_allele_n[[gene_base]] <<- base_max + 1L
    } else {
      next_allele_n[[gene_base]] <<- as.integer(next_allele_n[[gene_base]]) + 1L
    }
    sprintf("%02d", as.integer(next_allele_n[[gene_base]]))
  }

  # In ASC mode, use PIgLET new_allele directly (IGHVFx-Gy*01 style).
  # PIgLET has already removed exact duplicates (removed_duplicated=TRUE rows
  # are dropped here). No IMGT-name inference, no DECIPHER, no distance lookup.
  if (isTRUE(opt$use_asc)) {
    tbl_cust[, new_allele_final := normalise_name(new_allele)]

if (isTRUE(opt$trig_nc)) {
      # TRIG-NC FASTA naming is applied uniformly to the MERGED custom+reference
      # set in annotate_locus() via .trignc_rename_merged(), so that reference
      # and novel genes are size-ranked together. Here we only keep the ASC
      # new_allele in the annotation table for provenance; the authoritative
      # colon names are assigned downstream on the merged FASTA.
      cat(sprintf("  [TRIG-NC] %s naming deferred to merged rename (uniform gene sizing)\n",
                  if (nrow(tbl_cust)) as.character(tbl_cust$imgt_allele[1]) else ""))
    }

    tbl_cust[, new_allele       := new_allele_final]
    drop_cols <- intersect(c("new_allele_final","ref_rep","fam_rep"), names(tbl_cust))
    if (length(drop_cols)) tbl_cust[, (drop_cols) := NULL]
  } else {

  tbl_cust[, new_allele_final := new_allele]
  for (i in seq_len(nrow(tbl_cust))) {
    ref_rep <- tbl_cust$ref_rep[i]
    fam_rep <- if ("fam_rep" %in% names(tbl_cust)) tbl_cust$fam_rep[i] else NA_character_

    # Pick the best reference gene name available
    best_ref <- if (!is.na(ref_rep) && nchar(ref_rep) > 0L)  ref_rep
                else if (!is.na(fam_rep) && nchar(fam_rep) > 0L) fam_rep
                else NA_character_

    if (!is.na(best_ref)) {
      gene_base <- sub("[*].*$", "", best_ref)
      # Guard: ref allele had no * -> treat as allele 0
      if (!grepl("[*]", best_ref)) ref_allele_max[[gene_base]] <- 0L
      tbl_cust$new_allele_final[i] <-
        paste0(gene_base, "*", .next_allele(gene_base))
    } else {
      # No reference co-member at cluster or family level.
      # Use the PIgLET distance matrix to find the closest REFERENCE allele.
      # dist_mat[custom_name, ref_name] gives the similarity distance from
      # the joint clustering — exactly what PIgLET used to build clusters.
      cust_nm  <- tbl_cust$imgt_allele[i]
      best_gene <- NULL
      if (!is.null(dist_mat) && cust_nm %in% rownames(dist_mat)) {
        row       <- dist_mat[cust_nm, , drop = TRUE]
        # Only use non-_dup reference sequences as anchors
        ref_cols  <- names(row)[names(row) %in% ref_names_valid]
        if (length(ref_cols) > 0L) {
          # Prefer alleles with proper IMGT subtype (dash or S notation)
          subtyped_cols <- ref_cols[grepl("[-S][0-9]", sub("[*].*$","",ref_cols))]
          use_cols <- if (length(subtyped_cols) > 0L) subtyped_cols else ref_cols
          closest     <- use_cols[which.min(row[use_cols])]
          # Recover original name for gene_base (cleaned may have lost subtype)
          orig_closest <- if (exists("rev_ref_name_map") && closest %in% names(rev_ref_name_map))
                            rev_ref_name_map[[closest]] else closest
          best_gene <- sub("[*].*$", "", orig_closest)
        }
      }
      # Fallback: ungapped edit distance if dist_mat unavailable
      if (is.null(best_gene)) {
        seq_i     <- as.character(DECIPHER::RemoveGaps(
                       custom_gapped[cust_nm], removeGaps = "all"))
        # Exclude _dup sequences from adist reference pool
        ref_valid2 <- ref_for_joint[!grepl("_dup", names(ref_for_joint), fixed=TRUE)]
        ref_ug     <- as.character(DECIPHER::RemoveGaps(ref_valid2, removeGaps = "all"))
        dists      <- as.integer(adist(seq_i, ref_ug))
        best_cleaned  <- names(ref_valid2)[which.min(dists)]
        orig_closest  <- if (exists("rev_ref_name_map") && best_cleaned %in% names(rev_ref_name_map))
                           rev_ref_name_map[[best_cleaned]] else best_cleaned
        best_gene <- sub("[*].*$", "", orig_closest)
      }
      tbl_cust$new_allele_final[i] <-
        paste0(best_gene, "*", .next_allele(best_gene))
    }
  }
  tbl_cust[, new_allele := new_allele_final]

  # --- 5b. Sanitise any remaining _dup names ----------------------------------
  # _dup suffixes are internal deduplication markers that must NEVER appear in
  # the final output. Any allele whose new_allele still contains "_dup" has
  # missed the normal naming path — resolve it now via:
  #   (1) dist_mat closest reference (preferred — same distance PIgLET used)
  #   (2) ungapped edit distance fallback
  # Then assign the next available allele number on that gene base.
  dup_mask <- grepl("_dup", tbl_cust$new_allele, fixed = TRUE)
  if (any(dup_mask)) {
    message(sprintf(
      "  [FIX] %s: %d allele(s) still have _dup names — resolving via distance",
      label, sum(dup_mask)))
    for (i in which(dup_mask)) {
      cust_nm  <- tbl_cust$imgt_allele[i]
      best_gene <- NULL

      # (1) dist_mat lookup against reference sequences
      if (!is.null(dist_mat) && cust_nm %in% rownames(dist_mat)) {
        row      <- dist_mat[cust_nm, , drop = TRUE]
        # Exclude _dup sequences from candidate references
        ref_cols <- names(row)[names(row) %in% ref_names_valid]
        if (length(ref_cols) > 0L) {
          # Prefer alleles with a proper IMGT subtype name (dash or S notation)
          # e.g. IGHV1-2*01 or IGHV1S16*01 are preferred over bare IGHV1*68
          subtyped_cols <- ref_cols[grepl("[-S][0-9]", sub("[*].*$","",ref_cols))]
          use_cols <- if (length(subtyped_cols) > 0L) subtyped_cols else ref_cols
          closest      <- use_cols[which.min(row[use_cols])]
          orig_closest <- if (exists("rev_ref_name_map") &&
                              closest %in% names(rev_ref_name_map))
                            rev_ref_name_map[[closest]] else closest
          best_gene <- sub("[*].*$", "", orig_closest)
          message(sprintf(
            "  [FIX]   %s -> closest ref %s (dist=%.4f) -> gene base: %s",
            cust_nm, closest, row[closest], best_gene))
        }
      }

      # (2) ungapped edit distance fallback
      if (is.null(best_gene)) {
        seq_i    <- as.character(DECIPHER::RemoveGaps(custom_gapped[cust_nm], removeGaps="all"))
        # Exclude _dup sequences from adist reference pool
        ref_valid  <- ref_for_joint[!grepl("_dup", names(ref_for_joint), fixed=TRUE)]
        ref_ug     <- as.character(DECIPHER::RemoveGaps(ref_valid, removeGaps="all"))
        dists      <- as.integer(adist(seq_i, ref_ug))
        best_cl    <- names(ref_valid)[which.min(dists)]
        orig_cl    <- if (exists("rev_ref_name_map") && best_cl %in% names(rev_ref_name_map))
                        rev_ref_name_map[[best_cl]] else best_cl
        best_gene  <- sub("[*].*$", "", orig_cl)
        message(sprintf("  [FIX]   %s -> adist closest ref %s -> gene base: %s",
                        cust_nm, orig_cl, best_gene))
      }

      if (!is.null(best_gene)) {
        new_name <- paste0(best_gene, "*", .next_allele(best_gene))
        message(sprintf("  [FIX]   %s: '%s' -> '%s'",
                        label, tbl_cust$new_allele[i], new_name))
        tbl_cust$new_allele[i] <- new_name
      } else {
        message(sprintf("  [WARN] %s: could not resolve _dup name for '%s'",
                        label, cust_nm))
      }
    }
  }

  drop_cols <- intersect(c("new_allele_final","ref_rep","fam_rep"), names(tbl_cust))
  if (length(drop_cols)) tbl_cust[, (drop_cols) := NULL]

  } # end IMGT-naming else branch


  # --- 6. Rename DNAStringSet objects --------------------------------------
  annot_df <- data.frame(imgt_allele = tbl_cust$imgt_allele,
                         new_allele  = tbl_cust$new_allele,
                         stringsAsFactors = FALSE)
  renamed_custom <- rename_by_table(custom_gapped, annot_df, label = label)

  # --- 7. Assemble full annotation table -----------------------------------
  cols  <- intersect(c("imgt_allele","new_allele","cluster_id","family_cluster"),
                     names(tbl_ref))
  cols2 <- intersect(c("imgt_allele","new_allele","piglet_cluster",
                        "cluster_id","family_cluster"),
                     names(tbl_cust))
  full_annot <- rbindlist(list(
    tbl_ref[,  ..cols ][,  source := "reference"],
    tbl_cust[, ..cols2][, source := "custom"]
  ), fill = TRUE)

  # Strip internal _dup<N> disambiguation tags from all new_allele values.
  # These arise from identical bare-family reference sequences (e.g. several
  # ">IGHV1" records) and are never valid output names. In TRIG-NC mode the
  # merged rename re-derives names anyway; in other modes the base name is the
  # correct output. Strip rather than error so a clean reference _dup does not
  # abort the whole run.
  full_annot[, new_allele := sub("_dup[0-9]+$", "", new_allele)]

  # Guard only for CUSTOM sequences whose name still contains _dup after
  # stripping the numeric suffix (would indicate a genuine logic error).
  cust_bad <- full_annot[source == "custom" & grepl("_dup", new_allele, fixed = TRUE)]
  if (nrow(cust_bad) > 0L) {
    message(sprintf(
      "  [ERROR] %s: %d custom _dup name(s) survived: %s",
      label, nrow(cust_bad),
      paste(head(cust_bad$new_allele, 5L), collapse=", ")))
    stop(sprintf("custom _dup allele names in final output for locus %s.", label))
  }

  list(renamed_custom    = renamed_custom,
       annot_table       = full_annot,
       dropped           = pf_res$dropped,
       dist_mat          = dist_mat,
       ref_names_used    = ref_names_clean,
       custom_names_used = custom_names)
}

# =============================================================================
# SECTION 5 -- File discovery
# =============================================================================
find_custom <- function(segment, locus, custom_dir) {
  candidates <- c(
    file.path(custom_dir, "heavy", paste0(locus, segment, ".fasta")),
    file.path(custom_dir, "light", paste0(locus, segment, ".fasta")),
    file.path(custom_dir, paste0(locus, segment, ".fasta")),
    file.path(custom_dir, paste0("imgt_custom_", locus, segment, ".fasta"))
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) return(NULL)
  found[1L]
}
find_ref <- function(segment, locus, ref_dir, species) {
  candidates <- c(
    file.path(ref_dir, paste0("imgt_", species, "_", locus, segment, ".fasta")),
    file.path(ref_dir, paste0(locus, segment, ".fasta")),
    file.path(ref_dir, species, "vdj",
              paste0("imgt_", species, "_", locus, segment, ".fasta"))
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) return(NULL)
  found[1L]
}

cat("--- Locating input files ---\n")
loci <- list(
  IGHV = list(locus = "IGH", seg = "V"),
  IGHD = list(locus = "IGH", seg = "D"),
  IGHJ = list(locus = "IGH", seg = "J"),
  IGKV = list(locus = "IGK", seg = "V"),
  IGKJ = list(locus = "IGK", seg = "J"),
  IGLV = list(locus = "IGL", seg = "V"),
  IGLJ = list(locus = "IGL", seg = "J")
)
custom_paths <- list()
ref_paths    <- list()
for (nm in names(loci)) {
  cp <- find_custom(loci[[nm]]$seg, loci[[nm]]$locus, opt$custom_dir)
  rp <- find_ref(loci[[nm]]$seg,    loci[[nm]]$locus, opt$ref_dir, opt$species)
  custom_paths[[nm]] <- cp;  ref_paths[[nm]] <- rp
  cat(sprintf("  %-6s  custom: %-55s  ref: %s\n", nm,
              ifelse(is.null(cp),"(absent)",cp),
              ifelse(is.null(rp),"(absent)",rp)))
}

# =============================================================================
# SECTION 6 -- Load and normalise
# =============================================================================
cat("\n--- Loading and normalising sequences ---\n")
custom_norm <- setNames(
  lapply(names(loci), function(nm) read_fasta_safe(custom_paths[[nm]], nm, "custom")),
  names(loci))
ref_norm <- setNames(
  lapply(names(loci), function(nm) read_fasta_safe(ref_paths[[nm]], nm, "reference")),
  names(loci))

cseqs <- function(nm) if (is.null(custom_norm[[nm]])) NULL else custom_norm[[nm]]$seqs
rseqs <- function(nm) if (is.null(ref_norm[[nm]]))    NULL else ref_norm[[nm]]$seqs

# Collect header maps for annotation output
header_maps <- list()
for (nm in names(loci)) {
  for (tag in c("custom","reference")) {
    obj <- if (tag == "custom") custom_norm[[nm]] else ref_norm[[nm]]
    if (!is.null(obj) && nrow(obj$map) > 0L) {
      m <- copy(obj$map); m$locus <- nm; m$source <- tag
      header_maps[[paste0(nm,"_",tag)]] <- m
    }
  }
}

# =============================================================================
# SECTION 6b -- Write debug script immediately after path resolution
# (Written here so the file exists even if the pipeline crashes later)
# =============================================================================
cat("\n--- Writing debug R script ---\n")

debug_script_path <- file.path(opt$outdir, paste0(opt$prefix, "_debug_session.R"))

.write_debug_script <- function() {
  # Paths that are defined later (aux, annotations) are approximated here;
  # the debug script comments make clear which files appear after a successful run.
  ann_dir   <- file.path(opt$outdir, "annotations")
  gapped_dir_d <- file.path(opt$outdir, "germlines", "gapped")

  lines <- c(
    "# ================================================================",
    "# Debug / interactive session script",
    "# Auto-generated by piglet_annotate_and_build.R",
    "# Source this file in R to reproduce the pipeline interactively.",
    "# All input paths are pre-populated from the last pipeline run.",
    "# ================================================================",
    "",
    "suppressPackageStartupMessages({",
    "  library(data.table); library(DECIPHER)",
    "  library(piglet);     library(Biostrings)",
    "})",
    "",
    "# ---- Parameters ----",
    sprintf('custom_dir    <- "%s"', opt$custom_dir),
    sprintf('ref_dir       <- "%s"', opt$ref_dir),
    sprintf('species       <- "%s"', opt$species),
    sprintf('outdir        <- "%s"', opt$outdir),
    sprintf('igdata        <- "%s"', opt$igdata),
    sprintf('prefix        <- "%s"', opt$prefix),
    sprintf('v_trim3       <- %d',   opt$v_trim3prime),
    sprintf('j_trim3       <- %d',   opt$j_trim3prime),
    sprintf('fam_thresh    <- %g',   opt$family_threshold),
    sprintf('allele_thresh <- %g',   opt$allele_cluster_threshold),
    "",
    "# ---- Resolved input file paths ----"
  )

  for (nm in names(loci)) {
    cp <- custom_paths[[nm]]
    rp <- ref_paths[[nm]]
    lines <- c(lines,
      sprintf('custom_path_%s <- %s', nm,
              if (is.null(cp)) "NULL" else sprintf('"%s"', cp)),
      sprintf('ref_path_%s    <- %s', nm,
              if (is.null(rp)) "NULL" else sprintf('"%s"', rp))
    )
  }

  lines <- c(lines, "",
    "# ---- Load sequences (mirrors original working code) ----",
    "# Original pattern (from exploratory code):",
    "#   mrl_v <- readDNAStringSet(path)",
    "#   names(mrl_v) <- gsub(' ', '', names(mrl_v))",
    "",
    "read_seq <- function(path) {",
    "  if (is.null(path) || !file.exists(path)) return(NULL)",
    "  s <- readDNAStringSet(path)",
    "  names(s) <- gsub(' ', '', names(s))  # strip spaces from IMGT headers",
    "  s",
    "}",
    ""
  )

  for (nm in names(loci)) {
    lines <- c(lines,
      sprintf("custom_%s <- read_seq(custom_path_%s)", nm, nm),
      sprintf("ref_%s    <- read_seq(ref_path_%s)",    nm, nm)
    )
  }

  lines <- c(lines,
    "",
    "# ---- PIgLET: the critical named-vector fix ----",
    "# BUG:  as.character(DNAStringSet) drops names",
    "#       -> germ.dist subscript out of bounds",
    "# FIX:  setNames(as.character(seqs), names(seqs))",
    "",
    "run_piglet_debug <- function(seqs_gapped, trim3, label = '') {",
    "  # Remove gap characters so PIgLET receives ungapped sequences",
    "  seqs_ungapped <- DECIPHER::removeGaps(seqs_gapped, removeGaps = \"all\")",
    "  named_vec <- setNames(as.character(seqs_ungapped), names(seqs_ungapped))",
    "  cat(sprintf(\"Running PIgLET on %d seqs (%s)\\n\", length(named_vec), label))",
    "  piglet::inferAlleleClusters(",
    "    germline_set             = named_vec,",
    "    trim_3prime_side         = trim3,",
    "    mask_5prime_side         = 0L,",
    "    family_threshold         = fam_thresh,",
    "    allele_cluster_threshold = allele_thresh",
    "  )",
    "}",
    "",
    "# Minimal sanity check -- run PIgLET on custom IGHV alone:",
    "if (!is.null(custom_IGHV)) {",
    "  cat(\"Testing PIgLET on custom IGHV...\\n\")",
    "  asc_ighv <- run_piglet_debug(custom_IGHV, trim3 = v_trim3, label = \"IGHV\")",
    "  print(head(asc_ighv@alleleClusterTable))",
    "}",
    "",
    "# ---- Inspect annotation outputs (available after successful run) ----",
    sprintf('# header_map    <- fread("%s")',
            file.path(ann_dir, paste0(opt$prefix,"_header_normalisation_map.tsv"))),
    sprintf('# cluster_annot <- fread("%s")',
            file.path(ann_dir, paste0(opt$prefix,"_allele_cluster_annotation.tsv"))),
    sprintf('# provenance    <- fread("%s")',
            file.path(ann_dir, paste0(opt$prefix,"_full_provenance.tsv"))),
    "",
    "# Trace one sequence's journey from raw header to final cluster name:",
    "# provenance[raw_header %like% 'IGHV1-2']",
    "",
    "# ---- Load hybrid output FASTAs (available after successful run) ----"
  )

  for (nm in c("IGHV","IGKV","IGLV","IGHJ","IGKJ","IGLJ","IGHD")) {
    lines <- c(lines,
      sprintf('# hybrid_%s <- readDNAStringSet("%s")', nm,
              file.path(gapped_dir_d, paste0(imgt_file_prefix,"_",nm,".fasta")))
    )
  }

  lines <- c(lines, "")
  writeLines(lines, debug_script_path)
  cat(sprintf("  Wrote debug script: %s\n", debug_script_path))
}

.write_debug_script()

# =============================================================================
# SECTION 7 -- Per-locus annotation
# =============================================================================

# Helper: annotate a single V or J locus
annotate_locus <- function(nm, trim3) {
  cust <- cseqs(nm)
  ref  <- rseqs(nm)

  if (is.null(cust) && is.null(ref)) {
    message(sprintf("  [WARN] No sequences for %s", nm))
    return(list(custom_gapped = NULL, ref_gapped = NULL,
                hybrid_gapped = NULL, hybrid_ungapped = NULL,
                annot = .empty_annot(), dropped = data.table()))
  }
  if (is.null(cust)) {
    cat(sprintf("  %s: no custom seqs; reference only\n", nm))
    ref_use <- ref
    if (isTRUE(opt$trig_nc)) ref_use <- .trignc_rename_merged(ref, ref, nm)
    return(list(custom_gapped = NULL, ref_gapped = ref_use,
                hybrid_gapped = ref_use, hybrid_ungapped = ungap(ref_use),
                annot = .empty_annot(), dropped = data.table()))
  }
  if (is.null(ref)) {
    cat(sprintf("  %s: no reference; PIgLET de-novo only\n", nm))
    pf  <- run_piglet(cust, trim3 = trim3,
                      fam_thresh = opt$family_threshold,
                      allele_thresh = opt$allele_cluster_threshold, label = nm)
    cr  <- rename_by_table(cust, pf$tbl, label = nm)
    return(list(custom_gapped = cr, ref_gapped = NULL,
                hybrid_gapped = cr, hybrid_ungapped = ungap(cr),
                annot = pf$tbl[, source := "custom"],
                dropped = pf$dropped))
  }

  res <- annotate_custom_with_ref(cust, ref,
                                  trim3         = trim3,
                                  fam_thresh    = opt$family_threshold,
                                  allele_thresh = opt$allele_cluster_threshold,
                                  label = nm)
  # Strip internal _dup<N> tags from reference names before merging. Identical
  # bare-family reference records (e.g. many ">IGHV1") were disambiguated with
  # _dup suffixes at load time; those are not valid germline names. After
  # stripping, collapse any exact name+sequence duplicates that result.
  ref_clean <- ref
  rn <- sub("_dup[0-9]+$", "", names(ref_clean))
  names(ref_clean) <- rn
  # Drop entries that are now exact name duplicates (keep first)
  ref_clean <- ref_clean[!duplicated(names(ref_clean))]

  hg <- merge_with_priority(res$renamed_custom, ref_clean)
  # TRIG-NC: re-derive uniform colon names across the MERGED set so reference
  # and novel genes are size-ranked TOGETHER within each family (avoids gene
  # number collisions between IMGT ref genes and size-ranked novel genes).
  # Pass PIgLET's annotation table so truly-novel sequences are grouped into
  # families by the 75% family-clustering threshold, not one family each.
  if (isTRUE(opt$trig_nc)) {
    hg <- .trignc_rename_merged(hg, ref_clean, nm, annot = res$annot_table)
  }
  list(custom_gapped    = res$renamed_custom,
       ref_gapped       = ref,
       hybrid_gapped    = hg,
       hybrid_ungapped  = ungap(hg),
       annot            = res$annot_table,
       dropped          = res$dropped,
       dist_mat         = res$dist_mat,
       ref_names_used   = res$ref_names_used,
       custom_names_used = res$custom_names_used)
}

cat("\n--- V-gene family annotation (joint PIgLET clustering) ---\n")
# Accumulator for TRIG-NC before->after name mappings (populated by
# .trignc_rename_merged as each locus is processed).
assign(".trignc_map_acc",
       data.table(before = character(0L), after = character(0L), locus = character(0L)),
       envir = .GlobalEnv)
V_results <- list()
for (nm in c("IGHV","IGKV","IGLV")) {
  cat(sprintf("\nProcessing %s...\n", nm))
  V_results[[nm]] <- annotate_locus(nm, trim3 = opt$v_trim3prime)
}

cat("\n--- J-gene annotation (joint PIgLET clustering) ---\n")
J_results <- list()
for (nm in c("IGHJ","IGKJ","IGLJ")) {
  cat(sprintf("\nProcessing %s...\n", nm))
  J_results[[nm]] <- annotate_locus(nm, trim3 = opt$j_trim3prime)
}

cat("\n--- D genes (merge, custom priority) ---\n")
D_gapped <- merge_with_priority(cseqs("IGHD"), rseqs("IGHD"))
if (!is.null(D_gapped)) names(D_gapped) <- normalise_name(names(D_gapped))

# =============================================================================
# SECTION 8 -- Write annotation tables
# =============================================================================
cat("\n--- Writing annotation tables ---\n")

# 8a. Header normalisation map: raw_header -> normalised_name -> parse_flag
#     This is the primary provenance table linking every sequence back to its
#     original FASTA header (IMGT, OGRDB, pipe-format, etc.)
header_map_dt  <- if (length(header_maps) > 0L)
  rbindlist(header_maps, fill = TRUE) else data.table()
header_map_out <- file.path(opt$outdir, "annotations",
                             paste0(opt$prefix, "_header_normalisation_map.tsv"))
fwrite(header_map_dt, header_map_out, sep = "\t")
cat(sprintf("  Wrote header map        : %s\n", header_map_out))

# 8b. PIgLET allele cluster annotation: normalised_name -> final cluster name
all_annot <- list()
for (nm in c("IGHV","IGKV","IGLV")) {
  a <- V_results[[nm]]$annot
  if (!is.null(a) && nrow(a) > 0L) { a <- copy(a); a$locus <- nm; all_annot[[nm]] <- a }
}
for (nm in c("IGHJ","IGKJ","IGLJ")) {
  a <- J_results[[nm]]$annot
  if (!is.null(a) && nrow(a) > 0L) { a <- copy(a); a$locus <- nm; all_annot[[nm]] <- a }
}
combined_annot <- rbindlist(all_annot, fill = TRUE)
annot_out      <- file.path(opt$outdir, "annotations",
                             paste0(opt$prefix, "_allele_cluster_annotation.tsv"))
fwrite(combined_annot, annot_out, sep = "\t")
cat(sprintf("  Wrote cluster annotation: %s\n", annot_out))

# 8c. Nearest-reference distance lookup table
# For each novel custom allele, records which reference allele was closest
# in the PIgLET joint distance matrix and the distance value.
# This provides the full audit trail for every IMGT-name assignment decision.
dist_lookup_rows <- list()
for (.nm in c("IGHV","IGKV","IGLV","IGHJ","IGKJ","IGLJ")) {
  .res <- if (.nm %in% c("IGHV","IGKV","IGLV")) V_results[[.nm]]
           else J_results[[.nm]]
  .dm    <- .res$dist_mat
  .refs  <- .res$ref_names_used
  .custs <- .res$custom_names_used
  if (is.null(.dm) || is.null(.refs) || is.null(.custs)) next
  for (.cn in intersect(.custs, rownames(.dm))) {
    .row   <- .dm[.cn, , drop = TRUE]
    .rcols <- names(.row)[names(.row) %in% .refs]
    if (!length(.rcols)) next
    .best  <- .rcols[which.min(.row[.rcols])]
    .dist  <- min(.row[.rcols], na.rm = TRUE)
    dist_lookup_rows[[length(dist_lookup_rows)+1L]] <-
      data.table(locus             = .nm,
                 custom_imgt_allele = .cn,
                 closest_ref_imgt  = .best,
                 piglet_distance   = .dist)
  }
}
if (length(dist_lookup_rows) > 0L) {
  dist_lookup     <- rbindlist(dist_lookup_rows)
  dist_lookup_out <- file.path(opt$outdir, "annotations",
                               paste0(opt$prefix, "_nearest_ref_lookup.tsv"))
  fwrite(dist_lookup, dist_lookup_out, sep = "\t")
  cat(sprintf("  Wrote nearest-ref lookup : %s\n", dist_lookup_out))
}

# 8c. Full provenance: join header map + cluster annotation into one table
#     Columns present depend on whether clustering was performed (cluster_id
#     and family_cluster are NA for loci where custom == reference).
if (nrow(header_map_dt) > 0L && nrow(combined_annot) > 0L) {
  # Select only columns that actually exist in combined_annot
  annot_cols <- intersect(
    c("imgt_allele","new_allele","cluster_id","family_cluster","source","locus"),
    names(combined_annot)
  )
  provenance <- merge(
    header_map_dt,
    combined_annot[, annot_cols, with = FALSE],
    by.x = c("normalised_name","locus","source"),
    by.y = c("imgt_allele","locus","source"),
    all.x = TRUE
  )
} else {
  provenance <- copy(header_map_dt)
}

# Attach the FINAL FASTA name (TRIG-NC colon format when --trig_nc, else the
# IMGT/ASC name). Join through the TRIG-NC map (normalised IMGT -> colon), so
# every row — including reference rows with no clustering annotation — carries
# its final name. This is the authoritative original -> final liftover.
trignc_map2 <- if (exists(".trignc_map_acc", envir = .GlobalEnv))
  get(".trignc_map_acc", envir = .GlobalEnv) else
  data.table(before = character(0L), after = character(0L), locus = character(0L))

if (isTRUE(opt$trig_nc) && nrow(trignc_map2) > 0L) {
  tnc2 <- unique(trignc_map2[, .(normalised_name = before, final_name = after, locus)])
  provenance <- merge(provenance, tnc2, by = c("normalised_name","locus"), all.x = TRUE)
  # Second-chance join for names carrying an internal _dup<N> tag: strip it and
  # retry against the TRIG-NC map (which uses clean names).
  if (any(is.na(provenance$final_name))) {
    provenance[, .nn_clean := sub("_dup[0-9]+$", "", normalised_name)]
    tnc3 <- unique(trignc_map2[, .(.nn_clean = before, final2 = after, locus)])
    provenance <- merge(provenance, tnc3, by = c(".nn_clean","locus"), all.x = TRUE)
    provenance[is.na(final_name) & !is.na(final2), final_name := final2]
    provenance[, c(".nn_clean","final2") := NULL]
  }
  # IGHD is not TRIG-NC renamed (D genes keep IMGT names); and any remaining
  # unmatched rows fall back to their clean normalised name.
  provenance[is.na(final_name), final_name := sub("_dup[0-9]+$", "", normalised_name)]
} else {
  # Non-TRIG: final name is new_allele where present, else normalised_name
  if ("new_allele" %in% names(provenance)) {
    provenance[, final_name := new_allele]
    provenance[is.na(final_name) | final_name == "", final_name := normalised_name]
  } else {
    provenance[, final_name := normalised_name]
  }
}

setcolorder(provenance,
  intersect(c("raw_header","normalised_name","new_allele","final_name",
              "locus","source","parse_flag","cluster_id","family_cluster"),
            names(provenance)))

prov_out <- file.path(opt$outdir, "annotations",
                      paste0(opt$prefix, "_full_provenance.tsv"))
fwrite(provenance, prov_out, sep = "\t")
cat(sprintf("  Wrote full provenance    : %s  (%d rows, %d with final name)\n",
            prov_out, nrow(provenance),
            sum(!is.na(provenance$final_name) & provenance$final_name != "")))

# 8d. Pre-flight dropped sequences
all_dropped <- rbindlist(c(
  lapply(V_results, `[[`, "dropped"),
  lapply(J_results, `[[`, "dropped")
), fill = TRUE, idcol = "locus")
if (nrow(all_dropped) > 0L) {
  dropped_out <- file.path(opt$outdir, "annotations",
                            paste0(opt$prefix, "_dropped_sequences.tsv"))
  fwrite(all_dropped, dropped_out, sep = "\t")
  cat(sprintf("  Wrote dropped sequences  : %s  (%d seqs)\n",
              dropped_out, nrow(all_dropped)))
}


#' Post-normalise a DNAStringSet:
#'   1. Strip _dup<N> suffixes
#'   2. Infer subfamily for bare-family names (IGHV1 -> IGHV1-2 etc.) using
#'      DECIPHER alignment + distance matrix against same-family sequences
#'      that already have a subfamily annotation.
#'   3. Ensure every name has *NN allele number
#'   4. Re-number colliding names sequentially
#'   5. Deduplicate by ungapped content
#'
#' Returns list(seqs = DNAStringSet, name_map = data.table(before, after))
#' so downstream code can reconcile final FASTA names back to original IDs.
post_normalise_seqs <- function(seqs, label = "", use_asc = FALSE) {
  if (is.null(seqs) || length(seqs) == 0L)
    return(list(seqs = seqs, name_map = data.table(before=character(), after=character())))

  old_names <- names(seqs)

  # ---- Step 1: strip _dup<N> suffixes ----
  new_names <- sub("_dup\\d+$", "", old_names)

  # ---- Step 2: DECIPHER-based subfamily inference for bare-family names ----
  # A bare-family name has the locus prefix and a number but no dash-subtype:
  #   IGHV1  (no dash after family number)
  # A properly annotated name has: IGHV1-2  or  IGHV1-18  etc.
  # Strategy:
  #   - Find all sequences in this set that have a subfamily (have a dash after
  #     the family number, e.g. IGHV1-2).
  #   - For each bare-family sequence, restrict comparison to same-family
  #     sequences WITH a subfamily annotation.
  #   - Align with DECIPHER::AlignSeqs, compute DistanceMatrix.
  #   - Pick the closest annotated sequence -> adopt its gene name as base.
  #   - If no annotated same-family sequences exist, fall through to *01.

  # Identify bare-family indices: name matches ^IG[HKL][VDJ]<digits>$ (no dash, no *)
  bare_mask <- grepl("^IG[HKL][VDJ]\\d+$", new_names)

  if (any(bare_mask) && !isTRUE(use_asc)) {
    n_bare <- sum(bare_mask)
    message(sprintf("  [POST-NORM] %s: %d bare-family names -> DECIPHER subfamily inference (by family group)",
                    label, n_bare))

    # Helper: extract locus prefix and family number from any gene name
    # e.g. IGHV1-2*01 -> locus="IGHV", fam="1"
    #      IGHV1       -> locus="IGHV", fam="1"
    get_locus <- function(nm) sub("^(IG[HKL][VDJ]).*$", "\\1", nm)
    get_fam   <- function(nm) sub("^IG[HKL][VDJ](\\d+).*$", "\\1", nm)

    # Identify annotated sequences: have a dash after the family number
    # e.g. IGHV1-2*01, IGHV1-18*03 are annotated; IGHV1 is not
    annotated_mask <- grepl("^IG[HKL][VDJ]\\d+-", new_names)

    # Group bare sequences by (locus, family) -> build ONE matrix per group
    bare_idx <- which(bare_mask)
    fam_groups <- unique(paste0(get_locus(new_names[bare_idx]),
                                get_fam(new_names[bare_idx])))

    for (fg in fam_groups) {
      locus_pfx <- sub("^(IG[HKL][VDJ]).*$", "\\1", fg)
      fam_n     <- sub("^IG[HKL][VDJ]", "",         fg)

      # Indices of bare sequences in this family group
      grp_bare_idx <- bare_idx[
        get_locus(new_names[bare_idx]) == locus_pfx &
        get_fam(new_names[bare_idx])   == fam_n
      ]

      # Indices of annotated same-family sequences (IGHV1-<x>) anywhere in the set
      sub_pattern   <- paste0("^", locus_pfx, fam_n, "-")
      grp_annot_idx <- which(grepl(sub_pattern, new_names) & annotated_mask)

      if (length(grp_annot_idx) == 0L) {
        message(sprintf("  [POST-NORM] %s: family %s%s: no annotated references; leaving bare",
                        label, locus_pfx, fam_n))
        next
      }

      # Build one DNAStringSet: all bare members + all annotated same-family seqs.
      # De-duplicate names within the candidate set (identical names break row lookup).
      grp_all_idx <- c(grp_bare_idx, grp_annot_idx)
      grp_seqs    <- seqs[grp_all_idx]
      names(grp_seqs) <- new_names[grp_all_idx]
      # If duplicate names exist (e.g. multiple bare IGHV1), make them unique
      nms_tmp <- names(grp_seqs)
      dup_nms <- duplicated(nms_tmp)
      if (any(dup_nms))
        nms_tmp[dup_nms] <- paste0(nms_tmp[dup_nms], ".", seq_len(sum(dup_nms)))
      names(grp_seqs) <- nms_tmp

      # Parallel name vectors for bare members (using de-duplicated names)
      n_bare_grp   <- length(grp_bare_idx)
      bare_tmp_nms <- nms_tmp[seq_len(n_bare_grp)]
      annot_nms    <- nms_tmp[seq(n_bare_grp + 1L, length(nms_tmp))]

      tryCatch({
        grp_ung <- DECIPHER::RemoveGaps(grp_seqs, removeGaps = "all")
        aligned <- DECIPHER::AlignSeqs(grp_ung, verbose = FALSE)
        dm      <- DECIPHER::DistanceMatrix(aligned, verbose = FALSE)

        message(sprintf("  [POST-NORM] %s: family %s%s: %d bare vs %d annotated refs",
                        label, locus_pfx, fam_n, n_bare_grp, length(annot_nms)))

        # For each bare member, find its closest annotated sequence
        for (k in seq_len(n_bare_grp)) {
          orig_i    <- grp_bare_idx[k]
          bare_tmp  <- bare_tmp_nms[k]
          orig_bare <- new_names[orig_i]   # e.g. "IGHV1"

          if (!bare_tmp %in% rownames(dm)) next
          row_k         <- dm[bare_tmp, , drop = TRUE]
          cand_annot    <- annot_nms[annot_nms %in% colnames(dm)]
          if (length(cand_annot) == 0L) next

          closest       <- cand_annot[which.min(row_k[cand_annot])]
          # The closest name is a de-duplicated annotated name; strip allele -> gene base
          # e.g. "IGHV1-2*01" -> "IGHV1-2"
          new_gene_base <- sub("[*].*$", "", closest)
          new_names[orig_i] <- new_gene_base
          message(sprintf("  [POST-NORM] %s: '%s' -> '%s' (dist=%.4f)",
                          label, orig_bare, new_gene_base, row_k[closest]))
        }
      }, error = function(e) {
        message(sprintf("  [POST-NORM] %s: DECIPHER failed for family %s%s: %s",
                        label, locus_pfx, fam_n, conditionMessage(e)))
      })
    }
  }

  # ---- Step 3: ensure every name has *NN ----
  # TRIG-NC colon-format names (IGKJ:01:002:001) already encode the allele as
  # the final colon field, so they must NOT receive an appended *01.
  is_trignc <- grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ]):[0-9]", new_names)
  no_star <- !grepl("[*]", new_names) & !is_trignc
  if (any(no_star)) {
    message(sprintf("  [POST-NORM] %s: %d names still missing allele number -> *01",
                    label, sum(no_star)))
    new_names[no_star] <- paste0(new_names[no_star], "*01")
  }

  # ---- Step 4: re-number colliding gene bases sequentially ----
  # Skip TRIG-NC colon names — their allele numbering is already resolved by
  # the size-ranked assignment and the colon format has no * to split on.
  gene_allele_count <- list()
  for (i in seq_along(new_names)) {
    nm        <- new_names[i]
    if (grepl("^(IG[HKL][VDJ]|TR[ABGD][VDJ]):[0-9]", nm)) next
    gene_base <- sub("[*].*$", "", nm)
    cur_n     <- suppressWarnings(as.integer(sub(".*[*]", "", nm)))
    if (is.na(cur_n)) cur_n <- 1L
    cnt <- gene_allele_count[[gene_base]]
    if (is.null(cnt)) {
      gene_allele_count[[gene_base]] <- list(max_used = cur_n)
    } else {
      next_n <- as.integer(cnt$max_used) + 1L
      new_names[i] <- paste0(gene_base, "*", sprintf("%02d", next_n))
      gene_allele_count[[gene_base]]$max_used <- next_n
    }
  }

  # ---- Step 5: build name map BEFORE dedup (maps pre-dedup -> final name) ----
  name_map_pre <- data.table(before = old_names, after = new_names)

  names(seqs) <- new_names

  # ---- Step 6: deduplicate by ungapped content ----
  ug_content <- as.character(DECIPHER::RemoveGaps(seqs, removeGaps = "all"))
  keep       <- !duplicated(ug_content)
  seqs       <- seqs[keep]
  name_map_pre <- name_map_pre[keep]

  # ---- Step 7: final name-uniqueness safety check ----
  final_keep <- !duplicated(names(seqs))
  seqs       <- seqs[final_keep]
  name_map_pre <- name_map_pre[final_keep]

  list(seqs = seqs, name_map = name_map_pre)
}

# =============================================================================
# SECTION 9 -- Write FASTAs
# =============================================================================
write_fasta <- function(seqs, path) {
  if (!is.null(seqs) && length(seqs) > 0L) {
    writeXStringSet(seqs, path)
    cat(sprintf("  Wrote: %s  (%d seqs)\n", path, length(seqs)))
  } else {
    message(sprintf("  [SKIP] No sequences for: %s", path))
  }
}

gapped_dir   <- file.path(opt$outdir, "germlines", "gapped")
ungapped_dir <- file.path(opt$outdir, "germlines", "ungapped")

cat("\n--- Writing gapped germline FASTAs ---\n")
# Accumulate post-normalisation name maps for the final provenance table
all_post_norm_maps <- list()

for (nm in c("IGHV","IGKV","IGLV")) {
  pn <- post_normalise_seqs(V_results[[nm]]$hybrid_gapped, label = nm, use_asc = isTRUE(opt$use_asc))
  V_results[[nm]]$hybrid_gapped   <- pn$seqs
  V_results[[nm]]$hybrid_ungapped <- ungap(pn$seqs)
  if (nrow(pn$name_map) > 0L) {
    pn$name_map$locus <- nm
    all_post_norm_maps[[nm]] <- pn$name_map
  }
  write_fasta(pn$seqs, file.path(gapped_dir, paste0(imgt_file_prefix,"_",nm,".fasta")))
}
for (nm in c("IGHJ","IGKJ","IGLJ")) {
  pn <- post_normalise_seqs(J_results[[nm]]$hybrid_gapped, label = nm, use_asc = isTRUE(opt$use_asc))
  J_results[[nm]]$hybrid_gapped <- pn$seqs
  if (nrow(pn$name_map) > 0L) {
    pn$name_map$locus <- nm
    all_post_norm_maps[[nm]] <- pn$name_map
  }
  write_fasta(pn$seqs, file.path(gapped_dir, paste0(imgt_file_prefix,"_",nm,".fasta")))
}
{
  pn_d    <- post_normalise_seqs(D_gapped, label = "IGHD", use_asc = isTRUE(opt$use_asc))
  D_gapped <- pn_d$seqs
  if (nrow(pn_d$name_map) > 0L) {
    pn_d$name_map$locus <- "IGHD"
    all_post_norm_maps[["IGHD"]] <- pn_d$name_map
  }
  write_fasta(D_gapped, file.path(gapped_dir, paste0(imgt_file_prefix,"_IGHD.fasta")))
}

# ── Constant region FASTAs (ASC/PIgLET path) ────────────────────────────────
# Constant genes are not clustered or renamed — they are copied verbatim from
# the reference (name-deduped), exactly as in --as-is-ids mode. Without this,
# the ig_c database is built from an empty FASTA in ASC mode.
cat("\n--- Constant region sequences (from reference, name-deduped) ---\n")
# Resolve the reference constant directory. opt$ref_dir is the VDJ dir; the
# constant FASTAs live in a sibling constant/ dir. Try several layouts.
.const_dirs <- unique(c(
  sub("/vdj/?$", "/constant", opt$ref_dir),
  file.path(dirname(opt$ref_dir), "constant"),
  opt$ref_dir))
cat(sprintf("  Searching constant dirs: %s\n", paste(.const_dirs, collapse=", ")))
for (c_locus in c("IGHC","IGKC","IGLC")) {
  cand <- unlist(lapply(.const_dirs, function(d) c(
    file.path(d, paste0("imgt_", opt$species, "_", c_locus, ".fasta")),
    file.path(d, paste0(opt$species, "_", c_locus, ".fasta")),
    file.path(d, paste0(c_locus, ".fasta")))))
  c_ref <- Filter(file.exists, cand)[1L]
  if (is.na(c_ref)) {
    cat(sprintf("  [SKIP] %s: no reference FASTA found\n", c_locus)); next
  }
  ok_write <- tryCatch({
    c_seqs <- Biostrings::readDNAStringSet(c_ref)
    if (length(c_seqs) == 0L) { cat(sprintf("  [SKIP] %s: empty file\n", c_locus)); NULL }
    else {
      # Parse IMGT pipe headers to field 2 (gene name); fall back to raw name.
      raw <- names(c_seqs)
      parsed <- vapply(raw, function(h) {
        if (grepl("\\|", h)) {
          f <- strsplit(h, "\\|", fixed = FALSE)[[1]]
          if (length(f) >= 2L && nzchar(f[2])) f[2] else h
        } else sub("\\s.*$", "", h)
      }, character(1L), USE.NAMES = FALSE)
      names(c_seqs) <- parsed
      keep <- !duplicated(parsed)
      n_dup <- sum(!keep)
      c_seqs <- c_seqs[keep]
      out_fa <- file.path(gapped_dir, paste0(imgt_file_prefix, "_", c_locus, ".fasta"))
      write_fasta(c_seqs, out_fa)
      cat(sprintf("  %s: %d sequences%s -> %s\n", c_locus, length(c_seqs),
                  if (n_dup > 0L) sprintf(" (%d dup name(s) removed)", n_dup) else "",
                  basename(out_fa)))
      TRUE
    }
  }, error = function(e) {
    cat(sprintf("  [ERROR] %s: %s\n", c_locus, conditionMessage(e))); NULL
  })
}

# Write the final name map: original raw header -> final sequence name in FASTA.
# This joins: raw_header -> normalised_name -> post-norm name (after subfamily
# inference and sequential numbering) -- one row per sequence in the output FASTA.
if (length(all_post_norm_maps) > 0L) {
  post_norm_dt <- rbindlist(all_post_norm_maps, fill = TRUE)

  # In TRIG-NC mode the FASTA names were converted from IMGT/OGRDB to colon
  # format INSIDE annotate_locus (via .trignc_rename_merged), BEFORE
  # post_normalise_seqs ran. So post_norm_dt$before is already the colon name.
  # To map the ORIGINAL header through to the colon name we must insert the
  # TRIG-NC before(IMGT)->after(colon) step between normalised_name and the
  # post-norm map.
  trignc_map <- if (exists(".trignc_map_acc", envir = .GlobalEnv))
    get(".trignc_map_acc", envir = .GlobalEnv) else
    data.table(before = character(0L), after = character(0L), locus = character(0L))

  if (isTRUE(opt$trig_nc) && nrow(trignc_map) > 0L) {
    # Chain: raw_header -> normalised_name --(trignc)--> colon --(postnorm)--> final
    # Step A: normalised_name -> colon name (trignc_map: before=IMGT, after=colon)
    hdr <- header_map_dt[, .(raw_header, normalised_name, parse_flag, locus, source)]
    tnc <- unique(trignc_map[, .(normalised_name = before, colon_name = after, locus)])
    step_a <- merge(hdr, tnc, by = c("normalised_name","locus"), all.x = TRUE)
    # Where no TRIG-NC entry (shouldn't happen for V/J), colon = normalised
    step_a[is.na(colon_name), colon_name := normalised_name]
    # Step B: colon name -> post-norm final (post_norm_dt: before=colon, after=final)
    pnm <- post_norm_dt[, .(colon_name = before, final_fasta_name = after, locus)]
    final_name_map <- merge(step_a, pnm, by = c("colon_name","locus"), all.x = TRUE)
    final_name_map[is.na(final_fasta_name), final_fasta_name := colon_name]
    # Keep the intermediate colon name visible for provenance
    setnames(final_name_map, "colon_name", "trignc_name")
  } else {
    # Non-TRIG modes: original direct join normalised_name -> final
    final_name_map <- merge(
      header_map_dt[, .(raw_header, normalised_name, parse_flag, locus, source)],
      post_norm_dt[, .(normalised_name = before, final_fasta_name = after, locus)],
      by = c("normalised_name", "locus"),
      all.x = TRUE
    )
    final_name_map[is.na(final_fasta_name), final_fasta_name := normalised_name]
  }

  # For ASC / TRIG-NC analysis the map MUST let the user recover the ORIGINAL
  # input sequence id from the final FASTA name and vice versa. Provide explicit,
  # clearly-named columns:
  #   original_id      : the exact id from the user's input / reference FASTA
  #   normalised_name  : after header parsing (IMGT pipe -> allele, strain tag)
  #   final_fasta_name : the name written into the output FASTA (TRIG-NC colon
  #                      format when --trig_nc, else IMGT/ASC name)
  final_name_map[, original_id := raw_header]
  setcolorder(final_name_map,
    intersect(c("original_id","normalised_name","trignc_name","final_fasta_name",
                "locus","source","parse_flag","raw_header"),
              names(final_name_map)))
  final_map_out <- file.path(opt$outdir, "annotations",
                              paste0(opt$prefix, "_final_name_map.tsv"))
  fwrite(final_name_map[order(locus, source, original_id)],
         final_map_out, sep = "\t")
  cat(sprintf("  Wrote final name map     : %s  (%d sequences)\n",
              final_map_out, nrow(final_name_map)))
  cat(sprintf("    columns: original_id -> normalised_name -> final_fasta_name\n"))
  n_renamed <- sum(final_name_map$original_id != final_name_map$final_fasta_name, na.rm=TRUE)
  cat(sprintf("    %d of %d sequences renamed from their original id\n",
              n_renamed, nrow(final_name_map)))
}

combined_V <- do.call(c, Filter(Negate(is.null),
  lapply(c("IGHV","IGKV","IGLV"), function(nm) V_results[[nm]]$hybrid_gapped)))
write_fasta(combined_V, file.path(gapped_dir, paste0(imgt_file_prefix,"_ALL_V.fasta")))

combined_J <- do.call(c, Filter(Negate(is.null),
  lapply(c("IGHJ","IGKJ","IGLJ"), function(nm) J_results[[nm]]$hybrid_gapped)))
write_fasta(combined_J, file.path(gapped_dir, paste0(imgt_file_prefix,"_ALL_J.fasta")))

cat("\n--- Writing ungapped FASTAs ---\n")
for (nm in c("IGHV","IGKV","IGLV"))
  write_fasta(V_results[[nm]]$hybrid_ungapped,
              file.path(ungapped_dir, paste0(file_prefix,"_",nm,".fasta")))
for (nm in c("IGHJ","IGKJ","IGLJ")) {
  hg <- J_results[[nm]]$hybrid_gapped
  if (!is.null(hg))
    write_fasta(ungap(hg),
                file.path(ungapped_dir, paste0(file_prefix,"_",nm,".fasta")))
}
write_fasta(ungap(D_gapped),
            file.path(ungapped_dir, paste0(file_prefix,"_IGHD.fasta")))

# =============================================================================
# SECTION 10 -- IgBLAST auxiliary file (.aux)
#
# IgBLAST aux format (5 tab-separated fields, no header):
#   gene | frame_offset (0-based) | chain_type (JH/JK/JL) |
#   cdr3_stop (1-based nt pos of CDR3 end in J) | extra_bps
# Source: optional_file/mouse_gl.aux header:
#   "gene/allele name, first coding frame start position,
#    chain type, CDR3 stop, extra bps beyond J coding end"
# Internally we store anchor as 0-based; output converts to 1-based cdr3_stop.
#
# Inheritance priority per J gene:
#   1. Reference aux: exact name match
#   2. Reference aux: via header_normalisation_map (normalised_name -> raw)
#   3. Reference aux: via PIgLET annotation table (new_allele -> imgt_allele)
#   4. Conserved motif search in all 3 reading frames
#   5. Default 0 + warning
# =============================================================================
cat("\n--- Building IgBLAST auxiliary file ---\n")

ref_aux_dt <- NULL
for (cand in c(file.path(opt$igdata, "optional_file",
                          paste0(opt$species, "_gl.aux")),
               file.path(opt$igdata, "optional_file",
                          paste0(opt$species, ".aux")))) {
  if (file.exists(cand)) {
    cat(sprintf("  Loading reference aux: %s\n", cand))
    # IgBLAST aux format: 5 whitespace-separated fields, # comment lines, blanks.
    # The file uses MIXED separators (tab after gene name on some lines,
    # spaces on others) so fread(sep=...) always fails.  Parse with readLines.
    #   gene | frame (0-based) | chain_type (JH/JK/JL) | cdr3_stop (0-based) | extra_bps
    aux_raw   <- readLines(cand)
    aux_data  <- aux_raw[!grepl("^\\s*#", aux_raw) & nchar(trimws(aux_raw)) > 0L]
    aux_split <- strsplit(trimws(aux_data), "\\s+")
    aux_rows_ok <- aux_split[sapply(aux_split, length) == 5L]
    if (length(aux_rows_ok) == 0L) {
      message("  [WARN] Could not parse aux file ", cand)
    } else {
      ref_aux_dt <- as.data.table(do.call(rbind, aux_rows_ok))
      setnames(ref_aux_dt, c("gene","frame","chain_type","cdr3_stop","extra_bps"))
      ref_aux_dt[, frame     := as.integer(frame)]
      ref_aux_dt[, cdr3_stop := as.integer(cdr3_stop)]
      ref_aux_dt[, extra_bps := as.integer(extra_bps)]
      cat(sprintf("  Loaded %d entries from reference aux\n", nrow(ref_aux_dt)))
    }
    break
  }
}
if (is.null(ref_aux_dt))
  message("  [WARN] No reference aux found; all anchors from motif search")

build_aux_rows <- function(j_seqs_gapped, chain_type,
                           ref_aux_dt, hmap_dt, annot_dt, ref_seq_map = NULL) {
  if (is.null(j_seqs_gapped) || length(j_seqs_gapped) == 0L)
    return(data.table())
  j_ung <- ungap(j_seqs_gapped)
  rbindlist(lapply(seq_along(j_ung), function(i) {
    gene <- names(j_ung)[i]
    nt   <- as.character(j_ung[[i]])
    d <- .derive_j_aux(gene, nt, chain_type, ref_aux_dt, hmap_dt, annot_dt,
                       ref_seq_map = ref_seq_map)
    data.table(gene=gene, anchor=d$anchor, frame=d$frame, extra_bps=d$extra_bps,
               chain=chain_type, anchor_method=d$method,
               motif=if (is.na(d$motif)) NA_character_ else d$motif,
               motif_anchor=d$motif_anchor, anchor_agrees=d$agrees,
               seqmatch_stop=d$seqmatch_stop,
               seqmatch_gene=if (is.na(d$seqmatch_gene)) NA_character_ else d$seqmatch_gene,
               sequence=nt)
  }))
}

# Reference-J sequence -> curated anchor map for sequence-identity liftover
# validation (built from each locus's reference J sequences + the reference aux).
.asc_ref_jseq_map <- {
  rj <- Biostrings::DNAStringSet()
  for (jl in c("IGHJ","IGKJ","IGLJ")) {
    rg <- J_results[[jl]]$ref_gapped
    if (!is.null(rg) && length(rg) > 0L) rj <- c(rj, rg)
  }
  if (length(rj) > 0L) .build_ref_jseq_anchor_map(rj, ref_aux_dt) else list()
}

aux_rows <- rbindlist(list(
  build_aux_rows(J_results[["IGHJ"]]$hybrid_gapped, "IGH",
                 ref_aux_dt, header_map_dt, J_results[["IGHJ"]]$annot,
                 ref_seq_map = .asc_ref_jseq_map),
  build_aux_rows(J_results[["IGKJ"]]$hybrid_gapped, "IGK",
                 ref_aux_dt, header_map_dt, J_results[["IGKJ"]]$annot,
                 ref_seq_map = .asc_ref_jseq_map),
  build_aux_rows(J_results[["IGLJ"]]$hybrid_gapped, "IGL",
                 ref_aux_dt, header_map_dt, J_results[["IGLJ"]]$annot,
                 ref_seq_map = .asc_ref_jseq_map)
), fill = TRUE)

if (nrow(aux_rows) > 0L) {
  ms <- aux_rows[, .N, by = anchor_method][order(-N)]
  cat(sprintf("  Anchor summary (%d J genes):\n", nrow(aux_rows)))
  for (i in seq_len(nrow(ms)))
    cat(sprintf("    %-22s : %d\n", ms$anchor_method[i], ms$N[i]))

  # ── J-anchor validation report (sequence-identity ground truth) ─────────
  # Ground truth = the curated anchor of a reference J sequence IDENTICAL to the
  # gene's sequence (seqmatch_stop, works for novel names too), else the
  # name-lifted reference anchor. Compares BOTH our final anchor AND the
  # sequence-only motif inference against that ground truth.
  if ("motif_anchor" %in% names(aux_rows)) {
    gt <- copy(aux_rows)
    gt[, gt_stop := ifelse(!is.na(seqmatch_stop), seqmatch_stop,
                    ifelse(grepl("^reference_aux", anchor_method), anchor, NA_integer_))]
    gt[, gt_src := ifelse(!is.na(seqmatch_stop), "seq_identity",
                   ifelse(grepl("^reference_aux", anchor_method), "name_lift", NA_character_))]
    val <- gt[!is.na(gt_stop)]
    n_chk <- nrow(val)
    if (n_chk > 0L) {
      val[, delta_our   := anchor - gt_stop]
      val[, delta_motif := motif_anchor - gt_stop]
      n_seq  <- sum(val$gt_src == "seq_identity", na.rm = TRUE)
      n_ourok <- sum(abs(val$delta_our) <= 1L, na.rm = TRUE)
      n_ourex <- sum(val$delta_our == 0L, na.rm = TRUE)
      mv <- val$delta_motif[!is.na(val$delta_motif)]
      cat(sprintf("  J-anchor validation vs reference (%d genes; %d via exact sequence identity):\n",
                  n_chk, n_seq))
      cat(sprintf("    final anchor  : %d/%d within ±1nt, %d exact\n", n_ourok, n_chk, n_ourex))
      cat(sprintf("    motif-only    : %d/%d within ±1nt, %d exact\n",
                  sum(abs(mv) <= 1L), length(mv), sum(mv == 0L)))
      disc <- val[abs(delta_our) > 0L][order(-abs(delta_our))]
      if (nrow(disc) > 0L) {
        cat(sprintf("    %d final-anchor discrepancy(ies):\n", nrow(disc)))
        for (r in seq_len(min(nrow(disc), 20L)))
          cat(sprintf("      %-24s our=%d vs %s '%s'=%d Δ=%+d motif='%s'\n",
                      disc$gene[r], disc$anchor[r], disc$gt_src[r],
                      if (is.na(disc$seqmatch_gene[r])) disc$gene[r] else disc$seqmatch_gene[r],
                      disc$gt_stop[r], disc$delta_our[r],
                      if (is.na(disc$motif[r])) "?" else disc$motif[r]))
      }
      val_out <- file.path(opt$outdir, "annotations",
                           paste0(opt$prefix, "_jaux_validation.tsv"))
      dir.create(dirname(val_out), recursive = TRUE, showWarnings = FALSE)
      fwrite(val[order(-abs(delta_our)),
             .(gene, chain, anchor_method, ground_truth = gt_src,
               gt_gene = seqmatch_gene, reference_stop = gt_stop,
               our_stop = anchor, motif_stop = motif_anchor,
               delta_our, delta_motif,
               agrees_pm1 = abs(delta_our) <= 1L, motif, sequence)],
             val_out, sep = "\t")
      cat(sprintf("  Wrote J-anchor validation table: %s (%d rows)\n", val_out, n_chk))
    }

    # Novel J genes (motif_search) — motif-derived by definition, no reference
    # to validate against; flag nt-fallback ones as needing manual review.
    novel <- aux_rows[anchor_method == "motif_search"]
    n_novel <- nrow(novel)
    if (n_novel > 0L) {
      n_fallback <- sum(grepl("nt-fallback", novel$motif), na.rm = TRUE)
      cat(sprintf("  %d novel J gene(s) used motif-search anchors (no reference to lift)\n",
                  n_novel))
      if (n_fallback > 0L)
        cat(sprintf("    of which %d used the nucleotide-codon fallback (no clean [WF]G.G motif) — review recommended\n",
                    n_fallback))
    }
  }
}
aux_path <- file.path(opt$outdir, "auxiliary",
                      paste0(file_prefix,"_gl.aux"))
# Write in IgBLAST 5-field format:
# gene | frame_offset | chain_type | cdr3_stop (1-based) | extra_bps
chain_type_map <- c(IGH = "JH", IGK = "JK", IGL = "JL")
aux_out <- aux_rows[, .(
  gene      = gene,
  frame     = frame,
  chain_type = chain_type_map[chain],
  cdr3_stop  = anchor,         # 0-based position (matches reference aux format)
  extra_bps  = extra_bps       # inherited from the resolved reference row (build_aux_rows)
)]
# Write aux file matching the reference format:
#   1. Two comment lines
#   2. Legacy short-name rows (JH1, JK1 etc.) — copied from reference aux
#   3. IMGT allele-level rows for new/hybrid J genes
{
  con <- file(aux_path, open = "wt")
  writeLines(c(
    "#gene/allele name, first coding frame start position, chain type, CDR3 stop, extra bps beyond J coding end.",
    "#All positions are 0-based",
    ""
  ), con)

  # Write legacy short-name rows from reference (JH1, JK1, JL1 etc.)
  if (!is.null(ref_aux_dt) && nrow(ref_aux_dt) > 0L) {
    legacy <- ref_aux_dt[!grepl("\\*", gene) & grepl("^J[HKLA-Z]\\d", gene)]
    if (nrow(legacy) > 0L) {
      for (r in seq_len(nrow(legacy))) {
        writeLines(paste(legacy$gene[r], legacy$frame[r], legacy$chain_type[r],
                         legacy$cdr3_stop[r], legacy$extra_bps[r], sep = "\t"), con)
      }
      writeLines("", con)
    }
  }

  # Write allele-level rows for the hybrid J genes.
  # For strain-tagged names (e.g. IGHJ1*01_C57BL/6), also write a base-name
  # alias row (IGHJ1*01) with the same anchor so igblastn productivity
  # calculation finds the CDR3 anchor regardless of which name variant
  # appears in the alignment output.
  written_genes <- character(0L)
  for (r in seq_len(nrow(aux_out))) {
    nm      <- aux_out$gene[r]
    fr      <- aux_out$frame[r]
    ct      <- aux_out$chain_type[r]
    stop_   <- aux_out$cdr3_stop[r]
    extra_  <- aux_out$extra_bps[r]
    writeLines(paste(nm, fr, ct, stop_, extra_, sep = "\t"), con)
    written_genes <- c(written_genes, nm)
    # Write base-name alias if name has a strain tag (contains _ after the allele)
    # e.g. IGHJ1*01_C57BL/6 -> base = IGHJ1*01
    base_nm <- sub("_[^*_][^*]*$", "", nm)  # strip _anything that follows *xx
    # more precisely: strip the last _<tag> that is NOT part of the allele *xx
    base_nm2 <- sub("(\\*\\d+)_.*$", "\\1", nm)
    if (base_nm2 != nm && !base_nm2 %in% written_genes) {
      writeLines(paste(base_nm2, fr, ct, stop_, extra_, sep = "\t"), con)
      written_genes <- c(written_genes, base_nm2)
    }
  }
  close(con)
}
fwrite(aux_rows, paste0(aux_path, ".diagnostic"), sep = "\t", col.names = TRUE)
cat(sprintf("  Wrote aux file: %s\n", aux_path))


# =============================================================================
# SECTION 10b -- Generate ndm.imgt (V gene FWR/CDR boundary annotation)
# =============================================================================
# Format (13 tab-delimited columns, no header):
#   gene  fwr1s  fwr1e  cdr1s  cdr1e  fwr2s  fwr2e  cdr2s  cdr2e  fwr3s  fwr3e  chain_type  0
#
# CRITICAL: all positions are UNGAPPED nucleotide counts (i.e. count only
# real A/C/G/T characters, not gap dots or dashes).
# igblastn aligns the ungapped query sequence and looks up positions in the
# ndm.imgt by counting non-gap characters — gapped positions would cause
# FWR3 to end ~18-24 nt too late, truncating CDR3 extraction.
#
# Algorithm: scan the IMGT-gapped FASTA to find each boundary in gapped
# coordinates (using IMGT aa position × 3 = gapped nt position), then
# convert each gapped index to an ungapped count via .ung() before writing.
#
# IMGT aa region boundaries (gapped nt search windows):
#   FWR1: aa 1-25  -> gapped 1-75   (VH); 1-78 (VK/VL, aa 1-26)
#   CDR1: aa 27-38 -> gapped 76-114  (variable end — scan for last real nt)
#   FWR2: aa 39-55 -> gapped 115-165
#   CDR2: aa 56-65 -> gapped 166-195 (variable)
#   FWR3: aa 66-104-> gapped 196-312 (Cys104 = end of FWR3)
#
# Typical UNGAPPED positions for a VH with ~18-24 gaps:
#   fwr3_stop ~ 288 (= gapped 312 minus ~24 gap chars)
#
# Chain type codes: VH (heavy), VK (kappa), VL (lambda)
# =============================================================================
cat("\n--- Building ndm.imgt (V gene FWR/CDR annotation) ---\n")

# Return the 1-based gapped-nt position of the FIRST real nt at or after gapped_char_pos
first_nt_from <- function(chars, gapped_char_pos) {
  for (j in seq(gapped_char_pos, length(chars)))
    if (chars[j] != "." && chars[j] != "-") return(j)
  NA_integer_
}

# Return the 1-based gapped-nt position of the LAST real nt at or before gapped_char_pos
last_nt_to <- function(chars, gapped_char_pos) {
  gapped_char_pos <- min(gapped_char_pos, length(chars))
  for (j in seq(gapped_char_pos, 1L, -1L))
    if (chars[j] != "." && chars[j] != "-") return(j)
  NA_integer_
}

build_ndm_rows <- function(v_gapped, chain_type, label) {
  if (is.null(v_gapped) || length(v_gapped) == 0L) return(data.table())

  # One representative allele per gene base
  gene_bases <- sub("[*].*$", "", names(v_gapped))
  v_gapped   <- v_gapped[!duplicated(gene_bases)]

  # IMGT boundary in gapped-nt chars (1-based):
  # FWR1 always starts at char 1.
  # FWR1 ends at char 75 (VH: aa 1-25) or 78 (VK/VL: aa 1-26).
  fwr1_end_char <- if (chain_type == "VH") 75L else 78L
  # CDR1 starts right after FWR1; its end is variable (scan for last real nt).
  # CDR1 region spans up to aa 38 = char 114 (but many sequences are shorter).
  cdr1_max_char <- 114L
  # FWR2: from after CDR1 to aa 55 = char 165
  fwr2_end_char  <- 165L
  # CDR2: from after FWR2 to aa 65 = char 195 (variable)
  cdr2_max_char  <- 195L
  # FWR3: from after CDR2 to aa 104 = char 312
  fwr3_end_char  <- 312L

  # .ung(): convert a gapped character index (1-based, counting gap chars) to
  # the number of non-gap nucleotides up to and including that position.
  # igblastn ndm.imgt requires UNGAPPED positions — count only real nt chars.
  # Gap characters are: . (IMGT) - (dash) = (equals)
  .ung <- function(gapped_idx, chars) {
    if (is.na(gapped_idx) || gapped_idx < 1L) return(gapped_idx)
    sum(!chars[seq_len(min(gapped_idx, length(chars)))] %in% c(".","-","="))
  }

  clamp2 <- function(s_g, e_g, chars) {
    # Convert gapped positions to ungapped, then clamp
    s <- if (is.na(s_g) || s_g < 0L) -1L else as.integer(.ung(s_g, chars))
    e <- if (is.na(e_g) || e_g < 0L) -1L else as.integer(.ung(e_g, chars))
    if (s != -1L && e != -1L && s > e) { s <- -1L; e <- -1L }
    c(s, e)
  }

  rows <- lapply(seq_along(v_gapped), function(i) {
    chars <- strsplit(as.character(v_gapped[[i]]), "")[[1L]]
    n     <- length(chars)

    # Find boundaries in GAPPED coordinates, then convert to UNGAPPED
    # (igblastn ndm.imgt requires ungapped nucleotide positions)
    fwr1_g <- c(first_nt_from(chars, 1L),
                last_nt_to  (chars, min(fwr1_end_char, n)))
    fwr1   <- clamp2(fwr1_g[1L], fwr1_g[2L], chars)

    cdr1_g <- c(first_nt_from(chars, min(fwr1_end_char + 1L, n)),
                last_nt_to  (chars, min(cdr1_max_char, n)))
    cdr1   <- clamp2(cdr1_g[1L], cdr1_g[2L], chars)

    fwr2_search_start <- if (!is.na(cdr1_g[2L]) && cdr1_g[2L] > 0L)
                           cdr1_g[2L] + 1L else cdr1_max_char + 1L
    fwr2_g <- c(first_nt_from(chars, min(fwr2_search_start, n)),
                last_nt_to  (chars, min(fwr2_end_char, n)))
    fwr2   <- clamp2(fwr2_g[1L], fwr2_g[2L], chars)

    cdr2_search_start <- if (!is.na(fwr2_g[2L]) && fwr2_g[2L] > 0L)
                           fwr2_g[2L] + 1L else fwr2_end_char + 1L
    cdr2_g <- c(first_nt_from(chars, min(cdr2_search_start, n)),
                last_nt_to  (chars, min(cdr2_max_char, n)))
    cdr2   <- clamp2(cdr2_g[1L], cdr2_g[2L], chars)

    fwr3_search_start <- if (!is.na(cdr2_g[2L]) && cdr2_g[2L] > 0L)
                           cdr2_g[2L] + 1L else cdr2_max_char + 1L
    fwr3_g <- c(first_nt_from(chars, min(fwr3_search_start, n)),
                last_nt_to  (chars, min(fwr3_end_char, n)))
    fwr3   <- clamp2(fwr3_g[1L], fwr3_g[2L], chars)

    nm       <- names(v_gapped)[i]
    base_nm  <- sub("(\\*\\d+)_.*$", "\\1", nm)   # strip strain tag if present
    row_vals <- data.table(
      gene       = nm,
      fwr1_start = fwr1[1L], fwr1_stop = fwr1[2L],
      cdr1_start = cdr1[1L], cdr1_stop = cdr1[2L],
      fwr2_start = fwr2[1L], fwr2_stop = fwr2[2L],
      cdr2_start = cdr2[1L], cdr2_stop = cdr2[2L],
      fwr3_start = fwr3[1L], fwr3_stop = fwr3[2L],
      chain_type = chain_type,
      trailing   = 0L
    )
    # If strain-tagged, also add a base-name row so igblastn can find
    # FWR/CDR positions regardless of which name variant it reports
    if (base_nm != nm) {
      base_row <- copy(row_vals); base_row$gene <- base_nm
      rbindlist(list(row_vals, base_row))
    } else {
      row_vals
    }
  })
  rbindlist(rows)
}

ndm_rows <- rbindlist(list(
  build_ndm_rows(V_results[["IGHV"]]$hybrid_gapped, "VH", "IGHV"),
  build_ndm_rows(V_results[["IGKV"]]$hybrid_gapped, "VK", "IGKV"),
  build_ndm_rows(V_results[["IGLV"]]$hybrid_gapped, "VL", "IGLV")
), fill = TRUE)

ndm_path <- file.path(opt$outdir, "auxiliary",
                      paste0(file_prefix, ".ndm.imgt"))

# Define int_ndm path unconditionally (used in validation block below)
int_dir <- file.path(opt$outdir, "internal_data", opt$species)
int_ndm <- file.path(int_dir, paste0(file_prefix, ".ndm.imgt"))

if (nrow(ndm_rows) > 0L) {
  # Coerce position columns to integer; clamp any remaining NAs to -1
  pos_cols <- setdiff(names(ndm_rows), c("gene", "chain_type", "trailing"))
  for (col in pos_cols) {
    set(ndm_rows, NULL, col, as.integer(ndm_rows[[col]]))
    set(ndm_rows, which(is.na(ndm_rows[[col]])), col, -1L)
  }

  n_bad <- 0L
  con   <- file(ndm_path, open = "wt")
  for (r in seq_len(nrow(ndm_rows))) {
    pos_vals <- as.integer(unlist(ndm_rows[r, pos_cols, with = FALSE]))
    pos_vals[is.na(pos_vals)] <- -1L
    if (length(pos_vals) != 10L) { n_bad <- n_bad + 1L; next }
    writeLines(paste(c(ndm_rows$gene[r], pos_vals,
                       ndm_rows$chain_type[r], ndm_rows$trailing[r]),
                     collapse = "\t"), con)
  }
  close(con)
  if (n_bad > 0L)
    message(sprintf("  [NDM WARN] %d rows skipped (wrong field count)", n_bad))

  fwrite(ndm_rows,
         paste0(tools::file_path_sans_ext(ndm_path), "_annotated.tsv"),
         sep = "\t", col.names = TRUE)
  cat(sprintf("  Wrote ndm.imgt : %s  (%d V genes)\n", ndm_path, nrow(ndm_rows)))

  if (dir.exists(int_dir)) {
    file.copy(ndm_path, int_ndm, overwrite = TRUE)
    cat(sprintf("  Copied to      : %s\n", int_ndm))
  }

  # ---- Comprehensive validation ----
  # Re-read and check every row: 13 fields, 10 integer positions, valid chain type,
  # no inversions, no empty fields.
  ndm_raw <- readLines(ndm_path)
  ndm_raw <- ndm_raw[nchar(trimws(ndm_raw)) > 0L]

  valid_chains <- c("VH","VK","VL","VA","VB","VD","VG")

  diagnose_row <- function(row) {
    parts <- strsplit(row, "\t")[[1L]]
    if (length(parts) != 13L) return(sprintf("wrong_fields(%d)", length(parts)))
    issues <- character(0L)
    if (!nchar(trimws(parts[1L])))
      issues <- c(issues, "empty_gene")
    vals <- suppressWarnings(as.integer(parts[2:11]))
    na_idx <- which(is.na(vals))
    if (length(na_idx))
      issues <- c(issues, sprintf("non_int_cols(%s)", paste(na_idx+1L, collapse=",")))
    bad_rng <- which(!is.na(vals) & vals != -1L & (vals < 1L | vals > 400L))
    if (length(bad_rng))
      issues <- c(issues, sprintf("out_of_range_cols(%s)", paste(bad_rng+1L, collapse=",")))
    region_pairs <- list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10))
    region_names <- c("FWR1","CDR1","FWR2","CDR2","FWR3")
    for (k in seq_along(region_pairs)) {
      s <- vals[region_pairs[[k]][1]]; e <- vals[region_pairs[[k]][2]]
      if (!is.na(s) && !is.na(e) && s != -1L && e != -1L && s > e)
        issues <- c(issues, sprintf("%s_inverted(%d>%d)", region_names[k], s, e))
    }
    if (!trimws(parts[12L]) %in% valid_chains)
      issues <- c(issues, sprintf("bad_chain(%s)", parts[12L]))
    if (length(issues)) paste(issues, collapse="; ") else ""
  }

  diag     <- sapply(ndm_raw, diagnose_row, USE.NAMES = FALSE)
  bad_mask <- nchar(diag) > 0L
  n_bad2   <- sum(bad_mask)

  if (n_bad2 > 0L) {
    message(sprintf("  [NDM] %d/%d rows invalid -- repairing:", n_bad2, length(ndm_raw)))
    for (r in head(which(bad_mask), 10L))
      message(sprintf("    row %3d [%s]: %s", r, diag[r], substr(ndm_raw[r], 1L, 80L)))
    ndm_clean <- ndm_raw[!bad_mask]
    writeLines(ndm_clean, ndm_path)
    if (dir.exists(int_dir)) writeLines(ndm_clean, int_ndm)
    message(sprintf("  [NDM] Rewrote: %d clean rows", length(ndm_clean)))
  } else {
    writeLines(ndm_raw, ndm_path)
    if (dir.exists(int_dir)) writeLines(ndm_raw, int_ndm)
  }

  # Column stats: range of each position column
  clean_rows <- ndm_raw[!bad_mask]
  if (length(clean_rows) > 0L) {
    vals_mat  <- do.call(rbind, lapply(clean_rows, function(r)
      suppressWarnings(as.integer(strsplit(r, "\t")[[1L]][2:11]))))
    col_names <- c("fwr1s","fwr1e","cdr1s","cdr1e","fwr2s","fwr2e",
                   "cdr2s","cdr2e","fwr3s","fwr3e")
    ranges    <- apply(vals_mat, 2L, function(v) range(v[v != -1L], na.rm = TRUE))
    cat(sprintf("  Validated: %d/%d rows clean\n",
                length(clean_rows), length(ndm_raw)))
    cat("  Position ranges (excl. -1 sentinels):\n")
    for (k in seq_len(ncol(ranges)))
      cat(sprintf("    %-8s: %d - %d\n", col_names[k], ranges[1L,k], ranges[2L,k]))
  }
} else {
  cat("  [WARN] No V genes available for ndm.imgt\n")
}

# SECTION 11 -- Debug R script
# (Written early via write_debug_script(); see call after Section 6)
# =============================================================================
# (already written)

# =============================================================================
# SECTION 12 -- Manifest
# =============================================================================
manifest <- data.table(
  file_type = c(
    "header_normalisation_map", "allele_cluster_annotation",
    "full_provenance", "debug_script",
    "gapped_ALL_V",
    "gapped_IGHV","gapped_IGKV","gapped_IGLV",
    "gapped_ALL_J",
    "gapped_IGHJ","gapped_IGKJ","gapped_IGLJ","gapped_IGHD",
    "ungapped_IGHV","ungapped_IGKV","ungapped_IGLV",
    "ungapped_IGHJ","ungapped_IGKJ","ungapped_IGLJ","ungapped_IGHD",
    "aux_file","aux_diagnostic"
  ),
  path = c(
    header_map_out, annot_out, prov_out, debug_script_path,
    file.path(gapped_dir, paste0(imgt_file_prefix,"_ALL_V.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGHV.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGKV.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGLV.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_ALL_J.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGHJ.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGKJ.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGLJ.fasta")),
    file.path(gapped_dir, paste0(imgt_file_prefix,"_IGHD.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGHV.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGKV.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGLV.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGHJ.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGKJ.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGLJ.fasta")),
    file.path(ungapped_dir, paste0(file_prefix,"_IGHD.fasta")),
    aux_path, paste0(aux_path,".diagnostic")
  )
)
manifest$exists <- file.exists(manifest$path)
fwrite(manifest, file.path(opt$outdir, paste0(opt$prefix,"_manifest.tsv")), sep="\t")

cat("\n=== R annotation step complete ===\n")
cat(sprintf("Manifest    : %s\n", file.path(opt$outdir, paste0(opt$prefix,"_manifest.tsv"))))
cat(sprintf("Debug script: %s\n", debug_script_path))
cat(sprintf("Provenance  : %s\n", prov_out))