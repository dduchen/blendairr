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
  make_option("--organism",                 type = "character",    default = NULL,
              help = "Organism name used as output file prefix (default: <prefix>_<species>)")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("custom_dir","ref_dir","outdir","igdata"))
  if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))

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
cat(sprintf("Naming mode: %s\n", if (isTRUE(opt$use_asc)) "ASC (PIgLET cluster names)" else "IMGT (reference-based)"))

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
  raw_content <- as.character(RemoveGaps(seqs, removeGaps = "all"))
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
  # RemoveGaps requires XStringSet; for character vectors strip gap chars directly.
  if (is(seqs, "XStringSet"))
    return(RemoveGaps(seqs, removeGaps = "all"))
  if (is.character(seqs))
    return(gsub("[.-]", "", seqs))   # [.-] in a char class: literal dot and hyphen
  RemoveGaps(seqs, removeGaps = "all")
}

normalise_name <- function(nm) {
  nm <- gsub("-G-",   "-",  nm)
  nm <- gsub("-G\\*", "\\*", nm)
  nm <- sub("(\\*\\d+)[FP]$", "\\1", nm)
  nm
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
  tbl[, new_allele := normalise_name(new_allele)]
  if (!"imgt_allele" %in% names(tbl) && "allele" %in% names(tbl))
    setnames(tbl, "allele", "imgt_allele")

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
  custom_content <- as.character(RemoveGaps(custom_seqs, removeGaps = "all"))
  ref_content    <- as.character(RemoveGaps(ref_name_only, removeGaps = "all"))
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
  custom_seqs_chr <- as.character(RemoveGaps(custom_gapped, removeGaps = "all"))
  ref_seqs_chr    <- as.character(RemoveGaps(ref_gapped,    removeGaps = "all"))

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
  # Reverse map: cleaned name -> original IMGT name (preserves full subtype)
  # e.g. "IGHV1-2*01" -> "IGHV1-2*01" (unchanged if clean)
  # or   "IGHV1*02"   -> "IGHV1_dup2" -> original was "IGHV1_dup2" (bare family)
  # Used when recovering gene_base from dist_mat closest-ref column names.
  rev_ref_name_map <- setNames(ref_nms, ref_nms_clean)

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

  ref_ungapped_chr <- as.character(RemoveGaps(ref_gapped, removeGaps = "all"))

  for (dc in missing_custom) {
    dc_seq    <- as.character(RemoveGaps(custom_gapped[dc], removeGaps = "all"))

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
        extra[, imgt_allele   := dc]
        extra[, new_allele    := imgt_ref_name]
        extra[, piglet_cluster := if ("new_allele" %in% names(ref_row))
                                    ref_row$new_allele[1L] else NA_character_]
        tbl_custom <- rbind(tbl_custom, extra, fill = TRUE)
        message(sprintf("  [INFO] %s: '%s' identical to ref '%s' -> adopts IMGT name",
                        label, dc, imgt_ref_name))
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
      ref_ug_clean <- as.character(RemoveGaps(ref_for_joint, removeGaps = "all"))
      dists        <- as.integer(adist(dc_seq, ref_ug_clean))
      best_cleaned <- names(ref_for_joint)[which.min(dists)]
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
    cluster_ref_rep <- tbl_ref[
      order(-n_dashes, allele_num_sort),
      .(ref_rep = imgt_allele[1L]),
      by = cluster_id
    ]
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
        ref_cols  <- names(row)[names(row) %in% ref_names_clean]
        if (length(ref_cols) > 0L) {
          closest     <- ref_cols[which.min(row[ref_cols])]
          # Recover original name for gene_base (cleaned may have lost subtype)
          orig_closest <- if (exists("rev_ref_name_map") && closest %in% names(rev_ref_name_map))
                            rev_ref_name_map[[closest]] else closest
          best_gene <- sub("[*].*$", "", orig_closest)
        }
      }
      # Fallback: ungapped edit distance if dist_mat unavailable
      if (is.null(best_gene)) {
        seq_i     <- as.character(RemoveGaps(
                       custom_gapped[cust_nm], removeGaps = "all"))
        ref_ug    <- as.character(RemoveGaps(ref_for_joint, removeGaps = "all"))
        dists     <- as.integer(adist(seq_i, ref_ug))
        best_cleaned  <- names(ref_for_joint)[which.min(dists)]
        orig_closest  <- if (exists("rev_ref_name_map") && best_cleaned %in% names(rev_ref_name_map))
                           rev_ref_name_map[[best_cleaned]] else best_cleaned
        best_gene <- sub("[*].*$", "", orig_closest)
      }
      tbl_cust$new_allele_final[i] <-
        paste0(best_gene, "*", .next_allele(best_gene))
    }
  }
  tbl_cust[, new_allele := new_allele_final]
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
    "  seqs_ungapped <- DECIPHER::RemoveGaps(seqs_gapped, removeGaps = \"all\")",
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
    return(list(custom_gapped = NULL, ref_gapped = ref,
                hybrid_gapped = ref, hybrid_ungapped = ungap(ref),
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
  hg <- merge_with_priority(res$renamed_custom, ref)
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
prov_out <- file.path(opt$outdir, "annotations",
                      paste0(opt$prefix, "_full_provenance.tsv"))
fwrite(provenance, prov_out, sep = "\t")
cat(sprintf("  Wrote full provenance    : %s\n", prov_out))

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
        grp_ung <- RemoveGaps(grp_seqs, removeGaps = "all")
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
  no_star <- !grepl("[*]", new_names)
  if (any(no_star)) {
    message(sprintf("  [POST-NORM] %s: %d names still missing allele number -> *01",
                    label, sum(no_star)))
    new_names[no_star] <- paste0(new_names[no_star], "*01")
  }

  # ---- Step 4: re-number colliding gene bases sequentially ----
  gene_allele_count <- list()
  for (i in seq_along(new_names)) {
    nm        <- new_names[i]
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
  ug_content <- as.character(RemoveGaps(seqs, removeGaps = "all"))
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

# Write the final name map: original raw header -> final sequence name in FASTA.
# This joins: raw_header -> normalised_name -> post-norm name (after subfamily
# inference and sequential numbering) -- one row per sequence in the output FASTA.
if (length(all_post_norm_maps) > 0L) {
  post_norm_dt <- rbindlist(all_post_norm_maps, fill = TRUE)
  # Join with header_map_dt to get raw_header -> normalised_name -> final_name
  final_name_map <- merge(
    header_map_dt[, .(raw_header, normalised_name, parse_flag, locus, source)],
    post_norm_dt[, .(normalised_name = before, final_fasta_name = after, locus)],
    by = c("normalised_name", "locus"),
    all.x = TRUE
  )
  # Sequences unchanged by post_normalise_seqs have final_fasta_name == normalised_name
  final_name_map[is.na(final_fasta_name), final_fasta_name := normalised_name]
  final_map_out <- file.path(opt$outdir, "annotations",
                              paste0(opt$prefix, "_final_name_map.tsv"))
  fwrite(final_name_map[order(locus, source, raw_header)],
         final_map_out, sep = "\t")
  cat(sprintf("  Wrote final name map     : %s  (%d sequences)\n",
              final_map_out, nrow(final_name_map)))
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

find_j_anchor <- function(nt_seq, chain_type) {
  motifs <- list(
    IGH = c("WGQG","WGPG","WGRG","FGQG","FGAG","FGSG","WGAG","FGTG"),
    IGK = c("FGQG","FGPG","FGRG","WGQG","FGSG","FGTG","WGPG"),
    IGL = c("FGGG","FGSG","WGSG","FGAG","FGTG","WGGG")
  )
  chain_motifs <- motifs[[chain_type]] %||%
    unique(unlist(motifs, use.names = FALSE))
  seq_obj  <- DNAString(nt_seq)
  best_pos <- NA_integer_
  for (frame in 0L:2L) {
    sublen <- nchar(nt_seq) - frame
    sublen <- sublen - (sublen %% 3L)
    if (sublen < 3L) next
    aa <- as.character(translate(subseq(seq_obj, frame+1L, frame+sublen)))
    for (motif in chain_motifs) {
      pa <- regexpr(motif, aa, fixed = TRUE)
      if (pa > 0L) { best_pos <- (pa-1L)*3L+frame; break }
    }
    if (!is.na(best_pos)) break
  }
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
  # Try all candidate names derived from gene_name
  for (cand in .aux_candidates(gene_name)) {
    hit <- ref_aux_dt[gene == cand]
    if (nrow(hit) > 0L) return(hit[1L])
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

build_aux_rows <- function(j_seqs_gapped, chain_type,
                           ref_aux_dt, hmap_dt, annot_dt) {
  if (is.null(j_seqs_gapped) || length(j_seqs_gapped) == 0L)
    return(data.table())
  j_ung <- ungap(j_seqs_gapped)
  rbindlist(lapply(seq_along(j_ung), function(i) {
    gene   <- names(j_ung)[i]
    nt     <- as.character(j_ung[[i]])
    anchor <- NA_integer_; frame <- 0L; method <- "none"
    rh <- lookup_ref_anchor(gene, ref_aux_dt, hmap_dt, annot_dt)
    if (!is.null(rh)) {
      # cdr3_stop in aux is 0-based (per header comment "All positions are 0-based")
      anchor <- as.integer(rh$cdr3_stop)
      frame  <- as.integer(rh$frame)
      method <- "reference_aux"
    }
    if (is.na(anchor)) {
      anchor <- find_j_anchor(nt, chain_type)
      if (!is.na(anchor)) { frame <- anchor %% 3L; method <- "motif_search" }
    }
    if (is.na(anchor)) {
      message(sprintf("    [WARN] No anchor for %s (%s); defaulting 0", gene, chain_type))
      anchor <- 0L; method <- "default"
    }
    data.table(gene=gene, anchor=anchor, frame=frame,
               chain=chain_type, anchor_method=method)
  }))
}

aux_rows <- rbindlist(list(
  build_aux_rows(J_results[["IGHJ"]]$hybrid_gapped, "IGH",
                 ref_aux_dt, header_map_dt, J_results[["IGHJ"]]$annot),
  build_aux_rows(J_results[["IGKJ"]]$hybrid_gapped, "IGK",
                 ref_aux_dt, header_map_dt, J_results[["IGKJ"]]$annot),
  build_aux_rows(J_results[["IGLJ"]]$hybrid_gapped, "IGL",
                 ref_aux_dt, header_map_dt, J_results[["IGLJ"]]$annot)
), fill = TRUE)

if (nrow(aux_rows) > 0L) {
  ms <- aux_rows[, .N, by = anchor_method][order(-N)]
  cat(sprintf("  Anchor summary (%d J genes):\n", nrow(aux_rows)))
  for (i in seq_len(nrow(ms)))
    cat(sprintf("    %-22s : %d\n", ms$anchor_method[i], ms$N[i]))
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
  extra_bps  = 0L             # conservative default; inherited from ref when available
)]
# Where we inherited from reference, use its extra_bps directly
if (!is.null(ref_aux_dt) && nrow(ref_aux_dt) > 0L) {
  ref_extras <- ref_aux_dt[, .(gene, extra_bps)]
  aux_out <- merge(aux_out, ref_extras, by = "gene", all.x = TRUE,
                   suffixes = c("","_ref"))
  if ("extra_bps_ref" %in% names(aux_out))
    aux_out[!is.na(extra_bps_ref), extra_bps := extra_bps_ref]
  aux_out[, extra_bps_ref := NULL]
}
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

  # Write allele-level rows for the hybrid J genes
  for (r in seq_len(nrow(aux_out))) {
    writeLines(paste(aux_out$gene[r], aux_out$frame[r], aux_out$chain_type[r],
                     aux_out$cdr3_stop[r], aux_out$extra_bps[r], sep = "\t"), con)
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
# All positions are 1-based GAPPED nucleotide positions (each character in the
# IMGT-gapped sequence = 1 position, whether a real nucleotide or a gap dot).
# Positions are found by scanning the gapped sequence directly for the last/first
# real nucleotide within each IMGT region, in gapped-nt coordinates.
#
# IMGT aa boundaries (1-based aa positions, each aa = 3 gapped nt chars):
#   FWR1: aa 1-25   (gapped nt 1-75)   [VH uses 25 aa; VK/VL use 26 aa = 78 nt]
#   CDR1: variable length (starts right after FWR1)
#   FWR2: fixed end at aa 55 (gapped nt 165)
#   CDR2: variable
#   FWR3: ends at aa 104 (gapped nt 312)
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

  clamp2 <- function(s, e) {
    s <- if (is.na(s)) -1L else as.integer(s)
    e <- if (is.na(e)) -1L else as.integer(e)
    if (s != -1L && e != -1L && s > e) { s <- -1L; e <- -1L }
    c(s, e)
  }

  rows <- lapply(seq_along(v_gapped), function(i) {
    chars <- strsplit(as.character(v_gapped[[i]]), "")[[1L]]
    n     <- length(chars)

    # FWR1: chars 1 to fwr1_end_char
    fwr1 <- clamp2(first_nt_from(chars, 1L),
                   last_nt_to  (chars, min(fwr1_end_char, n)))

    # CDR1: chars fwr1_end_char+1 to cdr1_max_char
    cdr1 <- clamp2(first_nt_from(chars, min(fwr1_end_char + 1L, n)),
                   last_nt_to  (chars, min(cdr1_max_char, n)))

    # FWR2: first real nt after CDR1 region, to fwr2_end_char
    fwr2_search_start <- if (cdr1[2L] > 0L) cdr1[2L] + 1L else cdr1_max_char + 1L
    fwr2 <- clamp2(first_nt_from(chars, min(fwr2_search_start, n)),
                   last_nt_to  (chars, min(fwr2_end_char, n)))

    # CDR2: first real nt after FWR2, to cdr2_max_char
    cdr2_search_start <- if (fwr2[2L] > 0L) fwr2[2L] + 1L else fwr2_end_char + 1L
    cdr2 <- clamp2(first_nt_from(chars, min(cdr2_search_start, n)),
                   last_nt_to  (chars, min(cdr2_max_char, n)))

    # FWR3: first real nt after CDR2, to fwr3_end_char
    fwr3_search_start <- if (cdr2[2L] > 0L) cdr2[2L] + 1L else cdr2_max_char + 1L
    fwr3 <- clamp2(first_nt_from(chars, min(fwr3_search_start, n)),
                   last_nt_to  (chars, min(fwr3_end_char, n)))

    data.table(
      gene       = names(v_gapped)[i],
      fwr1_start = fwr1[1L], fwr1_stop = fwr1[2L],
      cdr1_start = cdr1[1L], cdr1_stop = cdr1[2L],
      fwr2_start = fwr2[1L], fwr2_stop = fwr2[2L],
      cdr2_start = cdr2[1L], cdr2_stop = cdr2[2L],
      fwr3_start = fwr3[1L], fwr3_stop = fwr3[2L],
      chain_type = chain_type,
      trailing   = 0L
    )
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