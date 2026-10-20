#!/usr/bin/env Rscript
# =============================================================================
# plot_asc_circos_interactive.R
#
# Interactive circular (circos-style) visualization of PIgLET ASC clustering and
# allele assignments, built from a blendAIRR cluster-annotation table. Rendered
# as concentric arc tracks with a radial clustering dendrogram in the centre, in
# the style of PIgLET's circular dendrogram output:
#
#   Outer ring  : FAMILY clusters (family_cluster) as coloured arcs
#   Middle ring : GENE clusters   (cluster_id)      as arcs within each family
#   Leaf ticks  : ALLELES as radial ticks with markers
#   Centre tree : radial dendrogram of the ASC hierarchy
#                 (allele -> gene cluster -> family cluster), with allele order
#                 inside each gene cluster set by a k-mer similarity tree
#
# Hovering any allele leaf shows its ORIGINAL annotation, the ASC name, and the
# mapped TR/IG NC colon name. Hovering a ring shows the cluster plus its TRIG-NC
# family / gene. The TR/IG NC names are read from the pipeline's
# *_final_name_map.tsv (trignc_name column), which carries the colon form for
# every clustered allele; the annotation table's new_allele holds the ASC name.
#
# INPUT  : *_allele_cluster_annotation.tsv (imgt_allele, new_allele,
#          piglet_cluster, cluster_id, family_cluster, source, locus, sequence)
#          plus, when present, the sibling *_final_name_map.tsv for TRIG-NC names
# OUTPUT : one interactive .html per locus (plotly).
#
# USAGE
#   Rscript plot_asc_circos_interactive.R \
#     --annot human_trig/annotations/hybrid_asc_allele_cluster_annotation.tsv \
#     --locus IGHV \
#     --out   human_trig/annotations/hybrid_asc_IGHV_circos.html
#
#   # all loci -> one file per locus, written next to --annot
#   Rscript plot_asc_circos_interactive.R --annot <file>
#
#   # point explicitly at the name map if it is not the sibling default
#   Rscript plot_asc_circos_interactive.R --annot <file> \
#     --name-map <dir>/hybrid_asc_final_name_map.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(plotly)
  library(htmlwidgets)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

opt <- parse_args(OptionParser(option_list = list(
  make_option("--annot",  type = "character", help = "cluster-annotation TSV"),
  make_option("--locus",  type = "character", default = NULL,
              help = "single locus (e.g. IGHV); omit to plot each locus"),
  make_option("--out",    type = "character", default = NULL,
              help = "output .html (single-locus); default derives from --annot"),
  make_option("--name-map", type = "character", default = NULL,
              help = "final_name_map.tsv with trignc_name; default is the sibling *_final_name_map.tsv"),
  make_option("--title",  type = "character", default = NULL),
  make_option("--gap-deg", type = "double", default = 2,
              help = "angular gap (degrees) between family wedges [2]"),
  make_option("--hide-dendrogram", action = "store_true", default = FALSE,
              help = "omit the central radial dendrogram"),
  make_option("--label-min", type = "integer", default = 0,
              help = "draw allele labels only when a gene cluster has <= this many alleles (0 = none)")
)))
if (is.null(opt$annot)) stop("--annot is required")

annot <- fread(opt$annot, sep = "\t", na.strings = c("", "NA"))
req <- c("imgt_allele", "cluster_id", "family_cluster", "locus")
miss <- setdiff(req, names(annot))
if (length(miss)) stop("annotation missing required column(s): ",
                       paste(miss, collapse = ", "))
# Ensure the string columns exist AND are character. fread types an all-empty
# column as logical (e.g. piglet_cluster when nothing was clustered), which then
# breaks fifelse() branches expecting character; coerce them here.
for (c0 in c("imgt_allele", "new_allele", "piglet_cluster", "cluster_id",
             "family_cluster", "cluster_status", "source", "sequence", "locus")) {
  if (!c0 %in% names(annot)) annot[[c0]] <- NA_character_
  else if (!is.character(annot[[c0]])) annot[[c0]] <- as.character(annot[[c0]])
}

# ── TR/IG NC name map ────────────────────────────────────────────────────────
# The colon-form names (IGHV:01:007:001) live in the pipeline's final_name_map,
# keyed by original_id, not in the annotation table. Load it (sibling by default)
# and keep only genuine colon-form names so unchanged custom D/J ids are not
# mistaken for TRIG-NC assignments.
if (is.null(opt[["name-map"]])) {
  cand <- sub("_allele_cluster_annotation\\.tsv$", "_final_name_map.tsv", opt$annot)
  if (file.exists(cand) && cand != opt$annot) opt[["name-map"]] <- cand
}
# Reason an allele is not in an ASC cluster, read from the IMGT pipe header
# (field 4 = F/(F)/P/ORF, plus a "partial in 3'/5'" note). Used only as a
# fallback when the annotation table has no cluster_status column.
imgt_reason <- function(h) {
  if (length(h) != 1L || is.na(h) || !grepl("\\|", h)) return(NA_character_)
  f  <- strsplit(h, "\\|", perl = TRUE)[[1]]
  fx <- if (length(f) >= 4L) toupper(gsub("[^A-Za-z]", "", f[4])) else ""
  parts <- character(0)
  if (grepl("partial in 3'", h, fixed = TRUE)) parts <- c(parts, "3' partial")
  if (grepl("partial in 5'", h, fixed = TRUE)) parts <- c(parts, "5' partial")
  if (fx == "P")   parts <- c(parts, "pseudogene")
  if (fx == "ORF") parts <- c(parts, "ORF")
  if (!length(parts)) NA_character_ else paste(parts, collapse = ", ")
}

trig_lk <- NULL
reason_lk <- character(0)
if (!is.null(opt[["name-map"]]) && file.exists(opt[["name-map"]])) {
  nm <- fread(opt[["name-map"]], sep = "\t", na.strings = c("", "NA"))
  if ("trignc_name" %in% names(nm)) {
    nmc <- nm[!is.na(trignc_name) & grepl(":", trignc_name, fixed = TRUE)]
    # Key the lookup by the CLEAN allele name (normalised_name) first, since that
    # is what the annotation's imgt_allele holds. original_id is the raw header
    # for reference-derived rows (e.g. "M77300|IGHV3-30*09|..."), so keying only
    # by it would miss every reference-only allele. Keep original_id as a
    # secondary key for safety.
    keys <- c(if ("normalised_name" %in% names(nmc)) nmc$normalised_name else character(0),
              if ("original_id"     %in% names(nmc)) nmc$original_id     else character(0))
    vals <- c(if ("normalised_name" %in% names(nmc)) nmc$trignc_name else character(0),
              if ("original_id"     %in% names(nmc)) nmc$trignc_name else character(0))
    ok <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
    trig_lk <- setNames(vals[ok], keys[ok])
    cat(sprintf("Loaded %d TR/IG NC name keys from %s\n",
                length(trig_lk), opt[["name-map"]]))
  }
  if (all(c("raw_header", "normalised_name") %in% names(nm))) {
    rf <- nm[source == "reference" & !is.na(raw_header)]
    if (nrow(rf)) {
      rr  <- vapply(rf$raw_header, imgt_reason, character(1), USE.NAMES = FALSE)
      keys <- c(rf$normalised_name,
                if ("trignc_name" %in% names(rf)) rf$trignc_name else character(0))
      vals <- c(rr, rr)
      ok   <- !is.na(vals) & nzchar(keys) & !duplicated(keys)
      reason_lk <- setNames(vals[ok], keys[ok])
    }
  }
}

# ── Geometry helpers ─────────────────────────────────────────────────────────
# Annular sector (arc band) a0..a1 (radians), radii r0..r1, as a closed polygon.
arc_band <- function(a0, a1, r0, r1, n = 64) {
  th_out <- seq(a0, a1, length.out = n)
  th_in  <- rev(th_out)
  x <- c(r1 * cos(th_out), r0 * cos(th_in))
  y <- c(r1 * sin(th_out), r0 * sin(th_in))
  list(x = c(x, x[1]), y = c(y, y[1]))
}

pal <- function(n) {
  base <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948",
            "#B07AA1","#FF9DA7","#9C755F","#BAB0AC","#1B9E77","#D95F02",
            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
  if (n <= length(base)) base[seq_len(n)]
  else grDevices::colorRampPalette(base)(n)
}

src_symbol_colour <- function(s) {
  s <- ifelse(is.na(s), "unknown", s)
  fifelse(grepl("custom",      s), "#D55E00",
  fifelse(grepl("reference",   s), "#0072B2",
  fifelse(grepl("unclustered", s), "#999999", "#111111")))
}

# k-mer similarity ordering (base R, no alignment / DECIPHER). Returns a position
# vector so that sorting by it places similar sequences next to each other.
kmer_rank <- function(seqs, k = 4L) {
  n <- length(seqs)
  if (n <= 2L) return(seq_len(n))
  s <- toupper(gsub("[^ACGT]", "", ifelse(is.na(seqs), "", seqs)))
  if (any(nchar(s) < k)) return(seq_len(n))
  kmers <- lapply(s, function(x)
    substring(x, seq_len(nchar(x) - k + 1L), seq(k, nchar(x))))
  lev <- unique(unlist(kmers))
  m <- t(vapply(kmers, function(kk)
    as.integer(table(factor(kk, levels = lev))), integer(length(lev))))
  ord <- tryCatch(hclust(dist(m), method = "average")$order,
                  error = function(e) seq_len(n))
  pos <- integer(n); pos[ord] <- seq_len(n); pos
}

# Append a transparent, still-hoverable marker grid across an arc band so the
# whole wedge reports its text on hover (filled-polygon fill-hover is unreliable
# across browsers and tends to show the raw fill colour instead).
hover_ring <- function(p, a0, a1, radii, text) {
  span <- a1 - a0
  nang <- max(6L, as.integer(round(span / (1.5 * pi / 180))))
  ang  <- seq(a0, a1, length.out = nang)
  g    <- expand.grid(a = ang, r = radii)
  add_trace(p, x = g$r * cos(g$a), y = g$r * sin(g$a),
            type = "scatter", mode = "markers",
            marker = list(size = 12, color = "rgba(0,0,0,0)"),
            hovertemplate = paste0(text, "<extra></extra>"),
            showlegend = FALSE)
}

# ── Build one locus plot ─────────────────────────────────────────────────────
build_locus_plot <- function(dt, locus_val, title = NULL) {
  d <- dt[locus == locus_val]
  if (nrow(d) == 0L) return(NULL)
  d <- copy(d)
  d[, family_cluster := as.character(family_cluster)]
  d[, cluster_id     := as.character(cluster_id)]
  if (!"cluster_status" %in% names(d)) d[, cluster_status := NA_character_]

  # One leaf per allele. The annotation carries a row per (allele, source), so an
  # allele present in both custom and reference appears twice; collapse those
  # exact duplicates to a single leaf, preferring the custom-tagged row for its
  # colour. Rows that genuinely differ (different cluster or ASC name) are kept.
  key_cols <- intersect(c("imgt_allele", "new_allele", "cluster_id",
                          "family_cluster", "cluster_status"), names(d))
  if (length(key_cols)) {
    n0 <- nrow(d)
    d[, .pref := fifelse(!is.na(source) & grepl("custom", source), 0L, 1L)]
    setorder(d, .pref)
    d <- unique(d, by = key_cols)
    d[, .pref := NULL]
    if (nrow(d) < n0)
      message(sprintf("  [%s] collapsed %d duplicate source row(s) to %d leaves",
                      locus_val, n0 - nrow(d), nrow(d)))
  }

  if (d[!is.na(family_cluster) & !is.na(cluster_id), .N] == 0L) return(NULL)

  # Which alleles PIgLET did not place in an ASC cluster, detected BEFORE the
  # empty cluster columns are filled. Prefer the builder's cluster_status; fall
  # back to empty cluster columns for older annotation tables.
  d[, is_unclust := if (any(!is.na(cluster_status)))
      (!is.na(cluster_status) & cluster_status != "clustered")
    else (is.na(family_cluster) | is.na(cluster_id))]

  # Names: ASC from piglet_cluster (or non-colon new_allele), TR/IG NC colon
  # name from the name map, falling back to a colon-form new_allele/imgt_allele.
  d[, asc_name  := fifelse(!is.na(piglet_cluster) & nzchar(piglet_cluster), piglet_cluster,
                   fifelse(!is.na(new_allele) & !grepl(":", new_allele, fixed = TRUE) &
                           nzchar(new_allele), new_allele, NA_character_))]
  d[, trig_name := if (!is.null(trig_lk)) unname(trig_lk[imgt_allele]) else NA_character_]
  d[is.na(trig_name) & !is.na(new_allele)  & grepl(":", new_allele,  fixed = TRUE),
    trig_name := new_allele]
  d[is.na(trig_name) & !is.na(imgt_allele) & grepl(":", imgt_allele, fixed = TRUE),
    trig_name := imgt_allele]
  d[, trig_fam  := fifelse(is.na(trig_name), NA_character_,
                           sub("^([^:]+):([^:]+):.*$", "\\1:\\2", trig_name))]
  d[, trig_gene := fifelse(is.na(trig_name), NA_character_,
                           sub("^([^:]+):([^:]+):([^:]+):.*$", "\\1:\\2:\\3", trig_name))]

  # In --trig_nc mode the authoritative grouping is the TRIG-NC family and gene,
  # not PIgLET's internal family_cluster. PIgLET's 75% family threshold can
  # collapse many closely-related families into a single ASC cluster (e.g. all
  # human IGHV land in family_cluster 1), which would draw the whole locus as one
  # family. Repoint the rings to the colon family/gene wherever a TRIG-NC name
  # exists; the ASC cluster numbers remain the fallback for anything without one.
  trig_mode <- any(!is.na(d$trig_name))
  if (trig_mode) {
    d[!is.na(trig_name), family_cluster := sub("^[^:]+:([^:]+):.*$", "\\1", trig_name)]
    d[!is.na(trig_name), cluster_id     := trig_gene]
  }

  # Reason for the unclustered ones: cluster_status (strip the prefix) if the
  # builder wrote it, else the IMGT-header reason from the name map, else a
  # generic note.
  d[, reason := NA_character_]
  d[is_unclust & !is.na(cluster_status),
    reason := sub("^not ASC-clustered:?\\s*", "", cluster_status)]
  if (length(reason_lk))
    d[is_unclust & (is.na(reason) | !nzchar(reason)),
      reason := { r <- unname(reason_lk[imgt_allele])
                  fifelse(!is.na(r), r, unname(reason_lk[new_allele])) }]
  d[is_unclust & (is.na(reason) | !nzchar(reason)),
    reason := "not placed in an ASC cluster"]

  # Place each unclustered allele with its real gene/family instead of a catch-all.
  #  1) by TRIG-NC gene  -> the ASC family+gene of a clustered sibling allele
  #  2) else by TRIG-NC family -> the ASC family, as its own thin gene arc
  #  3) else (no clustered member in the family at all) a small "unclustered" bucket
  # Most frequent non-NA value, always a length-1 character (never NULL, so the
  # data.table type probe is happy even when the grouped subset is empty, e.g.
  # when no name map is supplied and trig_gene is all NA).
  .safe_mode <- function(x) { x <- x[!is.na(x)]
    if (!length(x)) NA_character_ else names(which.max(table(x)))[1L] }
  cl <- d[is_unclust == FALSE & !is.na(cluster_id) & !is.na(family_cluster)]
  if (nrow(cl) && cl[!is.na(trig_gene), .N] > 0L) {
    plg <- cl[!is.na(trig_gene),
              .(pf = .safe_mode(family_cluster),
                pg = .safe_mode(cluster_id)), by = trig_gene]
    d <- merge(d, plg, by = "trig_gene", all.x = TRUE, sort = FALSE)
    d[is_unclust & is.na(family_cluster) & !is.na(pf), family_cluster := pf]
    d[is_unclust & is.na(cluster_id)     & !is.na(pg), cluster_id     := pg]
    d[, c("pf", "pg") := NULL]
  }
  if (nrow(cl) && cl[!is.na(trig_fam), .N] > 0L) {
    plf <- cl[!is.na(trig_fam),
              .(pf2 = .safe_mode(family_cluster)), by = trig_fam]
    d <- merge(d, plf, by = "trig_fam", all.x = TRUE, sort = FALSE)
    d[is_unclust & is.na(family_cluster) & !is.na(pf2), family_cluster := pf2]
    d[is_unclust & is.na(cluster_id) & !is.na(family_cluster) & !is.na(trig_gene),
      cluster_id := trig_gene]
    d[, pf2 := NULL]
  }
  d[is.na(family_cluster), family_cluster := "unclustered"]
  d[is.na(cluster_id),     cluster_id     := "unclustered"]

  # Leaf order: family (by size) -> gene cluster (by size) -> allele order set
  # by a k-mer similarity tree inside each gene cluster.
  fam_sz  <- d[, .N, by = family_cluster][order(-N)]
  gene_sz <- d[, .N, by = .(family_cluster, cluster_id)][order(-N)]
  d[, fam_rank  := match(family_cluster, fam_sz$family_cluster)]
  d[, gene_rank := match(paste(family_cluster, cluster_id),
                         paste(gene_sz$family_cluster, gene_sz$cluster_id))]
  setorder(d, fam_rank, gene_rank, imgt_allele)
  d[, kpos := kmer_rank(sequence), by = .(family_cluster, cluster_id)]
  setorder(d, fam_rank, gene_rank, kpos)
  d[, leaf := .I]

  # Family angular spans proportional to allele counts, separated by a gap.
  gap     <- opt[["gap-deg"]] * pi / 180
  fam_tab <- d[, .(n = .N, first = min(leaf)), by = family_cluster]
  setorder(fam_tab, first)
  n_fam   <- nrow(fam_tab)
  usable  <- 2*pi - gap * n_fam
  fam_tab[, span := usable * n / sum(n)]
  fam_tab[, a0 := cumsum(c(0, head(span, -1))) + (seq_len(.N) - 1) * gap]
  fam_tab[, a1 := a0 + span]

  d <- merge(d, fam_tab[, .(family_cluster, fa0 = a0, fa1 = a1, fn = n)],
             by = "family_cluster", all.x = TRUE)
  setorder(d, leaf)
  d[, k := seq_len(.N), by = family_cluster]
  d[, angle := fa0 + (fa1 - fa0) * (k - 0.5) / fn]

  # Track radii
  r_fam   <- c(1.34, 1.46)   # outer family ring
  r_gene  <- c(1.19, 1.31)   # middle gene ring
  r_leaf0 <- 1.06            # allele tick base
  r_leaf1 <- 1.15            # allele tick tip / marker
  r_lab   <- 1.52            # family number labels
  r_txt   <- 1.17            # optional allele labels
  # Central dendrogram radii
  r_leafin <- 1.05
  r_genen  <- 0.72
  r_famn   <- 0.40
  r_root   <- 0.14

  fam_col  <- setNames(pal(nrow(fam_tab)), fam_tab$family_cluster)
  gene_key <- unique(d[, .(family_cluster, cluster_id)])
  gene_col <- setNames(pal(nrow(gene_key)),
                       paste(gene_key$family_cluster, gene_key$cluster_id))

  # TR/IG NC roll-ups for the ring hovers.
  fam_trig_lk  <- { x <- d[!is.na(trig_fam),
                           .(tf = paste(sort(unique(trig_fam)), collapse = ", ")),
                           by = family_cluster]
                    setNames(x$tf, x$family_cluster) }
  gene_trig_lk <- { x <- d[!is.na(trig_gene),
                           .(tg = paste(sort(unique(trig_gene)), collapse = ", ")),
                           by = .(family_cluster, cluster_id)]
                    setNames(x$tg, paste(x$family_cluster, x$cluster_id)) }

  p <- plot_ly()

  # ── Central radial dendrogram (allele -> gene -> family -> root) ────────────
  if (!isTRUE(opt[["hide-dendrogram"]])) {
    geneN <- d[, .(a0 = min(angle), a1 = max(angle), gang = mean(angle)),
               by = .(family_cluster, cluster_id)]
    famN  <- geneN[, .(g0 = min(gang), g1 = max(gang), fang = mean(gang)),
                   by = family_cluster]

    spoke <- function(r1, r2, a)
      list(x = as.vector(rbind(r1 * cos(a), r2 * cos(a), NA)),
           y = as.vector(rbind(r1 * sin(a), r2 * sin(a), NA)))
    L <- spoke(r_leafin, r_genen, d$angle)       # allele spokes
    G <- spoke(r_genen,  r_famn,  geneN$gang)     # gene-node spokes
    F <- spoke(r_famn,   r_root,  famN$fang)      # family-node spokes
    p <- add_trace(p, x = c(L$x, G$x, F$x), y = c(L$y, G$y, F$y),
                   type = "scatter", mode = "lines",
                   line = list(color = "#9AA0A6", width = 0.6),
                   hoverinfo = "skip", showlegend = FALSE)

    ax <- numeric(0); ay <- numeric(0)
    addarc <- function(r, a0, a1) {
      n  <- max(2L, as.integer(round(abs(a1 - a0) / (pi / 90))))
      th <- seq(a0, a1, length.out = n)
      ax <<- c(ax, r * cos(th), NA); ay <<- c(ay, r * sin(th), NA)
    }
    for (i in seq_len(nrow(geneN))) addarc(r_genen, geneN$a0[i], geneN$a1[i])
    for (i in seq_len(nrow(famN)))  addarc(r_famn,  famN$g0[i],  famN$g1[i])
    if (nrow(famN) > 1L) addarc(r_root, min(famN$fang), max(famN$fang))
    p <- add_trace(p, x = ax, y = ay, type = "scatter", mode = "lines",
                   line = list(color = "#9AA0A6", width = 0.6),
                   hoverinfo = "skip", showlegend = FALSE)
  }

  # ── Family ring (outer) ─────────────────────────────────────────────────────
  for (i in seq_len(nrow(fam_tab))) {
    b  <- arc_band(fam_tab$a0[i], fam_tab$a1[i], r_fam[1], r_fam[2])
    ft <- fam_trig_lk[fam_tab$family_cluster[i]]
    htxt <- if (trig_mode && length(ft) && !is.na(ft))
      sprintf("family %s, %d allele(s)", ft, fam_tab$n[i])
    else
      sprintf("family cluster %s, %d allele(s)%s",
              fam_tab$family_cluster[i], fam_tab$n[i],
              if (length(ft) && !is.na(ft))
                sprintf("<br>TRIG-NC family: %s", ft) else "")
    p <- add_polygons(p, x = b$x, y = b$y,
                      line = list(width = 0.5, color = "#FFFFFF"),
                      fillcolor = unname(fam_col[fam_tab$family_cluster[i]]),
                      hoverinfo = "skip", showlegend = FALSE)
    p <- hover_ring(p, fam_tab$a0[i], fam_tab$a1[i],
                    c(1.37, 1.40, 1.43), htxt)
    mid <- (fam_tab$a0[i] + fam_tab$a1[i]) / 2
    p <- add_text(p, x = r_lab * cos(mid), y = r_lab * sin(mid),
                  text = fam_tab$family_cluster[i],
                  textfont = list(size = 11, color = "#222222"),
                  hoverinfo = "skip", showlegend = FALSE)
  }

  # ── Gene-cluster ring (middle) ──────────────────────────────────────────────
  half <- d[, .(step = (fa1 - fa0)[1] / (2 * fn[1])), by = family_cluster]
  gsp  <- d[, .(a0 = min(angle), a1 = max(angle)),
            by = .(family_cluster, cluster_id)]
  gsp  <- merge(gsp, half, by = "family_cluster")
  gsp  <- merge(gsp, fam_tab[, .(family_cluster, fa0 = a0, fa1 = a1)],
                by = "family_cluster")
  gsp[, a0 := pmax(a0 - step, fa0)]
  gsp[, a1 := pmin(a1 + step, fa1)]
  for (i in seq_len(nrow(gsp))) {
    b   <- arc_band(gsp$a0[i], gsp$a1[i], r_gene[1], r_gene[2])
    key <- paste(gsp$family_cluster[i], gsp$cluster_id[i])
    gt  <- gene_trig_lk[key]
    htxt <- if (trig_mode && length(gt) && !is.na(gt))
      sprintf("gene %s", gt)
    else
      sprintf("gene cluster %s (family %s)%s",
              gsp$cluster_id[i], gsp$family_cluster[i],
              if (length(gt) && !is.na(gt))
                sprintf("<br>TRIG-NC gene: %s", gt) else "")
    p <- add_polygons(p, x = b$x, y = b$y,
                      line = list(width = 0.4, color = "#FFFFFF"),
                      fillcolor = unname(gene_col[key]),
                      hoverinfo = "skip", showlegend = FALSE)
    p <- hover_ring(p, gsp$a0[i], gsp$a1[i], c(1.22, 1.25, 1.28), htxt)
  }

  # ── Allele leaves: ticks + markers with hover ───────────────────────────────
  seqlen <- ifelse(is.na(d$sequence), NA_integer_, nchar(d$sequence))
  d[, asc_line  := fifelse(is.na(asc_name),  "", sprintf("<br>ASC: <b>%s</b>", asc_name))]
  d[, trig_line := fifelse(is.na(trig_name), "", sprintf("<br>TRIG-NC: <b>%s</b>", trig_name))]
  d[, stat_line := fifelse(is_unclust & !is.na(reason),
                           sprintf("<br>ASC status: <b>not clustered (%s)</b>", reason), "")]
  d[, leaf_lab  := fifelse(!is.na(trig_name), trig_name,
                   fifelse(!is.na(asc_name), asc_name, imgt_allele))]
  d[, hover := sprintf(
      paste0("<b>%s</b><br><i>original annotation</i>",
             "%s%s%s",
             "<br>family cluster: %s   gene cluster: %s",
             "<br>source: %s%s"),
      imgt_allele, asc_line, trig_line, stat_line, family_cluster, cluster_id,
      ifelse(is.na(source), "n/a", source),
      ifelse(is.na(seqlen), "", sprintf("<br>length: %d nt", seqlen)))]
  d[, lcol := src_symbol_colour(source)]

  tick_x <- as.vector(rbind(r_leaf0 * cos(d$angle), r_leaf1 * cos(d$angle), NA))
  tick_y <- as.vector(rbind(r_leaf0 * sin(d$angle), r_leaf1 * sin(d$angle), NA))
  p <- add_trace(p, x = tick_x, y = tick_y, type = "scatter", mode = "lines",
                 line = list(color = "#BBBBBB", width = 0.6),
                 hoverinfo = "skip", showlegend = FALSE)
  # Clustered alleles as filled dots; ASC-unclustered alleles as hollow dots so
  # the two are distinguishable at a glance while sitting in their real gene.
  dc <- d[is_unclust == FALSE]; du <- d[is_unclust == TRUE]
  if (nrow(dc))
    p <- add_markers(p, x = r_leaf1 * cos(dc$angle), y = r_leaf1 * sin(dc$angle),
                     marker = list(size = 6, color = dc$lcol,
                                   line = list(width = 0.4, color = "#333333")),
                     hovertemplate = paste0(dc$hover, "<extra></extra>"),
                     showlegend = FALSE)
  if (nrow(du))
    p <- add_markers(p, x = r_leaf1 * cos(du$angle), y = r_leaf1 * sin(du$angle),
                     marker = list(size = 7, symbol = "circle-open",
                                   color = du$lcol,
                                   line = list(width = 1.4, color = du$lcol)),
                     hovertemplate = paste0(du$hover, "<extra></extra>"),
                     showlegend = FALSE)

  # Optional allele labels for sparse gene clusters (uses the updated label).
  if (opt[["label-min"]] > 0) {
    d[, gene_n := .N, by = .(family_cluster, cluster_id)]
    lab <- d[gene_n <= opt[["label-min"]]]
    for (j in seq_len(nrow(lab))) {
      deg  <- lab$angle[j] * 180 / pi
      pos  <- if (deg > 90 && deg < 270) "middle left" else "middle right"
      p <- add_text(p, x = r_txt * cos(lab$angle[j]), y = r_txt * sin(lab$angle[j]),
                    text = lab$leaf_lab[j], textposition = pos,
                    textfont = list(size = 8, color = "#444444"),
                    hoverinfo = "skip", showlegend = FALSE)
    }
  }

  # Legend: provenance (filled) plus the hollow marker for ASC-unclustered.
  # Empty-vector traces carry the swatch without plotting a point.
  leg <- data.table(lbl = c("custom", "reference"),
                    col = c("#D55E00", "#0072B2"),
                    sym = c("circle", "circle"))
  if (d[is_unclust == TRUE, .N] > 0L)
    leg <- rbind(leg, data.table(lbl = "not ASC-clustered",
                                 col = "#555555", sym = "circle-open"))
  for (k in seq_len(nrow(leg)))
    p <- add_trace(p, x = numeric(0), y = numeric(0),
                   type = "scatter", mode = "markers", name = leg$lbl[k],
                   marker = list(size = 9, symbol = leg$sym[k], color = leg$col[k],
                                 line = list(width = if (leg$sym[k] == "circle-open") 1.4 else 0.4,
                                             color = leg$col[k])),
                   showlegend = TRUE, hoverinfo = "skip")

  ttl <- title %||% sprintf(
    "%s ASC allele-similarity clusters (hover a leaf or ring for details)",
    locus_val)
  ax0 <- list(visible = FALSE, showgrid = FALSE, zeroline = FALSE,
              range = c(-1.62, 1.62))
  p |> layout(
    title = list(text = ttl, x = 0.02, font = list(size = 14)),
    xaxis = c(ax0, list(scaleanchor = "y", scaleratio = 1)),
    yaxis = ax0,
    legend = list(orientation = "h", x = 0.02, y = -0.02),
    margin = list(l = 10, r = 10, t = 50, b = 10),
    hoverlabel = list(bgcolor = "white", bordercolor = "#333333",
                      font = list(size = 12, color = "#111111"), align = "left"),
    plot_bgcolor = "white", paper_bgcolor = "white")
}

# ── Drive ────────────────────────────────────────────────────────────────────
out_for <- function(loc) {
  if (!is.null(opt$out)) return(opt$out)
  base <- sub("\\.tsv$", "", opt$annot)
  sprintf("%s_%s_circos.html", base, loc)
}

loci_to_plot <- if (!is.null(opt$locus)) opt$locus else sort(unique(annot$locus))
for (loc in loci_to_plot) {
  p <- build_locus_plot(annot, loc, title = opt$title)
  if (is.null(p)) {
    message("  [skip] ", loc, ": no ASC cluster structure (D/J loci are not ",
            "clustered, or no rows for this locus)")
    next
  }
  of <- out_for(loc)

  # selfcontained = TRUE needs pandoc to inline the dependencies. When pandoc
  # is absent saveWidget errors, so fall back to a plain .html plus a sibling
  # "<name>_files/" directory that holds its JS/CSS. The pair opens the same
  # way in a browser; keep them together when moving the file.
  saved <- tryCatch({
    saveWidget(p, of, selfcontained = TRUE, title = paste0(loc, " ASC circos"))
    TRUE
  }, error = function(e) FALSE)

  if (saved) {
    cat(sprintf("Wrote %s (self-contained)\n", of))
  } else {
    saveWidget(p, of, selfcontained = FALSE, title = paste0(loc, " ASC circos"))
    libdir <- sub("\\.html$", "_files", basename(of))
    cat(sprintf(paste0("Wrote %s (pandoc not found, so not self-contained)\n",
                       "  Keep the sibling directory '%s' next to the .html ",
                       "when moving or sharing it.\n"),
                of, libdir))
  }
  if (!is.null(opt$locus)) break
}
