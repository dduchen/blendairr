# ------------------------------------------------------------------------------
# Visualization functions for allele clustering
# Provides truncated tree plots, network plots, and silhouette optimization plots

#' @include piglet.R
#' @include allele_cluster.R
#' @include community_detection.R
NULL

# ------------------------------------------------------------------------------

#' Plot truncated tree visualization
#'
#' Creates a circular or dendrogram tree visualization collapsed to ASC subgroup level,
#' with optional heatmap annotations showing family assignments.
#'
#' @param x A GermlineCluster object from \code{\link{inferAlleleClusters}}
#' @param layout Tree layout: "circular" (default) or "dendrogram"
#' @param collapse_to Level to collapse tree: "asc_subgroup" (default, based on ASC names),
#'   "iuis_subgroup" (based on original IUIS gene names), or "family"
#' @param label_style Label style for tips: "asc" (default, show ASC names like IGHVF1-G1),
#'   "iuis" (show IUIS names with superscript markers if ASC splits IUIS group),
#'   or "both" (show both names)
#' @param show_threshold_line Logical. Show threshold line on tree. Default is TRUE.
#' @param threshold Threshold height for threshold line (0-1 scale). Default is 0.25.
#' @param tip_size_by Variable for tip point size: "n_alleles" (default), "fixed", or NULL
#' @param tip_color_by Variable for tip point color: "present" (default), "fraction_novel", or NULL
#' @param show_heatmap Logical. Show heatmap annotation for IUIS vs ASC families. Default is TRUE.
#' @param label_size Size of tip labels. Default is 7.
#' @param ... Additional arguments passed to ggtree
#'
#' @return A ggplot/ggtree object
#'
#' @details
#' This function creates a publication-quality tree visualization that:
#' \itemize{
#'   \item Renames tree tips from original allele names to ASC names (new_allele)
#'   \item Collapses alleles to ASC subgroup level (single representative per ASC group)
#'   \item Shows tip point size by number of alleles in cluster
#'   \item Adds optional heatmap track showing IUIS vs ASC family assignments
#'   \item Draws threshold line at specified height
#' }
#'
#' When using \code{label_style = "iuis"}, if multiple ASC groups split a single IUIS
#' subgroup, the labels are marked with superscript letters (e.g., IGHV1-2^A, IGHV1-2^B)
#' to distinguish them.
#'
#' Requires the ggtree package to be installed.
#'
#' @seealso \code{\link{inferAlleleClusters}}, \code{\link{plot.GermlineCluster}}
#'
#' @examples
#' \donttest{
#' data(HVGERM)
#' asc <- inferAlleleClusters(HVGERM[1:50])
#'
#' # Basic truncated tree with ASC labels
#' if (requireNamespace("ggtree", quietly = TRUE)) {
#'   plotTruncatedTree(asc, show_heatmap = FALSE)
#'
#'   # With IUIS labels (marked if ASC splits IUIS group)
#'   plotTruncatedTree(asc, label_style = "iuis", show_heatmap = FALSE)
#' }
#' }
#'
#' @export
plotTruncatedTree <- function(x,
                              layout = c("circular", "dendrogram"),
                              collapse_to = c("asc_subgroup", "iuis_subgroup", "family"),
                              label_style = c("asc", "iuis", "both"),
                              show_threshold_line = TRUE,
                              threshold = 0.25,
                              tip_size_by = "n_alleles",
                              tip_color_by = "present",
                              show_heatmap = TRUE,
                              label_size = 7,
                              ...) {

  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop("Package 'ggtree' is required for this function. ",
         "Install it from Bioconductor: BiocManager::install('ggtree')")
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for this function. ",
         "Install it with: install.packages('ape')")
  }

  if (!inherits(x, "GermlineCluster"))
    stop("Object must be of class GermlineCluster")

  if (is.null(x$hclustAlleleCluster))
    stop("Truncated tree plot requires hierarchical clustering. ",
         "Use clustering_method = 'hierarchical' in inferAlleleClusters().")

  layout <- match.arg(layout)
  collapse_to <- match.arg(collapse_to)
  label_style <- match.arg(label_style)

  ## get hclust object and allele table

  hc <- x$hclustAlleleCluster
  allele_table <- .compat_allele_table(as.data.frame(x$alleleClusterTable))

  ## extract ASC subgroup from new_allele (e.g., IGHVF1-G1 from IGHVF1-G1*01)
  allele_table$asc_subgroup <- sapply(allele_table$new_allele, function(a) {
    gsub("[*].*$", "", a)
  })

  ## get IUIS subgroup from iuis_allele
  allele_table$iuis_subgroup <- alakazam::getGene(
    allele_table$iuis_allele,
    strip_d = FALSE,
    omit_nl = FALSE
  )

  ## get IUIS family from iuis_allele
  allele_table$iuis_family <- alakazam::getFamily(
    allele_table$iuis_allele,
    strip_d = FALSE,
    omit_nl = FALSE
  )

  ## get ASC family from new_allele (e.g., IGHVF1 from IGHVF1-G1*01)
  allele_table$asc_family <- sapply(allele_table$new_allele, function(a) {
    parts <- strsplit(gsub("[*].*$", "", a), "-")[[1]]
    if (length(parts) > 0) parts[1] else a
  })

  ## Create mapping from iuis_allele to new_allele for renaming hclust labels
  imgt_to_asc <- stats::setNames(allele_table$new_allele, allele_table$iuis_allele)

  ## Rename hclust labels to ASC names
  hc_asc <- hc
  hc_asc$labels <- unname(imgt_to_asc[hc$labels])

  ## Convert to phylo for manipulation
  tree <- ape::as.phylo(hc_asc)

  ## Determine grouping for collapse
  if (collapse_to == "asc_subgroup") {
    ## Group by ASC subgroup (from new_allele)
    group_col <- "asc_subgroup"
  } else if (collapse_to == "iuis_subgroup") {
    ## Group by IUIS subgroup (from iuis_allele)
    group_col <- "iuis_subgroup"
  } else {
    ## Group by ASC family
    group_col <- "asc_family"
  }

  ## Create mapping from new_allele to group
  asc_to_group <- stats::setNames(allele_table[[group_col]], allele_table$new_allele)

  ## Get tips grouped by the chosen grouping
  asc_tips <- tapply(allele_table$new_allele, allele_table[[group_col]], function(a) a)

  ## collapse clades - keep first tip of each group
  collapse_clade <- function(tree_, tips) {
    if (length(tips) > 1) {
      ape::drop.tip(tree_, tip = tips[-1])
    } else {
      tree_
    }
  }

  collapsed_tree <- Reduce(function(tree_, tips) collapse_clade(tree_, tips),
                           asc_tips, init = tree)

  ## convert back to hclust for ggtree
  hc_collapsed <- ape::as.hclust.phylo(collapsed_tree)

  ## Get tip labels (these are now ASC names)
  tip_labels_asc <- collapsed_tree$tip.label

  ## Create tip metadata
  tip_metadata <- data.frame(
    label = tip_labels_asc,
    stringsAsFactors = FALSE
  )

  ## Map back to get info from allele_table
  ## For each tip (ASC name), find the corresponding row in allele_table
  tip_metadata$asc_subgroup <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$asc_subgroup[idx[1]] else gsub("[*].*$", "", lab)
  })

  tip_metadata$iuis_subgroup <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$iuis_subgroup[idx[1]] else NA
  })

  tip_metadata$iuis_family <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$iuis_family[idx[1]] else NA
  })

  tip_metadata$asc_family <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$asc_family[idx[1]] else NA
  })

  tip_metadata$family <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$family[idx[1]] else NA
  })

  tip_metadata$allele_cluster <- sapply(tip_labels_asc, function(lab) {
    idx <- which(allele_table$new_allele == lab)
    if (length(idx) > 0) allele_table$allele_cluster[idx[1]] else NA
  })

  ## count alleles per group
  if (collapse_to == "asc_subgroup") {
    allele_counts <- table(allele_table$asc_subgroup)
    tip_metadata$n_alleles <- as.integer(allele_counts[tip_metadata$asc_subgroup])
  } else if (collapse_to == "iuis_subgroup") {
    allele_counts <- table(allele_table$iuis_subgroup)
    tip_metadata$n_alleles <- as.integer(allele_counts[tip_metadata$iuis_subgroup])
  } else {
    allele_counts <- table(allele_table$asc_family)
    tip_metadata$n_alleles <- as.integer(allele_counts[tip_metadata$asc_family])
  }

  ## Create display labels based on label_style
  if (label_style == "asc") {
    ## Use ASC subgroup names directly
    tip_metadata$display_label <- tip_metadata$asc_subgroup
  } else if (label_style == "iuis") {
    ## Use IUIS names, but mark with superscripts if ASC splits IUIS group
    tip_metadata$display_label <- .create_iuis_labels_with_markers(
      tip_metadata$iuis_subgroup,
      tip_metadata$asc_subgroup
    )
  } else {
    ## Both: show "IUIS (ASC)"
    tip_metadata$display_label <- paste0(
      tip_metadata$iuis_subgroup, " (", tip_metadata$asc_subgroup, ")"
    )
  }

  ## create size bins
  if (!is.null(tip_size_by) && tip_size_by == "n_alleles") {
    max_alleles <- max(tip_metadata$n_alleles, na.rm = TRUE)
    if (max_alleles > 1) {
      breaks <- 2^(0:ceiling(log2(max_alleles)))
      tip_metadata$size_bin <- cut(tip_metadata$n_alleles,
                                   breaks = breaks,
                                   include.lowest = TRUE)
    } else {
      tip_metadata$size_bin <- factor(rep("1", nrow(tip_metadata)))
    }
  }

  ## create base tree
  if (layout == "circular") {
    p <- ggtree::ggtree(hc_collapsed, layout = "circular", linewidth = 1, ...)
  } else {
    p <- ggtree::ggtree(hc_collapsed, layout = "dendrogram", linewidth = 1, ...)
  }

  ## add metadata to tree
  p <- ggtree::`%<+%`(p, tip_metadata)

  ## add tip points with size mapping
  if (!is.null(tip_size_by) && tip_size_by == "n_alleles") {
    p <- p + ggtree::geom_tippoint(ggplot2::aes(size = !!rlang::sym("size_bin")),
                                   alpha = 0.7, color = "gray40")

    ## set size scale
    bin_labels <- levels(tip_metadata$size_bin)
    sizes <- seq(3, 10, length.out = length(bin_labels))
    p <- p + ggplot2::scale_size_manual(values = sizes, labels = bin_labels,
                                        name = "# Alleles")
  } else if (!is.null(tip_size_by) && tip_size_by == "fixed") {
    p <- p + ggtree::geom_tippoint(size = 3, alpha = 0.7, color = "gray40")
  }

  r_max <- max(p$data$x, na.rm = TRUE)
  hm_width <- if (show_heatmap) 0.10 else 0
  hm_thickness <- r_max * hm_width
  max_chars <- max(nchar(tip_metadata$display_label), na.rm = TRUE)
  char_to_radius <- 0.012
  label_room <- r_max * (max_chars * label_size * char_to_radius)
  heatmap_offset <- r_max * 0.02
  label_offset   <- heatmap_offset + hm_thickness + r_max * 0.03
  x_expand <- label_offset + label_room + r_max * 0.05

  ## add threshold line
  if (show_threshold_line && layout == "circular") {
    max_height <- attr(as.dendrogram(hc_collapsed), "height")
    radius_cutoff <- max_height * (1 - threshold)
    p <- p + ggplot2::geom_vline(xintercept = radius_cutoff,
                                 colour = "firebrick", linetype = "dashed")
  }

  ## add heatmap if requested
  if (show_heatmap && requireNamespace("ggtree", quietly = TRUE)) {
    ## Create heatmap data frame with IUIS family and ASC family
    ## Row names must match the tree tip labels exactly
    heatmap_df <- data.frame(
      IUIS = as.character(tip_metadata$iuis_family),
      ASC = as.character(tip_metadata$asc_family),
      stringsAsFactors = FALSE
    )
    rownames(heatmap_df) <- tip_metadata$label

    ## Try to add heatmap (requires gheatmap from ggtree)
    ## Wrap in tryCatch since gheatmap can be sensitive to tree structure
    if (exists("gheatmap", where = asNamespace("ggtree"), mode = "function")) {
      p <- tryCatch({
        ggtree::gheatmap(p, heatmap_df,
                         offset = heatmap_offset,
                         width = hm_width,
                         colnames_angle = 90,
                         colnames_offset_y = 0.5)
      }, error = function(e) {
        warning("Could not add heatmap annotation: ", conditionMessage(e))
        p
      })
    }
  }

  p <- p + ggtree::geom_tiplab(ggplot2::aes(label = !!rlang::sym("display_label")),
                               size = label_size, offset = label_offset)

  

  return(p)
}

#' Create IUIS labels with markers for split groups
#'
#' Internal function to create IUIS labels with superscript markers
#' when multiple ASC groups split a single IUIS subgroup.
#'
#' @param iuis_subgroups Vector of IUIS subgroup names
#' @param asc_subgroups Vector of corresponding ASC subgroup names
#'
#' @return Character vector of labels with markers
#'
#' @keywords internal
.create_iuis_labels_with_markers <- function(iuis_subgroups, asc_subgroups) {
  ## Find IUIS groups that are split into multiple ASC groups
  iuis_to_asc <- tapply(asc_subgroups, iuis_subgroups, function(x) unique(x))

  ## Create labels
  labels <- character(length(iuis_subgroups))

  for (i in seq_along(iuis_subgroups)) {
    iuis <- iuis_subgroups[i]
    asc <- asc_subgroups[i]

    if (is.na(iuis)) {
      labels[i] <- asc
      next
    }

    asc_groups <- iuis_to_asc[[iuis]]

    if (length(asc_groups) > 1) {
      ## This IUIS group is split - add superscript marker
      marker_idx <- which(asc_groups == asc)
      marker <- LETTERS[marker_idx]
      labels[i] <- paste0(iuis, "^", marker)
    } else {
      ## Single ASC group for this IUIS - no marker needed
      labels[i] <- iuis
    }
  }

  return(labels)
}

# ------------------------------------------------------------------------------

#' Plot community network
#'
#' Creates a network visualization of allele clusters from community detection.
#'
#' @param x A GermlineCluster object with Leiden clustering
#' @param layout Network layout: "fr" (Fruchterman-Reingold, default), "kk" (Kamada-Kawai), or "circle"
#' @param node_color Variable for node color: "cluster" (default), "family", or a color value
#' @param node_size Variable for node size: "degree" (default), "fixed", or a numeric value
#' @param edge_alpha Alpha transparency for edges. Default is 0.3.
#' @param show_labels Logical. Show node labels. Default is TRUE.
#' @param label_size Size of node labels. Default is 3.
#' @param ... Additional arguments
#'
#' @return A ggplot object
#'
#' @details
#' This function creates a network visualization showing:
#' \itemize{
#'   \item Nodes representing alleles, colored by cluster
#'   \item Edges weighted by sequence similarity
#'   \item Layout optimized by specified algorithm
#' }
#'
#' @seealso \code{\link{inferAlleleClusters}}, \code{\link{detect_communities_leiden}}
#'
#' @examples
#' \donttest{
#' data(HVGERM)
#' asc <- inferAlleleClusters(HVGERM[1:30],
#'                            clustering_method = "leiden",
#'                            target_clusters = 5)
#' plotCommunityNetwork(asc)
#' }
#'
#' @export
plotCommunityNetwork <- function(x,
                                 layout = c("fr", "kk", "circle"),
                                 node_color = "cluster",
                                 node_size = "degree",
                                 edge_alpha = 0.3,
                                 show_labels = TRUE,
                                 label_size = 3,
                                 ...) {

  if (!inherits(x, "GermlineCluster"))
    stop("Object must be of class GermlineCluster")

  if (is.null(x$graphObject))
    stop("Network plot requires Leiden clustering with graph object. ",
         "Use clustering_method = 'leiden' in inferAlleleClusters().")

  layout <- match.arg(layout)

  g <- x$graphObject
  allele_table <- .compat_allele_table(as.data.frame(x$alleleClusterTable))

  ## create layout
  if (layout == "fr") {
    layout_coords <- igraph::layout_with_fr(g)
  } else if (layout == "kk") {
    layout_coords <- igraph::layout_with_kk(g)
  } else {
    layout_coords <- igraph::layout_in_circle(g)
  }

  ## prepare node data
  node_data <- data.frame(
    name = igraph::V(g)$name,
    x = layout_coords[, 1],
    y = layout_coords[, 2],
    stringsAsFactors = FALSE
  )

  ## add cluster info
  cluster_map <- stats::setNames(allele_table$allele_cluster, allele_table$iuis_allele)
  node_data$cluster <- as.factor(cluster_map[node_data$name])

  ## add family info
  node_data$family <- sapply(node_data$name, function(n) {
    alakazam::getFamily(n, strip_d = FALSE, omit_nl = FALSE)
  })

  ## calculate degree for node size
  node_data$degree <- igraph::degree(g)

  ## prepare edge data
  edges <- igraph::as_data_frame(g, what = "edges")
  edges$x <- layout_coords[match(edges$from, igraph::V(g)$name), 1]
  edges$y <- layout_coords[match(edges$from, igraph::V(g)$name), 2]
  edges$xend <- layout_coords[match(edges$to, igraph::V(g)$name), 1]
  edges$yend <- layout_coords[match(edges$to, igraph::V(g)$name), 2]

  ## create plot
  p <- ggplot2::ggplot()

  ## add edges
  p <- p + ggplot2::geom_segment(
    data = edges,
    ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), xend = !!rlang::sym("xend"), yend = !!rlang::sym("yend")),
    alpha = edge_alpha,
    color = "gray60"
  )

  ## add nodes
  if (node_color == "cluster") {
    if (is.numeric(node_size) || node_size == "fixed") {
      size_val <- if (is.numeric(node_size)) node_size else 5
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = !!rlang::sym("cluster")),
        size = size_val
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = !!rlang::sym("cluster"), size = !!rlang::sym("degree"))
      )
    }
    p <- p + ggplot2::scale_color_discrete(name = "Cluster")
  } else if (node_color == "family") {
    if (is.numeric(node_size) || node_size == "fixed") {
      size_val <- if (is.numeric(node_size)) node_size else 5
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = !!rlang::sym("family")),
        size = size_val
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = !!rlang::sym("family"), size = !!rlang::sym("degree"))
      )
    }
    p <- p + ggplot2::scale_color_discrete(name = "Family")
  } else {
    if (is.numeric(node_size) || node_size == "fixed") {
      size_val <- if (is.numeric(node_size)) node_size else 5
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y")),
        color = node_color,
        size = size_val
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = node_data,
        ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), size = !!rlang::sym("degree")),
        color = node_color
      )
    }
  }

  ## add labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = node_data,
      ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"), label = !!rlang::sym("name")),
      size = label_size,
      vjust = -1
    )
  }

  ## theme
  p <- p + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")

  if (node_size == "degree") {
    p <- p + ggplot2::scale_size_continuous(name = "Degree", range = c(2, 8))
  }

  return(p)
}

# ------------------------------------------------------------------------------

#' Plot silhouette optimization results
#'
#' Creates a plot showing silhouette score and cluster count across resolution values.
#'
#' @param optimization_result Result from \code{\link{optimize_resolution}}
#' @param highlight_best Logical. Highlight optimal resolution. Default is TRUE.
#' @param ... Additional arguments
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{optimize_resolution}}, \code{\link{igClust}}
#'
#' @examples
#' \donttest{
#' data(HVGERM)
#' d <- igDistance(HVGERM[1:30], method = "hamming")
#' g <- distance_to_graph(d)
#' opt <- optimize_resolution(g, d, target_clusters = 5)
#' plotSilhouetteOptimization(opt)
#' }
#'
#' @export
plotSilhouetteOptimization <- function(optimization_result,
                                       highlight_best = TRUE,
                                       ...) {

  if (!is.list(optimization_result) || !"results" %in% names(optimization_result))
    stop("Input must be a result from optimize_resolution()")

  results <- optimization_result$results
  best_res <- optimization_result$best_resolution

  ## create plot
  p <- ggplot2::ggplot(results, ggplot2::aes(x = !!rlang::sym("Resolution")))

  ## add silhouette line
  p <- p + ggplot2::geom_line(
    ggplot2::aes(y = !!rlang::sym("Silhouette")),
    color = "steelblue",
    linewidth = 1
  )
  p <- p + ggplot2::geom_point(
    ggplot2::aes(y = !!rlang::sym("Silhouette")),
    color = "steelblue",
    size = 2
  )

  ## highlight best
  if (highlight_best) {
    best_row <- results[results$Resolution == best_res, ]
    p <- p + ggplot2::geom_vline(
      xintercept = best_res,
      linetype = "dashed",
      color = "red"
    )
    p <- p + ggplot2::geom_point(
      data = best_row,
      ggplot2::aes(x = !!rlang::sym("Resolution"), y = !!rlang::sym("Silhouette")),
      color = "red",
      size = 4
    )
  }

  ## add cluster count on secondary axis
  scale_factor <- max(results$Silhouette, na.rm = TRUE) /
                  max(results$ClusterCount, na.rm = TRUE)

  p <- p + ggplot2::geom_line(
    ggplot2::aes(y = !!rlang::sym("ClusterCount") * scale_factor),
    color = "darkorange",
    linetype = "dotted",
    linewidth = 0.8
  )

  ## labels and theme
  p <- p + ggplot2::scale_y_continuous(
    name = "Silhouette Score",
    sec.axis = ggplot2::sec_axis(
      ~ . / scale_factor,
      name = "Cluster Count"
    )
  )

  p <- p + ggplot2::labs(
    x = "Resolution Parameter",
    title = "Silhouette Optimization"
  )

  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y.left = ggplot2::element_text(color = "steelblue"),
      axis.text.y.left = ggplot2::element_text(color = "steelblue"),
      axis.title.y.right = ggplot2::element_text(color = "darkorange"),
      axis.text.y.right = ggplot2::element_text(color = "darkorange")
    )

  return(p)
}

# ------------------------------------------------------------------------------

#' Compare hierarchical and Leiden clustering
#'
#' Creates a comparison visualization showing cluster assignments from both methods.
#'
#' @param hierarchical_result GermlineCluster object from hierarchical clustering
#' @param leiden_result GermlineCluster object from Leiden clustering
#' @param ... Additional arguments
#'
#' @return A ggplot object showing cluster agreement
#'
#' @seealso \code{\link{inferAlleleClusters}}
#'
#' @export
plotClusterComparison <- function(hierarchical_result,
                                  leiden_result,
                                  ...) {

  if (!inherits(hierarchical_result, "GermlineCluster") ||
      !inherits(leiden_result, "GermlineCluster"))
    stop("Both inputs must be GermlineCluster objects")

  if (hierarchical_result$clusteringMethod != "hierarchical")
    stop("First argument must be from hierarchical clustering")

  if (leiden_result$clusteringMethod != "leiden")
    stop("Second argument must be from Leiden clustering")

  ## get cluster assignments
  hier_table <- .compat_allele_table(as.data.frame(hierarchical_result$alleleClusterTable))
  leid_table <- .compat_allele_table(as.data.frame(leiden_result$alleleClusterTable))

  ## merge by allele name
  comparison <- merge(
    hier_table[, c("iuis_allele", "allele_cluster")],
    leid_table[, c("iuis_allele", "allele_cluster")],
    by = "iuis_allele",
    suffixes = c("_hier", "_leid")
  )

  ## create confusion matrix
  confusion <- table(
    Hierarchical = comparison$allele_cluster_hier,
    Leiden = comparison$allele_cluster_leid
  )

  ## convert to data frame for plotting
  conf_df <- as.data.frame(confusion)

  ## create heatmap
  p <- ggplot2::ggplot(conf_df,
                       ggplot2::aes(
                        x = !!rlang::sym("Leiden"), 
                        y = !!rlang::sym("Hierarchical"), 
                        fill = !!rlang::sym("Freq")))

  p <- p + ggplot2::geom_tile()

  p <- p + ggplot2::geom_text(
    ggplot2::aes(label = !!rlang::sym("Freq")),
    color = "white"
  )

  p <- p + ggplot2::scale_fill_gradient(
    low = "gray90",
    high = "steelblue",
    name = "Count"
  )

  p <- p + ggplot2::labs(
    title = "Cluster Comparison: Hierarchical vs Leiden",
    x = "Leiden Clusters",
    y = "Hierarchical Clusters"
  )

  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  return(p)
}
