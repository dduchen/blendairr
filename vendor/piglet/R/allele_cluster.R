# ------------------------------------------------------------------------------
# Allele clusters functions
# The functions in this scripts are for generating the allele clusters from a given reference set

#' @include piglet.R
#' @include RcppExports.R
NULL

# ------------------------------------------------------------------------------

#### Classes ####

#' Create a GermlineCluster object
#'
#' \code{GermlineCluster} is an S3 class that stores the output of
#' \code{\link{inferAlleleClusters}}. It contains the allele cluster table,
#' clustering objects, and threshold parameters used for inference.
#'
#' @param germlineSet The original germline set provided.
#' @param alleleClusterSet The renamed germline set with allele clusters.
#' @param alleleClusterTable The allele cluster table.
#' @param threshold The threshold used for family and allele clusters.
#' @param hclustAlleleCluster A hierarchical clustering object for the germline set,
#'   or \code{NULL}.
#' @param clusteringMethod The clustering method used, either \code{"hierarchical"}
#'   or \code{"leiden"}.
#' @param communityObject A community detection object for Leiden clustering, or \code{NULL}.
#' @param graphObject An \pkg{igraph} graph object for Leiden clustering, or \code{NULL}.
#' @param distanceMatrix The distance matrix used for clustering, or \code{NULL}.
#' @param silhouetteScore The silhouette score for community detection.
#' @param resolutionParameter The resolution parameter used for Leiden clustering.
#' @param locus The locus identifier, for example \code{"IGHV"}, \code{"IGHD"}, \code{"IGHJ"}.
#'
#' @return An object of class \code{"GermlineCluster"}.
#'
#' @seealso \code{\link{inferAlleleClusters}}
#' @seealso \code{\link{GermlineCluster}}
#'
#' @export
#' @family GermlineCluster
new_germline_cluster <- function(germlineSet,
                                 alleleClusterSet,
                                 alleleClusterTable,
                                 threshold,
                                 hclustAlleleCluster = NULL,
                                 clusteringMethod = "hierarchical",
                                 communityObject = NULL,
                                 graphObject = NULL,
                                 distanceMatrix = NULL,
                                 silhouetteScore = NA_real_,
                                 resolutionParameter = NA_real_,
                                 locus = "IGHV",
                                 .familiesCut = NULL) {
  alleleClusterTable <- structure(
    as.data.frame(alleleClusterTable),
    class = c("AlleleClusterTable", "data.frame")
  )
  structure(
    list(
      germlineSet = germlineSet,
      alleleClusterSet = alleleClusterSet,
      alleleClusterTable = alleleClusterTable,
      threshold = threshold,
      hclustAlleleCluster = hclustAlleleCluster,
      clusteringMethod = clusteringMethod,
      communityObject = communityObject,
      graphObject = graphObject,
      distanceMatrix = distanceMatrix,
      silhouetteScore = silhouetteScore,
      resolutionParameter = resolutionParameter,
      locus = locus,
      .familiesCut = .familiesCut
    ),
    class = "GermlineCluster"
  )
}

#' Print method for GermlineCluster
#'
#' @param x A GermlineCluster object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.GermlineCluster <- function(x, ...) {
  cat("GermlineCluster Object\n")
  cat("======================\n")
  cat("Locus:", x$locus, "\n")
  cat("Clustering method:", x$clusteringMethod, "\n")

  cat("Number of alleles:", length(x$germlineSet), "\n")
  cat("Number of allele clusters:", length(unique(x$alleleClusterTable$allele_cluster)), "\n")
  cat("Number of family clusters:", length(unique(x$alleleClusterTable$family)), "\n")
  cat("Thresholds:\n")
  cat("  - Family:", x$threshold$family_threshold, "\n")
  cat("  - Allele cluster:", x$threshold$allele_cluster_threshold, "\n")
  if (!is.na(x$silhouetteScore)) {
    cat("Silhouette score:", round(x$silhouetteScore, 4), "\n")
  }
  if (!is.na(x$resolutionParameter)) {
    cat("Resolution parameter:", round(x$resolutionParameter, 4), "\n")
  }
  invisible(x)
}

#' Summary method for GermlineCluster
#'
#' @param object A GermlineCluster object
#' @param ... Additional arguments (ignored)
#'
#' @return A list with summary statistics
#' @export
summary.GermlineCluster <- function(object, ...) {
  list(
    locus = object$locus,
    clustering_method = object$clusteringMethod,
    n_alleles = length(object$germlineSet),
    n_clusters = length(unique(object$alleleClusterTable$allele_cluster)),
    n_families = length(unique(object$alleleClusterTable$family)),
    family_threshold = object$threshold$family_threshold,
    cluster_threshold = object$threshold$allele_cluster_threshold,
    silhouette_score = object$silhouetteScore,
    resolution_parameter = object$resolutionParameter
  )
}

# ------------------------------------------------------------------------------

#### AlleleClusterTable S3 class ####

#' Dollar accessor for AlleleClusterTable (backward-compatible deprecation)
#'
#' Intercepts access to the deprecated \code{imgt_allele} column name and
#' redirects it to \code{iuis_allele} with a deprecation warning.
#'
#' @param x An AlleleClusterTable object.
#' @param name Column name.
#' @keywords internal
#' @export
`$.AlleleClusterTable` <- function(x, name) {
  if (identical(name, "imgt_allele")) {
    .Deprecated(msg = "Column 'imgt_allele' has been renamed to 'iuis_allele'. Please update your code.")
    name <- "iuis_allele"
  } else if (identical(name, "Family")) {
    .Deprecated(msg = "Column 'Family' has been renamed to 'family'. Please update your code.")
    name <- "family"
  } else if (identical(name, "Allele_Cluster")) {
    .Deprecated(msg = "Column 'Allele_Cluster' has been renamed to 'allele_cluster'. Please update your code.")
    name <- "allele_cluster"
  }
  .subset2(x, name)
}

#' Double-bracket accessor for AlleleClusterTable (backward-compatible deprecation)
#'
#' Intercepts access to the deprecated \code{imgt_allele} column name and
#' redirects it to \code{iuis_allele} with a deprecation warning.
#'
#' @param x An AlleleClusterTable object.
#' @param i Index or column name.
#' @param ... Additional arguments.
#' @keywords internal
#' @export
`[[.AlleleClusterTable` <- function(x, i, ...) {
  if (is.character(i)) {
    if (identical(i, "imgt_allele")) {
      .Deprecated(msg = "Column 'imgt_allele' has been renamed to 'iuis_allele'. Please update your code.")
      i <- "iuis_allele"
    } else if (identical(i, "Family")) {
      .Deprecated(msg = "Column 'Family' has been renamed to 'family'. Please update your code.")
      i <- "family"
    } else if (identical(i, "Allele_Cluster")) {
      .Deprecated(msg = "Column 'Allele_Cluster' has been renamed to 'allele_cluster'. Please update your code.")
      i <- "allele_cluster"
    }
  }
  NextMethod()
}

#' Backward-compatible allele cluster table helper
#'
#' Renames the deprecated \code{imgt_allele} column to \code{iuis_allele} if present.
#' Called at the boundary of exported functions that accept an \code{alleleClusterTable}.
#'
#' @param tbl A data.frame that may contain an \code{imgt_allele} column.
#' @return The table with \code{imgt_allele} renamed to \code{iuis_allele} if needed.
#' @keywords internal
.compat_allele_table <- function(tbl) {
  if ("imgt_allele" %in% names(tbl) && !"iuis_allele" %in% names(tbl)) {
    .Deprecated(msg = "Column 'imgt_allele' has been renamed to 'iuis_allele'. Please update your alleleClusterTable.")
    names(tbl)[names(tbl) == "imgt_allele"] <- "iuis_allele"
  }
  if ("Family" %in% names(tbl) && !"family" %in% names(tbl)) {
    .Deprecated(msg = "Column 'Family' has been renamed to 'family'. Please update your alleleClusterTable.")
    names(tbl)[names(tbl) == "Family"] <- "family"
  }
  if ("Allele_Cluster" %in% names(tbl) && !"allele_cluster" %in% names(tbl)) {
    .Deprecated(msg = "Column 'Allele_Cluster' has been renamed to 'allele_cluster'. Please update your alleleClusterTable.")
    names(tbl)[names(tbl) == "Allele_Cluster"] <- "allele_cluster"
  }
  tbl
}

# ------------------------------------------------------------------------------

#' Alleles nucleotide position difference
#'
#' Compare the sequences of two alleles (reference and sample alleles)
#' and returns the differential nucleotide positions of the sample allele.
#'
#'
#' @param    reference_allele          The nucleotide sequence of the reference allele, character object.
#' @param    sample_allele             The nucleotide sequence of the sample allele, character object.
#' @param    position_threshold        A position from which to check for differential positions. If zero checks all position. Default to zero.
#' @param    snps                      If to return the SNP with the position (e.g., A2G where A is for the reference and G is for the sample.). If false returns just the positions. Default to True
#' @return
#' A \code{character} vector of the differential nucleotide positions of the sample allele.
#'
#' @examples
#' {
#' reference_allele = "AAGG"
#' sample_allele = "ATGA"
#'
#' # setting position_threshold = 0 will return all differences
#' diff <- allele_diff(reference_allele, sample_allele)
#' # "A2T", "G4A"
#' print(diff)
#'
#' # setting position_threshold = 3 will return the differences from position three onward
#' diff <- allele_diff(reference_allele, sample_allele, position_threshold = 3)
#' # "G4A"
#' print(diff)
#'
#' # setting snps = FALSE will return the differences as indices
#' diff <- allele_diff(reference_allele, sample_allele, snps = FALSE)
#' # 2, 4
#' print(diff)
#'
#' }
#' @details
#' The function utilizes c++ script to optimize the run time for large comparisons.
#'
#' @export
allele_diff <-
  function(reference_allele,
           sample_allele,
           position_threshold = 0,
           snps = TRUE) {
    germs <- c(reference_allele, sample_allele)
    
    if (snps)
      return(allele_diff_strings(germs, position_threshold))
    
    return(allele_diff_indices(germs, position_threshold))
    
  }

# ------------------------------------------------------------------------------

# alignSeqs <- function(germline_set){
#   if (!is.vector(germline_set, mode = "character"))
#     germline_set <- tigger::readIgFasta(germline_set)
#   
#   aligned_set <- msa::msa(germline_set, type = "dna")
#   xn <- as.character(unmasked(aligned_set))
#   setNames(unname(xn), names(xn))
# }


# ------------------------------------------------------------------------------
# beta cluster app
# run_piglet_app <- function() {
#   shinyAppDir(
#     system.file("shiny", package = "piglet"),
#     options = list(
#       width = "100%", height = 600
#     )
#   )
# }

# ------------------------------------------------------------------------------

#' Germline set alleles distance
#'
#' Calculates the distance between pairs of alleles based on their aligned germline sequences.
#' Supports multiple distance methods for different segment types.
#'
#' @param germline_set A character vector of aligned allele sequences. See details for curating options.
#' @param AA Logical (FALSE by default). If TRUE, calculate the distance based on amino acid sequences.
#' @param method Distance calculation method. One of:
#'   \itemize{
#'     \item \code{"decipher"}: Uses DECIPHER::DistanceMatrix (requires aligned sequences, best for V segments)
#'     \item \code{"hamming"}: Hamming distance (requires equal length, sequences padded if needed)
#'     \item \code{"lv"}: Levenshtein distance (handles variable length, good for D/J segments)
#'   }
#' @param trim_3prime Optional position to trim sequences from 3' end before distance calculation
#' @param return_type One of "matrix" (default) or "dist" to return a dist object
#' @param quiet Logical (TRUE by default). Suppress informative messages
#'
#' @return
#' A \code{matrix} or \code{dist} object of the computed distances between allele pairs.
#'
#' @details
#' The aligned IMGT IGHV allele germline set can be downloaded from the IMGT site
#' \url{https://www.imgt.org/} under the section genedb.
#'
#' For V segments, the "decipher" method is recommended as it handles alignment gaps properly.
#' For D and J segments which may have variable lengths, the "lv" (Levenshtein) method is appropriate.
#'
#' @seealso \code{\link{ighvDistance}} for backward compatibility wrapper
#'
#' @examples
#' data(HVGERM)
#' # Using DECIPHER method (default, for V segments)
#' d1 <- igDistance(HVGERM[1:10], method = "decipher")
#'
#' # Using Hamming distance
#' d2 <- igDistance(HVGERM[1:10], method = "hamming")
#'
#' # Using Levenshtein distance (good for D/J segments)
#' d3 <- igDistance(HVGERM[1:10], method = "lv")
#'
#' @export
igDistance <- function(germline_set,
                       AA = FALSE,
                       method = c("decipher", "hamming", "lv"),
                       trim_3prime = NULL,
                       return_type = c("matrix", "dist"),
                       quiet = TRUE) {

  method <- match.arg(method)
  return_type <- match.arg(return_type)

  if (!is.character(germline_set))
    stop("The input germline set is not in a character class.")

  if (length(germline_set) < 2)
    stop("At least two sequences are required for distance calculation.")

  ## optionally trim sequences
  if (!is.null(trim_3prime)) {
    stopifnot(is.numeric(trim_3prime), trim_3prime >= 1)
    germline_set <- substr(germline_set, 1, trim_3prime)
  }

  ## sort by name
  germline_set <- germline_set[order(names(germline_set))]

  if (method == "decipher") {
    ## DECIPHER method (original implementation)
    germline_set <-
      gsub("\\s", "N", format(germline_set, width = max(nchar(germline_set))))

    if (AA) {
      germline_set <- gsub("X", "-", germline_set)
      germline_set <- Biostrings::AAStringSet(germline_set)
    } else {
      germline_set <- gsub("[.]", "-", germline_set)
      germline_set <- Biostrings::DNAStringSet(germline_set)
    }

    germline_distance <-
      DECIPHER::DistanceMatrix(
        germline_set,
        includeTerminalGaps = FALSE,
        penalizeGapLetterMatches = TRUE,
        verbose = FALSE
      )

    if (return_type == "dist") {
      return(stats::as.dist(germline_distance))
    }
    return(germline_distance)

  } else {
    ## stringdist methods (hamming or levenshtein)

    if (method == "hamming") {
      ## pad sequences to equal length for Hamming
      if (!quiet) message(sprintf("Padding to max length %s.", max(nchar(germline_set))))
      germline_set <- gsub("\\s", "N", format(germline_set, width = max(nchar(germline_set))))
    }

    ## check if sequences have names
    useNames <- "names"
    if (is.null(names(germline_set)) || any(is.na(names(germline_set)))) {
      useNames <- "strings"
    }

    d <- stringdist::stringdistmatrix(germline_set, method = method, useNames = useNames)

    if (return_type == "dist") {
      return(stats::as.dist(d))
    }
    return(as.matrix(d))
  }
}

#' Germline set alleles distance (deprecated)
#'
#' This function is deprecated. Use \code{\link{igDistance}} instead.
#'
#' @param germline_set A character list of aligned IGHV allele sequences.
#' @param AA Logical (FALSE by default). If to calculate the distance based on amino acid sequences.
#'
#' @return A matrix of computed distances between allele pairs.
#'
#' @seealso \code{\link{igDistance}} for the current implementation
#'
#' @export
ighvDistance <- function(germline_set, AA = FALSE) {
  .Deprecated("igDistance")
  igDistance(germline_set, AA = AA, method = "decipher")
}

# ------------------------------------------------------------------------------

#' Allele similarity clustering
#'
#' Cluster the distance matrix to create allele clusters. Supports both
#' hierarchical clustering (default) and Leiden community detection.
#'
#' @param germline_distance A germline set distance matrix created by \code{\link{igDistance}}.
#' @param method Clustering method. One of "hierarchical" (default) or "leiden".
#' @param distance_method The distance method used to compute \code{germline_distance}.
#'   One of \code{"decipher"} (default), \code{"hamming"}, or \code{"lv"}.
#'   For \code{"decipher"}, thresholds are similarity percentages in \code{[0, 100]};
#'   for \code{"hamming"} and \code{"lv"}, thresholds are raw integer mismatch counts.
#' @param family_threshold The threshold for family-level grouping.
#'   For \code{distance_method = "decipher"}: similarity percentage (default 75).
#'   For \code{"hamming"} / \code{"lv"}: maximum number of mismatches (must be >= allele_cluster_threshold).
#' @param allele_cluster_threshold The threshold for allele-cluster-level grouping (hierarchical only).
#'   For \code{distance_method = "decipher"}: similarity percentage (default 95).
#'   For \code{"hamming"} / \code{"lv"}: maximum number of mismatches (must be <= family_threshold).
#' @param cluster_method The linkage method for hierarchical clustering (used for family assignment in both methods). Default is "complete".
#' @param resolution Resolution parameter for Leiden clustering. If NULL, will be optimized.
#' @param target_clusters Target number of clusters for Leiden optimization. Default is NULL.
#' @param optimize_silhouette Logical. Optimize resolution using silhouette score (Leiden only). Default is TRUE.
#' @param ncores Number of cores for parallel processing (Leiden only). Default is 1.
#' @param quiet Logical. Suppress messages. Default is FALSE.
#'
#' @return
#' A named list that includes:
#' \itemize{
#'   \item \code{alleleClusterTable}: data.frame of allele clusters
#'   \item \code{threshold}: list of threshold parameters
#'   \item \code{hclustAlleleCluster}: hierarchical clustering object (both methods)
#'   \item \code{communityObject}: community detection result (Leiden method)
#'   \item \code{graphObject}: igraph object (Leiden method)
#'   \item \code{silhouetteScore}: silhouette score (Leiden method)
#'   \item \code{resolutionParameter}: resolution used (Leiden method)
#' }
#'
#' @seealso \code{\link{igDistance}}, \code{\link{inferAlleleClusters}}
#'
#' @export
igClust <- function(germline_distance,
                    method = c("hierarchical", "leiden"),
                    distance_method = "decipher",
                    family_threshold = 75,
                    allele_cluster_threshold = 95,
                    cluster_method = "complete",
                    resolution = NULL,
                    target_clusters = NULL,
                    optimize_silhouette = TRUE,
                    ncores = 1,
                    quiet = FALSE) {

  method <- match.arg(method)

  ## check the class of the distance matrix
  if (!is.matrix(germline_distance) && !inherits(germline_distance, "dist"))
    stop("The input germline_distance must be a matrix or dist object.")

  if (is.matrix(germline_distance) && nrow(germline_distance) != ncol(germline_distance))
    stop("The input germline_distance must be a square matrix.")

  ## convert to matrix if needed for processing
  if (inherits(germline_distance, "dist")) {
    dist_obj <- germline_distance
    germline_distance <- as.matrix(germline_distance)
  } else {
    dist_obj <- stats::as.dist(germline_distance)
  }

  if (method == "hierarchical") {
    ## Original hierarchical clustering implementation
    return(.igClust_hierarchical(
      germline_distance = germline_distance,
      distance_method = distance_method,
      family_threshold = family_threshold,
      allele_cluster_threshold = allele_cluster_threshold,
      cluster_method = cluster_method
    ))
  } else {
    ## Leiden community detection
    return(.igClust_leiden(
      germline_distance = germline_distance,
      dist_obj = dist_obj,
      distance_method = distance_method,
      family_threshold = family_threshold,
      cluster_method = cluster_method,
      resolution = resolution,
      target_clusters = target_clusters,
      optimize_silhouette = optimize_silhouette,
      ncores = ncores,
      quiet = quiet
    ))
  }
}

## Internal function for hierarchical clustering
.igClust_hierarchical <- function(germline_distance,
                                  distance_method = "decipher",
                                  family_threshold,
                                  allele_cluster_threshold,
                                  cluster_method) {

  ## check that thresholds are valid
  if (distance_method == "decipher") {
    if (any(c(family_threshold, allele_cluster_threshold) < 0 |
            c(family_threshold, allele_cluster_threshold) > 100))
      stop("Thresholds must be between 0 and 100 for the 'decipher' distance method.")
  } else {
    if (any(c(family_threshold, allele_cluster_threshold) < 0))
      stop("Thresholds must be non-negative.")
  }

  ## for decipher (similarity %): higher = stricter, so family (75) must be <= cluster (95)
  ## for hamming/lv (mismatch counts): higher = looser, so family (e.g. 10) must be >= cluster (e.g. 2)
  wrong_order <- if (distance_method == "decipher") {
    family_threshold > allele_cluster_threshold
  } else {
    family_threshold < allele_cluster_threshold
  }
  if (wrong_order) {
    message("Switching thresholds so that family is more inclusive than allele cluster.")
    tmp                      <- family_threshold
    family_threshold         <- allele_cluster_threshold
    allele_cluster_threshold <- tmp
  }

  ## convert thresholds to the distance scale used by cutree
  if (distance_method == "decipher") {
    family_threshold_dist         <- 1 - family_threshold / 100
    allele_cluster_threshold_dist <- 1 - allele_cluster_threshold / 100
  } else {
    ## hamming / lv: distance values are already raw counts
    family_threshold_dist         <- family_threshold
    allele_cluster_threshold_dist <- allele_cluster_threshold
  }

  ## cluster the germline
  germline_cluster <- hclust(stats::as.dist(germline_distance), method = cluster_method)

  ## cut the groups based on the threshold
  families_cut <- data.frame(
    "family" = dendextend::cutree(
      as.dendrogram(germline_cluster, hang = -1),
      h = family_threshold_dist,
      order_clusters_as_data = FALSE
    )
  )
  families_cut$iuis_allele <- rownames(families_cut)

  allele_cluster_cut <- data.frame(
    "allele_cluster" = dendextend::cutree(
      as.dendrogram(germline_cluster, hang = -1),
      h = allele_cluster_threshold_dist,
      order_clusters_as_data = FALSE
    )
  )
  allele_cluster_cut$iuis_allele <- rownames(allele_cluster_cut)

  ## merge the tables
  alleleClusterTable <- merge(families_cut, allele_cluster_cut, by = "iuis_allele", all.y = TRUE)
  alleleClusterTable <- alleleClusterTable[, c("iuis_allele", "family", "allele_cluster")]

  list(
    alleleClusterTable = alleleClusterTable,
    threshold = list(
      family_threshold = family_threshold,
      allele_cluster_threshold = allele_cluster_threshold
    ),
    hclustAlleleCluster = germline_cluster,
    communityObject = NULL,
    graphObject = NULL,
    silhouetteScore = NA_real_,
    resolutionParameter = NA_real_,
    families_cut = families_cut
  )
}

## Internal function for Leiden community detection
.igClust_leiden <- function(germline_distance,
                            dist_obj,
                            distance_method = "decipher",
                            family_threshold = 75,
                            cluster_method = "complete",
                            resolution,
                            target_clusters,
                            optimize_silhouette,
                            ncores,
                            quiet) {

  ## convert distance to similarity graph
  g <- distance_to_graph(germline_distance)

  ## determine resolution
  if (is.null(resolution)) {
    if (optimize_silhouette && !is.null(target_clusters)) {
      if (!quiet) message("Optimizing resolution using silhouette score...")
      opt_result <- optimize_resolution(
        g = g,
        distance_matrix = dist_obj,
        target_clusters = target_clusters,
        ncores = ncores
      )
      resolution <- opt_result$best_resolution
      membership <- opt_result$best_partition
      silhouette_score <- opt_result$results$Silhouette[
        opt_result$results$Resolution == resolution
      ]
    } else if (!is.null(target_clusters)) {
      ## use binary search to find resolution
      if (!quiet) message("Finding resolution for target clusters...")
      tuned <- .getNClusters(g, target_clusters)
      resolution <- tuned$best_resolution
      membership <- igraph::membership(tuned$partition)
      silhouette_score <- NA_real_
    } else {
      ## use default resolution
      resolution <- 0.5
      comm <- detect_communities_leiden(g, resolution = resolution)
      membership <- igraph::membership(comm)
      silhouette_score <- NA_real_
    }
  } else {
    ## use provided resolution
    comm <- detect_communities_leiden(g, resolution = resolution)
    membership <- igraph::membership(comm)
    silhouette_score <- NA_real_
  }

  ## build allele cluster table
  allele_names <- names(membership)
  cluster_ids <- as.integer(membership)

  ## compute hierarchical families using the same distance matrix
  if (distance_method == "decipher") {
    family_threshold_dist <- 1 - family_threshold / 100
  } else {
    family_threshold_dist <- family_threshold
  }
  germline_cluster <- hclust(stats::as.dist(germline_distance), method = cluster_method)
  families_cut <- data.frame(
    "family" = dendextend::cutree(
      as.dendrogram(germline_cluster, hang = -1),
      h = family_threshold_dist,
      order_clusters_as_data = TRUE
    )
  )
  families_cut$iuis_allele <- rownames(families_cut)

  ## merge hierarchical families with Leiden allele clusters
  leiden_table <- data.frame(
    iuis_allele = allele_names,
    allele_cluster = cluster_ids,
    stringsAsFactors = FALSE
  )
  alleleClusterTable <- merge(families_cut, leiden_table, by = "iuis_allele", all.y = TRUE)
  alleleClusterTable <- alleleClusterTable[, c("iuis_allele", "family", "allele_cluster")]

  list(
    alleleClusterTable = alleleClusterTable,
    threshold = list(
      family_threshold = family_threshold,
      allele_cluster_threshold = NA_real_,
      resolution = resolution
    ),
    hclustAlleleCluster = germline_cluster,
    communityObject = if (exists("comm")) comm else NULL,
    graphObject = g,
    silhouetteScore = silhouette_score,
    resolutionParameter = resolution,
    families_cut = families_cut
  )
}

#' Allele similarity clustering (deprecated)
#'
#' This function is deprecated. Use \code{\link{igClust}} instead.
#'
#' @inheritParams igClust
#'
#' @return A named list with clustering results.
#'
#' @seealso \code{\link{igClust}} for the current implementation
#'
#' @export
ighvClust <- function(germline_distance,
                      family_threshold = 75,
                      allele_cluster_threshold = 95,
                      cluster_method = "complete") {
  .Deprecated("igClust")
  igClust(
    germline_distance = germline_distance,
    method = "hierarchical",
    family_threshold = family_threshold,
    allele_cluster_threshold = allele_cluster_threshold,
    cluster_method = cluster_method
  )
}

# ------------------------------------------------------------------------------

#' Allele similarity cluster naming scheme
#'
#' For a given cluster the function collapse similar sequences and renames the sequences based on the ASC name scheme
#'
#'
#' @param    cluster               A vector with the cluster identifier - the family and allele cluster number.
#' @param    allele.cluster.table  A data.frame with the list of all germline sequences and their clusters.
#' @param    germ.dist             A matrix with the germline distance between the germline set sequences.
#' @param    chain                 A character with the chain identifier: IGH/IGL/IGK/TRB/TRA... (Currently only IGH is supported)
#' @param    segment               A character with the segment identifier: IGHV/IGHD/IGHJ.... (Currently only IGHV is supported)
#' @param    family_prefix         Logical. If TRUE (default), prepend "F" to the family number in ASC names (e.g. IGHVF1-G1*01). If FALSE, omit the "F" (e.g. IGHV1-G1*01).
#'
#' @return
#'
#' A data.frame with the clusters renamed alleles based on the ASC scheme.
#' @export
alleleClusterNames <- function(cluster,
                               allele.cluster.table,
                               germ.dist,
                               chain,
                               segment,
                               family_prefix = TRUE) {
  # Extract cluster values
  family_cluster <- cluster[[1]]
  allele_cluster <- cluster[[2]]
  f_prefix <- if (family_prefix) "F" else ""
  
  # Filter allele cluster table
  allele.cluster.table <-
    allele.cluster.table[allele.cluster.table$family == family_cluster &
                           allele.cluster.table$allele_cluster == allele_cluster,]
  
  # Check the number of alleles
  if (nrow(allele.cluster.table) == 1) {
    allele.cluster.table$new_allele <-
      paste0(segment, f_prefix, family_cluster, "-G", allele_cluster, "*01")
    allele.cluster.table$removed_duplicated <- FALSE
    return(allele.cluster.table)
  }
  
  # Subset the distance matrix
  germ.dist <-
    germ.dist[allele.cluster.table$iuis_allele, allele.cluster.table$iuis_allele]
  diag(germ.dist) <- NA
  
  # Find similar alleles
  distances <-
    germ.dist[rowSums(germ.dist == 0, na.rm = TRUE) != 0, colSums(germ.dist == 0, na.rm = TRUE) != 0]
  similar <- which(distances == 0, arr.ind = TRUE)
  
  # Check for similar alleles
  if (length(similar) == 0) {
    n <- nrow(allele.cluster.table)
    allele.cluster.table$new_allele <-
      paste0(segment,
             f_prefix,
             family_cluster,
             "-G",
             allele_cluster,
             "*",
             sprintf("%02d", 1:n))
    allele.cluster.table$removed_duplicated <- FALSE
  } else {
    all_alleles <- rownames(germ.dist)
    gene <- strsplit(all_alleles[1], "[*]")[[1]][1]
    
    # Collapse similar alleles
    similar <- similar[order(similar[, 1], similar[, 2]),]
    similar_names <-
      apply(similar, 1, function(row)
        paste(sort(row), collapse = ","))
    unique_similar_names <- unique(similar_names)
    
    # Find alleles to keep and remove
    allele_keep <- character(0)
    allele_remove <- character(0)
    
    for (unique_name in unique_similar_names) {
      similar_indices <- which(similar_names == unique_name)
      alleles_idx <- unique(similar[similar_indices,])
      alleles_idx <- rownames(alleles_idx)
      alleles_gene <- strsplit(alleles_idx, "[*]")[[1]]
      
      if (all(alleles_gene != gene)) {
        allele_remove <- c(allele_remove, alleles_idx[2])
        allele_keep <- c(allele_keep, alleles_idx[1])
      } else {
        allele_remove <- c(allele_remove, alleles_idx[2])
        allele_keep <- c(allele_keep, alleles_idx[1])
      }
    }
    
    remove_allele_list <- setNames(allele_keep, allele_remove)
    keep_alleles_list <- setdiff(all_alleles, allele_remove)
    
    allele.cluster.table$new_allele <- ""
    for (i in seq_along(keep_alleles_list)) {
      a <- keep_alleles_list[i]
      allele.cluster.table$new_allele[allele.cluster.table$iuis_allele == a] <-
        paste0(segment,
               f_prefix,
               family_cluster,
               "-G",
               allele_cluster,
               "*",
               sprintf("%02d", i))
    }

    # Add the new names for the removed alleles
    if (length(remove_allele_list) != 0) {
      for (i in seq_along(remove_allele_list)) {
        a <- names(remove_allele_list)[i]
        a_new <- remove_allele_list[i]
        allele.cluster.table$new_allele[allele.cluster.table$iuis_allele == a] <-
          allele.cluster.table$new_allele[allele.cluster.table$iuis_allele == a_new]
      }
    }

    # Flag duplicated alleles
    allele.cluster.table$removed_duplicated <-
      allele.cluster.table$iuis_allele %in% allele_remove
  }
  
  return(allele.cluster.table)
}

#' Generate allele similarity reference set
#'
#' Generates the allele clusters reference set based on the clustering from \link{ighvClust}. The function collapse
#' similar alleles and assign them into their respective allele clusters and family clusters. See details for naming scheme
#'
#' @param    germline_distance     A germline set distance matrix created by \link{ighvDistance}.
#' @param    germline_set          A character list of the IMGT aligned IGHV allele sequences. See details for curating options.
#' @param    alleleClusterTable    A data.frame of the alleles and their clusters created by \link{ighvClust}.
#' @param    trim_3prime_side      If a 3' position trim is supplied, duplicated sequences will be checked for differential positions past the trim position. Default NULL; NULL will not activate the check. see @details
#' @param    family_prefix         Logical. If TRUE (default), prepend "F" to the family number in ASC names (e.g. IGHVF1-G1*01). If FALSE, omit the "F" (e.g. IGHV1-G1*01).
#'
#' @details
#' Each allele is named by this scheme:
#' IGHVF1-G1*01 - IGH = chain, V = region, F1 = family cluster numbering (the "F" prefix can be omitted by setting family_prefix = FALSE),
#' G1 - allele cluster numbering, and 01 = allele numbering (given by clustering order, no connection to the expression)
#'
#' In case there are alleles that are differentiated in a nucleotide position past the trimming position used for the clustering,
#' then the alleles are separated and are annotated with the differentiating position as so:
#' Say A1*01 and A1*02 are similar up to position 318, and thus collapsed in the clusters to G1*01.
#' Upon checking the sequences past the trim position (318), a differentiating nucleotide was seen in position 319,
#' A1*01 has a G, and A1*02 has a T.
#' Then the alleles will be separated, and the new annotation will be as so:
#' A1*01 = G1*01, and A1*02 = G1*01_G319T.
#' Where the first nucleotide indicate the base, the following number the position, and the last nucleotide the one the base changed into.
#'
#'
#' @return
#' A \code{list} with the re-named germline set, and a table of the allele clusters and thresholds.
#'

generateReferenceSet <-
  function(germline_distance,
           germline_set,
           alleleClusterTable,
           trim_3prime_side = NULL,
           family_prefix = TRUE) {
    # check the parameters
    ### check the class of the distance matrix

    if (!(is.matrix(germline_distance) &
          is.array(germline_distance)))
      stop("The input germline distance was not create with ighvDistance.")

    ### check the class of the germline set
    if (!is.character(germline_set))
      stop("The input germline set is not in a character class.")

    ### apply backward-compat shim before validation
    alleleClusterTable <- .compat_allele_table(alleleClusterTable)

    ### check the class and names in alleleClusterTable
    if (!(is.data.frame(alleleClusterTable) &
          all(
            names(alleleClusterTable) %chin% c("iuis_allele", "family", "allele_cluster")
          )))
      stop("alleleClusterTable does not match the ouput from ighvClust.")
    
    ## get the chain and the segment
    chain <- substr(names(germline_set)[1], 1, 3)
    segment <- substr(names(germline_set)[1], 1, 4)
    
    ## get the allele cluster
    clusters <-
      alleleClusterTable[!duplicated(alleleClusterTable[, c(2, 3)]), c(2, 3)]
    
    ## rename the alleles
    alleleClusterTable.tmp <- apply(
      clusters,
      1,
      alleleClusterNames,
      allele.cluster.table = alleleClusterTable,
      germ.dist = germline_distance,
      chain = chain,
      segment = segment,
      family_prefix = family_prefix
    )
    
    alleleClusterTable.tmp <-
      data.table::rbindlist(alleleClusterTable.tmp, fill = T)
    
    ## check if the alleles that were declared duplicated by the trimming
    ## have a differential position past the trimming.
    ## If so, adjust the name a the duplicated column.
    
    check_alleles <-
      function(alleles,
               germline_set,
               trim_3prime_side) {
        if (length(alleles) > 1) {
          reference_allele <- germline_set[alleles[1]]
          
          diffs <- sapply(alleles[-1], function(allele) {
            dif <-
              allele_diff(reference_allele, germline_set[allele], position_threshold = trim_3prime_side)
            paste0(dif, collapse = "_")
          })
          return(c("",diffs))
        } else{
          return("")
        }
      }
    
    if (!is.null(trim_3prime_side)) {
      alleleClusterTable.tmp[, "diff_pos_past_trim" := check_alleles(get("iuis_allele"), germline_set, trim_3prime_side), by = list(get("new_allele"))]
      
      alleleClusterTable.tmp$new_allele[alleleClusterTable.tmp$diff_pos_past_trim !=
                                          ""] <-
        paste0(
          alleleClusterTable.tmp$new_allele[alleleClusterTable.tmp$diff_pos_past_trim !=
                                              ""],
          "_",
          alleleClusterTable.tmp$diff_pos_past_trim[alleleClusterTable.tmp$diff_pos_past_trim !=
                                                      ""]
        )
      alleleClusterTable.tmp$removed_duplicated[alleleClusterTable.tmp$diff_pos_past_trim !=
                                                  ""] <- F
    }
    
    
    ## rename the germline set
    germline_set.tmp <- germline_set
    names(germline_set.tmp) <-
      sapply(names(germline_set.tmp), function(x) {
        alleleClusterTable.tmp$new_allele[alleleClusterTable.tmp$iuis_allele == x]
      })
    # remove the duplicated alleles from the germline set
    germline_set.tmp <-
      germline_set.tmp[!duplicated(germline_set.tmp)]
    
    alleleClusterTable.tmp$thresh <- 1e-04
    
    return(
      list(
        alleleClusterTable = alleleClusterTable.tmp,
        alleleClusterSet = germline_set.tmp
      )
    )
    
  }

# ------------------------------------------------------------------------------

frw1_primers <- setNames(
  c(
    'GGCCTCAGTGAAGGTCTCCTGCAAG',
    'GTCTGGTCCTACGCTGGTGAAACCC',
    'CTGGGGGGTCCCTGAGACTCTCCTG',
    'CTTCGGAGACCCTGTCCCTCACCTG',
    'CGGGGAGTCTCTGAAGATCTCCTGT',
    'TCGCAGACCCTCTCACTCACCTGTG'
  ),
  c(
    'VH1-FR1',
    'VH2-FR1',
    'VH3-FR1',
    'VH4-FR1',
    'VH5-FR1',
    'VH6-FR1'
  )
)

### the function gaps the sequences based on the germlines
### input: v sequence and reference; output: gapped sequence

togap <- function(vgap, vdj) {
  ##add in vdj gaps
  gapadd <- vdj
  for (i in which(unlist(strsplit(vgap, "", fixed = T)) == ".")) {
    gapadd <-
      paste0(substr(gapadd, 1, i - 1), ".", substr(gapadd, i, nchar(gapadd)))
  }
  return(gapadd)
}

#' FWR1 artificial dataset generator
#'
#' A function to artificially create an IGHV reference set with framework1 (FWR1) primers (see Details).
#'
#' @param    germline_set      A germline set distance matrix created by `ighvDistance`.
#' @param    mask_primer       Logical (TRUE by default). If to mask with Ns the region of the primer from the germline sequence
#' @param    trimm_primer      Logical (FALSE by default). If to trim the region of the primer from the germline sequence. If TRUE then, mask_primer is ignored.
#' @param    quite             Logical (FALSE by default). Do you want to suppress informative messages
#'
#' @details
#'
#' The FRW1 primers used in this function were taken from the BIOMED-2 protocol. For more information on the protocol and primer design go to:
#' van Dongen, J., Langerak, A., Brüggemann, M. et al. Design and standardization of PCR primers and protocols for detection of clonal immunoglobulin and
#' T-cell receptor gene recombinations in suspect lymphoproliferations: Report of the BIOMED-2 Concerted Action BMH4-CT98-3936.
#' Leukemia 17, 2257–2317 (2003). https://doi.org/10.1038/sj.leu.2403202Van Dongen, J. J. M., et al. "Design and standardization of PCR primers and protocols for detection of clonal immunoglobulin and T-cell
#' receptor gene recombinations in suspect lymphoproliferations: report of the BIOMED-2 Concerted Action BMH4-CT98-3936."
#' Leukemia 17.12 (2003): 2257-2317.
#'
#' @return
#' A \code{list} with the input germline set allele and the trimmed/masked sequences.
#'
#' @export
artificialFRW1Germline <-
  function(germline_set,
           mask_primer = TRUE,
           trimm_primer = FALSE,
           quite = FALSE) {
    counter <-
      setNames(rep(0L, length(frw1_primers)), names(frw1_primers))
    frw1_set <- c()
    for (seq_name in names(germline_set)) {
      seq <- germline_set[seq_name]
      
      imgt_family <-
        getFamily(
          seq_name,
          first = F,
          collapse = F,
          strip_d = F,
          omit_nl = F
        )
      
      id_primer <-
        grep(gsub("IGHV", "VH", imgt_family), names(frw1_primers))
      primer <- frw1_primers[id_primer]
      
      if (length(primer) != 0) {
        loc <-
          unlist(aregexec(
            text = gsub("[.]", "", seq),
            pattern = primer,
            max.distance = 4
          ))
        if (length(loc) != 0) {
          counter[names(frw1_primers)[id_primer]] = counter[names(frw1_primers)[id_primer]] + 1
          seq_n_gap <- gsub("[.]", "", seq)
          preceding <-
            substr(seq_n_gap, 1, (loc[1] - 1 + nchar(primer)))
          preceding <- gsub("[AGCT]", "N", preceding)
          frw1_seq <-
            substr(seq_n_gap, (loc[1] + nchar(primer)), nchar(seq_n_gap))
          seq_paste <- paste0(preceding, frw1_seq)
          seq_gapped <- togap(seq, seq_paste)
          frw1_set[seq_name] <-
            if (trimm_primer)
              substr(
                seq_gapped,
                start = aregexec(
                  text = seq_gapped,
                  pattern = "[ATGC]",
                  fixed = F
                )[[1]][1],
                nchar(seq_gapped)
              )
          else
            seq_gapped
        } else{
          frw1_set[seq_name] <- seq
        }
      } else{
        frw1_set[seq_name] <- seq
      }
    }
    
    if (!quite)
      cat(
        paste0(
          sum(counter),
          "/",
          length(germline_set),
          " germline sequences have passed"
        ),
        paste0("Counts by primers: "),
        paste(
          names(counter),
          counter,
          sep = ":",
          collapse = ","
        ),
        sep = "\n"
      )
    
    
    return(frw1_set)
  }

# ------------------------------------------------------------------------------
# @param    set_aligned                    If the germline set provided is aligned, if the set is not aligned an alignment with msa alignment is computed. Default is TRUE
#' Allele similarity cluster
#'
#' A wrapper function to infer the allele clusters. Supports both hierarchical
#' clustering (default) and Leiden community detection.
#'
#' @param germline_set A character vector of Ig sequence alleles (must be gapped by IMGT scheme for optimal results).
#' @param locus The locus type. One of "IGHV", "IGKV", "IGLV", "IGHD", "IGHJ", "IGKJ", "IGLJ".
#'   Default is NULL (auto-detected from sequence names).
#' @param clustering_method Clustering method. One of "hierarchical" (default) or "leiden".
#' @param distance_method Distance calculation method. One of "decipher" (default), "hamming", or "lv".
#' @param trim_3prime_side Position to trim sequences from 3' end. Default is 318; NULL uses full length.
#' @param mask_5prime_side Length to mask from 5' side. Default is 0.
#' @param family_threshold Similarity threshold for family level. Default is 75.
#' @param allele_cluster_threshold Similarity threshold for allele cluster level (hierarchical only). Default is 95.
#' @param cluster_method The linkage method for hierarchical clustering (used for family assignment in both methods). Default is "complete".
#' @param resolution Resolution parameter for Leiden clustering. Default is NULL (auto-optimized).
#' @param target_clusters Target number of clusters for Leiden optimization. Default is NULL.
#' @param optimize_silhouette Optimize resolution using silhouette score (Leiden only). Default is TRUE.
#' @param ncores Number of cores for parallel processing (Leiden only). Default is 1.
#' @param aa_set Logical. Is the sequence set amino acids? Default is FALSE.
#' @param quiet Logical. Suppress messages. Default is FALSE.
#' @param family_prefix Logical. If TRUE (default), prepend "F" to the family number in ASC names (e.g. IGHVF1-G1*01). If FALSE, omit the "F" (e.g. IGHV1-G1*01).
#'
#' @details
#' The distance between pairs of allele sequences is calculated, then the alleles are clustered.
#' For hierarchical clustering, two similarity thresholds define family and allele clusters.
#' For Leiden clustering, community detection identifies clusters at a specified resolution.
#'
#' The allele cluster names follow this scheme:
#' IGHVF1-G1*01 - IGH = chain, V = region, F1 = family cluster numbering,
#' G1 = allele cluster numbering, 01 = allele numbering (by clustering order)
#'
#' For V segments, the "decipher" distance method is recommended.
#' For D and J segments with variable lengths, "lv" (Levenshtein) is more appropriate.
#'
#' @return
#' An object of class \link{GermlineCluster} containing:
#' \itemize{
#'   \item germlineSet: Modified germline set (3' trimming and 5' masking)
#'   \item alleleClusterSet: Renamed germline set with ASC names
#'   \item alleleClusterTable: data.frame of allele similarity clusters
#'   \item threshold: List of threshold parameters
#'   \item hclustAlleleCluster: hclust object (both methods)
#'   \item clusteringMethod: Method used ("hierarchical" or "leiden")
#'   \item communityObject: Community object (Leiden method)
#'   \item graphObject: igraph object (Leiden method)
#'   \item silhouetteScore: Silhouette score (Leiden method)
#'   \item resolutionParameter: Resolution used (Leiden method)
#'   \item locus: Locus identifier
#' }
#'
#' @seealso \code{\link{igDistance}}, \code{\link{igClust}}, \code{\link{plot.GermlineCluster}}
#'
#' @examples
#' # load the initial germline set
#' \donttest{
#' data(HVGERM)
#'
#' germline <- HVGERM[!grepl("^[.]", HVGERM)]
#'
#' # Hierarchical clustering (default)
#' asc <- inferAlleleClusters(germline)
#'
#' # Leiden community detection
#' asc_leiden <- inferAlleleClusters(germline[1:50],
#'                                   clustering_method = "leiden",
#'                                   target_clusters = 10)
#'
#' ## plotting the clusters
#' plot(asc)
#' }
#' @export
inferAlleleClusters <- function(germline_set,
                                locus = NULL,
                                clustering_method = c("hierarchical", "leiden"),
                                distance_method = c("decipher", "hamming", "lv"),
                                trim_3prime_side = 318,
                                mask_5prime_side = 0,
                                family_threshold = 75,
                                allele_cluster_threshold = 95,
                                cluster_method = "complete",
                                resolution = NULL,
                                target_clusters = NULL,
                                optimize_silhouette = TRUE,
                                ncores = 1,
                                aa_set = FALSE,
                                quiet = FALSE,
                                family_prefix = TRUE) {

  clustering_method <- match.arg(clustering_method)
  distance_method <- match.arg(distance_method)

  if (!is.vector(germline_set, mode = "character"))
    stop("The input germline set must be a character vector.")

  if (length(germline_set) < 2)
    stop("At least two sequences are required.")

  ## detect locus from germline names if not provided
  if (is.null(locus)) {
    locus <- substr(names(germline_set)[1], 1, 4)
  }

  ## create a copy of the germline set to return (trimming/masking only for clustering)
  germline_set_copy <- germline_set

  ## apply 3' trimming
  if (!is.null(trim_3prime_side)) {
    germline_set <- substr(germline_set, 1, trim_3prime_side)
  }

  ## apply 5' masking
  if (mask_5prime_side != 0) {
    trailing <- substr(germline_set, mask_5prime_side + 1, stop = 1000000L)
    preceding <- paste0(rep("N", mask_5prime_side), collapse = "")
    germline_set <- paste0(preceding, trailing)
    names(germline_set) <- names(trailing)
  }

  ## calculate distance matrix
  germline_distance <- igDistance(
    germline_set,
    AA = aa_set,
    method = distance_method,
    quiet = quiet
  )

  if (clustering_method == "hierarchical") {
    ## check thresholds
    if (distance_method == "decipher") {
      if (any(c(family_threshold, allele_cluster_threshold) < 0 |
              c(family_threshold, allele_cluster_threshold) > 100))
        stop("Thresholds must be between 0-100.")
    } else {
      if (any(c(family_threshold, allele_cluster_threshold) < 0))
        stop("Thresholds must be non-negative.")
    }

    wrong_order <- if (distance_method == "decipher") {
      family_threshold > allele_cluster_threshold
    } else {
      family_threshold < allele_cluster_threshold
    }
    if (wrong_order) {
      if (!quiet) message("Switching family and allele cluster thresholds.")
      tmp <- family_threshold
      family_threshold <- allele_cluster_threshold
      allele_cluster_threshold <- tmp
    }

    ## hierarchical clustering
    cluster_results <- igClust(
      germline_distance,
      method = "hierarchical",
      distance_method = distance_method,
      family_threshold = family_threshold,
      allele_cluster_threshold = allele_cluster_threshold,
      cluster_method = cluster_method
    )

    cluster_renamed <- generateReferenceSet(
      germline_distance = germline_distance,
      germline_set = germline_set_copy,
      alleleClusterTable = cluster_results$alleleClusterTable,
      trim_3prime_side = trim_3prime_side,
      family_prefix = family_prefix
    )

    results <- new_germline_cluster(
      germlineSet = germline_set,
      alleleClusterSet = cluster_renamed$alleleClusterSet,
      alleleClusterTable = as.data.frame(cluster_renamed$alleleClusterTable),
      threshold = list(
        family_threshold = family_threshold,
        allele_cluster_threshold = allele_cluster_threshold
      ),
      hclustAlleleCluster = cluster_results$hclustAlleleCluster,
      clusteringMethod = "hierarchical",
      distanceMatrix = germline_distance,
      locus = locus,
      .familiesCut = cluster_results$families_cut
    )

  } else {
    ## Leiden community detection
    cluster_results <- igClust(
      germline_distance,
      method = "leiden",
      distance_method = distance_method,
      family_threshold = family_threshold,
      cluster_method = cluster_method,
      resolution = resolution,
      target_clusters = target_clusters,
      optimize_silhouette = optimize_silhouette,
      ncores = ncores,
      quiet = quiet
    )

    cluster_renamed <- generateReferenceSet(
      germline_distance = germline_distance,
      germline_set = germline_set_copy,
      alleleClusterTable = cluster_results$alleleClusterTable,
      trim_3prime_side = trim_3prime_side,
      family_prefix = family_prefix
    )

    results <- new_germline_cluster(
      germlineSet = germline_set,
      alleleClusterSet = cluster_renamed$alleleClusterSet,
      alleleClusterTable = as.data.frame(cluster_renamed$alleleClusterTable),
      threshold = list(
        family_threshold = family_threshold,
        allele_cluster_threshold = NA_real_,
        resolution = cluster_results$resolutionParameter
      ),
      hclustAlleleCluster = cluster_results$hclustAlleleCluster,
      clusteringMethod = "leiden",
      communityObject = cluster_results$communityObject,
      graphObject = cluster_results$graphObject,
      distanceMatrix = germline_distance,
      silhouetteScore = cluster_results$silhouetteScore,
      resolutionParameter = cluster_results$resolutionParameter,
      locus = locus,
      .familiesCut = cluster_results$families_cut
    )
  }

  return(results)
}

# ------------------------------------------------------------------------------

#' Plotting the dendrogram of the clusters
#'
#' @param    x                            The GermlineCluster object. See \link{inferAlleleClusters}
#' @param    y                            NULL. not in use.
#' @param    cex                          Controls the size of the allele label. Default is 1.
#' @param    seed                         Set a seed number for drawing the dendrogram. Default 9999.
#' @return
#' A plot of the allele clusters dendrogram
#'

plotAlleleCluster <- function(x,
                              y = NULL,
                              cex = 1,
                              seed = 9999) {
  ## check the class of the object

  if (!inherits(x, "GermlineCluster"))
    stop("Object is not of class GermlineCluster")

  alleleClusterTable <- x$alleleClusterTable
  hclustAlleleCluster <- x$hclustAlleleCluster

  fam_cut = 1 - x$threshold$family_threshold / 100
  group_cut = 1 - x$threshold$allele_cluster_threshold / 100
  
  ## get the colors
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category %chin% c('qual', 'div'),]
  col_vector = unlist(mapply(
    RColorBrewer::brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
  group_color_vector = setNames(c("#FFFFFF", col_vector), c(0, sort(
    unique(alleleClusterTable$allele_cluster)
  )))
  
  col_fun = circlize::colorRamp2(breaks = c(0:length(
    unique(alleleClusterTable$allele_cluster)
  )),
  colors = group_color_vector[1:(length(unique(alleleClusterTable$allele_cluster)) +
                                   1)])
  
  ## prepare the data
  db_sub <-
    alleleClusterTable[!alleleClusterTable$removed_duplicated,]
  chain <- substr(alleleClusterTable$iuis_allele[1], 1, 3)
  segment <- substr(alleleClusterTable$iuis_allele[1], 1, 4)
  iuis_to_asc <-
    setNames(alleleClusterTable$new_allele,
             alleleClusterTable$iuis_allele)

  group.df <- data.frame(
    A = labels(hclustAlleleCluster),
    Allele = iuis_to_asc[labels(hclustAlleleCluster)],
    fam = sapply(strsplit(iuis_to_asc[labels(hclustAlleleCluster)], "-"), "[[", 1),
    row.names = labels(hclustAlleleCluster)
  )
  
  
  
  group.df$Group <- 0
  for (i in 1:nrow(group.df)) {
    allele <- group.df$A[i]
    full_allele <- db_sub$new_allele[db_sub$iuis_allele == allele]
    
    if (length(full_allele) != 0) {
      # potential groups
      full_groups <-
        as.numeric(strsplit(gsub(
          paste0(segment, "F[0-9]+-G"), "", full_allele
        ), "[*]")[[1]][1])
      
      group.df$Group[i] <- as.numeric(full_groups)
      
    }
  }
  
  group.mat <- as.matrix(group.df$Group)
  colnames(group.mat) <- c("S1")
  rownames(group.mat) <- gsub(segment, "", group.df$A)
  useT.mat <- t(group.mat)
  mat_list = useT.mat
  dend_list = as.dendrogram(hclustAlleleCluster)
  
  
  circos.clear()
  circos.par(
    "start.degree" = 90,
    cell.padding = c(0, 0, 0, 0),
    gap.degree = 10,
    circle.margin = c(0.01, 0.01)
  )
  circos.initialize("a", xlim = c(0, nrow(alleleClusterTable))) #changed to 1 heatmap setting
  #adding new track for column label
  circos.track(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.05,
    panel.fun = function(x, y) {
      for (i in seq_len(ncol(mat_list))) {
        circos.text(
          i - 0.5,
          0,
          colnames(mat_list)[i],
          adj = c(0, 0.5),
          facing = "clockwise",
          niceFacing = F,
          cex = cex
        )
      }
    }
  )
  
  circos.track(
    ylim = c(0, 1),
    bg.border = NA,
    panel.fun = function(x, y) {
      m = mat_list
      dend = dend_list
      #changed variable for 1 heatmap setting
      m2 = m#[, order.dendrogram(dend)]
      col_mat = col_fun(m2)
      nr = nrow(m2)
      nc = ncol(m2)
      for (i in 1:nr) {
        circos.rect(
          1:nc - 1,
          rep(nr - i, nc),
          1:nc,
          rep(nr - i + 1, nc),
          border = col_mat[i, ],
          col = col_mat[i, ]
        )
        
        ## add lines for S2 and S3.
        draw_arrow <-
          function(x0,
                   x1,
                   y0,
                   y1,
                   col = "black",
                   lty = 1,
                   lwd = 1,
                   eps = 0.5) {
            circos.arrow(
              x1 = x0 ,
              #- eps,
              x2 = x1 - eps,
              y = y0,
              arrow.head.width = 0.1,
              width = 0,
              arrow.position = "start",
              arrow.head.length = 0.2,
              #y1 = y1,
              col = col
              #lty = lty,
              #lwd = lwd,
            )
          }
        
        draw_line <-
          function(x0,
                   x1,
                   y0,
                   y1,
                   col = "black",
                   lty = 1,
                   lwd = 1,
                   eps = 0.5) {
            circos.segments(
              x0 = x0 ,
              #- eps,
              x1 = x1,
              y0 = y0,
              y1 = y1,
              straight = T
            )
          }
        
        #if(i!=3){
        tt <- mat_list[i,]
        col_arrow <- c()
        zero_idx_start <- c()
        zero_idx_end <- c()
        end <- 1
        while (end <= length(tt)) {
          if (tt[end] == 0) {
            zero_idx_start <- c(zero_idx_start, end - 1)
            col_arrow <- c(col_arrow, tt[(end - 1)])
            start_i <- min(c(end + 1, length(tt)))
            for (col in start_i:length(tt)) {
              if (tt[col] != 0) {
                zero_idx_end <- c(zero_idx_end, col - 1)
                break
              }
            }
            end <- col
          }
          end <- end + 1
        }
        
        zero_idx_end <- zero_idx_end <- c(zero_idx_end, length(tt))
        for (ii in seq_along(zero_idx_start)) {
          x0 = zero_idx_start[ii]
          x1 = zero_idx_end[ii]
          #if(i==1) print(paste(x0,x1,col_arrow[ii], tt[x0], tt[x0+1]))
          ys <- c(rep(nr - i, nc), rep(nr - i + 1, nc))
          y0 = mean(ys)
          y1 = i
          draw_arrow(x0, x1, y0, y1, col = group_color_vector[as.character(col_arrow[ii])])
          draw_line(x1 - 0.5, x1 - 0.5, min(ys), max(ys))
        }
        
        #}
        
        
        
        #circos.points(1:nc, 0, colnames(useT.mat)[i])
        
        if (i == 1) {
          gs <- m2[i, ]
          gs <- gs[gs != 0]
          labs = unique(paste0("G", gs))
          groups_labs <-
            sapply(unique(gs), function(x)
              which(m2[i,] == x))
          y_pos = mean(c(rep(nr - i, nc), rep(nr - i + 1, nc)))#-ifelse(i==1, 0, 0.5)
          #adj = c(ifelse(i==1, 0, 0), 0.5)
          facing = ifelse(i == 1, "clockwise", "inside")
          for (ii in seq_len(length(groups_labs))) {
            x_pos = median(groups_labs[[ii]])
            if (m2[i, x_pos] == 0)
              x_pos = x_pos + 1
            circos.text(
              x_pos - 0.5,
              y_pos,
              labs[ii],
              #adj = adj,
              facing = facing,
              niceFacing = T,
              cex = 0.9,
              font = 2
            )
          }
        }
      }
      
      #adding row label
      # circos.text(rep(1, nr+1), 1:(nr+1),
      #             c(rev(rownames(useT.mat)),"IUIS"),
      #             facing = "downward", adj = c(3, 1.5), cex = 1)
      
      
    }
  )
  
  max_height = attr(dend_list, "height") #changed for 1 dendrogram setting
  set.seed(seed)
  
  fam_col <-
    c(
      "brown4",
      "darkblue",
      "darkorchid4",
      "darkgreen",
      "firebrick",
      "darkorange3",
      "deeppink4",
      "deepskyblue4",
      "darkslategrey",
      "green",
      "red",
      "blue",
      "purple",
      "orange",
      "pink",
      "gray"
    )
  
  max_fam <- length(unique(group.df$fam))
  fam_col <- c(fam_col, group_color_vector)[1:max_fam]
  
  circos.track(
    ylim = c(0, max_height),
    bg.border = NA,
    track.height = 0.65,
    panel.fun = function(x, y) {
      dend = dend_list
      circos.dendrogram.piglet(
        dend %>% dendextend::color_branches(dend, k = max_fam, col = fam_col),
        max_height = max_height,
        fam_cut = fam_cut,
        group_cut = group_cut
      )
    }
  )
  circos.clear()
  
  plot_circos <- recordPlot()
  
  
}

#' Plot method for GermlineCluster
#'
#' @param x GermlineCluster object
#' @param y Not used
#' @param cex Controls the size of the allele label. Default is 1.
#' @param seed Set a seed number for drawing the dendrogram. Default 9999.
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A plot of the allele clusters dendrogram
#' @export
plot.GermlineCluster <- function(x, y = NULL, cex = 1, seed = 9999, ...) {
  plotAlleleCluster(x, y = y, cex = cex, seed = seed)
}


### adapted from the circlize package.
circos.dendrogram.piglet <-
  function (dend,
            facing = c("outside", "inside"),
            max_height = NULL,
            use_x_attr = FALSE,
            sector.index = get.current.sector.index(),
            track.index = get.current.track.index(),
            lty_line = 1,
            fam_cut = 0.25,
            group_cut = 0.05,
            seed = 123456)
  {
    set.seed(seed)
    os = get.current.sector.index()
    ot = get.current.track.index()
    set.current.cell(sector.index, track.index)
    on.exit(set.current.cell(os, ot))
    facing = match.arg(facing)[1]
    if (is.null(max_height)) {
      max_height = attr(dend, "height")
    }
    is.leaf = function(object) {
      leaf = attr(object, "leaf")
      if (is.null(leaf)) {
        FALSE
      }
      else {
        leaf
      }
    }
    use_x_attr = use_x_attr
    lines_par = function(col = par("col"),
                         lty = par("lty"),
                         lwd = par("lwd"),
                         ...) {
      return(list(col = col, lty = lty, lwd = lwd))
    }
    points_par = function(col = par("col"),
                          pch = par("pch"),
                          cex = par("cex"),
                          ...) {
      return(list(col = col, pch = pch, cex = cex))
    }
    draw.d = function(dend,
                      max_height,
                      facing = "outside",
                      max_width = 0) {
      leaf = attr(dend, "leaf")
      height = attr(dend, "height")
      midpoint = attr(dend, "midpoint")
      n = length(dend)
      xl = numeric(n)
      yl = numeric(n)
      for (i in seq_len(n)) {
        if (use_x_attr) {
          xl[i] = attr(dend[[i]], "x")
        }
        else {
          if (is.leaf(dend[[i]])) {
            xl[i] = x[as.character(attr(dend[[i]], "label"))]
          }
          else {
            xl[i] = attr(dend[[i]], "midpoint") + x[as.character(labels(dend[[i]]))[1]]
          }
        }
        yl[i] = attr(dend[[i]], "height")
      }
      edge_par_lt = vector("list", n)
      for (i in seq_len(n)) {
        edge_par_lt[[i]] = do.call("lines_par", as.list(attr(dend[[i]],
                                                             "edgePar")))
      }
      node_par = attr(dend, "nodePar")
      if (!is.null(node_par))
        node_par = do.call("points_par", as.list(attr(dend,
                                                      "nodePar")))
      if (facing == "outside") {
        if (n == 1) {
          circos.lines(
            c(xl[1], xl[1]),
            max_height - c(yl[1],
                           height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd,
            straight = TRUE
          )
        }
        else {
          circos.lines(
            c(xl[1], xl[1]),
            max_height - c(yl[1],
                           height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd,
            straight = TRUE
          )
          circos.lines(
            c(xl[1], (xl[1] + xl[2]) / 2),
            max_height -
              c(height, height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd
          )
          if (n > 2) {
            for (i in seq(2, n - 1)) {
              circos.lines(
                c(xl[i], xl[i]),
                max_height -
                  c(yl[i], height),
                col = edge_par_lt[[i]]$col,
                lty = edge_par_lt[[i]]$lty,
                lwd = edge_par_lt[[i]]$lwd,
                straight = TRUE
              )
              circos.lines(
                c((xl[i - 1] + xl[i]) / 2, (xl[i] +
                                              xl[i + 1]) / 2),
                max_height - c(height, height),
                col = edge_par_lt[[i]]$col,
                lty = edge_par_lt[[i]]$lty,
                lwd = edge_par_lt[[i]]$lwd
              )
            }
          }
          circos.lines(
            c(xl[n], xl[n]),
            max_height - c(yl[n],
                           height),
            col = edge_par_lt[[n]]$col,
            lty = edge_par_lt[[n]]$lty,
            lwd = edge_par_lt[[n]]$lwd,
            straight = TRUE
          )
          circos.lines(
            c(xl[n], (xl[n] + xl[n - 1]) / 2),
            max_height - c(height, height),
            col = edge_par_lt[[n]]$col,
            lty = edge_par_lt[[n]]$lty,
            lwd = edge_par_lt[[n]]$lwd
          )
        }
        if (!is.null(node_par)) {
          circos.points(
            mean(xl) / 2,
            max_height - height,
            col = node_par$col,
            pch = node_par$pch,
            cex = node_par$cex
          )
        }
      }
      else if (facing == "inside") {
        if (n == 1) {
          circos.lines(
            c(xl[1], xl[1]),
            c(yl[1], height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd,
            straight = TRUE
          )
        }
        else {
          circos.lines(
            c(xl[1], xl[1]),
            c(yl[1], height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd,
            straight = TRUE
          )
          circos.lines(
            c(xl[1], (xl[1] + xl[2]) / 2),
            c(height,
              height),
            col = edge_par_lt[[1]]$col,
            lty = edge_par_lt[[1]]$lty,
            lwd = edge_par_lt[[1]]$lwd
          )
          if (n > 2) {
            for (i in seq(2, n - 1)) {
              circos.lines(
                c(xl[i], xl[i]),
                c(yl[i], height),
                col = edge_par_lt[[i]]$col,
                lty = edge_par_lt[[i]]$lty,
                lwd = edge_par_lt[[i]]$lwd,
                straight = TRUE
              )
              circos.lines(
                c((xl[i - 1] + xl[i]) / 2, (xl[i] +
                                              xl[i + 1]) / 2),
                c(height, height),
                col = edge_par_lt[[i]]$col,
                lty = edge_par_lt[[i]]$lty,
                lwd = edge_par_lt[[i]]$lwd
              )
            }
          }
          circos.lines(
            c(xl[n], xl[n]),
            c(yl[n], height),
            col = edge_par_lt[[n]]$col,
            lty = edge_par_lt[[n]]$lty,
            lwd = edge_par_lt[[n]]$lwd,
            straight = TRUE
          )
          circos.lines(
            c(xl[n], (xl[n] + xl[n - 1]) / 2),
            c(height, height),
            col = edge_par_lt[[n]]$col,
            lty = edge_par_lt[[n]]$lty,
            lwd = edge_par_lt[[n]]$lwd
          )
          
        }
        if (!is.null(node_par)) {
          circos.points(
            mean(xl) / 2,
            height,
            col = node_par$col,
            pch = node_par$pch,
            cex = node_par$cex
          )
        }
      }
      for (i in seq_len(n)) {
        if (is.leaf(dend[[i]])) {
          node_par = attr(dend[[i]], "nodePar")
          if (!is.null(node_par))
            node_par = do.call("points_par", as.list(attr(dend[[i]],
                                                          "nodePar")))
          if (facing == "outside") {
            if (!is.null(node_par)) {
              circos.points(
                xl[i],
                max_height,
                col = node_par$col,
                pch = node_par$pch,
                cex = node_par$cex
              )
            }
          }
          else if (facing == "inside") {
            if (!is.null(node_par)) {
              circos.points(
                xl[i],
                0,
                col = node_par$col,
                pch = node_par$pch,
                cex = node_par$cex
              )
            }
          }
        }
        else {
          draw.d(dend[[i]], max_height, facing, max_width)
        }
      }
    }
    
    dashed_line <- function(x0, x1, y0, y1, lty, lwd = 1, col) {
      idx <- seq(from = x0,
                 to = (x1 - 1),
                 by = 2)
      for (x in idx) {
        circos.segments(
          x0 = x,
          x1 = x + 1,
          y0 = y0,
          y1 = y1,
          col = col,
          lty = lty,
          lwd = lwd,
          straight = F
        )
      }
    }
    
    labels = as.character(labels(dend))
    x = seq_along(labels) - 0.5
    names(x) = labels
    n = length(labels)
    if (!is.leaf(dend)) {
      draw.d(dend, max_height, facing, max_width = n)
      dashed_line(
        x0 = 1,
        x1 = n,
        y0 = max_height - group_cut,
        y1 = max_height - group_cut,
        col =   "deepskyblue4",
        lty = lty_line,
        lwd = 3
      )
      
      dashed_line(
        x0 = 1,
        x1 = n,
        y0 = max_height - fam_cut,
        y1 = max_height - fam_cut,
        col =   "darkorange3",
        lty = lty_line,
        lwd = 3
      )
      
    }
    
    
  }
