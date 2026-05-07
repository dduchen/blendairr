# ------------------------------------------------------------------------------
# Community detection functions for allele clustering
# Provides Leiden community detection as an alternative to hierarchical clustering

#' @include piglet.R
#' @include allele_cluster.R
NULL

# ------------------------------------------------------------------------------

#' Convert distance matrix to weighted graph
#'
#' Converts a distance matrix to a weighted igraph object using a log transform
#' that spreads small distances and produces weights in \[0,1\].
#'
#' @param distance_matrix A distance matrix or dist object
#'
#' @return An igraph object with weighted edges
#'
#' @details
#' The transformation uses a log-based similarity measure:
#' 1. Normalize distances by max distance
#' 2. Apply -log transform to convert to similarity
#' 3. Normalize similarities to \[0,1\] range
#' 4. Create weighted undirected graph
#'
#' @seealso \code{\link{detect_communities_leiden}}, \code{\link{igClust}}
#'
#' @examples
#' data(HVGERM)
#' d <- igDistance(HVGERM[1:10], method = "hamming")
#' g <- distance_to_graph(d)
#'
#' @export
distance_to_graph <- function(distance_matrix) {
  if (inherits(distance_matrix, "dist")) {
    distance_matrix <- as.matrix(distance_matrix)
  }

  stopifnot(is.matrix(distance_matrix), nrow(distance_matrix) == ncol(distance_matrix))

  dmax <- max(distance_matrix)
  if (!is.finite(dmax) || dmax == 0) dmax <- 1

  ## normalize distances
  Dn <- distance_matrix / dmax

  ## log transform to similarity
  S <- -log(Dn)
  S[!is.finite(S)] <- 0

  ## normalize to [0,1]
  smax <- .max_finite(S)
  if (is.finite(smax) && smax > 0) S <- S / smax

  diag(S) <- 0
  S <- as.matrix(S)

  igraph::graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
}

## Safe helper for max that handles non-finite values
.max_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x, na.rm = TRUE)
}

## Safe helper for min that handles non-finite values
.min_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x, na.rm = TRUE)
}

# ------------------------------------------------------------------------------

#' Leiden community detection
#'
#' Performs community detection on a weighted graph using the Leiden algorithm
#' with CPM (Constant Potts Model) objective function.
#'
#' @param g An igraph graph object with weighted edges
#' @param resolution Resolution parameter for Leiden algorithm. Higher values
#'   produce more communities. Default is 1.0.
#'
#' @return An igraph communities object
#'
#' @details
#' The Leiden algorithm is a community detection method that optimizes a quality
#' function (here CPM). It guarantees connected communities and is generally
#' faster than Louvain while producing better quality partitions.
#'
#' @seealso \code{\link{distance_to_graph}}, \code{\link{optimize_resolution}}
#'
#' @examples
#' data(HVGERM)
#' d <- igDistance(HVGERM[1:10], method = "hamming")
#' g <- distance_to_graph(d)
#' comm <- detect_communities_leiden(g, resolution = 0.5)
#'
#' @export
detect_communities_leiden <- function(g, resolution = 1.0) {
  if (!igraph::is_igraph(g))
    stop("Input 'g' must be an igraph graph.")

  if (!"weight" %in% igraph::edge_attr_names(g))
    igraph::E(g)$weight <- 1

  igraph::cluster_leiden(
    g,
    objective_function = "CPM",
    resolution = resolution,
    weights = igraph::E(g)$weight
  )
}

# ------------------------------------------------------------------------------

#' Find resolution for target cluster count
#'
#' Uses binary search to find a resolution parameter that produces approximately
#' the target number of clusters.
#'
#' @param g An igraph graph object with weighted edges
#' @param n_cluster Target number of clusters
#' @param range_min Minimum resolution to search. Default is 0.
#' @param range_max Maximum resolution to search. Default is 6.
#' @param max_steps Maximum number of search iterations. Default is 20.
#' @param method Community detection method: "leiden" or "louvain". Default is "leiden".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{partition}: The community detection result
#'   \item \code{clusters}: Number of clusters found
#'   \item \code{best_resolution}: The resolution parameter used
#' }
#'
#' @keywords internal
.getNClusters <- function(g,
                          n_cluster,
                          range_min = 0,
                          range_max = 6,
                          max_steps = 20,
                          method = "leiden") {

  stopifnot(igraph::is_igraph(g), n_cluster >= 1, range_min <= range_max, max_steps >= 1)

  this_step <- 0
  this_min <- range_min
  this_max <- range_max
  pre_min_cluster <- 0
  pre_max_cluster <- NULL
  min_update <- FALSE
  max_update <- FALSE

  best_partition <- NULL
  closest_clusters <- NULL
  this_resolution <- (this_min + this_max) / 2

  if (!"weight" %in% igraph::edge_attr_names(g))
    igraph::E(g)$weight <- 1

  while (this_step < max_steps) {
    this_step <- this_step + 1

    if (this_step == 1) {
      this_resolution <- (this_min + this_max) / 2
    } else {
      if (max_update && !min_update) {
        this_resolution <- this_min + (this_max - this_min) *
          (n_cluster - pre_min_cluster) / (pre_max_cluster - pre_min_cluster)
      } else if (min_update && !max_update) {
        if (!is.null(pre_max_cluster)) {
          this_resolution <- this_min + (this_max - this_min) *
            (n_cluster - pre_min_cluster) / (pre_max_cluster - pre_min_cluster)
        } else {
          this_resolution <- (this_min + this_max) / 2
        }
      }
    }

    partition <- switch(
      method,
      louvain = igraph::cluster_louvain(g, resolution = this_resolution, weights = igraph::E(g)$weight),
      leiden = igraph::cluster_leiden(g, objective_function = "CPM",
                                      resolution = this_resolution, weights = igraph::E(g)$weight),
      stop("Error: Unsupported method. Choose 'louvain' or 'leiden'.")
    )

    this_clusters <- length(unique(igraph::membership(partition)))

    if (this_clusters > n_cluster) {
      this_max <- this_resolution
      pre_max_cluster <- this_clusters
      min_update <- FALSE
      max_update <- TRUE
    } else if (this_clusters < n_cluster) {
      this_min <- this_resolution
      pre_min_cluster <- this_clusters
      min_update <- TRUE
      max_update <- FALSE
    } else {
      closest_clusters <- this_clusters
      best_partition <- partition
      break
    }

    if (is.null(closest_clusters) || abs(this_clusters - n_cluster) < abs(closest_clusters - n_cluster)) {
      closest_clusters <- this_clusters
      best_partition <- partition
    }
  }

  list(partition = best_partition, clusters = closest_clusters, best_resolution = this_resolution)
}

# ------------------------------------------------------------------------------

#' Optimize resolution parameter using silhouette score
#'
#' Performs a grid search over resolution parameters and selects the one
#' that maximizes the silhouette score.
#'
#' @param g An igraph graph object with weighted edges
#' @param distance_matrix The distance matrix (as dist object) used for silhouette calculation
#' @param target_clusters Target number of clusters for initial tuning. Default is 80.
#' @param resolution_range_low Fractional range below tuned resolution. Default is 0.1.
#' @param resolution_range_high Fractional range above tuned resolution. Default is 0.5.
#' @param max_steps Maximum steps for initial tuning. Default is 20.
#' @param ncores Number of cores for parallel processing. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{results}: data.frame with Resolution, ClusterCount, Silhouette
#'   \item \code{partitions}: list of membership vectors for each resolution
#'   \item \code{best_resolution}: optimal resolution parameter
#'   \item \code{best_partition}: membership vector at optimal resolution
#'   \item \code{best_clusters}: number of clusters at optimal resolution
#' }
#'
#' @seealso \code{\link{detect_communities_leiden}}, \code{\link{igClust}}
#'
#' @export
optimize_resolution <- function(g,
                                distance_matrix,
                                target_clusters = 80,
                                resolution_range_low = 0.1,
                                resolution_range_high = 0.5,
                                max_steps = 20,
                                ncores = 1) {

  stopifnot(igraph::is_igraph(g))

  ## ensure distance_matrix is a dist object for silhouette
  D <- if (inherits(distance_matrix, "dist")) distance_matrix else stats::as.dist(distance_matrix)

  if (!"weight" %in% igraph::edge_attr_names(g))
    igraph::E(g)$weight <- 1

  ## tune a starting resolution to hit about target_clusters using Leiden
  tuned <- .getNClusters(
    g = g,
    n_cluster = target_clusters,
    range_min = 0,
    range_max = 1,
    max_steps = max_steps,
    method = "leiden"
  )
  base_res <- tuned$best_resolution

  ## local resolution sweep around the tuned value
  min_res <- base_res - base_res * resolution_range_low
  max_res <- base_res + base_res * resolution_range_high
  resolution_range <- seq(min_res, max_res, by = 0.05)
  if (!base_res %in% resolution_range)
    resolution_range <- sort(unique(c(resolution_range, base_res)))

  partitions <- vector("list", length(resolution_range))
  names(partitions) <- paste0("res-", resolution_range)

  results <- data.frame(
    Resolution = numeric(),
    ClusterCount = integer(),
    Silhouette = numeric()
  )

  for (resolution in resolution_range) {
    ## run Leiden at this resolution
    part <- igraph::cluster_leiden(
      g,
      objective_function = "CPM",
      resolution = resolution,
      weights = igraph::E(g)$weight
    )

    memb <- unlist(igraph::membership(part))
    n_clusters <- length(unique(memb))
    partitions[[paste0("res-", resolution)]] <- memb

    ## calculate silhouette score
    sil_score <- if (n_clusters > 1) {
      si <- cluster::silhouette(as.integer(memb), D)
      if (!any(is.na(si))) mean(si[, 3]) else NA_real_
    } else {
      NA_real_
    }

    results <- rbind(results, data.frame(
      Resolution = resolution,
      ClusterCount = n_clusters,
      Silhouette = sil_score
    ))
  }

  ## pick the resolution with the highest silhouette
  best_idx <- which.max(replace(results$Silhouette, is.na(results$Silhouette), -Inf))
  best_resolution <- results$Resolution[best_idx]
  best_partition <- partitions[[paste0("res-", best_resolution)]]
  best_clusters <- length(unique(best_partition))

  list(
    results = results,
    partitions = partitions,
    best_resolution = best_resolution,
    best_partition = best_partition,
    best_clusters = best_clusters
  )
}

# ------------------------------------------------------------------------------

#' Compute distance matrix
#'
#' Compute a pairwise distance matrix between sequences using stringdist.
#'
#' @param sequences A named character vector of sequences
#' @param method Distance method: "hamming" or "lv" (Levenshtein). Default is "hamming".
#' @param trim_3prime Optional position to trim sequences from 3' end
#' @param quiet Logical. Suppress messages. Default is TRUE.
#' @param return_type One of "dist" (default) or "matrix"
#'
#' @return A dist object or matrix of pairwise distances
#'
#' @seealso \code{\link{igDistance}} for more distance options
#'
#' @keywords internal
compute_distance <- function(sequences,
                             method = c("hamming", "lv"),
                             trim_3prime = NULL,
                             quiet = TRUE,
                             return_type = c("dist", "matrix")) {

  method <- match.arg(method)
  return_type <- match.arg(return_type)

  stopifnot(is.character(sequences), length(sequences) > 1)

  if (!is.null(trim_3prime)) {
    stopifnot(is.numeric(trim_3prime), trim_3prime >= 1)
    sequences <- substr(sequences, 1, trim_3prime)
  }

  if (method == "hamming") {
    ## pad to common length with N so Hamming is defined
    if (!quiet) message(sprintf("Padding to max length %s.", max(nchar(sequences))))
    sequences <- gsub("\\s", "N", format(sequences, width = max(nchar(sequences))))
  }

  ## check if sequences have names
  useNames <- "names"
  if (is.null(names(sequences)) || any(is.na(names(sequences)))) {
    useNames <- "strings"
  }

  d <- stringdist::stringdistmatrix(sequences, method = method, useNames = useNames)

  if (return_type == "dist") {
    return(stats::as.dist(d))
  }
  as.matrix(d)
}
