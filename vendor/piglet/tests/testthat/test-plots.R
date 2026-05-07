# Test visualization functions

test_that("plotCommunityNetwork doesn't error", {
  skip_on_cran()

  seqs <- c(
    A1 = "ACGTTGCAACGT",
    A2 = "ACGTTGCAACGA",
    B1 = "TTGCAAATTTGG"
  )

  d <- igDistance(seqs, method = "hamming")
  g <- distance_to_graph(d)

  obj <- new_germline_cluster(
    germlineSet = seqs,
    alleleClusterSet = seqs,
    alleleClusterTable = data.frame(
      iuis_allele = names(seqs),
      family =c(1, 1, 2),
      allele_cluster =c(1, 1, 2)
    ),
    threshold = list(family_threshold = NA, allele_cluster_threshold = NA),
    clusteringMethod = "leiden",
    graphObject = g
  )

  expect_error(plotCommunityNetwork(obj), NA)
})

test_that("plotCommunityNetwork validates input", {
  expect_error(plotCommunityNetwork("not an object"))

  # Object without graph
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c(A = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list()
  )

  expect_error(plotCommunityNetwork(obj), "graph object")
})

test_that("plotSilhouetteOptimization works", {
  opt_result <- list(
    results = data.frame(
      Resolution = c(0.1, 0.2, 0.3),
      ClusterCount = c(2, 3, 4),
      Silhouette = c(0.5, 0.7, 0.6)
    ),
    best_resolution = 0.2
  )

  expect_error(plotSilhouetteOptimization(opt_result), NA)
})

test_that("plotSilhouetteOptimization validates input", {
  expect_error(plotSilhouetteOptimization("not a list"))
  expect_error(plotSilhouetteOptimization(list(foo = 1)))
})

test_that("plotClusterComparison validates inputs", {
  obj1 <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c(A = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list(),
    clusteringMethod = "hierarchical"
  )

  obj2 <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c(A = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list(),
    clusteringMethod = "leiden"
  )

  # Should work with correct inputs
  expect_error(plotClusterComparison(obj1, obj2), NA)

  # Should error if both hierarchical
  expect_error(plotClusterComparison(obj1, obj1))
})
