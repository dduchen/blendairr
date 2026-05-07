# Test community detection functions

test_that("distance_to_graph creates valid igraph", {
  d <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), nrow = 3)
  rownames(d) <- colnames(d) <- c("A", "B", "C")

  g <- distance_to_graph(d)

  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), 3)
  expect_true("weight" %in% igraph::edge_attr_names(g))
})

test_that("distance_to_graph handles dist object", {
  d_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), nrow = 3)
  rownames(d_mat) <- colnames(d_mat) <- c("A", "B", "C")
  d_dist <- as.dist(d_mat)

  g <- distance_to_graph(d_dist)

  expect_true(igraph::is_igraph(g))
})

test_that("detect_communities_leiden returns communities object", {
  d <- matrix(c(0, 0.1, 0.5, 0.1, 0, 0.5, 0.5, 0.5, 0), nrow = 3)
  rownames(d) <- colnames(d) <- c("A", "B", "C")

  g <- distance_to_graph(d)
  comm <- detect_communities_leiden(g, resolution = 0.5)

  expect_s3_class(comm, "communities")
  expect_true(length(igraph::membership(comm)) == 3)
})

test_that("detect_communities_leiden requires igraph input", {
  expect_error(detect_communities_leiden("not a graph"))
})

test_that("optimize_resolution finds clusters", {
  skip_on_cran()

  # Create test data with clear cluster structure
  seqs <- c(
    A1 = "ACGTTGCAACGT",
    A2 = "ACGTTGCAACGA",
    A3 = "ACGTTGCAACGG",
    B1 = "TTGCAAATTTGG",
    B2 = "TTGCAAATTTGA"
  )

  d <- igDistance(seqs, method = "hamming", return_type = "dist")
  g <- distance_to_graph(as.matrix(d))

  opt <- optimize_resolution(g, d, target_clusters = 2, ncores = 1)

  expect_type(opt, "list")
  expect_true("best_resolution" %in% names(opt))
  expect_true("best_partition" %in% names(opt))
  expect_true("results" %in% names(opt))
  expect_true(is.data.frame(opt$results))
})

test_that("compute_distance works with both methods", {
  seqs <- c(A = "ACGT", B = "ACGA", C = "TTGC")

  d_ham <- compute_distance(seqs, method = "hamming", return_type = "matrix")
  d_lv <- compute_distance(seqs, method = "lv", return_type = "matrix")

  expect_true(is.matrix(d_ham))
  expect_true(is.matrix(d_lv))
  expect_equal(d_ham["A", "A"], 0)
  expect_equal(d_lv["A", "A"], 0)
})
