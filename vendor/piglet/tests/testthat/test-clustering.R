# Test igClust function and clustering methods

test_that("igClust works with hierarchical method", {
  # Create distance matrix
  d <- matrix(c(0, 0.1, 0.3, 0.1, 0, 0.25, 0.3, 0.25, 0), nrow = 3)
  rownames(d) <- colnames(d) <- c("A", "B", "C")

  result <- igClust(d, method = "hierarchical",
                    family_threshold = 75, allele_cluster_threshold = 95)

  expect_type(result, "list")
  expect_true("alleleClusterTable" %in% names(result))
  expect_true("hclustAlleleCluster" %in% names(result))
  expect_true(is.data.frame(result$alleleClusterTable))
  expect_equal(nrow(result$alleleClusterTable), 3)
})

test_that("igClust validates thresholds", {
  d <- matrix(c(0, 0.1, 0.1, 0), nrow = 2)
  rownames(d) <- colnames(d) <- c("A", "B")

  # Should error on invalid thresholds
  expect_error(igClust(d, method = "hierarchical", family_threshold = 150))
  expect_error(igClust(d, method = "hierarchical", family_threshold = -10))
})

test_that("igClust swaps thresholds when family > cluster", {
  d <- matrix(c(0, 0.1, 0.1, 0), nrow = 2)
  rownames(d) <- colnames(d) <- c("A", "B")

  # Should swap and produce message
  expect_message(
    result <- igClust(d, method = "hierarchical",
                      family_threshold = 95, allele_cluster_threshold = 75),
    "Switching"
  )
})

test_that("igClust works with leiden method", {
  skip_on_cran()

  # Create test distance matrix
  seqs <- c(
    A1 = "ACGTTGCAACGT",
    A2 = "ACGTTGCAACGA",
    B1 = "TTGCAAATTTGG",
    B2 = "TTGCAAATTTGA"
  )

  d <- igDistance(seqs, method = "hamming")

  result <- igClust(d, method = "leiden", resolution = 0.5)

  expect_type(result, "list")
  expect_true("alleleClusterTable" %in% names(result))
  expect_true("graphObject" %in% names(result))
  expect_equal(nrow(result$alleleClusterTable), 4)
})

test_that("igClust leiden produces hclust object and stores family_threshold", {
  skip_on_cran()

  seqs <- c(
    A1 = "ACGTTGCAACGT",
    A2 = "ACGTTGCAACGA",
    B1 = "TTGCAAATTTGG",
    B2 = "TTGCAAATTTGA"
  )

  d <- igDistance(seqs, method = "hamming")

  result <- igClust(d, method = "leiden", resolution = 0.5, family_threshold = 75)

  # hclust object should be populated (not NULL)
  expect_true(!is.null(result$hclustAlleleCluster))
  expect_s3_class(result$hclustAlleleCluster, "hclust")

  # family_threshold should be stored in threshold

  expect_equal(result$threshold$family_threshold, 75)

  # family and allele_cluster columns should both exist
  expect_true(all(c("family", "allele_cluster") %in% names(result$alleleClusterTable)))
})

test_that("igClust accepts dist object", {
  d_mat <- matrix(c(0, 0.1, 0.3, 0.1, 0, 0.25, 0.3, 0.25, 0), nrow = 3)
  rownames(d_mat) <- colnames(d_mat) <- c("A", "B", "C")
  d_dist <- as.dist(d_mat)

  result <- igClust(d_dist, method = "hierarchical")

  expect_true(is.data.frame(result$alleleClusterTable))
})

test_that("ighvClust is deprecated but works", {
  d <- matrix(c(0, 0.1, 0.1, 0), nrow = 2)
  rownames(d) <- colnames(d) <- c("A", "B")

  expect_warning(
    result <- ighvClust(d),
    "deprecated"
  )

  expect_true(is.data.frame(result$alleleClusterTable))
})

test_that("alleleClusterNames respects family_prefix = FALSE", {
  skip_on_cran()
  seqs <- c(
    "IGHV1*01" = "ACGTTGCAACGTACGTTGCAACGT",
    "IGHV1*02" = "ACGTTGCAACGAACGTTGCAACGA"
  )
  result_with <- inferAlleleClusters(seqs, trim_3prime_side = NULL, family_prefix = TRUE)
  result_without <- inferAlleleClusters(seqs, trim_3prime_side = NULL, family_prefix = FALSE)

  # With prefix: names contain "F"
  expect_true(any(grepl("F", result_with$alleleClusterTable$new_allele)))
  # Without prefix: no "F" immediately after the segment prefix
  expect_false(any(grepl("IGHVF", result_without$alleleClusterTable$new_allele)))
})
