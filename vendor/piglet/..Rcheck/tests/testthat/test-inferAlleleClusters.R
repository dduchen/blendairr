# Test inferAlleleClusters main function

test_that("inferAlleleClusters works with hierarchical method", {
  skip_on_cran()

  # Simple test sequences
  seqs <- c(
    "IGHV1*01" = "ACGTTGCAACGTACGTTGCAACGT",
    "IGHV1*02" = "ACGTTGCAACGAACGTTGCAACGA",
    "IGHV2*01" = "TTGCAAATTTGGTTGCAAATTTGG",
    "IGHV2*02" = "TTGCAAATTTGATTGCAAATTTGA"
  )

  result <- inferAlleleClusters(
    seqs,
    clustering_method = "hierarchical",
    trim_3prime_side = NULL,
    family_threshold = 50,
    allele_cluster_threshold = 90
  )

  expect_s3_class(result, "GermlineCluster")
  expect_equal(result$clusteringMethod, "hierarchical")
  expect_equal(nrow(result$alleleClusterTable), 4)
  expect_true(!is.null(result$hclustAlleleCluster))
})

test_that("inferAlleleClusters works with leiden method", {
  skip_on_cran()

  seqs <- c(
    "IGHV1*01" = "ACGTTGCAACGTACGTTGCAACGT",
    "IGHV1*02" = "ACGTTGCAACGAACGTTGCAACGA",
    "IGHV2*01" = "TTGCAAATTTGGTTGCAAATTTGG",
    "IGHV2*02" = "TTGCAAATTTGATTGCAAATTTGA"
  )

  result <- inferAlleleClusters(
    seqs,
    clustering_method = "leiden",
    distance_method = "hamming",
    trim_3prime_side = NULL,
    resolution = 0.5
  )

  expect_s3_class(result, "GermlineCluster")
  expect_equal(result$clusteringMethod, "leiden")
  expect_true(!is.null(result$graphObject))
})

test_that("inferAlleleClusters auto-detects locus", {
  seqs <- c(
    "IGHV1*01" = "ACGT",
    "IGHV1*02" = "ACGA"
  )

  result <- inferAlleleClusters(seqs, trim_3prime_side = NULL)

  expect_equal(result$locus, "IGHV")
})

test_that("inferAlleleClusters handles locus parameter", {
  seqs <- c(
    "IGHV1*01" = "ACGT",
    "IGHV1*02" = "ACGA"
  )

  result <- inferAlleleClusters(seqs, locus = "IGKV", trim_3prime_side = NULL)

  expect_equal(result$locus, "IGKV")
})

test_that("inferAlleleClusters validates input", {
  expect_error(inferAlleleClusters(123))  # Not character
  expect_error(inferAlleleClusters(c(A = "ACGT")))  # Only one sequence
})

test_that("inferAlleleClusters applies trimming", {
  seqs <- c(
    "IGHV1*01" = "ACGTACGTACGTACGT",  # 16 bp
    "IGHV1*02" = "ACGTACGTACGTACGA"
  )

  # With trimming, sequences are truncated before clustering
  result <- inferAlleleClusters(seqs, trim_3prime_side = 10)

  # The stored germlineSet should be trimmed
  expect_true(all(nchar(result$germlineSet) == 10))
})

test_that("inferAlleleClusters applies masking", {
  seqs <- c(
    "IGHV1*01" = "ACGTACGT",
    "IGHV1*02" = "ACGTACGA"
  )

  result <- inferAlleleClusters(seqs, mask_5prime_side = 2, trim_3prime_side = NULL)

  # First 2 positions should be masked with N
  expect_true(all(substr(result$germlineSet, 1, 2) == "NN"))
})

test_that("inferAlleleClusters handles different distance methods", {
  skip_on_cran()

  seqs <- c(
    "IGHV1*01" = "ACGTACGT",
    "IGHV1*02" = "ACGTACGA",
    "IGHV2*01" = "TTGCAAAT"
  )

  result_decipher <- inferAlleleClusters(seqs, distance_method = "decipher",
                                         trim_3prime_side = NULL)
  result_hamming <- inferAlleleClusters(seqs, distance_method = "hamming",
                                        trim_3prime_side = NULL)
  result_lv <- inferAlleleClusters(seqs, distance_method = "lv",
                                   trim_3prime_side = NULL)

  expect_s3_class(result_decipher, "GermlineCluster")
  expect_s3_class(result_hamming, "GermlineCluster")
  expect_s3_class(result_lv, "GermlineCluster")
})
