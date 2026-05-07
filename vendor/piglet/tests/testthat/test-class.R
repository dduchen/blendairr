# Test GermlineCluster S3 class

test_that("new_germline_cluster creates valid object", {
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c("G1*01" = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list(family_threshold = 75, allele_cluster_threshold = 95)
  )

  expect_s3_class(obj, "GermlineCluster")
  expect_true("germlineSet" %in% names(obj))
  expect_true("alleleClusterTable" %in% names(obj))
  expect_true("threshold" %in% names(obj))
})

test_that("GermlineCluster has default values for new slots", {
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c("G1*01" = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list(family_threshold = 75, allele_cluster_threshold = 95)
  )

  expect_null(obj$hclustAlleleCluster)
  expect_equal(obj$clusteringMethod, "hierarchical")
  expect_true(is.na(obj$silhouetteScore))
  expect_equal(obj$locus, "IGHV")
})

test_that("print.GermlineCluster works", {
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT", B = "ACGA"),
    alleleClusterSet = c("G1*01" = "ACGT", "G1*02" = "ACGA"),
    alleleClusterTable = data.frame(
      iuis_allele = c("A", "B"),
      family =c(1, 1),
      allele_cluster =c(1, 1)
    ),
    threshold = list(family_threshold = 75, allele_cluster_threshold = 95)
  )

  output <- capture.output(print(obj))

  expect_true(any(grepl("GermlineCluster Object", output)))
  expect_true(any(grepl("Number of alleles", output)))
})

test_that("summary.GermlineCluster returns list", {
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT", B = "ACGA"),
    alleleClusterSet = c("G1*01" = "ACGT", "G1*02" = "ACGA"),
    alleleClusterTable = data.frame(
      iuis_allele = c("A", "B"),
      family =c(1, 1),
      allele_cluster =c(1, 1)
    ),
    threshold = list(family_threshold = 75, allele_cluster_threshold = 95)
  )

  summ <- summary(obj)

  expect_type(summ, "list")
  expect_true("n_alleles" %in% names(summ))
  expect_equal(summ$n_alleles, 2)
})

test_that("plot.GermlineCluster exists", {
  expect_true(is.function(plot.GermlineCluster))
})

test_that("GermlineCluster slots can be accessed with $", {
  obj <- new_germline_cluster(
    germlineSet = c(A = "ACGT"),
    alleleClusterSet = c("G1*01" = "ACGT"),
    alleleClusterTable = data.frame(iuis_allele = "A", family =1, allele_cluster =1),
    threshold = list(family_threshold = 75, allele_cluster_threshold = 95),
    locus = "IGHV"
  )

  expect_equal(length(obj$germlineSet), 1)
  expect_equal(obj$threshold$family_threshold, 75)
  expect_equal(obj$locus, "IGHV")
})
