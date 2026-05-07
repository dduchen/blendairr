# Test igDistance function and distance methods

test_that("igDistance works with decipher method", {
  skip_on_cran()

  # Create simple test sequences
  seqs <- c(
    "IGHV1*01" = "ACGTTGCAACGT",
    "IGHV1*02" = "ACGTTGCAACGA",
    "IGHV2*01" = "TTGCAAATTTGG"
  )

  d <- igDistance(seqs, method = "decipher")

  expect_true(is.matrix(d))
  expect_equal(nrow(d), 3)
  expect_equal(ncol(d), 3)
  expect_equal(rownames(d), names(seqs))
  expect_equal(d[1,1], 0)  # Self-distance is 0
})

test_that("igDistance works with hamming method", {
  seqs <- c(
    A = "ACGTTGCA",
    B = "ACGTTGCG",
    C = "TTGCAAAT"
  )

  d <- igDistance(seqs, method = "hamming")

  expect_true(is.matrix(d))
  expect_equal(nrow(d), 3)
  expect_equal(d["A", "A"], 0)
  expect_equal(d["A", "B"], 1)  # One mismatch
})

test_that("igDistance works with levenshtein method", {
  seqs <- c(
    A = "ACGT",
    B = "ACGTT",   # One insertion
    C = "TTGC"
  )

  d <- igDistance(seqs, method = "lv")

  expect_true(is.matrix(d))
  expect_equal(d["A", "A"], 0)
  expect_equal(d["A", "B"], 1)  # One edit
})

test_that("igDistance handles trim_3prime parameter", {
  seqs <- c(
    A = "ACGTTGCAACGT",
    B = "ACGTTGCAACGA"
  )

  # Full length - 1 difference
  d_full <- igDistance(seqs, method = "hamming")

  # Trimmed - should have 0 difference if we trim before position 12
  d_trim <- igDistance(seqs, method = "hamming", trim_3prime = 11)

  expect_true(d_full["A", "B"] > 0)
  expect_equal(d_trim["A", "B"], 0)
})

test_that("igDistance returns dist object when requested", {
  seqs <- c(A = "ACGT", B = "ACGA", C = "TTGC")

  d <- igDistance(seqs, method = "hamming", return_type = "dist")

  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 3)
})

test_that("igDistance errors on invalid input", {
  expect_error(igDistance(123))  # Not character
  expect_error(igDistance(c(A = "ACGT")))  # Only one sequence
})

test_that("ighvDistance is deprecated but works", {
  seqs <- c(A = "ACGT", B = "ACGA")

  expect_warning(
    d <- ighvDistance(seqs),
    "deprecated"
  )

  expect_true(is.matrix(d))
})
