# ============================================================================
# Load the function to be tested
# ============================================================================

# Load pca_pSVD function from R/ directory
source(testthat::test_path("..", "..", "R", "pca_pSVD.R"))

# ============================================================================
# Cross-platform compatibility helpers
# ============================================================================

# Helper function to handle sign ambiguity in PCA
fix_sign <- function(A, B) {
  s <- sign(colSums(A * B))
  s[s == 0] <- 1
  sweep(A, 2, s, `*`)
}

# ============================================================================
# Basic functionality tests
# ============================================================================

test_that("pca_pSVD basic structure & dimensions", {
  set.seed(1)
  X <- matrix(rnorm(3000), 100, 30)
  out <- pca_pSVD(X, ncomp = 5, center = TRUE, scale = TRUE)

  expect_s3_class(out, "pca")
  expect_true(all(c("sdev", "loadings", "variates", "rotation", "x", "prop_expl_var", "cum.var", "var.tot") %in% names(out)))
  expect_equal(ncol(out$rotation), 5)
  expect_equal(ncol(out$x), 5)
  expect_equal(length(out$sdev), 5)

  # x ≈ Xs %*% rotation
  Xs <- scale(X, center = out$center, scale = out$scale)
  x2 <- Xs %*% out$rotation
  rownames(x2) <- rownames(out$x)
  expect_equal(out$x, x2, tolerance = 1e-7)
})

# ============================================================================
# Comparison with base R prcomp
# ============================================================================

test_that("pca_pSVD vs prcomp: sdev & variance explained close", {
  set.seed(42)
  X <- matrix(rnorm(5000), 100, 50)
  k <- 6
  out <- pca_pSVD(X, ncomp = k, center = TRUE, scale = TRUE)

  pr <- prcomp(X, center = TRUE, scale. = TRUE, retx = TRUE)

  # Compare singular values (sign-invariant)
  expect_equal(out$sdev, pr$sdev[1:k], tolerance = 1e-5)

  # Variance proportion
  var.tot <- sum(apply(scale(X, center = TRUE, scale = TRUE), 2, var))
  expect_equal(out$prop_expl_var$X, (out$sdev^2) / var.tot, tolerance = 1e-8)
})

# ============================================================================
# Comparison with mixOmics::pca (with sign correction)
# ============================================================================

test_that("pca_pSVD vs mixOmics::pca: sign-invariant loadings close", {
  skip_if_not_installed("mixOmics")
  set.seed(7)
  X <- matrix(rnorm(4000), 200, 20)
  k <- 5

  out_mx <- mixOmics::pca(X, ncomp = k, center = TRUE, scale = TRUE)
  out_psv <- pca_pSVD(X, ncomp = k, center = TRUE, scale = TRUE)

  # Variance proportion (sign-invariant)
  expect_equal(unname(out_psv$prop_expl_var$X), unname(out_mx$prop_expl_var$X[1:k]), tolerance = 1e-4)

  # Loading vectors: correct for sign flips before comparison
  L1 <- out_mx$loadings$X[, 1:k, drop = FALSE]
  L2 <- out_psv$loadings$X[, 1:k, drop = FALSE]

  # Verify both have same dimensions (mixOmics loadings don't have rownames)
  expect_equal(dim(L1), dim(L2))
  expect_null(rownames(L1))
  expect_null(rownames(L2))

  # Apply sign correction and compare directly (no subsetting needed)
  L2 <- fix_sign(L2, L1)
  expect_equal(L2, L1, tolerance = 1e-3)
})

# ============================================================================
# CLR transformation test (treat CLR as "none")
# ============================================================================

test_that("CLR logratio is treated as 'none'", {
  set.seed(99)
  X <- matrix(runif(2000, 0.1, 10), 100, 20) # Positive data for CLR
  k <- 5

  # CLR should be treated same as "none"
  out_clr <- pca_pSVD(X, ncomp = k, logratio = "CLR", center = TRUE, scale = TRUE)
  out_none <- pca_pSVD(X, ncomp = k, logratio = "none", center = TRUE, scale = TRUE)

  # Results should be identical
  expect_equal(out_clr$sdev, out_none$sdev, tolerance = 1e-10)
  expect_equal(out_clr$prop_expl_var$X, out_none$prop_expl_var$X, tolerance = 1e-10)
})

# ============================================================================
# Fallback mechanisms
# ============================================================================

test_that("fallback to svd works when svds unavailable", {
  set.seed(123)
  X <- matrix(rnorm(2000), 100, 20)
  # Use ncomp close to min(n,p) to potentially trigger full svd fallback
  k <- min(dim(X)) - 1
  out <- pca_pSVD(X, ncomp = k, center = TRUE, scale = FALSE)
  expect_s3_class(out, "pca")
  expect_equal(ncol(out$rotation), k)
})

test_that("missing values delegate to mixOmics (if installed)", {
  set.seed(1)
  X <- matrix(rnorm(1000), 100, 10)
  X[1, 1] <- NA
  if (requireNamespace("mixOmics", quietly = TRUE)) {
    out <- pca_pSVD(X, ncomp = 2)
    expect_s3_class(out, "pca")
  } else {
    expect_error(pca_pSVD(X, ncomp = 2), "mixOmics not installed")
  }
})

test_that("ILR transformation delegates to mixOmics (if installed)", {
  skip_if_not_installed("mixOmics")
  set.seed(2)
  X <- matrix(runif(1000, 0.1, 10), 100, 10)

  # ILR should delegate to mixOmics
  out <- pca_pSVD(X, ncomp = 3, logratio = "ILR")
  expect_s3_class(out, "pca")
})

# ============================================================================
# Cross-platform numerical stability tests
# ============================================================================

test_that("results are numerically stable across platforms", {
  set.seed(12345)
  X <- matrix(rnorm(3000, mean = 100, sd = 15), 150, 20)
  k <- 5

  out <- pca_pSVD(X, ncomp = k, center = TRUE, scale = TRUE)

  # Test basic properties that should hold on any platform
  # 1. Eigenvalues are non-negative and decreasing
  expect_true(all(out$sdev >= 0))
  expect_true(all(diff(out$sdev) <= 0))

  # 2. Loadings are orthonormal (within tolerance)
  loadings <- out$loadings$X
  gram <- t(loadings) %*% loadings
  expect_equal(unname(gram), diag(k), tolerance = 1e-6)

  # 3. Variance proportions sum to <= 1
  expect_true(sum(out$prop_expl_var$X) <= 1 + 1e-10)

  # 4. Cumulative variance is monotonically increasing
  expect_true(all(diff(out$cum.var) >= 0))
})

test_that("handles edge cases consistently", {
  set.seed(456)

  # Small matrix
  X_small <- matrix(rnorm(100), 10, 10)
  out_small <- pca_pSVD(X_small, ncomp = 3, center = TRUE, scale = TRUE)
  expect_s3_class(out_small, "pca")
  expect_equal(ncol(out_small$rotation), 3)

  # Wide matrix (n < p)
  X_wide <- matrix(rnorm(500), 10, 50)
  out_wide <- pca_pSVD(X_wide, ncomp = 5, center = TRUE, scale = TRUE)
  expect_s3_class(out_wide, "pca")

  # Tall matrix (n > p)
  X_tall <- matrix(rnorm(500), 50, 10)
  out_tall <- pca_pSVD(X_tall, ncomp = 5, center = TRUE, scale = TRUE)
  expect_s3_class(out_tall, "pca")
})
