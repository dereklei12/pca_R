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

  # Check multiple inheritance: class is c("pca_pSVD", "pca")
  expect_s3_class(out, "pca_pSVD")
  expect_s3_class(out, "pca")
  expect_equal(class(out), c("pca_pSVD", "pca"))

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

  # Verify both have same dimensions
  expect_equal(dim(L1), dim(L2))
  expect_null(rownames(L1))
  expect_null(rownames(L2))

  # Apply sign correction and compare directly
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
  k <- 5

  # Mock requireNamespace by unlocking and temporarily replacing it
  env <- asNamespace("base")
  old_requireNamespace <- base::requireNamespace

  # Unlock the binding, replace it, then lock it again on exit
  unlockBinding("requireNamespace", env)
  on.exit(
    {
      assign("requireNamespace", old_requireNamespace, envir = env)
      lockBinding("requireNamespace", env)
    },
    add = FALSE
  )

  assign("requireNamespace", function(pkg, quietly = FALSE) {
    if (pkg %in% c("RSpectra", "rARPACK")) {
      return(FALSE)
    }
    old_requireNamespace(pkg, quietly = quietly)
  }, envir = env)

  out <- expect_warning(
    pca_pSVD(X, ncomp = k, center = TRUE, scale = FALSE),
    "falling back to base::svd"
  )
  expect_s3_class(out, "pca_pSVD")
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
  expect_s3_class(out_small, "pca_pSVD")
  expect_s3_class(out_small, "pca")
  expect_equal(ncol(out_small$rotation), 3)

  # Wide matrix (n < p)
  X_wide <- matrix(rnorm(500), 10, 50)
  out_wide <- pca_pSVD(X_wide, ncomp = 5, center = TRUE, scale = TRUE)
  expect_s3_class(out_wide, "pca_pSVD")
  expect_s3_class(out_wide, "pca")

  # Tall matrix (n > p)
  X_tall <- matrix(rnorm(500), 50, 10)
  out_tall <- pca_pSVD(X_tall, ncomp = 5, center = TRUE, scale = TRUE)
  expect_s3_class(out_tall, "pca_pSVD")
  expect_s3_class(out_tall, "pca")
})

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("input validation catches invalid inputs", {
  # Non-numeric matrix
  expect_error(pca_pSVD(matrix(letters[1:20], 4, 5)), "numeric matrix")

  # Too few rows
  expect_error(pca_pSVD(matrix(1:5, 1, 5)), ">= 2 rows")

  # Too few columns
  expect_error(pca_pSVD(matrix(1:5, 5, 1)), ">= 2 columns")

  # Invalid ncomp
  expect_error(pca_pSVD(matrix(rnorm(100), 10, 10), ncomp = 0), "ncomp must be >= 1")
  expect_error(pca_pSVD(matrix(rnorm(100), 10, 10), ncomp = -1), "ncomp must be >= 1")

  # ncomp too large
  X <- matrix(rnorm(100), 10, 10)
  expect_warning(pca_pSVD(X, ncomp = 15), "ncomp truncated")

  # Missing rownames/colnames are auto-generated
  X <- matrix(rnorm(100), 10, 10)
  out <- pca_pSVD(X, ncomp = 3)
  expect_match(rownames(out$x)[1], "^sample_")
  expect_match(out$names$X[1], "^V")

  # Explicit data.frame input
  df <- data.frame(matrix(rnorm(100), 10, 10))
  out <- pca_pSVD(df, ncomp = 3)
  expect_s3_class(out, "pca_pSVD")
})

# ============================================================================
# Sparse Matrix Support Tests
# ============================================================================

test_that("sparse matrix handling works correctly", {
  skip_if_not_installed("Matrix")

  set.seed(42)
  X_dense <- matrix(rnorm(2000), 100, 20)
  X_dense[abs(X_dense) < 1] <- 0 # Make sparse
  X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)

  # dgCMatrix
  expect_true(inherits(X_sparse, "dgCMatrix"))
  out_sparse <- pca_pSVD(X_sparse, ncomp = 5, center = TRUE, scale = TRUE)
  out_dense <- pca_pSVD(X_dense, ncomp = 5, center = TRUE, scale = TRUE)

  expect_equal(out_sparse$sdev, out_dense$sdev, tolerance = 1e-6)

  # dgRMatrix
  X_rmat <- as(X_sparse, "RsparseMatrix")
  out_rmat <- pca_pSVD(X_rmat, ncomp = 5, center = TRUE, scale = TRUE)
  expect_equal(out_rmat$sdev, out_dense$sdev, tolerance = 1e-6)

  # Sparse with center=FALSE, scale=FALSE
  out_raw <- pca_pSVD(X_sparse, ncomp = 5, center = FALSE, scale = FALSE)
  expect_s3_class(out_raw, "pca_pSVD")
})

# ============================================================================
# Zero-Variance Column Handling Tests
# ============================================================================

test_that("zero-variance columns handled correctly", {
  set.seed(123)
  X <- matrix(rnorm(1000), 100, 10)
  X[, 3] <- 5 # Zero variance
  X[, 7] <- -2 # Zero variance

  # scale=FALSE should remove with warning
  expect_warning(
    out <- pca_pSVD(X, ncomp = 5, center = TRUE, scale = FALSE),
    "Removed 2 zero-variance columns"
  )
  expect_equal(ncol(out$rotation), 5)

  # scale=TRUE should protect near-zero variance columns (sparse matrices only)
  skip_if_not_installed("Matrix")
  set.seed(456)
  X2_dense <- matrix(rnorm(1000), 100, 10)
  # Add columns with extremely small but non-zero variance
  # The threshold is sqrt(.Machine$double.eps) ≈ 1.49e-08
  threshold <- sqrt(.Machine$double.eps)
  X2_dense[, 3] <- 100 + seq(0, threshold * 0.5, length.out = 100)
  X2_dense[, 7] <- -50 + seq(0, threshold * 0.5, length.out = 100)
  X2_sparse <- Matrix::Matrix(X2_dense, sparse = TRUE)

  # Warning only appears for sparse matrices (dense uses base::scale)
  expect_warning(
    out2 <- pca_pSVD(X2_sparse, ncomp = 5, center = TRUE, scale = TRUE),
    "near-zero-variance columns"
  )
})

# ============================================================================
# Center/Scale Combinations Tests
# ============================================================================

test_that("all center/scale combinations work", {
  set.seed(999)
  X <- matrix(rnorm(2000), 100, 20)

  # center=FALSE, scale=TRUE
  out1 <- pca_pSVD(X, ncomp = 5, center = FALSE, scale = TRUE)
  expect_s3_class(out1, "pca_pSVD")
  expect_equal(out1$center, FALSE)
  expect_equal(out1$scale, TRUE)

  # center=FALSE, scale=FALSE
  out2 <- pca_pSVD(X, ncomp = 5, center = FALSE, scale = FALSE)
  expect_s3_class(out2, "pca_pSVD")
  expect_equal(out2$center, FALSE)
  expect_equal(out2$scale, FALSE)

  # Numeric center vector
  cvec <- colMeans(X)
  out3 <- pca_pSVD(X, ncomp = 5, center = cvec, scale = FALSE)
  expect_equal(out3$center, cvec)

  # Numeric scale vector
  svec <- apply(X, 2, sd)
  out4 <- pca_pSVD(X, ncomp = 5, center = TRUE, scale = svec)
  expect_equal(out4$scale, svec)
})

# ============================================================================
# Multilevel Support Tests
# ============================================================================

test_that("multilevel parameter adds mlpca class and processes data", {
  skip_if_not_installed("mixOmics")

  set.seed(42)
  # Create data with repeated measures: 10 subjects, 10 measurements each
  n_subjects <- 10
  n_reps <- 10
  n_vars <- 20

  # Simulate data with subject effects
  subject_effects <- matrix(rnorm(n_subjects * n_vars, sd = 5), n_subjects, n_vars)
  X <- do.call(rbind, lapply(1:n_subjects, function(i) {
    matrix(rep(subject_effects[i, ], n_reps), n_reps, n_vars, byrow = TRUE) +
      matrix(rnorm(n_reps * n_vars, sd = 1), n_reps, n_vars)
  }))

  multilevel_design <- data.frame(subject = rep(1:n_subjects, each = n_reps))

  # Without multilevel
  out1 <- pca_pSVD(X, ncomp = 3, multilevel = NULL)
  expect_equal(class(out1), c("pca_pSVD", "pca"))
  expect_null(out1$Xw)
  expect_null(out1$design)

  # With multilevel
  out2 <- pca_pSVD(X, ncomp = 3, multilevel = multilevel_design)
  expect_equal(class(out2), c("mlpca", "pca_pSVD", "pca"))

  # Check multilevel-specific outputs
  expect_false(is.null(out2$Xw))
  expect_false(is.null(out2$design))
  expect_equal(dim(out2$Xw), dim(X))
  expect_equal(nrow(out2$design), nrow(X))
  expect_equal(ncol(out2$design), 1)
})

test_that("multilevel matches mixOmics::pca results", {
  skip_if_not_installed("mixOmics")

  set.seed(123)
  # Create repeated measures data
  n_subjects <- 8
  n_reps <- 5
  n_vars <- 15

  subject_effects <- matrix(rnorm(n_subjects * n_vars, sd = 3), n_subjects, n_vars)
  X <- do.call(rbind, lapply(1:n_subjects, function(i) {
    matrix(rep(subject_effects[i, ], n_reps), n_reps, n_vars, byrow = TRUE) +
      matrix(rnorm(n_reps * n_vars, sd = 1), n_reps, n_vars)
  }))

  multilevel_design <- data.frame(subject = rep(1:n_subjects, each = n_reps))
  k <- 3

  # Run both implementations
  out_mx <- mixOmics::pca(X, ncomp = k, center = TRUE, scale = FALSE, multilevel = multilevel_design)
  out_psv <- pca_pSVD(X, ncomp = k, center = TRUE, scale = FALSE, multilevel = multilevel_design)

  # Compare variance explained (sign-invariant, ignore names)
  expect_equal(unname(out_psv$prop_expl_var$X), unname(out_mx$prop_expl_var$X[1:k]), tolerance = 1e-4)

  # Compare standard deviations (ignore names)
  expect_equal(unname(out_psv$sdev), unname(out_mx$sdev[1:k]), tolerance = 1e-5)

  # Verify Xw matrices are identical
  expect_equal(out_psv$Xw, out_mx$Xw, tolerance = 1e-10)

  # Verify design is stored correctly
  expect_equal(out_psv$design, out_mx$design)
})

test_that("multilevel validation catches errors", {
  skip_if_not_installed("mixOmics")

  set.seed(42)
  X <- matrix(rnorm(1000), 100, 10)

  # Wrong number of rows
  multilevel_wrong_rows <- data.frame(subject = 1:50)
  expect_error(
    pca_pSVD(X, ncomp = 3, multilevel = multilevel_wrong_rows),
    "unequal number of rows"
  )

  # Multiple columns
  multilevel_multi_col <- data.frame(
    subject = rep(1:10, each = 10),
    condition = rep(1:2, 50)
  )
  expect_error(
    pca_pSVD(X, ncomp = 3, multilevel = multilevel_multi_col),
    "single column"
  )
})

test_that("multilevel without mixOmics throws error", {
  # Mock requireNamespace to simulate mixOmics unavailable
  env <- asNamespace("base")
  old_requireNamespace <- base::requireNamespace

  unlockBinding("requireNamespace", env)
  on.exit(
    {
      assign("requireNamespace", old_requireNamespace, envir = env)
      lockBinding("requireNamespace", env)
    },
    add = FALSE
  )

  assign("requireNamespace", function(pkg, quietly = FALSE) {
    if (pkg == "mixOmics") {
      return(FALSE)
    }
    old_requireNamespace(pkg, quietly = quietly)
  }, envir = env)

  X <- matrix(rnorm(1000), 100, 10)
  multilevel_design <- data.frame(subject = rep(1:10, each = 10))

  expect_error(
    pca_pSVD(X, ncomp = 3, multilevel = multilevel_design),
    "mixOmics package"
  )
})

# ============================================================================
# mixOmics Error Handling Tests
# ============================================================================

test_that("mixOmics absence causes appropriate errors", {
  # Mock requireNamespace to simulate mixOmics unavailable
  env <- asNamespace("base")
  old_requireNamespace <- base::requireNamespace

  unlockBinding("requireNamespace", env)
  on.exit(
    {
      assign("requireNamespace", old_requireNamespace, envir = env)
      lockBinding("requireNamespace", env)
    },
    add = FALSE
  )

  assign("requireNamespace", function(pkg, quietly = FALSE) {
    if (pkg == "mixOmics") {
      return(FALSE)
    }
    old_requireNamespace(pkg, quietly = quietly)
  }, envir = env)

  X <- matrix(rnorm(1000), 100, 10)
  X[1, 1] <- NA

  expect_error(
    pca_pSVD(X, ncomp = 2),
    "mixOmics not installed"
  )

  X_complete <- matrix(runif(1000, 0.1, 10), 100, 10)
  expect_error(
    pca_pSVD(X_complete, ncomp = 2, logratio = "ILR"),
    "mixOmics not installed"
  )
})

# ============================================================================
# SVD Fallback Chain Tests
# ============================================================================

test_that("rARPACK fallback on sparse warns about densification", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("rARPACK")

  set.seed(42)
  X_dense <- matrix(rnorm(1000), 100, 10)
  X_dense[abs(X_dense) < 1] <- 0
  X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)

  # Mock to disable RSpectra but keep rARPACK available
  env <- asNamespace("base")
  old_requireNamespace <- base::requireNamespace

  unlockBinding("requireNamespace", env)
  on.exit(
    {
      assign("requireNamespace", old_requireNamespace, envir = env)
      lockBinding("requireNamespace", env)
    },
    add = FALSE
  )

  assign("requireNamespace", function(pkg, quietly = FALSE) {
    if (pkg == "RSpectra") {
      return(FALSE)
    }
    old_requireNamespace(pkg, quietly = quietly)
  }, envir = env)

  expect_warning(
    pca_pSVD(X_sparse, ncomp = 3),
    "rARPACK fallback.*will densify"
  )
})
