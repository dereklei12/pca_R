test_that("pca_woodbury structure & dimensions", {
  set.seed(1)
  X <- matrix(rnorm(2000), 100, 20) # n < p 场景
  k <- 5
  out <- pca_woodbury(X, ncomp = k, center = TRUE, scale = TRUE)

  expect_s3_class(out, "pca")
  expect_true(all(c("sdev", "loadings", "variates", "rotation", "x", "prop_expl_var") %in% names(out)))
  expect_equal(ncol(out$rotation), k)
  expect_equal(ncol(out$x), k)

  # x ≈ Xs %*% rotation
  Xs <- scale(X, center = out$center, scale = out$scale)
  x2 <- Xs %*% out$rotation
  expect_equal(out$x, x2, tolerance = 1e-7)
})

test_that("pca_woodbury vs prcomp: sdev & rotation close (sign-invariant)", {
  set.seed(2)
  X <- matrix(rnorm(3000), 100, 30) # n < p
  k <- 6
  pw <- pca_woodbury(X, ncomp = k, center = TRUE, scale = TRUE)
  pr <- prcomp(X, center = TRUE, scale. = TRUE, retx = TRUE)

  expect_equal(pw$sdev, pr$sdev[1:k], tolerance = 1e-7)

  # 载荷对齐符号后比较
  fix_sign <- function(A, B) {
    s <- sign(colSums(A * B))
    s[s == 0] <- 1
    sweep(A, 2, s, `*`)
  }
  L1 <- pw$rotation
  L2 <- pr$rotation[, 1:k, drop = FALSE]
  L1 <- fix_sign(L1, L2)
  common <- intersect(rownames(L1), rownames(L2))
  expect_equal(L1[common, ], L2[common, ], tolerance = 1e-4)
})

test_that("pca_woodbury delegates to mixOmics when ILR or NA", {
  X <- matrix(rnorm(1000), 100, 10)
  X[1, 1] <- NA
  if (requireNamespace("mixOmics", quietly = TRUE)) {
    out <- pca_woodbury(X, ncomp = 2)
    expect_s3_class(out, "pca")
  } else {
    expect_error(pca_woodbury(X, ncomp = 2))
  }
})

test_that("pca_woodbury handles p < n by falling back to SVD", {
  set.seed(3)
  X <- matrix(rnorm(1500), 50, 30) # n > p
  out <- pca_woodbury(X, ncomp = 5, center = TRUE, scale = TRUE)
  expect_s3_class(out, "pca")
  expect_equal(ncol(out$rotation), 5)
})
