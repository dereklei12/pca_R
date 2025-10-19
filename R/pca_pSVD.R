## ---- helpers: sparse-aware standardization & stats ----

.is_sparse <- function(X) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    return(FALSE)
  }
  inherits(X, "dgCMatrix") || inherits(X, "dgRMatrix") || inherits(X, "dgcMatrix")
}

.col_means <- function(X, sparse) {
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    as.numeric(Matrix::colMeans(X))
  } else {
    as.numeric(colMeans(X))
  }
}

.col_sumsq <- function(X, sparse) {
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    as.numeric(Matrix::colSums(X^2))
  } else {
    as.numeric(colSums(X^2))
  }
}

# Compute column centering/scaling (standard deviation, n-1) and return data for dense/sparse paths
.prepare_center_scale <- function(X, center, scale) {
  n <- nrow(X)
  p <- ncol(X)
  sparse <- .is_sparse(X)

  if (!sparse) {
    # dense: directly use scale, simple and robust
    Xs <- base::scale(X,
      center = if (isTRUE(center)) TRUE else if (is.numeric(center)) center else FALSE,
      scale  = if (isTRUE(scale)) TRUE else if (is.numeric(scale)) scale else FALSE
    )
    list(
      sparse = FALSE,
      Xs = Xs,
      cvec = attr(Xs, "scaled:center"),
      svec = attr(Xs, "scaled:scale")
    )
  } else {
    # sparse: avoid densification, explicitly compute cvec/svec
    cvec <- if (isTRUE(center)) {
      .col_means(X, TRUE)
    } else if (is.numeric(center)) {
      if (length(center) != p) stop("length(center) must equal ncol(X).")
      as.numeric(center)
    } else {
      NULL
    }

    if (isTRUE(scale)) {
      if (is.null(cvec)) cvec <- .col_means(X, TRUE)
      css <- .col_sumsq(X, TRUE) # sum of squares by col
      ssd <- css - n * (cvec^2) # ∑(x-c)^2
      ssd[ssd < 0] <- 0
      svec <- sqrt(pmax(ssd / max(1, n - 1), .Machine$double.eps))
      # Protect near-zero variance columns
      near0 <- which(svec <= sqrt(.Machine$double.eps))
      if (length(near0)) {
        svec[near0] <- sqrt(.Machine$double.eps)
        warning(sprintf("Detected %d near-zero-variance columns; protected by epsilon scaling.", length(near0)))
      }
    } else if (is.numeric(scale)) {
      if (length(scale) != p) stop("length(scale) must equal ncol(X).")
      svec <- as.numeric(scale)
    } else {
      svec <- NULL
    }

    list(sparse = TRUE, Xs = NULL, cvec = cvec, svec = svec)
  }
}

# Compute scores for sparse path without materializing Xs, based on (X - 1 c^T) diag(1/s) %*% V
.scores_from_center_scale <- function(X, rotation, cvec, svec) {
  W <- if (is.null(svec)) rotation else sweep(rotation, 1L, svec, "/")
  x <- as.matrix(X %*% W) # N x k
  if (!is.null(cvec)) {
    cw <- as.numeric(cvec %*% W) # 1 x k
    x <- sweep(x, 2L, cw, "-")
  }
  x
}

# Compute total variance ||Xs||_F^2 / (n-1); if Xs not materialized, compute directly from column statistics
.var_tot_from_stats <- function(X, Xs, n, cvec, svec, sparse) {
  if (!is.null(Xs)) {
    return(sum(Xs * Xs) / max(1, n - 1))
  }
  css <- .col_sumsq(X, sparse)
  if (is.null(cvec)) cvec <- rep(0, ncol(X))
  ssd <- css - n * (cvec^2)
  ssd[ssd < 0] <- 0
  if (is.null(svec)) {
    sum(ssd) / max(1, n - 1)
  } else {
    sum(ssd / (svec^2)) / max(1, n - 1)
  }
}

#' PCA via partial SVD (RSpectra/rARPACK)
#'
#' @param X numeric matrix or data.frame
#' @param ncomp number of principal components to keep
#' @param center,scale same as base::scale
#' @param max.iter,tol iterative svds/svd convergence parameters
#' @param logratio one of c("none","CLR","ILR"); allows ignoring CLR/ILR
#' @param ilr.offset offset for ILR transformation (used when delegating to mixOmics)
#' @param V matrix for ILR transformation (used when delegating to mixOmics)
#' @param multilevel data.frame or matrix with a single column indicating repeated
#'   measurements on each individual (subject ID). When provided, performs
#'   within-subject variation decomposition to remove between-subject effects
#'   before PCA analysis. Requires mixOmics package.
#' @param verbose.call whether to save the call in the result
#' @return list with class "pca" (aligned with mixOmics::pca)
pca_pSVD <- function(X,
                     ncomp = 2,
                     center = TRUE,
                     scale = FALSE,
                     max.iter = 500, tol = 1e-9,
                     logratio = c("none", "CLR", "ILR"),
                     ilr.offset = 0.001, V = NULL,
                     multilevel = NULL,
                     verbose.call = FALSE) {
  ## -------- Input/checks --------
  logratio <- match.arg(logratio)

  if (is.data.frame(X)) X <- as.matrix(X)
  # Check if numeric or sparse matrix
  is_sparse <- .is_sparse(X)
  if (!is_sparse && !is.numeric(X)) stop("X must be a numeric matrix/data.frame.")

  n <- nrow(X)
  p <- ncol(X)
  if (is.null(n) || is.null(p) || min(n, p) < 2L) stop("X must have >= 2 rows and >= 2 columns.")
  if (!is.numeric(ncomp) || ncomp < 1) stop("ncomp must be >= 1")
  ncomp <- as.integer(ncomp)

  sample.names <- rownames(X)
  if (is.null(sample.names)) sample.names <- paste0("sample_", seq_len(n))
  var.names <- colnames(X)
  # Use "V" prefix to match mixOmics::pca naming convention
  if (is.null(var.names)) var.names <- paste0("V", seq_len(p))

  has_na <- anyNA(X)
  if (has_na || identical(logratio, "ILR")) {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("Missing values or ILR requested but mixOmics not installed. Install mixOmics or provide complete data / non-ILR logratio.")
    }
    return(
      mixOmics::pca(
        X,
        ncomp = ncomp, center = center, scale = scale,
        max.iter = max.iter, tol = tol,
        logratio = logratio, ilr.offset = ilr.offset, V = V,
        multilevel = multilevel, verbose.call = verbose.call
      )
    )
  }
  if (identical(logratio, "CLR")) logratio <- "none"

  ## -------- Multilevel handling --------
  Xw <- NULL # Store within-variation matrix for output
  if (!is.null(multilevel)) {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("multilevel analysis requires mixOmics package. Install mixOmics or set multilevel = NULL.")
    }

    # Convert to data.frame and validate
    multilevel <- data.frame(multilevel)

    if (nrow(X) != nrow(multilevel)) {
      stop("unequal number of rows in 'X' and 'multilevel'.")
    }

    if (ncol(multilevel) != 1) {
      stop("'multilevel' should have a single column for the repeated measurements.")
    }

    # Convert factor levels to numeric for subject IDs
    multilevel[, 1] <- as.numeric(factor(multilevel[, 1]))

    # Apply within-variation decomposition to remove between-subject effects
    Xw <- mixOmics::withinVariation(X, design = multilevel)
    X <- Xw

    # Update dimensions after within-variation processing
    n <- nrow(X)
    p <- ncol(X)
  }

  ## -------- Standardization --------
  prep <- .prepare_center_scale(X, center, scale)
  Xs <- prep$Xs
  cvec <- prep$cvec
  svec <- prep$svec
  sparse <- prep$sparse

  # For dense matrices, remove zero-variance columns when scale=FALSE
  if (!sparse && identical(scale, FALSE)) {
    col_var <- apply(Xs, 2L, var)
    keep <- which(col_var > 0)
    if (length(keep) < ncol(Xs)) {
      warning(sprintf("Removed %d zero-variance columns.", ncol(Xs) - length(keep)))
      Xs <- Xs[, keep, drop = FALSE]
      var.names <- var.names[keep]
    }
  }

  nn <- if (is.null(Xs)) n else nrow(Xs)
  pp <- if (is.null(Xs)) p else ncol(Xs)
  ncap <- min(nn, pp) - 1L
  if (ncomp > ncap) {
    warning(sprintf("ncomp truncated to %d (<= min(nrow, ncol) - 1).", ncap))
    ncomp <- ncap
  }
  if (ncomp < 1L) stop("ncomp became < 1 after dimension checks.")

  ## -------- Partial SVD --------
  svds_ok <- FALSE
  V <- NULL
  d <- NULL

  if (requireNamespace("RSpectra", quietly = TRUE)) {
    if (!sparse) {
      svd_res <- try(
        RSpectra::svds(Xs,
          k = ncomp, nu = 0, nv = ncomp,
          opts = list(maxitr = max.iter, tol = tol)
        ),
        silent = TRUE
      )
    } else {
      # Sparse: implicit center/scale
      svd_res <- try(
        RSpectra::svds(X,
          k = ncomp, nu = 0, nv = ncomp,
          opts = list(
            maxitr = max.iter, tol = tol,
            center = if (is.null(cvec)) FALSE else cvec,
            scale = if (is.null(svec)) FALSE else svec
          )
        ),
        silent = TRUE
      )
    }
    if (!inherits(svd_res, "try-error")) {
      V <- svd_res$v
      d <- svd_res$d
      svds_ok <- TRUE
    }
  }

  if (!svds_ok && requireNamespace("rARPACK", quietly = TRUE)) {
    # rARPACK lacks center/scale options; will densify sparse matrices
    if (sparse) {
      warning("rARPACK fallback on sparse input will densify; consider installing RSpectra.")
      Xs <- base::scale(as.matrix(X),
        center = if (is.null(cvec)) FALSE else cvec,
        scale  = if (is.null(svec)) FALSE else svec
      )
    }
    svd_res <- try(
      rARPACK::svds(Xs, k = ncomp, opts = list(maxitr = max.iter, tol = tol)),
      silent = TRUE
    )
    if (!inherits(svd_res, "try-error")) {
      V <- svd_res$v
      d <- svd_res$d
      svds_ok <- TRUE
    }
  }

  if (!svds_ok) {
    warning("svds unavailable or failed; falling back to base::svd (full decomposition).")
    if (is.null(Xs)) {
      Xs <- base::scale(as.matrix(X),
        center = if (is.null(cvec)) FALSE else cvec,
        scale  = if (is.null(svec)) FALSE else svec
      )
    }
    svd_full <- base::svd(Xs, nu = 0, nv = ncomp)
    V <- svd_full$v[, seq_len(ncomp), drop = FALSE]
    d <- svd_full$d[seq_len(ncomp)]
  }

  V <- matrix(V,
    nrow = if (is.null(Xs)) p else ncol(Xs), ncol = ncomp,
    dimnames = list(var.names, NULL)
  )
  colnames(V) <- paste0("PC", seq_len(ncomp))
  rotation <- V

  ## -------- Assemble result --------
  sdev <- as.numeric(d) / sqrt(max(1, nn - 1))

  if (!is.null(Xs)) {
    x <- Xs %*% rotation
  } else {
    x <- .scores_from_center_scale(X, rotation, cvec, svec)
  }
  rownames(x) <- sample.names

  var.tot <- .var_tot_from_stats(X, Xs, n, cvec, svec, sparse)
  prop <- (sdev^2) / var.tot
  names(prop) <- paste0("PC", seq_len(ncomp))
  cumvar <- cumsum(prop)

  # Create loadings without rownames to match mixOmics::pca behavior
  loadings_mat <- rotation
  rownames(loadings_mat) <- NULL
  loadings <- list(X = loadings_mat)
  variates <- list(X = x)

  # For mixOmics compatibility: dense path saves Xs; sparse path keeps original X with center/scale attributes
  X_field <- if (!is.null(Xs)) {
    structure(Xs,
      "scaled:center" = attr(Xs, "scaled:center"),
      "scaled:scale" = attr(Xs, "scaled:scale")
    )
  } else {
    structure(X,
      "scaled:center" = if (is.null(cvec)) FALSE else cvec,
      "scaled:scale" = if (is.null(svec)) FALSE else svec
    )
  }

  # Assemble result list
  res <- list(
    call = if (isTRUE(verbose.call)) match.call() else match.call(expand.dots = FALSE)[1],
    X = X_field,
    ncomp = ncomp,
    center = center,
    scale = scale,
    names = list(X = var.names, sample = sample.names),
    sdev = sdev,
    loadings = loadings,
    variates = variates,
    prop_expl_var = list(X = as.numeric(prop)),
    var.tot = var.tot,
    cum.var = as.numeric(cumvar),
    rotation = rotation,
    x = x
  )

  # Append multilevel-specific fields if applicable
  if (!is.null(multilevel)) {
    res <- c(res, list(Xw = Xw, design = multilevel))
  }

  # Set class hierarchy
  class(res) <- if (!is.null(multilevel)) {
    c("mlpca", "pca_pSVD", "pca")
  } else {
    c("pca_pSVD", "pca")
  }

  res
}
