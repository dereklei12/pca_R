#' Principal Components Analysis via EM-PPCA (Woodbury-optimized)
#'
#' A drop-in PCA function that mirrors the structure and outputs of the
#' existing `pca()` (mixOmics-style) but computes components using an
#' EM algorithm for Probabilistic PCA (PPCA). It leverages the
#' Sherman–Morrison–Woodbury identity implicitly by restricting all linear
#' solves to \code{k x k} systems (\code{k = ncomp}), using Cholesky-based
#' solves for speed and numerical stability. Missing values (\code{NA})
#' are handled in-iteration via model-based reconstruction, similar to
#' \code{pcaMethods::ppca}.
#'
#' The function keeps the same public signature and return structure as the
#' given `pca()` reference (including logratio transforms and multilevel
#' preprocessing). It returns the same components, so it should be usable as a
#' replacement in downstream code.
#'
#' @param X see mixOmics::pca
#' @param ncomp see mixOmics::pca
#' @param center see mixOmics::pca
#' @param scale see mixOmics::pca
#' @param max.iter Integer, maximum EM iterations (default 500)
#' @param tol Numeric, relative objective change tolerance (default 1e-9)
#' @param logratio see mixOmics::pca
#' @param ilr.offset see mixOmics::pca
#' @param V see mixOmics::pca (ILR back-transformation helper)
#' @param multilevel see mixOmics::pca
#' @param verbose.call see mixOmics::pca
#'
#' @return A list with class c("pca", "prcomp") like mixOmics::pca.
#' Names/structure are kept identical: call, X, ncomp, center, scale, names,
#' sdev, loadings (list X=...), variates (list X=...), prop_expl_var, var.tot,
#' cum.var, rotation, x; plus Xw / design if multilevel.
#'
#' @export
pca_em_woodbury <- function(X,
                            ncomp = 2,
                            center = TRUE,
                            scale = FALSE,
                            max.iter = 500,
                            tol = 1e-09,
                            logratio = c('none','CLR','ILR'),
                            ilr.offset = 0.001,
                            V = NULL,
                            multilevel = NULL,
                            verbose.call = FALSE) {
  
  # ----------- input checks (mirror mixOmics::pca style) -----------------
  arg.call <- match.call()
  
  X <- .check_numeric_matrix(X)
  
  X.names <- colnames(X); if (is.null(X.names)) X.names <- paste0("V", seq_len(ncol(X)))
  ind.names <- rownames(X); if (is.null(ind.names)) ind.names <- seq_len(nrow(X))
  
  if (is.null(ncomp)) ncomp <- min(nrow(X), ncol(X))
  if (!is.numeric(ncomp) || !is.finite(ncomp) || ncomp < 1) stop("invalid 'ncomp'")
  ncomp <- as.integer(round(ncomp))
  if (ncomp > min(nrow(X), ncol(X))) stop("use smaller 'ncomp'")
  
  logratio_choice <- c('CLR', 'ILR', 'none')
  logratio <- logratio_choice[pmatch(logratio, logratio_choice)]
  logratio <- match.arg(logratio)
  if (logratio != 'none' && any(X < 0)) stop("'X' contains negative values; cannot log-transform")
  
  if (!is.logical(center)) {
    if (!is.numeric(center) || length(center) != ncol(X))
      stop("'center' must be logical or a numeric vector of length ncol(X)")
  }
  if (!is.logical(scale)) {
    if (!is.numeric(scale) || length(scale) != ncol(X))
      stop("'scale' must be logical or a numeric vector of length ncol(X)")
  }
  
  if (is.null(max.iter) || !is.numeric(max.iter) || !is.finite(max.iter) || max.iter < 1)
    stop("invalid 'max.iter'")
  max.iter <- as.integer(round(max.iter))
  
  if (is.null(tol) || !is.numeric(tol) || !is.finite(tol) || tol < 0)
    stop("invalid 'tol'")
  
  # ---------------- logratio transform (external helpers expected) --------
  if (logratio != 'none') {
    if (!requireNamespace('mixOmics', quietly = TRUE)) {
      stop("logratio != 'none' but package 'mixOmics' is not available. Install it or set logratio='none'.")
    }
    if (is.null(V) && logratio == 'ILR') {
      V <- mixOmics::clr.backtransfo(X)
    }
    X <- mixOmics::logratio.transfo(X = X, logratio = logratio,
                                    offset = if (logratio == 'ILR') ilr.offset else 0)
  }
  if (ncomp > min(nrow(X), ncol(X))) stop("use smaller 'ncomp'")
  
  # ---------------- multilevel preprocessing ------------------------------
  if (!is.null(multilevel)) {
    multilevel <- data.frame(multilevel)
    if (nrow(X) != nrow(multilevel)) stop("unequal rows in 'X' and 'multilevel'")
    if (ncol(multilevel) != 1) stop("'multilevel' must have a single column")
    multilevel[,1] <- as.numeric(factor(multilevel[,1]))
    Xw <- withinVariation(X, design = multilevel)
    X <- Xw
  }
  
  .check_zero_var_columns(X, scale = scale)
  X <- scale(X, center = center, scale = scale)
  sc <- attr(X, 'scaled:scale'); cen <- attr(X, 'scaled:center')
  
  result <- list(call = match.call(), X = X, ncomp = ncomp,
                 scale = if (is.null(sc)) FALSE else sc,
                 center = if (is.null(cen)) FALSE else cen,
                 names = list(X = X.names, sample = ind.names))
  
  # ---------------- EM-PPCA core -----------------------------------------
  # Work on a copy M so that result$X preserves NA pattern
  N <- nrow(X); D <- ncol(X); k <- ncomp
  Obs <- !is.na(X)
  hidden_idx <- which(!Obs)
  missing <- length(hidden_idx)
  
  M <- X
  if (missing) M[hidden_idx] <- 0
  
  # initialization of loadings C (D x k) and noise variance ss
  # Start with random C and a least-squares X, like pcaMethods::ppca
  C <- matrix(rnorm(D * k), D, k)
  CtC <- .crossprod(C)                        # k x k
  # Initial scores via LS: X0 = M C (C^T C)^{-1}
  Xscores <- M %*% C %*% .solve_sympd(CtC)
  # Initial reconstruction and noise variance
  recon <- Xscores %*% t(C)
  if (missing) recon[hidden_idx] <- 0
  ss <- sum((recon - M)^2) / (N * D - missing + 1e-16)
  if (!is.finite(ss) || ss <= 0) ss <- 1.0
  
  # EM loop
  old_obj <- Inf
  iter <- 0L
  repeat {
    iter <- iter + 1L
    
    # E-step covariance in latent space: Sx = (I + C^T C / ss)^{-1}
    CtC <- .crossprod(C)
    Sx <- .solve_sympd(diag(k) + CtC / ss)   # k x k (Woodbury-equivalent)
    
    ss_old <- ss
    
    # If missing, impute with current projection (model-based)
    if (missing) {
      proj <- Xscores %*% t(C)               # N x D
      M[hidden_idx] <- proj[hidden_idx]
    }
    
    # E-step expected scores: X = M C Sx / ss
    Xscores <- (M %*% C) %*% (Sx / ss)       # N x k
    
    # M-step for C: C = (M^T X) (X^T X + N Sx)^{-1}
    SumXtX <- .crossprod(Xscores)            # k x k
    RHS <- crossprod(M, Xscores)             # D x k
    C <- RHS %*% .solve_sympd(SumXtX + N * Sx)
    
    # Update ss (noise variance)
    CtC <- .crossprod(C)
    # Compute reconstruction error term ||C X^T - M^T||_F^2 efficiently
    # (D x N) - we keep for clarity & robustness on moderate sizes
    E <- C %*% t(Xscores) - t(M)
    err_term <- sum(E * E)
    ss <- ( err_term + N * sum(diag(CtC %*% Sx)) + missing * ss_old ) / (N * D)
    
    # Objective monitor (relative change); use log|Sx| via chol for stability
    logdet_Sx <- tryCatch({ R <- chol(Sx); 2 * sum(log(diag(R))) }, error = function(e) log(det(Sx)))
    obj <- N * (D * log(ss) + sum(diag(Sx)) - logdet_Sx) + sum(diag(SumXtX)) - missing * log(ss_old)
    rel_ch <- abs(1 - obj / old_obj)
    old_obj <- obj
    
    if ((rel_ch < tol && iter > 5L) || iter >= max.iter) break
  }
  
  # Orthonormalize columns of C to get principal directions; rotate scores
  C <- .orth(C)                               # D x k, columns orthonormal
  tmp_scores <- M %*% C                       # N x k (in working space)
  # Rotate within the k-dim subspace to order by variance
  ev <- eigen(cov(tmp_scores), symmetric = TRUE)
  vals <- ev$values
  vecs <- ev$vectors
  C <- C %*% vecs                              # rotate loadings
  Xscores <- tmp_scores %*% vecs               # rotate scores
  
  # Prepare outputs
  # sdev should be singular values so that .add_... divides by sqrt(n-1)
  sdev <- sqrt(pmax(vals, 0)) * sqrt(max(1, N - 1))
  
  # Back-transform loadings for ILR case to clr-space (match mixOmics behavior)
  if (logratio == 'ILR') {
    if (is.null(V)) stop("ILR back-transformation matrix 'V' is required but missing")
    loadings <- V %*% C[, seq_len(ncomp), drop = FALSE]
    variates <- Xscores[, seq_len(ncomp), drop = FALSE]
  } else {
    loadings <- C[, seq_len(ncomp), drop = FALSE]
    variates <- Xscores[, seq_len(ncomp), drop = FALSE]
  }
  
  result <- .add_sdev_loadings_and_variates(result, sdev = sdev,
                                            loadings = loadings,
                                            variates = variates)
  result <- c(result, .get_var_stats(X = result$X, sdev = result$sdev))
  
  expected_output_names <- c("call", "X", "ncomp", "center", "scale", "names",
                             "sdev", "loadings", "variates", "prop_expl_var",
                             "var.tot", "cum.var", "rotation", "x")
  result <- result[expected_output_names]
  
  if (!is.null(multilevel)) result <- c(result, list(Xw = Xw, design = multilevel))
  
  if (verbose.call) {
    c.simple <- result$call
    result$call <- mget(names(formals()))
    result$call <- append(c.simple, result$call)
    names(result$call)[1] <- "simple.call"
  }
  
  class(result) <- c("pca", "prcomp")
  if (!is.null(multilevel)) class(result) <- c("mlpca", class(result))
  
  return(invisible(result))
}

# --------------------- helpers (self-contained) ------------------------------

# Numeric matrix check (placeholder – expected to exist in mixOmics env)
.check_numeric_matrix <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) stop("'X' must be a numeric matrix")
  storage.mode(X) <- "double"
  X
}

# Zero-variance check (placeholder – expected to exist in mixOmics env)
.check_zero_var_columns <- function(X, scale = FALSE) {
  v <- apply(X, 2L, function(z) var(z, na.rm = TRUE))
  if (isTRUE(scale) && any(v == 0, na.rm = TRUE))
    stop("'scale = TRUE' but some columns have zero variance; filter them first")
  invisible(NULL)
}

# crossprod shortcut
.crossprod <- function(A) {
  # returns t(A) %*% A for dense matrices
  crossprod(A)
}

# Solve symmetric positive-definite system via Cholesky; fallback to solve
.solve_sympd <- function(A) {
  out <- tryCatch({ R <- chol(A); chol2inv(R) }, error = function(e) NULL)
  if (is.null(out)) solve(A) else out
}

# Orthonormalize columns (like MATLAB orth) using economy SVD
.orth <- function(A) {
  s <- svd(A, nu = 0, nv = ncol(A))
  # Right singular vectors of A (in column space) give orthonormal basis after mapping
  # For column-orthonormalization, use QR for efficiency when D >> k
  q <- tryCatch({ qr.Q(qr(A))[, seq_len(ncol(A)), drop = FALSE] }, error = function(e) NULL)
  if (!is.null(q)) return(q)
  # fallback: normalize via SVD in column space
  u <- svd(A, nu = ncol(A), nv = 0)$u
  u
}

# Add sdev, loadings and components to result (copied from user-supplied code)
.add_sdev_loadings_and_variates <- function(result, sdev, loadings, variates = NULL) {
  ncomp <- result$ncomp
  pc_names <- paste0("PC", seq_len(ncomp))
  
  Xint <- result$X
  Xint[is.na(Xint)] <- 0
  
  sdev <- sdev[seq_len(ncomp)] / sqrt(max(1, nrow(Xint) - 1))
  loadings <- loadings[, seq_len(ncomp), drop = FALSE]
  if (is.null(variates)) variates <- Xint %*% loadings
  
  names(sdev) <- pc_names
  dimnames(loadings) <- list(colnames(Xint), pc_names)
  dimnames(variates) <- list(rownames(Xint), pc_names)
  
  result[c('sdev', 'loadings', 'variates', 'x', 'rotation')] <- list(
    sdev,
    list(X = loadings),
    list(X = variates),
    variates,
    loadings
  )
  result
}

# Variance stats (copied from user-supplied code)
.get_var_stats <- function(X, sdev) {
  var.tot <- sum(X^2, na.rm = TRUE) / max(1, nrow(X) - 1)
  prop_expl_var <- sdev^2 / var.tot
  cum.var <- cumsum(prop_expl_var)
  prop_expl_var <- list(X = prop_expl_var)
  list(prop_expl_var = prop_expl_var, var.tot = var.tot, cum.var = cum.var)
}
