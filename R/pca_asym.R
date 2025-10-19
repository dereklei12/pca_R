#' PCA with asymmetric-dimension optimization (full decompositions only)
#' - If p >> n: eigen on K = Xs %*% t(Xs) (n×n), then recover V
#' - If n >> p: eigen on C = t(Xs) %*% Xs (p×p), then recover U
#' Otherwise: fall back to base::svd(Xs) full SVD
#'
#' Fields align with mixOmics::pca: sdev, loadings$X, variates$X, rotation, x, etc.
pca_asym <- function(
    X,
    ncomp = 2,
    center = TRUE,
    scale = FALSE,
    logratio = c("none", "CLR", "ILR"),
    ilr.offset = 0.001, V = NULL,
    multilevel = NULL,
    tol = 1e-12,
    verbose.call = FALSE,
    asym_threshold = 2 # e.g., treat as asymmetric if larger/smaller dim >= 2×
    ) {
  logratio <- match.arg(logratio)

  # ---- 0) Input/basic checks ----
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (anyNA(X)) {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("Missing values detected: please install mixOmics or impute first.")
    }
    return(mixOmics::pca(X,
      ncomp = ncomp, center = center, scale = scale,
      logratio = logratio, ilr.offset = ilr.offset, V = V,
      multilevel = multilevel, verbose.call = verbose.call
    ))
  }
  if (identical(logratio, "ILR")) {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("ILR requested: please install mixOmics or use logratio='none'/'CLR'.")
    }
    return(mixOmics::pca(X,
      ncomp = ncomp, center = center, scale = scale,
      logratio = logratio, ilr.offset = ilr.offset, V = V,
      multilevel = multilevel, verbose.call = verbose.call
    ))
  }

  nr <- nrow(X)
  nc <- ncol(X)
  if (is.null(nr) || is.null(nc) || min(nr, nc) < 2) stop("X must have >=2 rows and cols.")
  ncomp <- as.integer(max(1, min(ncomp, min(nr, nc) - 1L)))

  # Names
  sample.names <- rownames(X)
  if (is.null(sample.names)) sample.names <- paste0("sample_", seq_len(nr))
  var.names <- colnames(X)
  if (is.null(var.names)) var.names <- paste0("var_", seq_len(nc))

  # ---- 1) Standardization (consistent with prcomp/mixOmics semantics) ----
  Xs <- base::scale(X, center = center, scale = scale)
  scaled.center <- attr(Xs, "scaled:center")
  scaled.scale <- attr(Xs, "scaled:scale")

  # Remove zero-variance columns (may occur only when scale=FALSE)
  if (identical(scale, FALSE)) {
    v <- apply(Xs, 2, var)
    keep <- which(v > 0)
    if (length(keep) < ncol(Xs)) {
      warning(sprintf("Removed %d zero-variance columns.", ncol(Xs) - length(keep)))
      Xs <- Xs[, keep, drop = FALSE]
      var.names <- var.names[keep]
      nc <- ncol(Xs)
      ncomp <- min(ncomp, min(nr, nc) - 1L)
    }
  }

  # ---- 2) Choose decomposition route (full decomposition, but on smaller dimension) ----
  # Convention: treat as "asymmetric" when dimension ratio >= asym_threshold
  route <- "svd_full_X" # Default: full svd on Xs directly
  if (nr * asym_threshold <= nc) route <- "eigen_K" # p >> n
  if (nc * asym_threshold <= nr) route <- "eigen_C" # n >> p

  # ---- 3) Full decomposition (on smaller dimension), then recover singular vectors on the other side ----
  U <- NULL
  V <- NULL
  d <- NULL

  if (route == "eigen_K") {
    # p >> n: decompose K = Xs Xs^T (n×n)
    K <- Xs %*% t(Xs) # n x n
    eg <- eigen(K, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    d_all <- sqrt(vals)
    U_all <- eg$vectors
    # Filter very small singular values (avoid division by zero)
    keep <- which(d_all > tol * max(d_all, 0))
    if (length(keep) == 0) stop("Numerical rank = 0.")
    # Right singular vectors: V = Xs^T U D^{-1}
    V_all <- t(Xs) %*% sweep(U_all[, keep, drop = FALSE], 2, d_all[keep], "/")
    # Numerical orthogonalization correction (optional): V_all <- qr.Q(qr(V_all))
    kk <- min(ncomp, length(keep))
    d <- d_all[keep][seq_len(kk)]
    U <- U_all[, keep, drop = FALSE][, seq_len(kk), drop = FALSE]
    V <- V_all[, seq_len(kk), drop = FALSE]
  } else if (route == "eigen_C") {
    # n >> p: decompose C = Xs^T Xs (p×p)
    C <- crossprod(Xs) # p x p
    eg <- eigen(C, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    d_all <- sqrt(vals)
    V_all <- eg$vectors
    keep <- which(d_all > tol * max(d_all, 0))
    if (length(keep) == 0) stop("Numerical rank = 0.")
    # Left singular vectors: U = Xs V D^{-1}
    U_all <- Xs %*% sweep(V_all[, keep, drop = FALSE], 2, d_all[keep], "/")
    kk <- min(ncomp, length(keep))
    d <- d_all[keep][seq_len(kk)]
    V <- V_all[, keep, drop = FALSE][, seq_len(kk), drop = FALSE]
    U <- U_all[, seq_len(kk), drop = FALSE]
  } else {
    # Dimensions not extreme: direct full svd(Xs)
    sv <- base::svd(Xs, nu = ncomp, nv = ncomp)
    U <- sv$u[, seq_len(ncomp), drop = FALSE]
    V <- sv$v[, seq_len(ncomp), drop = FALSE]
    d <- sv$d[seq_len(ncomp)]
  }

  # ---- 4) Assemble and align with mixOmics ----
  sdev <- as.numeric(d / sqrt(nrow(Xs) - 1))
  rotation <- matrix(V,
    nrow = ncol(Xs), ncol = length(d),
    dimnames = list(var.names, paste0("PC", seq_len(length(d))))
  )
  x <- Xs %*% rotation
  colnames(x) <- paste0("PC", seq_len(ncol(x)))
  rownames(x) <- sample.names

  loadings <- list(X = rotation)
  variates <- list(X = x)

  var.tot <- sum(apply(Xs, 2, var))
  prop <- (sdev^2) / var.tot
  cumvar <- cumsum(prop)

  res <- list(
    call = if (isTRUE(verbose.call)) match.call() else NULL,
    X = structure(Xs, "scaled:center" = scaled.center, "scaled:scale" = scaled.scale),
    ncomp = length(d),
    center = center,
    scale = scale,
    names = list(sample = sample.names, col = var.names),
    sdev = sdev,
    loadings = loadings,
    variates = variates,
    prop_expl_var = list(X = prop),
    var.tot = var.tot,
    cum.var = list(X = cumvar),
    rotation = rotation,
    x = x
  )
  wanted <- c(
    "call", "X", "ncomp", "center", "scale", "names", "sdev",
    "loadings", "variates", "prop_expl_var", "var.tot", "cum.var",
    "rotation", "x"
  )
  res <- res[wanted]
  class(res) <- "pca"
  if (!is.null(multilevel)) class(res) <- c("mlpca", class(res))
  return(res)
}
