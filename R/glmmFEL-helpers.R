##### glmmFEL-helpers.R ######################################################
# This file contains small internal helpers used across the matrix-only
# development branch of glmmFEL.
#
# Design goals for this branch:
#   * Keep the public API centered on (y, X, Z).
#   * Keep the random-effects covariance to a single variance component:
#       G = tau^2 I_q
#   * Retain the fully exponential (FE) mean and covariance correction machinery
#     (the point of the package), while removing the formula + structured-G
#     parsing layer to reduce complexity.
#
# Provenance notes (high level):
#   * FE Laplace EM machinery follows Karl, Yang, and Lohr (2014), and is
#     aligned with the fully exponential Laplace work used in joint models by
#     Rizopoulos, Verbeke, and Lesaffre (2009).
#   * The RSPL/MSPL pseudo-likelihood code paths are adapted from RealVAMS
#     (Broatch, Green, and Karl, 2018) which closely follows Wolfinger and
#     O'Connell (1993) working-response / working-weights linearization.
################################################################################

#' Resolve supported family specifications
#'
#' @description
#' Internal helper that normalizes user family input into one of the canonical
#' labels used by the internal FE/Laplace code:
#' \itemize{
#'   \item \code{"binomial_probit"}
#'   \item \code{"binomial_logit"}
#'   \item \code{"poisson_log"}
#' }
#'
#' Users may pass either the canonical character labels above, or common
#' `stats::family()` objects:
#' \itemize{
#'   \item \code{stats::binomial(link = "probit")}
#'   \item \code{stats::binomial(link = "logit")}
#'   \item \code{stats::poisson(link = "log")}
#' }
#'
#' @param family Character label or a `stats::family()` object.
#' @return A canonical family label.
#' @keywords internal
glmmfe_resolve_family <- function(family) {
  if (is.character(family)) {
    if (length(family) != 1L || is.na(family)) stop("family must be a single label.")
    fam <- tolower(family)
    if (fam %in% c("binomial_probit", "binomial-probit", "probit")) return("binomial_probit")
    if (fam %in% c("binomial_logit", "binomial-logit", "logit", "logistic")) return("binomial_logit")
    if (fam %in% c("poisson_log", "poisson-log", "poisson")) return("poisson_log")
    stop("Unsupported family label: ", family)
  }

  if (inherits(family, "family")) {
    if (identical(family$family, "binomial") && identical(family$link, "probit")) return("binomial_probit")
    if (identical(family$family, "binomial") && identical(family$link, "logit"))  return("binomial_logit")
    if (identical(family$family, "poisson")  && identical(family$link, "log"))    return("poisson_log")
    stop("Unsupported family object: family=", family$family, ", link=", family$link)
  }

  stop("family must be either a character label or a stats::family() object.")
}

#' Resolve approximation labels
#'
#' @description
#' Internal helper that normalizes approximation input to one of:
#' \itemize{
#'   \item \code{"Laplace"}
#'   \item \code{"FE_mean"}
#'   \item \code{"FE_full"} (alias \code{"FE"})
#'   \item \code{"RSPL"}
#'   \item \code{"MSPL"}
#' }
#'
#' @param approx Character label.
#' @return Canonical approximation label.
#' @keywords internal
glmmfe_resolve_approx <- function(approx) {
  if (length(approx) != 1L) approx <- approx[1L]
  if (!is.character(approx) || length(approx) != 1L || is.na(approx)) stop("approx must be a character string.")

  a <- tolower(approx)

  if (a %in% c("laplace", "l1")) return("Laplace")
  if (a %in% c("fe", "fe_full", "fefull")) return("FE_full")
  if (a %in% c("fe_mean", "femean", "fe-mean")) return("FE_mean")
  if (a %in% c("rspl")) return("RSPL")
  if (a %in% c("mspl")) return("MSPL")

  stop("Unsupported approx: ", approx,
       ". Expected one of Laplace, FE_mean, FE_full/FE, RSPL, MSPL.")
}

#' Coerce a fixed-effects design matrix to a numeric base matrix
#'
#' @param X Numeric matrix-like object.
#' @return A base numeric matrix.
#' @keywords internal
glmmfe_as_X <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  X
}

#' Coerce a random-effects design matrix to a sparse dgCMatrix
#'
#' @param Z Numeric matrix-like object (dense or sparse).
#' @param n Optional expected number of rows (used for defensive checks).
#' @return A sparse `dgCMatrix` with numeric storage.
#' @keywords internal
glmmfe_as_Z <- function(Z, n = NULL) {
  if (is.null(Z)) stop("Z must not be NULL.")

  Zs <- Matrix::Matrix(Z, sparse = TRUE)

  # Coerce to dgCMatrix for consistent sparse ops
  if (!inherits(Zs, "dgCMatrix")) {
    Zs <- methods::as(methods::as(methods::as(Zs, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }

  # Optional dimension check (used by PL engine)
  if (!is.null(n) && nrow(Zs) != n) {
    stop(sprintf("Z has %d rows but expected n = %d.", nrow(Zs), n))
  }

  # Ensure numeric storage
  if (!is.double(Zs@x)) storage.mode(Zs@x) <- "double"

  Zs
}

#' Construct a glmmFEL fitted-model object
#'
#' @description
#' Internal helper to standardize creation of the fitted-model object returned
#' by [glmmFEL()] and related engines.
#'
#' This branch stores only a single variance component `tau2` with
#' \eqn{G = \tau^2 I_q}. More structured covariance parameterizations and the
#' formula wrapper are intentionally removed to reduce complexity.
#'
#' @param y Numeric response vector of length `n`.
#' @param X Fixed-effects design matrix `n x p`.
#' @param Z Random-effects design matrix `n x q` (stored as sparse dgCMatrix).
#' @param beta Fixed-effect estimates (length `p`).
#' @param eta Random-effect predictions (length `q`).
#' @param tau2 Non-negative scalar variance component.
#' @param G Optional `q x q` covariance matrix (defaults to `tau2 * I_q`).
#' @param vcov_beta Optional `p x p` covariance matrix for `beta`.
#' @param vcov_eta Optional `q x q` covariance matrix for `eta`.
#' @param cov_beta_eta Optional `p x q` cross-covariance block.
#' @param var_eta Optional alias for prediction-error covariance of `eta`.
#' @param family Canonical family label.
#' @param approx Canonical approximation label.
#' @param control List of control settings used.
#' @param convergence List containing convergence information.
#' @param logLik Approximate log-likelihood/objective (may be NA).
#' @param call Captured match.call().
#' @param reml Logical flag (used by RSPL/MSPL only; may be NULL otherwise).
#'
#' @return A list of class `c("glmmFELMod", "glmmFEL")`.
#' @keywords internal
glmmfe_new_fit <- function(
    y,
    X,
    Z,
    beta,
    eta,
    tau2,
    G           = NULL,
    vcov_beta    = NULL,
    vcov_eta     = NULL,
    cov_beta_eta = NULL,
    var_eta      = NULL,
    family,
    approx,
    control,
    convergence,
    logLik   = NA_real_,
    call     = NULL,
    reml     = NULL
) {
  y <- as.numeric(y)
  X <- glmmfe_as_X(X)
  Z <- glmmfe_as_Z(Z)

  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)

  if (nrow(X) != n || nrow(Z) != n) {
    stop("Inconsistent dimensions in glmmfe_new_fit(): y, X, Z must have the same number of rows.")
  }

  beta <- as.numeric(beta)
  eta  <- as.numeric(eta)

  if (length(beta) != p) stop("Length of 'beta' (", length(beta), ") does not match ncol(X) = ", p, ".")
  if (length(eta)  != q) stop("Length of 'eta' (", length(eta),  ") does not match ncol(Z) = ", q, ".")

  tau2 <- as.numeric(tau2)
  if (length(tau2) != 1L || !is.finite(tau2) || tau2 < 0) {
    stop("'tau2' must be a single non-negative finite numeric scalar.")
  }

  if (is.null(G)) {
    G <- Matrix::Diagonal(q, x = rep.int(tau2, q))
  } else {
    if (!inherits(G, "Matrix") && !is.matrix(G)) stop("'G' must be a base matrix or a Matrix object if supplied.")
    if (!inherits(G, "Matrix")) G <- Matrix::Matrix(G, sparse = FALSE)
  }
  if (!all(dim(G) == c(q, q))) stop("Dimension of 'G' must be q x q, where q = ncol(Z).")

  if (!is.null(vcov_beta) && !all(dim(vcov_beta) == c(p, p))) stop("vcov_beta must be a p x p matrix, where p = ncol(X).")
  if (!is.null(vcov_eta)  && !all(dim(vcov_eta)  == c(q, q))) stop("vcov_eta must be a q x q matrix, where q = ncol(Z).")
  if (!is.null(cov_beta_eta) && !all(dim(cov_beta_eta) == c(p, q))) {
    stop("cov_beta_eta must be a p x q matrix (p = ncol(X), q = ncol(Z)).")
  }

  if (is.null(var_eta)) var_eta <- vcov_eta

  res <- list(
    y            = y,
    X            = X,
    Z            = Z,
    beta         = beta,
    eta          = eta,
    tau2         = tau2,
    G            = G,
    vcov_beta    = vcov_beta,
    vcov_eta     = vcov_eta,
    cov_beta_eta = cov_beta_eta,
    var_eta      = var_eta,
    family       = family,
    approx       = approx,
    control      = control,
    convergence  = convergence,
    logLik       = logLik,
    call         = call,
    reml         = if (!is.null(reml)) isTRUE(reml) else NULL
  )

  class(res) <- c("glmmFELMod", "glmmFEL")
  res
}

#' Internal infix helper for defaulting `NULL` values
#'
#' Returns `x` if it is not `NULL`, otherwise returns `y`.
#'
#' @param x,y Objects.
#' @return `x` if not `NULL`, else `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Fast trace of a matrix product
#'
#' @description
#' Internal helper to compute \eqn{\mathrm{tr}(A B)} without forming `A %*% B`:
#' \deqn{\mathrm{tr}(A B) = \sum (A \circ B^\top).}
#'
#' This identity is used repeatedly in the fully exponential (FE) trace
#' corrections (Karl, Yang, and Lohr, 2014, Appendix B), where a naive
#' `diag(A %*% B)` would allocate the full product.
#'
#' @param A,B Numeric matrices with identical dimensions.
#' @return A single numeric scalar equal to `tr(A %*% B)`.
#' @keywords internal
glmmfe_trAB <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)

  if (!all(dim(A) == dim(B))) {
    stop("glmmfe_trAB(): A and B must have the same dimensions.")
  }

  # tr(A %*% B) = sum(A * t(B))
  sum(A * t(B))
}


#' @noRd
glmmfe_frob <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  if (!all(dim(A) == dim(B))) stop("glmmfe_frob(): dims must match.")
  sum(A * B)
}


##### glmmFEL-helpers.R #######################################################
