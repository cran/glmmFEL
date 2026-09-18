#' Compute FE trace-corrections for diagonal G = tau2 * I_q
#'
#' @description
#' Computes fully-exponential (FE) trace corrections for the random-effects
#' posterior mean and (optionally) covariance under a **single variance
#' component** random-effects covariance structure
#' \deqn{G = \tau^2 I_q.}
#'
#' This implements the same algebra as the mvglmmRank binary/Poisson engines,
#' and matches Appendix B of Karl, Yang, and Lohr (2014):
#'  - Define \eqn{\Sigma} as the *negative Hessian* at the tilted mode, and
#'    \eqn{\Sigma^{-1}} as the Laplace covariance \eqn{\mathrm{Var}(\eta\mid y)}.
#'  - Let \eqn{D_k = -\partial\Sigma/\partial c_k} and
#'    \eqn{D_{kd}^{(2)} = -\partial^2\Sigma/\partial c_k\partial c_d}.
#'  - Then the mean and covariance corrections use trace terms involving
#'    \eqn{\mathrm{tr}(\Sigma^{-1} D_k)} and
#'    \eqn{\mathrm{tr}(\Sigma^{-1} D_{kk}^{(2)}) + \mathrm{tr}(\Sigma^{-1} D_k \Sigma^{-1} D_k)}.
#'
#' @param Z Random-effects design matrix (`dgCMatrix` recommended).
#' @param var_eta Laplace posterior covariance, i.e. \eqn{\Sigma^{-1}} (dense numeric `q x q`).
#' @param temp_trc_C Vector (length n) of the third-derivative term on the linear predictor scale.
#' @param temp_trc_D Vector (length n) of the fourth-derivative term on the linear predictor scale.
#' @param max_nq_mem Soft cap on forming intermediate `Svar = Z %*% var_eta` in memory.
#'
#' @return A list with:
#' * `trc_y1`: length-`q` numeric vector (mean correction term),
#' * `trc_y2`: `q x q` sparse diagonal matrix (covariance correction term).
#'
#' @noRd
#' @keywords internal
fe_trace_diagG <- function(Z, var_eta, temp_trc_C, temp_trc_D, max_nq_mem = 5e7) {

  n <- nrow(Z)
  q <- ncol(Z)

  # For disjoint clusters the posterior Hessian is diagonal. These identities
  # are the same trace formulas below, reduced to scalar sums per cluster.
  # Retain the general path for crossed/multiple-membership designs.
  if (all(var_eta[row(var_eta) != col(var_eta)] == 0) &&
      all(Matrix::rowSums(Z != 0) <= 1)) {
    v <- diag(var_eta)
    c3 <- as.numeric(Matrix::crossprod(Z^3, temp_trc_C))
    c4 <- as.numeric(Matrix::crossprod(Z^4, temp_trc_D))
    return(list(trc_y1 = v^2 * c3,
                trc_y2 = Matrix::Diagonal(q, x = v^3 * c4 + 2 * v^4 * c3^2)))
  }

  # Memory guard: forming Svar = Z %*% var_eta is n×q dense.
  use_full_Svar <- (as.double(n) * as.double(q) <= as.double(max_nq_mem))

  trc_y1     <- numeric(q)
  trc_y2_diag <- numeric(q)

  # Precompute Svar if allowed: Svar = Z %*% var_eta (n×q).
  Svar <- if (use_full_Svar) as.matrix(Z %*% var_eta) else NULL

  for (j in seq_len(q)) {

    # Svar_j = (Z %*% var_eta)[,j]
    Svar_j <- if (use_full_Svar) Svar[, j] else as.numeric(Z %*% var_eta[, j])

    # ----- First derivative term: D_j = - dSigma/dc_j (q×q) -----
    # mvglmmRank: dsig.dc.j = t(Z) %*% (Z * Svar[,j] * temp_trc.C)
    w1 <- as.numeric(Svar_j * temp_trc_C)      # length n
    dsig_dc_j <- as.matrix(Matrix::crossprod(Z, Z * w1))  # q×q

    # trc_y1[j] = tr( var_eta %*% dsig_dc_j )
    trc_y1[j] <- sum(diag(var_eta %*% dsig_dc_j))

    # ----- Second derivative term (diagonal only): D2_jj = - d^2Sigma/dc_j^2 -----
    # mvglmmRank diag case:
    # qv    = var_eta %*% (dsig_dc_j %*% var_eta[,j])
    # d2sig = t(Z) %*% ( Z * ( (Z%*%qv)*temp_trc.C + Svar_j^2*temp_trc.D ) )
    qv <- var_eta %*% (dsig_dc_j %*% var_eta[, j])

    w2 <- as.numeric(as.numeric(Z %*% qv) * temp_trc_C + (Svar_j * Svar_j) * temp_trc_D)
    d2sig <- as.matrix(Matrix::crossprod(Z, Z * w2))      # q×q

    tr1 <- sum(diag(var_eta %*% d2sig))

    # IMPORTANT: mvglmmRank uses a Frobenius inner product here:
    # sum( (t(dsig)%*%var_eta) * (var_eta%*%dsig) )
    A <- t(dsig_dc_j) %*% var_eta
    B <- var_eta %*% dsig_dc_j
    tr2 <- glmmfe_frob(A, B)

    trc_y2_diag[j] <- tr1 + tr2
  }

  list(
    trc_y1 = trc_y1,
    trc_y2 = Matrix::Diagonal(q, x = trc_y2_diag)
  )
}



#' Binomial probit family object for glmmFEL
#'
#' @description
#' Constructs the family/derivative bundle used by glmmFEL for the
#' binomial-probit GLMM case, including:
#'
#' * stable derivatives up to 4th order for FE trace terms,
#' * score components used by the Laplace Newton update for eta,
#' * initial beta estimates and input checks,
#' * precomputed FE trace inputs for a given (beta, eta).
#'
#' @param y Response vector.
#' @param X Fixed-effects design matrix.
#' @param Z Random-effects design matrix.
#' @return A list of functions used internally by glmmFEL.
#' @noRd
#' @keywords internal
glmmfe_family_binomial_probit <- function(y, X, Z) {

  derivs <- fe_derivatives_binomial_probit()

  check_y <- function() {
    if (!all(y %in% c(0, 1))) stop("For binomial_probit, y must be 0/1.")
    invisible(TRUE)
  }

  init_beta <- function() {
    # Use a simple GLM init on fixed effects only (eta=0).
    # This is not the final estimator; it's a stable starting point.
    fit0 <- suppressWarnings(stats::glm.fit(x = X, y = y, family = stats::binomial(link = "probit")))
    b <- fit0$coefficients
    b[!is.finite(b)] <- 0
    b
  }

  # Return E, R2, R3 used in the Laplace eta Newton update and FE beta score.
  # For probit, define q = s * eta_lin with s = (-1)^(1-y).
  E_R2_R3 <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)
    s <- ifelse(y == 1, 1, -1)
    q <- s * eta_lin

    d1 <- derivs$d1(q)
    d2 <- derivs$d2(q)
    d3 <- derivs$d3(q)

    # E corresponds to d/deta loglik contribution, mapped back to eta_lin scale:
    # d/d(eta_lin) log Phi(s * eta_lin) = s * d1(q)
    E  <- s * d1
    R2 <- d2
    R3 <- s * d3

    list(E = E, R2 = R2, R3 = R3)
  }

  # FE trace input vectors temp_trc_C and temp_trc_D (third/fourth derivatives on eta_lin scale)
  FE_trace_inputs <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)
    s <- ifelse(y == 1, 1, -1)
    q <- s * eta_lin

    d3 <- derivs$d3(q)
    d4 <- derivs$d4(q)

    temp_trc_C <- s * d3
    temp_trc_D <- d4

    list(temp_trc_C = temp_trc_C, temp_trc_D = temp_trc_D)
  }

  list(
    check_y = check_y,
    init_beta = init_beta,
    E_R2_R3 = E_R2_R3,
    FE_trace_inputs = FE_trace_inputs
  )
}


#' Binomial logit family object for glmmFEL
#'
#' @description
#' Constructs the family/derivative bundle used by glmmFEL for the
#' binomial-logit GLMM case.
#'
#' @param y Response vector.
#' @param X Fixed-effects design matrix.
#' @param Z Random-effects design matrix.
#' @return A list of functions used internally by glmmFEL.
#' @noRd
#' @keywords internal
glmmfe_family_binomial_logit <- function(y, X, Z) {

  derivs <- fe_derivatives_binomial_logit()

  check_y <- function() {
    if (!all(y %in% c(0, 1))) stop("For binomial_logit, y must be 0/1.")
    invisible(TRUE)
  }

  init_beta <- function() {
    # Use a simple GLM init on fixed effects only (eta=0).
    fit0 <- suppressWarnings(stats::glm.fit(x = X, y = y, family = stats::binomial(link = "logit")))
    b <- fit0$coefficients
    b[!is.finite(b)] <- 0
    b
  }

  # For logit, define q = s * eta_lin with s = (-1)^(1-y).
  E_R2_R3 <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)
    s <- ifelse(y == 1, 1, -1)
    q <- s * eta_lin

    d1 <- derivs$d1(q)
    d2 <- derivs$d2(q)
    d3 <- derivs$d3(q)

    # d/d(eta_lin) log plogis(s*eta_lin) = s * d1(q)
    E  <- s * d1
    R2 <- d2
    R3 <- s * d3

    list(E = E, R2 = R2, R3 = R3)
  }

  FE_trace_inputs <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)
    s <- ifelse(y == 1, 1, -1)
    q <- s * eta_lin

    d3 <- derivs$d3(q)
    d4 <- derivs$d4(q)

    temp_trc_C <- s * d3
    temp_trc_D <- d4

    list(temp_trc_C = temp_trc_C, temp_trc_D = temp_trc_D)
  }

  list(
    check_y = check_y,
    init_beta = init_beta,
    E_R2_R3 = E_R2_R3,
    FE_trace_inputs = FE_trace_inputs
  )
}


#' Poisson log-link family object for glmmFEL
#'
#' @description
#' Constructs the family/derivative bundle used by glmmFEL for the
#' Poisson-log GLMM case.
#'
#' @param y Response vector.
#' @param X Fixed-effects design matrix.
#' @param Z Random-effects design matrix.
#' @return A list of functions used internally by glmmFEL.
#' @noRd
#' @keywords internal
glmmfe_family_poisson_log <- function(y, X, Z) {

  derivs <- fe_derivatives_poisson_log()

  check_y <- function() {
    if (any(y < 0) || any(!is.finite(y))) stop("For poisson_log, y must be finite and nonnegative.")
    invisible(TRUE)
  }

  init_beta <- function() {
    # GLM init on fixed effects only (eta=0).
    fit0 <- suppressWarnings(stats::glm.fit(x = X, y = y, family = stats::poisson(link = "log")))
    b <- fit0$coefficients
    b[!is.finite(b)] <- 0
    b
  }

  E_R2_R3 <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)

    d1 <- derivs$d1(eta_lin, y)
    d2 <- derivs$d2(eta_lin)
    d3 <- derivs$d3(eta_lin)

    E  <- d1
    R2 <- d2
    R3 <- d3

    list(E = E, R2 = R2, R3 = R3)
  }

  FE_trace_inputs <- function(beta, eta) {
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)

    temp_trc_C <- derivs$d3(eta_lin)
    temp_trc_D <- derivs$d4(eta_lin)

    list(temp_trc_C = temp_trc_C, temp_trc_D = temp_trc_D)
  }

  list(
    check_y = check_y,
    init_beta = init_beta,
    E_R2_R3 = E_R2_R3,
    FE_trace_inputs = FE_trace_inputs
  )
}


#' Create a family bundle for glmmFEL
#'
#' @description
#' Internal constructor that returns the set of derivative and helper functions
#' needed by the Laplace and fully exponential Laplace algorithms, based on a
#' canonical family label.
#'
#' @param fam_name Canonical family label (see [glmmfe_resolve_family()]).
#' @param y Response vector.
#' @param X Fixed-effects design matrix.
#' @param Z Random-effects design matrix.
#'
#' @return A list with the family-specific functions used by glmmFEL.
#' @noRd
#' @keywords internal
glmmfe_make_family <- function(fam_name, y, X, Z) {

  if (identical(fam_name, "binomial_probit")) return(glmmfe_family_binomial_probit(y, X, Z))
  if (identical(fam_name, "binomial_logit"))  return(glmmfe_family_binomial_logit(y, X, Z))
  if (identical(fam_name, "poisson_log"))     return(glmmfe_family_poisson_log(y, X, Z))

  stop("Unsupported family label in glmmfe_make_family(): ", fam_name)
}
