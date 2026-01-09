#' Internal Gaussian inner fit for PL / weighted LMM with G = tau2 * I
#'
#' Model (conditional on working quantities):
#' \deqn{z_{\mathrm{work}} = X\beta + Z\eta + e,\qquad e \sim N(0, W^{-1})}
#' \deqn{\eta \sim N(0, \tau^2 I_q)}
#'
#' where \eqn{W = \mathrm{diag}(w_{\mathrm{num}})}.
#'
#' Returns \code{beta}, \code{eta} and covariance blocks needed by RSPL/MSPL updates.
#' Does NOT form \eqn{n \times n} matrices.
#'
#' @keywords internal
glmmfe_lmm_inner_fit <- function(
    z_work, w_num, X, Z, tau2,
    approx = c("RSPL", "MSPL"),
    vc_eps = 1e-12,
    ridge_init = 1e-8
) {
  approx <- glmmfe_resolve_approx(approx)
  if (!approx %in% c("RSPL", "MSPL")) stop("glmmfe_lmm_inner_fit only supports approx = RSPL or MSPL.")

  z_work <- as.numeric(z_work)
  w_num  <- as.numeric(w_num)

  X <- glmmfe_as_X(X)
  Z <- glmmfe_as_Z(Z, n = length(z_work))

  n <- length(z_work)
  p <- ncol(X)
  q <- ncol(Z)
  if (length(w_num) != n) stop("w_num must have length n.")

  inv_tau2 <- 1 / max(as.numeric(tau2), vc_eps)

  ## Build A = X' W X,  B = X' W Z  (small matrices)
  A <- as.matrix(crossprod(X, X * w_num))               # p x p
  B <- as.matrix(Matrix::crossprod(X, Z * w_num))       # p x q

  ## Build C_eta = Z' W Z + inv_tau2 I  (q x q)
  C_eta <- Matrix::crossprod(Z, Z * w_num) +
    Matrix::Diagonal(q, x = rep.int(inv_tau2, q))
  C_eta_mat <- as.matrix(C_eta)

  rhs1 <- as.numeric(crossprod(X, z_work * w_num))            # p
  rhs2 <- as.numeric(Matrix::crossprod(Z, z_work * w_num))    # q

  ## Full augmented system (vp_cp: C_inv) is (p+q) x (p+q)
  H_aug <- rbind(
    cbind(A, B),
    cbind(t(B), C_eta_mat)
  )
  H_aug <- 0.5 * (H_aug + t(H_aug))

  ## Solve for [beta; eta] with a small ridge fallback on A
  sol <- tryCatch(solve(H_aug, c(rhs1, rhs2)), error = function(e) NULL)
  if (is.null(sol)) {
    A2 <- A + diag(ridge_init, p)
    H2 <- rbind(cbind(A2, B), cbind(t(B), C_eta_mat))
    H2 <- 0.5 * (H2 + t(H2))
    sol <- solve(H2, c(rhs1, rhs2))
    H_aug <- H2
  }

  beta_new <- sol[1:p]
  eta_new  <- sol[(p + 1):(p + q)]

  ## var_eta_post = (Z'WZ + inv_tau2 I)^{-1}  (vp_cp: H.inv)
  C_sym <- 0.5 * (C_eta_mat + t(C_eta_mat))
  cholC <- tryCatch(chol(C_sym), error = function(e) NULL)
  var_eta_post <- if (!is.null(cholC)) chol2inv(cholC) else solve(C_sym)

  ## Full inverse for REML covariance blocks (vp_cp: cs = inverse(C_inv))
  cholH <- tryCatch(chol(H_aug), error = function(e) NULL)
  H_inv <- if (!is.null(cholH)) chol2inv(cholH) else solve(H_aug)

  vcov_beta    <- H_inv[1:p, 1:p, drop = FALSE]
  cov_beta_eta <- H_inv[1:p, (p + 1):(p + q), drop = FALSE]
  var_eta_reml <- H_inv[(p + 1):(p + q), (p + 1):(p + q), drop = FALSE]

  ## logdet(chol(.)) terms in vp_cp convention:
  ## - REML subtracts logdet(chol(C_inv))  -> here chol(H_aug)
  ## - ML   subtracts logdet(chol(H))      -> here chol(C_sym)
  logdet_chol_eta <- if (!is.null(cholC)) sum(log(diag(cholC))) else NA_real_
  logdet_chol_aug <- if (!is.null(cholH)) sum(log(diag(cholH))) else NA_real_

  list(
    beta = beta_new,
    eta  = eta_new,
    var_eta_post = var_eta_post,
    var_eta_reml = var_eta_reml,
    vcov_beta = vcov_beta,
    cov_beta_eta = cov_beta_eta,
    logdet_chol_eta = logdet_chol_eta,
    logdet_chol_aug = logdet_chol_aug
  )
}


#' vp_cp-style PL objective (includes constants)
#'
#' @keywords internal
glmmfe_pl_objective <- function(z_work, w_num, X, Z, beta, eta, tau2, inner, approx, vc_eps = 1e-12) {
  approx <- glmmfe_resolve_approx(approx)
  if (!approx %in% c("RSPL", "MSPL")) stop("glmmfe_pl_objective only supports approx = RSPL or MSPL.")

  z_work <- as.numeric(z_work)
  w_num  <- as.numeric(w_num)
  X <- glmmfe_as_X(X)
  Z <- glmmfe_as_Z(Z, n = length(z_work))

  n <- length(z_work)
  p <- ncol(X)
  q <- ncol(Z)

  cons.logLik <- if (identical(approx, "RSPL")) {
    0.5 * (q + p) * log(2 * pi)
  } else {
    0.5 * q * log(2 * pi)
  }

  resid <- as.numeric(z_work - (X %*% beta + Z %*% eta))

  log_p_y <- -(n / 2) * log(2 * pi) +
    0.5 * sum(log(w_num)) -
    0.5 * sum((resid^2) * w_num)

  log_p_eta <- -(q / 2) * log(2 * pi) -
    (q / 2) * log(max(tau2, vc_eps)) -
    0.5 * sum(eta^2) / max(tau2, vc_eps)

  det_term <- if (identical(approx, "RSPL")) inner$logdet_chol_aug else inner$logdet_chol_eta
  if (!is.finite(det_term)) return(NA_real_)

  as.numeric(cons.logLik + log_p_eta + log_p_y - det_term)
}
