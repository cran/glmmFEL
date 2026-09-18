#' Pseudo-likelihood engine for RSPL/MSPL (Wolfinger-style, simplified R = I)
#'
#' This follows the vp_cp / RealVAMS structure:
#'  - Outer PQL loop updates working response (z_work) and weights (w_num).
#'  - Inner EM loop (with z_work, w_num fixed) updates (beta, eta, tau2).
#'
#' RSPL vs MSPL:
#'  #'  - MSPL uses \eqn{\mathrm{Var}(\eta \mid \beta, y) = (Z^\top W Z + G^{-1})^{-1}}
#'    (called \code{H.inv} in \code{vp_cp}).
#'  - RSPL uses the eta-eta block of the inverse of the full augmented system
#'    (called C.mat in vp_cp) for the variance-component moment update.
#'
#' @keywords internal
glmmFEL_pl <- function(
    y, X, Z,
    family   = "binomial_probit",
    approx   = c("RSPL", "MSPL"),
    max_iter = 200L,
    tol      = 1e-6,
    control  = list()
) {
  fam_name <- glmmfe_resolve_family(family)
  approx   <- glmmfe_resolve_approx(approx)
  if (!approx %in% c("RSPL", "MSPL")) stop("glmmFEL_pl only supports approx = RSPL or MSPL.")

  ## Inputs
  X <- glmmfe_as_X(X)
  Z <- glmmfe_as_Z(Z, n = length(y))
  glmmfe_validate_data(y, X, Z, fam_name)
  y <- as.numeric(y)

  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  if (q == 0L) stop("Z must have at least one random-effects column.")

  ## Basic response checks (match family assumptions)
  if (fam_name %in% c("binomial_probit", "binomial_logit")) {
    if (!all(y %in% c(0, 1))) stop("For binomial families, y must be 0/1.")
  }
  if (fam_name == "poisson_log") {
    if (any(y < 0) || any(!is.finite(y))) stop("For poisson_log, y must be finite and nonnegative.")
  }

  ## Base family object for link, mu.eta, variance
  fam_obj <- switch(
    fam_name,
    binomial_probit = stats::binomial(link = "probit"),
    binomial_logit  = stats::binomial(link = "logit"),
    poisson_log     = stats::poisson(link = "log"),
    stop("Unsupported family in glmmFEL_pl.")
  )

  ## Controls
  ctrl <- list(
    pql_max_iter = max_iter,
    pql_tol      = tol,

    em_max_iter  = 50L,
    em_tol       = tol,

    vc_eps       = 1e-12,
    tau2_init    = 1,

    ## NEW safe control for inner solve fallback:
    lmm_ridge_init = 1e-8,

    verbose      = FALSE,
    trace        = FALSE
  )
  ctrl <- glmmfe_validate_control(control, ctrl)

  ## Basic init: fixed-only GLM start (stable)
  beta <- rep(0, p)
  fit0 <- suppressWarnings(
    try(stats::glm.fit(x = X, y = y, family = fam_obj), silent = TRUE)
  )
  if (!inherits(fit0, "try-error")) {
    b0 <- as.numeric(fit0$coefficients)
    b0[!is.finite(b0)] <- 0
    if (length(b0) == p) beta <- b0
  }
  if (!is.null(colnames(X))) names(beta) <- colnames(X)
  glmmfe_check_separation(y, X, beta, fam_name)

  eta  <- rep(0, q)

  tau2 <- as.numeric(ctrl$tau2_init)
  if (!is.finite(tau2) || tau2 <= 0) tau2 <- 1
  tau2 <- max(tau2, ctrl$vc_eps)

  ## Optional trace storage
  pl_trace <- if (isTRUE(ctrl$trace)) numeric(0) else NULL

  ## Outer PQL loop
  pql_converged <- FALSE
  iter_used     <- 0L
  em_iters_last <- 0L

  inner_last <- NULL
  w_last <- NULL
  z_last <- NULL

  for (it in seq_len(ctrl$pql_max_iter)) {
    iter_used <- it
    beta_old_outer <- beta
    tau2_old_outer <- tau2
    eta_old_outer <- eta

    ## Linear predictor and mean
    eta_lin <- as.numeric(X %*% beta + Z %*% eta)
    mu      <- fam_obj$linkinv(eta_lin)

    ## Working response and weights
    dmu <- fam_obj$mu.eta(eta_lin)
    vmu <- fam_obj$variance(mu)

    dmu <- pmax(dmu, 1e-12)
    vmu <- pmax(vmu, 1e-12)

    w_num  <- as.numeric((dmu^2) / vmu)
    z_work <- as.numeric(eta_lin + (y - mu) / dmu)

    ## Inner EM loop with z_work, w_num fixed (Gaussian core)
    em_converged <- FALSE
    for (em_it in seq_len(ctrl$em_max_iter)) {
      em_iters_last <- em_it
      beta_old <- beta
      tau2_old <- tau2

      inner <- glmmfe_lmm_inner_fit(
        z_work = z_work, w_num = w_num,
        X = X, Z = Z, tau2 = tau2,
        approx = approx,
        vc_eps = ctrl$vc_eps,
        ridge_init = ctrl$lmm_ridge_init
      )

      beta <- as.numeric(inner$beta)
      eta  <- as.numeric(inner$eta)

      var_eta_used <- if (identical(approx, "RSPL")) inner$var_eta_reml else inner$var_eta_post

      tau2 <- mean(eta^2 + diag(var_eta_used))
      tau2 <- max(tau2, ctrl$vc_eps)

      theta_old <- c(beta_old, tau2_old)
      theta_new <- c(beta, tau2)

      num <- max(abs(theta_new - theta_old))
      den <- max(1, max(abs(theta_old)))
      delta <- num / den

      if (delta < ctrl$em_tol) {
        em_converged <- TRUE
        break
      }
    }

    inner_last <- inner
    w_last <- w_num
    z_last <- z_work

    if (isTRUE(ctrl$trace)) {
      pl_val <- glmmfe_pl_objective(
        z_work = z_work, w_num = w_num,
        X = X, Z = Z,
        beta = beta, eta = eta, tau2 = tau2,
        inner = inner_last,
        approx = approx,
        vc_eps = ctrl$vc_eps
      )
      pl_trace <- c(pl_trace, pl_val)
    }

    ## Outer convergence on (beta, tau2, eta), requiring the inner solve as well.
    theta_old_outer <- c(beta_old_outer, tau2_old_outer, eta_old_outer)
    theta_new_outer <- c(beta, tau2, eta)

    num <- max(abs(theta_new_outer - theta_old_outer))
    den <- max(1, max(abs(theta_old_outer)))
    delta_outer <- num / den

    if (isTRUE(ctrl$verbose)) {
      cat("PQL iter:", it,
          " approx:", approx,
          " delta_outer:", sprintf("%.3e", delta_outer),
          " tau2:", sprintf("%.6g", tau2),
          " em_converged:", em_converged,
          " em_iters:", em_iters_last, "\n")
    }

    if (delta_outer < ctrl$pql_tol && em_converged) {
      pql_converged <- TRUE
      break
    }
  }

  # Refresh working quantities and covariance at the returned parameter values.
  final_linear <- as.numeric(X %*% beta + Z %*% eta)
  mu <- fam_obj$linkinv(final_linear)
  dmu <- pmax(fam_obj$mu.eta(final_linear), 1e-12)
  w_last <- dmu^2 / pmax(fam_obj$variance(mu), 1e-12)
  z_last <- final_linear + (y - mu) / dmu
  inner_last <- glmmfe_lmm_inner_fit(z_last, w_last, X, Z, tau2, approx,
                                    ctrl$vc_eps, ctrl$lmm_ridge_init)

  ## Final covariance outputs
  vcov_beta    <- if (!is.null(inner_last)) inner_last$vcov_beta else NULL
  cov_beta_eta <- if (!is.null(inner_last)) inner_last$cov_beta_eta else NULL

  var_eta_post <- if (!is.null(inner_last)) inner_last$var_eta_post else diag(NA_real_, q)

  var_eta_used <- if (!is.null(inner_last)) {
    if (identical(approx, "RSPL")) inner_last$var_eta_reml else inner_last$var_eta_post
  } else {
    var_eta_post
  }

  ## Retain the working Gaussian objective separately from a GLMM likelihood.
  logLik_val <- NA_real_
  if (!is.null(inner_last) && !is.null(w_last) && !is.null(z_last)) {
    logLik_val <- glmmfe_pl_objective(
      z_work = z_last, w_num = w_last,
      X = X, Z = Z,
      beta = beta, eta = eta, tau2 = tau2,
      inner = inner_last,
      approx = approx,
      vc_eps = ctrl$vc_eps
    )
  }

  ## Bookkeeping
  G <- Matrix::Diagonal(q, x = rep.int(tau2, q))

  fit <- glmmfe_new_fit(
    y = y, X = X, Z = Z,
    beta = beta, eta = eta, tau2 = tau2,
    G = G,
    vcov_beta = vcov_beta,
    vcov_eta  = var_eta_used,
    cov_beta_eta = cov_beta_eta,
    var_eta = var_eta_used,
    family = fam_name,
    approx = approx,
    control = ctrl,
    convergence = list(
      pql_converged = pql_converged,
      pql_iter = iter_used,
      em_iter_last = em_iters_last,
      em_converged = em_converged,
      reason = if (pql_converged) "converged" else "iteration_limit",
      pl_trace = pl_trace
    ),
    logLik = NA_real_,
    call = match.call(),
    reml = identical(approx, "RSPL")
  )

  fit$var_eta_post <- var_eta_post
  fit$var_eta_used <- var_eta_used
  fit$working_logLik <- logLik_val
  fit$logLik_type <- "working Gaussian objective only; GLMM likelihood unavailable"
  if (!pql_converged) warning("glmmFEL pseudo-likelihood did not converge: iteration_limit", call. = FALSE)

  fit
}
