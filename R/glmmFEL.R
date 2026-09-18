#' Fit GLMMs via Laplace and fully exponential Laplace (matrix interface)
#'
#' @description
#' `glmmFEL()` fits a generalized linear mixed model (GLMM) with multivariate
#' normal random effects using EM-type algorithms and likelihood approximations:
#'
#' * first-order Laplace (`approx = "Laplace"`),
#' * fully exponential corrections to the random-effects mean
#'   (`approx = "FE_mean"`),
#' * fully exponential corrections to both mean and variance diagonals
#'   (`approx = "FE_full"` / `"FE"`),
#' * pseudo-likelihood / PL linearization (`approx = "RSPL"` or `"MSPL"`)
#'   via [glmmFEL_pl()].
#'
#' The interface is matrix-based: provide the response `y`, fixed-effects
#' design matrix `X`, and random-effects design matrix `Z`.
#'
#' Random effects are assumed \eqn{\eta \sim N(0, G)} with a **single** variance
#' component
#' \deqn{G = \tau^2 I_q,}
#' allowing arbitrary (including multi-membership) `Z` while keeping the variance
#' update simple and stable.
#'
#' @name glmmFEL
#'
#' @param y
#'  Numeric response vector of length \eqn{n}. For
#'  `family = "binomial_probit"` / `binomial(link = "probit")` or
#'  `family = "binomial_logit"` / `binomial(link = "logit")`,
#'  values must be 0 or 1; grouped binomial counts and trial weights are not
#'  supported. Poisson responses must be nonnegative integer counts.
#' @param X
#'  Fixed-effects design matrix of dimension \eqn{n \times p}. May be a
#'  base R matrix or a matrix-like object; it is internally coerced to a
#'  base numeric matrix. Must have at least one column, finite entries, and
#'  full column rank.
#' @param Z
#'  Random-effects design matrix of dimension \eqn{n \times q}. May be a
#'  base R matrix or a \pkg{Matrix} object. Internally it is coerced to a
#'  sparse \code{"dgCMatrix"} where possible (to preserve sparsity).
#'  Must have at least one column (no purely fixed-effects models).
#' @param family
#'  Either a character string or a [stats::family] object indicating the
#'  model family. The argument is resolved via [glmmfe_resolve_family()].
#' @param approx
#'  Approximation type, resolved via [glmmfe_resolve_approx()]. Accepted
#'  values (case-insensitive) include:
#'
#'  * `"Laplace"` – first-order Laplace approximation,
#'  * `"FE_mean"` – staged algorithm: Laplace phase then FE mean corrections,
#'  * `"FE"` / `"FE_full"` – staged algorithm: Laplace phase, then FE mean,
#'    then FE covariance corrections (default),
#'  * `"RSPL"` – restricted pseudo-likelihood (REML-style) linearization,
#'  * `"MSPL"` – marginal pseudo-likelihood (ML-style) linearization.
#' @param max_iter
#'  Maximum number of EM iterations (outer iterations over \eqn{\beta}
#'  and \eqn{\tau^2}). Can be overridden by `control$em_max_iter`.
#' @param tol
#'  Baseline convergence tolerance for the EM algorithm. The staged thresholds
#'  default to:
#'  \itemize{
#'    \item Laplace stage: `tol_laplace = 10 * tol`,
#'    \item FE-mean stage: `tol_fe_mean = 3 * tol`,
#'    \item FE-full stage: `tol_fe_full = tol`.
#'  }
#'  You can override these via `control$tol_laplace`, `control$tol_fe_mean`,
#'  and `control$tol_fe_full`.
#' @param control
#'  List of optional control parameters. Recognized entries include:
#'  \itemize{
#'    \item `em_max_iter`, `em_tol`,
#'    \item `tol_laplace`, `tol_fe_mean`, `tol_fe_full`,
#'    \item `eta_max_iter`, `eta_tol_grad`,
#'    \item `beta_max_iter`, `beta_tol`,
#'    \item `tau2_init` (initial value for \eqn{\tau^2}),
#'    \item `vc_eps` (lower bound for \eqn{\tau^2}),
#'    \item `max_nq_mem` (memory guard for FE trace intermediates),
#'    \item `verbose` (logical),
#'    \item `beta_step_max` (max Newton step size for beta; default 2),
#'    \item `beta_ls_max_iter` (max line-search halvings; default 12),
#'    \item `beta_hess_ridge_init` (initial ridge for Hessian; default 1e-8),
#'    \item `beta_hess_ridge_max` (max ridge; default 1e2)
#'  }
#'
#' @details
#' `Laplace` denotes an EM algorithm using posterior modes in its expected
#' score and Gaussian posterior second moments. It is not direct maximization
#' of the first-order Laplace marginal likelihood (as in `lme4::glmer`).
#' `FE_mean` adds mean and expected-score corrections. `FE_full` also corrects
#' the posterior variance diagonals needed to update the single variance
#' component; off-diagonal entries of `var_eta` remain at their Laplace values.
#' Thus "full" distinguishes the two implemented EM stages, not a correction
#' of every entry of an arbitrary posterior covariance matrix.
#'
#' Returned moments are recomputed at the final parameter estimates. A failed
#' inner solve or exhausted iteration budget produces a warning and a false
#' convergence flag. Always inspect `convergence` before interpreting estimates.
#' Approximate posterior variances can become invalid in sparse/extreme settings;
#' such fits are flagged rather than reported as converged.
#' A strictly separating fixed-effect direction found at initialization is
#' rejected for binary data, since no finite estimate exists. This sufficient
#' check is not a complete detector of all forms of quasi-separation.
#'
#' `vcov_beta` is a joint-Hessian approximation conditional on the fitted
#' variance component; it is not a fully exponential observed-information
#' estimate and does not incorporate variance-component estimation uncertainty.
#' For EM fits, `logLik` evaluates the first-order Laplace marginal likelihood
#' at the returned estimates, which need not maximize it. For PL fits it is
#' unavailable (`NA`); `working_logLik` stores the working Gaussian objective.
#' These objectives should not be used for cross-method AIC comparisons.
#' Fitted response values plug the estimated random effects into the inverse
#' link; they are not posterior predictive averages over random effects.
#'
#' See [glmmFEL-benchmarks] for a reproducible simulation appraisal and its limits.
#'
#' @return
#' A fitted model object of class `glmmFELMod`, with estimates `beta`, `tau2`,
#' random-effect predictions `eta`, approximate covariance matrices, and
#' `convergence`. EM convergence diagnostics include `em_converged`, `reason`,
#' `mode_converged`, `mode_gradient`, `beta_converged`, and `moments_valid`.
#' PL fits report both `pql_converged` and `em_converged`.
#'
#' @template ref-doc
#' @examples
#' ## Example 1: Simulated probit random-intercept GLMM (matrix interface)
#' set.seed(1)
#' n_id <- 30
#' m_per_id <- 6
#' n <- n_id * m_per_id
#' id <- factor(rep(seq_len(n_id), each = m_per_id))
#' x  <- rnorm(n)
#' X  <- model.matrix(~ x)
#' Z  <- Matrix::sparseMatrix(i = seq_len(n),
#'                            j = as.integer(id),
#'                            x = 1,
#'                            dims = c(n, n_id))
#' beta_true <- c(0.2, 0.7)
#' tau2_true <- 0.5
#' eta_true  <- rnorm(n_id, sd = sqrt(tau2_true))
#' lp <- as.vector(X %*% beta_true + Z %*% eta_true)
#' y  <- rbinom(n, 1, pnorm(lp))
#'
#' fit <- glmmFEL(y, X, Z, family = "binomial_probit", approx = "Laplace")
#' fit$beta
#' fit$tau2
#'
#' ## Example 2: Get X, y, Z from an lme4 formula (without a glmmFEL formula wrapper)
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   dat <- data.frame(y = y, x = x, id = id)
#'   lf  <- lme4::lFormula(y ~ x + (1 | id), data = dat)
#'   X_lme4 <- lf$X
#'   Z_lme4 <- Matrix::t(lf$reTrms$Zt)
#'   y_lme4 <- lf$fr$y
#'
#'   fit2 <- glmmFEL(y_lme4, X_lme4, Z_lme4, family = "binomial_probit", approx = "Laplace")
#' }
#' }
#'
#' @export
glmmFEL <- function(
    y,
    X,
    Z,
    family   = stats::binomial(link = "probit"),
    approx   = c("FE", "Laplace", "FE_mean", "FE_full", "RSPL", "MSPL"),
    max_iter = 200,
    tol      = 1e-6,
    control  = list()
) {
  ## -----------------------------
  ## Resolve family / approx
  ## -----------------------------
  fam_name   <- glmmfe_resolve_family(family)
  approx_lab <- glmmfe_resolve_approx(approx)

  ## -----------------------------
  ## Basic argument processing
  ## -----------------------------
  X <- glmmfe_as_X(X)
  Z_M <- glmmfe_as_Z(Z, n = length(y))
  glmmfe_validate_data(y, X, Z_M, fam_name)
  y <- as.numeric(y)
  n <- length(y); p <- ncol(X); q <- ncol(Z_M)

  ## -----------------------------
  ## Dispatch to pseudo-likelihood engine (RSPL/MSPL)
  ## -----------------------------
  if (approx_lab %in% c("RSPL", "MSPL")) {
    fit <- glmmFEL_pl(
      y        = y,
      X        = X,
      Z        = Z_M,
      family   = fam_name,
      approx   = approx_lab,
      max_iter = max_iter,
      tol      = tol,
      control  = control
    )
    fit$call <- match.call()
    return(fit)
  }

  ## -----------------------------
  ## Control defaults
  ## -----------------------------
  ctrl <- list(
    em_max_iter   = max_iter,
    em_tol        = tol,
    tol_laplace   = NULL,
    tol_fe_mean   = NULL,
    tol_fe_full   = NULL,

    eta_max_iter  = 50L,
    eta_tol_grad  = tol,

    beta_max_iter = 50L,
    beta_tol      = tol,

    ## beta Newton robustness controls
    beta_step_max        = 2.0,     # cap on ||step||_inf after solving
    beta_ls_max_iter     = 12L,     # backtracking halvings
    beta_hess_ridge_init = 1e-8,    # initial ridge for Hessian; fallback only
    beta_hess_ridge_max  = 1e2,     # max ridge allowed

    tau2_init     = 1,
    vc_eps        = 1e-12,

    ## Memory guard for FE trace corrections and FE beta intermediates:
    ## if n*q is larger than this, helpers avoid forming dense n×q intermediates.
    max_nq_mem    = 2.5e7,

    verbose       = FALSE
  )
  ctrl <- glmmfe_validate_control(control, ctrl)
  if (is.null(ctrl$tol_laplace)) ctrl$tol_laplace <- 10 * ctrl$em_tol
  if (is.null(ctrl$tol_fe_mean)) ctrl$tol_fe_mean <- 3 * ctrl$em_tol
  if (is.null(ctrl$tol_fe_full)) ctrl$tol_fe_full <- ctrl$em_tol

  ## -----------------------------
  ## Family spec (derivatives & checks)
  ## -----------------------------
  fam_spec <- glmmfe_make_family(fam_name, y, X, Z_M)
  fam_spec$check_y()

  ## -----------------------------
  ## Initial values
  ## -----------------------------
  beta <- fam_spec$init_beta()
  glmmfe_check_separation(y, X, beta, fam_name)
  eta  <- rep(0, q)

  tau2 <- as.numeric(ctrl$tau2_init)
  if (!is.finite(tau2) || tau2 <= 0) tau2 <- 1
  tau2 <- max(tau2, ctrl$vc_eps)

  ## Initialize Var(eta|y) (overwritten after the first Laplace step)
  G       <- Matrix::Diagonal(q, x = rep.int(tau2, q))
  var_eta <- as.matrix(G)

  update_eta_laplace <- function(eta_init, beta, tau2) {
    glmmfe_mode(eta_init, beta, tau2, y, X, Z_M, fam_name, fam_spec, ctrl)
  }
  moments <- function(lap, beta, phase) {
    mean <- lap$eta; variance <- lap$var_eta
    if (phase > 1L && all(is.finite(variance))) {
      d <- fam_spec$FE_trace_inputs(beta, lap$eta)
      fe <- fe_trace_diagG(Z_M, variance, d$temp_trc_C, d$temp_trc_D,
                          max_nq_mem = ctrl$max_nq_mem)
      mean <- mean + 0.5 * fe$trc_y1
      if (phase == 3L) diag(variance) <- diag(variance) + 0.5 * Matrix::diag(fe$trc_y2)
    }
    list(mean = mean, variance = variance,
         valid = all(is.finite(mean)) && all(is.finite(variance)) && all(diag(variance) >= 0))
  }

  ## -----------------------------
  ## EM loop with staged approximation
  ## -----------------------------
  em_converged <- FALSE
  iter_used    <- 0L
  phase        <- 1L
  phase_used <- 1L
  failure_reason <- "iteration_limit"
  beta_info <- list(converged = FALSE, score = Inf, reason = "not_run")

  for (iter in seq_len(ctrl$em_max_iter)) {
    iter_used <- iter
    phase_used <- phase
    beta_old  <- beta
    tau2_old  <- tau2

    ## Laplace step: mode and covariance of eta given current beta,tau2
    lap <- update_eta_laplace(eta_init = eta, beta = beta, tau2 = tau2)
    eta0     <- lap$eta
    var_eta0 <- lap$var_eta

    if (!lap$converged) { failure_reason <- lap$reason; break }
    mom <- moments(lap, beta, phase)
    if (!mom$valid) { failure_reason <- "invalid_posterior_moments"; break }
    eta_hat <- mom$mean; var_eta_hat <- mom$variance

    ## Update tau2 using FE-corrected moments (or Laplace moments in phase 1):
    tau2 <- mean(diag(var_eta_hat) + eta_hat^2)
    tau2 <- max(tau2, ctrl$vc_eps)
    G    <- Matrix::Diagonal(q, x = rep.int(tau2, q))

    ## Update beta (robust Newton)
    beta_info <- glmmfe_beta_update(beta, eta0, var_eta0, phase, X, Z_M, fam_spec, ctrl)
    beta <- beta_info$beta
    if (!beta_info$converged) { failure_reason <- beta_info$reason; break }

    ## carry forward FE-corrected (or Laplace) eta moments
    eta     <- eta0  # Warm-start the next mode at the previous mode, not its FE shift.
    var_eta <- var_eta_hat

    ## convergence check on (beta, tau2)
    theta_old_full <- c(beta_old, tau2_old)
    theta_new_full <- c(beta, tau2)

    num   <- max(abs(theta_new_full - theta_old_full))
    den   <- max(1, max(abs(theta_old_full)))
    delta <- num / den

    if (isTRUE(ctrl$verbose)) {
      cat("EM iter:", iter,
          " phase:", phase,
          " approx:", approx_lab,
          " delta:", sprintf("%.3e", delta),
          " tau2:", sprintf("%.6g", tau2), "\n")
    }

    ## phase transitions
    if (phase == 1L) {
      if (delta < ctrl$tol_laplace) {
        if (approx_lab == "Laplace") {
          em_converged <- TRUE
          break
        }
        phase <- 2L
      }
    } else if (phase == 2L) {
      if (approx_lab == "FE_mean") {
        if (delta < ctrl$tol_fe_mean) {
          em_converged <- TRUE
          break
        }
      } else {
        if (delta < ctrl$tol_fe_mean) phase <- 3L
      }
    } else if (phase == 3L) {
      if (delta < ctrl$tol_fe_full) {
        em_converged <- TRUE
        break
      }
    }
  }

  ## -----------------------------
  ## Final Hessian-based covariance (Laplace-style) computed at the MODE
  ## even if eta stored in the fit is FE-corrected.
  ##
  ## CHANGE: Hessian not PD is now a SOFT FAILURE (vcovs NA, fit still returned).
  ## -----------------------------
  lap_final <- update_eta_laplace(eta_init = eta, beta = beta, tau2 = tau2)
  eta_mode  <- lap_final$eta
  final_mom <- moments(lap_final, beta, phase_used)
  eta <- final_mom$mean
  var_eta <- final_mom$variance
  em_converged <- em_converged && lap_final$converged && final_mom$valid
  if (em_converged) failure_reason <- "converged"
  if (!lap_final$converged) failure_reason <- lap_final$reason
  if (!final_mom$valid) failure_reason <- "invalid_posterior_moments"

  comp  <- fam_spec$E_R2_R3(beta, eta_mode)
  R2    <- comp$R2
  temp  <- -R2

  inv_tau2 <- 1 / max(tau2, ctrl$vc_eps)
  G_inv    <- diag(inv_tau2, q)
  logdet_G <- q * log(max(tau2, ctrl$vc_eps))

  H11 <- crossprod(X, X * temp)
  H12 <- Matrix::crossprod(X, Z_M * temp)
  H21 <- t(as.matrix(H12))
  H22 <- G_inv + as.matrix(Matrix::crossprod(Z_M, Z_M * temp))

  chol_H22 <- tryCatch(chol(H22), error = function(e) NULL)
  logdet_H <- if (is.null(chol_H22)) NA_real_ else 2 * sum(log(diag(chol_H22)))

  H_joint <- rbind(
    cbind(H11, as.matrix(H12)),
    cbind(H21, H22)
  )
  H_joint <- 0.5 * (H_joint + t(H_joint))

  chol_Hjoint <- tryCatch(chol(H_joint), error = function(e) e)
  cov_pd <- !inherits(chol_Hjoint, "error")
  cov_err <- if (!cov_pd) conditionMessage(chol_Hjoint) else NA_character_
  if (!cov_pd) {
    em_converged <- FALSE
    failure_reason <- "fixed_effects_information_not_positive_definite"
  }

  if (cov_pd) {
    vcov_joint <- chol2inv(chol_Hjoint)
  } else {
    vcov_joint <- matrix(NA_real_, nrow = p + q, ncol = p + q)
  }

  vcov_beta <- vcov_joint[1:p, 1:p, drop = FALSE]
  vcov_eta  <- vcov_joint[(p + 1):(p + q), (p + 1):(p + q), drop = FALSE]
  cov_beta_eta_block <- vcov_joint[1:p, (p + 1):(p + q), drop = FALSE]

  se_beta <- sqrt(diag(vcov_beta))
  se_eta  <- sqrt(diag(vcov_eta))

  convergence <- list(
    em_converged = em_converged,
    em_iter      = iter_used,
    phase        = phase_used,
    reason       = failure_reason,
    mode_converged = lap_final$converged,
    mode_gradient = lap_final$gradient,
    beta_converged = beta_info$converged,
    beta_score = beta_info$score,
    moments_valid = final_mom$valid,
    tol_laplace  = ctrl$tol_laplace,
    tol_fe_mean  = ctrl$tol_fe_mean,
    tol_fe_full  = ctrl$tol_fe_full,
    cov_pd       = cov_pd,
    cov_err      = cov_err
  )

  ## Laplace-approximate marginal log-likelihood at the final parameter values
  logLik_val <- NA_real_
  if (is.finite(logdet_G) && is.finite(logdet_H)) {
    eta_lin <- as.numeric(X %*% beta + Z_M %*% eta_mode)

    if (fam_name == "binomial_probit") {
      y1 <- (y == 1)
      ll_y <- sum(stats::pnorm(eta_lin[y1], log.p = TRUE)) +
        sum(stats::pnorm(eta_lin[!y1], lower.tail = FALSE, log.p = TRUE))
    } else if (fam_name == "binomial_logit") {
      y1 <- (y == 1)
      ll_y <- sum(stats::plogis(eta_lin[y1], log.p = TRUE)) +
        sum(stats::plogis(eta_lin[!y1], lower.tail = FALSE, log.p = TRUE))
    } else if (fam_name == "poisson_log") {
      mu   <- exp(eta_lin)
      ll_y <- sum(y * eta_lin - mu - lgamma(y + 1))
    } else {
      ll_y <- NA_real_
    }

    if (is.finite(ll_y)) {
      quad <- sum(eta_mode^2) * inv_tau2
      logLik_val <- ll_y - 0.5 * (logdet_G + quad + logdet_H)
    }
  }

  fit <- glmmfe_new_fit(
    y            = y,
    X            = X,
    Z            = Z_M,
    beta         = beta,
    eta          = eta,        # FE-corrected mean (or mode if Laplace)
    tau2         = tau2,
    G            = G,
    vcov_beta    = vcov_beta,
    vcov_eta     = vcov_eta,
    cov_beta_eta = cov_beta_eta_block,
    var_eta      = var_eta,    # EM/FE prediction-error covariance
    family       = fam_name,
    approx       = approx_lab,
    control      = ctrl,
    convergence  = convergence,
    logLik       = logLik_val,
    call         = match.call(),
    reml         = NULL
  )

  ## extras (safe additions)
  fit$vcov_joint <- vcov_joint
  fit$se_beta    <- se_beta
  fit$se_eta     <- se_eta

  ## store Laplace mode used for covariance/logLik
  fit$eta_mode     <- eta_mode
  fit$var_eta_mode <- lap_final$var_eta
  fit$logLik_type <- "Laplace marginal log-likelihood evaluated at the EM estimates"
  if (!em_converged) warning("glmmFEL did not converge: ", failure_reason, call. = FALSE)

  fit
}
