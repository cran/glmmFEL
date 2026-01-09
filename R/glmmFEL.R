#' Fit GLMMs via Laplace and fully exponential Laplace (matrix interface)
#'
#' @description
#' `glmmFEL()` fits a generalized linear mixed model (GLMM) with multivariate
#' normal random effects using EM-type algorithms and likelihood approximations:
#'
#' * first-order Laplace (`approx = "Laplace"`),
#' * fully exponential corrections to the random-effects mean
#'   (`approx = "FE_mean"`),
#' * fully exponential corrections to both mean and covariance
#'   (`approx = "FE_full"` / `"FE"`),
#' * pseudo-likelihood / PL linearization (`approx = "RSPL"` or `"MSPL"`)
#'   via [glmmFEL_pl()].
#'
#' This **development branch is matrix-only**: you provide the response `y`,
#' fixed-effects design matrix `X`, and random-effects design matrix `Z`.
#' A formula interface (via optional 'lme4' helpers) and structured \eqn{G}
#' parameterizations are in development.
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
#'  values must be 0 or 1.
#' @param X
#'  Fixed-effects design matrix of dimension \eqn{n \times p}. May be a
#'  base R matrix or a matrix-like object; it is internally coerced to a
#'  base numeric matrix. Must have full column rank.
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
#' @return
#'  A fitted model object of class `glmmFELMod`.
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
  y   <- as.numeric(y)
  X   <- glmmfe_as_X(X)
  Z_M <- glmmfe_as_Z(Z, n = length(y))

  n <- length(y)
  if (nrow(X) != n || nrow(Z_M) != n) {
    stop("y, X, and Z must have the same number of rows/observations.")
  }

  p <- ncol(X)
  q <- ncol(Z_M)
  if (q == 0L) stop("Z must have at least one random-effects column.")

  rankX <- qr(X)$rank
  if (rankX < p) stop("Fixed-effects design matrix X is not full rank; remove collinear columns from X.")

  ## -----------------------------
  ## Dispatch to pseudo-likelihood engine (RSPL/MSPL)
  ## -----------------------------
  if (approx_lab %in% c("RSPL", "MSPL")) {
    return(glmmFEL_pl(
      y        = y,
      X        = X,
      Z        = Z_M,
      family   = fam_name,
      approx   = approx_lab,
      max_iter = max_iter,
      tol      = tol,
      control  = control
    ))
  }

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for glmmFEL(). Please install it.")
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
  if (length(control) > 0L) {
    for (nm in names(control)) ctrl[[nm]] <- control[[nm]]
  }
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
  eta  <- rep(0, q)

  tau2 <- as.numeric(ctrl$tau2_init)
  if (!is.finite(tau2) || tau2 <= 0) tau2 <- 1
  tau2 <- max(tau2, ctrl$vc_eps)

  ## Initialize Var(eta|y) (overwritten after the first Laplace step)
  G       <- Matrix::Diagonal(q, x = rep.int(tau2, q))
  var_eta <- as.matrix(G)

  ## -----------------------------
  ## Newton helper for eta (Laplace mode + Laplace var_eta = Sigma^{-1})
  ##
  ## Paper mapping (Karl et al.):
  ##   L(eta)   is the gradient of NEG log-posterior:
  ##     L(eta) = -Z' * (d/deta log f(y|eta)) + G^{-1} eta
  ##   Sigma    is the NEG Hessian at the mode.
  ## -----------------------------
  update_eta_laplace <- function(eta_init, beta, tau2) {
    eta <- as.numeric(eta_init)

    inv_tau2 <- 1 / max(tau2, ctrl$vc_eps)
    G_inv    <- diag(inv_tau2, q)

    var_eta_loc <- diag(1, q)

    for (it in seq_len(ctrl$eta_max_iter)) {
      comp <- fam_spec$E_R2_R3(beta, eta)
      E    <- comp$E
      R2   <- comp$R2

      ## grad of NEG log-posterior:  grad = -Z' E + G^{-1} eta
      grad_r <- as.numeric(Matrix::crossprod(Z_M, E))
      grad_p <- - as.numeric(G_inv %*% eta)
      grad   <- - (grad_r + grad_p)

      ## Sigma = G^{-1} + Z' diag(-R2) Z
      temp <- -R2
      H <- G_inv + as.matrix(Matrix::crossprod(Z_M, Z_M * temp))
      H <- 0.5 * (H + t(H))

      H_chol <- tryCatch(chol(H), error = function(e) NULL)
      if (!is.null(H_chol)) {
        var_eta_loc <- chol2inv(H_chol)
      } else {
        var_eta_loc <- tryCatch(solve(H), error = function(e) diag(1, q))
      }

      step <- var_eta_loc %*% grad

      ## Newton convergence criterion aligned with grad' * step
      crit <- as.numeric(crossprod(grad, step))
      if (!is.finite(crit) || crit <= ctrl$eta_tol_grad) break

      eta <- eta - as.numeric(step)
    }

    var_eta_loc <- 0.5 * (var_eta_loc + t(var_eta_loc))
    list(eta = eta, var_eta = var_eta_loc)
  }

  ## -----------------------------
  ## Newton helper for beta (robust + sparse-friendly)
  ##
  ## Safe changes implemented:
  ##  (1) keep Z sparse (no densification)
  ##  (2) avoid forming dense n×q intermediates Svar2/temp.C
  ##  (3) PRECOMPUTE row_qf = diag(Z var_eta0 Z') ONCE per update_beta()
  ##      (depends on var_eta0, not on b), so numDeriv Jacobian is much cheaper.
  ##  (4) ridge/line-search are "fallback only" to preserve equivalence as much as possible:
  ##      - first try plain solve(H, score)
  ##      - use ridge only if solve fails
  ##      - line-search only if full Newton step fails to reduce score norm
  ## -----------------------------
  update_beta <- function(beta, eta0, var_eta0, phase) {

    if (phase == 1L) {
      score_fun <- function(b) {
        comp <- fam_spec$E_R2_R3(b, eta0)
        E    <- comp$E
        as.numeric(crossprod(X, E))
      }
    } else {

      ## dgCMatrix slot access for O(nnz) accumulation
      Zp <- Z_M@p
      Zi <- Z_M@i
      Zx <- Z_M@x

      use_full_M <- (as.double(n) * as.double(q) <= as.double(ctrl$max_nq_mem))

      ## PRECOMPUTE row_qf[i] = z_i' var_eta0 z_i, once per update_beta() call
      row_qf <- numeric(n)

      if (use_full_M) {
        ## M = Z %*% var_eta0 (dense n×q), formed once per beta update
        M <- as.matrix(Z_M %*% var_eta0)

        ## row_qf[rows] += Z_rc * M_rc across nnz(Z)
        for (col in seq_len(q)) {
          lo <- Zp[col] + 1L
          hi <- Zp[col + 1L]
          if (hi >= lo) {
            idx  <- lo:hi
            rows <- Zi[idx] + 1L
            vals <- Zx[idx]
            row_qf[rows] <- row_qf[rows] + vals * M[rows, col]
          }
        }
      } else {
        ## memory-guarded: compute M[,col] on the fly per column; still only once per beta update
        for (col in seq_len(q)) {
          lo <- Zp[col] + 1L
          hi <- Zp[col + 1L]
          if (hi >= lo) {
            Mj <- as.numeric(Z_M %*% var_eta0[, col])  # length n
            idx  <- lo:hi
            rows <- Zi[idx] + 1L
            vals <- Zx[idx]
            row_qf[rows] <- row_qf[rows] + vals * Mj[rows]
          }
        }
      }

      score_fun <- function(b) {
        comp <- fam_spec$E_R2_R3(b, eta0)
        E    <- comp$E
        R2   <- comp$R2
        R3   <- comp$R3

        ## term1 = (z_i' var_eta0 z_i) * R3_i
        term1 <- row_qf * as.numeric(R3)

        ## bvec = (Z var_eta0)' term1 = var_eta0 %*% (Z' term1)
        zTa  <- as.numeric(Matrix::crossprod(Z_M, term1))  # q
        bvec <- as.numeric(var_eta0 %*% zTa)               # q

        ## term2 = R2 ∘ (Z bvec)
        term2 <- as.numeric(R2) * as.numeric(Z_M %*% bvec) # n

        trc_beta <- term1 + term2

        as.numeric(crossprod(X, E + 0.5 * trc_beta))
      }
    }

    b <- as.numeric(beta)

    for (it in seq_len(ctrl$beta_max_iter)) {
      score <- score_fun(b)
      if (!all(is.finite(score))) break

      smax <- max(abs(score))
      if (!is.finite(smax) || smax < ctrl$beta_tol) break

      ## Numerical Jacobian of score
      H <- tryCatch(numDeriv::jacobian(score_fun, b),
                    error = function(e) NULL)
      if (is.null(H) || any(!is.finite(H))) {
        if (isTRUE(ctrl$verbose)) message("beta Newton Jacobian failed; leaving beta unchanged for this EM iteration.")
        break
      }
      H <- 0.5 * (H + t(H))

      ## Try plain Newton solve first (most equivalent to previous behavior)
      step <- tryCatch(solve(H, score), error = function(e) NULL)

      ## If plain solve fails, use ridge escalation (fallback only)
      if (is.null(step) || any(!is.finite(step))) {
        ridge <- as.numeric(ctrl$beta_hess_ridge_init)
        step <- NULL

        repeat {
          Hr <- H + diag(ridge, length(b))
          step_try <- tryCatch(solve(Hr, score), error = function(e) NULL)

          if (!is.null(step_try) && all(is.finite(step_try))) {
            step <- step_try
            break
          }

          ridge <- ridge * 10
          if (!is.finite(ridge) || ridge > ctrl$beta_hess_ridge_max) break
        }

        if (is.null(step) || any(!is.finite(step))) {
          if (isTRUE(ctrl$verbose)) message("beta Newton solve failed (even with ridge); leaving beta unchanged for this EM iteration.")
          break
        }
      }

      ## Cap step size (inf norm) (safe numerical guard)
      step_inf <- max(abs(step))
      if (is.finite(step_inf) && step_inf > ctrl$beta_step_max) {
        step <- step * (ctrl$beta_step_max / step_inf)
      }

      ## Accept full step if it reduces score norm; otherwise do backtracking (fallback only)
      score_norm0 <- max(abs(score))

      b_full <- b - as.numeric(step)
      score_full <- score_fun(b_full)
      score_norm_full <- if (all(is.finite(score_full))) max(abs(score_full)) else Inf

      if (is.finite(score_norm_full) && score_norm_full <= (1 - 1e-4) * score_norm0) {
        b <- b_full
        next
      }

      ## Backtracking line search
      alpha <- 1.0
      accepted <- FALSE

      for (ls in seq_len(as.integer(ctrl$beta_ls_max_iter))) {
        b_new <- b - alpha * as.numeric(step)
        score_new <- score_fun(b_new)

        if (all(is.finite(score_new))) {
          score_norm1 <- max(abs(score_new))
          if (is.finite(score_norm1) && score_norm1 <= (1 - 1e-4 * alpha) * score_norm0) {
            b <- b_new
            accepted <- TRUE
            break
          }
        }

        alpha <- alpha / 2
      }

      if (!accepted) {
        ## If line search fails, take a very small step if it improves at all; else bail.
        b_try <- b - 1e-3 * as.numeric(step)
        score_try <- score_fun(b_try)
        if (all(is.finite(score_try)) && max(abs(score_try)) < score_norm0) {
          b <- b_try
        } else {
          if (isTRUE(ctrl$verbose)) message("beta line search failed; leaving beta unchanged for this EM iteration.")
          break
        }
      }
    }

    b
  }

  ## -----------------------------
  ## EM loop with staged approximation
  ## -----------------------------
  em_converged <- FALSE
  iter_used    <- 0L
  phase        <- 1L

  for (iter in seq_len(ctrl$em_max_iter)) {
    iter_used <- iter
    beta_old  <- beta
    tau2_old  <- tau2

    ## Laplace step: mode and covariance of eta given current beta,tau2
    lap <- update_eta_laplace(eta_init = eta, beta = beta, tau2 = tau2)
    eta0     <- lap$eta
    var_eta0 <- lap$var_eta

    ## FE corrections (if phase > 1)
    if (phase == 1L) {
      eta_hat     <- eta0
      var_eta_hat <- var_eta0
    } else {
      fe_inp <- fam_spec$FE_trace_inputs(beta, eta0)

      fe <- fe_trace_diagG(
        Z           = Z_M,
        var_eta     = var_eta0,
        temp_trc_C  = fe_inp$temp_trc_C,
        temp_trc_D  = fe_inp$temp_trc_D,
        max_nq_mem  = ctrl$max_nq_mem
      )

      trc_y1 <- fe$trc_y1
      trc_y2 <- fe$trc_y2

      ## Mean correction: eta_hat = eta_mode + 0.5 * trc_y1
      eta_hat <- eta0 + 0.5 * trc_y1

      if (approx_lab == "FE_mean") {
        var_eta_hat <- var_eta0
      } else {
        ## FE_full: stage covariance correction only in phase 3 (after FE-mean stabilizes)
        if (phase < 3L) {
          var_eta_hat <- var_eta0
        } else {
          ## trc_y2 is a Matrix::Diagonal (S4). Apply diagonal-only correction safely.
          var_eta_hat <- as.matrix(var_eta0)
          diag(var_eta_hat) <- diag(var_eta_hat) + 0.5 * Matrix::diag(trc_y2)
        }
      }

      ## symmetry guard (keeps type as base matrix)
      var_eta_hat <- 0.5 * (var_eta_hat + t(var_eta_hat))
    }

    ## Update tau2 using FE-corrected moments (or Laplace moments in phase 1):
    tau2 <- mean(diag(var_eta_hat) + eta_hat^2)
    tau2 <- max(tau2, ctrl$vc_eps)
    G    <- Matrix::Diagonal(q, x = rep.int(tau2, q))

    ## Update beta (robust Newton)
    beta <- update_beta(beta = beta, eta0 = eta0, var_eta0 = var_eta0, phase = phase)

    ## carry forward FE-corrected (or Laplace) eta moments
    eta     <- eta_hat
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
    phase        = phase,
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

  fit
}
