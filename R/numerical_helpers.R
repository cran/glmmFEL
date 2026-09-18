# Internal numerical and validation helpers. Statistical identities and numerical
# checks are documented in project_review/NUMERICAL_AUDIT.md (development only).

glmmfe_validate_data <- function(y, X, Z, family) {
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y) || any(!is.finite(y)))
    stop("y must be a nonempty finite numeric response vector.")
  if (nrow(X) != length(y) || nrow(Z) != length(y))
    stop("y, X, and Z must have the same number of rows/observations.")
  if (!ncol(X) || any(!is.finite(X))) stop("X must have at least one column and finite entries.")
  if (!ncol(Z) || any(!is.finite(Z@x))) stop("Z must have at least one column and finite entries.")
  if (!any(Z@x != 0)) stop("Z must contain nonzero entries to identify random effects.")
  if (qr(X)$rank < ncol(X)) stop("Fixed-effects design matrix X is not full rank.")
  if (family == "poisson_log") {
    if (any(y < 0 | y != floor(y))) stop("For poisson_log, y must be nonnegative integer counts.")
  } else if (any(!y %in% c(0, 1))) stop("For binomial families, y must be 0/1 (Bernoulli responses).")
  no_finite_intercept <- (family == "poisson_log" && all(y == 0)) ||
    (family != "poisson_log" && length(unique(y)) == 1L)
  if (no_finite_intercept && qr(cbind(X, 1))$rank == ncol(X))
    stop("Constant response has no finite intercept estimate for this model.")
  invisible(TRUE)
}

glmmfe_check_separation <- function(y, X, beta, family) {
  if (family == "poisson_log") return(invisible(NULL))
  margin <- (2 * y - 1) * as.numeric(X %*% beta)
  scale <- max(abs(margin))
  # A strictly positive signed margin is a certificate of complete separation:
  # multiplying beta by a growing constant increases every success probability.
  # Require a relative margin to guard against floating-point roundoff. This is
  # deliberately a sufficient check, not a complete quasi-separation detector.
  if (is.finite(scale) && scale > 0 && all(margin > 1e-8 * scale))
    stop("Complete separation in X: binary data have no finite fixed-effect estimate.")
  invisible(NULL)
}

glmmfe_validate_control <- function(control, defaults) {
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("control must be a list with unique, nonempty names.")
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control setting: ", paste(unknown, collapse = ", "))
  for (name in names(control)) {
    if (is.null(control[[name]]) && !is.null(defaults[[name]]))
      stop(name, " cannot be NULL.")
  }
  result <- utils::modifyList(defaults, control)
  for (name in names(result)) {
    x <- result[[name]]
    if (is.null(x)) next
    if (name %in% c("verbose", "trace")) {
      if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(name, " must be TRUE or FALSE.")
    } else {
      if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
        stop(name, " must be a positive finite numeric scalar.")
      if (grepl("_iter$", name) && (x != floor(x) || x > .Machine$integer.max))
        stop(name, " must be a positive integer.")
    }
  }
  result
}

glmmfe_log_response <- function(y, linear, family) {
  if (family == "poisson_log") return(sum(stats::dpois(y, lambda = exp(linear), log = TRUE)))
  signed <- (2 * y - 1) * linear
  if (family == "binomial_probit") sum(stats::pnorm(signed, log.p = TRUE)) else
    sum(stats::plogis(signed, log.p = TRUE))
}

glmmfe_row_variance <- function(Z, V, max_nq_mem) {
  if (all(V[row(V) != col(V)] == 0))
    return(as.numeric((Z^2) %*% diag(V)))
  if (as.double(nrow(Z)) * ncol(Z) <= max_nq_mem)
    return(as.numeric(Matrix::rowSums(Z * as.matrix(Z %*% V))))
  out <- numeric(nrow(Z))
  for (j in seq_len(ncol(Z))) {
    out <- out + as.numeric(Z[, j]) * as.numeric(Z %*% V[, j])
  }
  out
}

glmmfe_mode <- function(start, beta, tau2, y, X, Z, family, bundle, control) {
  eta <- as.numeric(start)
  q <- ncol(Z)
  disjoint <- all(Matrix::rowSums(Z != 0) <= 1)
  Z_squared <- if (disjoint) Z^2 else NULL
  objective <- function(b) -glmmfe_log_response(y, as.numeric(X %*% beta + Z %*% b), family) +
    sum(b^2) / (2 * tau2)
  evaluate <- function(b) {
    d <- bundle$E_R2_R3(beta, b)
    gradient <- as.numeric(-Matrix::crossprod(Z, d$E)) + b / tau2
    hessian <- if (disjoint) as.numeric(Matrix::crossprod(Z_squared, -d$R2)) + 1 / tau2 else
      as.matrix(Matrix::crossprod(Z, Z * (-d$R2))) + diag(1 / tau2, q)
    list(gradient = gradient, hessian = hessian)
  }
  converged <- FALSE
  reason <- "iteration_limit"
  for (iteration in seq_len(control$eta_max_iter)) {
    state <- evaluate(eta)
    if (any(!is.finite(state$gradient)) || any(!is.finite(state$hessian))) {
      reason <- "nonfinite_mode_derivatives"; break
    }
    if (max(abs(state$gradient)) <= control$eta_tol_grad) {
      converged <- TRUE; reason <- "converged"; break
    }
    if (disjoint) {
      if (any(state$hessian <= 0)) { reason <- "mode_hessian_not_positive_definite"; break }
      step <- state$gradient / state$hessian
    } else {
      ch <- tryCatch(chol(state$hessian), error = function(e) NULL)
      if (is.null(ch)) { reason <- "mode_hessian_not_positive_definite"; break }
      step <- backsolve(ch, forwardsolve(t(ch), state$gradient))
    }
    value <- objective(eta)
    accepted <- FALSE
    for (halving in 0:30) {
      scale <- 2^(-halving)
      candidate <- eta - scale * step
      next_value <- objective(candidate)
      allowance <- 10 * .Machine$double.eps * max(1, abs(value))
      if (is.finite(next_value) && next_value <= value -
          1e-4 * scale * sum(state$gradient * step) + allowance) {
        eta <- candidate; accepted <- TRUE; break
      }
    }
    if (!accepted) { reason <- "mode_line_search_failed"; break }
  }
  # Always recompute at the returned mode, including iteration-limit exits.
  state <- evaluate(eta)
  covariance <- if (disjoint) {
    if (all(is.finite(state$hessian)) && all(state$hessian > 0)) diag(1 / state$hessian, q) else NULL
  } else tryCatch(chol2inv(chol(state$hessian)), error = function(e) NULL)
  gradient_norm <- max(abs(state$gradient))
  converged <- !is.null(covariance) && is.finite(gradient_norm) && gradient_norm <= control$eta_tol_grad
  if (converged) reason <- "converged"
  list(eta = eta, var_eta = if (is.null(covariance)) matrix(NA_real_, q, q) else covariance,
       converged = converged, gradient = gradient_norm, reason = reason, iterations = iteration)
}

glmmfe_beta_update <- function(beta, mode, covariance, phase, X, Z, bundle, control) {
  n <- nrow(X)
  rowvar <- shift <- numeric(n)
  if (phase > 1L) {
    rowvar <- glmmfe_row_variance(Z, covariance, control$max_nq_mem)
    # The posterior in the E step is fixed at the OLD parameters. In particular,
    # its skew correction must not be recomputed using each candidate beta.
    old_third <- bundle$E_R2_R3(beta, mode)$R3
    shift <- 0.5 * as.numeric(Z %*% (covariance %*% Matrix::crossprod(Z, rowvar * old_third)))
  }
  evaluate <- function(b) {
    d <- bundle$E_R2_R3(b, mode)
    score <- as.numeric(crossprod(X, d$E + d$R2 * shift + 0.5 * d$R3 * rowvar))
    fourth <- if (phase > 1L) bundle$FE_trace_inputs(b, mode)$temp_trc_D else numeric(n)
    jacobian <- crossprod(X, X * (d$R2 + d$R3 * shift + 0.5 * fourth * rowvar))
    list(score = score, jacobian = jacobian)
  }
  b <- as.numeric(beta)
  reason <- "iteration_limit"
  for (iteration in seq_len(control$beta_max_iter)) {
    state <- evaluate(b)
    if (any(!is.finite(state$score)) || any(!is.finite(state$jacobian))) {
      reason <- "nonfinite_beta_derivatives"; break
    }
    norm <- max(abs(state$score))
    if (norm <= control$beta_tol) { reason <- "converged"; break }
    step <- tryCatch(solve(state$jacobian, state$score), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      ridge <- control$beta_hess_ridge_init
      while (ridge <= control$beta_hess_ridge_max) {
        step <- tryCatch(solve(state$jacobian - diag(ridge, length(b)), state$score),
                         error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) break
        ridge <- ridge * 10
      }
    }
    if (is.null(step) || any(!is.finite(step))) { reason <- "beta_solve_failed"; break }
    step <- step / max(1, max(abs(step)) / control$beta_step_max)
    accepted <- FALSE
    for (halving in 0:control$beta_ls_max_iter) {
      candidate <- b - 2^(-halving) * step
      new_score <- evaluate(candidate)$score
      if (all(is.finite(new_score)) && max(abs(new_score)) < norm) {
        b <- candidate; accepted <- TRUE; break
      }
    }
    if (!accepted) { reason <- "beta_line_search_failed"; break }
  }
  score <- evaluate(b)$score
  norm <- max(abs(score))
  converged <- is.finite(norm) && norm <= control$beta_tol
  list(beta = b, converged = converged, score = norm, reason = if (converged) "converged" else reason)
}
