# ======================================================================
# tests/testthat/test-glmmFEL-matrix-only.R
#
# Tests for the *simplified* glmmFEL development setting:
#   - Matrix-only interface: glmmFEL(y, X, Z, ...)
#   - Single variance component: G = tau2 * I_q (one tau2 scalar)
#   - Keep FE mean/full machinery + RSPL/MSPL (PQL) machinery
#   - Convert dense X/Z to sparse immediately
#
# Coverage:
#   (1) Derivative equivalence vs mvglmmrank rder*/jder* (moderate q-range)
#   (2) FE trace corrections for diag-G match finite differences (tilting)
#   (3) Smoke fits for Laplace / FE_mean / FE_full and RSPL/MSPL
#   (4) Multi-membership Z smoke (2 nonzeros per row)
#   (5) Dense->sparse coercion
#   (6) Smoke fits for binomial_logit across Laplace/FE/PL
# ======================================================================

testthat::local_edition(3)
testthat::skip_on_cran() # Extended simulation smoke tests; fast independent checks run on CRAN.

# ---- helpers ----

make_Z_oneway <- function(id, n_id = max(as.integer(id))) {
  n <- length(id)
  Matrix::sparseMatrix(
    i = seq_len(n),
    j = as.integer(id),
    x = 1,
    dims = c(n, n_id)
  )
}

simulate_oneway <- function(family = c("binomial_probit", "binomial_logit", "poisson_log"),
                            n_id = 30L,
                            m_per_id = 6L,
                            beta = c(0.2, 0.7),
                            tau2 = 0.5,
                            seed = NULL) {
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)

  id <- factor(rep(seq_len(n_id), each = m_per_id))
  n  <- length(id)

  x <- stats::rnorm(n)
  X <- stats::model.matrix(~ x)

  Z <- make_Z_oneway(id, n_id)

  eta <- stats::rnorm(n_id, sd = sqrt(tau2))
  lp  <- as.vector(X %*% beta + Z %*% eta)

  y <- switch(
    family,
    binomial_probit = stats::rbinom(n, 1, stats::pnorm(lp)),
    binomial_logit  = stats::rbinom(n, 1, stats::plogis(lp)),
    poisson_log     = stats::rpois(n, exp(lp))
  )

  list(y = y, X = X, Z = Z, id = id, x = x, truth = list(beta = beta, tau2 = tau2))
}

simulate_multimember_2perrow <- function(family = c("binomial_probit", "binomial_logit", "poisson_log"),
                                         n_re = 12L,
                                         n_obs = 60L,
                                         beta = 0.0,
                                         tau2 = 0.7,
                                         seed = NULL) {
  # Multi-membership: each row has +1 for "home", -1 for "away"
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)

  teams <- seq_len(n_re)
  home <- sample(teams, n_obs, replace = TRUE)
  away <- sample(teams, n_obs, replace = TRUE)
  # enforce home != away
  same <- home == away
  while (any(same)) {
    away[same] <- sample(teams, sum(same), replace = TRUE)
    same <- home == away
  }

  i <- rep(seq_len(n_obs), times = 2)
  j <- c(home, away)
  x <- c(rep(1, n_obs), rep(-1, n_obs))

  Z <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_obs, n_re))
  X <- matrix(1, n_obs, 1)
  colnames(X) <- "(Intercept)"

  eta <- stats::rnorm(n_re, sd = sqrt(tau2))
  lp  <- as.vector(X %*% beta + Z %*% eta)

  y <- switch(
    family,
    binomial_probit = stats::rbinom(n_obs, 1, stats::pnorm(lp)),
    binomial_logit  = stats::rbinom(n_obs, 1, stats::plogis(lp)),
    poisson_log     = stats::rpois(n_obs, exp(lp))
  )

  list(y = y, X = X, Z = Z, truth = list(beta = beta, tau2 = tau2))
}



# ---- (2) FE trace diag-G finite-difference validation ----

testthat::test_that("fe_trace_diagG matches finite-difference tilting derivatives (sign check)", {
  testthat::skip_on_cran()  # safe: this is a numeric finite-difference validation

  set.seed(123)
  n <- 120
  q <- 25
  tau2 <- 0.8

  # Build an arbitrary sparse-ish Z
  Z <- Matrix::rsparsematrix(n, q, density = 0.08)
  Z <- Matrix::drop0(Z)

  # One fixed effect (intercept) to keep it simple here
  X <- matrix(1, n, 1)
  beta <- 0.0

  # Simulate eta, then y from probit
  eta_true <- stats::rnorm(q, sd = sqrt(tau2))
  lin_true <- as.vector(X %*% beta + Z %*% eta_true)
  y <- stats::rbinom(n, 1, stats::pnorm(lin_true))

  # === Helpers needed for finite-difference tilting ===
  derivs <- glmmFEL:::fe_derivatives_binomial_probit()
  d1 <- derivs$d1; d2 <- derivs$d2; d3 <- derivs$d3; d4 <- derivs$d4

  Ginv <- Matrix::Diagonal(q, x = rep(1 / tau2, q))

  # Negative log posterior up to constants (tilted by c):
  # f_c(eta) = -sum log Phi(s * (Xb + Z eta)) + 0.5 eta'Ginv eta - c'eta
  # where s = +1 for y=1, s=-1 for y=0
  f_c <- function(eta, c) {
    s <- ifelse(y == 1, 1, -1)
    lin <- as.vector(X %*% beta + Z %*% eta)
    qv  <- s * lin
    nll <- -sum(stats::pnorm(qv, log.p = TRUE))
    prior <- 0.5 * as.numeric(crossprod(eta, as.matrix(Ginv %*% eta)))
    tilt  <- -as.numeric(crossprod(c, eta))
    nll + prior + tilt
  }

  grad_c <- function(eta, c) {
    s <- ifelse(y == 1, 1, -1)
    lin <- as.vector(X %*% beta + Z %*% eta)
    qv  <- s * lin
    # d/deta [-log Phi(qv)] = -Z' [ s * d1(qv) ]
    g_like <- -Matrix::crossprod(Z, s * d1(qv))
    g_prior <- as.matrix(Ginv %*% eta)
    g_tilt  <- -c
    as.vector(g_like + g_prior + g_tilt)
  }

  hess_c <- function(eta, c) {
    s <- ifelse(y == 1, 1, -1)
    lin <- as.vector(X %*% beta + Z %*% eta)
    qv  <- s * lin
    # Second derivative of -log Phi(qv) w.r.t. linear predictor is (-d2(qv)) >= 0
    W <- -d2(qv)
    # Hessian: Z' diag(W) Z + Ginv
    Matrix::crossprod(Z, W * Z) + Ginv
  }

  solve_mode <- function(c, eta_init = rep(0, q), maxit = 60, tol = 1e-10) {
    eta <- eta_init
    for (it in seq_len(maxit)) {
      g <- grad_c(eta, c)
      H <- hess_c(eta, c)
      step <- as.vector(solve(as.matrix(H), g))
      eta_new <- eta - step
      if (max(abs(step)) < tol) break
      eta <- eta_new
    }
    eta
  }

  # c=0 Laplace mode and covariance
  c0 <- rep(0, q)
  eta0 <- solve_mode(c0)
  H0 <- hess_c(eta0, c0)
  var_eta <- solve(as.matrix(H0))

  # Build FE ingredients at c=0
  s <- ifelse(y == 1, 1, -1)
  lin0 <- as.vector(X %*% beta + Z %*% eta0)
  q0 <- s * lin0
  temp_trc_C <- s * d3(q0)          # matches Appendix-B form via s * d3(s*lin)
  temp_trc_D <- d4(q0)              # d4 is already w.r.t q0

  fe <- glmmFEL:::fe_trace_diagG(
    Z = Z,
    var_eta = var_eta,
    temp_trc_C = temp_trc_C,
    temp_trc_D = temp_trc_D
  )
  trc_y1 <- fe$trc_y1
  trc_y2 <- fe$trc_y2

  testthat::expect_equal(length(trc_y1), q)
  testthat::expect_equal(dim(trc_y2), c(q, q))
  testthat::expect_true(all(abs(trc_y2[row(trc_y2) != col(trc_y2)]) < 1e-12))

  # --- Finite differences for one component j ---
  j <- 7L
  eps <- 1e-4
  ej <- rep(0, q); ej[j] <- 1

  # dsig/dc_j and d2sig/dc_j^2 at c=0 by finite differences:
  # Sigma(c) = Hessian at tilted mode eta_hat(c)
  Sigma_at <- function(c) {
    eta_c <- solve_mode(c, eta_init = eta0)
    hess_c(eta_c, c)
  }

  Sigma_p <- Sigma_at(+eps * ej)
  Sigma_m <- Sigma_at(-eps * ej)
  Sigma_0 <- Sigma_at(c0)

  dSig_dc_fd  <- (Sigma_p - Sigma_m) / (2 * eps)
  d2Sig_dc2_fd <- (Sigma_p - 2 * Sigma_0 + Sigma_m) / (eps^2)

  # Our convention inside FE code is D = -dSigma/dc and D2 = -d2Sigma/dc2
  dsig_fd  <- -as.matrix(dSig_dc_fd)
  d2sig_fd <- -as.matrix(d2Sig_dc2_fd)

  trc_y1_fd <- sum(diag(var_eta %*% dsig_fd))
  term1_fd  <- sum(diag(var_eta %*% d2sig_fd))
  A_fd <- var_eta %*% dsig_fd
  term2_fd  <- sum(A_fd * t(A_fd))  # trace(A^2) when var_eta is symmetric

  trc_y2_fd <- term1_fd + term2_fd

  # Compare signs + magnitude (relative tolerance)
  testthat::expect_equal(trc_y1[j], trc_y1_fd, tolerance = 5e-3)
  testthat::expect_equal(trc_y2[j, j], trc_y2_fd, tolerance = 1e-2)
})

# ---- (3) Smoke fits for core methods ----

testthat::test_that("glmmFEL runs for Laplace / FE_mean / FE_full (one-way intercept, probit)", {
  dat <- simulate_oneway(
    family = "binomial_probit",
    n_id = 25L, m_per_id = 8L,
    beta = c(-0.2, 0.6),
    tau2 = 0.6,
    seed = 101
  )

  fits <- list(
    Laplace = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_probit",
                               approx = "Laplace", max_iter = 120, tol = 1e-5,
                               control = list(verbose = FALSE)),
    FE_mean = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_probit",
                               approx = "FE_mean", max_iter = 160, tol = 1e-5,
                               control = list(verbose = FALSE)),
    FE_full = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_probit",
                               approx = "FE_full", max_iter = 200, tol = 1e-5,
                               control = list(verbose = FALSE))
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    testthat::expect_true(inherits(fit, "glmmFELMod"))
    testthat::expect_equal(length(fit$beta), ncol(dat$X))
    testthat::expect_equal(length(fit$eta), ncol(dat$Z))
    testthat::expect_true(is.finite(fit$logLik))
    testthat::expect_true(is.finite(fit$tau2))
    testthat::expect_true(fit$tau2 >= 0)
  }
})

testthat::test_that("glmmFEL runs for RSPL and MSPL (one-way intercept, poisson)", {
  dat <- simulate_oneway(
    family = "poisson_log",
    n_id = 20L, m_per_id = 10L,
    beta = c(log(2), 0.2),
    tau2 = 0.4,
    seed = 202
  )

  fits <- list(
    RSPL = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "poisson_log",
                            approx = "RSPL", max_iter = 80, tol = 1e-5,
                            control = list(verbose = FALSE)),
    MSPL = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "poisson_log",
                            approx = "MSPL", max_iter = 80, tol = 1e-5,
                            control = list(verbose = FALSE))
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    testthat::expect_true(inherits(fit, "glmmFELMod"))
    testthat::expect_equal(length(fit$beta), ncol(dat$X))
    testthat::expect_equal(length(fit$eta), ncol(dat$Z))
    testthat::expect_true(is.finite(fit$tau2))
    testthat::expect_true(fit$tau2 >= 0)
  }
})

# ---- (4) Multi-membership smoke ----

testthat::test_that("glmmFEL runs with multi-membership Z (2 nonzeros per row)", {
  dat <- simulate_multimember_2perrow(
    family = "binomial_probit",
    n_re = 12L, n_obs = 50L,
    beta = 0.0, tau2 = 0.8,
    seed = 303
  )

  fit <- glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_probit",
                          approx = "Laplace", max_iter = 150, tol = 1e-5,
                          control = list(verbose = FALSE))

  testthat::expect_true(inherits(fit, "glmmFELMod"))
  testthat::expect_equal(length(fit$eta), ncol(dat$Z))
  testthat::expect_true(is.finite(fit$tau2))
})

# ---- (5) Dense->sparse coercion ----

testthat::test_that("dense Z is converted to sparse internally", {
  dat <- simulate_oneway(
    family = "binomial_probit",
    n_id = 10L, m_per_id = 6L,
    beta = c(0.0, 0.3),
    tau2 = 0.5,
    seed = 404
  )

  # Make a dense Z on purpose
  Z_dense <- as.matrix(dat$Z)

  fit <- glmmFEL::glmmFEL(dat$y, dat$X, Z_dense, family = "binomial_probit",
                          approx = "Laplace", max_iter = 100, tol = 1e-5,
                          control = list(verbose = FALSE))

  testthat::expect_true(inherits(fit$Z, "dgCMatrix"))
  testthat::expect_equal(dim(fit$Z), dim(dat$Z))
})

# ---- (6) binomial_logit smoke fits ----

testthat::test_that("glmmFEL runs for Laplace / FE_mean / FE_full (one-way intercept, logit)", {
  dat <- simulate_oneway(
    family = "binomial_logit",
    n_id = 25L, m_per_id = 8L,
    beta = c(-0.2, 0.6),
    tau2 = 0.6,
    seed = 505
  )

  fits <- list(
    Laplace = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                               approx = "Laplace", max_iter = 120, tol = 1e-5,
                               control = list(verbose = FALSE)),
    FE_mean = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                               approx = "FE_mean", max_iter = 160, tol = 1e-5,
                               control = list(verbose = FALSE)),
    FE_full = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                               approx = "FE_full", max_iter = 200, tol = 1e-5,
                               control = list(verbose = FALSE))
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    testthat::expect_true(inherits(fit, "glmmFELMod"))
    testthat::expect_equal(length(fit$beta), ncol(dat$X))
    testthat::expect_equal(length(fit$eta), ncol(dat$Z))
    testthat::expect_true(is.finite(fit$logLik))
    testthat::expect_true(is.finite(fit$tau2))
    testthat::expect_true(fit$tau2 >= 0)
  }
})

testthat::test_that("glmmFEL runs for RSPL and MSPL (one-way intercept, logit)", {
  dat <- simulate_oneway(
    family = "binomial_logit",
    n_id = 20L, m_per_id = 10L,
    beta = c(-0.1, 0.2),
    tau2 = 0.4,
    seed = 606
  )

  fits <- list(
    RSPL = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                            approx = "RSPL", max_iter = 80, tol = 1e-5,
                            control = list(verbose = FALSE)),
    MSPL = glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                            approx = "MSPL", max_iter = 80, tol = 1e-5,
                            control = list(verbose = FALSE))
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    testthat::expect_true(inherits(fit, "glmmFELMod"))
    testthat::expect_equal(length(fit$beta), ncol(dat$X))
    testthat::expect_equal(length(fit$eta), ncol(dat$Z))
    testthat::expect_true(is.finite(fit$tau2))
    testthat::expect_true(fit$tau2 >= 0)
  }
})

testthat::test_that("glmmFEL runs with multi-membership Z (2 nonzeros per row, logit)", {
  dat <- simulate_multimember_2perrow(
    family = "binomial_logit",
    n_re = 12L, n_obs = 50L,
    beta = 0.0, tau2 = 0.8,
    seed = 707
  )

  fit <- glmmFEL::glmmFEL(dat$y, dat$X, dat$Z, family = "binomial_logit",
                          approx = "Laplace", max_iter = 150, tol = 1e-5,
                          control = list(verbose = FALSE))

  testthat::expect_true(inherits(fit, "glmmFELMod"))
  testthat::expect_equal(length(fit$eta), ncol(dat$Z))
  testthat::expect_true(is.finite(fit$tau2))
})


test_that("glmmFEL input checks fire", {
  dat <- simulate_oneway(seed = 1)

  testthat::expect_error(
    glmmFEL(dat$y[-1], dat$X, dat$Z),
    "Z has .* rows but expected n"
  )
  expect_error(glmmFEL(dat$y, dat$X[,1,drop=FALSE], matrix(0, nrow(dat$X), 0)), "at least one")
  expect_error(glmmFEL(c(0,2,1,0), dat$X[1:4,], dat$Z[1:4,], family="binomial_probit"), "0/1")
  expect_error(glmmFEL(c(0,-1,2,3), dat$X[1:4,], dat$Z[1:4,], family="poisson_log"), "nonnegative")
})


test_that("permuting Z columns leaves beta/tau2 unchanged", {
  dat <- simulate_oneway(family="binomial_logit", seed=2)
  perm <- sample(seq_len(ncol(dat$Z)))
  Zp <- dat$Z[, perm]

  expect_warning(fit1 <- glmmFEL(dat$y, dat$X, dat$Z, family="binomial_logit", approx="Laplace",
                  max_iter=3, tol=1e-5), "did not converge")
  expect_warning(fit2 <- glmmFEL(dat$y, dat$X, Zp, family="binomial_logit", approx="Laplace",
                  max_iter=3, tol=1e-5), "did not converge")

  expect_equal(fit1$beta, fit2$beta, tolerance=1e-6)
  expect_equal(fit1$tau2, fit2$tau2, tolerance=1e-6)
})


test_that("Laplace fit roughly recovers parameters (sanity)", {
  dat <- simulate_oneway(family="binomial_probit", n_id=40, m_per_id=8,
                         beta=c(-0.3, 0.8), tau2=0.7, seed=3)
  fit <- glmmFEL(dat$y, dat$X, dat$Z, family="binomial_probit", approx="Laplace",
                 max_iter=150, tol=1e-6)

  expect_true(sign(fit$beta[2]) == sign(dat$truth$beta[2]))
  expect_true(fit$tau2 > 0.05 && fit$tau2 < 5)   # wide bounds
})


test_that("tau2_init extremes remain stable", {
  dat <- simulate_oneway(family="poisson_log", seed=4)

  fit_small <- glmmFEL(dat$y, dat$X, dat$Z, family="poisson_log", approx="Laplace",
                       control=list(tau2_init=1e-10), max_iter=80, tol=1e-5)
  fit_big <- glmmFEL(dat$y, dat$X, dat$Z, family="poisson_log", approx="Laplace",
                     control=list(tau2_init=1e3), max_iter=80, tol=1e-5)

  expect_true(is.finite(fit_small$tau2) && fit_small$tau2 >= 0)
  expect_true(is.finite(fit_big$tau2) && fit_big$tau2 >= 0)
})


test_that("FE_full var_eta is base matrix and symmetric", {
  dat <- simulate_oneway(family="poisson_log", n_id=6, m_per_id=20, seed=5)
  fit <- glmmFEL(dat$y, dat$X, dat$Z, family="poisson_log", approx="FE_full",
                 max_iter=300, tol=1e-5)

  expect_true(is.matrix(fit$var_eta))
  expect_equal(fit$var_eta, t(fit$var_eta), tolerance=1e-10)
  expect_true(all(is.finite(diag(fit$var_eta))))
})
