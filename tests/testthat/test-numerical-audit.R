test_that("probit derivatives agree with independent high-precision values in the tails", {
  # Generated with mpmath 1.3.0 at 70 decimal digits by differentiating log(erfc).
  reference <- read.csv(test_path("fixtures", "probit-reference.csv"))
  d <- glmmFEL:::fe_derivatives_binomial_probit()
  for (k in 0:4) {
    name <- paste0("d", k)
    expect_lt(max(abs(d[[name]](reference$q) - reference[[name]]) / pmax(1, abs(reference[[name]]))), 1e-9)
  }
})

test_that("logit and Poisson derivatives agree with numerical differentiation", {
  skip_if_not_installed("numDeriv")
  logit <- glmmFEL:::fe_derivatives_binomial_logit()
  poisson <- glmmFEL:::fe_derivatives_poisson_log()
  points <- c(-4, -0.3, 0, 1.4, 4)
  for (k in 1:4) {
    expect_equal(logit[[paste0("d", k)]](points),
                 vapply(points, function(x) numDeriv::grad(logit[[paste0("d", k - 1)]], x), 0),
                 tolerance = 1e-7)
    previous <- if (k <= 2) function(x) poisson[[paste0("d", k - 1)]](x, y = 3) else
      poisson[[paste0("d", k - 1)]]
    actual <- if (k == 1) poisson$d1(points, 3) else poisson[[paste0("d", k)]](points)
    expect_equal(actual, vapply(points, function(x) numDeriv::grad(previous, x), 0), tolerance = 1e-7)
  }
})

test_that("FE corrections agree with independently tilted Laplace integrals", {
  # A dense, nonnested, two-dimensional design exercises off-diagonal Hessians.
  Z <- Matrix::Matrix(matrix(c(1, 0.2, 1, -0.5, 0.3, 1, -0.2, 1, 1, 0.4, 0.6, 1),
                              ncol = 2, byrow = TRUE), sparse = TRUE)
  X <- matrix(1, 6, 1); beta <- 0.4; tau2 <- 0.7; y <- c(2, 1, 3, 2, 4, 1)
  tilted <- function(c) {
    objective <- function(b) -sum(dpois(y, exp(beta + as.numeric(Z %*% b)), log = TRUE)) +
      sum(b^2) / (2 * tau2) - sum(c * b)
    gradient <- function(b) as.numeric(Matrix::crossprod(Z, exp(beta + as.numeric(Z %*% b)) - y)) +
      b / tau2 - c
    opt <- optim(c(0, 0), objective, gradient, method = "BFGS",
                 control = list(reltol = 1e-13, maxit = 500))
    b <- opt$par
    # Refine to machine accuracy without relying on glmmFEL's mode routine.
    for (i in 1:6) {
      mu <- exp(beta + as.numeric(Z %*% b))
      H <- as.matrix(Matrix::crossprod(Z, Z * mu)) + diag(1 / tau2, 2)
      b <- b - solve(H, gradient(b))
    }
    mu <- exp(beta + as.numeric(Z %*% b))
    H <- as.matrix(Matrix::crossprod(Z, Z * mu)) + diag(1 / tau2, 2)
    list(value = -objective(b) - as.numeric(determinant(H, logarithm = TRUE)$modulus) / 2,
         mode = b, variance = solve(H), mu = mu)
  }
  base <- tilted(c(0, 0))
  fe <- glmmFEL:::fe_trace_diagG(Z, base$variance, -base$mu, -base$mu)
  small_memory <- glmmFEL:::fe_trace_diagG(Z, base$variance, -base$mu, -base$mu, max_nq_mem = 1)
  expect_equal(fe, small_memory, tolerance = 1e-12)
  h <- 0.001
  for (j in 1:2) {
    c <- c(0, 0); c[j] <- h
    plus <- tilted(c)$value; minus <- tilted(-c)$value
    expect_equal(base$mode[j] + fe$trc_y1[j] / 2, (plus - minus) / (2 * h), tolerance = 2e-6)
    expect_equal(base$variance[j, j] + fe$trc_y2[j, j] / 2,
                 (plus - 2 * base$value + minus) / h^2, tolerance = 2e-6)
  }
})

test_that("disjoint-cluster FE fast path agrees with scalar quadrature", {
  y <- c(3, 4, 7, 2, 5, 6, 3, 5, 4, 6)
  xbeta <- log(3); tau2 <- 0.5
  logposterior <- function(b) sum(dpois(y, exp(xbeta + b), log = TRUE)) + dnorm(b, sd = sqrt(tau2), log = TRUE)
  mode <- optimize(function(b) -logposterior(b), c(-3, 3), tol = 1e-12)$minimum
  density <- function(b) vapply(b, function(x) exp(logposterior(x) - logposterior(mode)), 0)
  norm <- integrate(density, -5, 5, rel.tol = 1e-11)$value
  mean <- integrate(function(b) b * density(b), -5, 5, rel.tol = 1e-11)$value / norm
  variance <- integrate(function(b) (b - mean)^2 * density(b), -5, 5, rel.tol = 1e-11)$value / norm
  mu <- rep(exp(xbeta + mode), length(y))
  V <- matrix(1 / (sum(mu) + 1 / tau2))
  Z <- Matrix::Matrix(matrix(1, length(y), 1), sparse = TRUE)
  fe <- glmmFEL:::fe_trace_diagG(Z, V, -mu, -mu)
  corrected_mean <- mode + fe$trc_y1 / 2
  corrected_variance <- V[1] + fe$trc_y2[1, 1] / 2
  expect_lt(abs(corrected_mean - mean), abs(mode - mean))
  expect_lt(abs(corrected_variance - variance), abs(V[1] - variance))
  expect_lt(abs(corrected_mean - mean), 0.0002)
  expect_lt(abs(corrected_variance - variance), 0.0002)
})

audit_data <- function() {
  x <- rep(c(-1, 0, 1), 4)
  list(y = c(1, 2, 4, 2, 5, 6, 1, 1, 2, 4, 5, 9), X = cbind(1, x),
       Z = Matrix::sparseMatrix(i = 1:12, j = rep(1:4, each = 3), x = 1))
}

test_that("both binary links fit and predict with FEL and PL", {
  id <- rep(1:4, each = 16)
  x <- rep(seq(-1, 1, length.out = 16), 4)
  X <- cbind(1, x); Z <- Matrix::sparseMatrix(i = seq_along(id), j = id, x = 1)
  linear <- -0.3 + 0.5 * x + c(-1.5, -0.5, 0.5, 1.5)[id]
  for (link in c("logit", "probit")) {
    set.seed(214)
    y <- rbinom(length(id), 1, if (link == "logit") plogis(linear) else pnorm(linear))
    for (method in c("FE_full", "MSPL")) {
      fit <- glmmFEL(y, X, Z, family = binomial(link), approx = method, max_iter = 400)
      expect_true(isTRUE(fit$convergence$em_converged))
      if (method == "MSPL") expect_true(fit$convergence$pql_converged)
      expect_true(all(fitted(fit) > 0 & fitted(fit) < 1))
    }
  }
})

test_that("invalid inputs and numerical controls fail clearly", {
  d <- audit_data()
  fit <- function(y = d$y, X = d$X, Z = d$Z, ...) glmmFEL(y, X, Z, family = "poisson_log", ...)
  expect_error(fit(y = replace(d$y, 1, NA)), "finite")
  expect_error(fit(y = d$y + 0.5), "integer")
  expect_error(fit(y = rep(0, 12)), "finite intercept")
  expect_error(fit(X = replace(d$X, 1, Inf)), "finite")
  expect_error(fit(Z = matrix(NA_real_, 12, 1)), "finite")
  expect_error(fit(Z = matrix(0, 12, 1)), "nonzero")
  expect_error(fit(max_iter = 0), "positive")
  expect_error(fit(max_iter = 1.5), "integer")
  expect_error(fit(tol = -1), "positive")
  expect_error(fit(control = list(typo = 1)), "Unknown")
  expect_error(fit(control = list(2)), "names")
  expect_error(fit(control = list(tau2_init = NA_real_)), "positive")
  expect_error(fit(control = list(beta_max_iter = NULL)), "cannot be NULL")
  expect_error(fit(approx = "RSPL", control = list(em_max_iter = 0)), "positive")
  expect_error(fit(approx = "Laplace", control = list(eta_max_iter = 0)), "positive")
})

test_that("iteration limits are visible and final moments match final parameters", {
  d <- audit_data()
  expect_warning(short <- glmmFEL(d$y, d$X, d$Z, family = "poisson_log", max_iter = 1), "did not converge")
  expect_false(short$convergence$em_converged)
  expect_warning(inner <- glmmFEL(d$y, d$X, d$Z, family = "poisson_log",
                                   control = list(eta_max_iter = 1)), "did not converge")
  expect_false(inner$convergence$em_converged)
  fit <- glmmFEL(d$y, d$X, d$Z, family = "poisson_log", approx = "FE_full", max_iter = 400)
  expect_true(fit$convergence$em_converged)
  mu <- exp(as.numeric(d$X %*% fit$beta + d$Z %*% fit$eta_mode))
  score <- as.numeric(Matrix::crossprod(d$Z, d$y - mu)) - fit$eta_mode / fit$tau2
  expect_lt(max(abs(score)), 1e-6)
  expect_equal(fit$var_eta_mode, solve(as.matrix(Matrix::crossprod(d$Z, d$Z * mu)) +
                                      diag(1 / fit$tau2, 4)), tolerance = 1e-10)
  fe <- glmmFEL:::fe_trace_diagG(d$Z, fit$var_eta_mode, -mu, -mu)
  expect_equal(fit$eta, fit$eta_mode + fe$trc_y1 / 2, tolerance = 1e-10)
  expect_equal(diag(fit$var_eta), diag(fit$var_eta_mode) + Matrix::diag(fe$trc_y2) / 2, tolerance = 1e-10)
  expect_equal(as.numeric(fitted(fit, "link")), as.numeric(d$X %*% coef(fit) + d$Z %*% fit$eta))
  expect_equal(fitted(fit), exp(fitted(fit, "link")))
  expect_equal(predict(fit), fitted(fit))
  expect_identical(dimnames(vcov(fit)), list(colnames(d$X), colnames(d$X)))
  expect_error(predict(fit, newdata = data.frame(x = 0)), "not supported")
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(attr(logLik(fit), "df"), 3L)
})

test_that("a separating fixed-effects direction cannot be called converged", {
  y <- rep(0:1, each = 10); X <- cbind(1, x = rep(c(-1, 1), each = 10))
  Z <- Matrix::sparseMatrix(i = 1:20, j = rep(1:4, 5), x = 1)
  for (family in c("binomial_logit", "binomial_probit")) {
    for (method in c("FE_full", "MSPL"))
      expect_error(glmmFEL(y, X, Z, family = family, approx = method), "Complete separation")
  }
})

test_that("PL reports both inner and outer convergence and no GLMM likelihood", {
  d <- audit_data()
  for (method in c("MSPL", "RSPL")) {
    expect_warning(short <- glmmFEL(d$y, d$X, d$Z, family = "poisson_log", approx = method,
                                     control = list(pql_max_iter = 1, em_max_iter = 1)), "did not converge")
    expect_false(short$convergence$pql_converged)
    fit <- glmmFEL(d$y, d$X, d$Z, family = "poisson_log", approx = method)
    expect_true(fit$convergence$pql_converged)
    expect_true(fit$convergence$em_converged)
    expect_true(is.na(as.numeric(logLik(fit))))
    expect_true(is.finite(fit$working_logLik))
  }
})

test_that("Gaussian PL core agrees with an independent marginal GLS calculation", {
  d <- audit_data(); Z <- as.matrix(d$Z); X <- d$X
  weights <- seq(0.7, 2, length.out = length(d$y)); tau2 <- 0.6
  V <- diag(1 / weights) + tau2 * tcrossprod(Z)
  Vi <- solve(V); info <- crossprod(X, Vi %*% X)
  beta <- solve(info, crossprod(X, Vi %*% d$y))
  residual <- d$y - X %*% beta
  b <- tau2 * crossprod(Z, Vi %*% residual)
  posterior <- tau2 * diag(ncol(Z)) - tau2^2 * crossprod(Z, Vi %*% Z)
  adjust <- tau2^2 * crossprod(Z, Vi %*% X) %*% solve(info) %*% crossprod(X, Vi %*% Z)
  loglik <- -length(d$y) * log(2 * pi) / 2 - as.numeric(determinant(V, logarithm = TRUE)$modulus) / 2 -
    as.numeric(crossprod(residual, Vi %*% residual)) / 2
  for (method in c("MSPL", "RSPL")) {
    inner <- glmmFEL:::glmmfe_lmm_inner_fit(d$y, weights, X, d$Z, tau2, method)
    expect_equal(unname(inner$beta), as.numeric(beta), tolerance = 1e-10)
    expect_equal(unname(inner$eta), as.numeric(b), tolerance = 1e-10)
    expect_equal(inner$var_eta_post, posterior, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(inner$var_eta_reml, posterior + adjust, tolerance = 1e-10, ignore_attr = TRUE)
    objective <- glmmFEL:::glmmfe_pl_objective(d$y, weights, X, d$Z, inner$beta,
                                              inner$eta, tau2, inner, method)
    reference <- if (method == "MSPL") loglik else
      loglik + ncol(X) * log(2 * pi) / 2 - as.numeric(determinant(info, logarithm = TRUE)$modulus) / 2
    expect_equal(objective, reference, tolerance = 1e-10)
  }
})
