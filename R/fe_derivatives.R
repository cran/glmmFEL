#' Derivative functions for the binomial–probit FE family
#'
#' @description
#' Internal helper that returns a list of derivative functions of the
#' log-likelihood contribution for a probit GLMM, expressed in terms of
#' the transformed linear predictor
#' \deqn{q = (-1)^{1-y} \eta,}
#' where \eqn{\eta} is the linear predictor and \eqn{y \in \{0,1\}}.
#'
#' For a single observation, the contribution to the log-likelihood can be
#' written as
#' \deqn{
#'   \ell(q) = \log \Phi(q),
#' }
#' where \eqn{\Phi} is the standard normal distribution function.
#'
#' @details
#' The returned list has components:
#'
#' * `d0(q)` – \eqn{\log \Phi(q)};
#' * `d1(q)` – first derivative;
#' * `d2(q)` – second derivative w.r.t. \eqn{q};
#' * `d3(q)` – third derivative w.r.t. \eqn{q};
#' * `d4(q)` – fourth derivative w.r.t. \eqn{q}.
#'
#' Implementation notes:
#' * Uses `pnorm(..., log.p=TRUE)` and `dnorm(..., log=TRUE)` to form
#'   the Mills ratio \eqn{\phi(q)/\Phi(q)} stably where possible.
#' * For large negative `q`, differentiates a continued fraction for the
#'   Mills ratio to avoid cancellation in higher derivatives.
#'
#' @return
#' A list with components `d0`, `d1`, `d2`, `d3`, and `d4`, each a function
#' of a numeric argument `q`.
#'
#' @noRd
#' @keywords internal
fe_derivatives_binomial_probit <- function(cutoff = -8) {
  # For t = -q > 0, phi(t)/Phi(-t) has the continued fraction
  # t + 1/(t + 2/(t + 3/(...))). Differentiate the fraction itself;
  # subtracting powers of a large Mills ratio destroys d3/d4 precision.
  tail_derivatives <- function(q) {
    t <- -q
    v <- v1 <- v2 <- v3 <- numeric(length(t))
    for (k in 40:1) {
      g <- t + v; g1 <- 1 + v1; g2 <- v2; g3 <- v3
      v3 <- k * (-6 * g1^3 / g^4 + 6 * g1 * g2 / g^3 - g3 / g^2)
      v2 <- k * (2 * g1^2 / g^3 - g2 / g^2)
      v1 <- -k * g1 / g^2
      v <- k / g
    }
    list(t + v, -(1 + v1), v2, -v3)
  }
  derivative <- function(q, order) {
    q <- as.numeric(q)
    out <- numeric(length(q))
    left <- q < cutoff
    if (any(left)) out[left] <- tail_derivatives(q[left])[[order]]
    if (any(!left)) {
      x <- q[!left]
      m <- exp(stats::dnorm(x, log = TRUE) - stats::pnorm(x, log.p = TRUE))
      out[!left] <- switch(order,
        m,
        -x * m - m^2,
        (x^2 - 1) * m + 3 * x * m^2 + 2 * m^3,
        -(x^3 - 3 * x) * m - (7 * x^2 - 4) * m^2 - 12 * x * m^3 - 6 * m^4)
    }
    out
  }
  list(d0 = function(q) stats::pnorm(q, log.p = TRUE),
       d1 = function(q) derivative(q, 1L),
       d2 = function(q) derivative(q, 2L),
       d3 = function(q) derivative(q, 3L),
       d4 = function(q) derivative(q, 4L))
}

#' Derivative functions for the binomial–logit FE family
#'
#' @description
#' Internal helper that returns a list of derivative functions of the
#' log-likelihood contribution for a logit (logistic) GLMM, expressed in terms of
#' the transformed linear predictor
#' \deqn{q = (-1)^{1-y} \eta,}
#' where \eqn{\eta} is the linear predictor and \eqn{y \in \{0,1\}}.
#'
#' For a single observation, the contribution to the log-likelihood can be
#' written as
#' \deqn{
#'   \ell(q) = \log \operatorname{logit}^{-1}(q) = \log \left(\frac{1}{1 + e^{-q}}\right).
#' }
#'
#' @details
#' The returned list has components:
#'
#' * `d0(q)` – \eqn{\log \operatorname{logit}^{-1}(q)};
#' * `d1(q)` – first derivative;
#' * `d2(q)` – second derivative w.r.t. \eqn{q};
#' * `d3(q)` – third derivative w.r.t. \eqn{q};
#' * `d4(q)` – fourth derivative w.r.t. \eqn{q}.
#'
#' Implementation notes:
#' * Uses `plogis(q)` and `plogis(-q)` to avoid cancellation in the tails.
#' * Uses `plogis(q, log.p = TRUE)` for a stable `d0`.
#'
#' @return
#' A list with components `d0`, `d1`, `d2`, `d3`, and `d4`, each a function
#' of a numeric argument `q`.
#'
#' @noRd
#' @keywords internal
fe_derivatives_binomial_logit <- function() {

  d0 <- function(q) {
    q <- as.numeric(q)
    stats::plogis(q, log.p = TRUE)
  }

  d1 <- function(q) {
    q <- as.numeric(q)
    # d/dq log(plogis(q)) = plogis(-q)
    stats::plogis(-q)
  }

  d2 <- function(q) {
    q <- as.numeric(q)
    p <- stats::plogis(q)
    r <- stats::plogis(-q)
    # -p(1-p) but computed as -p*r for stability
    -(p * r)
  }

  d3 <- function(q) {
    q <- as.numeric(q)
    p <- stats::plogis(q)
    r <- stats::plogis(-q)
    # -p(1-p)(1-2p)
    -(p * r) * (1 - 2 * p)
  }

  d4 <- function(q) {
    q <- as.numeric(q)
    p <- stats::plogis(q)
    r <- stats::plogis(-q)
    # -p(1-p)(1 - 6p + 6p^2)
    -(p * r) * (1 - 6 * p + 6 * p^2)
  }

  list(d0 = d0, d1 = d1, d2 = d2, d3 = d3, d4 = d4)
}



#' Derivative functions for the Poisson–log FE family
#'
#' @description
#' Internal helper that returns a list of derivative functions of the
#' log-likelihood contribution for a Poisson GLMM with log link:
#' \deqn{
#'   \ell(\eta) = y \eta - \exp(\eta) - \log(y!).
#' }
#'
#' @return
#' A list with components `d0`, `d1`, `d2`, `d3`, and `d4`, each a function
#' of a numeric argument `eta` and (for `d0`/`d1`) also the response `y`.
#'
#' @noRd
#' @keywords internal
fe_derivatives_poisson_log <- function() {

  safe_exp <- function(x) {
    x <- as.numeric(x)
    exp(pmin(x, log(.Machine$double.xmax)))
  }

  d0 <- function(eta, y) {
    eta <- as.numeric(eta)
    y <- as.numeric(y)
    y * eta - safe_exp(eta) - lgamma(y + 1)
  }

  d1 <- function(eta, y) {
    eta <- as.numeric(eta)
    y <- as.numeric(y)
    y - safe_exp(eta)
  }

  d2 <- function(eta) {
    eta <- as.numeric(eta)
    -safe_exp(eta)
  }

  d3 <- function(eta) {
    eta <- as.numeric(eta)
    -safe_exp(eta)
  }

  d4 <- function(eta) {
    eta <- as.numeric(eta)
    -safe_exp(eta)
  }

  list(d0 = d0, d1 = d1, d2 = d2, d3 = d3, d4 = d4)
}
