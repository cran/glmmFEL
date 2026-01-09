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
#' * For large negative `q`, uses a Mills-ratio asymptotic expansion to
#'   avoid overflow/underflow.
#' * Higher derivatives are computed via stable recurrences in terms of
#'   the Mills ratio.
#'
#' @return
#' A list with components `d0`, `d1`, `d2`, `d3`, and `d4`, each a function
#' of a numeric argument `q`.
#'
#' @noRd
#' @keywords internal
fe_derivatives_binomial_probit <- function(cutoff = -10) {

  # asymptotic Mills ratio for q -> -Inf:
  # m(q) = phi(q)/Phi(q) ~ -q - 1/q - 2/q^3 - 5/q^5 - 14/q^7 - 42/q^9
  mills_asym <- function(q) {
    q <- as.numeric(q)
    inv  <- 1 / q
    inv3 <- inv^3
    inv5 <- inv^5
    inv7 <- inv^7
    inv9 <- inv^9
    m <- (-q) - inv - 2 * inv3 - 5 * inv5 - 14 * inv7 - 42 * inv9
    pmax(m, 1e-300)  # keep positive, avoid log(0)
  }

  log_mills <- function(q) {
    q <- as.numeric(q)

    # If far in the left tail, use asymptotic directly (avoids any exp overflow downstream)
    out <- rep(NA_real_, length(q))
    left <- (q < cutoff)

    if (any(!left)) {
      qq <- q[!left]
      logphi <- stats::dnorm(qq, log = TRUE)
      logPhi <- stats::pnorm(qq, log.p = TRUE)

      out0 <- logphi - logPhi

      # Guard against:
      #  - non-finite out0 (logPhi = -Inf)
      #  - exp(out0) overflow later (out0 too large)
      too_big <- out0 > log(.Machine$double.xmax)
      bad <- !is.finite(out0) | too_big

      if (any(bad)) {
        qb <- qq[bad]
        out0[bad] <- log(mills_asym(qb))
      }

      out[!left] <- out0
    }

    if (any(left)) {
      out[left] <- log(mills_asym(q[left]))
    }

    out
  }

  mills <- function(q) {
    lm <- log_mills(q)

    # Prevent exp overflow turning into Inf silently; if it would overflow, return +Inf.
    lm_cap <- log(.Machine$double.xmax)
    out <- exp(pmin(lm, lm_cap))
    out
  }

  # Stable logPhi: use pnorm(log.p=TRUE) where possible; otherwise use logphi - log(m)
  logPhi_stable <- function(q) {
    q <- as.numeric(q)
    logPhi <- stats::pnorm(q, log.p = TRUE)

    bad <- !is.finite(logPhi) | (q < cutoff)
    if (any(bad)) {
      qb <- q[bad]
      logphi_b <- stats::dnorm(qb, log = TRUE)
      logm_b   <- log(mills_asym(qb))
      logPhi[bad] <- logphi_b - logm_b
    }

    logPhi
  }

  d0 <- function(q) {
    logPhi_stable(q)
  }

  d1 <- function(q) {
    mills(q)
  }

  # Recurrences in terms of m = mills(q)
  d2 <- function(q) {
    q <- as.numeric(q)
    m <- mills(q)
    -q * m - m^2
  }

  d3 <- function(q) {
    q <- as.numeric(q)
    m <- mills(q)
    (q^2 - 1) * m + 3 * q * m^2 + 2 * m^3
  }

  d4 <- function(q) {
    q <- as.numeric(q)
    m <- mills(q)
    -(q^3 - 3 * q) * m - (7 * q^2 - 4) * m^2 - 12 * q * m^3 - 6 * m^4
  }

  list(d0 = d0, d1 = d1, d2 = d2, d3 = d3, d4 = d4)
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
