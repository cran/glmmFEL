#' glmmFEL: Generalized Linear Mixed Models via Fully Exponential Laplace in EM
#'
#' @description
#' **glmmFEL** fits generalized linear mixed models (GLMMs) with
#' normal random effects using a matrix-based interface where users supply
#' \eqn{(y, X, Z)} directly. Model fitting is performed with an EM algorithm whose
#' E-step can be approximated using first-order Laplace or fully exponential Laplace
#' (mean-only or mean + covariance corrections), and includes pseudo-likelihood
#' alternatives based on working-response / working-weights linearization.
#'
#' This development branch is intentionally **matrix-first** and supports a flexible
#' (including multiple-membership) random-effects design matrix \eqn{Z}; see the
#' package NOTE for details on current covariance support and planned extensions.
#'
#' Supported families in this branch:
#' \itemize{
#'   \item \code{family = stats::binomial(link = "probit")} (binary probit),
#'   \item \code{family = stats::binomial(link = "logit")} (binary logit),
#'   \item \code{family = stats::poisson(link = "log")} (Poisson log-link).
#' }
#'
#' @template ref-doc
#'
#' @section Approximations:
#' [glmmFEL()] supports:
#' \itemize{
#'   \item \code{"Laplace"}: first-order Laplace approximation,
#'   \item \code{"FE_mean"}: fully exponential Laplace corrections to \eqn{\widehat\eta} only,
#'   \item \code{"FE_full"} (or \code{"FE"}): fully exponential Laplace corrections to both
#'         \eqn{\widehat\eta} and \eqn{\widehat{\mathrm{Var}}(\eta\mid y)},
#'   \item \code{"RSPL"} / \code{"MSPL"}: restricted/marginal pseudo-likelihood
#'         (working response / working weights).
#' }
#'
#' @section Output:
#' [glmmFEL()] returns an object of class \code{"glmmFELMod"} containing:
#' \itemize{
#'   \item \code{beta}: fixed-effect estimates,
#'   \item \code{eta}: empirical Bayes predictions of random effects,
#'   \item \code{tau2}: the scalar variance component,
#'   \item \code{G}: \eqn{q\times q} covariance matrix (diagonal in this branch),
#'   \item \code{var_eta}: prediction-error covariance for \code{eta} (approx.),
#'   \item \code{vcov_beta}: approximate covariance of \code{beta} when available,
#'   \item \code{convergence}: iteration counts and flags.
#' }
#'
#' @seealso
#'  [glmmFEL_pl()] for
#' the pseudo-likelihood engines.
#'
#' @docType package
#' @name glmmFEL-package
"_PACKAGE"
