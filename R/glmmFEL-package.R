#' glmmFEL: Generalized Linear Mixed Models via Fully Exponential Laplace in EM
#'
#' @description
#' **glmmFEL** fits generalized linear mixed models (GLMMs) with
#' normal random effects using a matrix-based interface where users supply
#' \eqn{(y, X, Z)} directly. Model fitting is performed with an EM algorithm whose
#' E-step can be approximated using first-order Laplace or fully exponential Laplace
#' (mean-only or mean and variance-diagonal corrections), and includes pseudo-likelihood
#' alternatives based on working-response / working-weights linearization.
#'
#' The matrix interface supports multiple-membership designs with a single
#' variance component. It does not implement all covariance structures in the
#' cited papers. See [glmmFEL()] for the precise approximation and output limits.
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
#'   \item \code{"FE_full"} (or \code{"FE"}): mean corrections plus the posterior
#'         variance-diagonal corrections needed for the scalar variance update,
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
#'   \item \code{var_eta}: approximate posterior covariance, with diagonal-only
#'         FEL covariance corrections in \code{FE_full},
#'   \item \code{vcov_beta}: approximate covariance of \code{beta} when available,
#'   \item \code{convergence}: iteration counts and flags.
#' }
#'
#' @seealso
#'  [glmmFEL_pl()] for
#' the pseudo-likelihood engines; [glmmFEL-benchmarks] for the simulation appraisal.
#'
#' @docType package
#' @name glmmFEL-package
"_PACKAGE"
