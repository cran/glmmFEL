##### glmmFEL-methods.R ######################################################
# S3 methods for the simplified, matrix-only glmmFEL development branch.
#
# This file intentionally avoids any dependency on lme4.  The object returned
# by glmmFEL() is a plain list with class "glmmFELMod", and we provide a small
# set of base-compatible methods (print/summary/coef/vcov/fitted/predict/logLik).
################################################################################

#' Print a glmmFEL model object
#'
#' @param x A `glmmFELMod` object.
#' @param ... Unused.
#' @return Returns `x` invisibly (a `glmmFELMod` object), called for its side
#'   effect of printing model information to the console.
#' @export
print.glmmFELMod <- function(x, ...) {
  cat("glmmFEL model fit\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }

  cat("\nFamily:", x$family, "\n")
  cat("Approximation:", x$approx, "\n")

  if (is.finite(x$logLik)) {
    cat("logLik (approx):", format(x$logLik, digits = 6), "\n")
  }

  cat("\nFixed effects (beta):\n")
  b <- x$beta
  if (!is.null(colnames(x$X))) names(b) <- colnames(x$X)
  print(b)

  cat("\nRandom-effects variance component:\n")
  cat("  tau2 =", format(x$tau2, digits = 6), "\n")

  invisible(x)
}

#' Summary for a glmmFEL model object
#'
#' @param object A `glmmFELMod` object.
#' @param ... Unused.
#' @return An object of class `summary.glmmFELMod` (a list) containing summary
#'   information for the fitted model, including the call, family,
#'   approximation label, approximate log-likelihood, estimated variance
#'   component `tau2`, a coefficient table for fixed effects (and standard
#'   errors when available), number of observations, and convergence
#'   information.
#' @export
summary.glmmFELMod <- function(object, ...) {
  b <- object$beta
  if (!is.null(colnames(object$X))) names(b) <- colnames(object$X)

  se <- NULL
  if (!is.null(object$vcov_beta)) {
    se <- sqrt(diag(object$vcov_beta))
    names(se) <- names(b)
  }

  tab <- NULL
  if (!is.null(se)) {
    tab <- cbind(Estimate = b, `Std. Error` = se)
  } else {
    tab <- cbind(Estimate = b)
  }

  res <- list(
    call    = object$call,
    family  = object$family,
    approx  = object$approx,
    logLik  = object$logLik,
    tau2    = object$tau2,
    coef    = tab,
    nobs    = length(object$y),
    converg = object$convergence
  )

  class(res) <- "summary.glmmFELMod"
  res
}

#' Print a summary.glmmFELMod object
#'
#' @param x A `summary.glmmFELMod` object.
#' @param ... Unused.
#' @return Returns `x` invisibly (a `summary.glmmFELMod` object), called for its
#'   side effect of printing the summary to the console.
#' @export
print.summary.glmmFELMod <- function(x, ...) {
  cat("glmmFEL model summary\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }

  cat("\nFamily:", x$family, "\n")
  cat("Approximation:", x$approx, "\n")

  if (is.finite(x$logLik)) {
    cat("logLik (approx):", format(x$logLik, digits = 6), "\n")
  }

  cat("\nFixed effects:\n")
  print(x$coef)

  cat("\nRandom-effects variance component:\n")
  cat("  tau2 =", format(x$tau2, digits = 6), "\n")

  if (!is.null(x$converg)) {
    cat("\nConvergence:\n")
    print(x$converg)
  }

  invisible(x)
}



#' Extract model coefficients (fixed effects)
#'
#' @param object A `glmmFELMod` object.
#' @param ... Unused.
#' @return A named numeric vector of estimated fixed-effect regression
#'   coefficients (on the linear predictor scale). Names correspond to columns
#'   of the fixed-effects design matrix `X` when available.
#' @export
coef.glmmFELMod <- function(object, ...) {
  b <- object$beta
  if (!is.null(colnames(object$X))) names(b) <- colnames(object$X)
  b
}

#' Extract the covariance matrix of the fixed effects
#'
#' @param object A `glmmFELMod` object.
#' @param ... Unused.
#' @return A variance-covariance matrix for the estimated fixed-effect
#'   regression coefficients. Row and column names correspond to coefficient
#'   names when available.
#' @export
vcov.glmmFELMod <- function(object, ...) {
  object$vcov_beta
}

#' Extract fitted values
#'
#' @param object A `glmmFELMod` object.
#' @param type Either `"link"` (linear predictor) or `"response"`.
#' @param ... Unused.
#' @return A numeric vector of fitted values. If `type = "link"`, the fitted
#'   linear predictor values are returned. If `type = "response"`, the fitted
#'   mean response values on the response scale are returned. The length equals
#'   the number of observations in the fitted object.
#' @export
fitted.glmmFELMod <- function(object, type = c("response", "link"), ...) {
  type <- match.arg(type)

  eta_lin <- as.vector(object$X %*% object$beta + object$Z %*% object$eta)

  if (type == "link") return(eta_lin)

  if (identical(object$family, "binomial_probit")) {
    return(stats::pnorm(eta_lin))
  }
  if (identical(object$family, "binomial_logit")) {
    return(stats::plogis(eta_lin))
  }
  if (identical(object$family, "poisson_log")) {
    return(exp(eta_lin))
  }

  stop("Unsupported family in fitted.glmmFELMod().")
}

#' Predict from a fitted glmmFEL model
#'
#' @param object A `glmmFELMod` object.
#' @param newdata Not supported in this branch (matrix interface only).
#' @param type Either `"link"` or `"response"`.
#' @param ... Unused.
#' @return A numeric vector of predictions. If `type = "link"`, predictions are
#'   returned on the linear predictor scale. If `type = "response"`, predictions
#'   are returned on the response scale. The length equals the number of
#'   observations used to fit the model. Supplying `newdata` triggers an error
#'   in this matrix-only branch.
#' @export
predict.glmmFELMod <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  if (!is.null(newdata)) {
    stop("predict(..., newdata=) is not supported in the matrix-only branch. ",
         "Supply X/Z directly and refit if you need out-of-sample prediction.")
  }
  fitted.glmmFELMod(object, type = type)
}

#' Extract log-likelihood (approximate)
#'
#' @param object A `glmmFELMod` object.
#' @param ... Unused.
#' @return An object of class `"logLik"` giving the approximate log-likelihood
#'   stored in `object$logLik`, with attributes `"df"` (effective number of
#'   parameters, taken as `length(beta) + 1` for the variance component) and
#'   `"nobs"` (number of observations).
#' @export
logLik.glmmFELMod <- function(object, ...) {
  ll <- as.numeric(object$logLik)
  attr(ll, "df") <- length(object$beta) + 1L
  attr(ll, "nobs") <- length(object$y)
  class(ll) <- "logLik"
  ll
}

##### glmmFEL-methods.R #######################################################
