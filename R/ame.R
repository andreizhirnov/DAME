#' @title AME function
#' @description
#' \code{ame} computes the average marginal effects of variable \code{x} at the specified values of \code{at} variables.
#' @param x a character string representing the name of the main variable of interest. Marginal effects will be computed for this variable.
#' @param model fitted model object. The package works best with GLM objects and will extract the formula, dataset, family, coefficients, and
#' the QR components of the design matrix if arguments \code{formula}, \code{data}, \code{link}, \code{coefficients}, and/or
#' \code{vcov} are not explicitly specified.
#' @param data the dataset to be used to compute marginal effects (if not specified, it is extracted from the fitted model object).
#' @param formula the formula used in estimation (if not specified, it is extracted from the fitted model object).
#' @param link the name of the link function used in estimation (if not specified, it is extracted from the fitted model object).
#' @param coefficients the named vector of coefficients produced during the estimation (if not specified, it is extracted from the fitted model object).
#' @param vcov the variance-covariance matrix to be used for computing standard errors (if not specified, it is extracted from the fitted model object).
#' @param discrete A logical variable. If TRUE, the function will compute the effect of a discrete change in \code{x}. If FALSE, the function will compute the partial derivative of \code{x}.
#' @param discrete_step The size of a discrete change in \code{x} used in computations (used only if \code{discrete=TRUE}).
#' @param at an optional named list of values of independent variables. These variables will be set to these value before computations. 
#' The remaining numeric variables (except \code{x} and \code{over}) will be set to their means. The remaining factor variables will be set
#' to their modes.
#' @param mc logical. If TRUE, the standard errors and confidence intervals will be computed using simulations. 
#' If FALSE (default), the delta method will be used.
#' @param iter the number of interations used in Monte-Carlo simulations. Default = 1,000. 
#' @param pct a numeric vector with the quantiles to be output with the DAME estimates. Default = \code{c(2.5,97.5)}. 
#' @param weights an optional vector of sampling weights.
#' @return \code{ame} returns a data frame with the estimates of the average marginal effects, standard errors, confidence intervals, 
#' and the used values of the independent variables.
#' @export

ame <- function(x, model = NULL, data = NULL, formula = NULL, link = NULL,
               coefficients = NULL, vcov = NULL,
               discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
               pct = c(2.5, 97.5), iter = 1000, weights = NULL) {


# compute the derivatives
  link <- link[1]
  if (is.null(link)) link <- eval(model)[["family"]][["link"]]
  check.required("link","character")
  if (!(link %in% c("logit","probit","cauchit","cloglog","identity","log","sqrt","1/mu^2","inverse"))) {
    stop("Invalid link name. Valid links include 'logit','probit','cauchit','cloglog','identity','log','sqrt','1/mu^2','inverse'", call. = FALSE)
  }
  calc <- make.dydm(link=link)

# make a data frame specific to AME
  obj <- list(data=data)
  if (is.null(obj[["data"]])) obj[["data"]] <- eval(model)[["data"]]
  check.required("data","data.frame", list=obj)

  obj[["bin_id"]] <- rep(1, nrow(obj[["data"]]))

  calc[["formula"]] <- formula
  if (is.null(calc[["formula"]])) calc[["formula"]] <- eval(model)[["formula"]]
  calc[["formula"]][[2L]] <- NULL
  check.required("formula","formula", list=calc)

  obj[["allvars"]] <- all.vars(calc[["formula"]])

  calc[["x"]] <- x
  check.required("x","character", list=calc)
  outside.formula <- setdiff(c(calc[["x"]],names(at)),obj[["allvars"]])
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x is included in the data
  outside.data <- setdiff(calc[["x"]],names(obj[["data"]]))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,sep="\n"), call. = FALSE)

  if (length(at)>0) obj[["at"]] <- as.list(at)

  obj[["weights"]] <- eval(weights)
  if (!is.null(obj[["weights"]]) && !is.numeric(obj[["weights"]])) stop("'weights' must be a numeric vector", call. = FALSE)
  if (!is.null(obj[["weights"]]) && length(obj[["weights"]]) != nrow(obj[["data"]])) stop("'weights' must have the same length as the dataset", call. = FALSE)
  if (is.null(obj[["weights"]])) obj[["weights"]] <- rep(1, nrow(obj[["data"]]))

  calc <- c(calc, do.call("makeframes.dame", obj))

  # computation
  calc[["discrete"]] <- discrete
  calc[["discrete_step"]] <- discrete_step
  calc[["coefficients"]] <- coefficients
  if (is.null(calc[["coefficients"]])) calc[["coefficients"]] <- stats::coef(model)
  check.required("coefficients", "numeric", list=calc)

  calc[["vcov"]] <- vcov
  if (is.null(calc[["vcov"]])) calc[["vcov"]] <- stats::vcov(model)
  check.required("vcov", "matrix", list=calc)

  calc[["pct"]] <- pct
  check.required("pct", "numeric", list=calc)
  names(calc[["pct"]]) <- paste0("p",pct)
  if (any(calc[["pct"]] > 100) || any(calc[["pct"]] <0)) stop("Error: 'pct' must be between 0 and 100", call. = FALSE)

  if (mc) {
    calc[["iter"]] <- as.integer(iter)
    if (calc[["iter"]] < 1) stop("Error: 'iter' must be positive.", call. = FALSE)
    effects <- do.call("simulated.me", calc)
  } else {
    effects <- do.call("analytical.me", calc)
  }
  # merge with other variables
  if (nrow(calc[["grid"]]) > 0) effects <- cbind(effects, calc[["grid"]])
  effects[["bin_id"]] <- NULL
  rownames(effects) <- c()
  return(effects)
}
