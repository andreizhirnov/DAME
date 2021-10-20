#' @title AME function
#' @description
#' \code{ame} computes the average marginal effects of variable \code{x} at the specified values of \code{at} variables.
#' @param x A character string representing the name of a variable. ME will be computed for this variable.
#' @param model Fitted model object. The package works best with GLM objects and will extract the formula, dataset, family, coefficients, and
#' the QR components of the design matrix if arguments \code{formula}, \code{data}, \code{link}, \code{coefficients}, and/or
#' \code{variance} are not explicitly specified.
#' @param data The dataset to be used to compute marginal effects. DAME provides most meaningful results if the original dataset is used for computing marginal effects.
#' @param formula The formula used in estimation (if not specified, it is extracted from the fitted model object).
#' @param link The link function used in estimation (if not specified, it is extracted from the fitted model object).
#' @param coefficients The named vector of coefficients produced during the estimation (if not specified, it is extracted from the fitted model object).
#' @param variance The variance-covariance matrix to be used for computing standard errors (if not specified, it is extracted from the fitted model object).
#' @param discrete A logical variable. If TRUE, the function will compute the effect of a discrete change in \code{x}. If FALSE, the function will compute the partial derivative of \code{x}.
#' @param discrete_step The size of a discrete change in \code{x} used in computations (used only if \code{discrete=TRUE}).
#' @param at A named list of values of independent variables. These variables will be set to these value before computations. All other quantitative variables (except \code{x} and \code{over}) will be set to their means. All other factor variables will be set to their modes.
#' @param mc A logical variable. If TRUE, the function will compute standard errors and sampling quantiles using Monte-Carlo simulations. If FALSE, the function will use the delta method.
#' @param pct A numeric vector with the sampling quantiles to be output with the DAME estimates. Default = \code{c(2.5,97.5)}.
#' @param iter Number of iterations used in Monte-Carlo simulations. Default = 1,000.
#' @param weights an optional vector of sampling weights.
#' @return \code{ame} returns a data frame with the estimates of the average marginal effects, standard errors, confidence intervals, and the used values of the independent variables.
#' @export

ame <- function(x, model = NULL, data = NULL, formula = NULL, link = NULL,
               coefficients = NULL, variance = NULL,
               discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
               pct = c(2.5, 97.5), iter = 1000, weights = NULL) {

# extract arguments from the call
  args <- as.list(match.call())
  if (!("formula" %in% names(args))) args[["formula"]] <- eval(args[["model"]])[["formula"]]
  if (!("data" %in% names(args))) args[["data"]] <- eval(args[["model"]])[["data"]]
  if (!("link" %in% names(args))) args[["link"]] <- eval(args[["model"]])[["family"]][["link"]]
  if (!("coefficients" %in% names(args))) args[["coefficients"]] <- stats::coef(eval(args[["model"]]))
  if (!("variance" %in% names(args))) args[["variance"]] <- stats::vcov(eval(args[["model"]]))
  if (is.null(weights)) weights <- rep(1, nrow(data))
  # check the required arguments and coerce the specified arguments into a proper class
  checks <- list(
    required=c("x","formula","data","link","coefficients","variance"),
    types = list(x ="character", data = "data.frame", link = "character", formula  = "formula",
                 coefficients= "numeric", variance = "matrix", discrete = "logical", discrete_step = "numeric",
                 at= "list", mc = "logical", pct = "numeric", iter = "integer"),
    lengths = list(x = 1L, link = 1L, discrete = 1L, discrete_step = 1L, mc = 1L, iter = 1L)
  )
  check.args(args=args, checks=checks)

  # check if x and at variables are included in the formula
  updform <- formula
  updform[[2L]] <- NULL
  allvars <- all.vars(updform)
  outside.formula <- setdiff(c(x,names(at)),allvars)
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x variable is included in the data
  outside.data <- setdiff(x,names(data))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,sep="\n"), call. = FALSE)

  # misc checks
  if (!is.null(pct)) {
    for (p in pct) {
      if (p > 100 | p <= 0) stop("Error: 'pct' must be between 0 and 100", call. = FALSE)
    }
  }
  names(pct) <- paste0("p",pct)
  if (!is.null(link)) {
    if (!(link %in% c("logit","probit","cauchit","cloglog","identity","log","sqrt","1/mu^2","inverse"))) {
      stop("Invalid link name. Valid links include 'logit','probit','cauchit','cloglog','identity','log','sqrt','1/mu^2','inverse'", call. = FALSE)
    }
  }
  dyli <- make.dydm(link=link)
## make a data frame specific to AME
  mfli <- makeframes.dame(data=data,allvars=allvars,at=at,bin_id=rep(1,nrow(data)))
##  computation
  if (mc) {
    effects <- simulated.me(discrete=discrete, discrete_step=discrete_step, iter=iter, coefficients=args[["coefficients"]], variance=args[["variance"]],
                              data=mfli[["data.compressed"]], x = x, formula=updform, ym=dyli$ym, mx=mx, dydm=dyli$dydm, wmat=mfli$wmat, pct=pct)
  } else {
    effects <- analytical.me(discrete=discrete, discrete_step=discrete_step, coefficients=args[["coefficients"]], variance=args[["variance"]],
                               data=mfli[["data.compressed"]], x = x, formula=updform, ym=dyli$ym, mx=mx, dydm=dyli$dydm, d2ydm2=dyli$d2ydm2, wmat=mfli$wmat, pct=pct)
  }
  ## merge in other variables
  if (length(at) > 0) {
    mfli[["grid"]][["bin_id"]] <- NULL
    effects <- cbind(effects,mfli$grid)
    }
  rownames(effects) <- c()
  return(effects)
}
