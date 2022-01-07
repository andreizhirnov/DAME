#' @title Marginal effects at means
#' @description
#' \code{mem} computes the marginal effects of variable \code{x} at the specified values of \code{at} variables and the mean values of all the other variables (including \code{x}).
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
#' @param pct a named numeric vector with the sampling quantiles to be output with the DAME estimates (the names are used as the new variable names).
#' Default = \code{c(lb=2.5,ub=97.5)}.
#' @param weights an optional vector of sampling weights.
#' @return \code{mem} returns a data frame with the estimates of the marginal effects at means, along with standard errors, confidence intervals,
#' and the used values of the independent variables. All quantitative variable not included in \code{at} are set to their means,
#'  and all qualitative variables (except those listed in \code{at}) are converted to factors and set to their modes.
#' @examples
#' ##poisson regression with 2 variables and an interaction between them
#' #fit the regression first
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5),
#' x1 = rpois(10000, 3), w=c("a","b","c","d"))
#' y <- glm(y ~ x1*x2 + w, data = data, family = "poisson")
#' #compute marginal effects at means
#' mem(model = y, x = "x1",
#'      at=list(x2=seq(min(data$x2), max(data$x2), length.out=5)))
#' \dontrun{
#' ## logit
#' m <- glm(any_dispute ~ flows.ln*polity2 + gdp_pc, data=strikes, family="binomial")
#' summary(m)
#' ## marginal effects at means with a robust (heteroscedasticity-consistent)
#' variance-covariance matrix
#' library(sandwich)
#' mem(model=m, x="flows.ln", vcov=vcovHC(m), at=list(polity2=c(-10,0,10)))
#' }
#' @export


## ! is it jsut a wrapper of me() when over=NULL and x is set to its mean?
mem <- function(x, model = NULL, data = NULL, formula = NULL, link = NULL,
                coefficients = NULL, vcov = NULL,
                discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
                pct = c(lb=2.5, ub=97.5), iter = 1000, weights = NULL) {

  # compute the derivatives
  link <- link[1]
  if (is.null(link)) link <- eval(model)[["family"]][["link"]]
  check.required("link","character")

  if (!(link %in% c("logit","probit","cauchit","cloglog","identity","log","sqrt","1/mu^2","inverse"))) {
    stop("Invalid link name. Valid links include 'logit','probit','cauchit','cloglog','identity','log','sqrt','1/mu^2','inverse'", call. = FALSE)
  }
  calc <- make.dydm(link=link)

  # make a data frame specific to MEM
  if (is.null(data)) data <- eval(model)[["data"]]
  check.required("data","data.frame")

  calc[["formula"]] <- formula
  if (is.null(calc[["formula"]])) calc[["formula"]] <- stats::formula(model)
  calc[["formula"]][[2L]] <- NULL
  check.required("formula","formula", list=calc)

  allvars <- all.vars(calc[["formula"]])
  tomeans <- setdiff(allvars, names(at))
  names(tomeans) <- tomeans

  at <- as.list(at)
  if (length(at)>0) {
    for (v in names(at)) {
      if (is.character(at[[v]]) & !is.factor(at[[v]])) {
        xle <-  model[["xlevels"]][[v]]
        if (is.null(xle)) xle <- sort(unique(data[[v]]))
        if (is.null(xle)) {
          stop("Please convert the character variables in the 'at' list into factors", call. = FALSE)
        }
        if (any(!at[[v]] %in% xle)) {
          stop(paste0("Could not find all listed values of ",v," in the model"), call. = FALSE)
        }
        at[[v]] <- factor(at[[v]], levels=xle)
      }
    }
  }

  if (length(tomeans)>0) at <- c(at, lapply(tomeans, find.central, data=data, weights=weights))

  calc[["data"]] <- makeframes.mem(at)

  ## calculations
  calc[["x"]] <- x
  check.required("x","character", list=calc)
  outside.formula <- setdiff(c(calc[["x"]], names(at)),allvars)
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x and over variables are included in the data
  outside.data <- setdiff(calc[["x"]],names(data))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,sep="\n"), call. = FALSE)

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
  if (is.null(names(calc[["pct"]]))) {
	names(calc[["pct"]]) <- paste0("p",pct)
	} else {
      names(calc[["pct"]]) <- make.names(names(calc[["pct"]]))
	}
  if (any(calc[["pct"]] > 100) || any(calc[["pct"]] <0)) stop("Error: 'pct' must be between 0 and 100", call. = FALSE)

  if (mc) {
    calc[["iter"]] <- as.integer(iter)
    if (calc[["iter"]] < 1) stop("Error: 'iter' must be positive.", call. = FALSE)
    effects <- do.call("simulated.me", calc)
  } else {
    effects <- do.call("analytical.me", calc)
  }
  # merge with other variables
  if (nrow(calc[["data"]]) > 0) effects <- cbind(effects, calc[["data"]])
  rownames(effects) <- c()
  return(effects)
}
