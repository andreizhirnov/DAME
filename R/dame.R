#' @title Distribution-weighted average marginal effect
#' @description
#' \code{dame} computes the distribution-weighted average marginal effects (DAME) as described in Zhirnov, Moral, and Sedashov (2021).
#' @param x a character string representing the name of the main variable of interest. Marginal effects will be computed for this variable.
#' @param over a character string representing the name of the conditionning variable. DAME will be computed for the bins long the range of this variable.
#' @param model fitted model object. The package works best with GLM objects and will extract the formula, dataset, family, coefficients, and
#' the QR components of the design matrix if arguments \code{formula}, \code{data}, \code{link}, \code{coefficients}, and/or
#' \code{vcov} are not explicitly specified.
#' @param data the dataset to be used to compute marginal effects (if not specified, it is extracted from the fitted model object).
#' @param formula the formula used in estimation (if not specified, it is extracted from the fitted model object).
#' @param link the name of the link function used in estimation (if not specified, it is extracted from the fitted model object).
#' @param coefficients the named vector of coefficients produced during the estimation (if not specified, it is extracted from the fitted model object).
#' @param vcov the variance-covariance matrix to be used for computing standard errors (if not specified, it is extracted from the fitted model object).
#' @param nbins the number of bins to be used for aggregating marginal effects; the default is 10 bins of equal size; ignored if \code{bin_id} is specified or
#' \code{use_distinct_values} is TRUE.
#' @param bin_id a numeric vector identifying the bins used for aggregating marginal effects (if not specified and \code{use_distinct_values=FALSE},
#' the function uses \code{nbins} bins with roughly equal number of observations; if not specified and \code{use_distinct_values=TRUE}, the function
#' uses all unique values of the \code{over} variable).
#' @param use_distinct_values logical; if TRUE, the function uses all unique values of the \code{over} variable; ignored if \code{bin_id} is specified.
#' @param discrete logical. If TRUE, the function will compute the effect of a discrete change in \code{x}. If FALSE, the function will compute the partial derivative of \code{x}.
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
#' @author Function \code{dame} is an implementation of a procedure described in Zhirnov, Moral, and  Sedashov (2021).
#' Standard errors are computed using either the delta method (Greene 2012) for more details) or Monte-Carlo simulations (King, Tomz, and Wittenberg 2000).
#' @references
#' Greene, William. 2012. \emph{Econometric Analysis, 7 ed.} Pearson Education Limited.
#'
#' King, Garry, Michael Tomz, and Jason Wittenberg. 2000. ``Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' \emph{American Journal of Political Science} 44(2): 341-355.
#'
#' Zhirnov, Andrei, Mert Moral, and Evgeny Sedashov (2021). ``Taking Distributions Seriously: On the Interpretation
#' of the Estimates of Interactive Nonlinear Models.'' Working paper.
#' @return \code{dame} returns a data frame with the estimates of the distribution-weighted average marginal effects, standard errors, confidence intervals,
#' the corresponding bin IDs (by default, the median value of the conditioning variable within the bin), and the targeted combinations of
#' \code{at} values if specified.
#' @examples
#' ##poisson regression with 2 variables and an interaction between them
#' #fit the regression first
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5),
#' x1 = rpois(10000, 3), w=c("a","b","c","d"))
#' y <- glm(y ~ x1*x2 + w, data = data, family = "poisson")
#' #compute DAME
#' dame(model = y, x = "x1", over = "x2")
#' \dontrun{
#' ## logit
#' m <- glm(any_dispute ~ flows.ln*polity2 + gdp_pc, data=strikes, family="binomial")
#' summary(m)
#' ## DAME with a robust (heteroscedasticity-consistent) variance-covariance matrix and 4 bins
#' library(sandwich)
#' dame(model=m, x="flows.ln", over="polity2", nbins=4, vcov=vcovHC(m))
#'}
#' @export

dame <- function(x, over = NULL, model = NULL,
                 data = NULL, formula = NULL, link = NULL,
                 coefficients = NULL, vcov = NULL,
                 nbins = 10, bin_id = NULL, use_distinct_values = FALSE,
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

# make a data frame specific to DAME
  obj <- list(data=data)
  if (is.null(obj[["data"]])) obj[["data"]] <- eval(model)[["data"]]
  check.required("data","data.frame", list=obj)

  obj[["bin_id"]] <- eval(bin_id)
  if (!is.null(obj[["bin_id"]]) && !is.numeric(obj[["bin_id"]])) stop("'bin_id' must be a numeric vector", call. = FALSE)
  if (!is.null(obj[["bin_id"]]) && length(obj[["bin_id"]]) != nrow(obj[["data"]])) stop("'bin_id' must have the same length as the dataset", call. = FALSE)
  if (is.null(obj[["bin_id"]])) {
    if (inherits(over, "character") && use_distinct_values) {
      obj[["bin_id"]] <- obj[["data"]][[over]]
    } else if (inherits(over, "character") && !is.null(obj[["data"]][[over]])) {
      obj[["bin_id"]] <- make.bins(obj[["data"]][[over]], nbins)
    } else {
      obj[["bin_id"]] <- rep(1, nrow(obj[["data"]]))
    }
  }
  calc[["formula"]] <- formula
  if (is.null(calc[["formula"]])) calc[["formula"]] <- eval(model)[["formula"]]
  calc[["formula"]][[2L]] <- NULL
  check.required("formula","formula", list=calc)

  obj[["allvars"]] <- all.vars(calc[["formula"]])

  calc[["x"]] <- x
  check.required("x","character", list=calc)
  outside.formula <- setdiff(c(calc[["x"]],over,names(at)),obj[["allvars"]])
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x and over variables are included in the data
  outside.data <- setdiff(c(calc[["x"]],over),names(obj[["data"]]))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,sep="\n"), call. = FALSE)

  if (length(at)>0) {
    obj[["at"]] <- as.list(at)
    for (v in names(obj[["at"]])) {
      if (is.character(obj[["at"]][[v]]) & !is.factor(obj[["at"]][[v]])) {
        xle <-  model[["xlevels"]][[v]]
        if (is.null(xle)) xle <- sort(unique(obj[["data"]][[v]]))
        if (is.null(xle)) {
          stop("Please convert the character variables in the 'at' list into factors", call. = FALSE)
        }
        if (any(!obj[["at"]][[v]] %in% xle)) {
          stop(paste0("Could not find all listed values of ",v," in the model"), call. = FALSE)
        }
        obj[["at"]][[v]] <- factor(obj[["at"]][[v]], levels=xle)
      }
    }
  }

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
  if (nrow(calc[["grid"]]) > 0) effects <- cbind(effects, calc[["grid"]])
  rownames(effects) <- c()
  return(effects)
}
