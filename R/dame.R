#' @title DAME function
#' @description
#' \code{dame} computes the distribution-weighted average marginal effects (DAME) as described in Moral, Sedashov, and Zhirnov (2017).
#' @param x A character string representing the name of a variable. DAME will be computed for this variable.
#' @param over A character string representing the name of a variable. DAME will be computed for each of the values of this variable.
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
#' @param iter Number of interations used in Monte-Carlo simulations. Default = 1,000.
#' @author Function \code{dame} comes as a complimentary open-source implementation of procedures described in Moral, Sedashov, and Zhirnov (2017).
#' Standard errors are computed using either delta method (Greene 2012) for more details) or Monte-Carlo simulations (King, Tomz, and Wittenberg 2000).
#' @references
#' Greene, William. 2012. \emph{Econometric Analysis, 7 ed.} Pearson Education Limited.
#'
#' King, Garry, Michael Tomz, and Jason Wittenberg. 2000. ``Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' \emph{American Journal of Political Science} 44(2): 341-355.

#' Moral, Mert, Evgeny Sedashov, and Andrei Zhirnov. (2017) ``Taking Distributions Seriously Interpreting the Effects of Constitutive Variables in Nonlinear Models with Interactions.'' Working paper.
#' @return A list of the following:
#' \itemize{
#' \item\code{dame} A data frame with the DAME estimates, standard errors, quantiles of the sampling distribution, and the values of the independent variables used for computing DAME.
#' \item\code{execute_time} Execution time
#' }
#' @examples
#' ##poisson regression with 2 variables and an interaction between them
#' #fit the regression first
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
#' y <- glm(y ~ x1*x2, data = data, family = "poisson")
#' #compute DAME
#' dame(model = y, x = "x1", over = "x2")
#' @export

dame <- function(x, over = NULL, model = NULL, data = NULL, formula = NULL, link = NULL,
                 coefficients = NULL, variance = NULL,
                 discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
                 pct = c(2.5, 97.5), iter = 1000) {
  start_time <- Sys.time()

  # extract arguments from the call
  args <- as.list(match.call())
  if (!("formula" %in% names(args))) args[["formula"]] <- eval(args[["model"]])[["formula"]]
  if (!("data" %in% names(args))) args[["data"]] <- eval(args[["model"]])[["data"]]
  if (!("link" %in% names(args))) args[["link"]] <- eval(args[["model"]])[["family"]][["link"]]
 # if (!("link" %in% names(args))) args[["link"]] <- quote(eval(args[["model"]])[["family"]][["link"]])
  if (!("coefficients" %in% names(args))) args[["coefficients"]] <- stats::coef(eval(args[["model"]]))
  if (!("variance" %in% names(args))) args[["variance"]] <- stats::vcov(eval(args[["model"]]))

# check the required arguments and coerce the specified arguments into a proper class
  checks <- list(
    required=c("x","formula","data","link","coefficients","variance"),
    types = list(x ="character", over = "character", data = "data.frame", link = "character", formula  = "formula",
                 coefficients= "numeric", variance = "matrix", discrete = "logical", discrete_step = "numeric",
                 at= "list", mc = "logical", pct = "numeric", iter = "integer"),
    lengths = list(x = 1L, link = 1L, discrete = 1L, discrete_step = 1L, mc = 1L, iter = 1L)
    )
  check.args(args=args, checks=checks)

  # check if x, over, at variables are included in the formula
  updform <- args[["formula"]]
  updform[[2L]] <- NULL
  allvars <- all.vars(updform)
  outside.formula <- setdiff(c(x,over,names(at)),allvars)
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x and over variables are included in the data
    outside.data <- setdiff(c(x,over),names(data))
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
# preliminaries
  dyli <- make.dydm(link=link)
# make a data frame specific to DAME
  mfli <- makeframes.dame(data=args[["data"]],allvars=allvars,at=args[["at"]],over=args[["over"]],x=args[["x"]])
# computation
  if (mc) {
  to_insert <- simulated.me(discrete=discrete, discrete_step=discrete_step, iter=iter, coefficients=args[["coefficients"]], variance=args[["variance"]],
                           data=mfli[["data.compressed"]], x = x, formula=updform, ym=dyli[["ym"]], mx=mx, dydm=dyli[["dydm"]], wmat=mfli[["wmat"]], pct=pct)
  } else {
  to_insert <- analytical.me(discrete=discrete, discrete_step=discrete_step, coefficients=args[["coefficients"]], variance=args[["variance"]],
                                 data=mfli[["data.compressed"]], x = x, formula=updform, ym=dyli[["ym"]], mx=mx, dydm=dyli[["dydm"]], d2ydm2=dyli[["d2ydm2"]], wmat=mfli[["wmat"]], pct=pct)
  }
 # merge in other variables
  if (nrow(mfli$grid) == 0) {
    effects <- as.data.frame(to_insert)
  } else {
    effects <- data.frame(to_insert,mfli$grid)
  }
  rownames(effects) <- c()
  return(list("dame" = effects, "execute_time" = Sys.time() - start_time))
}
