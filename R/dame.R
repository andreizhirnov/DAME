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
#' the corresponding bin IDs (by default, the mean value of the conditioning variable within the bin), and the targeted combinations of
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
                 nbins = 10, bin_id = NULL, use_distinct_values = TRUE,
                 discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
                 pct = c(lb=2.5, ub=97.5), iter = 1000, weights = NULL) {
## extract arguments
  link_id <- check_link(link = link, model = model)
  data <- check_data(data=data, model=model)
  f <- check_formula(formula=formula, model=model)
  bnames <- stats::model.matrix(f, data[0L,]) |> colnames()
  at <- check_at(at=at, data=data)
  weights <- check_weights(weights=weights, data=data)
  bins <- make_bins(bin_id=bin_id, over=over, data=data, 
                    use_distinct_values=use_distinct_values, nbins=nbins)
  coefficients <- check_coefs(bnames, coefficients, model)
  vcov <- check_vcov(bnames, vcov, model)
  probs <- make_bounds(pct)
  
## data pieces
  obj <- makeframes(data=data, f=f, bins=bins, at=at, weights=weights)

## model matrix
  mmat <- make_mmat(f, obj[["samples"]])
  offset <- make_offset(f, obj[["samples"]])

## adjust for an offset
  if (!is.null(offset)) {
    mmat <- cbind(offset, mmat)
    coefficients <- c(1, coefficients)
    vcov <- rbind(0, cbind(0,vcov))
  }
  
## second matrix
  if (discrete) {
    ## model matrix with a shift
    mmat_p <- make_mmat_p(f, obj[["samples"]], x=x, discrete_step = discrete_step)
    if (!is.null(offset)) mmat_p <- cbind(offset, mmat_p)    
    } else { 
    ## a matrix with cross partial derivatives of the linear prediction
    xpdm <- make_d2mdxdb(f, obj[["samples"]], x=x)
    if (!is.null(offset)) xpdm <- cbind(0, xpdm)
  }

## calculations   
  if (discrete && mc) { 
    effects <- get_ddx_mc(coefficients, vcov,
                        mmat, mmat_p,
                        obj[["wei_locs"]], obj[["wei_vals"]],
                        probs, link_id, iter)
      
    } else if (discrete) {
      effects <- get_ddx_delta(coefficients, vcov,
                           mmat, mmat_p,
                           obj[["wei_locs"]], obj[["wei_vals"]],
                           probs, link_id)
    } else if (mc) {
      effects <- get_dydx_mc(coefficients, vcov,
                  mmat, xpdm,
                  obj[["wei_locs"]], obj[["wei_vals"]],
                  probs, link_id, iter)
    } else {
      effects <- get_dydx_delta(coefficients, vcov,
                            mmat, xpdm,
                            obj[["wei_locs"]], obj[["wei_vals"]],
                            probs, link_id)
    }
  colnames(effects) <- c("est","se", names(probs))
  if (nrow(obj[["grid"]]) > 0) effects <- data.frame(effects, obj[["grid"]])
  rownames(effects) <- c() 
  return(effects)
}
