#' @title Marginal effects at unique combinations of independent variables
#' @description
#' \code{me} computes the marginal effects of variable \code{x} for the values of variables \code{x} and \code{over}.
#' @param x a character string representing the name of the main variable of interest. Marginal effects will be computed for this variable.
#' @param over a character string with the name of the conditioning variable. DAME will be computed for the bins long the range of this variable.
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
#' @return \code{me} returns a data frame with the estimates of the marginal effects for each combination of the variables specified in \code{x} and \code{over},
#' along with standard errors, confidence intervals, and the used values of the independent variables. All quantitative variable not included in
#' \code{at}, \code{x} and \code{over} are set to their means, and all qualitative variables (except those listed in \code{at}, \code{x} and \code{over}) are
#' converted to factors and set to their modes.
#' @examples
#' ##poisson regression with 2 variables and interaction between them
#' #fit the regression first
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
#' y <- glm(y ~ x1 + x2 + x1*x2, data = data, family = "poisson")
#' me(model = y, x = "x1", over = "x2")
#' @export

me <- function(x, over = NULL, model = NULL, data = NULL, formula = NULL, link = NULL,
               coefficients = NULL, vcov = NULL,
               discrete = FALSE, discrete_step = 1, at = NULL, mc = FALSE,
               pct = c(lb=2.5, ub=97.5), iter = 1000, weights = NULL) {

  ## extract arguments
  link_id <- check_link(link = link, model = model)
  data <- check_data(data=data, model=model)
  f <- check_formula(formula=formula, model=model)
  bnames <- stats::model.matrix(f, data[0L,]) |> colnames()
  weights <- check_weights(weights=weights, data=data)
  at <- check_at(at=at, data=data)
  #### check the offset
  tt <- stats::terms(f)
  mt <- attr(tt, "variables")
  if (!is.null(attr(tt, "offset"))) {
    mt[[1L+attr(tt, "offset")]] <- NULL
  }
  ### add the central stats for the remaining variables
  tomeans <- setdiff(all.vars(mt), c(names(at), x, over))
  for (v in tomeans) {
    at[[v]] <- find_central(x=v, data=data, weights=weights)
  }
  bins <- make_bins(bin_id=1:nrow(data), data = data)
  
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
  if (nrow(obj[["grid"]]) > 0) effects <- data.frame(effects, obj[["samples"]])
  rownames(effects) <- c() 
  return(effects)
}

