#' @title Visualizing Marginal Effects
#' @description
#' \code{plot_me} produces a heatmap of the marginal effects of variable \code{x} plotted against the combinations of \code{x} and \code{over} and
#' adds a scatterplot representing the joint distribution these two variables in the given sample. The size of the markers represents the number of observations.
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
#' @param discrete logical. If TRUE, the function will compute the effect of a discrete change in \code{x}. If FALSE, the function will compute the partial derivative of \code{x}.
#' @param discrete_step The size of a discrete change in \code{x} used in computations (used only if \code{discrete=TRUE}).
#' @param at an optional named list of values of independent variables. These variables will be set to these value before computations.
#' The remaining numeric variables (except \code{x} and \code{over}) will be set to their means. The remaining factor variables will be set
#' to their modes.
#' @param mc logical. If TRUE, the standard errors and confidence intervals will be computed using simulations.
#' If FALSE (default), the delta method will be used.
#' @param iter the number of interations used in Monte-Carlo simulations. Default = 1,000.
#' @param weights an optional vector of sampling weights.
#' @param heatmap_dim a numeric vector containing the number of rows and columns used for drawing the heatmap. Default = 100 each.
#' @param p the singificance level for the marginal effects. Default = 0.05.
#' @author \code{plot_me} visualizes ME procedure described in Zhirnov, Moral, and Sedashov (2021) using the tools from \code{ggplot2} package.
#' @references
#' Zhirnov, Andrei, Mert Moral, and Evgeny Sedashov (2021). ``Taking Distributions Seriously: On the Interpretation
#' of the Estimates of Interactive Nonlinear Models.'' Working paper.
#' @details \code{plot_me} provides a convenient way to interpret two-way interactions using heatmaps.
#' It returns a ggplot object, which allows users to customize the plot using functions and layers from the \code{ggplot2} package.
#' @examples
#' ##Poisson regression with 2 variables and interaction between them
#' \dontrun{
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
#' y <- glm(y ~ x1 + x2 + x1*x2, data = data, family = "poisson")
#' ## A contour-plot with 4 areas
#' library(ggplot2)
#' plot_me(model = y, data = data, x = "x1", over = "x2") +
#'     scale_fill_steps(low="yellow", high="red", n.breaks=4)
#' ## A heatmap with smooth transition of colors
#' plot_me(model = y, data = data, x = "x1", over = "x2") +
#'     scale_fill_gradient(low="yellow", high="red")
#' ## A heatmap with histograms at the edges
#' library(ggExtra)
#' g <- plot_me(model = y, data = data, x = "x1", over = "x2")
#' gt <- g + theme(legend.position="left") + scale_fill_gradient(low="yellow", high="red")
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x)
#' ## if more control over the histograms needed:
#' nbins <- sapply(data[c("x1","x2")], grDevices::nclass.FD)
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x,
#'    xparams=list(bins=nbins['x2']), yparams=list(bins=nbins['x1']))
#' }
#' \dontrun{
#' ## logit
#' m <- glm(any_dispute ~ flows.ln*polity2 + gdp_pc, data=strikes, family="binomial")
#' summary(m)
#' plot_me(model = m, x = "flows.ln", over = "polity2") +
#'     scale_fill_gradient(low="yellow", high="red") +
#'     labs(x="Polity", y="ln(FDI flows)")
#'}
#' @export

plot_me <- function(x, over, model = NULL, data = NULL,
                    link = NULL, formula = NULL, coefficients = NULL, vcov = NULL,
                    discrete = FALSE, discrete_step = 1,
                    at = NULL, mc = FALSE, iter = 1000,
                    heatmap_dim = c(100,100),
                    p = 0.05, weights = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",call. = FALSE)
  }
 
## extract arguments
  link_id <- check_link(link = link, model = model) 
  data <- check_data(data=data, model=model)
  f <- check_formula(formula=formula, model=model)
  bnames <- stats::model.matrix(f, data[0L,]) |> colnames()
  weights <- check_weights(weights=weights, data=data)
  at <- check_at(at=at, data=data)
  if (any(sapply(at, length)>1)) {
    at <- lapply(at, `[`, 1L)
  }
  ### add the central stats for the remaining variables
  tomeans <- setdiff(all.vars(f), c(names(at), x, over))
  for (v in tomeans) {
    at[[v]] <- find_central(x=v, data=data, weights=weights)
  }
  probs <- setNames(100*c(p/2, (1-p/2)), c("lb","ub")) |> make_bounds()
  coefficients <- check_coefs(bnames, coefficients, model)
  vcov <- check_vcov(bnames, vcov, model) 
  
  ### make a compressed dataset by grid and send a copy to obj
  grid.li <- list(
    x = seq(from = min(data[[x]], na.rm=TRUE), to = max(data[[x]], na.rm=TRUE), length.out = heatmap_dim[1L]),
    over = seq(from = min(data[[over]], na.rm=TRUE), to = max(data[[over]], na.rm=TRUE), length.out = heatmap_dim[2L])
  )
  gdata <- expand.grid(grid.li)
  colnames(gdata) <- c(x,over)
  bins <- make_bins(bin_id=1:nrow(gdata), data = gdata)

  ### calculate
  obj <- makeframes(data=gdata, f=f, bins=bins, at=at, weights=NULL)
  
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
    plotdata <- get_ddx_mc(coefficients, vcov,
                          mmat, mmat_p,
                          obj[["wei_locs"]], obj[["wei_vals"]],
                          probs, link_id, iter)
    
  } else if (discrete) {
    plotdata <- get_ddx_delta(coefficients, vcov,
                             mmat, mmat_p,
                             obj[["wei_locs"]], obj[["wei_vals"]],
                             probs, link_id)
  } else if (mc) {
    plotdata <- get_dydx_mc(coefficients, vcov,
                           mmat, xpdm,
                           obj[["wei_locs"]], obj[["wei_vals"]],
                           probs, link_id, iter)
  } else {
    plotdata <- get_dydx_delta(coefficients, vcov,
                              mmat, xpdm,
                              obj[["wei_locs"]], obj[["wei_vals"]],
                              probs, link_id)
  }
  
  colnames(plotdata) <- c("est","se", "lb", "ub")
  plotdata <- data.frame(plotdata, gdata)
  
  ### find the number of observations by grid bin
  if (is.null(weights)) weights <- 1.0
  plotdata$nobs <- count_nearest(as.matrix(data[,c(x, over),drop=FALSE]),
                        as.matrix(plotdata[,c(x, over),drop=FALSE]),
                        weights) |> as.vector()
# data for heatmaps
  plotdata[["sig"]] <- factor(rowSums(plotdata[c("lb","ub")]>0) %% 2, 
                              levels=c(0,1), 
                              labels=paste0(c("p<","p>"),p))
  
# plot
  ggplot2::ggplot(data = plotdata, ggplot2::aes(x = .data[[over]], y = .data[[x]])) +
    ggplot2::geom_raster(ggplot2::aes(fill = est), interpolate=FALSE) +
    ggplot2::geom_point(ggplot2::aes(size = nobs, shape = sig), 
                        color = "black", 
                        data=plotdata[which(plotdata$nobs>0),]) +
    ggplot2::scale_shape_manual(values=c(1L,16L), drop=FALSE) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(order = 1L), 
                    shape = ggplot2::guide_legend(order = 2L), 
                    size = "none") +
    ggplot2::labs(fill="Effect Size", 
                  shape=NULL, 
                  x=over, y=x) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
           panel.grid.major = ggplot2::element_blank(), 
           panel.border = ggplot2::element_rect(colour = "black"),
           aspect.ratio = 1, legend.position = "right")
}

