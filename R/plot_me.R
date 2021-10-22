#' @title Vizualizing Marginal Effects
#' @description
#' \code{plot_me} combines a heatmap of the marginal effects of a variable with the elements representing
#' its joint distribution with another independent variable. The heatmaps is overlaid with markers whose size
#' represents the number of observations. On the sides, we append histograms of the independent variables
#' under consideration.
#' @param model Fitted model object. The package works best with GLM objects and will extract the formula, dataset, family, coefficients, and
#' the QR components of the design matrix if arguments \code{formula}, \code{data}, \code{link}, \code{coefficients}, and/or
#' \code{variance} are not explicitly specified.
#' @param data The dataframe used in estimation.
#' @param formula The formula used in estimation (if not specified, it is extracted from the fitted model object).
#' @param link The link function used in estimation (if not specified, it is extracted from the fitted model object).
#' @param coefficients The named vector of coefficients produced during the estimation (if not specified, it is extracted from the fitted model object).
#' @param variance The variance-covariance matrix to be used for computing standard errors (if not specified, it is extracted from the fitted model object).
#' @param x The independent variable whose effect is being assessed. This variable is mapped along the vertical axis.
#' @param over The second independent variable. This variable is mapped along the horizontal axis.
#' @param at A named list of the values for other indepedent variables. Unless listed, other continuous variables
#' are fixed at their means, and all factor variables are set to their modes.
#' @param mc A logical variable. If TRUE, the function will compute standard errors and sampling quantiles using Monte-Carlo simulations. If FALSE, the function will use the delta method.
#' @param iter Number of interations used in Monte-Carlo simulations. Default = 1,000.
#' @param heatmap_dim Numeric vector containin the number of rows and columns used for drawing the heatmap. Default = 100 each.
#' @param p Singificance level for the marginal effects. Default = 0.05.
#' @param gradient String vector representing color gradient for the contour plot. First element
#' is the color of the lowest value, and second element is the color of the highest value.
#' @param breaks_hor Numeric vector of breaks for the horizontal axis ("over" variable) \code{waiver()} by default.
#' @param breaks_ver Numeric vector of breaks for the vertical axis ("x" variable)  \code{waiver()} by default.
#' @param breaks_hb Numeric vector specifying probability axis for the bottom histogram.
#' \code{waiver()} by default.
#' @param breaks_hl Numeric vector specifying probability axis for the left histogram.
#' \code{waiver()} by default.
#' @param theme_hl The theme for the left panel. Function \code{theme} from \code{ggplot2} package.
#' @param theme_hb The theme for the bottom panel. Function \code{theme} from \code{ggplot2} package.
#' @param theme_mp The theme for the main panel. Function \code{theme} from \code{ggplot2} package.
#' @author \code{plot_me} visualizes ME procedure described in Moral, Sedashov, and Zhirnov (2017)
#' using the tools from \code{ggplot2} package.
#' @references
#' Moral, Mert, Evgeny Sedashov, and Andrei Zhirnov. (2017) ``Taking Distributions Seriously
#' Interpreting the Effects of Constitutive Variables in Nonlinear Models with Interactions.'' Working paper.
#' @details \code{plot_me} provides a convenient way to interpret two-way interactions using heatmaps.
#' It returns a ggplot object, which allows users to customize the plot using functions and layers from the \code{ggplot2} package.
#' @examples
#' ##Poisson regression with 2 variables and interaction between them
#' \donttest{
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
#' y <- glm(y ~ x1 + x2 + x1*x2, data = data, family = "poisson")
#' ## A heatmap with 4 areas (the default)
#' plot_me(model = y, data = data, x = "x1", over = "x2")
#' ## A heatmap with smooth transition of colors
#' plot_me(model = y, data = data, x = "x1", over = "x2", smooth=TRUE)
#' ## A heatmap with histograms at the edges
#' library(ggExtra)
#' gt <- g + theme(legend.position="left")
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x)
#' ## if more control over the histograms needed:
#' nbins <- sapply(data[c("x","z")], grDevices::nclass.FD)
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x, xparams=list(bins=nbins['z']), yparams=list(bins=nbins['x']))
#' }
#' @export

plot_me <- function(x, over, model = NULL, data = NULL,
                    link = NULL, formula = NULL, coefficients = NULL, variance = NULL,
                    at = NULL, mc = FALSE, iter = 1000,
                    heatmap_dim = c(100,100), nsteps=4, smooth=FALSE, gradient=c("#f2e6e7", "#f70429"),
                    p = 0.05, weights = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",call. = FALSE)
  }

  args <- as.list(match.call())
  obj <- lapply(args[intersect(names(formals(me)), names(args))], eval)

  if (is.null(formula)) formula <- model[["formula"]]
  formula[[2L]] <- NULL
  check.required("formula","formula")
  allvars <- all.vars(formula)

  check.required("x","character")
  check.required("over","character")

  if (is.null(data)) data <- eval(model)[["data"]]
  check.required("data","data.frame")

  outside.formula <- setdiff(c(x,over,names(at)),allvars)
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,collapse="\n"), call. = FALSE)
  outside.data <- setdiff(c(x,over),names(data))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,collapse="\n"), call. = FALSE)

  if (length(weights) != nrow(data)) weights <- rep(1,nrow(data))

  tomeans <- setdiff(allvars, c(x,over,names(at)))
  names(tomeans) <- tomeans

  obj[["at"]] <- as.list(at)
  if (length(tomeans)>0) obj[["at"]] <- c(obj[["at"]], lapply(tomeans, find.central, data=data, weights=weights))

  obj[["pct"]] <- 100*c(p/2, (1-p/2))


# data for heatmaps
  grid.li <- list(
    x = seq(from = min(data[[x]], na.rm=TRUE), to = max(data[[x]], na.rm=TRUE), length.out = heatmap_dim[1]),
    over = seq(from = min(data[[over]], na.rm=TRUE), to = max(data[[over]], na.rm=TRUE), length.out = heatmap_dim[2])
    )
  obj[["data"]] <- grid <- expand.grid(grid.li)
  colnames(obj[["data"]]) <- c(x,over)
  plotdata.hm <- do.call("me",obj)
  plotdata.hm <- merge(grid, plotdata.hm, by.x=c("x","over"), by.y=c(x,over))
  plotdata.hm[["sig"]] <- factor(rowSums(plotdata.hm[paste0("p", obj[["pct"]])]>0) %% 2, levels=c(0,1), labels=paste0(c("p<","p>"),p))

# data for scatterplots
  temp <- aggregate(list("nobs" = weights), by = list(x=data[[x]], over=data[[over]]), FUN = sum, na.action=NULL, na.rm=TRUE)
  for (j in c("x","over")) {
    bm <- data.frame(y.n=grid.li[[j]][seq_len(length(grid.li[[j]])-1)],
                     y.x=grid.li[[j]][seq_len(length(grid.li[[j]])-1)+1])
    bm[1,"y.n"] <- -Inf
    bm[nrow(bm),"y.x"] <- Inf
    temp <- merge(temp, bm, by=NULL)
    temp <- temp[temp[[j]]>temp$y.n & temp[[j]]<=temp$y.x,]
    d.n <- abs(temp$y.n-temp[[j]])
    d.x <- abs(temp$y.x-temp[[j]])
    temp[[j]] <- temp$y.x
    temp[[j]][which(d.x > d.n)] <- temp$y.n[which(d.x > d.n)]
    temp <- temp[,c("nobs","x","over"), drop=FALSE]
  }
  temp <- aggregate(nobs ~ x + over, data=temp, FUN = sum, na.action=NULL, na.rm=TRUE)
  plotdata <- merge(plotdata.hm, temp, by=c("x","over"), all=TRUE)

# plot
  nsteps <- round(nsteps,0)
  if (abs(min(plotdata$est)) > max(plotdata$est)) gradient <- gradient[c(2,1)]
  if (!smooth && nsteps >1) {
    shading <- ggplot2::scale_fill_steps(low = gradient[1], high = gradient[2], n.breaks=nsteps[1])
  } else {
    shading <- ggplot2::scale_fill_gradient(low = gradient[1], high = gradient[2])
  }
  ggplot2::ggplot(data = plotdata, aes_string(x = "over", y = "x")) +
    ggplot2::geom_raster(aes_string(fill = "est"), interpolate=FALSE) +
    ggplot2::geom_point(aes_string(size = "nobs", shape = "sig"), color = "black", data=subset(plotdata, !is.na(nobs))) +
    ggplot2::scale_shape_manual(values=c(16L,1L), drop=FALSE) +
    shading +
    ggplot2::guides(fill = guide_colourbar(order = 1L), shape = guide_legend(order = 2L), size = "none") +
    ggplot2::labs(fill="Effect Size", shape=element_blank(), x=over, y=x) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = element_blank(),
           panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"),
           aspect.ratio = 1, legend.position = "right")
}

