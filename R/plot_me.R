#' @title Vizualizing Marginal Effects
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
#' @param at an optional named list of values of independent variables. These variables will be set to these value before computations.
#' The remaining numeric variables (except \code{x} and \code{over}) will be set to their means. The remaining factor variables will be set
#' to their modes.
#' @param mc logical. If TRUE, the standard errors and confidence intervals will be computed using simulations.
#' If FALSE (default), the delta method will be used.
#' @param iter the number of interations used in Monte-Carlo simulations. Default = 1,000.
#' @param weights an optional vector of sampling weights.
#' @param heatmap_dim a numeric vector containing the number of rows and columns used for drawing the heatmap. Default = 100 each.
#' @param smooth logical. If FALSE, the values of marginal effects are broken down into bins. The number of bins is specified using \code{nlevels}.
#' @param nlevels the number of bins used for displaying marginal effects on the heatmap. Default = 4.
#' @param gradient a string vector representing the color gradient for the heatamp. The first element
#' is the color assigned to the lowest value of the marginal effect if the marginal effects are predominantly positive, and second element is the color of the
#' assigned to its highest value. If marginal effects are predominantly negative, the color assignments are swapped.
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
#' ## A heatmap with 4 areas (the default)
#' library(ggplot2)
#' plot_me(model = y, data = data, x = "x1", over = "x2")
#' ## A heatmap with smooth transition of colors
#' plot_me(model = y, data = data, x = "x1", over = "x2", smooth=TRUE)
#' ## A heatmap with histograms at the edges
#' library(ggExtra)
#' gt <- g + theme(legend.position="left")
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x)
#' ## if more control over the histograms needed:
#' nbins <- sapply(data[c("x","z")], grDevices::nclass.FD)
#' ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x,
#'    xparams=list(bins=nbins['z']), yparams=list(bins=nbins['x']))
#' }
#' \dontrun{
#' ## logit
#' m <- glm(any_dispute ~ flows.ln*polity2 + gdp_pc, data=strikes, family="binomial")
#' summary(m)
#' dame(model=m, x="flows.ln", over="polity2")
#'}
#' @export

plot_me <- function(x, over, model = NULL, data = NULL,
                    link = NULL, formula = NULL, coefficients = NULL, vcov = NULL,
                    at = NULL, mc = FALSE, iter = 1000,
                    heatmap_dim = c(100,100), smooth=FALSE, nlevels=4, gradient=c("#f2e6e7", "#f70429"),
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
  if (length(at)>0) {
    for (v in names(obj[["at"]])) {
      if (is.character(obj[["at"]][[v]]) & !is.factor(obj[["at"]][[v]])) {
        xle <-  model[["xlevels"]][[v]]
        if (is.null(xle)) xle <- sort(unique(data[[v]]))
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
  if (length(tomeans)>0) obj[["at"]] <- c(obj[["at"]], lapply(tomeans, find.central, data=data, weights=weights))

  obj[["pct"]] <- 100*c(p/2, (1-p/2))
  names(obj[["pct"]]) <- c("lb","ub")

# data for heatmaps
  grid.li <- list(
    x = seq(from = min(data[[x]], na.rm=TRUE), to = max(data[[x]], na.rm=TRUE), length.out = heatmap_dim[1]),
    over = seq(from = min(data[[over]], na.rm=TRUE), to = max(data[[over]], na.rm=TRUE), length.out = heatmap_dim[2])
    )
  obj[["data"]] <- grid <- expand.grid(grid.li)
  colnames(obj[["data"]]) <- c(x,over)
  plotdata.hm <- do.call("me",obj)
  plotdata.hm <- merge(grid, plotdata.hm, by.x=c("x","over"), by.y=c(x,over))
  plotdata.hm[["sig"]] <- factor(rowSums(plotdata.hm[c("lb","ub")]>0) %% 2, levels=c(0,1), labels=paste0(c("p<","p>"),p))

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
  nlevels <- round(nlevels,0)
  if (abs(min(plotdata$est, na.rm=TRUE)) > max(plotdata$est, na.rm=TRUE)) gradient <- gradient[c(2,1)]
  if (!smooth && nlevels >1) {
    shading <- ggplot2::scale_fill_steps(low = gradient[1], high = gradient[2], n.breaks=nlevels[1])
  } else {
    shading <- ggplot2::scale_fill_gradient(low = gradient[1], high = gradient[2])
  }
  ggplot2::ggplot(data = plotdata, ggplot2::aes_string(x = "over", y = "x")) +
    ggplot2::geom_raster(ggplot2::aes_string(fill = "est"), interpolate=FALSE) +
    ggplot2::geom_point(ggplot2::aes_string(size = "nobs", shape = "sig"), color = "black", data=plotdata[which(!is.na(plotdata$nobs)),]) +
    ggplot2::scale_shape_manual(values=c(16L,1L), drop=FALSE) +
    shading +
    ggplot2::guides(fill = ggplot2::guide_colourbar(order = 1L), shape = ggplot2::guide_legend(order = 2L), size = "none") +
    ggplot2::labs(fill="Effect Size", shape=ggplot2::element_blank(), x=over, y=x) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
           panel.grid.major = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black"),
           aspect.ratio = 1, legend.position = "right")
}

