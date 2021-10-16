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
#' @param geom_hl Specifies ggplot layer of the type \code{geom_histogram} for the left histogram.
#' If NULL, default settings are used.
#' @param geom_hb Specifies ggplot layer of the type \code{geom_histogram} for the bottom histogram.
#' If NULL, default settings are used.
#' @param theme_hl The theme for the left panel. Function \code{theme} from \code{ggplot2} package.
#' @param theme_hb The theme for the bottom panel. Function \code{theme} from \code{ggplot2} package.
#' @param theme_mp The theme for the main panel. Function \code{theme} from \code{ggplot2} package.
#' @author \code{plot_me} visualizes ME procedure described in Moral, Sedashov, and Zhirnov (2017)
#' using the tools from \code{ggplot2} package.
#'
#' @references
#' Moral, Mert, Evgeny Sedashov, and Andrei Zhirnov. (2017) ``Taking Distributions Seriously
#' Interpreting the Effects of Constitutive Variables in Nonlinear Models with Interactions.'' Working paper.
#' @details \code{plot_me} provides a convenient way to interpret two-way interactions using
#' a heatmap with the histograms of the independent variables on x and y axis. Procedure allows users to customize the
#' plot using functions and layers from the \code{ggplot2} package. For instance, the look
#' of contour plot and histograms is fully customizable, but the overall style (histograms on the
#' left and on the bottom with contourplot in the center) is preserved for all plots.
#' @examples
#' ##Poisson regression with 2 variables and interaction between them
#' data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
#' y <- glm(y ~ x1 + x2 + x1*x2, data = data, family = "poisson")
#' plot_me(model = y, data = data, x = "x1", over = "x2")
#' @export
#' @import ggplot2
#' @import grid

plot_me <- function(x, over, model = NULL, data = NULL,
                    link = NULL, formula = NULL, coefficients = NULL, variance = NULL,
                    at = NULL, mc = FALSE, iter = 1000,
                    heatmap_dim = c(100,100),
                    p = 0.05,
                    gradient=c("#f2e6e7", "#f70429"),
                    breaks_hor=waiver(), breaks_ver=waiver(),
                    breaks_hb=waiver(), breaks_hl=waiver(),
                    geom_hl = geom_histogram(aes(y=stat(count)/sum(stat(count)))),
                    geom_hb = geom_histogram(aes(y=stat(count)/sum(stat(count)))),
                    theme_mp=theme_bw() + theme(panel.grid.minor = element_blank(),
                                                panel.grid.major = element_blank(),
                                                panel.border = element_rect(colour = "black"),
                                                axis.text.y = element_text(angle = 90),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                plot.margin = unit(c(0.5, 0.5, 0.01, 0.01), "lines"),
                                                aspect.ratio = 1,
                                                legend.position = "right"),
                    theme_hl=theme_bw() + theme(panel.border = element_rect(colour = NA),
                                                panel.grid.minor = element_blank(),
                                                panel.grid.major = element_blank(),
                                                axis.ticks.length = unit(0.1, "cm"),
                                                plot.margin = unit(c(0.5, 0, 0, 0.5), "lines"),
                                                axis.title.y = element_text(face = "bold"),
                                                axis.title.x = element_blank(),
                                                aspect.ratio = 5),
                    theme_hb=theme_bw() + theme(panel.border = element_rect(colour = NA),
                                                panel.grid.minor = element_blank(),
                                                panel.grid.major = element_blank(),
                                                axis.ticks.length = unit(0.1, "cm"),
                                                plot.margin= unit(c(0, 0.5, 0, 0), "lines"),
                                                axis.text.y = element_text(angle=90),
                                                axis.title.x = element_text(face = "bold"),
                                                axis.title.y = element_blank(),
                                                aspect.ratio = 1/5)) {
  # extract arguments from the call
  args <- as.list(match.call())
  if (!("formula" %in% names(args))) args[["formula"]] <- quote(eval(args[["model"]])[["formula"]])
  if (!("data" %in% names(args))) args[["data"]] <- quote(eval(args[["model"]])[["data"]])
  if (!("link" %in% names(args))) args[["link"]] <- quote(eval(args[["model"]])[["family"]][["link"]])
  if (!("coefficients" %in% names(args))) args[["coefficients"]] <- quote(stats::coef(eval(args[["model"]])))
  if (!("variance" %in% names(args))) args[["variance"]] <- quote(stats::vcov(eval(args[["model"]])))

  # check the required arguments and coerce the specified arguments into a proper class
  checks <- list(
    required=c("x","over","formula","data","link","coefficients","variance"),
    types = list(x = "character", over = "character",  data = "data.frame", link = "character", formula = "formula", coefficients = "numeric", variance = "matrix", at = "list",
                 mc = "logical", iter = "integer", heatmap_dim = "integer", p = "numeric", gradient = "character", breaks_hor = "numeric", breaks_ver = "numeric", breaks_hb = "numeric", breaks_hl = "numeric",
                 geom_hl = "Layer", geom_hb = "Layer", theme_mp = "theme", theme_hl = "theme", theme_hb = "theme"),
    lengths = list(x = 1L, over = 1L, link = 1L, mc = 1L, iter = 1L, heatmap_dim = 2L, p = 1L, gradient = 2L)
  )
  check.args(args=args, checks=checks)

  # check if x, over, at variables are included in the formula
  updform <- formula
  updform[[2L]] <- NULL
  allvars <- all.vars(updform)
  outside.formula <- setdiff(c(x,over,names(at)),allvars)
  if (length(outside.formula)>0) stop(paste("Failed to find the following variables in the formula:",outside.formula,sep="\n"), call. = FALSE)
  # check if x and over variables are included in the data
  outside.data <- setdiff(c(x,over),names(data))
  if (length(outside.data)>0) stop(paste("Failed to find the following variables in the dataset:",outside.data,sep="\n"), call. = FALSE)
  # misc checks
  if (!is.null(at)) {
    for (i in 1:length(at)) {
      if (length(at[[i]])!=1) stop("Error: Each 'at' variable must have one value", call. = FALSE)
    }
  }
  if (!is.null(p)) {
    if (p > 1 | p <= 0) stop("Error: 'p' must be between 0 and 1", call. = FALSE)
  }

  if (!is.null(link)) {
    if (!(link %in% c("logit","probit","cauchit","cloglog","identity","log","sqrt","1/mu^2","inverse"))) {
      stop("Invalid link name. Valid links include 'logit','probit','cauchit','cloglog','identity','log','sqrt','1/mu^2','inverse'", call. = FALSE)
    }
  }
  # preliminaries
  condvar <- setdiff(allvars,x)
  make.dydm(link=link)
  colsum <- c(at,lapply(data[intersect(allvars,setdiff(names(data),c(names(at),x,over)))], function(x) if (is.numeric(x)) return(mean(x,na.rm=TRUE)) else return(find.mode(x))))
  # prepare dataframe for heatmap
  data.hm <- list()
  data.hm[[x]] <- seq(from = min(data[[x]]), to = max(data[[x]]), length.out = heatmap_dim[2])
  data.hm[[over]] <- seq(from = min(data[[over]]), to = max(data[[over]]), length.out = heatmap_dim[1])
  data.hm <- data.frame(c(expand.grid(data.hm),colsum))
  # prepare dataframe for scatterplots
  newdata <- na.omit(data)
  data.sp <- as.data.frame(c(aggregate(list("nobs" = 1L:nrow(newdata)), by = newdata[c(x,over)], FUN = length),colsum))
  data.hb <- aggregate(list("proportion" = 1L:nrow(newdata)), by = newdata[over], FUN = length)
  data.hb[["proportion"]] <- data.hb[["proportion"]]/sum(data.hb[["proportion"]])
  data.hl <- aggregate(list("proportion" = 1L:nrow(newdata)), by = newdata[x], FUN = length)
  data.hl[["proportion"]] <- data.hl[["proportion"]]/sum(data.hl[["proportion"]])
  # compute ME
  mf.hm <- model.frame(formula = updform, data = data.hm)
  mf.sp <- model.frame(formula = updform, data = data.sp)
  mmat.hm <- model.matrix(object = updform, data = mf.hm)
  mmat.sp <- model.matrix(object = updform, data = mf.sp)
  beta.names <- intersect(names(coefficients),colnames(mmat.hm))
  if (length(beta.names) == 0) stop("Failed to associate the coefficients with the columns names in the dataset", call. = FALSE)
  coef_vect <- matrix(coefficients[beta.names],nrow=1L)
  var_covar <- variance[beta.names,beta.names]
  mmat.hm <- mmat.hm[,beta.names,drop=FALSE]
  mmat.sp <- mmat.sp[,beta.names,drop=FALSE]
  make.dmdx(formula=formula,bnames=beta.names,xvarname=x)
  if (mc) {
    coef_matrix <- MASS::mvrnorm(n = iter, mu = coef_vect, Sigma = var_covar, empirical = TRUE)
    me_matrix.hm <- dydm(mx(mmat = mmat.hm, coefficients = coef_matrix)) * dmdx(mmat = mmat.hm, data=data.hm, coefficients = coef_matrix)
    me_matrix.sp <- dydm(mx(mmat = mmat.sp, coefficients = coef_matrix)) * dmdx(mmat = mmat.sp, data=data.sp, coefficients = coef_matrix)
    estimate.hm <- colMeans(me_matrix.hm)
    plotdata.sp <- data.frame("estimate" = colMeans(me_matrix.sp), "se" = apply(me_matrix.sp, 2L, sd))
  } else {
    estimate.hm <- as.vector(dydm(mx(mmat = mmat.hm, coefficients = coef_vect)) * dmdx(mmat = mmat.hm, data=data.hm, coefficients = coef_vect))
    plotdata.sp <- data.frame("estimate" = as.vector(dydm(mx(mmat = mmat.sp, coefficients = coef_vect)) * dmdx(mmat = mmat.sp, data=data.sp, coefficients = coef_vect)))
    m_0 <- as.vector(mx(mmat = mmat.sp, coefficients = coef_vect))
    L <- d2ydm2(m_0)*as.vector(dmdx(mmat = mmat.sp, data=data.sp, coefficients = coef_vect))*dmdb(mmat = mmat.sp, data=data.sp) + dydm(m_0)*d2mdxdb(mmat=mmat.sp, data = data.sp)
    var <- as.vector(apply(L,1L, function(x) (x %*% var_covar) %*% x))
    var[var < 0] <- 0
    plotdata.sp[["se"]] <- sqrt(var)
  }
  plotdata.sp[["vertical"]] <- data.sp[[x]]
  plotdata.sp[["horizontal"]] <- data.sp[[over]]
  plotdata.sp[["size"]] <- data.sp[["nobs"]]/sum(data.sp[["nobs"]])
  plotdata.sp[["prob0"]] <- with(plotdata.sp, stats::pnorm(q = 0, mean = estimate, sd = se))
  plotdata.sp[["significance"]] <- with(plotdata.sp, ifelse((prob0 < p/2 | prob0 > 1-p/2),"Significant","Not Significant"))
  #limits to the axes
  hor_range <- max(data[[over]]) - min(data[[over]])
  ver_range <- max(data[[x]]) - min(data[[x]])
  hor_limits <- c((min(data[[over]])-0.05*hor_range),(max(data[[over]])+0.05*hor_range))
  ver_limits <- c((min(data[[x]])-0.05*ver_range),(max(data[[x]])+0.05*ver_range))
  #plotting by part
  mainpanel <- ggplot(data = data.hm, aes_string(x = over, y = x)) +
    geom_raster(aes(fill = estimate.hm), interpolate=TRUE) +
    geom_point(data = plotdata.sp,
               aes(x = horizontal, y = vertical, size = size, shape = significance),
               color = "black") +
    scale_shape_manual(values=c(1,16))+
    scale_fill_gradient(low = gradient[1], high = gradient[2])+
    scale_x_continuous(limits = hor_limits, position = "top", breaks = breaks_hor, expand = c(0,0))+
    scale_y_continuous(limits = ver_limits, position = "right", breaks = breaks_ver, expand = c(0,0))+
    guides(fill = guide_colourbar(order = 1),
           shape = guide_legend(order = 2),
           size = FALSE)+
    labs(fill="Marginal Effect", shape=paste0("Significance\n(p=",p,")"))+
    theme_mp
  hb <- ggplot(data = data.hb, aes_string(x = over, y="proportion")) +
    geom_col() +
    scale_x_continuous(limits = hor_limits, position = "bottom", breaks = breaks_hor, expand = c(0,0)) +
    scale_y_reverse(position = "right", breaks = breaks_hb) +
    theme_hb
  hl <- ggplot(data = data.hl, aes_string(x = x, y="proportion")) +
    geom_col() +
    scale_x_continuous(limits = ver_limits, position = "bottom", breaks = breaks_ver, expand = c(0,0)) +
    scale_y_reverse(position = "right", breaks = breaks_hl) +
    theme_hl +  coord_flip()
  gt1 <- ggplot_gtable(ggplot_build(mainpanel))
  gt2 <- ggplot_gtable(ggplot_build(hb))
  gt3 <- ggplot_gtable(ggplot_build(hl))
  # adjust the grids if the number of columns is different between the plots
  findnull <- function(what) {
    for (i in 1:length(what)) {
      if (is.null(attr(what[[i]],"unit"))) { next }
      else if (attr(what[[i]],"unit") == "null") {
        return(i)
        break
      }
    }
    return(0)
  }
  adj <- findnull(gt1[["widths"]]) - findnull(gt2[["widths"]])
  if (adj>0) {
    gt2 <- gtable::gtable_add_cols(x = gt2, widths = unit(rep(0,adj),"cm"), pos = 0)
  } else if (adj<0) {
    gt1 <- gtable::gtable_add_cols(x = gt1, widths = unit(rep(0,-adj),"cm"), pos = 0)
  }
  adj <- length(gt1[["widths"]]) - length(gt2[["widths"]])
  if (adj>0) {
    gt2 <- gtable::gtable_add_cols(x = gt2, widths = unit(rep(0,adj),"cm"), pos = -1)
  } else if (adj<0) {
    gt1 <- gtable::gtable_add_cols(x = gt1, widths = unit(rep(0,adj<0),"cm"), pos = -1)
  }
  gt <- rbind(gt1, gt2, size="first")
  adj <- findnull(gt[["heights"]]) - findnull(gt3[["heights"]])
  if (adj>0) {
    gt3 <- gtable::gtable_add_rows(x = gt3, heights = unit(rep(0,adj),"cm"), pos = 0)
  } else if (adj<0) {
    gt <- gtable::gtable_add_rows(x = gt, heights = unit(rep(0,-adj),"cm"), pos = 0)
  }
  adj <- length(gt[["heights"]]) - length(gt3[["heights"]])
  if (adj>0) {
    gt3 <- gtable::gtable_add_rows(x = gt3, heights = unit(rep(0,adj),"cm"), pos = -1)
  } else if (adj<0) {
    gt <- gtable::gtable_add_rows(x = gt, heights = unit(rep(0,-adj),"cm"), pos = -1)
  }
  mid1 <- length(gt[["heights"]])
  mid2 <- length(gt3[["heights"]])
  gt <- cbind(gt3, gt, size="last")
  ar_l <- ifelse(is.null(theme_hl[["aspect.ratio"]]),1/5,theme_hl[["aspect.ratio"]])
  ar_m <- ifelse(is.null(theme_mp[["aspect.ratio"]]),1,theme_mp[["aspect.ratio"]])
  ar_b <- ifelse(is.null(theme_hb[["aspect.ratio"]]),5,theme_hb[["aspect.ratio"]])
  c_widths <- c(1L,ar_l/ar_m)
  c_heights <- c(ar_m/ar_b,1L)
  # adjust the widths and heights according to the aspect ratios
  adjust_dim <- function(what,value){
    result <- what
    for (i in 1L:length(result)) {
      if (is.null(attr(result[[i]],"unit"))) { next }
      else if (attr(result[[i]],"unit") == "null") {result[[i]][[1L]] <- value}
    }
    result
  }
  range_l <- 1L:length(gt3[["widths"]])
  range_r <- (length(gt3[["widths"]])+1L):length(gt[["widths"]])
  range_t <- 1L:length(gt1[["heights"]])
  range_b <- (length(gt1[["heights"]])+1L):length(gt[["heights"]])
  gt[["widths"]][range_l] <- adjust_dim(gt[["widths"]][range_l],c_widths[1L])
  gt[["widths"]][range_r] <- adjust_dim(gt[["widths"]][range_r],c_widths[2L])
  gt[["heights"]][range_t] <- adjust_dim(gt[["heights"]][range_t],c_heights[1L])
  gt[["heights"]][range_b] <- adjust_dim(gt[["heights"]][range_b],c_heights[2L])
  grid.newpage()
  grid.draw(gt)
}
