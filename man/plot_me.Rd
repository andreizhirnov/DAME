% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_me.R
\name{plot_me}
\alias{plot_me}
\title{Visualizing Marginal Effects}
\usage{
plot_me(
  x,
  over,
  model = NULL,
  data = NULL,
  link = NULL,
  formula = NULL,
  coefficients = NULL,
  vcov = NULL,
  discrete = FALSE,
  discrete_step = 1,
  at = NULL,
  mc = FALSE,
  iter = 1000,
  heatmap_dim = c(100, 100),
  p = 0.05,
  weights = NULL
)
}
\arguments{
\item{x}{a character string representing the name of the main variable of interest. Marginal effects will be computed for this variable.}

\item{over}{a character string representing the name of the conditionning variable. DAME will be computed for the bins long the range of this variable.}

\item{model}{fitted model object. The package works best with GLM objects and will extract the formula, dataset, family, coefficients, and
the QR components of the design matrix if arguments \code{formula}, \code{data}, \code{link}, \code{coefficients}, and/or
\code{vcov} are not explicitly specified.}

\item{data}{the dataset to be used to compute marginal effects (if not specified, it is extracted from the fitted model object).}

\item{link}{the name of the link function used in estimation (if not specified, it is extracted from the fitted model object).}

\item{formula}{the formula used in estimation (if not specified, it is extracted from the fitted model object).}

\item{coefficients}{the named vector of coefficients produced during the estimation (if not specified, it is extracted from the fitted model object).}

\item{vcov}{the variance-covariance matrix to be used for computing standard errors (if not specified, it is extracted from the fitted model object).}

\item{discrete}{logical. If TRUE, the function will compute the effect of a discrete change in \code{x}. If FALSE, the function will compute the partial derivative of \code{x}.}

\item{discrete_step}{The size of a discrete change in \code{x} used in computations (used only if \code{discrete=TRUE}).}

\item{at}{an optional named list of values of independent variables. These variables will be set to these value before computations.
The remaining numeric variables (except \code{x} and \code{over}) will be set to their means. The remaining factor variables will be set
to their modes.}

\item{mc}{logical. If TRUE, the standard errors and confidence intervals will be computed using simulations.
If FALSE (default), the delta method will be used.}

\item{iter}{the number of interations used in Monte-Carlo simulations. Default = 1,000.}

\item{heatmap_dim}{a numeric vector containing the number of rows and columns used for drawing the heatmap. Default = 100 each.}

\item{p}{the singificance level for the marginal effects. Default = 0.05.}

\item{weights}{an optional vector of sampling weights.}
}
\description{
\code{plot_me} produces a heatmap of the marginal effects of variable \code{x} plotted against the combinations of \code{x} and \code{over} and
adds a scatterplot representing the joint distribution these two variables in the given sample. The size of the markers represents the number of observations.
}
\details{
\code{plot_me} provides a convenient way to interpret two-way interactions using heatmaps.
It returns a ggplot object, which allows users to customize the plot using functions and layers from the \code{ggplot2} package.
}
\examples{
##Poisson regression with 2 variables and interaction between them
\dontrun{
data <- data.frame(y = rpois(10000, 10), x2 = rpois(10000, 5), x1 = rpois(10000, 3))
y <- glm(y ~ x1 + x2 + x1*x2, data = data, family = "poisson")
## A contour-plot with 4 areas
library(ggplot2)
plot_me(model = y, data = data, x = "x1", over = "x2") +
    scale_fill_steps(low="yellow", high="red", n.breaks=4)
## A heatmap with smooth transition of colors
plot_me(model = y, data = data, x = "x1", over = "x2") +
    scale_fill_gradient(low="yellow", high="red")
## A heatmap with histograms at the edges
library(ggExtra)
g <- plot_me(model = y, data = data, x = "x1", over = "x2")
gt <- g + theme(legend.position="left") + scale_fill_gradient(low="yellow", high="red")
ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x)
## if more control over the histograms needed:
nbins <- sapply(data[c("x1","x2")], grDevices::nclass.FD)
ggExtra::ggMarginal(gt, type="histogram", data=data, x=z, y=x,
   xparams=list(bins=nbins['x2']), yparams=list(bins=nbins['x1']))
}
\dontrun{
## logit
m <- glm(any_dispute ~ flows.ln*polity2 + gdp_pc, data=strikes, family="binomial")
summary(m)
plot_me(model = m, x = "flows.ln", over = "polity2") +
    scale_fill_gradient(low="yellow", high="red") +
    labs(x="Polity", y="ln(FDI flows)")
}
}
\references{
Zhirnov, Andrei, Mert Moral, and Evgeny Sedashov (2021). ``Taking Distributions Seriously: On the Interpretation
of the Estimates of Interactive Nonlinear Models.'' Working paper.
}
\author{
\code{plot_me} visualizes ME procedure described in Zhirnov, Moral, and Sedashov (2021) using the tools from \code{ggplot2} package.
}
