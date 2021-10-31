## "DAME" package: Tools for Data-Conscious Interpretation of Nonlinear Interactive Models in R

To install:
```
library(devtools)
install_github("andreizhirnov/DAME")
```

To load:
```
library(DAME)
```

Exported functions:

*dame()* computes the distribution-weighted average marginal effects (DAME) as described in Zhirnov, Moral, and Sedashov (2021).

*me()* computes the marginal effects of a variable at the unique combinations of the values of listed variables.

*mem()* computes the marginal effects of a variable at the specified values of select independent variables and the means/modes of all the other variables.

*ame()* computes the average marginal effects of a variable.

*plot_me()* produces a heatmap of the marginal effects of a variable plotted against the combinations of the values of two variables and adds a scatterplot representing the joint distribution these two variables in the given sample.

Please Cite:  Zhirnov, Andrei, Mert Moral, and Evgeny Sedashov. 2021. "Taking Distributions Seriously: On the Interpretation of the Estimates of Interactive Nonlinear Models." Working paper.
