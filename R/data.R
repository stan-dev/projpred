#' Binomial toy example
#'
#' @format A simulated classification dataset containing 100 observations.
#' \describe{
#'   \item{y}{response, 0 or 1.}
#'   \item{x}{predictors, 30 in total.}
#' }
#' @source <https://web.stanford.edu/~hastie/glmnet/glmnetData/BNExample.RData>
"df_binom"

#' Gaussian toy example
#'
#' @format A simulated regression dataset containing 100 observations.
#' \describe{
#'   \item{y}{response, real-valued.}
#'   \item{x}{predictors, 20 in total. Mean and SD are approximately 0 and 1,
#'   respectively.}
#' }
#' @source <https://web.stanford.edu/~hastie/glmnet/glmnetData/QSExample.RData>
"df_gaussian"

#' Mesquite data set
#'
#' The mesquite bushes yields dataset from Gelman and Hill (2006)
#' (<http://www.stat.columbia.edu/~gelman/arm/>).
#'
#' @format The response variable is the total weight (in grams) of
#'   photosynthetic material as derived from actual harvesting of the bush. The
#'   predictor variables are:
#' \describe{
#'   \item{diam1}{diameter of the canopy (the leafy area of the bush) in meters,
#'   measured along the longer axis of the bush.}
#'   \item{diam2}{canopy diameter measured along the shorter axis.}
#'   \item{canopy height}{height of the canopy.}
#'   \item{total height}{total height of the bush.}
#'   \item{density}{plant unit density (# of primary stems per plant unit).}
#'   \item{group}{group of measurements (0 for the first group, 1 for the second
#'   group).}
#' }
#' @references Gelman, A. and Hill, J. (2006). *Data Analysis Using Regression
#'   and Multilevel/Hierarchical Models*. Cambridge University Press. DOI:
#'   [10.1017/CBO9780511790942](https://doi.org/10.1017/CBO9780511790942).
#' @source <http://www.stat.columbia.edu/~gelman/arm/examples/mesquite/mesquite.dat>
"mesquite"
