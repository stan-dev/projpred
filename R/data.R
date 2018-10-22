#' Binomial toy example.
#'
#' @format A simulated classification dataset containing 100 observations.
#' \describe{
#'   \item{y}{target, 0 or 1.}
#'   \item{x}{features, 30 in total.}
#' }
#' @source \url{http://web.stanford.edu/~hastie/glmnet/glmnetData/BNExample.RData}
"df_binom"

#' Gaussian toy example.
#'
#' @format A simulated regression dataset containing 100 observations.
#' \describe{
#'   \item{y}{target, real-valued.}
#'   \item{x}{features, 20 in total. Mean and sd approximately 0 and 1.}
#' }
#' @source \url{http://web.stanford.edu/~hastie/glmnet/glmnetData/QSExample.RData}
"df_gaussian"

#' Mesquite data set.
#' 
#' The mesquite bushes yields data set from Gelman
#' and Hill (2007) (\url{http://www.stat.columbia.edu/~gelman/arm/}).
#'
#' @format  The outcome variable is the total weight (in grams) of photosynthetic
#' material as derived from actual harvesting of the bush. The predictor
#' variables are:
#' \describe{
#' \item{diam1}{diameter of the canopy (the leafy area of the bush)
#' in meters, measured along the longer axis of the bush.}
#' \item{diam2}{canopy diameter measured along the shorter axis}
#' \item{canopy height}{height of the canopy.}
#' \item{total height}{total height of the bush.}
#' \item{density}{plant unit density (# of primary stems per plant unit).}
#' \item{group}{group of measurements (0 for the first group, 1 for the second group)}
#' }
#'
#' @source \url{http://www.stat.columbia.edu/~gelman/arm/examples/}
"mesquite"

