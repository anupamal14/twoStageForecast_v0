#' @title  Classical decomposition method for seasonality estimation
#'
#' @description Estimates the seasonality at different levels by
#'     subtracting the mean for additive seasonality
#'     (dividing for multiplicative) and then averaging across
#'     the entire data.
#'
#' @param d the actual data.
#'
#' @param yhat the forecasted vector.
#'
#'
#' @examples
#'

twoStage.classical <- function(d, period) {
  n <- length(d)
  agg_vec <- rep(1:period,ceiling(n/period))
  agg <- aggregate(d,list(vec = agg_vec[1:n]),mean)
  return(agg[,2])

}
