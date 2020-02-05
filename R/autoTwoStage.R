#' @title  Two-stage methods for forecasting time series with mulitple
#'    levels of seasonality
#'
#' @description Estimates the seasonality at different levels by
#'     subtracting the mean for additive seasonality
#'     (dividing for multiplicative) and then averaging across
#'     the entire data.
#'
#' @param d the actual data.
#'
#' @param seas_periods periods associated with the seasonality levels.
#'
#'
#' @examples
#'

auto.twoStage <- function(d, n = NULL, p = NULL, seas_periods = NULL, dReg = NULL,
                          plotFlag = FALSE) {
  # number of levels
  levels <- length(seas_periods)
  len <- length(d)
  if (is.null(n) & is.null(p)) {
    p <- seas_periods[levels()]
    n <- len - p
  }

  meanBy <- as.vector(t(matrix(rep(c(1:ceiling(len/seas_periods[3])),seas_periods[3]),
                               nrow = ceiling(len/seas_periods[3]))))
  d2 <- aggregate(d,by = list(meanBy[1:len]),mean)
  ylow <- d2$x

  # First Stage
  stage1 <- firstStage(ylow, n, p, seas_periods, regMat, plotFlag)
  # Second Stage
  stage2 <- secondStage(d, seas_periods)
  # Combine the output of the two stages

}
