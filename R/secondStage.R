#' @title  Seasonality estimation for high frequency time series
#'
#' @description Estimates the seasonality using different methods like
#'     classical decomposition method based, polynomial fitting,
#'     trigonometric methods and spline.
#'
#' @param d the actual data.
#'
#' @param seas_periods periods associated with the levels of seasonality
#'
#'
#' @examples
#'


secondStage <- function(d, seas_periods) {
  s <- sort(seas_periods)
  nrSeas <- length(seas_periods)

  seas <- list()
  seas$level <- list()

  for (k in 1:nrSeas) {
    # lowest level
    d_1 <- colMeans(matrix(data, nrow = s[k]))
    # Replicate this 12 times and reshape into a vector so that
    # it is the same length as d
    d_1_rep <- as.vector(t(matrix(t(rep(t(d_1), s[k])), ncol = s[k])))


    # Get de-mean signal
    multInput <- d/d_1_rep
    addInput <- d - d_1_rep

    seas$level[k]$trigMult <- twoStage.trig(multInput)
    seas$level[k]$trigAdd <- twoStage.trig(addInput)
    seas$level[k]$splineMult <- twoStage.spline(multInput)
    seas$level[k]$splineAdd <- twoStage.spline(addInput)
    #seas$level[k]$polyMult <- twoStage.poly(multInput)
    #seas$level[k]$polyAdd <- twoStage.poly(addInput)
    seas$level[k]$classicalMult <- twoStage.classical(multInput, s[min(k+1, nrSeas)])
    seas$level[k]$classicalAdd <- twoStage.classical(addInput, s[min(k+1, nrSeas)])
  }

  return(seas)
}
