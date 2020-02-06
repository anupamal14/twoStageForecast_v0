#' @title  Find the accuracy of a forecast
#'
#' @description Calculates the Mean Square Error (MSE), Mean Absolute
#'     Deviation (MAD) and Mean Absolute Percentage Error (MAPE) of
#'     a forecast.
#'
#' @param y the actual data.
#'
#' @param yhat the forecasted vector.
#'
#'
#' @examples
#'


twoStage.accuracy <- function(y, yhat) {
  err <- y - yhat
  # MSE
  MSE <- mean(err^2)
  # MAD
  MAD <- mean(abs(err))
  # MAPE
  MAPE <- mean(abs(err/y)) * 100

  ErrVec <- cbind(MSE, MAD, paste(toString(MAPE),"%"))
  colnames(ErrVec) <- c("MSE","MAD","MAPE")
  return(ErrVec)
}
