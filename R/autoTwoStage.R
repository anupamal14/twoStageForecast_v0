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
  nrLevels <- length(seas_periods)
  len <- length(d)
  if (is.null(n) & is.null(p)) {
    p <- seas_periods[nrLevels]
    n <- len - p
  }

  if (nrLevels == 1) {
    final <- firstStage(d, n, p, seas_periods, plotFlag = plotFlag)

    # Prepare output
    modNames <- c('TBATS', 'ARIMA')
    matchNames <- c('tbats','arima')
    psiMat <- array(length(final))
    psiMat[1] <- twoStage.penalizedLikelihood(n,
                                              as.double(final$tbats$accuracy[1]),
                                              length(final$tbats$model$parameters$vect) +
                                                length(final$tbats$model$seed.states))
    psiMat[2] <- twoStage.penalizedLikelihood(n,
                                              as.double(final$arima$accuracy[1]),
                                              length(final$arima$model$coef) + 1)
    myCall <- gsub("language", "Call: ",str(sys.call()))
    cat(paste(myCall,
                    "Seasonality level selected: Not applicable",
                    paste("First Stage method: ", modNames[which.min(psiMat)]),
                    "Second Stage method: Not applicable",
                    " ",
                    paste("psi: ",psiMat[which.min(psiMat)]),
              " ",
                    "Error measures:",
              "MSE                                    MAD                                MAPE",
                    toString(eval(parse(text=gsub(" ","",paste("final$",matchNames[which.min(psiMat)],"$accuracy"))))),
              sep="\n"))
  } else {
  meanBy <- as.vector(t(matrix(rep(c(1:ceiling(len/seas_periods[2])),seas_periods[2]),
                               nrow = ceiling(len/seas_periods[2]))))
  d2 <- aggregate(d,by = list(meanBy[1:len]),mean)
  ylow <- d2$x

  # First Stage
  stage1 <- firstStage(ylow, n/seas_periods[2], p/seas_periods[2], seas_periods[1], regMat, plotFlag)
  # Second Stage
  stage2 <- secondStage(d, seas_periods)
  # Combine the output of the two stages
  }

}
