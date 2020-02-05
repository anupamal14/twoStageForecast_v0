#' @title  First stage of the two-stage method
#'
#' @description Fits time series models for the averaged (Y_low)
#'     data. It computes and returns five possible models and
#'     their preidctions.
#'
#' @param ylow averaged data.
#'
#' @param n Number of samples in the model period.
#'
#' @param p Number of samples in the holdout period.
#'
#' @param seas_periods periods corresponding to levels of
#'     seasonality in the data
#'
#' @param regMat covariates to be used for regression.
#'
#' @param plotFlag plot results. default is TRUE.
#'
#' @examples
#'
#' @importFrom forecast tbats msts
#' @importFrom stats arima predict
#'
#' @export
#'
#'

firstStage <- function(ylow, n, p, seas_periods, regMat = NULL,
                       plotFlag = FALSE){

  yStg1 <- msts(ylow, seasonal.periods = seas_periods)

  minPsi <- 1000000 # large value initialization
  for (d in 0:2)
    for (D in 0:2) {
      fit <- auto.arima(ylow[1:n], d = d, D = D, max.p = 6, max.q = 6, max.P = 6, max.Q = 6,
                        start.p = 0, start.q = 0, start.P = 0, start.Q = 0,
                        seasonal = TRUE)
      psi <- twoStage.penalizedLikelihood(length(fit$residuals),
                                          mean(fit$residuals^2), (length(fit$coef) + 1))
      if (psi < minPsi) {
        fit_ARIMA <- fit
        minPsi <- psi
      }
    }

  xhat_ARIMA <- yStg1[1:n] - fit_ARIMA$residuals
  yhat_ARIMA <- predict(fit_ARIMA, n.ahead = p)$pred

  fit_TBATS <- tbats(yStg1[1:n])
  yhat_TBATS <- predict(fit_TBATS, h = p)$mean

  xAct <- yStg1[1:n]
  yAct <- yStg1[n + (1:p)]

  stage1 <- list()
  stage1$tbats$model <- fit_TBATS
  stage1$tbats$accuracy <- cbind(twoStage.accuracy(fit_TBATS$fitted.values,xAct),
                                 twoStage.accuracy(yhat_TBATS,yAct))
  stage1$arima$model <- fit_ARIMA
  stage1$arima$accuracy <- cbind(twoStage.accuracy(xhat_ARIMA,xAct),
                                 twoStage.accuracy(yhat_ARIMA,yAct))

  if (~is.null(regMat)){
    matrixForReg <- cbind(ylow, regMat)
    fit_Reg <- lm(ylow ~ ., matrixForReg[1:n,])
    xhat_Reg <- fit_Reg$fitted.values
    yhat_Reg <- predict(fit_Reg, newdata = regMat[n+(1:p),])

    xregres <- msts(fit_Reg$residuals, seasonal.periods = seas_periods)
    fit_Reg_TBATS <- tbats(xregres)
    xhat_Reg_TBATS <- fit_Reg$fitted.values + fit_Reg_TBATS$fitted.values
    yhat_Reg_TBATS <- predict(fit_Reg_TBATS, h = p)$mean + yhat_Reg

    ARIMAorder <- c(fit_ARIMA$arma[1],fit_ARIMA$arma[6], fit_ARIMA$arma[2])
    sARIMAorder <- c(fit_ARIMA$arma[3],fit_ARIMA$arma[7], fit_ARIMA$arma[4])
    if (sum(sARIMAorder) == 0){
      fit_reg_ARIMA <- arima(fit_Reg$residuals, order = ARIMAorder)
    } else {
      fit_reg_ARIMA <- arima(fit_Reg$residuals, order = ARIMAorder,
                             seasonal = list(order = sARIMAorder,
                                             period = fit_ARIMA$arma[5]))
    }
    xhat_Reg_ARIMA <- fit_Reg$fitted.values - fit_reg_ARIMA$residuals # Regression + ARIMA
    yhat_Reg_ARIMA <-  yhat_Reg - predict(fit_reg_ARIMA, n.ahead = p)$pred

    stage1$reg$model <- fit_Reg
    stage1$reg$accuracy <- cbind(twoStage.accuracy(xhat_Reg,xAct),
                                twoStage.accuracy(yhat_Reg,yAct))
    stage1$regResTBATS$model <- fit_Reg_TBATS
    stage1$regResTBATS$accuracy <- cbind(twoStage.accuracy(xhat_Reg_TBATS,xAct),
                                         twoStage.accuracy(yhat_Reg_TBATS,yAct))

    stage1$regResARIMA$model <- fit_reg_ARIMA
    stage1$regResARIMA$accuracy <- cbind(twoStage.accuracy(xhat_Reg_ARIMA,xAct),
                             twoStage.accuracy(yhat_Reg_ARIMA,yAct))


    if(plotFlag){

      # Stage 1 plot
      plot(1:(n+p),ylow , type = 'l', lwd = 3, col = "blue")
      lines(1:n, ylow[1:n], lwd = 3)
      lines(c(xhat_Reg[1:n],yhat_Reg), col= 'red')
      lines(c(fit_TBATS$fitted.values[1:n],yhat_TBATS), col = 'cyan')
      lines(c(xhat_ARIMA[1:n],yhat_ARIMA), col = 'orange')
      lines(c(xhat_Reg_TBATS[1:n], yhat_Reg_TBATS), col = 'magenta')
      lines(c(xhat_Reg_ARIMA[1:n], yhat_Reg_ARIMA), col = 'green')

      legend("bottomright",c("Actual Model Period","Actual Holdout",
                             "Regression", "TBATS", "ARIMA", "TBATS on Reg. Residuals",
                             "ARIMA on Reg. Residuals"),
             col = c('black','blue','red','cyan', 'orange', 'magenta', 'green'),
             lwd = c(2,2,1,1,1,1,1),
             bty = 'n')
    }

  } else {

    if(plotFlag){

      # Stage 1 plot
      plot(1:(n+p),ylow , type = 'l', lwd = 3, col = "blue")
      lines(1:n, ylow[1:n], lwd = 3)
      lines(c(fit_TBATS$fitted.values[1:n],yhat_TBATS), col = 'cyan')
      lines(c(xhat_ARIMA[1:n],yhat_ARIMA), col = 'orange')

      legend("bottomright",c("Actual Model Period","Actual Holdout",
                              "TBATS", "ARIMA"),
             col = c('black','blue','cyan', 'orange'),
             lwd = c(2,2,1,1),
             bty = 'n')
    }

  }
  return(stage1)
}
