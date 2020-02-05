#' @title  Polynomial fit for seasonality estimation
#'
#' @description Estimates seasonality using polynomial fitting. Currently
#'     this works only for daily and weekly level of seasonality.
#'
#' @param d data whose seasonality is to be estimated.
#'
#' @param freq period of the seasonality to be estimated.
#'
#' @param nrDays number of days in a week in the input data.
#'
#' @examples
#'


twoStage.poly <- function(d, freq, nrDays){
  l <- length(d)
  dayVec <- rep(1:nrDays, ceiling(l/nrDays))
  dayNumber <- dayVec[1:l]

  meanBy <- as.vector(t(matrix(rep(c(1:(nrow(d)/freq)),freq),
                               nrow = nrow(d)/freq)))
  davg <- aggregate(d,by = list(meanBy),mean)$x

  davg_rep <- as.vector(matrix(rep(davg,freq),nrow=freq,byrow=TRUE))
  dMult <- d/davg_rep
  dAdd <- d - davg_rep

  cpoly_season_add <- matrix(0,freq,nrDays)
  cpoly_season_mult <- matrix(1,freq,nrDays)


  ### Sunday to Saturday
  for (i in 1:nrDays) {
    dayIdx <- which(dayNumber == i)
    dday_mult <- dMult[dayIdx]
    dday_add <- dAdd[dayIdx]

    nday_model <- length(dday_mult)/freq		# no of days irrespective of multiplicative or additive
    xday <- rep((1:freq)/freq,nday_model)
    u1 <- xday - xday^4
    u2 <- xday^2 - xday^4
    u3 <- xday^3 - xday^4

    multFit <- lm(dday_mult ~ u1+u2+u3 )
    addFit <- lm(dday_add ~ u1+u2+u3)

    cpoly_season_mult[,i] <- multFit$fitted.values[1:freq]
    cpoly_season_add[,i] <- addFit$fitted.values[1:freq]
  }

  seas <- list()
  seas$mult <- cpoly_season_mult
  seas$add <- cpoly_season_add
  return(seas)
}
