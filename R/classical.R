#' @title  Classical decomposition method for seasonality estimation
#'
#' @description Estimates the seasonality at different levels by
#'     subtracting the mean for additive seasonality
#'     (dividing for multiplicative) and then averaging across
#'     the entire data.
#'
#' @param y the actual data.
#'
#' @param yhat the forecasted vector.
#'
#'
#' @examples
#'

twoStage.classical <- function(data) {
  BlockFreq   <- 12
  NrHrsInADay <- 24
  DayFreq     <- 288
  NrDays      <- 7

  dayNames <- c("sun","mon","tue","wed","thu","fri","sat")

  n <- nrow(data)

  #------------------------------------------------------
  # Block seasonality
  #------------------------------------------------------
  d_hour <- colMeans(matrix(data[,1], nrow = BlockFreq))
  # Replicate this 12 times and reshape into a vector so that
  # it is the same length as d
  d_hour_rep <- as.vector(t(matrix(t(rep(t(d_hour), BlockFreq)), ncol = BlockFreq)))
  # Get block seasonal index
  if (seasType == "m") {
    b.ratio <- data[,1]/d_hour_rep
  } else {
    b.ratio <- data[,1] - d_hour_rep
  }
  b.day <- rowSums(data[,10:16]*matrix(rep(c(1:7),n),nrow=n,byrow=TRUE))
  tmp <- rep(c(1:DayFreq),ceiling(n/DayFreq))
  b.block <- tmp[1:n]

  a <- aggregate(b.ratio,list(day = b.day, block = b.block),mean)
  block <- matrix(a[,3],nrow = DayFreq, ncol = NrDays, byrow=TRUE)
  colnames(block) <- dayNames

  #--------------------------------------------------------------------------
  # Hourly seasonality indices
  #--------------------------------------------------------------------------
  # First get daily averages and then divide hourly averages by daily averages
  d_day <- colMeans(matrix(data[,1], nrow = DayFreq))
  # Replicate by 24 so that we have the same value for each hour in a day
  d_day_rep <- as.vector(t(matrix(t(rep(t(d_day), NrHrsInADay)), ncol = NrHrsInADay)))
  # Get hourly seasonal index
  if (seasType == "m") {
    b2.ratio <- d_hour/d_day_rep
  } else {
    b2.ratio <- d_hour - d_day_rep
  }
  tmp <- rep(1:NrHrsInADay,ceiling(length(d_hour)/NrHrsInADay))#as.vector(t(matrix(rep(1:NrHrsInADay,ceiling(length(d_hour)/NrHrsInADay)), nrow=NrHrsInADay)))
  b2.hour <- tmp[1:length(d_hour)]
  b2.day <- b.day[seq(1,length(b.day),DayFreq/NrHrsInADay)]
  a2 <- aggregate(b2.ratio,list(day = b2.day, hour = b2.hour),mean)
  hourly <- matrix(a2[,3],nrow = NrHrsInADay, byrow = TRUE)
  colnames(hourly) <- dayNames


  # Initialize list to return the values
  SI <- list()
  SI$block  <- block
  SI$hourly <- hourly
  #SI$daily  <- daily
  return(SI)

}
