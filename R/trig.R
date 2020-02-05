#' @title  Trigonometric estimation of seasonality
#'
#' @description Estimate the seasonality using trigonometric functions.
#'
#' @param x data whose seasonality is to be estimated.
#'
#'
#' @examples
#'

twoStage.trig <- function(x){
  N <- length(x)
  X <- fft(x)/N
  halfX <- abs(X[1:(N/2)]) # + X[seq(n,(n/2+1),-1)]
  idx <- order(halfX, decreasing = TRUE)

  nrCos <- 3
  #apprCos <- matrix(0, nrow = nrCos, ncol = N)
  apprCos <- rep(0,N)
  psi <- rep(0,nrCos)

  #apprCos[1,] <- 2*halfX[idx[1]]*cos(2*pi*((idx[1]-1)/N)*(0:(N-1)) +
  #                                     atan2(Im(X[idx[1]]),Re(X[idx[1]])))
  #apprCos[2,] <- apprCos[1,] + 2*halfX[idx[2]]*cos(2*pi*((idx[2]-1)/N)*(0:(N-1)) +
  #                                     atan2(Im(X[idx[2]]),Re(X[idx[2]])))
  for (i in 1:nrCos) {
    #apprCos[i,] <- colSums(apprCos[1:(i-1),]) + 2*halfX[idx[i]]*cos(2*pi*((idx[i]-1)/N)*(0:(N-1)) +
    #                                                                  atan2(Im(X[idx[i]]),Re(X[idx[i]])))
    apprCos <- apprCos + 2*halfX[idx[i]]*cos(2*pi*((idx[i]-1)/N)*(0:(N-1)) +
                                               atan2(Im(X[idx[i]]),Re(X[idx[i]])))

    #sse <- sum((x - apprCos[i])^2)
    #psi[i] <- N*log(sse) + 2*(3*i)
  }
  #minPsiIdx <- which.min(psi) # Pick the sum of cosines which gives min psi

  #return(apprCos[minPsiIdx,])
  return(apprCos)
}
