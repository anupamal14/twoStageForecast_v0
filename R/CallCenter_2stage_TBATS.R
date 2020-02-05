library(forecast)
# setwd("C:/Anupama/FPM/RA/Forecasting/TwoStage")

Acc <- function(Yhat,Y) {
  err <- Y - Yhat
  MSE <- round(mean(err^2), 4)
  MAD <- round(mean(abs(err)), 4)
  MAPE <- round(mean(abs(err/Y)) * 100, 4)
  ErrVec <- cbind(MSE, MAD, paste(toString(MAPE),"%"))
  colnames(ErrVec) <- c("MSE","MAD","MAPE")
  return(ErrVec)
}

classicalSeasonality <- function(data, seasType) {
  
  BlockFreq   <- 13
  NrHrsInADay <- 13
  DayFreq     <- 169
  NrDays      <- 5
  
  n1 <- length(data)
  n <- floor(n1/(169*5))*169*5
  data <- data[1:n]
  
  #------------------------------------------------------
  # Block seasonality
  #------------------------------------------------------
  d_hour <- colMeans(matrix(data, nrow = BlockFreq))
  # Replicate this 12 times and reshape into a vector so that 
  # it is the same length as d
  d_hour_rep <- as.vector(t(matrix(t(rep(t(d_hour), BlockFreq)), ncol = BlockFreq)))
  # Get block seasonal index
  if (seasType == "m") {
    b.ratio <- data/d_hour_rep
  } else {
    b.ratio <- data - d_hour_rep
  }
  
  
  b.day <- as.vector(matrix(rep(c(1:NrDays),n/NrDays),nrow=n/NrDays,byrow=TRUE))
  tmp <- rep(c(1:DayFreq),ceiling(n/DayFreq))
  b.block <- tmp[1:n]
  
  a <- aggregate(b.ratio,list(day = b.day, block = b.block),mean)
  block <- matrix(a[,3],nrow = DayFreq, ncol = NrDays, byrow=TRUE)
  #colnames(block) <- dayNames
  
  #--------------------------------------------------------------------------
  # Hourly seasonality indices
  #--------------------------------------------------------------------------
  # First get daily averages and then divide hourly averages by daily averages
  d_day <- colMeans(matrix(data, nrow = DayFreq))
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
  #colnames(hourly) <- dayNames
  
  
  # Initialize list to return the values
  SI <- list()
  SI$block  <- block
  SI$hourly <- hourly
  #SI$daily  <- daily
  return(SI)
  
}

trigSeas <- function(x)
{
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


# (b) Number of calls handled on weekdays between 7:00 am and 9:05 pm
# Five-minute call volume from March 3, 2003, to May 23, 2003
# in a large North American commercial bank.
#
# Periodicity comes from 14 hours - 7 am to 9 pm = 14*12 (5 min) = 168 + 1 (9:05 pm record)
# Livera paper indicates only 11 weeks and 4 days of data but we see 33 weeks!
#
calls <- unlist(read.csv("callcenter.txt",
                         header=TRUE,sep="\t"))
calls <- calls[170:length(calls)] # Something seems to be wrong with the first day, dropping it
TotNr <- length(calls)
freq1 <- 169

# Forecast period 20 days (1 month) => 20*169
n <- TotNr - 30*freq1

# Seperation fixed for now - hourly data for Stage 1
#---------------------------------------------------

meanBy <- as.vector(t(matrix(rep(c(1:(TotNr/freq1)),freq1),
                             nrow = TotNr/freq1)))
data2 <- aggregate(calls,by = list(meanBy),mean)
callsDaily <- data2$x


x <- msts(calls[1:n], start=2003 + (31+28+2)/365.25, 
          seasonal.periods = c(freq1, freq1*5)) 
fitcalls <- tbats(x)
tbatsModelErr <- Acc(fitcalls$fitted.values, calls[1:n])

tbatsModelPsi <- n*log(as.numeric(tbatsModelErr[1])) + 2*length(fitcalls$parameters$vect) 


yhat_fulltbats <- predict(fitcalls, h = (TotNr - n))$mean
tbatsHoldErr <- Acc(yhat_fulltbats, calls[(n+1):TotNr])
tbatsHoldPsi <- (TotNr - n)*log(as.numeric(tbatsHoldErr[1])) + 2*length(fitcalls$parameters$vect)

# Two-stage
xdailyts <- msts(callsDaily[1:(n/freq1)],start = 2003+ (31+28+2)/365.25,
                  seasonal.periods = 5)
dayWeek <- rep(1:5, ceiling(TotNr/freq1))
matrixForReg <- data.frame(NrCalls = callsDaily, DNo = c(1:(TotNr/freq1)), 
                           dayWeek = as.factor(dayWeek[1:(TotNr/freq1)]))
colnames(matrixForReg) <- c("NrCalls","DNo","dayWeek")
fit_Reg <- lm(NrCalls ~ DNo + dayWeek, matrixForReg[c(1:(n/freq1)),])
xhat_Reg <- fit_Reg$fitted.values
xhat_Reg_rep <- as.vector(matrix(rep(xhat_Reg,freq1),nrow=freq1,byrow=TRUE))
newMatrix <- data.frame(DNo = matrixForReg$DNo[c(((n/freq1)+1):(TotNr/freq1))],
                        dayWeek = matrixForReg$dayWeek[c(((n/freq1)+1):(TotNr/freq1))])
yhat_Reg <- predict(fit_Reg, newdata = newMatrix)
yhat_Reg_rep <- as.vector(matrix(rep(yhat_Reg,freq1),nrow=freq1,byrow=TRUE))

# aic <- matrix(NA, nrow = 6*2, ncol = 6)
# for (d in c(1,0)){
# 
# for (p in 0:5) {
#   for (q in 0:5) {
#       model <- arima(xdailyts, order = c(p,d,q))
#       aic[(d*6)+p+1,q+1] <- model$aic
#     }
#   }
# }
# modelIdx <- which.min(aic)
#
# Best model from above model is (4,0,5)
#

fit_ARIMA <- arima(xdailyts, order = c(4,0,5))
xhat_ARIMA <- xdailyts - fit_ARIMA$residuals
xhat_ARIMA_rep <- as.vector(matrix(rep(xhat_ARIMA,freq1),nrow=freq1,byrow=TRUE))
yhat_ARIMA <- predict(fit_ARIMA, n.ahead = (TotNr - n)/freq1)$pred
yhat_ARIMA_rep <- as.vector(matrix(rep(yhat_ARIMA,freq1),nrow=freq1,byrow=TRUE))

fit_TBATS <- tbats(xdailyts)
xhat_TBATS_rep <- as.vector(matrix(rep(fit_TBATS$fitted.values,freq1),nrow=freq1,byrow=TRUE))

yhat_TBATS <- predict(fit_TBATS, h = (TotNr - n)/freq1)$mean
yhat_TBATS_rep <- as.vector(matrix(rep(yhat_TBATS,freq1),nrow=freq1,byrow=TRUE))

xregres <- msts(fit_Reg$residuals,start = 2003+ (31+28+2)/365.25,
                 seasonal.periods = 5)
fit_Reg_TBATS <- tbats(xregres)
xhat_Reg_TBATS <- fit_Reg$fitted.values + fit_Reg_TBATS$fitted.values
yhat_Reg_TBATS <- predict(fit_Reg_TBATS, h = (TotNr - n)/freq1)$mean + yhat_Reg

fit_reg_ARIMA <- arima(fit_Reg$residuals, order = c(4,0,5))
xhat_Reg_ARIMA <- fit_Reg$fitted.values + fit_Reg$residuals - fit_reg_ARIMA$residuals # Regression + ARIMA
yhat_Reg_ARIMA <-  yhat_Reg - predict(fit_reg_ARIMA, n.ahead = (TotNr - n)/freq1)$pred

#------------------------------------------------------
#CLASSICAL
#------------------------------------------------------
forecastLen <- TotNr - n

si <- classicalSeasonality(calls[1:n],"m")
# Now unwrap the indices for the entire week and then replicate to length of predicted period.
testStartInd <- ((n/freq1)%%5)*freq1 
SImult_vec1 <- as.vector(si$block)
repNr <- ceiling(length(xhat_Reg_rep)/length(SImult_vec1))
SImultBlk_vec1 <- rep(SImult_vec1, repNr)
SImultBlk_vec <- SImultBlk_vec1[(1:forecastLen)+testStartInd]

SImult_vec1 <- as.vector(si$hourly)
SImulttmp <- as.vector(t(matrix(rep(SImult_vec1,13), ncol=13)))
SImultHr_vec1 <- rep(SImulttmp, repNr)
SImultHr_vec <- SImultHr_vec1[(1:forecastLen)+testStartInd]

# Apply the seasonality corrections
yhat_TBATS_final_mc <- yhat_TBATS_rep*SImultBlk_vec*SImultHr_vec
xhat_TBATS_final_mc <- xhat_TBATS_rep*SImultBlk_vec1[1:n]*SImultHr_vec1[1:n]

yhat_Reg_final_mc <- yhat_Reg_rep*SImultBlk_vec*SImultHr_vec
xhat_Reg_final_mc <- xhat_Reg_rep*SImultBlk_vec1[1:n]*SImultHr_vec1[1:n]

yhat_ARIMA_final_mc <- yhat_ARIMA_rep*SImultBlk_vec*SImultHr_vec
xhat_ARIMA_final_mc <- xhat_ARIMA_rep*SImultBlk_vec1[1:n]*SImultHr_vec1[1:n]

# ADDITIVE
si <- classicalSeasonality(calls[1:n],"a")
# Now unwrap the indices for the entire week and then replicate to length of predicted period.
SImult_vec1 <- as.vector(si$block)
repNr <- ceiling(length(xhat_Reg_rep)/length(SImult_vec1))
SImultBlk_vec1 <- rep(SImult_vec1, repNr)
SImultBlk_vec <- SImultBlk_vec1[(1:forecastLen)+testStartInd]

SImult_vec1 <- as.vector(si$hourly)
SImulttmp <- as.vector(t(matrix(rep(SImult_vec1,13), ncol=13)))
SImultHr_vec1 <- rep(SImulttmp, repNr)
SImultHr_vec <- SImultHr_vec1[(1:forecastLen)+testStartInd]

# Apply the seasonality corrections
yhat_TBATS_final_ac <- yhat_TBATS_rep + SImultBlk_vec + SImultHr_vec
xhat_TBATS_final_ac <- xhat_TBATS_rep + SImultBlk_vec1[1:n] + SImultHr_vec1[1:n]

yhat_Reg_final_ac <- yhat_Reg_rep + SImultBlk_vec + SImultHr_vec
xhat_Reg_final_ac <- xhat_Reg_rep + SImultBlk_vec1[1:n] + SImultHr_vec1[1:n]

yhat_ARIMA_final_ac <- yhat_ARIMA_rep + SImultBlk_vec + SImultHr_vec
xhat_ARIMA_final_ac <- xhat_ARIMA_rep + SImultBlk_vec1[1:n] + SImultHr_vec1[1:n]

#-----------------------------------------------------
# TRIGONOMETRIC SEASONALITY
#-----------------------------------------------------
# Repeat daily average for every 5min
x_rep <- as.vector(matrix(rep(callsDaily[1:(n/freq1)],freq1),nrow=freq1,byrow=TRUE))

# Additive
inputAdd <- calls[1:n] - x_rep
addTrSeas <- trigSeas(inputAdd)
# Multiplicative
inputMult <- log(calls[1:n]) - log(x_rep)
multTrSeas <- trigSeas(inputMult)


yhat_TBATS_at <- yhat_TBATS_rep + addTrSeas[(1:forecastLen)+testStartInd]
xhat_TBATS_at <- xhat_TBATS_rep + addTrSeas

yhat_TBATS_mt <- exp(log(yhat_TBATS_rep) + multTrSeas[(1:forecastLen)+testStartInd])
xhat_TBATS_mt <- exp(log(xhat_TBATS_rep) + multTrSeas) 



yhat_ARIMA_at <- yhat_ARIMA_rep + addTrSeas[(1:forecastLen)+testStartInd]
xhat_ARIMA_at <- xhat_ARIMA_rep + addTrSeas

yhat_ARIMA_mt <- exp(log(yhat_ARIMA_rep) + multTrSeas[(1:forecastLen)+testStartInd])
xhat_ARIMA_mt <- exp(log(xhat_ARIMA_rep) + multTrSeas) 



yhat_Reg_at <- yhat_Reg_rep + addTrSeas[(1:forecastLen)+testStartInd]
xhat_Reg_at <- xhat_Reg_rep + addTrSeas

yhat_Reg_mt <- exp(log(yhat_Reg_rep) + multTrSeas[(1:forecastLen)+testStartInd])
xhat_Reg_mt <- exp(log(xhat_Reg_rep) + multTrSeas) 


y <- calls[(n+1):TotNr]
errMatHold <- rbind(Acc(yhat_Reg_final_ac,y),Acc(yhat_Reg_final_mc,y),
                    Acc(yhat_Reg_at,y),Acc(yhat_Reg_mt,y),
                    Acc(yhat_TBATS_final_ac,y),Acc(yhat_TBATS_final_mc,y),
                    Acc(yhat_TBATS_at,y),Acc(yhat_TBATS_mt,y),
                    Acc(yhat_ARIMA_final_ac,y),Acc(yhat_ARIMA_final_mc,y),
                    Acc(yhat_ARIMA_at,y),Acc(yhat_ARIMA_mt,y),
                    tbatsHoldErr)
                    
errMatModel <- rbind(Acc(xhat_Reg_final_ac,calls[1:n]),Acc(xhat_Reg_final_mc,calls[1:n]),
                    Acc(xhat_Reg_at,calls[1:n]),Acc(xhat_Reg_mt,calls[1:n]),
                    Acc(xhat_TBATS_final_ac,calls[1:n]),Acc(xhat_TBATS_final_mc,calls[1:n]),
                    Acc(xhat_TBATS_at,calls[1:n]),Acc(xhat_TBATS_mt,calls[1:n]),
                    Acc(xhat_ARIMA_final_ac,calls[1:n]),Acc(xhat_ARIMA_final_mc,calls[1:n]),
                    Acc(xhat_ARIMA_at,calls[1:n]),Acc(xhat_ARIMA_mt,calls[1:n]),
                    tbatsModelErr)


nrParReg <- length(fit_Reg$coefficients) + 1 
nrParTBATS <- length(fit_TBATS$parameters$vect) + 1 
nrParARIMA <- length(fit_ARIMA$parameters$vect) + 1

nrParClassSeas <- 5*freq1 + 5*13
nrParTrigSeas <- 12

nrParVec <- cbind(nrParReg + nrParClassSeas, nrParReg + nrParClassSeas, 
                  nrParReg + nrParTrigSeas, nrParReg + nrParTrigSeas,
                  nrParTBATS + nrParClassSeas, nrParTBATS + nrParClassSeas,
                  nrParTBATS + nrParTrigSeas, nrParTBATS + nrParTrigSeas,
                  nrParARIMA + nrParClassSeas, nrParARIMA + nrParClassSeas, 
                  nrParARIMA + nrParTrigSeas, nrParARIMA + nrParTrigSeas)

psiModel <- n*log(as.numeric(errMatModel[1:12,1])) + 2*nrParVec
psiModel <- cbind(psiModel,tbatsModelPsi)
psiHold <- forecastLen*log(as.numeric(errMatHold[1:12,1])) + 2*nrParVec
psiHold <- cbind(psiHold,tbatsHoldPsi)
psiOutput <- cbind(t(psiModel), t(psiHold))
colnames(psiOutput) <- c('Model','Hold out')


errorMatColNames <- c("Regression; Additive; Classical","Regression; Multiplicative; Classical",
                      "Regression; Additive; Trig","Regression; Multiplicative; Trig",
                      "TBATS with weekly seas.; Additive; Classical","TBATS with weekly seas.; Multiplicative; Classical",
                      "TBATS with weekly seas.; Additive; Trig","TBATS with weekly seas.; Multiplicative; Trig",
                      "ARIMA; Additive; Classical","ARIMA; Multiplicative; Classical",
                      "ARIMA; Additive; Trig","ARIMA; Multiplicative; Trig",
                      "TBATS with 2 levels of Seasonality")
                      
rownames(errMatHold) <- errorMatColNames
# write.csv(errMatHold,"CallCenterForecastErrors2.csv")
# 
rownames(errMatModel) <- errorMatColNames
# write.csv(errMatModel,"CallCenterModelErrors2.csv")
# 
# rownames(psiOutput) <- errorMatColNames
# write.csv(psiOutput,"CallCenterPsi2.csv")
