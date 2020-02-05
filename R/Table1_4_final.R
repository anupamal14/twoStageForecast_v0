library(forecast)
library(lubridate)
setwd("C:/Anupama/FPM/RA/Forecasting/Paper1TwoStage")

Acc <- function(Yhat,Y) {
  err <- Y - Yhat
  MSE <- round(mean(err^2), 1)
  MAD <- round(mean(abs(err)), 1)
  MAPE <- round(mean(abs(err/Y)) * 100, 2)
  ErrVec <- cbind(MSE, MAD, paste(toString(MAPE),"%"))
  colnames(ErrVec) <- c("MSE","MAD","MAPE")
  return(ErrVec)
}

classicalSeasonality <- function(data, seasType) {
  
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

polySeasEst <- function(d) {
  freq <- 288 #for daily average
  meanBy <- as.vector(t(matrix(rep(c(1:(nrow(d)/freq)),freq),
                               nrow = nrow(d)/freq)))
  davg <- aggregate(d$load,by = list(meanBy),mean)$x
  
  davg_rep <- as.vector(matrix(rep(davg,freq),nrow=freq,byrow=TRUE)) 
  dMult <- d$load/davg_rep
  dAdd <- d$load - davg_rep
  
  cpoly_season_add <- matrix(0,freq,7)
  cpoly_season_mult <- matrix(1,freq,7)
  
  
  ### Sunday to Saturday
  for (i in 1:7) {
    dayIdx <- which(d$dayweek == i)
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



#d2 <- read.csv("C:/Users/lanupama/Dropbox/Anupama and SD/Forecasting/NYC_final.csv",header = TRUE)
d2 <- read.csv("cleanNYdataFinal.csv", header = TRUE)
d2[,1] <- na.interp(d2[,1])
# Remove 29th Feb
feb29ind <- grep("29-Feb",d2[,2])
d <- d2[-feb29ind,]

#d <- cbind(d,d[,c(19:22)]^2)
#colnames(d) <- c("load","date","time","datetime","nobs","yeardate","weekyear","dayweek",
#                 "DNo","sun","mon","tue","wed","thu","fri","sat",
#                 "weekday","weekend","tmax","tmin","bmax","bmin","tmaxSq","tminSq","bmaxSq","bminSq")
d <- cbind(d,d[,c(19:22)]^2, d[,c(19:22)]^3)
monthCol <- d$dayweek
monthCol[grep("Jan",d$date)] <- 1
monthCol[grep("Feb",d$date)] <- 2
monthCol[grep("Mar",d$date)] <- 3
monthCol[grep("Apr",d$date)] <- 4
monthCol[grep("May",d$date)] <- 5
monthCol[grep("Jun",d$date)] <- 6
monthCol[grep("Jul",d$date)] <- 7
monthCol[grep("Aug",d$date)] <- 8
monthCol[grep("Sep",d$date)] <- 9
monthCol[grep("Oct",d$date)] <- 10
monthCol[grep("Nov",d$date)] <- 11
monthCol[grep("Dec",d$date)] <- 12
d <- cbind(d, as.factor(monthCol))
d$dayweek <- as.factor(d$dayweek)

colnames(d) <- c("load","date","time","datetime","nobs","yeardate","weekyear","dayweek",
                 "DNo","sun","mon","tue","wed","thu","fri","sat",
                 "weekday","weekend","tmax","tmin","bmax","bmin","tmaxSq","tminSq","bmaxSq","bminSq",
                 "tmax3","tmin3","bmax3","bmin3","month")

# Seperation fixed for now - daily data for Stage 1
#---------------------------------------------------
freq1 <- 288

meanBy <- as.vector(t(matrix(rep(c(1:(nrow(d)/freq1)),freq1),
                             nrow = nrow(d)/freq1)))
dataAvg <- d[seq(1,nrow(d),freq1),]
data2 <- aggregate(d$load,by = list(meanBy),mean)
dataAvg$load <- data2$x
dataAvg$DNo <- c(1:nrow(dataAvg))




# Separation into model and hold out periods
ind <- which(dataAvg$date == "28-02-2015") #"31-Dec-14")
n <- ind[length(ind)]
# Get end period for the 5 min dataset
ind2 <- which(d$date == "28-02-2015") #"31-Dec-14")
n2 <- ind2[length(ind2)]

DL <- dataAvg[1:n,]


# Wet and dry bulb temperature forcasts using NYISO forecasts
tempForecast <- read.csv("weatherForecastSep08ToDec15.csv", header = TRUE)

fcstDate2 <- parse_date_time(tempForecast[,1],c("dmy","mdy"))
fcstDate <- gsub(" UTC", "", fcstDate2)

feb29fcst <-  which((month(fcstDate) == 2) &
                      (day(fcstDate) == 29))
tempForecast <- tempForecast[-feb29fcst,]
fcstDate <- fcstDate[-feb29fcst]
fcstDateDay <- fcstDate[seq(1,length(fcstDate), 288)]

# Find the start day index of the test period
testPeriodDayStartInd <- d$dayweek[n2+1]

# F2- model
fitF2 <- lm(load ~ DNo + tmax + tmin + bmax + bmin +
              tmaxSq + tminSq + bmaxSq + bminSq + weekday + sat, data = DL)

fitF2minus <- lm(load ~ DNo + tmax + tmin + bmax + bmin +
                   tmaxSq + tminSq + bminSq + weekday + sat, data = DL)

regCoeff <- round(t(summary(fitF2minus)$coefficients[,c(1,2,4)]),2)
rownames(regCoeff)[3] <- "Sig."
#write.csv(regCoeff,"Table1_1.csv")

fitBReg <- lm(load ~ DNo + tmaxSq + tminSq + tmax3 + tmin3 + bmax3 + bmin3 + dayweek + month
              + tmaxSq*month + tminSq*month + tmin3*month + tmax3*month + bmin3*month + bmax3*month, data = DL)
fitF2minus <- fitBReg
forecastLen <- 365
newDATA <- dataAvg[c((n+1):(n+forecastLen)),c(9,23,24,27,28,29,30,8,31)]
yhat1 <- predict(fitF2minus,newdata=newDATA)

# For TBATS  
xMSTS <- msts(fitF2minus$residuals, start=8+248/365, seasonal.periods=c(7,365))

xMSTS2 <- msts(DL$load, start=8+248/365, seasonal.periods=c(7,365))


fit_Reg_TBATS <- tbats(xMSTS)
xhat_Reg_TBATS <- fitF2minus$fitted.values + fit_Reg_TBATS$fitted.values
# Replicate to get 288 values for each day
xhat_Reg_TBATS_full <- as.vector(matrix(rep(xhat_Reg_TBATS,freq1),nrow=freq1,byrow=TRUE)) # Regression + TBATS
xhat_Reg_full <- as.vector(matrix(rep(fitF2minus$fitted.values,freq1),nrow=freq1,byrow=TRUE)) # Only reg

# (p,0,q) = (5,0,5) -- in terms of AIC
fit_Reg_ARIMA505 <- arima(fitF2minus$residuals, order = c(5,0,5)) 
xhat_Reg_ARIMA505 <- fitF2minus$fitted.values + fitF2minus$residuals - fit_Reg_ARIMA505$residuals # Regression + ARIMA
xhat_Reg_ARIMA505_full <- as.vector(matrix(rep(xhat_Reg_ARIMA505,freq1),nrow=freq1,byrow=TRUE)) # Regression + ARIMA
# (p,0,q) = (3,0,1) -- in terms of BIC
fit_Reg_ARIMA301 <- arima(fitF2minus$residuals, order = c(3,0,1)) 
xhat_Reg_ARIMA301 <- fitF2minus$fitted.values + fitF2minus$residuals - fit_Reg_ARIMA301$residuals # Regression + ARIMA
xhat_Reg_ARIMA301_full <- as.vector(matrix(rep(xhat_Reg_ARIMA301,freq1),nrow=freq1,byrow=TRUE)) # Regression + ARIMA

fit_ARIMA505 <- arima(DL$load, order = c(5,0,5)) 
xhat_ARIMA505 <- fitted(fit_ARIMA505)# Regression + ARIMA
yhat_ARIMA505 <- predict(fit_ARIMA505, n.ahead = forecastLenSt1)$pred

fitTBATSwk_yr <- tbats(xMSTS2)
xhat_TBATSwk_yr_full <- as.vector(matrix(rep(fitTBATSwk_yr$fitted.values,freq1),nrow=freq1,byrow=TRUE))

#fit_ARIMA_xreg <- arima(DL$load, order = c(5,0,5), seasonal = list(order = c(0,0,0), period = 365),
#              xreg = cbind(DL$DNo, DL$tmax, DL$tmin, DL$bmax, DL$bmin, DL$tmaxSq, DL$tminSq, DL$bminSq,
#                           DL$weekday, DL$sat)) 
fit_ARIMA_xreg <- arima(DL$load, order = c(5,0,5), seasonal = list(order = c(0,0,0), period = 365),
                        xreg = cbind(DL$DNo, DL$tmaxSq, DL$tminSq, 
                                     DL$tmax3, DL$tmin3, DL$bmax3, DL$bmin3)) 

xhat_ARIMAxreg <- DL$load - fit_ARIMA_xreg$residuals
xhat_ARIMAxreg_full <- as.vector(matrix(rep(xhat_ARIMAxreg,freq1),nrow=freq1,byrow=TRUE))


#arima(DL$load, 	order = c(3, 1, 5), seasonal = list(order = c(0, 1, 0), 
# xreg= cbind(1:n, DL$weekday, DL$sat, DL$tmax, DL$tmin, DL$bmax,
# DL$bmin, DL$tmaxSq, DL$tminSq, DL$bminSq), period = 365))


#----------------------------
# Check for best ARIMA model
#----------------------------
#  options(warn=1)
 #  aic <- matrix(NA,6,6)
 #  bic <- matrix(NA,6,6)
 # for (p in 0:5) {
 #   for (q in 0:5) {
 #    model <- arima(DL$load, 	order = c(p, 0, q), seasonal = list(order = c(0, 1, 0),
 #                            xreg = cbind(DL$DNo, DL$tmaxSq, DL$tminSq, 
 #                            DL$tmax3, DL$tmin3, DL$bmax3, DL$bmin3, DL$month,
 #                            DL$dayweek),period=365)) 
 #                   #xreg= cbind(1:n, DL$weekday, DL$sat, DL$tmax, DL$tmin, DL$bmax,
 #                   #DL$bmin, DL$tmaxSq, DL$tminSq, DL$bminSq), period = 365))
 #    #models[[p+1]][[q+1]] <- model
 #     aic[p+1,q+1] <- model$aic
 #     bic[p+1,q+1] <- AIC(model, k= log(n))
 #     print(cat("p=",p,", q=",q))
 # 
 #   }
 # }
# write.csv(cbind(aic,bic),"ARIMAxreg_AIC_BIC.csv")
#
#-----------------------------------------------------------------------------------------

# Forecast - 1 day and 1 week
forecastLenSt1 <- 365
rowsIndices2 <- c(n2+(1:(forecastLenSt1*freq1)))
y <- d[rowsIndices2,1]
rowsIndices <- c(n+(1:forecastLenSt1))

forecastMat <- dataAvg[rowsIndices,c(9,23,24,27:30)]
yhat_Reg<- predict(fitF2minus, newdata = newDATA)
yhat_Reg_TBATS <- yhat_Reg + as.vector(predict(fit_Reg_TBATS, h=forecastLenSt1)$mean)

yhat_Reg_ARIMA505 <-  yhat_Reg + predict(fit_Reg_ARIMA505, n.ahead = forecastLenSt1)$pred
yhat_Reg_ARIMA301 <-  yhat_Reg + predict(fit_Reg_ARIMA301, n.ahead = forecastLenSt1)$pred

yhatTBATS_wk_yr <- as.vector(predict(fitTBATSwk_yr,h=forecastLenSt1)$mean)

yhat_ARIMAxreg <- predict(fit_ARIMA_xreg, newxreg = forecastMat)$pred

# Replicate to get 288 values for each day
yhat_Reg_TBATS_full <- as.vector(matrix(rep(yhat_Reg_TBATS,freq1),nrow=freq1,byrow=TRUE))
yhat_Reg_ARIMA505_full <- as.vector(matrix(rep(yhat_Reg_ARIMA505,freq1),nrow=freq1,byrow=TRUE))
yhat_Reg_ARIMA301_full <- as.vector(matrix(rep(yhat_Reg_ARIMA301,freq1),nrow=freq1,byrow=TRUE))
yhat_Reg_full <- as.vector(matrix(rep(yhat_Reg,freq1),nrow=freq1,byrow=TRUE))
yhatTBATS_wk_yr_full <- as.vector(matrix(rep(yhatTBATS_wk_yr,freq1),nrow=freq1,byrow=TRUE))
yhat_ARIMAxreg_full <- as.vector(matrix(rep(yhat_ARIMAxreg,freq1),nrow=freq1,byrow=TRUE))


# TABLE 2
#--------------
xMSTS3 <- msts(DL$load, start=8+248/365, seasonal.periods=7)
fitTBATS_wk <- tbats(xMSTS3)
xMSTS4 <- msts(DL$load, start=8+248/365, seasonal.periods=365)
fitTBATS_yr <- tbats(xMSTS4)
ySt1 <- dataAvg$load[(n+1):(n+forecastLenSt1)]
errStage1 <- rbind(cbind(Acc(fitF2minus$fitted.values,DL$load),Acc(yhat_Reg,ySt1)),
                   cbind(Acc(fitTBATS_wk$fitted.values,DL$load),Acc(as.vector(predict(fitTBATS_wk,h=forecastLenSt1)$mean),ySt1)),
                   cbind(Acc(fitTBATS_yr$fitted.values,DL$load),Acc(as.vector(predict(fitTBATS_yr,h=forecastLenSt1)$mean),ySt1)),
                   cbind(Acc(fitTBATSwk_yr$fitted.values,DL$load),Acc(yhatTBATS_wk_yr,ySt1)),
                   cbind(Acc(xhat_Reg_TBATS,DL$load),Acc(yhat_Reg_TBATS,ySt1)),
                   cbind(Acc(xhat_Reg_ARIMA505,DL$load),Acc(yhat_Reg_ARIMA505,ySt1)),
                   cbind(Acc(xhat_Reg_ARIMA301,DL$load),Acc(yhat_Reg_ARIMA301,ySt1)),
                   cbind(Acc(xhat_ARIMAxreg,DL$load),Acc(yhat_ARIMAxreg,ySt1)))


rownames(errStage1) <- c("Regression","TBATS weekly",
                         "TBATS annual","TBATS weekly & annual",
                         "Reg + TBATS (res)",
                         "Reg + ARIMA on res (5,0,5)","Reg + ARIMA on res (3,0,1)",
                         "ARIMA w xreg")
write.csv(errStage1,"Table2_rev2.csv")

#------------------------------------------------------------------------------------------------
# MULTIPLICATIVE SEASONALITY
# Get the seasonal corrections for block and hour because we are not using them in the TBATS call
# and apply to the predicted values

siMult <- classicalSeasonality(d[1:n2,],"m")

forecastLen <- 365*freq1

# Construct a matrix, with rows having indices 1:288 repeated for number of days and columns as
# the dayweek indicator
modelBlkSeasAppl <- cbind(rep(1:freq1, n), d$dayweek[1:n2])
holdBlkSeasAppl <- cbind(rep(1:freq1, 365), d$dayweek[(1:forecastLen)+n2])
tmpHrIdx <- as.vector(matrix(rep(1:24, 12),nrow=12,byrow=TRUE))
modelHrSeasAppl <- cbind(rep(tmpHrIdx, n), d$dayweek[1:n2])
holdHrSeasAppl <- cbind(rep(tmpHrIdx, 365), d$dayweek[(1:forecastLen)+n2])


# Apply the seasonality corrections
yhat_Reg_TBATS_final_mc <- yhat_Reg_TBATS_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_Reg_TBATS_fullmc <- xhat_Reg_TBATS_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]
  

yhat_Reg_ARIMA505_final_mc <- yhat_Reg_ARIMA505_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_Reg_ARIMA505_fullmc <- xhat_Reg_ARIMA505_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]


yhat_Reg_ARIMA301_final_mc <- yhat_Reg_ARIMA301_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_Reg_ARIMA301_fullmc <- xhat_Reg_ARIMA301_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]


yhat_Reg_final_mc <- yhat_Reg_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_Reg_fullmc <- xhat_Reg_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]


yhatTBATS_wk_yr_final_mc <- yhatTBATS_wk_yr_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_TBATSwk_yr_fullmc <- xhat_TBATSwk_yr_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]


yhat_ARIMAxreg_final_mc <- yhat_ARIMAxreg_full*siMult$block[holdBlkSeasAppl]*siMult$hourly[holdHrSeasAppl]
xhat_ARIMAxreg_fullmc <- xhat_ARIMAxreg_full*siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]

if (0) {
# For plot for paper
par(mai=c(.82,0.82,0.32,0.42))
plot(siMult$hourly[,7],xlab='Hours in a day', ylab='Seasonality estimated',
     type='b',ylim = c(0.7,1.2),lwd=2,pch=0)
lines(siMult$hourly[,1],col='red',pch=1,type='b',lwd=2)
lines(siMult$hourly[,2],col='green',pch=2,type='b',lwd=2)
lines(siMult$hourly[,3],col='blue',pch=3,type='b',lwd=2)
lines(siMult$hourly[,4],col='magenta',pch=4,type='b',lwd=2)
lines(siMult$hourly[,5],col='yellow',pch=5,type='b',lwd=2)
lines(siMult$hourly[,6],col='cyan',pch=6,type='b',lwd=2)
legend("bottomright",c('Sat','Sun','Mon','Tue','Wed','Thu','Fri'),
       pch=0:6,col=c('black','red','green','blue','magenta','yellow','cyan'),
       lty=rep(1,7),lwd=rep(2,7))

par(mai=c(.82,0.82,0.32,0.42))
plot(siMult$block[,7],xlab='5 minute intervals', ylab='Seasonality estimated',
     type='b',ylim = c(0.93,1.05),lwd=2,pch=0,xaxt='n')
lines(siMult$block[,1],col='red',pch=1,type='b',lwd=2)
lines(siMult$block[,2],col='green',pch=2,type='b',lwd=2)
lines(siMult$block[,3],col='blue',pch=3,type='b',lwd=2)
lines(siMult$block[,4],col='magenta',pch=4,type='b',lwd=2)
lines(siMult$block[,5],col='yellow',pch=5,type='b',lwd=2)
lines(siMult$block[,6],col='cyan',pch=6,type='b',lwd=2)
legend("bottomright",c('Sat','Sun','Mon','Tue','Wed','Thu','Fri'),
       pch=0:6,col=c('black','red','green','blue','magenta','yellow','cyan'),
       lty=rep(1,7),lwd=rep(2,7))
axis(1,at=seq(40,288,40), labels = c('2:40','6:00','9:20','12:40','16:00','19:20','22:40'))

par(mai=c(.82,0.82,0.32,0.42))
plot(seas$mult[seq(1,288,8),7],xlab='5 minute intervals', ylab='Seasonality estimated',
     type='b',ylim = c(0.75,1.2),lwd=2,pch=0,xaxt='n')
lines(seas$mult[seq(1,288,8),1],col='red',pch=1,type='b',lwd=2)
lines(seas$mult[seq(1,288,8),2],col='green',pch=2,type='b',lwd=2)
lines(seas$mult[seq(1,288,8),3],col='blue',pch=3,type='b',lwd=2)
lines(seas$mult[seq(1,288,8),4],col='magenta',pch=4,type='b',lwd=2)
lines(seas$mult[seq(1,288,8),5],col='yellow',pch=5,type='b',lwd=2)
lines(seas$mult[seq(1,288,8),6],col='cyan',pch=6,type='b',lwd=2)
legend("bottomright",c('Sat','Sun','Mon','Tue','Wed','Thu','Fri'),
       pch=0:6,col=c('black','red','green','blue','magenta','yellow','cyan'),
       lty=rep(1,7),lwd=rep(2,7))
axis(1,at=seq(5,36,5), labels = c('2:40','6:00','9:20','12:40','16:00','19:20','22:40'))

}
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
# ADDITIVE SEASONALITY
# Get the seasonal corrections for block and hour because we are not using them in the TBATS call
# and apply to the predicted values
siAdd <- classicalSeasonality(d[1:n2,],"a")

# Apply the seasonality corrections
yhat_Reg_TBATS_final_ac <- yhat_Reg_TBATS_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_Reg_TBATS_fullac <- xhat_Reg_TBATS_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


yhat_Reg_ARIMA505_final_ac <- yhat_Reg_ARIMA505_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_Reg_ARIMA505_fullac <- xhat_Reg_ARIMA505_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


yhat_Reg_ARIMA301_final_ac <- yhat_Reg_ARIMA301_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_Reg_ARIMA301_fullac <- xhat_Reg_ARIMA301_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


yhat_Reg_final_ac <- yhat_Reg_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_Reg_fullac <- xhat_Reg_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


yhatTBATS_wk_yr_final_ac <- yhatTBATS_wk_yr_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_TBATSwk_yr_fullac <- xhat_TBATSwk_yr_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


yhat_ARIMAxreg_final_ac <- yhat_ARIMAxreg_full + siAdd$block[holdBlkSeasAppl] + siAdd$hourly[holdHrSeasAppl]
xhat_ARIMAxreg_fullac <- xhat_ARIMAxreg_full + siAdd$block[modelBlkSeasAppl] + siAdd$hourly[modelHrSeasAppl]


#------------------------------------------------------------------------------------------------

#-----------------------------------------------------
# POLYNOMIAL SEASONALITY
#-----------------------------------------------------
# Repeat daily average for every 5min
x_rep <- as.vector(matrix(rep(DL$load,freq1),nrow=freq1,byrow=TRUE))

seas <- polySeasEst(d[1:n2,])



# Construct a matrix, with rows having indices 1:288 repeated for number of days and columns as
# the dayweek indicator
modelBlkSeasApplP <- cbind(rep(1:freq1, n), d$dayweek[1:n2])
holdBlkSeasApplP <- cbind(rep(1:freq1, 365), d$dayweek[(1:forecastLen)+n2])


yhat_Reg_TBATS_final_at <- yhat_Reg_TBATS_full + seas$add[holdBlkSeasApplP]
xhat_Reg_TBATS_final_at <- xhat_Reg_TBATS_full + seas$add[modelBlkSeasApplP]


yhat_Reg_TBATS_final_mt <- yhat_Reg_TBATS_full * seas$mult[holdBlkSeasApplP]
xhat_Reg_TBATS_final_mt <- xhat_Reg_TBATS_full * seas$mult[modelBlkSeasApplP]


yhat_Reg_ARIMA505_final_at <- yhat_Reg_ARIMA505_full + seas$add[holdBlkSeasApplP]
xhat_Reg_ARIMA505_final_at <- xhat_Reg_ARIMA505_full + seas$add[modelBlkSeasApplP]


yhat_Reg_ARIMA301_final_at <- yhat_Reg_ARIMA301_full + seas$add[holdBlkSeasApplP]
xhat_Reg_ARIMA301_final_at <- xhat_Reg_ARIMA301_full + seas$add[modelBlkSeasApplP]


yhat_Reg_ARIMA505_final_mt <- yhat_Reg_ARIMA505_full * seas$mult[holdBlkSeasApplP]
xhat_Reg_ARIMA505_final_mt <- xhat_Reg_ARIMA505_full * seas$mult[modelBlkSeasApplP]



yhat_Reg_ARIMA301_final_mt <- yhat_Reg_ARIMA301_full * seas$mult[holdBlkSeasApplP]
xhat_Reg_ARIMA301_final_mt <- xhat_Reg_ARIMA301_full * seas$mult[modelBlkSeasApplP]



yhat_Reg_final_at <- yhat_Reg_full + seas$add[holdBlkSeasApplP]
xhat_Reg_final_at <- xhat_Reg_full + seas$add[modelBlkSeasApplP]



yhat_Reg_final_mt <- yhat_Reg_full * seas$mult[holdBlkSeasApplP]
xhat_Reg_final_mt <- xhat_Reg_full * seas$mult[modelBlkSeasApplP]



yhatTBATS_wk_yr_final_at <- yhatTBATS_wk_yr_full + seas$add[holdBlkSeasApplP]
xhat_TBATSwk_yr_final_at <- xhat_TBATSwk_yr_full + seas$add[modelBlkSeasApplP]



yhatTBATS_wk_yr_final_mt <- yhatTBATS_wk_yr_full * seas$mult[holdBlkSeasApplP]
xhat_TBATSwk_yr_final_mt <- xhat_TBATSwk_yr_full * seas$mult[modelBlkSeasApplP]



yhat_ARIMAxreg_final_at <- yhat_ARIMAxreg_full + seas$add[holdBlkSeasApplP]
xhat_ARIMAxreg_final_at <- xhat_ARIMAxreg_full + seas$add[modelBlkSeasApplP]



yhat_ARIMAxreg_final_mt <- yhat_ARIMAxreg_full * seas$mult[holdBlkSeasApplP]
xhat_ARIMAxreg_final_mt <- xhat_ARIMAxreg_full * seas$mult[modelBlkSeasApplP]


#------------------------------------------------------


errMatHold <- rbind(Acc(yhat_Reg_final_ac,y),Acc(yhat_Reg_final_mc,y),
                    Acc(yhat_Reg_final_at,y),Acc(yhat_Reg_final_mt,y),
                    Acc(yhatTBATS_wk_yr_final_ac,y),Acc(yhatTBATS_wk_yr_final_mc,y),
                    Acc(yhatTBATS_wk_yr_final_at,y),Acc(yhatTBATS_wk_yr_final_mt,y),
                    Acc(yhat_Reg_TBATS_final_ac,y),Acc(yhat_Reg_TBATS_final_mc,y),
                    Acc(yhat_Reg_TBATS_final_at,y),Acc(yhat_Reg_TBATS_final_mt,y),
                    Acc(yhat_Reg_ARIMA505_final_ac,y),Acc(yhat_Reg_ARIMA505_final_mc,y),
                    Acc(yhat_Reg_ARIMA505_final_at,y),Acc(yhat_Reg_ARIMA505_final_mt,y),
                    Acc(yhat_Reg_ARIMA301_final_ac,y),Acc(yhat_Reg_ARIMA301_final_mc,y),
                    Acc(yhat_Reg_ARIMA301_final_at,y),Acc(yhat_Reg_ARIMA301_final_mt,y),
                    Acc(yhat_ARIMAxreg_final_ac,y),Acc(yhat_ARIMAxreg_final_mc,y),
                    Acc(yhat_ARIMAxreg_final_at,y),Acc(yhat_ARIMAxreg_final_mt,y))

#------------------------------------------------------------
xactual <- d[1:n2,1] 

errMatModel <- rbind(Acc(xhat_Reg_fullac,xactual),Acc(xhat_Reg_fullmc,xactual),
                     Acc(xhat_Reg_final_at,xactual),Acc(xhat_Reg_final_mt,xactual),
                     Acc(xhat_TBATSwk_yr_fullac,xactual),Acc(xhat_TBATSwk_yr_fullmc,xactual),
                     Acc(xhat_TBATSwk_yr_final_at,xactual),Acc(xhat_TBATSwk_yr_final_mt,xactual),
                     Acc(xhat_Reg_TBATS_fullac,xactual),Acc(xhat_Reg_TBATS_fullmc,xactual),
                     Acc(xhat_Reg_TBATS_final_at,xactual),Acc(xhat_Reg_TBATS_final_mt,xactual),
                     Acc(xhat_Reg_ARIMA505_fullac,xactual),Acc(xhat_Reg_ARIMA505_fullmc,xactual),
                     Acc(xhat_Reg_ARIMA505_final_at,xactual),Acc(xhat_Reg_ARIMA505_final_mt,xactual),
                     Acc(xhat_Reg_ARIMA301_fullac,xactual),Acc(xhat_Reg_ARIMA301_fullmc,xactual),
                     Acc(xhat_Reg_ARIMA301_final_at,xactual),Acc(xhat_Reg_ARIMA301_final_mt,xactual),
                     Acc(xhat_ARIMAxreg_fullac,xactual),Acc(xhat_ARIMAxreg_fullmc,xactual),
                     Acc(xhat_ARIMAxreg_final_at,xactual),Acc(xhat_ARIMAxreg_final_mt,xactual))

nrParReg <- length(fitF2minus$coefficients) + 1 
nrParTBATSwkYr <- length(fitTBATSwk_yr$parameters$vect) + 1 
nrParRegTBATS <- length(fit_Reg_TBATS$parameters$vect) + 1
nrParARIMA505 <- length(fit_Reg_ARIMA505$coef) + 1
nrParARIMA301 <- length(fit_Reg_ARIMA301$coef) + 1
nrParARIMAxreg <- length(fit_ARIMA_xreg$coef) + 1

nrParClassSeas <- 7*freq1 + 7*24
nrParTrigSeas <- 9

nrParVec <- cbind(nrParReg + nrParClassSeas, nrParReg + nrParClassSeas, 
                  nrParReg + nrParTrigSeas, nrParReg + nrParTrigSeas,
                  nrParTBATSwkYr + nrParClassSeas, nrParTBATSwkYr + nrParClassSeas,
                  nrParTBATSwkYr + nrParTrigSeas, nrParTBATSwkYr + nrParTrigSeas,
                  nrParRegTBATS + nrParClassSeas, nrParRegTBATS + nrParClassSeas, 
                  nrParRegTBATS + nrParTrigSeas, nrParRegTBATS + nrParTrigSeas,
                  nrParARIMA505 + nrParClassSeas, nrParARIMA505 + nrParClassSeas,
                  nrParARIMA505 + nrParTrigSeas, nrParARIMA505 + nrParTrigSeas,
                  nrParARIMA301 + nrParClassSeas, nrParARIMA301 + nrParClassSeas, 
                  nrParARIMA301 + nrParTrigSeas, nrParARIMA301 + nrParTrigSeas,
                  nrParARIMAxreg + nrParClassSeas, nrParARIMAxreg + nrParClassSeas,
                  nrParARIMAxreg + nrParTrigSeas, nrParARIMAxreg + nrParTrigSeas)
              
psiModel <- n2*log(as.numeric(errMatModel[,1])) + 2*nrParVec
psiHold <- forecastLen*log(as.numeric(errMatHold[,1])) + 2*nrParVec
psiOutput <- cbind(t(psiModel), t(psiHold))
colnames(psiOutput) <- c('Model','Hold out')


errorMatColNames <- c("Regression; Additive; Classical","Regression; Multiplicative; Classical",
                      "Regression; Additive; Polynomial","Regression; Multiplicative; Polynomial",
                      "TBATS with weekly, annual seas.; Additive; Classical","TBATS with weekly, annual seas.; Multiplicative; Classical",
                      "TBATS with weekly, annual seas.; Additive; Polynomial","TBATS with weekly, annual seas.; Multiplicative; Polynomial",
                      "Reg. and TBATS on residuals; Additive; Classical","Reg. and TBATS on residuals; Multiplicative; Classical",
                      "Reg. and TBATS on residuals; Additive; Polynomial","Reg. and TBATS on residuals; Multiplicative; Polynomial",
                      "Reg. and ARIMA (5,0,5) on resid.; Additive; Classical","Reg and ARIMA (5,0,5) on resid.; Multiplicative; Classical",
                      "Reg. and ARIMA (5,0,5) on resid.; Additive; Polynomial","Reg. and ARIMA (5,0,5) on resid.; Multiplicative; Polynomial",
                      "Reg. and ARIMA on res (3,0,1); Additive; Classical","Reg. and ARIMA on res (3,0,1); Multiplicative; Classical",
                      "Reg. and ARIMA on res (3,0,1); Additive; Polynomial","Reg + ARIMA on res (3,0,1); Multiplicative; Polynomial",
                      "ARIMA with xreg; Additive; Classical","ARIMA with xreg; Multiplicative; Classical",
                      "ARIMA with xreg; Additive; Polynomial","ARIMA with xreg; Multiplicative; Polynomial")

rownames(errMatHold) <- errorMatColNames
#write.csv(errMatHold,"Table4_1.csv")

rownames(errMatModel) <- errorMatColNames
#write.csv(errMatModel,"Table3_1.csv")
##write.csv(cbind(errMatModel,errMatHold), "Table_3_4_rev2.csv")

rownames(psiOutput) <- errorMatColNames
write.csv(psiOutput,"psiTables3_4.csv")

if (0) {

#---------------------------------------------------------------------------------------------------------
# FORECAST COMBINATIONS
#---------------------------------------------------------------------------------------------------------
library(ForecastCombinations)

fhat <- cbind(xhat_Reg_fullmc, xhat_Reg_TBATS_fullmc, xhat_TBATSwk_yr_fullmc, xhat_Reg_ARIMA505_fullmc,
              xhat_ARIMAxreg_fullac)
fhatNew <- cbind(yhat_Reg_final_mc, yhat_Reg_TBATS_final_mc, yhatTBATS_wk_yr_final_mc , yhat_Reg_ARIMA505_final_mc,
              yhat_ARIMAxreg_final_mc)
fhat <- cbind(xhat_Reg_fullmc, xhat_Reg_TBATS_fullmc,xhat_TBATSwk_yr_fullmc)
fhatNew <- cbind(yhat_Reg_final_mc, yhat_Reg_TBATS_final_mc,yhatTBATS_wk_yr_final_mc)
comb_simple <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'simple')
comb_ols <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'ols')
comb_var <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'variance based')
comb_robust <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'robust')

comb_err <- rbind(cbind(Acc(comb_simple$fitted,xactual), Acc(comb_simple$pred,y)),
                  cbind(Acc(comb_ols$fitted,xactual), Acc(comb_ols$pred,y)),
                  cbind(Acc(comb_var$fitted,xactual), Acc(comb_var$pred,y)),
                  cbind(Acc(comb_robust$fitted,xactual), Acc(comb_robust$pred,y)))
rownames(comb_err) <- c("Simple","OLS","Variance based","RObust Reg")

fhat <- cbind(xhat_Reg_final_mt, xhat_Reg_TBATS_final_mt, xhat_TBATSwk_yr_final_mt, xhat_Reg_ARIMA505_final_mt,
              xhat_ARIMAxreg_final_mt)
fhatNew <- cbind(yhat_Reg_final_mt, yhat_Reg_TBATS_final_mt, yhatTBATS_wk_yr_final_mt, yhat_Reg_ARIMA505_final_mt,
                 yhat_ARIMAxreg_final_mt)
comb_simple <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'simple')
comb_ols <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'ols')
comb_var <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'variance based')
comb_robust <- Forecast_comb(xactual, fhat, fhat_new = fhatNew, Averaging_scheme = 'robust')


comb_err_trig <- rbind(cbind(Acc(comb_simple$fitted,xactual), Acc(comb_simple$pred,y)),
                  cbind(Acc(comb_ols$fitted,xactual), Acc(comb_ols$pred,y)),
                  cbind(Acc(comb_var$fitted,xactual), Acc(comb_var$pred,y)),
                  cbind(Acc(comb_robust$fitted,xactual), Acc(comb_robust$pred,y)))
rownames(comb_err_trig) <- c("Simple","OLS","Variance based","Robust Reg")

print(rbind(comb_err, comb_err_trig))
write.csv(rbind(comb_err, comb_err_trig), "forecastCombErrors.csv")
}






# #------------SEasonality Plot--------------------------
# freq <- 288 #for daily average
# meanBy <- as.vector(t(matrix(rep(c(1:(nrow(d)/freq)),freq),
#                              nrow = nrow(d)/freq)))
# davg <- aggregate(d$load,by = list(meanBy),mean)$x
# davg_rep <- as.vector(matrix(rep(davg,freq),nrow=freq,byrow=TRUE)) 
# dMult <- d$load/davg_rep
# par(mai=c(.82,0.82,0.32,0.42))
# stIdx <- which(d$date == "02-Jan-12")[1]
# plot(dMult[(1:(7*288))+stIdx-1],type='l',xlab='Jan 2, 2012 - Jan 8, 2012',
#      ylab='Load / Average Daily Load', lwd=2,xaxt='n')
# seasClass <- siMult$block[modelBlkSeasAppl]*siMult$hourly[modelHrSeasAppl]
# lines(as.vector(seasClass[(1:(7*288))+(3*288)]),col='red',lty = 2,lwd=2)
# lines(as.vector(seas$mult[,c(2:7,1)]),col='blue',lty=3,lwd=2)
# legend('bottomright',ncol =3,c('Data','Classical','Poly'),col=c('black','red','blue'),lwd = c(2,2,2),lty=c(1,2,3))
# axis(1,at=seq(144,2000,288),labels = c('Mon','Tue','Wed','Thu','Fri','Sat','Sun'))
