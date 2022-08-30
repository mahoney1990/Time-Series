library(vars)
library(mFilter)
library(tseries)
library(TSstudio)
library(forecast)
library(tidyverse)
library(readxl)
agg=FALSE

#Time series VAR models. Key assumption: data is stationary. We'll have some tests below.

#Load Federal Reserve data
setwd("C:/Users/mahon/Tapering Project")

data=read_xlsx("Rule Data.xlsx")

#Read date variable as POSIX date
data$Date=as.POSIXct(data$Date)

#Index of observations between 2007 and 2020
idx_lower=which(data$Date>"2007-1-1")
idx_upper=which(data$Date<"2020-1-1")
idx=intersect(idx_lower,idx_upper)

#Subset data
data=data[idx,]
N=length(data$Date)
int_index=seq(from=3,to=N,by=4)

#Generate first difference of output gap variable (i.e. diff GDP and GDP*)
for(i in 1:(length(int_index)-1)){
 l=data$OutputGap[int_index[i]]
 u=data$OutputGap[int_index[i+1]]
 int_values=approx(c(l,u),n=5)$y[2:4]
 data$OutputGap[(int_index[i]+1):(int_index[i+1]-1)]=int_values
}

#Option one: quarterly aggregation
#Aggregate to quarterly averages -- Fed funds rate, quanitative easing variable, output gap, inflation
if(agg==TRUE){
FFR=ts(aggregate(data$FFR, by=list(qtr=data$Index), FUN="mean"))[,2]
qe=ts(aggregate(data$Monthly_Secs, by=list(qtr=data$Index), FUN="mean"))[,2]
X_hat=ts(aggregate(data$OutputGap, by=list(qtr=data$Index), FUN="mean"))[,2]
pi=ts(aggregate(data$PCE, by=list(qtr=data$Index), FUN="mean"))[,2]

#Convert to first differences on key variables
FFR=ts(diff(FFR,frequency = 4))
qe=ts(diff(qe,frequency = 4))
X_hat=ts(X_hat,frequency = 4)
pi=ts(pi,frequency = 4)}

#Option Two: Monthly observations
if(agg==FALSE){
FFR=ts(diff((data$FFR),frequency = 12))
qe=ts(diff(data$Monthly_Secs,frequency = 12))
X_hat=ts(data$OutputGap,frequency = 12)
pi=ts(data$PCE,frequency = 12)}

X_hat=ts(diff(X_hat,frequency = 12))

#Okay now, lets smooth out inflation data with a Hodrick-Prescott filter
#Then we can form some plots

hp<-hpfilter(pi,(1000))
pi=hp$cycle

#Form some time-series plots
ts.plot(pi) #inflation
ts.plot(FFR) #Fed funds rate
ts.plot(qe) #Fed bond portfolio size
ts.plot(X_hat) #Output gap
 


#***This is our test for a unit root+*** We could alos use a DF test
#If the p-value of this test is low (less than .10), we reject the null of
# non stationary data

pp.test(pi)
pp.test(FFR)
pp.test(qe)
pp.test(X_hat)


#Okay, now that we've established the stationary of our data, lets
#Build a VAR forecasting model. In this first model, the only two
#time-series we use are inflation and the output gap
n=length(FFR)

#MODEL 1 -- two variable VAR
#VAR requires our data to be formatted. Lets build a new data frame to talk to the VAR commmand
v1 <- cbind(FFR[-1], qe[-1])
colnames(v1) <- cbind("FFR","qe")
pi=pi[-c(n-1,n)]
X_hat=X_hat[-n]

exo <- cbind(X_hat, pi)
colnames(exo)<-cbind('X_hat','pi')

#Lag select command will choose the best lag order for our VAR (between 1 and 10 lags)
lagselect <- VARselect(v1, lag.max = 10, type = "const")
lagselect$selection

#Now we put it all together. Here we estimate the VAR at optimal lag order by OLS
Model1 <- VAR(v1, p = 2, type = "const", season = NULL) 
summary(Model1)

#Lets project out a year
fore=predict(Model1,n.ahead = 12, ci = 0.95)

#Forecast shows slowing FFR and a declining bond portfolio -- which was correct before COVID...
ts.plot(fore$fcst$FFR[,1:3])
ts.plot(fore$fcst$qe[,1:3])

#MODEL 2 -- 4 variable VAR
#Lets add inflation and output gap to the VAR and build another model
v2 <- cbind(FFR[-1], qe[-1],pi[-1],X_hat[-1])
colnames(v2) <- cbind("FFR","qe","inflation","output")

#Lag select command will choose the best lag order for our VAR (between 1 and 10 lags)
lagselect <- VARselect(v1, lag.max = 10, type = "const")
lagselect$selection

#Now we put it all together. Here we estimate the VAR at optimal lag order by OLS
Model2 <- VAR(v2, p = 2, type = "const", season = NULL) 
summary(Mode2)

#Lets project out a year
fore=predict(Model2,n.ahead = 12, ci = 0.95)

ts.plot(fore$fcst$FFR[,1:3])
ts.plot(fore$fcst$qe[,1:3])
ts.plot(fore$fcst$pi[,1:3])
ts.plot(fore$fcst$X_hat[,1:3])

#MODEL 3 -- 2 variable VAR with level realtionship on inflation and output variables

data=read_xlsx("Rule Data.xlsx")
data$Date=as.POSIXct(data$Date)

idx_lower=which(data$Date>"2020-1-1")
idx_upper=which(data$Date<"2022-1-1")
idx=intersect(idx_lower,idx_upper)
data=data[idx,]

hp<-hpfilter(data$PCE,(1000))
pi_fut=hp$cycle
ts.plot(pi_fut)
pp.test(pi_fut)
pi_fut=ts(pi_fut[-22])

N=length(data$Date)
int_index=seq(from=3,to=N,by=4)

for(i in 1:(length(int_index)-1)){
  l=data$OutputGap[int_index[i]]
  u=data$OutputGap[int_index[i+1]]
  int_values=approx(c(l,u),n=5)$y[2:4]
  data$OutputGap[(int_index[i]+1):(int_index[i+1]-1)]=int_values
}

X_hat_fut=ts(diff(data$OutputGap))

exog_future=cbind(X_hat_fut,pi_fut)
colnames(exog_future)=cbind('X_hat','pi')


#Miscellaneous Tests
#Want the residuals to exhibit no autocorrelation
# 
Serial1 <- serial.test(Model1, lags.pt = 5, type = "PT.asymptotic")
Serial1

#Hetroskdasticity
Arch1 <- arch.test(Model1, lags.multi = 5, multivariate.only = TRUE)
Arch1

#Cointegration
library(aTSA)
coint.test(data$FFR,data$Monthly_Secs)

#VECM

library(urca) # Load package

# Estimate
vec <- ca.jo(v1, ecdet = "none", type = "trace",
             K = 2, spec = "transitory", season = 12)

summary(vec)


coin_reg=glm(data$Monthly_Secs~data$FFR)
summary(coin_reg)

beta_reg<-glm(lead_beta ~ Phase+State+time,data=panel)



