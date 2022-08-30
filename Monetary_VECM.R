rm(list=ls()) #Removes all items in Environment!
library(tseries) # for `adf.test()`
library(dynlm) #for function `dynlm()`
library(vars) # for function `VAR()`
library(nlWaldTest) # for the `nlWaldtest()` function
library(lmtest) #for `coeftest()` and `bptest()`.
library(broom) #for `glance(`) and `tidy()`
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(knitr) #for `kable()`
library(forecast)
library(readxl)
library(mFilter)

#Load in federal reserve data
setwd("C:/Users/mahon/Tapering Project")
data=read_xlsx("Rule Data.xlsx")

#Convert to date format
data$Date=as.POSIXct(data$Date)

#Subset data
idx_lower=which(data$Date>"2008-1-1")
idx_upper=which(data$Date<"2020-1-1")

idx=intersect(idx_lower,idx_upper)
data=data[idx,]

N=length(data$Date)
int_index=seq(from=3,to=N,by=4)

#Build first differnces to make our data stationary -- basically none of this stuff
#works if our data is non stationary
for(i in 1:(length(int_index)-1)){
  l=data$OutputGap[int_index[i]]
  u=data$OutputGap[int_index[i+1]]
  int_values=approx(c(l,u),n=5)$y[2:4]
  data$OutputGap[(int_index[i]+1):(int_index[i+1]-1)]=int_values
}


#Cointegration -- basically, we need to determine if there is a relationship in levels
#between our time series variables. We can then account for this levels relation by
#using a Vector Error Correction model (VECM).

coin_reg=glm(data$Monthly_Secs~data$FFR)
summary(coin_reg)

cint1.dyn <- dynlm(FFR~Monthly_Secs, data=data)
kable(tidy(cint1.dyn), digits=3,
      caption="The results of the cointegration equation 'cint1.dyn'")

#Test for Cointegration -- 
ehat <- resid(cint1.dyn)
diff_ehat=diff(ehat)
lag_ehat=lag(ehat)[-1]

cint2.dyn <- dynlm(diff_ehat~lag_ehat-1)
summary(cint2.dyn)


#Lets generate some additional time series variables from the data frame
X_hat=ts(data$OutputGap,frequency = 12)
pi=ts(data$PCE,frequency = 12)

X_hat=ts(diff(X_hat,frequency = 12))

#Smooth out the inflation time series
hp<-hpfilter(pi,(1000))
pi=hp$cycle

pi.dyn <- dynlm(PCE~OutputGap, data=data)
summary(pi.dyn)


#Sure looks like there might be a level type relationship here. So lets make that VECM
#Estimate VECM

vec_FFR<- dynlm(diff(FFR)~lag(ehat)[-1]+pi[-N]+X_hat[-N], data=data)
vec_sec <- dynlm(diff(Monthly_Secs)~lag(ehat)[-1]+pi[-N]+X_hat[-N], data=data)

summary(vec_FFR)
summary(vec_sec)

ts.plot(diff(data$FFR))
ts.plot(predict(vec_FFR))

ts.plot(diff(data$Monthly_Secs))
ts.plot(predict(vec_sec))

data$Monthly_Secs-predict(vec_sec)

#####Predict into 2020

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

#Generate exogenous dataframe. We allow the qe and FFR variables to depend on these
#variables in levels. I.e. we believe there is a defined, equilbrium relationship between t
#the Feds decision variables and economic conditions.

exog_future=cbind(X_hat_fut,pi_fut)
colnames(exog_future)=cbind('X_hat','pi')

#Generate bivariate time-series
v1 <- cbind(data$FFR[-1], data$Monthly_Secs[-1])
colnames(v1) <- cbind("FFR","qe")
v1=data.frame(v1)

#Generate forecasts
fore_FFR=predict(vec_FFR,newdata=v1,dumvar = exog_future,n.ahead = 12, ci = 0.95)
fore=predict(vec_sec,newdata=v1,dumvar = exog_future,n.ahead = 12, ci = 0.95)

ehat=predict(coin_reg,v1)

ehat=(v1$qe-(coin_reg$coefficients[1]+v1$FFR*coin_reg$coefficients[2]))

int_ffr=vec_FFR$coefficients[1]
alpha_ehat=vec_FFR$coefficients[2]
alpha_pi=vec_FFR$coefficients[3]
alpha_xhat=vec_FFR$coefficients[4]

delta_FFR=int_ffr+alpha_ehat*ehat[1]+alpha_pi*pi[1]+alpha_xhat*xhat[1]


int_sec=vec_sec$coefficients[1]
phi_ehat=vec_sec$coefficients[2]
phi_pi=vec_sec$coefficients[3]
phi_xhat=vec_sec$coefficients[4]

delta_sec=int_sec+phi_ehat*ehat[1]+phi_pi*pi[1]+phi_xhat*xhat[1]




