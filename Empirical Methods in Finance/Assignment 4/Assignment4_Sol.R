#Empirical Method HW4 
setwd("~/Acads/MFE 2016/2nd Quarter/Empirical Methods in Finance/HW4")

install.packages("vrtest")
install.packages("fGarch")

library(lubridate)
library(xts)
library(ggplot2)
library(dplyr)
library(vrtest)
library(fUnitRoots)
library(zoo)
library(fGarch)
library(xlsx)
# question 1 part b
# get monthly exchange rate from 1970 to 2008
usdi<-read.xlsx("dollar_broadindex_tradeweighted.xlsx",'dollar_broadindex_tradeweighted')
usdi$Date<-as.Date(usdi$Date,format("%m/%d/%Y"))
usdi<-usdi[which(usdi$Date=="1971-01-31"):which(usdi$Date=="2008-11-30"),]
usdi<-select(usdi,Date,Close)
usdi$log<-log(usdi$Close)
#calculate and plot variance ratio for log exchange rate
# Changes made from row 24 to 29, to replace the two plots before
sigma_sq=var(diff(usdi$log))
par(mfrow=c(1,1)) 
vr_usdi<- sapply(1:100,function(x){
  var(diff(usdi$log,lag=x))/(x*sigma_sq)
}) 
plot(vr_usdi,main="Log Exchange Rate Variance Ratio", xlab="", ylab="") 

# variance ratio approximate 0 as lags goes larger, so the series may not follow a random walk, 
# but a time trend.

# question 1 part c
y <- usdi$log 
nob <- length(y) 
r <- diff(y,lag=1) 
Auto.VR(r) 

kvec <- c(2,5,10,30,50,100,150) 
VR.plot(r,kvec) 
# variance ratio increase above 1 first, and decrease below 1 as holding period gets gets longer.
# It also shows that the series is not unit root.

# question 2 part a
gdpq<-read.csv("fredgraph.csv")
colnames (gdpq) <-c("Date","GDP")
gdpq$Date<-as.yearqtr(gdpq$Date,format = "%Y-%m-%d")
gdpq$Date <-as.Date(gdpq$Date, format="%Y%q")
gdpq$loggdp<-log(gdpq$GDP)
lggdpq<-gdpq[,c("Date","loggdp")]
# get difference for quarterly gdp
diff_lggdp<-diff(gdpq$loggdp)
diff_lggdp<-data.frame(gdpq$Date[2:length(gdpq$Date)],diff_lggdp)
colnames (diff_lggdp)<-c("Date","diff_loggdp")
# plot Log GDP and its changes
p1<-ggplot(lggdpq)+geom_line(aes(lggdpq$Date,lggdpq$loggdp)) + xlab("Quarter") + ylab("Log GDP Per Capita") + ggtitle("Log GDP Per Capita (Quarterly)")
p1
p2<-ggplot(diff_lggdp)+geom_line(aes(diff_lggdp$Date,diff_lggdp$diff_loggdp)) + xlab("Quarter") + ylab("Difference of Log GDP Per Capita") + ggtitle("Difference of Log GDP Per Capita (Quarterly)")
p2

par(mfrow=c(2,1))
acf(lggdpq$loggdp)
pacf(lggdpq$loggdp)

par(mfrow=c(2,1))
acf(diff_lggdp$diff_loggdp)
pacf(diff_lggdp$diff_loggdp)

m1=ar(diff_lggdp$diff_loggdp,method = 'mle')
m1$order

# The plot of log GDP shows an upward trend. The first difference of log GDP seems to vary around a 
# fixed mean leve, although the variabilty appears to be smaller in recent years.
# AR model of difference series is significant at lag 3, and the sample PACF of the differenced 
# series shown significant lag of 12. Therefore, Lag =13 that will take this into account in the 
# augmented Dicky-Fuller unit root test.

# With growth of GDP over time, if the process has a unit root, it is expected to have a drift.
# If it is trend stationary, it is reasonable to expect a time trend.

# question 2 part b

# Therefore, we will do two augmented Dicky-Fuller test, one with a constant/drift, 
# and the other one with a constant and a time trend.

# ADF test with constant
# H0: phi_1=1 H1:|phi_1|<1
# null model: loggdp_t=phi_0+loggdp_t-1+error_t
# alternative model: loggdp_t=phi_0+phi_1*loggdp_t-1+error_t

# ADF test with a time trend and a constant
# H0: phi_1=1 H1:|phi_1|<1
# null model: loggdp_t = phi_0 + loggdp_t-1 + omega*t + error_t
# alternative model: loggdp_t = phi_0 + phi_1*loggdp_t-1 + omega*t + error_t

# question 2 part c
# ADF test with constant
adfTest (lggdpq$loggdp,lags=13,type=c("c"))
# ADF test with a time trend and a constant
adfTest (lggdpq$loggdp,lags=13,type=c("ct"))

# Both test fail to reject the null hypothesis of unit root in log GDP. It indicates that log GDP 
# is non-stationary process. 

# question 3 part a
data <-read.xlsx("currency_fund_prices.xlsx","currency_fund_prices")
data$Date <- as.Date(data$Date,format("%m/%d/%Y"))
# calculate daily log return
ccy<-select(data,Date,Adj.Close)
ccy$ret[2:length(ccy$Adj.Close)] <-diff(log(ccy$Adj.Close))
rets <- select(ccy,Date,ret)
rets <- rets[2:length(ccy$ret),]
par(mfrow=c(1,1))
plot(rets,type="l",xlab = 'Year', ylab = 'ln-return', main ='Time plot of daily log return')
par(mfrow=c(2,1))
acf(rets$ret, lag=30)
pacf(rets$ret, lag=30)

# seems log return is a constant plus innovation. abstract conditional mean from variation
rets$innovation <-rets$ret-mean(rets$ret)
Box.test(rets$innovation^2,lag=30,type='Ljung')
# LB Q-test rejects the null of no autocorrelation in volatility
par(mfrow=c(2,1))
acf(rets$innovation^2, lag=30)
pacf(rets$innovation^2, lag=30)

# innovation squared shows significant PACF for 1 to 12 lags as well as lag 18 and 30.
# one can employ more parsimonious GARCH model
# Try GARCH (1,1) with Gaussian innovations
m1 <- garchFit(~1+garch(1,1),data=rets$ret,trace=F)
summary(m1)
# GARCH(1,1) model: r_t = 1.521e-04 + a_t, a_t = sigma_t * e_t
# sigma_t^2 = 5.952e-07 + 1.377e-01 * a_(t-1)^2 + 8.664e-01 * sigma_(t-1)^2
# AIC = -7.161273
v1 = volatility(m1) # obtain volatility
resi <- residuals(m1,standardize=T) # standard residuals
vol <- data.frame(rets$Date,v1)
res <- data.frame(rets$Date,resi)
par(mfrow=c(2,1)) # show volatility and residuals
plot (vol,xlab='year',ylab='volatility', type='l')
plot (res,xlab='year',ylab='st. resi', type='l')
par(mfrow=c(2,2)) # obtain ACF & PACF
acf(resi, lag=24)
pacf(resi, lag=24)
acf(resi^2, lag=24)
pacf(resi^2, lag=24)
# no significant ACF/PACF for GARCH residuals
# model is good fit for the data

# question 3 part b
# prediction of log return level for next 20 periods
pred <- predict(m1,20) 
# log return volatility over 20 days, assuming iid daily log return, variance over 20 day period 
# equals sum of standard error squared on each day.
std <- sqrt(sum(pred$standardDeviation^2))
# the 20-trading-day return volatility (standard deviation) is 0.03525

# question 3 part c
# VaR at 5% on $2bn long position in currency is $115.95 million,
# m1 model predict constant mean is 1.521e-04, but it is not significant different from 0, given
# its p-value is 0.168>.05. So we assume mean of 20-day log return approximates 0.
2000*abs(qnorm(.05))*std


