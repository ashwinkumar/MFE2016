# MGMT237E HW6
# Question1

library(lubridate)
library(xts)
library(dplyr)

##Import Data
#Read 48 Industry Data/3 Factor Model Data
ind48<-read.csv("48_Industry_Portfolios.CSV",header=T)
fac3<-read.csv("F-F_Research_Data_Factors.CSV",header=T)
#Rename Date Column
colnames(ind48)[1]<-"Date"
colnames(fac3)[1]<-"Date"
#Reformat Date Column
ind48$Date<-parse_date_time(ind48$Date,"%y%m")
fac3$Date<-parse_date_time(fac3$Date,"%y%m")
ind48$Date<-as.Date(ind48$Date)
fac3$Date<-as.Date(fac3$Date)
#Select data from 1960 to 2015
ind48<-ind48[ind48$Date>=as.Date("1960-01-01"),]
fac3<-fac3[fac3$Date>=as.Date("1960-01-01"),]
fac3<-fac3[fac3$Date<as.Date("2015-12-31"),]

#Redefine NA convention
ind48[ind48==-99.99]=NA
#Remove columns with NA
ind48<-ind48[,colSums(is.na(ind48))==0]

#Excess return matrix: subrtract risk free rate from portfolio returns
xsret<-ind48[,2:length(ind48)]-fac3$RF

# number of periods
T=length(fac3$Date)
# number of portfolios/industries
N=dim(ind48)[2]-1
# number of factors
K=3 

# run the time series regression
beta=matrix(0,K+1,N)
predxsret=matrix(0,1,N)


# X is T*4 maxtrix for factors scale overtime
# beta is a 4*N matrix: factor loading for each industry with first row constant
X=cbind(1,fac3[,2:4])
X=as.matrix(X)
for (i in 1:N){
  out=lm(xsret[,i]~fac3$Mkt.RF+fac3$SMB+fac3$HML)
  beta[,i]=out$coefficients
}
beta=as.matrix(beta)
pred=X%*%beta
#regression intercept are the pricing errors.
alpha=as.matrix(beta[1,])

# difference between actual excess return and predicted excess return
error=xsret-pred

# plot predicted mean excess return vs. realized mean excess return
Min=min(colMeans(xsret))
Max=max(colMeans(xsret))
plot((colMeans(pred)-beta[1,])~colMeans(xsret),xlim=range(Min,Max),ylim=range(Min,Max),xlab="actual mean excess return", ylab="predicted mean excess return", main="FF Model on industry portfolios", col="blue")
fit=lm(colMeans(xsret)~(colMeans(pred)-beta[1,]))
abline(fit$coef,col="red",lwd=2)
summary(fit)

# covariance matrix across industry portfolios
sigma=cov(error)
facmean=as.matrix(colMeans(fac3[,2:4]))
facsigma=cov(fac3[,2:4])
# calculate F-statistic for pricing error, which follow F(N,T-N-K))
Fstat=(T-N-K)/N*(1+t(facmean)%*%facsigma%*%facmean)^(-1)*(t(alpha)%*%chol2inv(chol(sigma))%*%alpha)
pf(Fstat,df1=N,df2=T-N-K)
# Reject the null, pricing error are jointly deviated from zero.


