#Empirical method HW1 Problem 2
#part 1
library(xts)
Fund <- read.xls("DBV.xlsx",  na.strings = FALSE )
SP <- read.xls("GSPC.xlsx",  na.strings = FALSE )


Fund[,1]<-as.Date(Fund[,1],format("%m/%d/%Y"))
#create xts object
Funddata<-as.xts(Fund[,7],order.by=Fund[,1])
Fundret<-diff(log(Funddata),lag=1)
Fundret<-Fundret[-1,]
head(Fundret)

SP[,1]<-as.Date(SP[,1],format("%m/%d/%Y"))
SPdata<-as.xts(SP[,7],order.by=SP[,1])
SPret<-diff(log(SPdata),lag=1)
SPret<-SPret[-1,]

#time-series plot of daily log-returns
plot(Fundret,main="Time series plot of DBV log return")
plot(SPret,main="Time series plot of GSPC log return")

hist(Fundret,main="DBV log return histogram",breaks=40)
hist(SPret,main="GSPC log return histogram",breaks=40)

#part 2
library(fBasics)
basicStats(Fundret)
t3<-function(ret){
  S3<-skewness(ret)
  T<-length(ret)
  t3<-S3/sqrt(6/T)
  t3
}

t4<-function(ret){
  K4<-kurtosis(ret)
  T<-length(ret)
  t4<-K4/(sqrt(24/T))
  t4
}

t3(Fundret)
t4(Fundret)
normalTest(Fundret,method='jb')

criticalVal=abs(qnorm(0.05/2))
#abs(t3(Fundret))>criticalVal
#reject the null hypothesis, log return of DBV is negatively skewed at 5% significant level
#t4(Fundret)>criticalVal
#reject the null hypothesis, log return of DBV has heavy tails at 5% significant level
#Jarque-Bera test Asymptotic p Value: < 2.2e-16 
#reject the null hypothesis at 5% significant level

basicStats(SPret)
t3(SPret)
t4(SPret)
normalTest(SPret,method='jb')
normalTest(SPret,method='jb')@test$p.value

#abs(t3(SPret))>criticalVal
#reject the null hypothesis, log return of GSPC is negatively skewed at 5% significant level
#t4(Fundret)>criticalVal
#reject the null hypothesis, log return of GSPC has heavy tails at 5% significant level
#Jarque-Bera test Asymptotic p Value: < 2.2e-16 
#reject the null hypothesis at 5% significant level


#part 3
stats<-data.frame(cbind(basicStats((Fundret)),basicStats((SPret))))
table<-rbind(stats["Mean",],stats["Variance",],stats["Stdev",],stats["Skewness",],stats["Kurtosis",])
table["skewness test Tvalue",]<-cbind(t3(Fundret),t3(SPret))
table["kurosis test Tvalue",]<-cbind(t4(Fundret),t4(SPret))
table["Jarque-Bera test Pvalue",]<-cbind(normalTest(Fundret,method='jb')@test$p.value,normalTest(SPret,method='jb')@test$p.value)
colnames(table)=c("DBV","GSPC")
table

