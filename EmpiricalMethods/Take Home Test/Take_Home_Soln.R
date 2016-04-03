setwd("~/Acads/MFE 2016/2nd Quarter/Empirical Methods in Finance/Week 5")
library(lubridate)
library(xlsx)
library(xts)
library(ggplot2)
data_level <- read.xlsx('open_book_data.xlsx','Sheet1')
t <- 1:length(data_level[,1])
colnames(data_level)<- "Levels"
reg <- lm(data_level[,1]~t)
c1 <- reg$coefficients[1]
m1 <- reg$coefficients[2]
logdata <- log(data_level)
reg2 <- lm(logdata[,1]~t)
c2 <- reg2$coefficients[1]
m2 <- reg2$coefficients[2]
y_diff <- diff(logdata[,1], lag=1)
diffdata <- data.frame(y_diff)
colnames(diffdata) <- "DiffLog"
t2 <- 1:length(diffdata[,1])
p1 <- ggplot(data_level,aes(x=t, y=data_level[,1])) + geom_line(color= '#9999FF', size=1.5) + ggtitle("Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue")) 
p1 + geom_abline(intercept = c1 , slope = m1)
par(mfrow =c(1,2))
acf(data_level,lag.max = 12)
pacf(data_level,lag.max = 12)


p2 <- ggplot(logdata,aes(x=t, y=logdata[,1])) + geom_line(color= '#9999FF', size=1.5) + ggtitle("Log Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue")) 
p2+ geom_abline(intercept = c2 , slope = m2)
par(mfrow = c(1,2))
acf(logdata,lag.max = 12)
pacf(logdata)

ggplot(diffdata,aes(x=t2, y=diffdata[,1])) + geom_line(color= '#9999FF', size=1.5) + ggtitle("Diff Log Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue")) 
par(mfrow = c(1,2))
acf(diffdata, lag.max = 12)
pacf(diffdata)

y1 <- diff(data_level[,1],lag=12)
t3 <- 1:length(y1)
diff_12 <- data.frame(y1)

ggplot(diff_12,aes(x=t3, y=diff_12[,1])) + geom_line(color= '#9999FF', size=1.5) + ggtitle("12th Diff Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue")) 
y<- diff(logdata[,1],lag=12)
logdiff_12 <- data.frame(y)
ggplot(logdiff_12,aes(x=t3, y=logdiff_12[,1])) + geom_line(color= '#9999FF', size=1.5) + ggtitle("12th Log Diff Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue")) 
par(mfrow = c(1,2))
acf(logdiff_12)
pacf(logdiff_12)

myModel1<-arima(logdiff_12, order=c(2,0,0))
myModel1Coef<-coef(myModel1)
myModel1Res<-residuals(myModel1)
par(mfrow = c(1,2))
acf(myModel1Res,lag.max = 12)
pacf(myModel1Res,lag.max = 12)

predition_t <- 24
y_predict <- rep(0,predition_t+2)
y_predict[1] <- logdiff_12[length(t3)-1,1] 
y_predict[2] <- logdiff_12[length(t3),1] 
for(i in 1:24){
    f <- y_predict[i+1]*as.numeric(myModel1Coef[1]) + y_predict[i]*as.numeric(myModel1Coef[2]) + as.numeric(myModel1Coef[3])
    y_predict[i+2] <- f
}
t0 <- length(logdata[,1])
logdata_prediction <- logdata[,1]
for(i in 1: predition_t){
    logdata_prediction[t0+i] <- logdata_prediction[t0+i-12] + y_predict[t+2]
}

data_predict <- data.frame(exp(logdata_prediction))
t4 <- 1:length(data_predict[,1])
ggplot(data_predict,aes(x=t4, y=data_predict[,1])) + geom_line(color= ifelse(t4<=359, '#9999FF','red'), size=1.5) + ggtitle(" Predicted Time Series Plot") + 
    xlab("T (in Months)") + ylab("Industrial Output Levels") + 
    theme(plot.title =element_text(family="Times", face="bold.italic", size=18) ) + theme(axis.title = element_text(color = "blue"))
