# Clear my environment
rm(list=ls())

# set working directory
setwd("C:/Users/saurabh/Desktop/Time Series Analysis/Assignment_3")

# upload the following libraries 
library(TSA)
library(fUnitRoots)
library(lmtest)
library(FitAR)
library(forecast)

# source Utility functions
source('TSHandy.r')
source('residual.analysis.R')
source('sort.score.R')

# reading in the raw data
sunshine = read.csv("bom_timeseries.csv",header=TRUE)
sunshine
  
#Converting data to time series object 
# parse different date formats to one format
#sunshine$Date <-parse_date_time(sunshine$Date, orders = c('dmy', 'ymd'))

# convert the date column to datetime object
#sunshine$Date <- as.POSIXct(sunshine$Date, format="%y/%m/%d")

sunshine.ts <- ts(as.vector(sunshine$Sunshine_hrs),  frequency=7)
class(sunshine.ts)

##time plot
par(mfrow=c(1,1))
plot(sunshine.ts,ylab='Sunshine hrs',main = "Time series plot for daily sunshine hrs")

# summary statistics
summary(sunshine.ts)
# No negative values 
# No missing values. 
# sample size is good

# correlation
Originaldata = sunshine.ts              
Firstlag     = zlag(sunshine.ts)           
index        = 2:length(Firstlag)                   
cor(Originaldata[index],Firstlag[index]) 
# Positive correlation. 
# 0.35

# Time series Plot
par(mfrow=c(1,1))
plot(sunshine.ts, xlab="Year", ylab = "Daily Sunshine hrs",main ="Time series plot of Sunshine (hrs)") 
# There is seasonality and hence it is hard to predict the behaviour of the series
# Trend                                   - yes
# Changing Variance                       - yes
# Auto regressive behaviour - yes and MA  - no- MA yes
# Intervention Point                      - no
# Seasonality                             - yes

## scatter plot
par(mfrow=c(1,1))
plot(y=sunshine.ts,x=Firstlag,ylab='Solar Exposure', xlab='Previous year Sunshine',main = "Scatter plot of Daily Sunshine hrs")
# Correlation with previous years

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Fit deterministic Model
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##linear trend model
model1.linear = lm(sunshine.ts~time(sunshine.ts))
summary(model1.linear)

# Add the fitted least squares line from the linear model
par(mfrow=c(1,1))
plot(sunshine.ts,ylab='y', main = "Fitted linear model to Sunshine series.")
abline(model1.linear) 
##parameters are insignificant

##Quadratic trend model
t = time(sunshine.ts)
t2 = t^2
model2.quadratic = lm(sunshine.ts~t+t2)  
summary(model2.quadratic)
##parameters are insignificant

##plot
plot(ts(fitted(model2.quadratic)), ylim = c(min(c(fitted(model2.quadratic),as.vector(sunshine.ts))), max(c(fitted(model2.quadratic),as.vector(sunshine.ts)))),ylab='y' ,
     main = "Fitted quadratic curve to the the Sunshine series.")
lines(as.vector(sunshine.ts))

## lets fit a harmonic model
## sine/cosine one 
har.=harmonic(sunshine.ts,1)
model3.harmonic=lm(sunshine.ts~har.)
summary(model3.harmonic) 

plot(ts(fitted(model3.harmonic)), ylim = c(min(c(fitted(model3.harmonic),
                                                         as.vector(sunshine.ts))), max(c(fitted(model3.harmonic),as.vector(sunshine.ts)))),
     ylab='y' , main = "Fitted quadratic curve to random walk data", type="l",lty=2,col="red")
lines(as.vector(sunshine.ts))
## it does not captures the period behaviour 
 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Fit Stochastic Trend Models
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##approach
## First eliminating the effect of any kind of trend.
## And then proceeding with the usual analysis (modeling theautocorrelation structure in the series)

# ACF/PACF Plot
par(mfrow=c(2,1))
a=acf(sunshine.ts)
a$acf
pacf(sunshine.ts)
# ACF slowly decaying pattern and repeating patterns. Seasonality
# Non-stationary series.
 
# Unit root test
ar(diff(sunshine.ts))   
adfTest(sunshine.ts, lags = 22,  title = NULL,description = NULL)
# series nonstationary

# check Normality
par(mfrow=c(1,1))
qqnorm(sunshine.ts,main="QQ plot of solar exposure")
qqline(sunshine.ts, col = 2)
shapiro.test(sunshine.ts)
# data sample is not normal

## lets stabilise the variance and make the sample closer to normal distribution
options(warn=-1)
BC.sunshine=BoxCox.ar(sunshine.ts+0.01)
 
title(main = " Log-likelihood versus the values of lambda for Egg deposition series.")
# not a smooth curve for box cox
# Let's use MoM methof 

##BC.sunshine = BoxCox.ar(sunshine.ts+0.01, method = "yule-walker")
##title(main = " Log-likelihood versus the values of lambda for Solar exposure series.")
BC.sunshine$ci   ## confidence intervals close to 1
lambda = (0.7+0.8)/2
lambda
BC.sunshine = (sunshine.ts^lambda-1)/lambda
 

# Time series plot
par(mfrow=c(2,1))
plot(sunshine.ts,ylab='solar exposure series',main='solar exposure series')
plot(BC.sunshine,ylab='Box cox transformed solar exposure series',main='Box cox transformed solar exposure series')
# variance much better in boxcox transformed series
summary(BC.sunshine)




# QQ plot
par(mfrow=c(2,1))
qqnorm(sunshine.ts,main="normality of solar exposure series")
qqline(sunshine.ts, col = 2)

qqnorm(BC.sunshine,main="normality of BoxCox Transformed solar exposure series")
qqline(BC.sunshine, col = 2)

##histogram
par(mfrow=c(2,1))
hist(sunshine.ts,main="Histogram of sunshine series")
hist(BC.sunshine,main="Histogram of BoxCox transformed sunshine series")
##sample large  >30
## we keep the box cox transformed series
  
# Let's detrend using the boxcox transformed series
diff.BC.sunshine = diff(BC.sunshine,differences = 1)
par(mfrow=c(1,1))
plot(diff.BC.sunshine,ylab='First difference',main="Time series plot of first differenced series")

# ##trend exists
# diff.BC.sunshine = diff(BC.sunshine,differences = 2)
# par(mfrow=c(1,1))
# plot(diff.BC.sunshine,type='o',ylab='Second difference',main="Time series plot of Second differenced series")
# # trend removed

# unit root test
ar(diff(diff.BC.sunshine))
adfTest(diff.BC.sunshine, lags = 18,  title = NULL,description = NULL)
## suggests stationary

# Let's check the Model
# check ACF/PACF
par(mfrow=c(2,1))
a=acf(diff.BC.sunshine)
a$acf
pacf(diff.BC.sunshine)

## alternative bounds
par(mfrow=c(1,1))
acf(diff.BC.sunshine,ci.type='ma',xaxp=c(0,10,10),main="Figure 2. Alternative bounds for the ACF.")
# ACF 2 significant lags
# Pattern in PACF hence p=0
# ARIMA(0,2,2)

# EACF
eacf(diff.BC.sunshine, ar.max = 10, ma.max =10)
#EACF suggests ARIMA (0,2,2),ARIMA(0,2,3),ARIMA (1,2,1),ARIMA(1,2,2),ARIMA(1,2,2),ARIMA(1,2,3)

#BIC Table
par(mfrow=c(1,1))
res = armasubsets(y=diff.BC.sunshine,nar=10,nma=10,y.name='test',ar.method='ols')
plot(res)
# Shaded columns are AR(1),AR(2),AR(3),MA(1),MA(2),MA(4)
# BIC Suggests ARIMA (1,2,1),ARIMA(1,2,2),ARIMA(1,2,4),ARIMA(2,2,1),ARIMA(2,2,2),ARIMA(2,2,4),ARIMA(3,2,1),ARIMA(3,2,2),ARIMA(3,2,4)

##CAndidate Models
# ARIMA(0,2,2) , ARIMA(0,2,3),ARIMA (1,2,1) , ARIMA(1,2,2), ARIMA(1,2,3),ARIMA(1,2,4),ARIMA(2,2,1),ARIMA(2,2,2),ARIMA(2,2,4),ARIMA(3,2,1),ARIMA(3,2,2),ARIMA(3,2,4)

# Model Fitting - MOM
modelList <- list(c(0,2,2),c(0,2,3),c(1,2,1),c(1,2,2),c(1,2,3),c(1,2,4),c(2,2,1),c(2,2,2),c(2,2,4),c(3,2,1),c(3,2,2),c(3,2,4))
modelEstimation <- myCandidate(BC.sunshine, orderList = modelList, methodType = "CSS")
modelEstimation$significanceTest
 
# Model Fitting - MLE
modelList <- list(c(0,2,2),c(0,2,3),c(1,2,1),c(1,2,2),c(1,2,3),c(1,2,4),c(2,2,1),c(2,2,2),c(2,2,4),c(3,2,1),c(3,2,2),c(3,2,4))
modelEstimation <- myCandidate(BC.sunshine, orderList = modelList, methodType = "ML")
modelEstimation$significanceTest

#AIC/BIC
modelEstimation$IC
# AIC and BIC suggests ARIMA(0,2,3)


#OVERFITTING
# Model Fitting - MOM
modelList <- list(c(0,2,3),c(0,2,4),c(1,2,3),c(1,2,4))
modelEstimation <- myCandidate(BC.sunshine, orderList = modelList, methodType = "CSS")
modelEstimation$significanceTest

# Model Fitting - MLE
modelList <- list(c(0,2,3),c(0,2,4),c(1,2,3),c(1,2,4))
modelEstimation <- myCandidate(BC.sunshine, orderList = modelList, methodType = "ML")
modelEstimation$significanceTest

#AIC/BIC
modelEstimation$IC
#ARIMA(0,2,3) ## winner


# residual analysis
residual.analysis(modelEstimation$model[[1]], std = FALSE,start = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Fit Stochastic  seasonal  trends
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##time plot
par(mfrow=c(1,1))
plot(sunshine.ts,ylab='Sunshine hrs',main = "Time series plot for daily sunshine hrs")

par(mfrow=c(1,2))
acf( sunshine.ts,  lag.max = 36,main="The sample  ACF of sunshine series")
pacf(sunshine.ts,  lag.max = 36,main="The sample PACF of sunshine series")
## wave pattern seasonality. seasonal periods there is decreasing pattern
## We identify seasonal trend and hence perform D=1
 
## we will follow the residuals approach
## we will go step by step
## At every step we will fit one model
## D=1 and fit and see how it goes
# Seasonality and existence of trend are obvious from the ACF and PACF plots

# First fit a plain model with only the first seasonal difference with order D = 1 
# and see if we can get rid of the seasonal trend effect
# by inspecting the autocorrelation structure of the residuals.

#                                       p,d,q                        P,D,Q
m1.sunshine = arima(sunshine.ts,order=c(0,0,0),seasonal=list(order=c(0,1,0), period=7))
res.m1 = residuals(m1.sunshine);  
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
# TS plot: we can conclude that we got rid of the trend but there is variance

## important in ACF nd PACF we just check the significant lags on the Period 1 that is 12 months
acf(res.m1, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m1, lag.max = 36, main = "The sample PACF of the residuals")
# Look at ACF and PACF plots over the lags corresponding to the periods namely 7, 14, 18, etc.
# (Intermediate lags are used to determine p and q of the ordinary trend)
# We have one significant correlation at the first seasonal lag in  ACF and PACF there is a pattern. 
# This indicates P=0, Q=1.

## based on the above ACF pacf include P and Q in the model


#So, we will add the SARMA(0,1) component and see if we get rid of seasonal component.
m2.sunshine = arima(sunshine.ts,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=7))
res.m2 = residuals(m2.sunshine);  


par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")

par(mfrow=c(1,2))
acf(res.m2, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m2, lag.max = 36, main = "The sample PACF of the residuals")
# There are no significant correlation on the periods. we can conclude that the seasonality is filtered out.
# Even if there is one then If i handle it suitable in a good model then in the end it will be white noise residual

# Now, we will specify the orders of ARIMA component. 
# First of all: There was change in variation in the TS plot of the original series. We need to deal with this.

## In ACF/PACF lets look at the ordinary part before the seasonal lag.
## one significant lag and then after some period again many lags then it is sign of changing variance

## lets stabilise the variance and make the sample closer to normal distribution
options(warn=-1)
BC.sunshine=BoxCox.ar(sunshine.ts+0.01)

title(main = " Log-likelihood versus the values of lambda for Sunshine series.")
# not a smooth curve for box cox
# Let's use MoM methof 
##BC.sunshine = BoxCox.ar(sunshine.ts+0.01, method = "yule-walker")
##title(main = " Log-likelihood versus the values of lambda for Solar exposure series.")
BC.sunshine$ci   ## confidence intervals close to 1
lambda = (0.7+0.8)/2
lambda
BC.sunshine = (sunshine.ts^lambda-1)/lambda

##time plot
par(mfrow=c(1,1))
plot(BC.sunshine,ylab='Box Cox Transformed  Series',xlab='Days', main = "TS plot of Box Cox Transformed sunshine series.")
# # variance stabilized

## unit root test

# and now see if we can see the ordinary trend more clearly.
m3.sunshine = arima(BC.sunshine,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=7))
res.m3 = residuals(m3.sunshine);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m3, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m3, lag.max = 36, main = "The sample PACF of the residuals")
# We have a very high correlation at the first lag of PACF  and ACF. 



# There is an ordinary trend. Lets detrend it.hence d=0
#                                       p d q
m4.sunshine = arima(BC.sunshine,order=c(0,1,0),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m4.sunshine)
res.m4 = residuals(m4.sunshine);  
par(mfrow=c(1,1))
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m4, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m4, lag.max = 36, main = "The sample PACF of the residuals")
# From ACF: q=2 or q=5 
# From PACF p=0 or p=9
# d=1

eacf(res.m4)
# SARIMA(0,1,2)x(0,1,1)_7 
# SARIMA(0,1,3)x(0,1,1)_7
# SARIMA(1,1,2)x(0,1,1)_7
# SARIMA(1,1,1)x(0,1,1)_7
# SARIMA(2,1,0)x(0,1,1)_7


# The tentative models are specified as 
# SARIMA(0,1,2)x(0,1,1)_7 by ACF and PACF
# SARIMA(0,1,5)x(0,1,1)_7 by ACF and PACF  
# SARIMA(9,1,2)x(0,1,1)_7 by ACF and PACF
# SARIMA(9,1,5)x(0,1,1)_7 by ACF and PACF
# SARIMA(0,1,2)x(0,1,1)_7 by EACF
# SARIMA(0,1,3)x(0,1,1)_7 by EACF
# SARIMA(1,1,2)x(0,1,1)_7 by EACF
# SARIMA(1,1,1)x(0,1,1)_7 by EACF
# SARIMA(2,1,0)x(0,1,1)_7 by EACF


## lets fit these models and see how it goes in the ordinary bit without the s...look for significance
# check if residuals are white noise (i.e. uncorrelated) in the orddinary parameter not seasonal parameter
##ljung box test will actually tell if the significance is really true or not..it will catch it properly.

# SARIMA(0,1,2)x(0,1,1)_7 by ACF and PACF
m4_012.sunshine = arima(BC.sunshine,order=c(0,1,2),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m4_012.sunshine)
res.m4 = residuals(m4_012.sunshine);  
par(mfrow=c(1,1))
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m4, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m4, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(0,1,5)x(0,1,1)_7 by ACF and PACF 
m5_015.sunshine = arima(BC.sunshine,order=c(0,1,5),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m5_015.sunshine)
res.m5 = residuals(m5_015.sunshine);  
par(mfrow=c(1,1))
plot(res.m5,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m5, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m5, lag.max = 36, main = "The sample PACF of the residuals")
 

# SARIMA(9,1,2)x(0,1,1)_7 by ACF and PACF 
m6_912.sunshine = arima(BC.sunshine,order=c(9,1,2),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m6_912.sunshine)
res.m6 = residuals(m6_912.sunshine);  
par(mfrow=c(1,1))
plot(res.m6,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m6, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m6, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(9,1,5)x(0,1,1)_7 by ACF and PACF
m7_915.sunshine = arima(BC.sunshine,order=c(9,1,5),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m7_915.sunshine)
res.m7 = residuals(m7_915.sunshine);  
par(mfrow=c(1,1))
plot(res.m7,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m7, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m7, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(0,1,2)x(0,1,1)_7 by EACF
m8_012.sunshine = arima(BC.sunshine,order=c(0,1,2),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m8_012.sunshine)
res.m8 = residuals(m8_012.sunshine);  
par(mfrow=c(1,1))
plot(res.m8,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m8, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m8, lag.max = 36, main = "The sample PACF of the residuals")


# SARIMA(0,1,3)x(0,1,1)_7 by EACF
m9_013.sunshine = arima(BC.sunshine,order=c(0,1,3),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m9_013.sunshine)
res.m9 = residuals(m9_013.sunshine);  
par(mfrow=c(1,1))
plot(res.m9,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m9, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m9, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(1,1,2)x(0,1,1)_7 by EACF
m10_112.sunshine = arima(BC.sunshine,order=c(1,1,2),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m10_112.sunshine)
res.m10 = residuals(m10_112.sunshine);  
par(mfrow=c(1,1))
plot(res.m10,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m10, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m10, lag.max = 36, main = "The sample PACF of the residuals")


# SARIMA(1,1,1)x(0,1,1)_7 by EACF
m11_111.sunshine = arima(BC.sunshine,order=c(1,1,1),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m11_111.sunshine)
res.m11 = residuals(m11_111.sunshine);  
par(mfrow=c(1,1))
plot(res.m11,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m11, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m11, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(2,1,0)x(0,1,1)_7 by EACF
m12_210.sunshine = arima(BC.sunshine,order=c(2,1,0),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m12_210.sunshine)
res.m12 = residuals(m12_210.sunshine);  
par(mfrow=c(1,1))
plot(res.m12,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m12, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m12, lag.max = 36, main = "The sample PACF of the residuals")

##AIC
sc.AIC=AIC(m4_012.sunshine,m5_015.sunshine,m6_912.sunshine,m7_915.sunshine,m8_012.sunshine,m9_013.sunshine,m10_112.sunshine,m11_111.sunshine,m12_210.sunshine)
sort.score(sc.AIC, score = "aic")

## computes BIC
AIC(m4_012.sunshine, k = log(length(sunshine.ts)))
AIC(m5_015.sunshine, k = log(length(sunshine.ts)))
AIC(m6_912.sunshine, k = log(length(sunshine.ts)))
AIC(m7_915.sunshine, k = log(length(sunshine.ts)))
AIC(m8_012.sunshine, k = log(length(sunshine.ts)))
AIC(m9_013.sunshine, k = log(length(sunshine.ts)))
AIC(m10_112.sunshine, k = log(length(sunshine.ts)))
AIC(m11_111.sunshine, k = log(length(sunshine.ts)))
AIC(m12_210.sunshine, k = log(length(sunshine.ts)))
 

## AIC and BIC m4_012.sunshine
residual.analysis(model = m4_012.sunshine)


## OVERFIT
# SARIMA(0,1,3)x(0,1,1)_7 
m13_013.sunshine = arima(BC.sunshine,order=c(0,1,3),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m13_013.sunshine)
res.m13 = residuals(m13_013.sunshine);  
par(mfrow=c(1,1))
plot(res.m13,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m13, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m13, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(1,1,2)x(0,1,1)_7 
m14_112.sunshine = arima(BC.sunshine,order=c(1,1,2),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m14_112.sunshine)
res.m14 = residuals(m14_112.sunshine);  
par(mfrow=c(1,1))
plot(res.m14,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m14, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m14, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(1,1,3)x(0,1,1)_7 
m15_113.sunshine = arima(BC.sunshine,order=c(1,1,3),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m15_113.sunshine)
res.m15 = residuals(m15_113.sunshine);  
par(mfrow=c(1,1))
plot(res.m15,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m15, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m15, lag.max = 36, main = "The sample PACF of the residuals")

# SARIMA(0,1,1)x(0,1,1)_7 
m16_011.sunshine = arima(BC.sunshine,order=c(0,1,1),seasonal=list(order=c(0,1,1), period=7),method = "ML")
coeftest(m16_011.sunshine)
res.m16 = residuals(m16_011.sunshine);  
par(mfrow=c(1,1))
plot(res.m11,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m16, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m16, lag.max = 36, main = "The sample PACF of the residuals")

##AIC
sc.AIC=AIC(m4_012.sunshine, m13_013.sunshine, m14_112.sunshine , m15_113.sunshine,m16_011.sunshine)
sort.score(sc.AIC, score = "aic")

## computes BIC
AIC(m4_012.sunshine,  k = log(length(sunshine.ts)))
AIC(m13_013.sunshine,  k = log(length(sunshine.ts)))
AIC(m14_112.sunshine, k = log(length(sunshine.ts)))
AIC(m15_113.sunshine, k = log(length(sunshine.ts)))
AIC(m16_011.sunshine, k = log(length(sunshine.ts)))
##  m4_012.sunshine  winner

residual.analysis(model = m4_012.sunshine)

par(mfrow=c(1,1))
m4.sunshine = Arima((sunshine.ts),order=c(0,1,2),seasonal=list(order=c(0,1,1), period=7))
future = forecast(m4.sunshine, h = 60)
plot(future)
accuracy(m4.sunshine)

#-------------------------------------------------------------------

# Seasonal Decomposition
decompose(sunshine.ts)

plot(decompose(sunshine.ts))

# Using the stl method
#Seasonal Decomposition of Time Series by Loess
plot(stl(sunshine.ts, s.window = 7))

# stl forecasting
#Forecasting using stl objects
plot(stlf(sunshine.ts, method = "ets"))

# comparison with a standard ets forecast
plot(forecast(ets(sunshine.ts), h = 24))

# using autoplot
library(ggplot2)
autoplot(stlf(sunshine.ts, method = "ets"))

## Seasonal Arima (package forecast)
auto.arima(sunshine.ts, stepwise = T, 
           approximation = F, trace = T)

# Getting an object
sunshinearima = auto.arima(sunshine.ts, 
                             stepwise = T, 
                             approximation = F, 
                             trace = T)

# Forecast
forec = forecast(sunshinearima, h = 60)
plot(forec)

## Exponential Smoothing with ets
# Auto gemerated
ets(sunshine.ts)
# Forecast plot
sunshine.ts.ets = ets(sunshine.ts)

plot(forecast(sunshine.ts, h = 60))

# Comparison with seasonal Holt Winters model
hw <-hw(sunshine.ts, h = 60)
plot(hw(sunshine.ts, h = 60))


