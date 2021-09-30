
#######################################
#######################################
######### Time Series Project #########
########## Felipe Montenegro ##########
#######################################
#######################################



# 1. Importing Data and Packages
##################################################

# setting path and importing data
setwd("data")
getwd()
datafile <- "valeurs_mensuelles.csv"
data <- read.csv(datafile, sep=";")


# cheking data
head(data, n=10)


# importing time series packages
require(zoo)
require(tseries)
require(fUnitRoots)



# 2. Treating Data
##################################################

# changing dates type
typeof(data$Dates) # integer
dates_char <- as.character(data$Dates)
typeof(dates_char) # character


# checking size of sample
dates_char[1]; dates_char[length(dates_char)]


# changing date format
dates <- as.yearmon(seq(from=2000+7/12, to=2013+8/12, by=1/12))


# Creating time series Serie (double)
X_t <- zoo(data$Indice, order.by=dates)
typeof(X_t)



# 3. Time series, Log and Log Residues
##################################################

# time series plot
# try to find a tendence or saisonnality
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = "Indice brut de la production industrielle") 
acf(X_t, main = "")


# time series log plot
W_t = log(X_t)
plot(W_t, xlab = "Années", ylab = "Logarithme de l'indice") 
acf(W_t, main = "")


# tendence correction
lt <- lm(W_t ~ dates) # log regression
summary(lt)


# log regression residues
Z_t <- lt$residuals
acf(Z_t, main = ""); pacf(Z_t, main = "")



# 4. Unit root test on the logarithm of the series
##################################################


# Dickey-Fuller Test
adf <- adfTest(W_t, lag=0, type="ct")
adf


# autocorrelation test function
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients))


# Phillips-Perron Test
pp.test(W_t)



# 5. Comparison with the original Time series
##################################################

# plots before and after the transformation
plot(X_t, xlab = "Années", ylab = expression(X[t])) ; plot(Z_t, xlab = "Années", ylab = expression(Z[t]))



# 6. Modeling
##################################################

# Discover p* and q* values
acf(Z_t); pacf(Z_t)


# Validation Test
# Null Coefficients Test

signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}


# p = 0 and q = 1
arima001 <- arima(Z_t,c(0,0,1))
Qtests(arima001$residuals, 24, 1)
signif(arima001)

# p = 0 and q = 2
arima002 <- arima(Z_t,c(0,0,2))
Qtests(arima002$residuals, 24, 2)
signif(arima002)
ma2 <- arima002

# p = 1 and q = 1
arima101 <- arima(Z_t,c(1,0,1))
Qtests(arima101$residuals, 24, 2)
signif(arima101)

# p = 1 and q = 0
arima100 <- arima(Z_t,c(1,0,0))
Qtests(arima100$residuals, 24, 1)
signif(arima100)
ar1 <- arima100

# p = 1 and q = 2
arima102 <- arima(Z_t,c(1,0,2))
Qtests(arima102$residuals, 24, 3)
signif(arima102)                    


# AIC and BIC
models <- c("ma2", "ar1"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))



# 7. Prediction
##################################################

# prediction
prediction <- predict(ma2,2)
prediction


# sigma's noise
var_prev_1 <- ma2$sigma2
var_prev_2 <- ma2$sigma2*(1+(ma2$coef[1] + ma2$coef[2])**2)
ma2$sigma2


# theoretical formula
Bound_sup <- rbind(prediction$pred[1] +1.96*sqrt(var_prev_1),prediction$pred[2] +1.96*sqrt(var_prev_2))
Bound_inf <- rbind(prediction$pred[1] -1.96*sqrt(var_prev_1),prediction$pred[2] -1.96*sqrt(var_prev_2))


# sigma2 gives us the MLE of the innovations variance.
date_pred <- c(2013+9/12, 2013+10/12) 
coef_reg <- lt$coefficients


par(mfrow = c(1,1))
plot(NULL,NULL,xlim=c(2011+7/12, 2013+10/12),ylim=c(400,800),
     xlab="Années",ylab=expression(X[t]))
lines(X_t, type = "o", pch = 16)


polygon(x=c(date_pred, rev(date_pred)),
        y=c(c(exp(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1]), 
              exp(Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1])),
            rev(c(exp(Bound_sup[1] + coef_reg[2]*date_pred[1] + coef_reg[1]), 
                  exp(Bound_sup[2] + coef_reg[2]*date_pred[1] + coef_reg[1])))),col="grey",border=NA)
lines(c(X_t,exp(prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1])) )
points(date_pred, exp(prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1]), col="red", pch=16)



# 8. Report's Annexe
##################################################

par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = expression(X[t])); acf(X_t, main = "")
plot(W_t, xlab = "Années", ylab = expression(W[t])); acf(W_t, main = "")
plot(Z_t, xlab = "Années", ylab = expression(Z[t])); acf(Z_t, main = "")

