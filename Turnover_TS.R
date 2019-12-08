df <- read.csv("C:\\Users\\DanBi Han\\Documents\\Documents\\Data Science\\WUDAC\\Penn_HR_Data\\PennHR\\Ts_Penn_HR.csv")
other.data <- read.csv("C:\\Users\\DanBi Han\\Documents\\Documents\\UPenn\\STAT 435\\ImportsfromJapan.txt")

df <- df[-c(1)]   # drop Pandas index 


# Prepping the data
Time <- seq(1: length(df$Turnover))
month <- rep(seq(1:12), 12)[7:138]
fMonth <- as.factor(month)
Year <- c(rep(2008, 6), rep(2009, 12), rep(2010, 12), rep(2011, 12),
          rep(2012, 12), rep(2013, 12), rep(2014, 12), rep(2015, 12),
          rep(2016, 12), rep(2017, 12), rep(2018, 12), rep(2019, 6))

log.data <- log(df$Turnover)


# Export calendar variables from Prof. Shaman's data.
c220 <- other.data$c220[7:138]
s220 <- other.data$s220[7:138]
c348 <- other.data$c348[7:138]
s348 <- other.data$s348[7:138]
c432 <- other.data$c432[7:138]
s432 <- other.data$s432[7:138]


# Attach the new variables to the df.
df <- data.frame(df, log.data, Time, Year, month, fMonth, 
                 c220, s220, c348, s348, c432, s432)


### ----------- Plot ts -------------
plot(ts(df$Turnover, start = c(2008, 7), end = c(2019, 6), freq = 12),
     xlab = 'Year',
     ylab = 'Rate (%)',
     main = 'Turnover vs. Time')


plot(ts(df$log.data, start = c(2008, 7), end = c(2019, 6), freq = 12),
     xlab = 'Year',
     ylab = 'Rate (%)', ylim = c(-5, 1),
     main = 'Logged Turnover vs. Time')

   # Log transformation takes care of volatility.
   # Notice: No trend.


# ------------- Observations ---------------
# Due to recession, 2008 data is unusable.
# Start arond July, 2009? Jan, 2011?
# (Recession hit rock bottom in June, 2009; took a while till recovery)


# ----------- Fit seasonal ARIMA -----------
# Chop off Recession. Start at July, 2009.


# 1) Plot the spectrum to get a sense of the structures in the data.
spectrum(df$log.data[13:132], span = 6)
   # seems like no calendar effect...


# 2) Fit regression to eliminate calendar effect + outliers
### Let's first see if there is any calendar effect.
model.log <- lm(log.data ~ fMonth + c348 + s348,
             data = df[13:132, ])

summary(model.log)

### Partial F-test to see if 0.348 pair should be kept
model.without <- lm(log.data ~ fMonth,
                    data = df[13:132, ])

anova(model.log, model.without)
   # Yes, keep the pair.


### Are there any outliers?
norm.plot <- qqnorm(resid(model.log))
qqline(resid(model.log))
identify(norm.plot)   # obs114

obs114 <- c(rep(0, length(df$log.data)))
obs114[114] <- 1
obs114 <- obs114[13:132]

model1 <- lm(log.data ~ fMonth + c348 + s348 + obs114,
             data = df[13:132, ])

summary(model1)   # no significant outliers!


### Regression to remove calendar effect
model.reg <- lm(log.data ~ c348 + s348, data = df[13:132, ])


# 3) Seasonal differencing to kill static seasonal
residuals.reg <- resid(model.reg)

residuals.diff <- diff(residuals.reg, 12)


# 4) Inspect ACF and PACF to determine ARIMA orders.
acf(residuals.diff, 40)   # lag 12 => dynamic seasonal remaining
pacf(residuals.diff, 40)   # lags 12, 24 => yep, dynamic seasonal


# 5) Fit seasonal ARIMA to regression residuals
model.ar <- arima(residuals.reg, order = c(0,0,0),
                  seasonal = list(order = c(2,1,0), period = 12))


### Which model is better? 
# Comparing the metrics below, 
# ARIMA(0,0,0)(2,1,0)12 did better than ARIMA(0,0,0)(1,1,0)12.

library('lmtest')
coeftest(model.ar)

model.ar   # aic = -9.46


# 6) Evaluate the ARIMA model.
### Remove the first 12 residuals since seasonally differenced.
sel <- 1:12
residuals.ar <- resid(model.ar)[-sel]


# Residual analysis
### Plot residuals
plot(ts(residuals.ar, start = c(2009, 7), freq = 12),
     xlab = 'Year',
     ylab = 'Residuals',
     main = 'Residuals vs. Year')


### ACF and PACF
acf(residuals.ar, 40)
pacf(residuals.ar, 40)


### Spectrum
spectrum(residuals.ar, span = 12)


### Bartlett's Kolmogorov-Smirnov Test
library('hwwntest')
bartlettB.test(ts(residuals.ar))   # p-value = 0.9341


### Clear reduction to white noise!


# 7) Let's get static and dynamic seasonal indices.

### Static seasonal indices using regression.
coefficients <- coef(model.log)

intercept <- coefficients[1]
month.dummies <- coefficients[2:12] + intercept

seas.indices <- c(intercept, month.dummies)
seasonal <- exp(seas.indices - mean(seas.indices))
print(seasonal)

plot(ts(seasonal),
     ylab = "Seasonal Indices",
     xlab = "Months",
     main = 'Static Seasonal Indices from Regression')


### ----- Using ARIMA -------
arimapred <- residuals.reg[-sel] - residuals.ar # already removed first 12 residuals
arimapred2 <- resid(lm(arimapred ~ poly(df$Time[1:length(arimapred)],4)))

monthly.means <- tapply(arimapred2, df$month[13:132][-sel], mean)
static.seas <- exp(monthly.means - mean(monthly.means))

print(static.seas)

# Plot the static estimates from reg vs. ARIMA to compare.
plot(ts(static.seas),
     ylab = "Seasonal Indices",
     xlab = "Months",
     main = 'Static Seasonal Indices',
     lty=1, lwd=2, col="red")

lines(ts(seasonal),
      lty=2, lwd=2, col="blue")

legend("topleft", legend = c("ARIMA", "Regression"),
       col = c("red", "blue"), lty = 1, cex = 0.8)
###-----------------


### Dynamic seasonal indices -- Throw out.
arimapred2.ts <- ts(arimapred2, start = c(2009, 7), freq = 12)

boxplot(arimapred2 ~ cycle(arimapred2.ts))


# Calculating dynamic indices
y <- arimapred2
seas.m <- matrix(rep(0,108), ncol = 9)
j <- -11

for (i in 1:9){
   j <- j + 12
   j1 <- j + 11
   
   seas.m[,i] <- exp(y[j:j1] - mean(y[j:j1]))
}

year <- seq(2011, 2019)
seas.m2 <- matrix(rep(static.seas, 9), ncol = 9)

name <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July',
          'Aug', 'Sept', 'Oct', 'Nov', 'Dec')

par(mfrow = c(1,1))

for (i in 1:12){   # change this range to get other months 
   plot(year, seas.m[i,], xlab = 'Year', ylab = 'Indices',
        main = name[i], type = 'l', lwd = 2, col = 'red',
        ylim = c(0.30, 2.3))
   
   lines(year, seas.m2[i,], lty = 1, lwd = 2, col = 'blue')
}


# 8) ARIMA forecasts
### Create new df, containing new data.
Time <- (1:12)
month <- (1:12)
fMonth <- as.factor(month)
Year <- c(rep(2019, 6), rep(2020, 6))


# Export calendar variables from Prof. Shaman's data.
c220 <- other.data$c220[139:150]
s220 <- other.data$s220[139:150]
c348 <- other.data$c348[139:150]
s348 <- other.data$s348[139:150]
c432 <- other.data$c432[139:150]
s432 <- other.data$s432[139:150]


# Attach the new variables to the df.
df.new <- data.frame(Time, Year, month, fMonth, 
                     c220, s220, c348, s348, c432, s432)

# Calculate forecasts for 1 year ahead.
f.reg <- predict(model.reg, newdata = df.new)
f.arima <- predict(model.ar, n.ahead = 12)

final.forecast <- exp(f.reg + f.arima$pred)

# Tabulate and plot the forecasts
cbind(c('July','Aug', 'Sept', 'Oct', 'Nov', 'Dec',
              'Jan', 'Feb', 'Mar', 'Apr', 'May', 'June'),
            final.forecast)

plot(ts(final.forecast, start = c(2019, 7), freq = 12),
     xlab = 'Time',
     ylab = 'Turnover %',
     main = 'Turnover Forecasts for FY 2020')
