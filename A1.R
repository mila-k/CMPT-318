# Assignment 1

library(ggplot2)  
library('psych') 
library(corrplot)
 
DataDf <- read.table("Group_Assignment_1_Dataset.txt", header = TRUE, sep = ",")
DataDf$Date <- as.POSIXlt(DataDf$Date, format = "%d/%m/%Y")

# Extract week 2
week2 <- with(DataDf, DataDf[Date >= as.POSIXlt("8/1/2007", format = "%d/%m/%Y") & Date <= as.POSIXlt("14/1/2007", format = "%d/%m/%Y"), ])
week2$TimeModified <- as.POSIXlt(week2$Time, format = "%H:%M:%S")

# Modify Time
week2 <- subset(week2, select = -c(Time))
Time <- subset(week2$TimeModified, week2$Date >= as.POSIXlt("8/1/2007",format = "%d/%m/%Y"))
week2 <- cbind(week2, Time)


## 1

getmode <- function(v) {  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
} # https://www.tutorialspoint.com/r/r_mean_median_mode.htm

wday <- c(as.POSIXlt(week2$Date, format = "%d/%m/%Y")$wday)
week2 <- cbind(week2, wday)
weekdaysDayHoursDf <- with(week2, week2[(wday != 0 & wday != 6) & TimeModified >= as.POSIXlt("06:00:00", format = "%H:%M:%S") & TimeModified < as.POSIXlt("18:00:00", format = "%H:%M:%S"), ])
weekdaysNightHoursDf <- with(week2, week2[(wday != 0 & wday != 6) & (TimeModified < as.POSIXlt("06:00:00", format = "%H:%M:%S") | TimeModified >= as.POSIXlt("18:00:00", format = "%H:%M:%S")), ])
weekendDayHoursDf <- with(week2, week2[(wday == 0 | wday == 6) & TimeModified >= as.POSIXlt("06:00:00", format = "%H:%M:%S") & TimeModified < as.POSIXlt("18:00:00", format = "%H:%M:%S"), ])
weekendNightHoursDf <- with(week2, week2[(wday == 0 | wday == 6) & (TimeModified < as.POSIXlt("06:00:00", format = "%H:%M:%S") | TimeModified >= as.POSIXlt("18:00:00", format = "%H:%M:%S")), ])

# Global_active_power (A)
globalActivePowerArithMean <- mean(week2$Global_active_power)
globalActivePowerGeoMean <- geometric.mean(week2$Global_active_power)
globalActivePowerMedian <- median(week2$Global_active_power, na.rm = TRUE)
globalActivePowerMode <- getmode(week2$Global_active_power)
globalActivePowerStdDev <- sd(week2$Global_active_power)

globalActivePowerWeekdayDayMax <- max(weekdaysDayHoursDf$Global_active_power)
globalActivePowerWeekdayDayMin <- min(weekdaysDayHoursDf$Global_active_power)
globalActivePowerWeekdayNightMax <- max(weekdaysNightHoursDf$Global_active_power)
globalActivePowerWeekdayNightMin <- min(weekdaysNightHoursDf$Global_active_power)

globalActivePowerWeekendDayMax <- max(weekendDayHoursDf$Global_active_power)
globalActivePowerWeekendDayMin <- min(weekendDayHoursDf$Global_active_power)
globalActivePowerWeekendNightMax <- max(weekendNightHoursDf$Global_active_power)
globalActivePowerWeekendNightMin <- min(weekendNightHoursDf$Global_active_power)

# Global_reactive_power (B)
globalReactivePowerArithMean <- mean(week2$Global_reactive_power, na.rm = TRUE)
globalReactivePowerGeoMean <- geometric.mean(week2$Global_reactive_power)
globalReactivePowerMedian <- median(week2$Global_reactive_power, na.rm = TRUE)
globalReactivePowerMode <- getmode(week2$Global_reactive_power)
globalReactivePowerStdDev <- sd(week2$Global_reactive_power, na.rm = TRUE)

globalReactivePowerWeekdayDayMax <- max(weekdaysDayHoursDf$Global_reactive_power)
globalReactivePowerWeekdayDayMin <- min(weekdaysDayHoursDf$Global_reactive_power)
globalReactivePowerWeekdayNightMax <- max(weekdaysNightHoursDf$Global_reactive_power)
globalReactivePowerWeekdayNightMin <- min(weekdaysNightHoursDf$Global_reactive_power)

globalReactivePowerWeekendDayMax <- max(weekendDayHoursDf$Global_reactive_power)
globalReactivePowerWeekendDayMin <- min(weekendDayHoursDf$Global_reactive_power)
globalReactivePowerWeekendNightMax <- max(weekendNightHoursDf$Global_reactive_power, na.rm = TRUE)
globalReactivePowerWeekendNightMin <- min(weekendNightHoursDf$Global_reactive_power, na.rm = TRUE)

# Voltage (C)
VoltageArithMean <- mean(week2$Voltage, na.rm = TRUE)
VoltageGeoMean <- geometric.mean(week2$Voltage, na.rm = TRUE)
VoltageMedian <- median(week2$Voltage, na.rm = TRUE)
VoltageMode <- getmode(week2$Voltage)
VoltageStdDev <- sd(week2$Voltage, na.rm = TRUE)




## 2

# Correlation coefficients
AB <- cor(week2$Global_active_power, week2$Global_reactive_power, method = "pearson", use = "complete.obs")
AC <- cor(week2$Global_active_power, week2$Voltage, method = "pearson", use = "complete.obs")
AD <- cor(week2$Global_active_power, week2$Global_intensity, method = "pearson", use = "complete.obs")
AE <- cor(week2$Global_active_power, week2$Sub_metering_1, method = "pearson", use = "complete.obs")
AF <- cor(week2$Global_active_power, week2$Sub_metering_2, method = "pearson", use = "complete.obs")
AG <- cor(week2$Global_active_power, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

BC <- cor(week2$Global_reactive_power, week2$Voltage, method = "pearson", use = "complete.obs")
BD <- cor(week2$Global_reactive_power, week2$Global_intensity, method = "pearson", use = "complete.obs")
BE <- cor(week2$Global_reactive_power, week2$Sub_metering_1, method = "pearson", use = "complete.obs")
BF <- cor(week2$Global_reactive_power, week2$Sub_metering_2, method = "pearson", use = "complete.obs")
BG <- cor(week2$Global_reactive_power, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

CD <- cor(week2$Voltage, week2$Global_intensity, method = "pearson", use = "complete.obs")
CE <- cor(week2$Voltage, week2$Sub_metering_1, method = "pearson", use = "complete.obs")
CF <- cor(week2$Voltage, week2$Sub_metering_2, method = "pearson", use = "complete.obs")
CG <- cor(week2$Voltage, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

DE <- cor(week2$Global_intensity, week2$Sub_metering_1, method = "pearson", use = "complete.obs")
DF <- cor(week2$Global_intensity, week2$Sub_metering_2, method = "pearson", use = "complete.obs")
DG <- cor(week2$Global_intensity, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

EF <- cor(week2$Sub_metering_1, week2$Sub_metering_2, method = "pearson", use = "complete.obs")
EG <- cor(week2$Sub_metering_1, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

FG <- cor(week2$Sub_metering_2, week2$Sub_metering_3, method = "pearson", use = "complete.obs")

# correlation matrix
keeps <- c("Global_active_power","Global_reactive_power", "Voltage", "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
variablesDf = week2[keeps]
corrCoeffs <- cor(variablesDf, method = "pearson", use = "complete.obs")

matrix1 <- corrplot(corrCoeffs, method="color")
matrix2 <- corrplot(corrCoeffs, method="number")





## 3

dayHoursAvgDf <- with(week2, week2[(wday == 1) & TimeModified >= as.POSIXlt("06:00:00", format = "%H:%M:%S") & TimeModified < as.POSIXlt("18:00:00", format = "%H:%M:%S"), ])
dayHoursAvgDf <- subset(dayHoursAvgDf, select=c("TimeModified", "Time"))
nightHoursAvgDf <- with(week2, week2[(wday == 0) & (TimeModified < as.POSIXlt("06:00:00", format = "%H:%M:%S") | TimeModified >= as.POSIXlt("18:00:00", format = "%H:%M:%S")), ])
nightHoursAvgDf <- subset(nightHoursAvgDf, select=c("TimeModified", "Time"))

weekdayDayGlobalIntensity <- vector(length=nrow(dayHoursAvgDf))
weekdayNightGlobalIntensity <- vector(length=nrow(nightHoursAvgDf))
weekendDayGlobalIntensity <- vector(length=nrow(dayHoursAvgDf))
weekendNightGlobalIntensity <- vector(length=nrow(nightHoursAvgDf))


# average for weekday day hours
for(i in 1:nrow(dayHoursAvgDf)) 
{
 timeDf <- with(weekdaysDayHoursDf, weekdaysDayHoursDf[TimeModified == dayHoursAvgDf[i,1], ])
 average <- mean(timeDf$Global_intensity, na.rm = TRUE)
 weekdayDayGlobalIntensity[i] <- average
}
dayHoursAvgDf$weekday_average <- weekdayDayGlobalIntensity

# average for weekday night hours
for(i in 1:nrow(nightHoursAvgDf)) 
{
  timeDf <- with(weekdaysNightHoursDf, weekdaysNightHoursDf[TimeModified == nightHoursAvgDf[i,1], ])
  average <- mean(timeDf$Global_intensity, na.rm = TRUE)
  weekdayNightGlobalIntensity[i] <- average
}
nightHoursAvgDf$weekday_average <- weekdayNightGlobalIntensity

# average for weekend day hours
for(i in 1:nrow(dayHoursAvgDf)) 
{
  timeDf <- with(weekendDayHoursDf, weekendDayHoursDf[TimeModified == dayHoursAvgDf[i,1], ])
  average <- mean(timeDf$Global_intensity, na.rm = TRUE)
  weekendDayGlobalIntensity[i] <- average
}
dayHoursAvgDf$weekend_average <- weekendDayGlobalIntensity

# average for weekend night hours
for(i in 1:nrow(nightHoursAvgDf)) 
{
  timeDf <- with(weekendNightHoursDf, weekendNightHoursDf[TimeModified == nightHoursAvgDf[i,1], ])
  average <- mean(timeDf$Global_intensity, na.rm = TRUE)
  weekendNightGlobalIntensity[i] <- average
}
nightHoursAvgDf$weekend_average <- weekendNightGlobalIntensity

dayHoursAvgDf <- subset(dayHoursAvgDf, select = -c(TimeModified))
nightHoursAvgDf <- subset(nightHoursAvgDf, select = -c(TimeModified))


# Linear Regression

fit <- lm(formula=weekday_average ~ Time, data=dayHoursAvgDf)
weekday_predicted <- predict(fit, dayHoursAvgDf)
dayHoursAvgDf = cbind(dayHoursAvgDf, weekday_predicted)

fit <- lm(formula=weekend_average ~ Time, data=dayHoursAvgDf)
weekend_predicted <- predict(fit, dayHoursAvgDf)
dayHoursAvgDf = cbind(dayHoursAvgDf, weekend_predicted)

fit <- lm(formula=weekday_average ~ Time, data=nightHoursAvgDf)
weekday_predicted <- predict(fit, nightHoursAvgDf)
nightHoursAvgDf = cbind(nightHoursAvgDf, weekday_predicted)

fit <- lm(formula=weekend_average ~ Time, data=nightHoursAvgDf)
weekend_predicted <- predict(fit, nightHoursAvgDf)
nightHoursAvgDf = cbind(nightHoursAvgDf, weekend_predicted)


# Plot Linear Regression

ggplot(dayHoursAvgDf, aes(x=Time)) + 
  geom_line(aes(y = weekday_predicted), color = "darkred") + 
  geom_point(aes(y = weekday_average), color="darkred") +
  geom_line(aes(y = weekend_predicted), color="steelblue", linetype="twodash") +
  geom_point(aes(y = weekend_average), color="steelblue")

