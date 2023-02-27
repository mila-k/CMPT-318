### Assignment 3

library(depmixS4)
library(zoo)
library(ggplot2)

DataDf <- read.table("Group_Assignment_3_Dataset.txt", header = TRUE, sep = ",")

# linear interpolation on NA values
DataDf$Global_active_power <- na.approx(DataDf$Global_active_power)
DataDf$Global_reactive_power <- na.approx(DataDf$Global_reactive_power)
DataDf$Voltage <- na.approx(DataDf$Voltage)
DataDf$Global_intensity <- na.approx(DataDf$Global_intensity)
DataDf$Sub_metering_1 <- na.approx(DataDf$Sub_metering_1)
DataDf$Sub_metering_2 <- na.approx(DataDf$Sub_metering_2)
DataDf$Sub_metering_3 <- na.approx(DataDf$Sub_metering_3)




### 1

# Scale the dataset
DataDf$Global_active_power <- scale(DataDf$Global_active_power, center = TRUE, scale = TRUE)
DataDf$Global_reactive_power <- scale(DataDf$Global_reactive_power, center = TRUE, scale = TRUE)
DataDf$Global_intensity <- scale(DataDf$Global_intensity, center = TRUE, scale = TRUE)
DataDf$Voltage <- scale(DataDf$Voltage, center = TRUE, scale = TRUE)
DataDf$Sub_metering_1 <- scale(DataDf$Sub_metering_1, center = TRUE, scale = TRUE)
DataDf$Sub_metering_2 <- scale(DataDf$Sub_metering_2, center = TRUE, scale = TRUE)
DataDf$Sub_metering_3 <- scale(DataDf$Sub_metering_3, center = TRUE, scale = TRUE)


# Determine a time window for a specific weekday 
DataDf$week_num <- strftime(as.Date(DataDf$Date, format = "%d/%m/%Y"), format = "%V")
DataDf <- subset(DataDf, DataDf$Date != "31/12/2007")
WeeksDf <- data.frame(matrix(nrow = 10080, ncol = 0))
for (i in 1:52) 
{
  week <- subset(DataDf, DataDf$week_num == i | DataDf$week_num == paste0("0", i))
  WeeksDf <- cbind(WeeksDf, week$Global_intensity)
  colnames(WeeksDf)[i] <- paste0("week_", i)
}
WeeksDf$average <- rowMeans(WeeksDf)
WeeksDf$Datetime <- with(DataDf, as.POSIXct(paste(Date, Time), format="%d/%m/%Y %H:%M:%S"))[1:10080]

ggplot(WeeksDf[1440:2880, ]) +
  geom_line(aes(x = Datetime, y= average, colour="Average"), linewidth=0.5) # Choose Tuesday 12:00 - 14:59


# Extract the same time window (Tuesday 12:00 - 14:59) for each week
DataDf$weekday <- c(as.POSIXlt(DataDf$Date, format = "%d/%m/%Y")$wday)
DataDf$Time <- as.POSIXlt(DataDf$Time, format = "%H:%M:%S")
trainingDf <- subset(DataDf, DataDf$weekday == 2 & DataDf$Time >= as.POSIXlt("12:00:00", format = "%H:%M:%S") & DataDf$Time < as.POSIXlt("14:00:00", format = "%H:%M:%S"))



### 2

logLik_global_intensity <- vector(length=14)
bic_global_intensity <- vector(length=14)
n_vector <- rep(180, 52)

for (i in 3:16)
{
  model <- depmix(Global_intensity~1, data = trainingDf, nstates = i, ntimes = n_vector)
  fitModel <- fit(model)
  # summary <- summary(fitModel)
  bic <- BIC(fitModel)
  logLik <- logLik(fitModel)
  logLik_global_intensity[i-2] <- logLik
  bic_global_intensity[i-2] <- bic
}

logLik_active_power <- vector(length=14)
bic_active_power <- vector(length=14)

for (i in 3:16)
{
  model <- depmix(Global_active_power~1, data = trainingDf, nstates = i, ntimes = n_vector)
  fitModel <- fit(model)
  # summary <- summary(fitModel)
  bic <- BIC(fitModel)
  logLik <- logLik(fitModel)
  logLik_active_power[i-2] <- logLik
  bic_active_power[i-2] <- bic
}

logLik_reactive_power <- vector(length=14)
bic_reactive_power <- vector(length=14)

for (i in 3:16)
{
  model <- depmix(Global_reactive_power~1, data = trainingDf, nstates = i, ntimes = n_vector)
  fitModel <- fit(model)
  # summary <- summary(fitModel)
  bic <- BIC(fitModel)
  logLik <- logLik(fitModel)
  logLik_reactive_power[i-2] <- logLik
  bic_reactive_power[i-2] <- bic
}

HMM_Global_Intensity <- data.frame(bic = bic_global_intensity, logLik = logLik_global_intensity)
HMM_Active_Power <- data.frame(bic = bic_active_power, logLik = logLik_active_power)
HMM_Reactive_Power <- data.frame(bic = bic_reactive_power, logLik = logLik_reactive_power)

ggplot(NULL, aes(x=(3:16))) + 
  geom_line(data = HMM_Global_Intensity, aes(y = bic), color = "darkred") + 
  geom_line(data = HMM_Global_Intensity,aes(y = logLik), color="blue")

ggplot(NULL, aes(x=(3:16))) + 
  geom_line(data = HMM_Active_Power, aes(y = bic), color = "darkred") + 
  geom_line(data = HMM_Active_Power,aes(y = logLik), color="blue")

ggplot(NULL, aes(x=(3:16))) + 
  geom_line(data = HMM_Reactive_Power, aes(y = bic), color = "darkred") + 
  geom_line(data = HMM_Reactive_Power,aes(y = logLik), color="blue")
