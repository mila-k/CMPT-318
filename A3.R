### Assignment 3

library(depmixS4)
library(zoo)
library(ggplot2)

DataDf <- read.table("Group_Assignment_3_Dataset.txt", header = TRUE, sep = ",")

# linear interpolation on NA values
DataDf$Global_active_power <- na.approx(DataDf$Global_active_power)
DataDf$Global_reactive_power <- na.approx(DataDf$Global_reactive_power)
DataDf$Global_intensity <- na.approx(DataDf$Global_intensity)




### 1

# Scale the dataset
DataDf$Global_active_power <- scale(DataDf$Global_active_power, center = TRUE, scale = TRUE)
DataDf$Global_reactive_power <- scale(DataDf$Global_reactive_power, center = TRUE, scale = TRUE)
DataDf$Global_intensity <- scale(DataDf$Global_intensity, center = TRUE, scale = TRUE)



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
WeeksDf$Average <- rowMeans(WeeksDf)
WeeksDf$Datetime <- with(DataDf, as.POSIXct(paste(Date, Time), format="%d/%m/%Y %H:%M:%S"))[1:10080]

# plot averages
ggplot(WeeksDf[1440:2880, ]) +
  geom_line(aes(x = Datetime, y= Average, color="Average"), linewidth=0.5) + 
  ggtitle("Average Global Intensity") +
  theme(plot.title = element_text(hjust = 0.5))

# plot individual weeks
ggplot(WeeksDf[1440:2880, ]) +
  geom_line(aes(x = Datetime, y= week_30, color="Average"), linewidth=0.5) + # Choose Tuesday 18-21
  labs(y="Global_Intensity") +
  ggtitle("Global Intensity for Week 2") +
  theme(plot.title = element_text(hjust = 0.5))

# Extract the same time window (Tuesday 18-21) for each week
DataDf$weekday <- c(as.POSIXlt(DataDf$Date, format = "%d/%m/%Y")$wday)
DataDf$Time <- as.POSIXlt(DataDf$Time, format = "%H:%M:%S")
trainingDf <- subset(DataDf, DataDf$weekday == 2 & DataDf$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & DataDf$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))



### 2

logLik_global_intensity <- vector(length=14)
bic_global_intensity <- vector(length=14)
n_vector <- rep(180, 52)

# HMM for global intensity
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

HMM_Global_Intensity <- data.frame(bic = bic_global_intensity, logLik = logLik_global_intensity)

ggplot(NULL, aes(x=(3:16))) + 
  geom_line(data = HMM_Global_Intensity, aes(y = bic, colour = "BIC")) + 
  geom_line(data = HMM_Global_Intensity,aes(y = logLik, colour="LogLik")) +
  labs(x="nstates", y="values") +
  ggtitle("LogLik and BIC comparison for n states of HMM (Global Intensity)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual("", values = c("BIC"="darkred", 
                                     "LogLik"="steelblue")) 
