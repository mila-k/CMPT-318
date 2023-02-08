## Assignment 2 Part 2

library(TTR)
library(moments)
library(dplyr)
library(ggplot2)  
library(zoo)

DataDf <- read.table("Group_Assignment_2_Dataset.txt", header = TRUE, sep = ",")
DataDf$week_num <- strftime(as.Date(DataDf$Date, format = "%d/%m/%Y"), format = "%V")
DataDf <- subset(DataDf, DataDf$Date != "31/12/2007") # drop 31/12/2007 (Monday) to get complete 52 weeks 

# linear interpolation on NA values
DataDf$Global_intensity <- na.approx(DataDf$Global_intensity)


# Moving Averages
smoothenedWeeksDf <- data.frame(matrix(nrow = 10080, ncol = 0))
for (i in 1:52) 
{
  week <- subset(DataDf, DataDf$week_num == i | DataDf$week_num == paste0("0", i))
  moving_average <- SMA(week$Global_intensity, 10)
  smoothenedWeeksDf <- cbind(smoothenedWeeksDf, moving_average)
  colnames(smoothenedWeeksDf)[i] <- paste0("smoothened_week_", i)
}

# average smoothened week
smoothenedWeeksDf <- tail(smoothenedWeeksDf, -9) # remove the first 9 NA values produced from SMA
smoothenedWeeksDf$average_smoothened_week <- rowMeans(smoothenedWeeksDf)


# Euclidean distance (to measure similarity between the smoothened week and average)
euclidean <- function(a, b) sqrt(sum((a - b)^2))


# find euclidean distance for every week
euclideanDistances <- vector(length=52)
for (i in 1:52) 
{
  col <- smoothenedWeeksDf %>% select(paste0("smoothened_week_", i))
  euclideanDistances[i] <- euclidean(col, smoothenedWeeksDf$average_smoothened_week)
}


# Most anomalous week (Max euclidean distance)
mostAnomalous <- which.max(euclideanDistances)   # week 52

# Least anomalous week (Min euclidean distance)
leastAnomalous <- which.min(euclideanDistances)   # week 30


# store scores in a dataframe
scoresDf <- data.frame(week = 1:52, score = euclideanDistances)

# add datetime column for the x axis
smoothenedWeeksDf$Datetime <- with(DataDf, as.POSIXct(paste(Date, Time), format="%d/%m/%Y %H:%M:%S"))[10:10080]


# plot the smoothened versions of the most and the least anomalous weeks against the average smoothened week

ggplot(smoothenedWeeksDf) +
  geom_point(aes(x = Datetime, y= smoothened_week_30, colour="Least Anomalous"), size=0.5) +
  geom_line(aes(x = Datetime, y= average_smoothened_week, colour="Average"), linewidth=1.0) +
  scale_colour_manual("", values = c("Least Anomalous"="darkred", 
                                     "Average"="yellow")) +
  labs(x="Time", y="Global Intensity") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Least Anomalous Week")

ggplot(smoothenedWeeksDf) +
  geom_point(aes(x = Datetime, y= smoothened_week_52, colour="Most Anomalous"), size=0.5) +
  geom_line(aes(x = Datetime, y= average_smoothened_week, colour="Average"), linewidth=1.0) +
  scale_colour_manual("", values = c("Most Anomalous"="steelblue",
                                     "Average"="yellow")) +
  labs(x="Time", y="Global Intensity") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Most Anomalous Week")
