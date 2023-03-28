### Term Project

library(ggbiplot)
library(zoo)
library(stats)
library(factoextra)
library(depmixS4)

#load("termProject.RData")

DataDf <- read.table("Term_Project_Dataset.txt", header = TRUE, sep = ",")

# linear interpolation on NA values
DataDf <- DataDf[2:nrow(DataDf),]
DataDf[,3:9] <- na.approx(DataDf[,3:9])

# scale train data
DataDf[,3:9] <- scale(DataDf[,3:9])

# split data in train data (80% approx) and test data (20% approx)
TrainDf <- subset(DataDf, as.POSIXlt(DataDf$Date, format = "%d/%m/%Y") < as.POSIXlt("29/4/2009", format = "%d/%m/%Y"))
TestDf <- DataDf[(nrow(TrainDf)+1):nrow(DataDf),]





######### Part 1 #########


# Find the best selection of features using PCA on the train data
TrainPCA <- prcomp(TrainDf[,3:9])
TrainPCA
summary(TrainPCA)

# plot the results (PCs)
fviz_eig(TrainPCA, geom="bar", addlabels=T, width=0.5) +
  ggtitle("Scree Plot") +
  theme(plot.title = element_text(hjust = 0.5))

fviz_pca_var(TrainPCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +     # Avoid text overlapping
            xlab("PC1") +
            ylab("PC2")





######### Part 2 #########

# Divide data into weeks to find time window
TrainDf$week_num <- strftime(as.Date(TrainDf$Date, format = "%d/%m/%Y"), format = "%V")
TrainDf$year <- strftime(as.Date(TrainDf$Date, format = "%d/%m/%Y"), format = "%Y")
WeeksDf <- data.frame(matrix(nrow = 10080, ncol = 0))
for (i in 1:52) 
{
  week <- subset(TrainDf, (TrainDf$week_num == i | TrainDf$week_num == paste0("0", i)) & TrainDf$year == "2007" & TrainDf$Date != "31/12/2007")
  WeeksDf <- cbind(WeeksDf, week$Global_intensity)
  colnames(WeeksDf)[i] <- paste0("week_", i)
}
WeeksDf$Average <- rowMeans(WeeksDf)
WeeksDf$Datetime <- with(TrainDf, as.POSIXct(paste(Date, Time), format="%d/%m/%Y %H:%M:%S"))[21997:32076]

# plot averages to find time window
ggplot(WeeksDf[1440:2880, ]) +
  geom_line(aes(x = Datetime, y = Average, color="Average"), linewidth=0.5) + 
  ggtitle("Average Global Intensity") +
  theme(plot.title = element_text(hjust = 0.5))

# plot individual weeks to find time window
ggplot(WeeksDf[1440:2880, ]) +
  geom_line(aes(x = Datetime, y= week_30, color="Average"), linewidth=0.5) + # Choose Tuesday 18-21
  labs(y="Global_Intensity") +
  ggtitle("Global Intensity for Week 2") +
  theme(plot.title = element_text(hjust = 0.5))




# Extract the same time window (Tuesday 18-21) for each week
TrainDf$weekday <- c(as.POSIXlt(TrainDf$Date, format = "%d/%m/%Y")$wday)
TrainDf$Time <- as.POSIXlt(TrainDf$Time, format = "%H:%M:%S")
TrainDf_window <- subset(TrainDf, TrainDf$weekday == 2 & TrainDf$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & TrainDf$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))

TestDf$weekday <- c(as.POSIXlt(TestDf$Date, format = "%d/%m/%Y")$wday)
TestDf$Time <- as.POSIXlt(TestDf$Time, format = "%H:%M:%S")
TestDf_window <- subset(TestDf, TestDf$weekday == 2 & TestDf$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & TestDf$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))

logLik_train <- vector(length=11)
bic_train <- vector(length=11)
fit_list <- list()

# Train HMM for states 4 to 24
train_vector <- rep(180, 124) # 124 weeks in TrainDf
index <- 1
for (i in seq(4, 24, 2)) { # i increments by 2
  modTrain <- depmix(list(Global_active_power~1, Global_intensity~1), 
                     data=TrainDf_window, nstates = i, ntimes = train_vector, 
                     family=list(gaussian(), gaussian())) 
  fmodTrain <- fit(modTrain)
  fit_list[[i]] <- fmodTrain
  bic <- BIC(fmodTrain)
  logLik <- logLik(fmodTrain)
  logLik_train[index] <- logLik
  bic_train[index] <- bic
  index <- index + 1
}


HMM_train <- data.frame(bic = bic_train, logLik = logLik_train)

save(HMM_train, fit_list, file = "HMM.RData")
# load("HMM.RData")

ggplot(NULL, aes(x=(4:24))) + 
  geom_line(data = HMM_train, aes(x = seq(4,24,2), y = bic, colour = "BIC")) + 
  geom_line(data = HMM_train,aes(x = seq(4,24,2), y = logLik, colour="LogLik")) +
  labs(x="nstates", y="values") +
  ggtitle("LogLik and BIC comparison for n states of HMM") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual("", values = c("BIC"="darkred", 
                                     "LogLik"="steelblue")) 


# Calculate the log-likelihood of the test data for the selected models

test_vector <- rep(180, 30)

# 4 states
params_4states <- getpars(fit_list[[4]])
modTest_4states <- depmix(list(Global_active_power~1, Global_intensity~1), 
                           data=TestDf_window, nstates = 4, ntimes = test_vector, 
                           family=list(gaussian(), gaussian())) 
modTest_4states <- setpars(modTest_4states, params_4states)
test_loglikelihood_4states <- forwardbackward(modTest_4states, TestDf_window, return.all=TRUE, useC=FALSE)$logLike

# 8 states
params_8states <- getpars(fit_list[[8]])
modTest_8states <- depmix(list(Global_active_power~1, Global_intensity~1), 
                           data=TestDf_window, nstates = 8, ntimes = test_vector, 
                           family=list(gaussian(), gaussian())) 
modTest_8states <- setpars(modTest_8states, params_8states)
test_loglikelihood_8states <- forwardbackward(modTest_8states, TestDf_window, return.all=TRUE, useC=FALSE)$logLike

# 10 states
params_10states <- getpars(fit_list[[10]])
modTest_10states <- depmix(list(Global_active_power~1, Global_intensity~1), 
                   data=TestDf_window, nstates = 10, ntimes = test_vector, 
                   family=list(gaussian(), gaussian())) 
modTest_10states <- setpars(modTest_10states, params_10states)
test_loglikelihood_10states <- forwardbackward(modTest_10states, TestDf_window, return.all=TRUE, useC=FALSE)$logLike


# normalize the loglikelihood
logLik_train_normalized <- logLik_train / nrow(TrainDf_window)
test_loglikelihood_4states_normalized <- test_loglikelihood_4states / nrow(TestDf_window)
test_loglikelihood_8states_normalized <- test_loglikelihood_8states / nrow(TestDf_window)
test_loglikelihood_10states_normalized <- test_loglikelihood_10states / nrow(TestDf_window)

ggplot(NULL, aes(x=(4:24))) + 
  geom_line(data = HMM_train, aes(x = seq(4,24,2), y = logLik_train_normalized, colour = "Train logLik")) + 
  geom_point(aes(x = 4, y = test_loglikelihood_4states_normalized, colour="Test logLik")) +     # underfit
  geom_point(aes(x = 8, y = test_loglikelihood_8states_normalized, colour="Test logLik")) + 
  geom_point(aes(x = 10, y = test_loglikelihood_10states_normalized, colour="Test logLik")) +   # overfit
  labs(x="nstates", y="values") +
  ggtitle("Train LogLik and Test LogLik comparison for n states of HMM") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual("", values = c("Train logLik"="darkred", 
                                     "Test logLik"="steelblue")) 


######### Part 3 #########

DataWithAnomaliesDf1 <- read.table("Data_with_Anomalies/Dataset_with_Anomalies_1.txt", header = TRUE, sep = ",")
DataWithAnomaliesDf2 <- read.table("Data_with_Anomalies/Dataset_with_Anomalies_2.txt", header = TRUE, sep = ",")
DataWithAnomaliesDf3 <- read.table("Data_with_Anomalies/Dataset_with_Anomalies_3.txt", header = TRUE, sep = ",")

# linear interpolation on NA values
DataWithAnomaliesDf1[,3:9] <- na.approx(DataWithAnomaliesDf1[,3:9])
DataWithAnomaliesDf2[,3:9] <- na.approx(DataWithAnomaliesDf2[,3:9])
DataWithAnomaliesDf3[,3:9] <- na.approx(DataWithAnomaliesDf3[,3:9])

# scale train data
DataWithAnomaliesDf1[,3:9] <- scale(DataWithAnomaliesDf1[,3:9])
DataWithAnomaliesDf2[,3:9] <- scale(DataWithAnomaliesDf2[,3:9])
DataWithAnomaliesDf3[,3:9] <- scale(DataWithAnomaliesDf3[,3:9])

# Extract the same time window (Tuesday 18-21) for each week
DataWithAnomaliesDf1$weekday <- c(as.POSIXlt(DataWithAnomaliesDf1$Date, format = "%d/%m/%Y")$wday)
DataWithAnomaliesDf1$Time <- as.POSIXlt(DataWithAnomaliesDf1$Time, format = "%H:%M:%S")
DataWithAnomaliesDf1 <- subset(DataWithAnomaliesDf1, DataWithAnomaliesDf1$weekday == 2 & DataWithAnomaliesDf1$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & DataWithAnomaliesDf1$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))

DataWithAnomaliesDf2$weekday <- c(as.POSIXlt(DataWithAnomaliesDf2$Date, format = "%d/%m/%Y")$wday)
DataWithAnomaliesDf2$Time <- as.POSIXlt(DataWithAnomaliesDf2$Time, format = "%H:%M:%S")
DataWithAnomaliesDf2 <- subset(DataWithAnomaliesDf2, DataWithAnomaliesDf2$weekday == 2 & DataWithAnomaliesDf2$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & DataWithAnomaliesDf2$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))

DataWithAnomaliesDf3$weekday <- c(as.POSIXlt(DataWithAnomaliesDf3$Date, format = "%d/%m/%Y")$wday)
DataWithAnomaliesDf3$Time <- as.POSIXlt(DataWithAnomaliesDf3$Time, format = "%H:%M:%S")
DataWithAnomaliesDf3 <- subset(DataWithAnomaliesDf3, DataWithAnomaliesDf3$weekday == 2 & DataWithAnomaliesDf3$Time >= as.POSIXlt("18:00:00", format = "%H:%M:%S") & DataWithAnomaliesDf3$Time < as.POSIXlt("21:00:00", format = "%H:%M:%S"))


anomaly_data_vector = rep(180, 52)
anomaly_data_params <- getpars(modTest_10states)

# compute the log-likelihood for the data with anomalies
# Data with Anomalies 1
modAnomalyData1 <- depmix(list(Global_active_power~1, Global_intensity~1), 
                          data=DataWithAnomaliesDf1, nstates = 10, ntimes = anomaly_data_vector, 
                          family=list(gaussian(), gaussian())) 
modAnomalyData1 <- setpars(modAnomalyData1, anomaly_data_params)
AnomalyData_loglikelihood1 <- forwardbackward(modAnomalyData1, DataWithAnomaliesDf1, return.all=TRUE, useC=FALSE)$logLike

# Data with Anomalies 2
modAnomalyData2 <- depmix(list(Global_active_power~1, Global_intensity~1), 
                          data=DataWithAnomaliesDf2, nstates = 10, ntimes = anomaly_data_vector, 
                          family=list(gaussian(), gaussian())) 
modAnomalyData2 <- setpars(modAnomalyData2, anomaly_data_params)
AnomalyData_loglikelihood2 <- forwardbackward(modAnomalyData2, DataWithAnomaliesDf2, return.all=TRUE, useC=FALSE)$logLike

# Data with Anomalies 3
modAnomalyData3 <- depmix(list(Global_active_power~1, Global_intensity~1), 
                          data=DataWithAnomaliesDf3, nstates = 10, ntimes = anomaly_data_vector, 
                          family=list(gaussian(), gaussian())) 
modAnomalyData3 <- setpars(modAnomalyData3, anomaly_data_params)
AnomalyData_loglikelihood3 <- forwardbackward(modAnomalyData3, DataWithAnomaliesDf3, return.all=TRUE, useC=FALSE)$logLike

# normalize the loglikelihoods

AnomalyData_loglikelihood1_normalized <- AnomalyData_loglikelihood1 / nrow(DataWithAnomaliesDf1)
AnomalyData_loglikelihood2_normalized <- AnomalyData_loglikelihood2 / nrow(DataWithAnomaliesDf2)
AnomalyData_loglikelihood3_normalized <- AnomalyData_loglikelihood3 / nrow(DataWithAnomaliesDf3)

