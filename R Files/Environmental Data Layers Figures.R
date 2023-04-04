#import necessary libraries
library(raster)
library(reshape)
library(ggplot2)
library(dplyr)
library(viridis)
source("R Files/Helper Functions.R")

#grab list of desired BNAM files
full_file_list <- list.files("C:/Users/natha/Research/Data/BNAM/BNAM_1990-2019_monthly", full.names = TRUE, recursive = TRUE)
BtmTempFiles <- full_file_list[grepl("BtmTemp",full_file_list)==TRUE&grepl(".asc",full_file_list)==TRUE]
BtmSalinityFiles <- full_file_list[grepl("BtmSalinity",full_file_list)==TRUE&grepl(".asc",full_file_list)==TRUE]
BtmStressFiles <- full_file_list[grepl("BtmStress",full_file_list)==TRUE&grepl(".asc",full_file_list)==TRUE]

########################################TIME SERIES OF CERTAIN LOCATIONS#######################

######TEMPERATURE########

#extract points of interest for bottom temperature
BNAM_df <- extract(stack(BtmTempFiles), 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)

#set up a dataframe of the time series
BanquereauBank <- order_BNAM_values(BNAM_df[1,])
BrownsBank <- order_BNAM_values(BNAM_df[2,])
SableBank <- order_BNAM_values(BNAM_df[3,])
NorthOfCapeBreton <- order_BNAM_values(BNAM_df[4,])
monthyears <- BNAM_month_year_order()
data <- cbind.data.frame(monthyears, BanquereauBank, BrownsBank, 
                         SableBank, NorthOfCapeBreton)
data$monthyears <- as.Date(data$monthyears, format = "%d %b %Y")
data <- data.frame(data[1], stack(data[2:ncol(data)]))

#make temperature plots
ggplot(data=data, mapping=aes(x = monthyears, y=values, col = ind)) +
  xlab("Date") + ylab("Bottom Temperature (°C)") +
  geom_point() + geom_line(size = 1) + scale_x_date(date_labels = "%b %Y", breaks = data$monthyears[seq(1, 240, by = 6)]) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Temperature Time Series (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location")

#######SALINITY######

#extract points of interest for bottom salinity
BNAM_df <- extract(stack(BtmSalinityFiles), 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)

#set up a dataframe of the time series
BanquereauBank <- order_BNAM_values(BNAM_df[1,])
BrownsBank <- order_BNAM_values(BNAM_df[2,])
SableBank <- order_BNAM_values(BNAM_df[3,])
NorthOfCapeBreton <- order_BNAM_values(BNAM_df[4,])
monthyears <- BNAM_month_year_order()
data <- cbind.data.frame(monthyears, BanquereauBank, BrownsBank, 
                         SableBank, NorthOfCapeBreton)
data$monthyears <- as.Date(data$monthyears, format = "%d %b %Y")
data <- data.frame(data[1], stack(data[2:ncol(data)]))

#make salinity plots
ggplot(data=data, mapping=aes(x = monthyears, y=values, col = ind)) +
  xlab("Date") + ylab("Bottom Salinity (psu)") +
  geom_point() + geom_line(size = 1) + scale_x_date(date_labels = "%b %Y", breaks = data$monthyears[seq(1, 240, by = 6)]) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Salinity Time Series (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location")

#######STRESS######

#extract points of interest for bottom stress
BNAM_df <- extract(stack(BtmStressFiles), 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)

#set up a dataframe of the time series
BanquereauBank <- order_BNAM_values(BNAM_df[1,])
BrownsBank <- order_BNAM_values(BNAM_df[2,])
SableBank <- order_BNAM_values(BNAM_df[3,])
NorthOfCapeBreton <- order_BNAM_values(BNAM_df[4,])
monthyears <- BNAM_month_year_order()
data <- cbind.data.frame(monthyears, BanquereauBank, BrownsBank, 
                         SableBank, NorthOfCapeBreton)
data$monthyears <- as.Date(data$monthyears, format = "%d %b %Y")
data <- data.frame(data[1], stack(data[2:ncol(data)]))

#make stress plots
ggplot(data=data, mapping=aes(x = monthyears, y=values, col = ind)) +
  xlab("Date") + ylab(expression("Bottom Stress (kg*m"^(-1)*"s"^(-2)*")")) +
  geom_point() + geom_line(size = 1) + scale_x_date(date_labels = "%b %Y", breaks = data$monthyears[seq(1, 240, by = 6)]) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Bottom Stress Time Series (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location")

###########################MONTHLY AVG/VARIABILITY PLOTS##########################################

#plot the mean time series for btmtemp
stack <- monthly_average_stack(file_subset = BtmTempFiles)
stack_df <- as.data.frame(stack[[1]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "",
                                                colours = viridis(100), limits = c(-5,18)) + coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle("Mean BNAM Bottom Temperature (°C) (2000-2019)") + ylab("Latitude") + xlab("Longitude")
#plot the SD time series for btm salinity
stack_df <- as.data.frame(stack[[2]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "",
                                                colours = viridis(100)) + coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle("Standard Deviation for BNAM Bottom Temperature (°C) (2000-2019)") + ylab("Latitude") + xlab("Longitude")

#plot the mean time series for btmsalinity
stack <- monthly_average_stack(file_subset = BtmSalinityFiles)
stack_df <- as.data.frame(stack[[1]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(26,36)) +
  coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle("Mean BNAM Bottom Salinity (psu) (2000-2019)") + ylab("Latitude") + xlab("Longitude")
#plot the SD time series for btmsalinity
stack_df <- as.data.frame(stack[[2]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,2)) +
  coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle("Standard Deviation for BNAM Bottom Salinity (psu) (2000-2019)") + ylab("Latitude") + xlab("Longitude")

#plot the mean time series for btmstress
stack <- monthly_average_stack(file_subset = BtmStressFiles)
stack_df <- as.data.frame(stack[[1]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0, 0.15)) +
  coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle(expression("Mean BNAM Bottom Stress (kg*m"^(-1)*"s"^(-2)*") (2000-2019)")) + ylab("Latitude") + xlab("Longitude")
#plot the SD time series for btmstress
stack_df <- as.data.frame(stack[[2]], xy = TRUE) %>% melt(id.vars = c('x','y'))
levels(stack_df$variable) <- month.name[1:12]
ggplot() + geom_raster(data = stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,0.05)) +
  coord_equal() + xlim(-70,-56) + ylim(42,48) +
  ggtitle(expression("Standard Deviation for BNAM Bottom Stress (kg*m"^(-1)*"s"^(-2)*") (2000-2019)")) + 
  ylab("Latitude") + xlab("Longitude")

###################################################SEASONAL TREND PLOTS FOR FOUR LOCATIONS#########################

######TEMPERATURE######

#extract the average for each month at each location and put into dataframe
stack <- monthly_average_stack(file_subset = BtmTempFiles)
BNAM_df <- extract(stack[[1]], 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)
Locations <- rep(c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton"), 12)
BNAM_df <- data.frame(Locations, stack(BNAM_df[1:ncol(BNAM_df)]))
levels(BNAM_df$ind) <- month.name
#also attach the standard deviations
sd_df <- extract(stack[[2]], 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
sd_df <- as.data.frame(sd_df)
BNAM_df$sd <- stack(sd_df[1:ncol(sd_df)])[,1]

#make temperature plots
ggplot(data=BNAM_df, mapping=aes(x = ind, y=values, col = Locations)) +
  xlab("Month") + ylab("Value (°C)") +
  geom_point(size = 2) + geom_line(size = 1, group = 1) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Bottom Temperature - Mean and Standard Deviation (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location") + geom_errorbar(mapping = aes(x = ind, ymin = values - sd, ymax = values + sd)) +
  facet_wrap(~Locations) + theme(legend.position = "none")

########SALINITY #########

#extract the average for each month at each location and put into dataframe
stack <- monthly_average_stack(file_subset = BtmSalinityFiles)
BNAM_df <- extract(stack[[1]], 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)
Locations <- rep(c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton"), 12)
BNAM_df <- data.frame(Locations, stack(BNAM_df[1:ncol(BNAM_df)]))
levels(BNAM_df$ind) <- month.name
#also attach the standard deviations
sd_df <- extract(stack[[2]], 
                 rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                 method = "bilinear")
sd_df <- as.data.frame(sd_df)
BNAM_df$sd <- stack(sd_df[1:ncol(sd_df)])[,1]

#make salinity plots
ggplot(data=BNAM_df, mapping=aes(x = ind, y=values, col = Locations)) +
  xlab("Month") + ylab("Value (psu)") +
  geom_point(size = 2) + geom_line(size = 1, group = 1) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Bottom Salinity - Mean and Standard Deviation (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location") + geom_errorbar(mapping = aes(x = ind, ymin = values - sd, ymax = values + sd)) +
  facet_wrap(~Locations) + theme(legend.position = "none")

########STRESS #########

#extract the average for each month at each location and put into dataframe
stack <- monthly_average_stack(file_subset = BtmStressFiles)
BNAM_df <- extract(stack[[1]], 
                   rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                   method = "bilinear")
BNAM_df <- as.data.frame(BNAM_df)
Locations <- rep(c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton"), 12)
BNAM_df <- data.frame(Locations, stack(BNAM_df[1:ncol(BNAM_df)]))
levels(BNAM_df$ind) <- month.name
#also attach the standard deviations
sd_df <- extract(stack[[2]], 
                 rbind(c(-58,44.5),c(-66, 42.65), c(-60.95, 44), c(-60, 46.5)), 
                 method = "bilinear")
sd_df <- as.data.frame(sd_df)
BNAM_df$sd <- stack(sd_df[1:ncol(sd_df)])[,1]

#make stress plots
ggplot(data=BNAM_df, mapping=aes(x = ind, y=values, col = Locations)) +
  xlab("Month") + ylab(expression("Value (kg*m"^(-1)*"s"^(-2)*")")) +
  geom_point(size = 2) + geom_line(size = 1, group = 1) +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 16), 
        text=element_text(family="serif", size = 14))  + 
  ggtitle ("BNAM Bottom Stress - Mean and Standard Deviation (2000-2019)") + 
  scale_color_discrete(labels = c("Banquereau Bank","Browns Bank","Sable Bank","North of Cape Breton")) + 
  labs(col = "Location") + geom_errorbar(mapping = aes(x = ind, ymin = values - sd, ymax = values + sd)) +
  facet_wrap(~Locations) + theme(legend.position = "none")

##############################Plot environmental statistics (min, max, range, etc.)##########

#create stacks for 2000-2019 bottom temperature, stress and salinity
BtmTempStack <- stack(BtmTempFiles[121:360])
BtmSalinityStack <- stack(BtmSalinityFiles[121:360])
BtmStressStack <- stack(BtmStressFiles[121:360])
#crop them to spatial domain
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")
BtmTempStack <- mask(crop(BtmTempStack, shp), shp)
BtmSalinityStack <- mask(crop(BtmSalinityStack, shp), shp)
BtmStressStack <- mask(crop(BtmStressStack, shp), shp)

#calculate the mean, min, max, and range for bottom temperature
AverageTemp <- calc(BtmTempStack, fun = mean)
MinTemp <- calc(BtmTempStack, fun = min)
MaxTemp <- calc(BtmTempStack, fun = max)
RangeTemp <- MaxTemp - MinTemp
#calculate the mean, min, max, and range for bottom salinity
AverageSalinity <- calc(BtmSalinityStack, fun = mean)
MinSalinity <- calc(BtmSalinityStack, fun = min)
MaxSalinity <- calc(BtmSalinityStack, fun = max)
RangeSalinity <- MaxSalinity - MinSalinity
#calculate the mean, min, max, and range for bottom stress
AverageStress <- calc(BtmStressStack, fun = mean)
MinStress <- calc(BtmStressStack, fun = min)
MaxStress <- calc(BtmStressStack, fun = max)
RangeStress <- MaxStress - MinStress

#grab the average min, max, and range for bottom temperature
years <- seq(2000, 2019, by = 1)
average_min_max_range_temp <- average_min_max_range_calc(BtmTempFiles[121:360], years)
AverageMinTemp <- mask(crop(average_min_max_range_temp[[1]], shp), shp)
AverageMaxTemp <- mask(crop(average_min_max_range_temp[[2]], shp), shp)
AverageRangeTemp <- mask(crop(average_min_max_range_temp[[3]], shp), shp)

#grab the average min, max, and range for bottom salinity
average_min_max_range_salinity <- average_min_max_range_calc(BtmSalinityFiles[121:360], years)
AverageMinSalinity <- mask(crop(average_min_max_range_salinity[[1]], shp), shp)
AverageMaxSalinity <- mask(crop(average_min_max_range_salinity[[2]], shp), shp)
AverageRangeSalinity <- mask(crop(average_min_max_range_salinity[[3]], shp), shp)

#grab the average min, max, and range for bottom stress
average_min_max_range_stress <- average_min_max_range_calc(BtmStressFiles[121:360], years)
AverageMinStress <- mask(crop(average_min_max_range_stress[[1]], shp), shp)
AverageMaxStress <- mask(crop(average_min_max_range_stress[[2]], shp), shp)
AverageRangeStress <- mask(crop(average_min_max_range_stress[[3]], shp), shp)

#plot overall means for temperature, salinity, and stress
par(mfrow = c(2,2))
plot(AverageTemp, col = viridis(100), main = "Average Bottom Temperature (°C) (2000-2019)",
     family = "serif", cex.axis = 1.5, cex.lab = 1.5, 
     xlab = "Lon", ylab = "Lat", cex.main = 1.5, axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageSalinity, col = viridis(100), main = "Average Bottom Salinity (psu) (2000-2019)",
     family = "serif", cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5, xlab = "Lon", ylab = "Lat"
     , axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageStress, col = viridis(100), main = substitute(paste(bold("Average Bottom Stress (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))),
     family = "serif", cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5, xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))

#plot temp min, max, range (and averages of them)
par(mfrow = c(2,3))
plot(MinTemp, main = "Bottom Temperature Minimum  (°C) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100), xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(MaxTemp, main = "Bottom Temperature Maximum (°C) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(RangeTemp, main = "Bottom Temperature Range (°C) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMinTemp, main = "Bottom Temperature Average\nMinimum (°C) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMaxTemp, main = "Bottom Temperature Average\nMaximum (°C) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageRangeTemp, main = "Bottom Temperature Average\nRange (°C) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))

#plot salinity min, max, range (and averages of them)
par(mfrow = c(2,3))
plot(MinSalinity, main = "Bottom Salinity Minimum  (psu) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100), xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(MaxSalinity, main = "Bottom Salinity Maximum (psu) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(RangeSalinity, main = "Bottom Salinity Range (psu) (2000-2019)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMinSalinity, main = "Bottom Salinity Average\nMinimum (psu) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMaxSalinity, main = "Bottom Salinity Average\nMaximum (psu) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageRangeSalinity, main = "Bottom Salinity Average\nRange (psu) (2000-2019)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))

#plot stress min, max, range (and averages of them)
par(mfrow = c(2,3))
plot(MinStress, main = substitute(paste(bold("Bottom Stress Minimum (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100), xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(MaxStress, main = substitute(paste(bold("Bottom Stress Maximum (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(RangeStress, main = substitute(paste(bold("Bottom Stress Range (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMinStress, main = substitute(paste(bold("Bottom Stress Average Minimum (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageMaxStress, main = substitute(paste(bold("Bottom Stress Average Maximum (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(AverageRangeStress, main = substitute(paste(bold("Bottom Stress Average Range (kg*m"^(-1)*"s"^(-2)*") (2000-2019)"))), family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))

################################################QUICK DEM PLOTS#######################################

###LOAD IN DEM DATA (AND DERIVED LAYERS)###
DEM = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
Easterness = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/east_bathych/w001001.adf")
Northerness = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/nort_bathych/w001001.adf")
Slope = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/slop_bathych/w001001.adf")
RDMV = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/RDMV_bathych/w001001.adf")

#crop to spatial domain
DEM <- mask(crop(DEM, shp), shp)
Easterness <- mask(crop(Easterness, shp), shp)
Northerness <- mask(crop(Northerness, shp), shp)
Slope <- mask(crop(Slope, shp), shp)
RDMV <- mask(crop(RDMV, shp), shp)

#plot DEM layers
par(mfrow = c(2,3))
plot(log(-DEM), main = "Log(Depth [m])", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100), xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(sqrt(Slope), main = "sqrt(Slope [°])", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(Easterness, main = "Easterness (Unitless)", family = "serif", cex.axis = 1.5,
     cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat", 
     axis.args=list(family = "serif", cex.axis=1.5))
plot(Northerness, main = "Northerness (Unitless)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))
plot(RDMV, main = "Relative Deviation from Mean Value (RDMV) (Unitless)", family = "serif", 
     cex.axis = 1.5, cex.main = 1.7, cex.lab = 1.5, col = viridis(100),  xlab = "Lon", ylab = "Lat",
     axis.args=list(family = "serif", cex.axis=1.5))