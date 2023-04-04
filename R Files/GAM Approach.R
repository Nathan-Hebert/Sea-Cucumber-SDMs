#load necessary libraries
library(mgcv)
library(ggplot2)
library(knitr)
library(dplyr)
library(sp)
library(reshape)
library(ggpubr)
library(gtools)
library(raster)
library(gratia)
library(viridis)
source("R Files/Helper Functions.R")

#increase memory limit
memory.limit(size = 50000)

#load in data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))
#if RDMV is NA, set to zero
survey_data_combined$DEM_RDMV <- ifelse(is.na(survey_data_combined$DEM_RDMV), 0, survey_data_combined$DEM_RDMV)
#transform some covariates
covariates <- c("BtmTempBNAM","BtmSalinityBNAM","BtmStressBNAMLog","RangeTemp","RangeStressLog","DEM_log","sqrt_DEM_Slope","DEM_Easterness",
                "DEM_Northerness","DEM_RDMV")
survey_data_combined$RangeStressLog <- log(survey_data_combined$RangeStress)
survey_data_combined$BtmStressBNAMLog <- log(survey_data_combined$BtmStressBNAM)
survey_data_combined$sqrt_DEM_Slope <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$DEM_log <- log(-survey_data_combined$DEM)
#scale covariates using all data
survey_data_combined[paste(covariates,"Scaled", sep="")] <- scale(survey_data_combined[covariates])

#grab NS, New Brunswick, PEI, NFLD maps for later plots
map <- raster::getData(country = "CAN", level = 1)
map <- map[which(map$NAME_1 == "Nova Scotia"|map$NAME_1 == "Prince Edward Island"|map$NAME_1 == "New Brunswick"|map$NAME_1 == "Newfoundland and Labrador"),]
map <- spTransform(map, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#grab Maine map for later plots
map1 <- raster::getData(country = "USA", level = 1)
map1 <- map1[which(map1$NAME_1 == "Maine"),]
map1 <- spTransform(map1, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

#survey_data_combined <- survey_data_combined[-which(survey_data_combined$season_winter == 1),] #remove RV winter survey tows

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")
shp <- spTransform(shp, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

##########################################PRESENCE MODEL###########################################

#fit the chosen model on all data
start_time <- Sys.time()
presence_model <- gam(presence~te(UTMX, UTMY, year, k = c(400, 5), 
                                  bs = c("ds","tp"), m = list(c(1,0.5),NA), d = c(2,1)) + 
                        s(DEM_logScaled, k = 3)+s(BtmTempBNAMScaled, k = 3) + s(RangeTempScaled, k = 3) + s(BtmSalinityBNAMScaled, k = 3) + 
                        s(BtmStressBNAMLogScaled, k = 3) + s(RangeStressLogScaled, k = 3) + s(sqrt_DEM_SlopeScaled, k = 3) + 
                        s(DEM_EasternessScaled, DEM_NorthernessScaled, k = 9) + s(DEM_RDMVScaled, k = 3) + as.numeric(snowcrab), 
                      family = "binomial", 
                      data = survey_data_combined, 
                      method = "REML") 
end_time <- Sys.time()
fit_time <- end_time - start_time #get time to fit

#plot the spatio-temporal smooth 
draw(presence_model, select = "te(UTMX,UTMY,year)", n_3d = 20, contour = F)+
  geom_polygon(data = fortify(shp), aes(group = group, x = long, y = lat), fill = NA, col = "white") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  coord_cartesian(ylim = c(4700, 5200), xlim = c(175, 925)) + 
  scale_fill_gradientn(colours = viridis(100)) + ggtitle("f(Eastings, Northings, Year)") + 
  xlab("Eastings") + ylab("Northings") + labs(fill = "Partial\nEffect")
#plot the other model smooths
par(mfrow = c(3,3), mar = c(5,5,4,2))
plot(presence_model, select = 2, shade = TRUE, unconditional = T, se = 1.96, rug = T, main = "p-value < 0.001 ***", cex.main = 2,
     shade.col = "light blue", xlab = "Log(Depth)", ylab = "f(Log(Depth))", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-6,1.2))
plot(presence_model, select = 3, shade = TRUE, unconditional = T, se = 1.96, rug = T, ylim = c(-3,1.2), 
     shade.col = "light blue", xlab = "Bottom Temperature", main = "p-value = 0.023 ***", cex.main = 2,
     ylab = "f(Bottom Temperature)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(presence_model, select = 4, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-3,1.2), rug = T, main = "p-value = 0.020 ***", cex.main = 2,
     shade.col = "light blue", xlab = "Btm. Temperature Range", 
     ylab = "f(Btm. Temperature Range)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(presence_model, select = 5, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-3,1.2), rug = T, main = "p-value = 0.023 ***", cex.main = 2,
     shade.col = "light blue", xlab = "Bottom Salinity", 
     ylab = "f(Bottom Salinity)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(presence_model, select = 6, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-3,1.2), rug = T,
     shade.col = "light blue", xlab = "Log(Bottom Stress)", main = "p-value = 0.002 ***", cex.main = 2,
     ylab = "f(Log(Bottom Stress))", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(presence_model, select = 7, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-3,1.2), rug = T,
     shade.col = "light blue", xlab = "Log(Btm. Stress Range)", main = "p-value = 0.035 ***", cex.main = 2,
     ylab = "f(Log(Btm. Stress Range))", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(presence_model, select = 8, shade = TRUE, unconditional = T, se = 1.96, rug = T, main = "p-value = 0.067", cex.main = 2,
     shade.col = "light blue", xlab = "Square Root of Slope", ylab = "f(Square Root of Slope)", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1.2))
plot(presence_model, select = 9, shade = TRUE, unconditional = T, se = 1.96, rug = T, main = "p-value = 0.062", 
     family = "serif", cex.main = 2, cex.axis = 2, cex.lab = 1.85, scheme = 2, xlab = "Easterness", ylab = "Northerness")
plot(presence_model, select = 10, shade = TRUE, unconditional = T, se = 1.96, rug = T, main = "p-value = 0.381", cex.main = 2,
     shade.col = "light blue", xlab = "RDMV", ylab = "f(RDMV)", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1.2))

#ready information to put into gam summary table for presence model
a <- data.frame(format(round(summary(presence_model)$p.table, 3), digits = 2))
colnames(a) <- c("Estimate", "Std. Error", "Z","P-value")
rownames(a) <- c("Intercept", "Snow Crab Survey")
b <- data.frame(format(round(summary(presence_model)$s.table, 3), digits = 3))
colnames(b) <- c("EDF","Ref. DF","Chi-sq.","P-value")
rownames(b) <- c("f(Eastings, Northings, Year)", "f(Log(Depth))", "f(Bottom Temperature)", "f(Btm. Temperature Range)", 
                 "f(Bottom Salinity)", "f(Log(Bottom Stress))", "f(Log(Btm. Stress Range))", "f(Sqrt(Slope))",
                 "f(Easterness, Northerness)", "f(RDMV)")
#make the table
a <- kable(list(a,b), format = "latex", caption = "Presence Model", align = "rc")

#compute the residuals
survey_data_combined["residuals_presence"] <- residuals(presence_model, type = "response")

#geographically plot the training data's residuals on the initial presence model
a <- ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_presence),
                           data = survey_data_combined) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Response\nResidual\n") + theme(text=element_text(size=16, family="serif")) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot the training data's residuals on the final presence model
b <- ggplot(aes(x = year, y = residuals_presence), data = survey_data_combined) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Response Residual") + xlab("Year") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,0.5)) #plot together

#################################################CPUE MODEL###########################################

#grab presence tows and set up response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),]
survey_data_combined_CPUE$std.WGT_log <- log(survey_data_combined_CPUE$std.WGT)

#fit the chosen model with all the training data
start_time <- Sys.time()
CPUE_model <- gam(log(std.WGT) ~ te(UTMX, UTMY, year, k = c(300, 5), 
                                      bs = c("ds","tp"), m = list(c(1,0.5),NA), d = c(2,1)) + 
                      s(DEM_logScaled, k = 3)+s(BtmTempBNAMScaled, k = 3) + s(RangeTempScaled, k = 3) + s(BtmSalinityBNAMScaled, k = 3) + 
                      s(BtmStressBNAMLogScaled, k = 3) + s(RangeStressLogScaled, k = 3) + s(sqrt_DEM_SlopeScaled, k = 3) + 
                      s(DEM_EasternessScaled, DEM_NorthernessScaled, k = 9) + s(DEM_RDMVScaled, k = 3) + as.numeric(snowcrab), 
                    family = "gaussian", 
                    data = survey_data_combined_CPUE, 
                    method = "REML") 
end_time <- Sys.time()
fit_time3 <- end_time - start_time #get time to fit

#plot the spatio-temporal smooth 
draw(CPUE_model, select = "te(UTMX,UTMY,year)", n_3d = 20, contour = F)+
  geom_polygon(data = fortify(shp), aes(group = group, x = long, y = lat), fill = NA, col = "white") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  coord_cartesian(ylim = c(4750, 5150), xlim = c(250, 900)) + 
  scale_fill_gradientn(colours = viridis(100)) + ggtitle("f(Eastings, Northings, Year)") +
  xlab("Eastings") + ylab("Northings") + labs(fill = "Partial\nEffect")
#plot the other model smooths
par(mfrow = c(3,3), mar = c(5,5,4,2))
plot(CPUE_model, select = 2, shade = TRUE, unconditional = T, se = 1.96, main = "p-value < 0.001 ***", cex.main = 2,
     shade.col = "light blue", xlab = "Log(Depth)", ylab = "f(Log(Depth))", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-6,1.2))
plot(CPUE_model, select = 3, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-2,1.2), main = "p-value = 0.207", cex.main = 2,
     shade.col = "light blue", xlab = "Bottom Temperature", 
     ylab = "f(Bottom Temperature)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(CPUE_model, select = 4, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-2,1.2), 
     shade.col = "light blue", xlab = "Btm. Temperature Range", main = "p-value = 0.049 ***", cex.main = 2,
     ylab = "f(Btm. Temperature Range)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(CPUE_model, select = 5, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-2,1.2), 
     shade.col = "light blue", xlab = "Bottom Salinity", main = "p-value = 0.489", cex.main = 2,
     ylab = "f(Bottom Salinity)", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(CPUE_model, select = 6, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-2,1.2), 
     shade.col = "light blue", xlab =  "Log(Bottom Stress)", main = "p-value = 0.252", cex.main = 2,
     ylab = "f(Log(Bottom Stress))", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(CPUE_model, select = 7, shade = TRUE, unconditional = T, se = 1.96, ylim = c(-2,1.2), 
     shade.col = "light blue", xlab = "Log(Btm. Stress Range)", main = "p-value = 0.429", cex.main = 2,
     ylab = "f(Log(Btm. Stress Range))", family = "serif", cex.lab = 1.85, cex.axis = 2)
plot(CPUE_model, select = 8, shade = TRUE, unconditional = T, se = 1.96, 
     shade.col = "light blue", xlab = "Square Root of Slope", ylab = "f(Square Root of Slope)", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1.2), main = "p-value = 0.098 ***", cex.main = 2)
plot(CPUE_model, select = 9, shade = TRUE, unconditional = T, se = 1.96, main = "p-value = 0.751", 
     family = "serif", cex.main = 2, cex.axis = 2, cex.lab = 1.85, scheme = 2, xlab = "Easterness", ylab = "Northerness")
plot(CPUE_model, select = 10, shade = TRUE, unconditional = T, se = 1.96, main = "p-value = 0.966", cex.main = 2,
     shade.col = "light blue", xlab = "RDMV", ylab = "f(RDMV)", 
     family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1.2))

#ready information to put into gam summary table for CPUE model
a <- data.frame(format(round(summary(CPUE_model)$p.table, 3), digits = 2))
colnames(a) <- c("Estimate", "Std. Error", "Z","P-value")
rownames(a) <- c("Intercept", "Snow Crab Survey")
b <- data.frame(format(round(summary(CPUE_model)$s.table, 3), digits = 3))
colnames(b) <- c("EDF","Ref. DF","Chi-sq.","P-value")
rownames(b) <- c("f(Eastings, Northings, Year)", "f(Log(Depth))", "f(Bottom Temperature)", "f(Btm. Temperature Range)", 
                 "f(Bottom Salinity)", "f(Log(Bottom Stress))", "f(Log(Btm. Stress Range))", "f(Sqrt(Slope))",
                 "f(Easterness, Northerness)", "f(RDMV)")
#make the table
a <- kable(list(a,b), format = "latex", caption = "CPUE Model", align = "rc")

#compute the residuals
survey_data_combined_CPUE[,"residuals_CPUE"] <- residuals(CPUE_model, type = "response")

#geographically plot the training data's residuals on the CPUE model
a <- ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_CPUE),
                           data = survey_data_combined_CPUE) + 
  scale_colour_gradientn(colours = viridis(100), limits = c(-5,5)) +
  labs(col = "Residual") + theme(text=element_text(size=20, family="serif")) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(75, 1050)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot the training data's residuals on the final CPUE model
b <- ggplot(aes(x = year, y = residuals_CPUE), data = survey_data_combined_CPUE) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Residual") + xlab("Year") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,0.5)) #plot together

#residual vs fits
fits <- fitted(CPUE_model)
a <- ggplot() + geom_point(aes(y = residuals_CPUE, x = fits), data = survey_data_combined_CPUE, alpha =.5) + 
  theme(text=element_text(size=16, family="serif")) + 
  ggtitle("Residuals vs. Fits Plot") + xlab("Fitted Value") + ylab("Residual")
#qqplot
b <- ggplot(data = survey_data_combined_CPUE, aes(sample = residuals_CPUE)) + 
  stat_qq(size = I(0.2)) + 
  stat_qq_line() +
  ggtitle("Normal Q-Q Plot of Residuals") + ylab("Residual Quantile") + xlab("Theoretical Quantile") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,1)) #plot together

############################################PREDICTIONS##################################################

#load in data... change to raster_data_df_2019 file if want 2019 predictions
raster_data_df <- read.csv(paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))

#transform some covariates
raster_data_df$RangeStressLog <- log(raster_data_df$RangeStress)
raster_data_df$BtmStressBNAMLog <- log(raster_data_df$BtmStressBNAM)
raster_data_df$sqrt_DEM_Slope <- sqrt(raster_data_df$DEM_Slope)
raster_data_df$DEM_log <- log(-raster_data_df$DEM)
#scale covariates using all training data
raster_data_df[paste(covariates,"Scaled", sep="")] <- mapply(x = raster_data_df[covariates], y = survey_data_combined[covariates], 
                                                       function(x, y) scale(x, center = mean(y), scale = sd(y)))

#get predictions...
stack_predict_CPUE <- predict(CPUE_model, newdata = raster_data_df, se.fit = TRUE)
stack_predict_zero <- predict(presence_model, newdata = raster_data_df, se.fit = TRUE)
stack_predict_zero$fit <- 1/(1+exp(-stack_predict_zero$fit)) #convert log odds to probability
#combine to get final abundance predictions
stack_predict_combined <- stack_predict_zero
stack_predict_combined$fit <- log(stack_predict_zero$fit*exp(stack_predict_CPUE$fit+0.5*CPUE_model$scale))
#get combined standard errors
stack_predict_combined$se.fit <- sqrt(stack_predict_zero$se.fit^2*(1-stack_predict_zero$fit)^2+
                                        stack_predict_CPUE$se.fit^2)

#plot the presence predictions
ggplot() + geom_raster(data = raster_data_df, aes(x = UTMX, y = UTMY, fill = stack_predict_zero$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Probability of C. frondosa Presence") + facet_wrap(~year, ncol = 4)

#plot the conditional CPUE predictions
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_CPUE$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Log(CPUE [kg/square kilometre])\n[Conditional on C. frondosa Presence]") + facet_wrap(~year, ncol = 4)

#plot the unconditional CPUE predictions
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Log(Predicted CPUE [kg/square kilometre])") + facet_wrap(~year, ncol = 4)

#get uncertainty in final predictions
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$se.fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Log(Predicted CPUE [kg/square kilometre]) - Standard Errors") + facet_wrap(~year, ncol = 4)

#get 2001,2010,2019 predictions
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year%in%c(2001,2010,2019)),], 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$fit[which(raster_data_df$year%in%c(2001,2010,2019))])) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("GAM") + facet_wrap(~year, ncol = 4) + theme(legend.position = "none")

#get 2001, 2010, 2019 prediction uncertainty
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year%in%c(2001,2010,2019)),], 
aes(x = UTMX, y = UTMY, fill = stack_predict_combined$se.fit[which(raster_data_df$year%in%c(2001,2010,2019))])) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("GAM") + facet_wrap(~year, ncol = 4) + theme(legend.position = "none")

####################2019 HIGH RES PREDICTIONS#################

#load shape file for reserves and get the units into km
fishing_area_reserves <- shapefile("Data/All Fishing Areas with Reserves/Fishing_Areas_2018.shp")
fishing_area_reserves <- spTransform(fishing_area_reserves, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#subset it to our spatial domain
fishing_area_reserves <- fishing_area_reserves[-which(fishing_area_reserves$Region == "SWNB"|fishing_area_reserves$Region=="4X"),]
#split up fishing zones and reserves
fishing_area <- fishing_area_reserves[which(fishing_area_reserves$Type=="Fishing Area"),]
reserve <- fishing_area_reserves[which(fishing_area_reserves$Type=="Reserve"),]

#load the WEBCA boundary and fix formatting
WEBCA <- shapefile("Data/WEBCA/DFO_OECM_MPO_AMCEZ.shp")
WEBCA <- WEBCA[which(WEBCA$NAME_E == "Western/Emerald Banks Conservation Area (restricted fisheries zone)"),]
WEBCA <- spTransform(WEBCA, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

#plot the presence predictions
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_zero$fit[which(raster_data_df$year==2019)]), 
                       interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  coord_equal() +
  ggtitle("2019 - Predicted Probability of Presence") + ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))

#plot the conditional CPUE predictions
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_CPUE$fit[which(raster_data_df$year==2019)]), 
                       interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  coord_equal() +
  ggtitle("2019 - Predicted Log(CPUE [kg/square kilometre])\n[Conditional on Presence]") + 
  ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (0.1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))

#plot the final combined predictions
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$fit[which(raster_data_df$year==2019)]), 
                       interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  coord_equal() + ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])") + ylab("Northings") + 
  xlab("Eastings") + geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))

#final combined predictions... zoomed in with management regions
ggplot() + geom_raster(data = as.data.frame(raster_data_df[which(raster_data_df$year==2019),]), 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$fit[which(raster_data_df$year==2019)]), 
                       interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  coord_equal() + ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])\n- GAM") + ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c(1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4800, 5035), xlim = c(607.5, 950)) +
  geom_polygon(data = rbind(fortify(WEBCA)), aes(group = group, x = long, 
                                                 y = lat), col = "brown", fill = NA, size = 1.35) +
  geom_polygon(data = rbind(fortify(reserve)), aes(group = group, x = long, 
                                                   y = lat), col = "white", fill = NA, size = 1.35) +
  geom_polygon(data = fortify(fishing_area), aes(group = group, x = long, 
                                                 y = lat), col = "black", fill = NA, size = 1.35)

#uncertainty in final predictions for 2019
ggplot() + geom_raster(data = raster_data_df[which(raster_data_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$se.fit[which(raster_data_df$year==2019)]), 
                       interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  coord_equal() + ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])\n- Standard Errors") + ylab("Northings") + 
  xlab("Eastings") + geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (0.1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))