#load necessary libraries
library(sdmTMB)
library(INLA) #Use SPDE functionality
library(raster)
library(knitr)
library(sp)
library(ggplot2)
library(viridis)
library(ggpubr)
source("R Files/Helper Functions.R")

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

#reorder the data by year
survey_data_combined <- survey_data_combined[order(survey_data_combined$year, decreasing = FALSE),]

#survey_data_combined <- survey_data_combined[-which(survey_data_combined$season_winter == 1),] #remove RV winter survey tows

#########################################PRESENCE MODEL#########################################

#set up response
survey_data_combined$presence <- as.numeric(survey_data_combined$presence)

#set up the mesh for the model fitted to all data
non_convex_bdry <- inla.nonconvex.hull(cbind(survey_data_combined$UTMX, survey_data_combined$UTMY), -0.03, resolution = c(200, 50))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(25, 200), 
                      offset = c(100), cutoff = 10, loc = cbind(survey_data_combined$UTMX, survey_data_combined$UTMY))
#plot the mesh
plot(mesh, family = "serif", cex.main = 2, main = "")
points(cbind(survey_data_combined$UTMX, survey_data_combined$UTMY), col = "orange", cex = 0.4)

#fit the presence model in sdmTMB using the mesh
mesh_model <- make_mesh(survey_data_combined, xy_cols = c("UTMX", "UTMY"), cutoff = 10, mesh = mesh)
start_time <- Sys.time()
presence_model <- sdmTMB(presence ~ DEM_logScaled+I(DEM_logScaled^2)+BtmTempBNAMScaled+I(BtmTempBNAMScaled^2)+RangeTempScaled+
                           I(RangeTempScaled^2)+BtmSalinityBNAMScaled+I(BtmSalinityBNAMScaled^2)+BtmStressBNAMLogScaled+RangeStressLogScaled+
                           sqrt_DEM_SlopeScaled+DEM_NorthernessScaled+DEM_EasternessScaled+DEM_RDMVScaled+
                           as.numeric(snowcrab),
  data = survey_data_combined,
  mesh = mesh_model,
  family = binomial(),
  spatiotemporal = "ar1", spatial = "off", time = "year", control = sdmTMBcontrol(newton_loops = 1))
end_time <- Sys.time()
fit_time1 <- end_time - start_time #get time to fit
#sanity check
sanity(presence_model)

#create table to describe fixed effects, etc.
fixed_effects_table <- round(tidy(presence_model)[2:3],3)
rownames(fixed_effects_table) <- c("Intercept","Log(Depth)","Log(Depth) Squared","Bottom Temperature", "Bottom Temperature Squared",
                                   "Btm. Temperature Range","Btm. Temperature Range Squared","Bottom Salinity", "Bottom Salinity Squared", 
                                   "Log(Bottom Stress)", "Log(Btm. Stress Range)", "Sqrt(Slope)","Northerness","Easterness","RDMV",
                                   "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatiotemporal_parameters_table <- round(tidy(presence_model, effects = "ran_pars", conf.int = FALSE)[2],3)
se <- as.list(presence_model$sd_report, "Std. Error", report = TRUE) #grab SEs in natural space
spatiotemporal_parameters_table <- cbind(spatiotemporal_parameters_table, round(c(se$range[1], se$sigma_E[1], se$rho),3))
colnames(spatiotemporal_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatiotemporal_parameters_table) <- c("Range", "Marginal Spatial SD","AR(1) Correlation")
kable(list(fixed_effects_table, spatiotemporal_parameters_table), format = "latex")

#plot the model terms
par(mfrow = c(3,4), mar = c(5,5,4,2))
effects_plot(presence_model, 2, 3, original_vector = survey_data_combined$DEM_logScaled,
             xlab = "Log(Depth)", ylab = "Prediction Contribution", main = "***", cex.main = 2, 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-6,1.2))
effects_plot(presence_model, 4, 5, original_vector = survey_data_combined$BtmTempBNAMScaled,
             xlab = "Bottom Temperature", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 6,7, original_vector = survey_data_combined$RangeTempScaled,
             xlab = "Btm. Temperature Range", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 8,9, original_vector = survey_data_combined$BtmSalinityBNAMScaled,
             xlab = "Bottom Salinity", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 10, original_vector = survey_data_combined$BtmStressBNAMLogScaled,
             xlab = "Log(Bottom Stress)", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 11, original_vector = survey_data_combined$RangeStressLogScaled,
             xlab = "Log(Btm. Stress Range)", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 12, original_vector = survey_data_combined$sqrt_DEM_SlopeScaled,
             xlab = "Square Root of Slope", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 13, original_vector = survey_data_combined$DEM_NorthernessScaled,
             xlab = "Northerness", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 14, original_vector = survey_data_combined$DEM_EasternessScaled,
             xlab = "Easterness", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))
effects_plot(presence_model, 15, original_vector = survey_data_combined$DEM_RDMVScaled,
             xlab = "RDMV", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-3,1))

#predict for training data, get out of log link, get residuals
predict_train <- predict(presence_model)
predict_train$response <- exp(predict_train$est)/(1+exp(predict_train$est))
predict_train$residuals <- as.numeric(predict_train$presence)-predict_train$response

#spatio-temporally plot residuals
a <- ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = predict_train$residuals),
                           data = survey_data_combined) + 
  scale_colour_gradientn(colours = viridis(100), limits = c(-1,1)) +
  labs(col = "Conditional\nResponse\nResidual\n") + theme(text=element_text(size=16, family="serif")) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot residuals
b <- ggplot(aes(x = year, y = predict_train$residuals), data = survey_data_combined) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Conditional Response Residual") + xlab("Year") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,0.5)) #plot together

#########################################CPUE MODEL###########################################

#grab presence tows and set up the response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),]
survey_data_combined_CPUE$std.WGT_log <- log(survey_data_combined_CPUE$std.WGT)

#set up the mesh for the model fitted to all data
non_convex_bdry <- inla.nonconvex.hull(cbind(survey_data_combined_CPUE$UTMX, survey_data_combined_CPUE$UTMY), -0.05, resolution = c(200, 50))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(25, 200), 
                     offset = c(100), cutoff = 10, loc = cbind(survey_data_combined_CPUE$UTMX, survey_data_combined_CPUE$UTMY))
#plot the mesh
plot(mesh, family = "serif", cex.main = 2, main = "")
points(cbind(survey_data_combined_CPUE$UTMX, survey_data_combined_CPUE$UTMY), col = "orange", cex = 0.4)

#fit the CPUE model in sdmTMB using the mesh
mesh_model <- make_mesh(survey_data_combined_CPUE, xy_cols = c("UTMX", "UTMY"), cutoff = 10, mesh = mesh)
start_time <- Sys.time()
CPUE_model <- sdmTMB(std.WGT_log ~ DEM_logScaled+I(DEM_logScaled^2)+BtmTempBNAMScaled+RangeTempScaled+BtmSalinityBNAMScaled+
                         BtmStressBNAMLogScaled+RangeStressLogScaled+ sqrt_DEM_SlopeScaled+DEM_NorthernessScaled+DEM_EasternessScaled+
                         DEM_RDMVScaled+as.numeric(snowcrab),
                     data = survey_data_combined_CPUE,
                     mesh = mesh_model,
                     family = gaussian(),
                     spatiotemporal = "ar1", spatial = "off", time = "year", control = sdmTMBcontrol(newton_loops = 1))
end_time <- Sys.time()
fit_time3 <- end_time - start_time #get time to fit
#sanity check
sanity(CPUE_model)

#create table to describe fixed effects, etc.
fixed_effects_table <- round(tidy(CPUE_model)[2:3],3)
rownames(fixed_effects_table) <- c("Intercept","Log(Depth)","Log(Depth) Squared","Bottom Temperature",
                                   "Btm. Temperature Range","Bottom Salinity", 
                                   "Log(Bottom Stress)", "Log(Btm. Stress Range)", "Sqrt(Slope)","Northerness","Easterness","RDMV",
                                   "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatiotemporal_parameters_table <- round(tidy(CPUE_model, effects = "ran_pars", conf.int = FALSE)[-2,2],3)
se <- as.list(CPUE_model$sd_report, "Std. Error", report = TRUE) #grab SEs in natural space
spatiotemporal_parameters_table <- cbind(spatiotemporal_parameters_table, round(c(se$range[1], se$sigma_E[1], se$rho),3))
colnames(spatiotemporal_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatiotemporal_parameters_table) <- c("Range", "Marginal Spatial SD", "AR(1) Correlation")
kable(list(fixed_effects_table, spatiotemporal_parameters_table), format = "latex")

#plot the model terms
par(mfrow = c(3,4), mar = c(5,5,4,2))
effects_plot(CPUE_model, 2, 3, original_vector = survey_data_combined_CPUE$DEM_logScaled,
             xlab = "Log(Depth)", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-6,1.2))
effects_plot(CPUE_model, 4, original_vector = survey_data_combined_CPUE$BtmTempBNAMScaled,
             xlab = "Bottom Temperature", ylab = "Prediction Contribution", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 5, original_vector = survey_data_combined_CPUE$RangeTempScaled,
             xlab = "Btm. Temperature Range", ylab = "Prediction Contribution", main = "***", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 6, original_vector = survey_data_combined_CPUE$BtmSalinityBNAMScaled,
             xlab = "Bottom Salinity", ylab = "Prediction Contribution", main = "***", cex.main = 2,
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 7, original_vector = survey_data_combined_CPUE$BtmStressBNAMLogScaled,
             xlab = "Log(Bottom Stress)", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 8, original_vector = survey_data_combined_CPUE$RangeStressLogScaled,
             xlab = "Log(Btm. Stress Range)", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 9, original_vector = survey_data_combined_CPUE$sqrt_DEM_SlopeScaled,
             xlab = "Square Root of Slope", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 10, original_vector = survey_data_combined_CPUE$DEM_NorthernessScaled,
             xlab = "Northerness", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 11, original_vector = survey_data_combined_CPUE$DEM_EasternessScaled,
             xlab = "Easterness", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))
effects_plot(CPUE_model, 12, original_vector = survey_data_combined_CPUE$DEM_RDMVScaled,
             xlab = "RDMV", ylab = "Prediction Contribution", 
             family = "serif", cex.axis = 2, cex.lab = 1.85, ylim = c(-2,1))

#predict for training data, get residuals
predict_train <- predict(CPUE_model)
predict_train$residuals <- predict_train$std.WGT_log-predict_train$est

#spatio-temporally plot residuals
a <- ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = predict_train$residuals),
                           data = survey_data_combined_CPUE) + 
  scale_colour_gradientn(colours = viridis(100), limits = c(-5,5)) +
  labs(col = "Conditional\nResidual\n") + theme(text=element_text(size=16, family="serif")) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot residuals
b <- ggplot(aes(x = year, y = predict_train$residuals), data = survey_data_combined_CPUE) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Conditional Residual") + xlab("Year") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,0.5)) #plot together

#residual vs fits
a <- ggplot() + geom_point(aes(y = predict_train$residuals, x = predict_train$est), 
                           data = survey_data_combined_CPUE, alpha =.5) + 
  theme(text=element_text(size=16, family="serif")) + 
  ggtitle("Conditional Residuals vs. Fits Plot") + xlab("Fitted Value") + ylab("Conditional Residual")
#qqplot
b <- ggplot(data = survey_data_combined_CPUE, aes(sample = predict_train$residuals)) + 
  stat_qq(size = I(0.2)) + 
  stat_qq_line() +
  ggtitle("Normal Q-Q Plot of Conditional Residuals") + ylab("Residual Quantile") + xlab("Theoretical Quantile") + 
  theme(text=element_text(size=16, family="serif"))
ggarrange(a, NULL, b, ncol = 1, heights = c(1,0.05,1)) #plot together

###########################################PREDICTIONS##########################################

#increase memory limit
memory.limit(size = 400000)

#load in data... change to raster_data_df_2019 file if want 2019 predictions
raster_data_df <- read.csv(paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))

#transform some covariates
raster_data_df$RangeStressLog <- log(raster_data_df$RangeStress)
raster_data_df$BtmStressBNAMLog <- log(raster_data_df$BtmStressBNAM)
raster_data_df$sqrt_DEM_Slope <- sqrt(raster_data_df$DEM_Slope)
raster_data_df$DEM_log <- log(-raster_data_df$DEM)
#scale covariates using training data
raster_data_df[paste(covariates,"Scaled", sep="")] <- mapply(x = raster_data_df[covariates], y = survey_data_combined[covariates], 
                                                       function(x, y) scale(x, center = mean(y), scale = sd(y)))

#make predictions
presence_prediction_df <- predict(presence_model, newdata = raster_data_df, type = "response")
CPUE_prediction_df <- predict(CPUE_model, newdata = raster_data_df)
#get standard errors
presence_se <- predict(presence_model, newdata = raster_data_df, nsim = 500)
CPUE_se <- predict(CPUE_model, newdata = raster_data_df, nsim = 500)
if(nrow(presence_se)>nrow(raster_data_df))#remove possible extra rows... sdmTMB quirk
{
  presence_se <- presence_se[-(nrow(raster_data_df)+1:nrow(presence_se)),]
  CPUE_se <- CPUE_se[-(nrow(raster_data_df)+1:nrow(CPUE_se)),]
}
presence_prediction_df$se <- apply(presence_se, 1, sd)
CPUE_prediction_df$se <- apply(CPUE_se, 1, sd)
#create final combined predictions, standard errors, and add to presence_prediction_df
presence_prediction_df$combined_predict <- log(exp(CPUE_prediction_df$est+0.5*CPUE_model$sd_report$value[2]^2)*presence_prediction_df$est)
presence_prediction_df$combined_se <- sqrt((presence_prediction_df$se)^2*(1-presence_prediction_df$est)^2+(CPUE_prediction_df$se)^2)

#plot random field for presence
ggplot() + geom_raster(data = presence_prediction_df, aes(x = UTMX, y = UTMY, fill = est_rf)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Presence Model Random Effects") + facet_wrap(~year, ncol = 4)

#plot random field for conditional CPUE
ggplot() + geom_raster(data = CPUE_prediction_df, aes(x = UTMX, y = UTMY, fill = est_rf)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("CPUE Model Random Effects") + facet_wrap(~year, ncol = 4)

#plot the presence predictions
ggplot() + geom_raster(data = presence_prediction_df, aes(x = UTMX, y = UTMY, fill = est)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Probability of C. frondosa Presence") + facet_wrap(~year, ncol = 4)

#plot the conditional CPUE predictions
ggplot() + geom_raster(data = CPUE_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = est)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Log(CPUE [kg/square kilometre])\n[Conditional on C. frondosa Presence]") + facet_wrap(~year, ncol = 4)

#plot the unconditional CPUE predictions
ggplot() + geom_raster(data = presence_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = combined_predict)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12))  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Log(Predicted CPUE [kg/square kilometre])") + facet_wrap(~year, ncol = 4)

#get uncertainty in final predictions
ggplot() + geom_raster(data = presence_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = combined_se)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6))  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Log(Predicted CPUE [kg/square kilometre]) - Standard Errors") + facet_wrap(~year, ncol = 4)

#get 2001, 2010, 2019 predictions
ggplot() + geom_raster(data = presence_prediction_df[which(raster_data_df$year%in%c(2001,2010,2019)),], 
                       aes(x = UTMX, y = UTMY, fill = combined_predict)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("GMRF Model") + facet_wrap(~year, ncol = 4) + theme(legend.position = "none")

#get 2001, 2010, 2019 prediction uncertainty
ggplot() + geom_raster(data = presence_prediction_df[which(raster_data_df$year%in%c(2001,2010,2019)),], 
                       aes(x = UTMX, y = UTMY, fill = combined_se)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16, family="serif")) + coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("GMRF Model") + facet_wrap(~year, ncol = 4) + theme(legend.position = "none")

####################2019 HIGH-RES PREDICTIONS#################

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
ggplot() + geom_raster(data = presence_prediction_df[which(presence_prediction_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = est), interpolate = TRUE) + 
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
ggplot() + geom_raster(data = CPUE_prediction_df[which(CPUE_prediction_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = est), interpolate = TRUE) + 
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
ggplot() + geom_raster(data = presence_prediction_df[which(presence_prediction_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = combined_predict), interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  coord_equal() +
  ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])") + ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))

#final combined predictions... zoomed in with management regions
ggplot() + geom_raster(data = presence_prediction_df[which(presence_prediction_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = combined_predict), interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-15,12)) + 
  coord_equal() + ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])\n- GMRF Model") + ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c(1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4800, 5035), xlim = c(607.5, 950)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  geom_polygon(data = rbind(fortify(WEBCA)), aes(group = group, x = long, 
                                                                   y = lat), col = "brown", fill = NA, size = 1.35) +
  geom_polygon(data = rbind(fortify(reserve)), aes(group = group, x = long, 
                                                   y = lat), col = "white", fill = NA, size = 1.35) +
  geom_polygon(data = fortify(fishing_area), aes(group = group, x = long, 
                                                 y = lat), col = "black", fill = NA, size = 1.35)

#plot the uncertainty for 2019 predictions
ggplot() + geom_raster(data = presence_prediction_df[which(presence_prediction_df$year==2019),], 
                       aes(x = UTMX, y = UTMY, fill = combined_se), interpolate = TRUE) + 
  scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  coord_equal() +
  ggtitle("2019 - Log(Predicted CPUE [kg/square kilometre])\n- Standard Errors") + ylab("Northings") + xlab("Eastings") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat))+
  theme(text=element_text(size=20, family="serif"), legend.position = "bottom", 
        legend.key.size = unit(1.25, 'cm'), legend.text = element_text(size=18), 
        plot.margin=unit (c (0.1,2,2,1.5), 'cm')) + 
  coord_cartesian(ylim = c(4600, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))