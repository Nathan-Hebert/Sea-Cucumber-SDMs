#############################PRESENCE SUB-MODEL PLOTS FOR SPATIAL AND NON-SPATIAL GAM############################

#need to compare smooths
library(gratia)
library(mgcv)

#load in spatio-temporal model
load("~/R/Sea Cucumber Project/WorkspaceFiles/GAMlognormal.RData")

#scale covariates using all data
survey_data_combined[paste(covariates,"Scaled", sep="")] <- scale(survey_data_combined[covariates])
#fit the model without spatio-temporal component
start_time <- Sys.time()
presence_model_nonST <- gam(presence~s(DEM_logScaled, k = 3)+s(BtmTempBNAMScaled, k = 3) + s(RangeTempScaled, k = 3) + s(BtmSalinityBNAMScaled, k = 3) + 
                        s(BtmStressBNAMLogScaled, k = 3) + s(RangeStressLogScaled, k = 3) + s(sqrt_DEM_SlopeScaled, k = 3) + 
                        s(DEM_EasternessScaled, DEM_NorthernessScaled, k = 9) + s(DEM_RDMVScaled, k = 3) + as.numeric(snowcrab), 
                      family = "binomial", 
                      data = survey_data_combined, 
                      method = "REML") 
end_time <- Sys.time()
fit_time <- end_time - start_time #get time to fit

#fit image margins
par(mar = c(5,5,4,2))

#plot depth smooth for spatio-temporal presence sub-model
plot(presence_model, select = 2, shade = TRUE, rug = T, ylim = c(-9,5), 
     shade.col = "light blue", xlab = "Log(Depth)", cex.main = 2,
     ylab = "f(Log(Depth))", cex.lab = 2, cex.axis = 2, lwd = 2, col = "blue", se = 1.96)
rug(survey_data_combined$DEM_logScaled, lwd = 2)
#store information from other smooth
non_spatial_smooth <- compare_smooths(presence_model, presence_model_nonST, smooths = "s(DEM_logScaled)",
                                      overall_uncertainty = F)
x <- non_spatial_smooth["data"][[1]][2][[1]]$DEM_logScaled
y <- non_spatial_smooth["data"][[1]][2][[1]]$est
se <- non_spatial_smooth["data"][[1]][2][[1]]$se
#add non-ST model's depth smooth onto same panel
lines(x, y, lwd = 2, col = "goldenrod4")
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("gold", alpha.f=0.3), lty = 0)

#plot btm stress smooth for spatio-temporal presence sub-model
plot(presence_model, select = 6, shade = TRUE, rug = T, ylim = c(-1.6,1.5), 
     shade.col = "light blue", xlab = "Log(Bottom Stress)", cex.main = 2,
     ylab = "f(Log(Bottom Stress))", cex.lab = 2, cex.axis = 2, lwd = 2, col = "blue", se = 1.96)
rug(survey_data_combined$BtmStressBNAMLogScaled, lwd = 2)
#store information from other smooth
non_spatial_smooth <- compare_smooths(presence_model, presence_model_nonST, smooths = "s(BtmStressBNAMLogScaled)",
                                      overall_uncertainty = F)
x <- non_spatial_smooth["data"][[1]][2][[1]]$BtmStressBNAMLogScaled
y <- non_spatial_smooth["data"][[1]][2][[1]]$est
se <- non_spatial_smooth["data"][[1]][2][[1]]$se
#add non-ST model's depth smooth onto same panel
lines(x, y, lwd = 2, col = "goldenrod4")
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("gold", alpha.f=0.3), lty = 0)

#plot btm temp range smooth for spatio-temporal presence sub-model
plot(presence_model, select = 4, shade = TRUE, rug = T, ylim = c(-4,3), 
     shade.col = "light blue", xlab = "Bottom Temperature Range", cex.main = 2,
     ylab = "f(Bottom Temperature Range)", cex.lab = 2, cex.axis = 2, lwd = 2, col = "blue", se = 1.96)
rug(survey_data_combined$RangeTempScaled, lwd = 2)
#store information from other smooth
non_spatial_smooth <- compare_smooths(presence_model, presence_model_nonST, smooths = "s(RangeTempScaled)",
                                      overall_uncertainty = F)
x <- non_spatial_smooth["data"][[1]][2][[1]]$RangeTempScaled
y <- non_spatial_smooth["data"][[1]][2][[1]]$est
se <- non_spatial_smooth["data"][[1]][2][[1]]$se
#add non-ST model's depth smooth onto same panel
lines(x, y, lwd = 2, col = "goldenrod4")
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("gold", alpha.f=0.3), lty = 0)

#####################################PRESENCE SUB-MODEL PLOTS FOR RANDOM EFFECT MODELS############################

#DEPTH
load("~/R/Sea Cucumber Project/WorkspaceFiles/GMRFlognormal.RData")
source("R Files/Helper Functions.R")
#plot depth from GMRF
par(mar = c(5,5,4,2), mfrow = c(1,1))
effects_plot(presence_model, 2, 3, original_vector = survey_data_combined$DEM_logScaled,
             xlab = "Log(Depth)", ylab = "Prediction Contribution", cex.main = 2, 
             cex.axis = 2, cex.lab = 2, ylim = c(-9,5), line_color = "red", ci.color = adjustcolor("#F88279", alpha.f=0.3), lwd = 2)
#add starve lines
load("~/R/Sea Cucumber Project/WorkspaceFiles/starvelognormal.RData")
x <- seq(min(survey_data_combined$DEM_logScaled), max(survey_data_combined$DEM_logScaled), length.out = 1000)
y <- presence_model@TMB_out@sdr$par.fixed[6]*x+presence_model@TMB_out@sdr$par.fixed[7]*x^2
se <- sqrt(x^2*presence_model@tracing@parameter_covariance[6,6]+x^4*presence_model@tracing@parameter_covariance[7,7]+2*x^3*presence_model@tracing@parameter_covariance[6,7]) 
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("#000000", alpha.f=0.3), lty = 0)
lines(x,y, lwd = 2)

#BTM STRESS
load("~/R/Sea Cucumber Project/WorkspaceFiles/GMRFlognormal.RData")
source("R Files/Helper Functions.R")
#plot stress line from GMRF
par(mar = c(5,5,4,2), mfrow = c(1,1))
effects_plot(presence_model, 10, original_vector = survey_data_combined$BtmStressBNAMLogScaled,
             xlab = "Log(Bottom Stress)", ylab = "Prediction Contribution", cex.main = 2, 
             cex.axis = 2, cex.lab = 2, ylim = c(-1.6,1.5), line_color = "red", ci.color = adjustcolor("#F88279", alpha.f=0.3), lwd = 2)
#add starve lines
load("~/R/Sea Cucumber Project/WorkspaceFiles/starvelognormal.RData")
x <- seq(min(survey_data_combined$BtmStressBNAMLogScaled), max(survey_data_combined$BtmStressBNAMLogScaled), length.out = 1000)
y <- presence_model@TMB_out@sdr$par.fixed[14]*x
se <- x*sqrt(presence_model@tracing@parameter_covariance[14,14]) 
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("#000000", alpha.f=0.3), lty = 0)
lines(x,y, lwd = 2)

#BTM TEMP RANGE
load("~/R/Sea Cucumber Project/WorkspaceFiles/GMRFlognormal.RData")
source("R Files/Helper Functions.R")
#plot stress line from GMRF
par(mar = c(5,5,4,2), mfrow = c(1,1))
effects_plot(presence_model, 6,7, original_vector = survey_data_combined$RangeTempScaled,
             xlab = "Bottom Temperature Range", ylab = "Prediction Contribution", cex.main = 2, 
             cex.axis = 2, cex.lab = 2, ylim = c(-4,3), line_color = "red", ci.color = adjustcolor("#F88279", alpha.f=0.3), lwd = 2)
#add starve lines
load("~/R/Sea Cucumber Project/WorkspaceFiles/starvelognormal.RData")
x <- seq(min(survey_data_combined$RangeTempScaled), max(survey_data_combined$RangeTempScaled), length.out = 1000)
y <- presence_model@TMB_out@sdr$par.fixed[10]*x+presence_model@TMB_out@sdr$par.fixed[11]*x^2
se <- sqrt(x^2*presence_model@tracing@parameter_covariance[10,10]+x^4*presence_model@tracing@parameter_covariance[11,11]+2*x^3*presence_model@tracing@parameter_covariance[10,11]) 
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("#000000", alpha.f=0.3), lty = 0)
lines(x,y, lwd = 2)