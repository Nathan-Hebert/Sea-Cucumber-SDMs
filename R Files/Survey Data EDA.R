#load necessary libraries
library(ggplot2)
library(fields)
library(dplyr)
library(sp)
library(ggpubr)
library(gtools)
library(raster)
library(corrplot)
library(car)
library(knitr)

#load in data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))
#if RDMV is NA, set to zero
survey_data_combined$DEM_RDMV <- ifelse(is.na(survey_data_combined$DEM_RDMV), 0, survey_data_combined$DEM_RDMV)

#get new metric
survey_data_combined$weight_per_unit <- survey_data_combined$std.WGT/survey_data_combined$std.No

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")

##########################################EXAMINE TOWS############################################

#examine the number of tows per Survey Year (separated by month)
tows_per_year <- as.data.frame(table(survey_data_combined$year, survey_data_combined$month))
ggplot(data = tows_per_year, aes(x = Var1, y = Freq, fill = as.factor(Var2)))+geom_bar(stat = "identity") +
  labs(fill = "Month") + 
  scale_x_discrete(name = "Year", labels = factor(ggplot2:::interleave(seq(2000,2020,by=2), ""))) + 
  theme_classic() +  scale_y_continuous(breaks = c(0,50,100,150,200,250,300,350,400,450)) +
  ylab("Number of tows") + scale_fill_manual(values = tol(12), 
                                             labels=c("January", "February", "March", "April", "May","June","July","August","September","October","November","December"))

#examine distribution of surveys per year
snowcrab_per_year <- as.data.frame(prop.table(table(survey_data_combined$year, survey_data_combined$snowcrab),1))
snowcrab_per_year <- snowcrab_per_year[which(snowcrab_per_year$Var2==TRUE),]
ggplot(data = snowcrab_per_year, aes(x = Var1, y = Freq))+geom_bar(stat = "identity") +
  scale_x_discrete(name = "Year", labels = factor(ggplot2:::interleave(seq(2000,2019,by=2), ""))) + 
  ylab("Proportion of tows that were from the snow crab survey")
#create table
kable(table(survey_data_combined$year, survey_data_combined$season_winter, useNA = "always")[1:20,], format = "latex")

#examine each year's proportion of tows with sea cucumbers present
presence_per_year <- as.data.frame(prop.table(table(survey_data_combined$year, survey_data_combined$presence),1))
presence_per_year <- presence_per_year[which(presence_per_year$Var2==TRUE),]
ggplot(data = presence_per_year, aes(x = Var1, y = Freq))+geom_bar(stat = "identity") +
  scale_x_discrete(name = "Year", labels = factor(ggplot2:::interleave(seq(2000,2020,by=2), ""))) + 
  ylab("Proportion of tows with sea cucumbers present")

#compare weights from the two surveys with histograms
ggplot()+geom_histogram(data=survey_data_combined[which(survey_data_combined$snowcrab == TRUE&is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
                        mapping=aes(x=std.WGT, y = stat(density)), fill = "orange", col = "black", bins = 10000) +
  coord_cartesian(xlim = c(0,1000)) + xlab("CPUE (kg/squared kilometre)") + ylab("Density") +
  theme(text=element_text(size=16, family="serif")) +
  scale_x_continuous(breaks = seq(0, 2000, 100)) +
  geom_histogram(data=survey_data_combined[which(survey_data_combined$snowcrab == FALSE&is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
                 mapping=aes(x=std.WGT, y = stat(density), fill = "blue"), fill = "blue", col = "black", bins = 10000, alpha = 0.3) + 
  annotate("text", x = 130, y = 0.01, label = "RV Survey", col = "blue", family = "serif", size = 6) + 
  annotate("text", x = 280, y = 0.0025, label = "Snow Crab Survey", col = "brown1", family = "serif", size = 6)
#compare log of weights from two surveys
ggplot()+geom_histogram(data=survey_data_combined[which(survey_data_combined$snowcrab == TRUE&is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
                        mapping=aes(x=log(std.WGT), y = stat(density)), fill = "orange", col = "black", bins = 50) + 
  xlab("log(CPUE [kg/squared kilometre])") + ylab("Density") +
  theme(text=element_text(size=16, family="serif")) +
  geom_histogram(data=survey_data_combined[which(survey_data_combined$snowcrab == FALSE&is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
                 mapping=aes(x=log(std.WGT), y = stat(density), fill = "blue"), fill = "blue", col = "black", bins = 50, alpha = 0.3) + 
  annotate("text", x = -1, y = 0.1, label = "RV Survey", col = "blue", family = "serif", size = 6) + 
  annotate("text", x = 10, y = 0.2, label = "Snow Crab Survey", col = "brown1", family = "serif", size = 6)

#load in shape file with NAFO divisions
shp1 <- shapefile("Data/NAFO/NAFO_Divisions_SHP/NAFO_Divisions_2021_poly_not_clipped.shp")
#grab coordinates from tows and make sure they are comparable to the shape file
coordinates <- SpatialPoints(cbind(survey_data_combined$mid.lon, survey_data_combined$mid.lat))
proj4string(coordinates) <- "+proj=longlat +datum=NAD83 +no_defs"
#use the coordinates and shape file to create a NAFO division column in the RV data frame
survey_data_combined$NAFO <- over(coordinates, as(shp1[,"Division"],"SpatialPolygons"))
survey_data_combined$NAFO <- ifelse(survey_data_combined$NAFO == 27|survey_data_combined$NAFO == 29, 28, survey_data_combined$NAFO)
survey_data_combined$NAFO <- as.factor(survey_data_combined$NAFO)
levels(survey_data_combined$NAFO) <- c("4Vn", "4Vs", "4W","4X","5Y","5Ze")

#compare weights from the two surveys faceted by NAFO division and year
ggplot(data=survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$NAFO == "4Vn"),], 
       mapping=aes(x = as.factor(snowcrab), y=std.WGT))+geom_boxplot(fill = "grey") +
  coord_cartesian(ylim = c(0,700)) + ylab("CPUE (kg/squared kilometre)") + xlab("") + 
  scale_x_discrete(labels = c("RV", "Snow Crab")) + facet_wrap(~year) + ggtitle("NAFO Division 4Vn") +
  theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))
ggplot(data=survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$NAFO == "4Vs"),], 
       mapping=aes(x = as.factor(snowcrab), y=std.WGT))+geom_boxplot(fill = "grey") +
  coord_cartesian(ylim = c(0,3500)) + ylab("CPUE (kg/squared kilometre)") + xlab("") + 
  scale_x_discrete(labels = c("RV", "Snow Crab")) + facet_wrap(~year) + ggtitle("NAFO Division 4Vs")
theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))  
ggplot(data=survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$NAFO == "4W"),], 
       mapping=aes(x = as.factor(snowcrab), y=std.WGT))+geom_boxplot(fill = "grey") +
  coord_cartesian(ylim = c(0,1250)) + ylab("CPUE (kg/squared kilometre)") + xlab("") + 
  scale_x_discrete(labels = c("RV", "Snow Crab")) + facet_wrap(~year) + ggtitle("NAFO Division 4W")
theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90)) 
ggplot(data=survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$NAFO == "4X"),], 
       mapping=aes(x = as.factor(snowcrab), y=std.WGT))+geom_boxplot(fill = "grey") +
  coord_cartesian(ylim = c(0,500)) + ylab("CPUE (kg/squared kilometre)") + xlab("") + 
  scale_x_discrete(labels = c("RV", "Snow Crab")) + facet_wrap(~year) + ggtitle("NAFO Division 4X")
theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))  

world <- map_data('world') #get world data
n_am <- world[world$region %in% c('Canada','USA'), ] #subset North America
#plot where tows conducted
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-70,-56.5), ylim = c(40,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat, col = snowcrab), 
             data = survey_data_combined, alpha =.5) + facet_wrap(year ~ ., ncol = 5)+
  theme(text=element_text(size=11, family="serif"))
#plot RV tow locations
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-70,-56.5), ylim = c(40,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat, col = as.factor(season_winter)), 
             data = survey_data_combined[which(survey_data_combined$snowcrab==FALSE),], alpha =.33) + facet_wrap(year ~ ., ncol = 4)+
  theme(text=element_text(size=14, family="serif")) + scale_color_discrete(type=c("cornflowerblue","orange"), 
                                                                           name = "Survey Season", labels = c("Summer", "Winter"))
#plot snow crab tow locations
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-70,-56.5), ylim = c(40,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat), 
             data = survey_data_combined[which(survey_data_combined$snowcrab==TRUE),], alpha =.33, col = "cornflowerblue") + 
  facet_wrap(year ~ ., ncol = 4)+
  theme(text=element_text(size=14, family="serif"))
#plot amounts where sea cucumbers caught for RV
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-68,-57), ylim = c(42,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat, col = log(std.WGT)), 
             data = survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$snowcrab==FALSE),], alpha =.5) + 
  facet_wrap(year ~ ., ncol = 4)+
  labs(col ="log(CPUE\n[kg/square\nkilometre])") +
  theme(text=element_text(size=16, family="serif")) + 
  scale_colour_gradientn(colours = viridis(100))
#plot amounts where sea cucumbers caught for snowcrab
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-68,-57), ylim = c(42,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat, col = log(std.WGT)), 
             data = survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0&survey_data_combined$snowcrab==TRUE),], alpha =.5) + 
  facet_wrap(year ~ ., ncol = 3)+
  labs(col ="log(CPUE\n[kg/square\nkilometre])") +
  theme(text=element_text(size=16, family="serif")) + 
  scale_colour_gradientn(colours = viridis(100))

################################EXAMINE ENVIRONMENTAL COVARIATES#######################################

#correlation plot for temperature layers alone
cor <- cor(survey_data_combined[grep("Temp",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Temp.", "Temp. Overall Average","Temp. Minimum", 
                   "Temp. Maximum", "Temp. Range", "Temp. Average Minimum", "Temp. Average Maximum"
                   , "Temp. Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")
#correlation plot for salinity layers alone
cor <- cor(survey_data_combined[grep("Salinity",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Salinity", "Salinity Overall Average","Salinity Minimum", 
                   "Salinity Maximum", "Salinity Range", "Salinity Average Minimum", "Salinity Average Maximum"
                   , "Salinity Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")
#correlation plot for stress layers alone
cor <- cor(survey_data_combined[grep("Stress",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Stress", "Stress Overall Average","Stress Minimum", 
                   "Stress Maximum", "Stress Range", "Stress Average Minimum", "Stress Average Maximum"
                   , "Stress Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")

#correlation plot with considered predictors
survey_data_combined$RangeStressLog <- log(survey_data_combined$RangeStress)
survey_data_combined$BtmStressBNAMLog <- log(survey_data_combined$BtmStressBNAM)
survey_data_combined$sqrt_DEM_Slope <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$sqrt_DEM_StandardDeviation <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$DEM_log <- log(-survey_data_combined$DEM)
cor <- cor(survey_data_combined[c("BtmTempBNAM","BtmSalinityBNAM","BtmStressBNAMLog","RangeTemp","RangeSalinity",
                                "RangeStressLog","DEM_log","sqrt_DEM_Slope","DEM_StandardDeviation","DEM_RDMV", "DEM_Easterness","DEM_Northerness")])
rownames(cor) <- c("Btm. Temperature", "Btm. Salinity","Log(Btm. Stress)", 
                   "Btm. Temperature Range", "Btm. Salinity Range", "Log(Btm. Stress Range)", 
                   "Log(Depth)", "Square Root of Slope", "Square Root of Rugosity", "RDMV","Easterness", "Northerness")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper", tl.col = "black")

#check VIFs for presence model
b <- glm(presence~BtmTempBNAM+
           BtmSalinityBNAM+
           Log(BtmStressBNAM)+
           Log(RangeStress)+
           RangeTemp+
           RangeSalinity+
           log(-DEM)+
           sqrt(DEM_Slope)+
           DEM_RDMV+
           DEM_Easterness+
           DEM_Northerness+
           snowcrab, family = "binomial", data = survey_data_combined)
vif(b)

#check VIFs for CPUE model
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),]
training <- which(!(survey_data_combined_CPUE$Cluster %in% sampled_cluster))
testing <- which(survey_data_combined_CPUE$Cluster %in% sampled_cluster)
b <- glm(log(std.WGT)~BtmTempBNAM+
           BtmSalinityBNAM+
           Log(BtmStressBNAM)+
           Log(RangeStress)+
           RangeTemp+
           RangeSalinity+
           log(-DEM)+
           sqrt(DEM_Slope)+
           DEM_RDMV+
           DEM_Easterness+
           DEM_Northerness+
           snowcrab, family = "gaussian", data = survey_data_combined[which(survey_data_combined$year<2016&
                                  survey_data_combined$std.WGT>0&is.na(survey_data_combined$std.WGT)==FALSE),])
vif(b)

#load-in and crop DEM raster
DEM = raster("C:/Users/natha/Research/Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
#plot of log-transformed depth and where tows conducted
imagePlot(log(-mask(crop(DEM, shp), shp)), breaks = c(2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8, 8.5, 9), 
          col = c(paletteer::paletteer_c("ggthemes::Classic Blue", 14)), xlab = "Longitude", 
          ylab = "Latitude", main = "Log(Depth [in metres])", family = "serif", cex.lab = 1.5, 
          cex.axis = 1.5, cex.main = 1.5)
points(survey_data_combined$mid.lon, survey_data_combined$mid.lat, col = "orange", cex = 0.065)
#add where sea cucumbers are
points(survey_data_combined$mid.lon[which(survey_data_combined$presence == TRUE)], 
       survey_data_combined$mid.lat[which(survey_data_combined$presence == TRUE)], 
       col = "white", cex = 0.065)
shp <- spTransform(shp, crs("+datum=WGS84 +proj=utm +zone=20T +units=km")) #get it into UTM for later plot use

#plot the bottom salinities geographically (using middle) (and separate by year)
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + 
  coord_map(xlim = c(-68,-57), ylim = c(42,48)) + 
  geom_point(size = 0.9, aes(x = mid.lon, y = mid.lat, col = BtmSalinityBNAM),
             data = survey_data_combined, alpha =.5) + 
  scale_colour_gradientn(colours = viridis(100), limits = c(30,36)) + facet_wrap(year ~ ., ncol = 4)+
  labs(col = "Bottom\nSalinity (psu)") + theme(text=element_text(size=16, family="serif"))

#plot the bottom temp. geographically (using middle) (and separate by year)
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + 
  coord_map(xlim = c(-68,-57), ylim = c(42,48)) + 
  geom_point(size = 0.9, aes(x = mid.lon, y = mid.lat, col = BtmTempBNAM),
             data = survey_data_combined, alpha =.5) + 
  scale_colour_gradientn(colours = viridis(100)) + facet_wrap(year ~ ., ncol = 4)+
  labs(col = "Bottom\nTemperature\n(°C)") + theme(text=element_text(size=16, family="serif"))

#plot log(Bottom Stress)
ggplot() + geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group)) + coord_map(xlim = c(-70,-56.5), ylim = c(40,48)) + 
  geom_point(aes(x = mid.lon, y = mid.lat, col = log(BtmStressBNAM)), 
             data = survey_data_combined, alpha =.5) +scale_colour_gradientn(colours = viridis(100)) + 
  labs(col = "Log(Bottom Stress)") + geom_polygon(data = shp, mapping = aes(x = long, y = lat, group = group), fill = NA, col = "white") + 
  theme_dark() + theme(text=element_text(size=16, family="serif"))

#########################################WEIGHT PER UNIT EXAMINATION#######################################

#load shape file for reserves
fishing_area_reserves <- shapefile("Data/All Fishing Areas with Reserves/Fishing_Areas_2018.shp")
fishing_area_reserves <- spTransform(fishing_area_reserves, 
                                     crs("+proj=longlat +datum=NAD83 +no_defs"))
#subset it to our spatial domain
fishing_area_reserves <- fishing_area_reserves[-which(fishing_area_reserves$Region == "SWNB"|fishing_area_reserves$Region=="4X"),]
#grab coordinates of tows
coordinates <- SpatialPoints(cbind(survey_data_combined$mid.lon, survey_data_combined$mid.lat))
proj4string(coordinates) <- "+proj=longlat +datum=NAD83 +no_defs"
#use the coordinates and shape file to create a  division column in the RV data frame
survey_data_combined$Type <- unlist(over(coordinates, as(fishing_area_reserves[,"Type"],"SpatialPolygonsDataFrame")))

#function for adding n to each boxplot
f <- function(x){
  return(data.frame(y = 1.12, label = paste0("n = ",length(x))))
}

###plot this information temporally only with NAFO divisions separated (ignoring RV survey tows)
ggplot(data=survey_data_combined[which(survey_data_combined$mid.lon>-63&survey_data_combined$snowcrab == TRUE&
                                       is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
       mapping=aes(x = as.factor(year), y=weight_per_unit))+geom_boxplot(fill = "white") +
  coord_cartesian(ylim = c(0,1.125)) + ylab("kg/sea cucumber caught") + xlab("year") +
  theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))+
  geom_smooth(method = "loess", aes(group = 1)) + 
  stat_summary(fun=mean, geom="point", shape = 3) +
  stat_summary(fun.data=f, geom="text", col="blue", size = 5, family = "serif", angle = 90)+
  facet_wrap(~NAFO)

#plot this information temporally with all divisions combined (ignoring RV survey tows)
ggplot(data=survey_data_combined[which(survey_data_combined$mid.lon>-63&survey_data_combined$snowcrab == TRUE&
                                       is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
       mapping=aes(x = as.factor(year), y=weight_per_unit))+geom_boxplot(fill = "white") +
  coord_cartesian(ylim = c(0,1.15)) + ylab("kg/sea cucumber caught") + xlab("year") +
  theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))+
  geom_smooth(method = "loess", aes(group = 1)) + 
  stat_summary(fun=mean, geom="point", shape = 3) +
  stat_summary(fun.data=f, geom="text", col="blue", size = 5, family = "serif", angle = 90)

###plot this information temporally (ignoring RV survey tows) and consider fishing zones
survey_data_combined$FishingArea <- ifelse(survey_data_combined$Type == "Fishing Area"&is.na(survey_data_combined$Type)==FALSE, "Fishing Areas","Elsewhere")
survey_data_combined$FishingArea <- factor(survey_data_combined$FishingArea, levels = c("Fishing Areas","Elsewhere"))
ggplot(data=survey_data_combined[which(survey_data_combined$mid.lon>-63&survey_data_combined$snowcrab == TRUE&
                                       is.na(survey_data_combined$std.WGT)==FALSE&survey_data_combined$std.WGT!=0),], 
       mapping=aes(x = as.factor(year), y=weight_per_unit))+geom_boxplot(fill = "white") +
  coord_cartesian(ylim = c(0,1.125)) + ylab("kg/sea cucumber caught") + xlab("year") +
  theme(text=element_text(size=16, family="serif"), axis.text.x = element_text(angle = 90))+
  geom_smooth(method = "loess", aes(group = 1)) + 
  stat_summary(fun=mean, geom="point", shape = 3) +
  stat_summary(fun.data=f, geom="text", col="blue", size = 5, family = "serif", angle = 90)+
  facet_wrap(~FishingArea)