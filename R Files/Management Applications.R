#load necessary libraries in
library(raster)
library(ggplot2)
library(viridis)
library(dplyr)
library(reshape)
library(sf)
source("R Files/Helper Functions.R")

#set starting year to consider (i.e., will go from this year to 2019)
year_start = 2012

#grab NS, New Brunswick, PEI, NFLD maps for later plots
map <- raster::getData(country = "CAN", level = 1)
map <- map[which(map$NAME_1 == "Nova Scotia"|map$NAME_1 == "Prince Edward Island"|map$NAME_1 == "New Brunswick"|map$NAME_1 == "Newfoundland and Labrador"),]
map <- spTransform(map, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#grab Maine map for later plots
map1 <- raster::getData(country = "USA", level = 1)
map1 <- map1[which(map1$NAME_1 == "Maine"),]
map1 <- spTransform(map1, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")
shp <- spTransform(shp, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
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

#load in commercial data and toss out tows that aren't in one of the six groups (or are beyond 2019)
Seacuke_fishery_data_wgroupings <- read.csv("D:/Seacuke_fishery_data_wgroupings.csv")
Seacuke_fishery_data_wgroupings <- Seacuke_fishery_data_wgroupings[-which(is.na(Seacuke_fishery_data_wgroupings$group)|Seacuke_fishery_data_wgroupings$yr>2019),]

#extract the locations and convert to UTM-coordinates
XY = data.frame(Seacuke_fishery_data_wgroupings$lon, Seacuke_fishery_data_wgroupings$lat)
names(XY) = c("X","Y")
latLon = SpatialPoints(XY, proj4string = CRS(as.character(NA)))
proj4string(latLon) <- CRS("+proj=longlat +datum=WGS84") 
Seacuke_fishery_data_wgroupings[c("UTMX","UTMY")] <- as.data.frame(spTransform(latLon,paste("+proj=utm +zone=20T +units=km",sep = "")))

#create a time series of rasters (i.e., a raster stack) from our predictions
presence_prediction_df <- st_set_geometry(presence_prediction_df, NULL) #remove starve geometry
prediction_stack <- raster_create(presence_prediction_df, "combined_predict", year_start = year_start)
#spatio-temporally (by calendar year) intersect predictions with commercial data
for (i in 1:nlayers(prediction_stack))
{
  years_tows <- which(Seacuke_fishery_data_wgroupings$yr==(year_start=1+i))
  Seacuke_fishery_data_wgroupings[years_tows,"model"] <- extract(prediction_stack[[i]], 
                                                                 Seacuke_fishery_data_wgroupings[years_tows, c("UTMX","UTMY")])
}

#histograms for spatio-temporally intersected model predictions... by calendar year
ggplot(aes(x = model), data = Seacuke_fishery_data_wgroupings) + geom_histogram() + facet_wrap(~yr) + 
  xlab("Model predicted log(survey index [in kilograms/kilometre squared])") + theme(text=element_text(size=14, family="serif"))
#histograms for spatio-temporally intersected model predictions... collapsed as yearly trend wasn't interesting
ggplot(aes(x = model), data = Seacuke_fishery_data_wgroupings) + geom_histogram() + 
  xlab("Model predicted log(survey index [in kilograms/kilometre squared])") + theme(text=element_text(size=14, family="serif"))
#95% of the spatially-intersected prediction values are above 4.6...
#"analysis of fishery data suggested 4.6 to be a quality cutoff for defining quality habitat"
quantile(Seacuke_fishery_data_wgroupings$model, na.rm = TRUE, prob = c(0.05))

#apply threshold for quality habitat
prediction_stack_quality <- mask(prediction_stack, prediction_stack<4.6, maskvalue = T)
prediction_stack_quality_plot <- na.omit(as.data.frame(prediction_stack_quality, xy = TRUE)) %>% melt(id.vars = c('x','y'))
levels(prediction_stack_quality_plot$variable) <- seq(year_start, 2019, by = 1)
#intersect quality habitat with fishing zones
prediction_stack_quality_fishing <- mask(prediction_stack_quality, fishing_area)
prediction_stack_quality_plot_fishing <- na.omit(as.data.frame(prediction_stack_quality_fishing, xy = TRUE)) %>% melt(id.vars = c('x','y'))
levels(prediction_stack_quality_plot_fishing$variable) <- seq(year_start, 2019, by = 1)
#remove fishing zones from quality habitat
prediction_stack_quality_no_fishing <- mask(prediction_stack_quality, fishing_area, inverse = T)
prediction_stack_quality_plot_no_fishing <- na.omit(as.data.frame(prediction_stack_quality_no_fishing, xy = TRUE)) %>% melt(id.vars = c('x','y'))
levels(prediction_stack_quality_plot_no_fishing$variable) <- seq(year_start, 2019, by = 1)

#plot result of applying that threshold spatially
ggplot() + geom_polygon(data = fortify(shp), aes(group = group, x = long, y = lat), col = "black", fill = "grey") + 
  geom_raster(data = prediction_stack_quality_plot, aes(x = x, y = y, fill = value)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + facet_wrap(~variable)  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=18, family="serif")) + coord_cartesian(xlim = c(550, 1050), ylim = c(4800, 5150)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted log(survey index)")
#violin plots to summarize quality habitat by year
ggplot() + geom_violin(data = prediction_stack_quality_plot, aes(x = as.factor(variable), y = value)) + xlab("Year") + 
  ylab("Model log(survey index) prediction") + ggtitle("Quality habitat") + 
  theme(text=element_text(size=18, family="serif"), plot.title = element_text(size=26, vjust = 3), 
        axis.title.y = element_text(vjust=5), axis.title.x = element_text(vjust=-2), plot.margin=unit (c (1,1,1,1), 'cm'))
#plots to further summarize quality habitat by year
prediction_stack_quality_plot %>% group_by(variable) %>% summarize(Mean = mean(value), sd = sd(value)) %>%
  ggplot() + geom_point(aes(x = variable, y = Mean, group = 1), size = 2)+
  theme(text=element_text(size=18, family="serif")) + ylab("Mean log(survey index) prediction for quality habitat") + 
  xlab("Year") + geom_linerange(mapping = aes(x = variable, ymin = Mean - sd, ymax = Mean + sd), col = "black", size = .75) + 
  coord_cartesian(ylim = c(0,9)) + 
  theme(axis.title.y = element_text(vjust=5), axis.title.x = element_text(vjust=-2), plot.margin=unit (c (1,1,1,1), 'cm'))

#violin plots to summarize quality habitat in fished areas
ggplot() + geom_violin(data = prediction_stack_quality_plot_fishing, aes(x = as.factor(variable), y = value)) + xlab("Year") + 
  ylab("Model log(survey index) prediction") + ggtitle("Quality habitat in fished areas") + 
  theme(text=element_text(size=18, family="serif"), plot.title = element_text(size=26, vjust = 3), 
        axis.title.y = element_text(vjust=5), axis.title.x = element_text(vjust=-2), plot.margin=unit (c (1,1,1,1), 'cm'))
#violin plots to summarize quality habitat in non-fished areas
ggplot() + geom_violin(data = prediction_stack_quality_plot_no_fishing, aes(x = as.factor(variable), y = value)) + xlab("Year") + 
  ylab("Model log(survey index) prediction") + ggtitle("Quality habitat in non-fished areas") + 
  theme(text=element_text(size=18, family="serif"), plot.title = element_text(size=26, vjust = 3), 
        axis.title.y = element_text(vjust=5), axis.title.x = element_text(vjust=-2), plot.margin=unit (c (1,1,1,1), 'cm'))

#calculate total quality habitat area for each year
good_area <- cellStats(prediction_stack>4.6, sum)*res(prediction_stack)[1]*res(prediction_stack)[2]
#calculate total quality habitat area open to fishery
good_area_fishery <- cellStats(mask(prediction_stack>4.6, fishing_area), sum)*res(prediction_stack)[1]*res(prediction_stack)[2]
#plot proportion of total quality habitat open to the fishery
data.frame(cbind(seq(year_start,2019, by = 1), good_area_fishery/good_area)) %>%
  ggplot() +  geom_point(aes(x = X1, y = X2, group = 1), size = 1.5) + 
  geom_line(aes(x = X1, y = X2, group = 1), size = 1)+
  theme(text=element_text(size=18, family="serif")) + ylab("Proportion of quality habitat open to the fishery") + 
  xlab("Year") + coord_cartesian(ylim = c(0,0.25)) + 
  theme(axis.title.y = element_text(vjust=5), axis.title.x = element_text(vjust=-2), plot.margin=unit (c (1,1,1,1), 'cm'))