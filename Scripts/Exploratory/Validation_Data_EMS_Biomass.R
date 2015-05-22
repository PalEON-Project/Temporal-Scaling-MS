# ---------------------------------------
# Extracting, and organizing empirical data for Model Valdiation
# Valdiation Data Set 2: Harvard Forest LTER Datasets
# PalEON MIP: Temporal Scaling Paper
# Christy Rollinson, crollinson@gmail.com
# 27 March 2015
# ---------------------------------------
# Data Sources:
# EMS Biomass
# ---------------------------------------

# Set directories
setwd("~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dir.valid <- "~/Dropbox/PalEON CR/PalEON_MIP_Site/Validation Data/HarvardForest_LTER"

# load some packages
library(ggplot2)

sites <- data.frame(Site=c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL", "PMB" ), 
					Lon =c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83) + 360,
					Lat =c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61)
					)
# coordinates(sites) <- sites[,c("Lon", "Lat")]
# proj4string(sites) <- CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# # sites <- spTransform(sites, CRS('+init=epsg:3175'))
#plot(sites, pch=19)

sec2yr <- 1/(60*60*24*365)


ems.biomass <- read.csv(file.path(dir.valid, "EMS_Biomass", "hf069-13-annual-summ.csv"))
ems.biomass$AGB <- ems.biomass$agwb + (ems.biomass$anpp - ems.biomass$agwi)

for(y in 2000:2010)
summary(ems.biomass)

ems.lai <- read.csv(file.path(dir.valid, "EMS_Biomass", "hf069-02-LAI-site.csv"))
ems.lai$time <- ems.lai$year + ems.lai$doy/365
summary(ems.lai)



ems.swc <- read.csv(file.path(dir.valid, "EMS_Biomass", "hf069-14-swc.csv"))
ems.swc$year.day <- ems.swc$year + (ems.swc$day+ems.swc$time)/365
summary(ems.swc)


ggplot(data=ems.biomass[ems.biomass$site=="ems",]) + 
	geom_line(aes(y=anpp, x=year), size=2) +
	ggtitle("EMS Aboveground Net Primary Productivity") +
	ylab("Mg/Ha/yr")

ggplot(data=ems.biomass[ems.biomass$site=="ems",]) + 
	geom_line(aes(y=agwi, x=year), size=2) +
	ggtitle("EMS Aboveground Woody Biomass Increment") +
	ylab("Mg/Ha/yr")

ggplot(data=ems.biomass[ems.biomass$site=="ems",]) + 
	geom_line(aes(y=agwb, x=year), size=2) +
	ggtitle("EMS Aboveground Woody Biomass") +
	ylab("Biomass Mg/Ha")

ggplot(data=ems.lai[ems.lai$site=="ems" & ems.lai$all.plots==1,]) +
	geom_ribbon(aes(x=time, ymin=all.mean-all.stdev, ymax=all.mean+all.stdev), alpha=0.5) +
	geom_line(aes(x=time, y=all.mean), size=2) +
	ggtitle("EMS LAI") +
	ylab("LAI m2/m2")


ggplot(data=ems.swc) + facet_grid(.~pit) + 
	geom_line(aes(x=year.day, y=swc.15cm), col="orange3") +
	geom_line(aes(x=year.day, y=swc.40cm), col="red3") + 	
	geom_line(aes(x=year.day, y=swc.50cm), col="green3") +
	geom_line(aes(x=year.day, y=swc.90cm), col="black")
	