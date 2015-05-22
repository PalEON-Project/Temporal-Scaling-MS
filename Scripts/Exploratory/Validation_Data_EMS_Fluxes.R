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


ems.fluxes <- read.csv(file.path(dir.valid, "EMS_fluxes", "hf004-02-filled.csv"))
summary(ems.fluxes)

ggplot(data=ems.fluxes) + 
	geom_line(aes(y=nee, x=seq.time.90), size=0.1)

ems.fluxes.mo <- aggregate(ems.fluxes$nee, by=list(ems.fluxes$month, ems.fluxes$year), FUN=mean)
names(ems.fluxes.mo) <- c("month", "year", "nee")
summary(ems.fluxes.mo)

ems.fluxes.yr <- aggregate(ems.fluxes.mo$nee, by=list(ems.fluxes.mo$year), FUN=mean)
names(ems.fluxes.yr) <- c("year", "nee")
summary(ems.fluxes.yr)

plot(nee ~ year, data=ems.fluxes.yr, type="l")


#---------------------------------
flux.year <- read.csv(file.path(dir.valid, "..", "Flux_All_Year.csv"))
summary(flux.year[flux.year$Site=="HarvardForest1",])
