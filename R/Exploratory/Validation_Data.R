# ---------------------------------------
# Extracting, and organizing empirical data for Model Valdiation
# PalEON MIP: Temporal Scaling Paper
# Christy Rollinson, crollinson@gmail.com
# 27 March 2015
# ---------------------------------------
# Data Sources:
# 1) Paleon Biomass 
# 2) Tree Rings
# 3) Flux Records
# 4) Modern Biomass -- NACP NBCD_2000
# ---------------------------------------

# Set directories
setwd("~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dir.valid <- "~/Desktop/PalEON CR/PalEON_MIP_Site/Validation Data"

# load some packages
library(raster)
library(rgdal)
library(ncdf4)


sites <- data.frame(Site=c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL", "PMB" ), 
					Lon =c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83) + 360,
					Lat =c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61)
					)
coordinates(sites) <- sites[,c("Lon", "Lat")]
proj4string(sites) <- CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# sites <- spTransform(sites, CRS('+init=epsg:3175'))
#plot(sites, pch=19)

sec2yr <- 1/(60*60*24*365)

###################################################################

# ---------------------------------------
# Paleon Settlement Veg Raw Biomass Product
# Variable: Aboveground Biomass
# Spatial Scale: gridded Upper Midwest, (extent, resolution)
# Temporal Scale: single snapshot circa ??? (extent, resolution)
# Data Authors: Simon Goring
# Data Link: https://paleon.geography.wisc.edu/doku.php/data_and_products;settlement_vegetation_biomass
# Working Group Link: https://paleon.geography.wisc.edu/doku.php/working_groups;settlement_vegetation_biomass
# ---------------------------------------
# Read in data
set.veg <- read.csv(file.path(dir.valid, "paleon_biomass", "plss_biomass_alb_v0.9-3.csv"), row.names=1)
# set.veg <- read.csv(file.path(dir.valid, "paleon_biomass", "biomass.output_v0.9-1.csv"), row.names=1)
set.veg$Total <- rowSums(set.veg[,4:ncol(set.veg)])
summary(set.veg)

# Make data spatially explicit
set.veg.grid <- SpatialPixelsDataFrame(points=set.veg[,c("x", "y")], data=set.veg[,4:ncol(set.veg)], proj4string=CRS('+init=epsg:3175'))
gridded(set.veg.grid) <- TRUE
set.veg.grid <- stack(set.veg.grid)
set.veg.grid
# set.veg.grid[set.veg.grid>500] <- NA

plot(set.veg.grid[["Total"]], cex=0.1)
plot(sites, add=T, cex=2, pch=19)
# hist(set.veg[set.veg$Total<500,"Total"])
# ---------------------------------------

###################################################################

# ---------------------------------------
# Paleon Settlement Veg Statistical Biomass Product
# Variable: Aboveground Biomass
# Spatial Scale: gridded Upper Midwest, (extent, resolution)
# Temporal Scale: single snapshot circa ??? (extent, resolution)
# Data Authors: Jun Zhu, Xioping Feng, Chris Paciorek (stats modelers)
# Data Link: https://paleon.geography.wisc.edu/doku.php/data_and_products;settlement_vegetation_biomass
# Working Group Link: https://paleon.geography.wisc.edu/doku.php/working_groups;settlement_vegetation_biomass
# ---------------------------------------
set.veg.stat <- read.csv(file.path(dir.valid, "paleon_biomass", "biomass_prediction-2.csv"))
set.veg.stat$Total <- set.veg.stat$TNE + set.veg.stat$TBD
summary(set.veg.stat)
dim(set.veg.stat)
dim(set.veg)

set.veg.old <- read.csv(file.path(dir.valid, "paleon_biomass", "biomass.output_v0.9-1.csv"), row.names=1)
set.veg.old$Total <- rowSums(set.veg.old[,4:ncol(set.veg.old)])
dim(set.veg.old)
summary(set.veg.old)

plot(set.veg.stat$Total[set.veg.old$Total<500] ~ set.veg.old$Total[set.veg.old$Total<500])

set.veg.stat$x <- set.veg.old$x
set.veg.stat$y <- set.veg.old$y
summary(set.veg.stat)
quantile(set.veg.stat$Total, c(0.025, 0.975))


set.veg.stat.grid <- SpatialPixelsDataFrame(points=set.veg.stat[,c("x", "y")], data=set.veg.stat[,4:ncol(set.veg.stat)], proj4string=CRS('+init=epsg:3175'))
gridded(set.veg.stat.grid) <- TRUE
set.veg.stat.grid <- stack(set.veg.stat.grid)
set.veg.stat.grid

plot(set.veg.stat.grid[["Total"]])
plot(sites, add=T, cex=2, pch=19)
length(set.veg.stat$Total[set.veg.stat$Total>200])
length(set.veg.stat$Total)

sites$Biomass.SetVegStat <- extract(set.veg.stat.grid[["Total"]], sites, method="bilinear")
data.frame(sites)
write.csv(sites, file.path(dir.valid, "SetVegBiomass_MIP_Sites.csv"), row.names=F)
# ---------------------------------------

###################################################################

# ---------------------------------------
# Tree Ring Records
# Variable:
# Spatial Scale: (extent, resolution)
# Temporal Scale: (extent, resolution)
# Authors: Neil Pederson
# ---------------------------------------
library(dlpR)

rwl.dir <- "~/Desktop/PalEON CR/paleon_mip_site/Validation Data/TreeRings/Lyford_Data_13m/RW/Combined"
acru <- read.rwl(file.path(rwl.dir, "Lyford_ACRU.rw"))
beal <- read.rwl(file.path(rwl.dir, "Lyford_BEAL.rw"))
fagr <- read.rwl(file.path(rwl.dir, "Lyford_FAGR.rw"))
havi <- read.rwl(file.path(rwl.dir, "Lyford_HAVI.rw"))
pist <- read.rwl(file.path(rwl.dir, "Lyford_PIST.rw"))
quru <- read.rwl(file.path(rwl.dir, "Lyford_QURU.rw"))
tsca <- read.rwl(file.path(rwl.dir, "Lyford_TSCA.rw"))

quru.detrend <- detrend(quru, method="Spline", make.plot=F)
summary(quru.detrend)

quru.chron <- chron(quru.detrend, biweight=T)
summary(quru.chron)
crn.plot(quru.chron)
write.csv(quru.chron, "Lyford_QURU_Chronology.csv", row.names=T)

# ---------------------------------------

###################################################################

# ---------------------------------------
# Flux Towers
# Variables: H (sensible heat flux), LE (latent heat flux), VPD, SWC1, 
#			 NEE, RE, GPP, CO2top, 
# Spatial Scale: (extent, resolution)
# Temporal Scale: (extent, resolution)
# Authors: 
# ---------------------------------------
# ----------------------------------------------------------------
# Looping through all sites
# ----------------------------------------------------------------
flux.sites <- c("Harvard_Hemlock", "HarvardForest1", "Howland_East_Harvest", "Howland_Main", "Howland_West", "Sylvania", "Willow_Creek")
for(s in 1:length(flux.sites)){
	if(flux.sites[s]=="Harvard_Hemlock"){
		dir.flux <-  file.path(dir.valid, "Flux_Towers", flux.sites[s], "with_gaps")
	} else {
		dir.flux <-  file.path(dir.valid, "Flux_Towers", flux.sites[s], "gap_filled")
	}	
	
	# Not all sites had .nc files, so let just use CSV for simplicity
	files.flux <- dir(dir.flux, ".csv")
	for(i in 1:length(files.flux)){
		flux.file.dummy <- read.csv(file.path(dir.flux, files.flux[i]), header=T, skip=17, nrows=1)
		flux.file <- read.csv(file.path(dir.flux, files.flux[i]), header=F, skip=19, na.strings=c("-9999", "-6999"))
		names(flux.file) <- names(flux.file.dummy)
		summary(flux.file)
			
		temp <- data.frame(Site=flux.sites[s], Year=flux.file$YEAR, DOY=as.factor(flux.file$DOY), HRMIN=flux.file$HRMIN, NEE=flux.file$NEE, GPP=flux.file$GPP, RE=flux.file$RE, CO2=flux.file$CO2, Qh=flux.file$H, Qle=flux.file$LE, SWC=flux.file$SWC1)

	if(s==1 & i==1) flux.raw <- temp else flux.raw <- rbind(flux.raw, temp)	
	}

}

flux.raw <- flux.raw[!flux.raw$Year=="100" & !flux.raw$DOY==366,]
summary(flux.raw)
# summary(flux.raw[flux.raw$Site==flux.sites[1],])
write.csv(flux.raw, file=file.path(dir.valid, "Flux_All_Raw.csv"), row.names=F)

flux.day <- aggregate(flux.raw[,5:ncol(flux.raw)], by=list(flux.raw$Site, flux.raw$Year, flux.raw$DOY), FUN=mean)
names(flux.day)[1:3] <- c("Site", "Year", "Day")
summary(flux.day)

flux.year <- aggregate(flux.day[,4:ncol(flux.day)], by=list(flux.day$Site, flux.day$Year), FUN=mean)
names(flux.year)[1:2] <- c("Site","Year")
summary(flux.year)
write.csv(flux.year, file=file.path(dir.valid, "Flux_All_Year.csv"), row.names=F)


# -------------------
flux.year <- read.csv(file.path(dir.valid, "Flux_All_Year.csv"))
summary(flux.year)
summary(flux.year[flux.year$Site=="HarvardForest1", "GPP"]*sec2yr)

flux.sites
tower.colors <- c("green3", "darkgreen", "skyblue3", "blue2", "turquoise3", "sienna3", "maroon3")

pdf("~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/FluxTower_NEE_Annual.pdf")
plot(NEE ~ Year, data=flux.year, type="n", main="Net Ecosystem Exchange", ylab="NEE, umol/m2/s")
	for(s in 1:length(flux.sites)){
		lines(NEE ~ Year, data=flux.year[flux.year$Site==flux.sites[s],], lwd=2, col=tower.colors[s])
	}
	abline(h=0, col="red", lty="dashed")
	legend("topright", legend=flux.sites, col=tower.colors, lwd=2, bty="n")
dev.off()

pdf("~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/FluxTower_GPP_Annual.pdf")
plot(GPP ~ Year, data=flux.year, type="n", lwd=2, main="Gross Primary Productivity", ylab="GPP, umol/m2/s")
	for(s in 1:length(flux.sites)){
		lines(GPP ~ Year, data=flux.year[flux.year$Site==flux.sites[s],], lwd=2, col=tower.colors[s])
	}
	abline(h=0, col="red", lty="dashed")
	legend("bottomright", legend=flux.sites, col=tower.colors, lwd=2, bty="n")
dev.off()

pdf("~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/FluxTower_RE_Annual.pdf")
plot(RE ~ Year, data=flux.year, type="n", lwd=2, main="Total Ecosystem Respiration", ylab="Total Resp, umol/m2/s")
	for(s in 1:length(flux.sites)){
		lines(RE ~ Year, data=flux.year[flux.year$Site==flux.sites[s],], lwd=2, col=tower.colors[s])
	}
	abline(h=0, col="red", lty="dashed")
	legend("topleft", legend=flux.sites, col=tower.colors, lwd=2, bty="n")
dev.off()
# ---------------------------------------


###################################################################


# ---------------------------------------
# NACP NBCD Biomass Product
# Variables: AGB
# Spatial Scale: (extent, resolution = 30m)
# Temporal Scale: (extent, resolution)
# Authors: 
# ---------------------------------------
biomass.wd <- "~/Desktop/NBCD2000_V2" 
biomass.ne <- raster(file.path(biomass.wd,"PaleonDomain_Merge"))

sites <- spTransform(sites, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

plot(biomass.ne)
plot(sites, add=T, pch=19, col="blue")

sites$nbcd.mean <- extract(biomass.ne, sites, buffer=1000, fun=mean)
sites$nbcd.sd <- extract(biomass.ne, sites, buffer=1000, fun=sd)

summary(sites)
write.csv(sites, file.path(dir.valid, "NBCD_MIP_Sites.csv"), row.names=F)
# ---------------------------------------
