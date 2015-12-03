# ----------------------------------------
# Climate Sensitivity & Scale Paper
# Non-linear driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 28 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Script Objective: Create baseline model sensitivity to temperature, 
#    precipitation and CO2.  This becomes the basis
# -------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(raster)
library(ncdf4)
library(dplR)
library(ggplot2); library(grid)
library(car)
library(zoo)
# ----------------------------------------

# ----------------------------------------
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/"
fig.base="Figures/"

# Directory for PRISM meterorology
# NOTE NOTE NOTE: Currently it's annul means; need to change it to growing season
met.dir = "~/Dropbox/PalEON_CR/PRISM_PalEON/annual"
co2.dir = "~/Dropbox/PalEON_CR/phase1a_env_drivers/phase1a_env_drivers_v5/paleon_co2"

# Tree Ring Directories:
in.base = "raw_inputs"
dir.lyford  <- file.path(in.base, "Lyford_Data_20m_Final")
dir.harvard <- file.path(in.base, "HarvardForest")
dir.howland <- file.path(in.base, "HowlandTreeRings")
dir.itrdb   <- file.path(in.base, "ITRDB")



# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- file.path(fig.base)
dat.dir <- file.path(dat.base, "TreeRings")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Extract Met data for the plot locations
# ----------------------------------------
# ---------------
# Extract CO2
# ---------------
file.co2 <- nc_open(file.path(co2.dir, "paleon_annual_co2.nc"))
co2 <- data.frame(Year=850:2010, CO2=ncvar_get(file.co2, "co2"))
nc_close(file.co2)

summary(co2)
# ---------------

# ---------------
# Load PRISM Data
# ---------------
# Precipitation
ppt.files <- dir(file.path(met.dir, "ppt"))
ppt.yrs <- substr(ppt.files, nchar(ppt.files)-6, nchar(ppt.files)-3)

ppt <- stack(file.path(met.dir, "ppt", ppt.files))
names(ppt) <- ppt.yrs
ppt

# Temperature
tmean.files <- dir(file.path(met.dir, "tmean"))
tmean.yrs <- substr(tmean.files, nchar(tmean.files)-6, nchar(tmean.files)-3)

tmean <- stack(file.path(met.dir, "tmean", tmean.files))
names(tmean) <- tmean.yrs
tmean
# ---------------

# ---------------
# Load plot locations and extract info
# ---------------
plots <- read.csv(file.path(in.base, "TreeRing_Locations.csv"))
coordinates(plots) <- c("Lon..W.", "Lat..N.")
summary(plots)

# plot(tmean[[1]])
# plot(plots, add=T, pch=19)

plots.clim           <- stack(data.frame(extract(ppt, plots)))
names(plots.clim)    <- c("ppt.ann", "Year")
plots.clim$tmean.ann <- stack(data.frame(extract(tmean, plots)))[,1]
plots.clim$Year <- as.numeric(substr(plots.clim$Year, 2, 5))
plots.clim$PlotID <- plots$Plot
plots.clim <- merge(plots.clim, co2, all.x=T, all.y=F)
summary(plots.clim)
# ---------------

# ----------------------------------------


# -------------------------------------------------
# Formatting the tree data for the plot-based sampling
# -------------------------------------------------
# ----------------------------
# A) Tree Data
#  Note: lack of consistency among sites mean we have to do this piece-meal rather than smoothly all together;
#        this includes just going ahead and manually naming columns in the read.csv stage to make life a lot easier
# ----------------------------

# ------------
# A.1) Lyford Plot
# ------------
trees.lyford <- read.csv(file.path(dir.lyford, "LyfordAllPlots.csv"), skip=3, col.names=c("PlotID", "Tree", "Species", "Canopy", "Status", "DBH", "Distance", "Azimuth", "Tag.Harvard", "TopDist", "TopAz", "Decay"))
summary(trees.lyford)

# Doing some recoding to make things consistent across data sets
trees.lyford$Site        <- as.factor("PHA")
trees.lyford$Site2       <- as.factor(substr(trees.lyford$PlotID,1,2))
trees.lyford$Tree        <- as.factor(ifelse(nchar(trees.lyford$Tree)==1, paste0("0", trees.lyford$Tree), trees.lyford$Tree))
trees.lyford$TreeID      <- as.factor(paste0(trees.lyford$PlotID, trees.lyford$Tree))
trees.lyford$Canopy      <- recode(trees.lyford$Canopy, "'codominant'='C'; 'dominant'='D'; 'intermediate'='I'; 'supressed'='S'; 'suppressed'='S'")
trees.lyford$Status       <- recode(trees.lyford$Status, "'alive'='Li'; 'dead'='Sn'")
trees.lyford$Decay        <- as.factor(trees.lyford $Decay)
summary(trees.lyford)

# ------------

# ------------
# A.2) Havard (Tower Plot)
# ------------
trees.harvard <- read.csv(file.path(dir.harvard, "TowerAllPlots.csv"), skip=3, col.names=c("PlotID", "Tree", "Species", "Canopy", "Status", "DBH", "Distance", "Azimuth", "Tag.Harvard", "Distance.Top", "Azimuth.Top", "Decay"))
summary(trees.harvard)

# Doing some recoding to make things consistent across data sets
trees.harvard$Site        <- as.factor("PHA")
trees.harvard$Site2       <- as.factor(substr(trees.harvard$PlotID,1,2))
trees.harvard$Tree        <-  as.factor(ifelse(nchar(paste(trees.harvard$Tree))==1, paste0("00", trees.harvard$Tree),
  	                                    ifelse(nchar(paste(trees.harvard$Tree))==2, paste0("0", trees.harvard$Tree), 
  	                                    trees.harvard$Tree)))
trees.harvard$TreeID      <- as.factor(paste0(trees.harvard$PlotID, trees.harvard$Tree))
trees.harvard$Canopy      <- recode(trees.harvard$Canopy, "'codominant'='C'; 'dominant'='D'; 'intermediate'='I'; 'supressed'='S'; 'suppressed'='S'")
trees.harvard$Status       <- recode(trees.harvard$Status, "'dead'='Sn'")
trees.harvard$Decay        <- as.factor(trees.harvard$Decay)
summary(trees.harvard)
# ------------


# ------------
# A.3) Howland
# ------------
trees.howland <- read.csv(file.path(dir.howland, "HOWall_FieldData.csv"), col.names=c("PlotID", "Tree", "Species", "DBH", "Distance", "Azimuth", "Canopy", "Status", "Decay", "Lat", "Lon"), na.strings="")
summary(trees.howland)

# Doing some recoding to make things consistent across data sets
trees.howland$Site        <- as.factor("PHO")
trees.howland$Site2       <- as.factor(substr(trees.howland$PlotID,1,3))
trees.howland$Tree        <-  as.factor(ifelse(nchar(paste(trees.howland$Tree))==1, paste0("00", trees.howland$Tree),
  	                                    ifelse(nchar(paste(trees.howland$Tree))==2, paste0("0", trees.howland$Tree), 
  	                                    trees.howland$Tree)))
trees.howland$TreeID      <- as.factor(paste0(trees.howland$PlotID, trees.howland$Tree))
trees.howland$Canopy      <- recode(trees.howland$Canopy, "'S-I'='I'")
trees.howland$Species     <- recode(trees.howland$Species, "'PCRU'='PIRU'; 'PCRu'='PIRU'")
trees.howland$Decay       <- as.factor(trees.howland$Decay)
summary(trees.howland)
# ------------

# ------------
# Merging Datasets together
# ------------
tree.data <- merge(trees.lyford, trees.harvard, all.x=T, all.y=T)
tree.data <- merge(tree.data, trees.howland, all.x=T, all.y=T)
summary(tree.data)

# Figuring out the basal area of each tree (for relative importance calculation later on)
tree.data$BA.tree <- pi*(tree.data$DBH/2)^2
summary(tree.data)
# ------------
# ----------------------------
# -------------------------------------------------

# -------------------------------------------------
# Reading in and formatting the ring widths
# See Peters et al 2015 for Detrending Info
# -------------------------------------------------
# ------------------
# Plot-Level Tree Rings
# ------------------
# PHA
rings.harvard <- read.rwl(file.path(dir.harvard, "RW", "TP_All.rwl"))
rings.lyford <- read.rwl(file.path(dir.lyford, "RW", "combined", "LF_All.rwl"))

# PHO
rings.howland <- read.rwl(file.path(dir.howland, "RawRW_Data_by_sp", "Howland_All.rwl"))
# ------------------

# ------------------
# ITRDB Tree Rings
# ------------------
# PDL
dir(file.path(dir.itrdb, "PDL"), ".rwl")
mn008 <- read.rwl(file.path(dir.itrdb, "PDL", "mn008.rwl.txt"))
mn009 <- read.rwl(file.path(dir.itrdb, "PDL", "mn009.rwl.txt"))
mn010 <- read.rwl(file.path(dir.itrdb, "PDL", "mn010.rwl.txt"))
mn017 <- read.rwl(file.path(dir.itrdb, "PDL", "mn017.rwl.txt"))
mn019 <- read.rwl(file.path(dir.itrdb, "PDL", "mn019.rwl.txt"))

# PHO
dir(file.path(dir.itrdb, "PHO"), ".rwl")
me029 <- read.rwl(file.path(dir.itrdb, "PHO", "me029.rwl.txt"))

# PUN
dir(file.path(dir.itrdb, "PUN"), ".rwl")
mi006 <- read.rwl(file.path(dir.itrdb, "PUN", "mi006.rwl.txt"))
wi002 <- read.rwl(file.path(dir.itrdb, "PUN", "wi002.rwl.txt"))
# ------------------

# ------------------
# Detrending -- Spline (Conservative)
# ------------------
# Plot-Level Data
harvard.detrend <- detrend(rings.harvard, method="Spline")
lyford.detrend  <- detrend(rings.lyford, method="Spline")
howland.detrend <- detrend(rings.howland, method="Spline")

# ITRDB data
## PDL
mn008.detrend   <- detrend(mn008, method="Spline")
mn009.detrend   <- detrend(mn009, method="Spline")
mn010.detrend   <- detrend(mn010, method="Spline")
mn017.detrend   <- detrend(mn017, method="Spline")
mn019.detrend   <- detrend(mn019, method="Spline")

## PHO
me029.detrend   <- detrend(me029, method="Spline")

# PUN
mi006.detrend   <- detrend(mi006, method="Spline")
wi002.detrend   <- detrend(wi002, method="Spline")
# ------------------

# -------------------------------------------------

# -------------------------------------------------
# Putting tree rings & met data together
# -------------------------------------------------
# -----------------
# Binding the tree rings together
# -----------------
# PBL

# PDL
mn008.stack <- stack(mn008.detrend)[,c(2,1)]
names(mn008.stack) <- c("CoreID", "RWI")
mn008.stack$RW     <- stack(mn008)[,1]
mn008.stack$Year   <- as.numeric(row.names(mn008.detrend))
mn008.stack$TreeID <- as.factor(substr(mn008.stack$CoreID, 1, 5))
mn008.stack$PlotID <- as.factor("MN008")
mn008.stack$Site   <- as.factor("PDL")
summary(mn008.stack)

mn009.stack <- stack(mn009.detrend)[,c(2,1)]
names(mn009.stack) <- c("CoreID", "RWI")
mn009.stack$RW     <- stack(mn009)[,1]
mn009.stack$Year   <- as.numeric(row.names(mn009.detrend))
mn009.stack$TreeID <- as.factor(substr(mn009.stack$CoreID, 1, 5))
mn009.stack$PlotID <- as.factor("MN009")
mn009.stack$Site   <- as.factor("PDL")
summary(mn009.stack)

mn010.stack <- stack(mn010.detrend)[,c(2,1)]
names(mn010.stack) <- c("CoreID", "RWI")
mn010.stack$RW     <- stack(mn010)[,1]
mn010.stack$Year   <- as.numeric(row.names(mn010.detrend))
mn010.stack$TreeID <- as.factor(substr(mn010.stack$CoreID, 1, 5))
mn010.stack$PlotID <- as.factor("MN010")
mn010.stack$Site   <- as.factor("PDL")
summary(mn010.stack)

mn017.stack <- stack(mn017.detrend)[,c(2,1)]
names(mn017.stack) <- c("CoreID", "RWI")
mn017.stack$RW     <- stack(mn017)[,1]
mn017.stack$Year   <- as.numeric(row.names(mn017.detrend))
mn017.stack$TreeID <- as.factor(substr(mn017.stack$CoreID, 1, 5))
mn017.stack$PlotID <- as.factor("MN017")
mn017.stack$Site   <- as.factor("PDL")
summary(mn017.stack)

mn019.stack <- stack(mn019.detrend)[,c(2,1)]
names(mn019.stack) <- c("CoreID", "RWI")
mn019.stack$RW     <- stack(mn019)[,1]
mn019.stack$Year   <- as.numeric(row.names(mn019.detrend))
mn019.stack$TreeID <- as.factor(substr(mn019.stack$CoreID, 1, 5))
mn019.stack$PlotID <- as.factor("MN019")
mn019.stack$Site   <- as.factor("PDL")
summary(mn019.stack)

pdl.rwl <- rbind(mn008.stack, mn009.stack, mn010.stack, mn017.stack, mn019.stack)

# PHA
harvard.stack <- stack(harvard.detrend)[,c(2,1)]
names(harvard.stack) <- c("CoreID", "RWI")
harvard.stack$RW     <- stack(rings.harvard)[,1]
harvard.stack$Year   <- as.numeric(row.names(harvard.detrend))
harvard.stack$TreeID <- as.factor(substr(harvard.stack$CoreID, 1, 6))
harvard.stack$PlotID <- as.factor(substr(harvard.stack$CoreID, 1, 3))
harvard.stack$Site   <- as.factor("PHA")
summary(harvard.stack)

lyford.stack <- stack(lyford.detrend)[,c(2,1)]
names(lyford.stack) <- c("CoreID", "RWI")
lyford.stack$RW     <- stack(rings.lyford)[,1]
lyford.stack$Year   <- as.numeric(row.names(lyford.detrend))
lyford.stack$TreeID <- as.factor(substr(lyford.stack$CoreID, 1, 5))
lyford.stack$PlotID <- as.factor(substr(lyford.stack$CoreID, 1, 3))
lyford.stack$Site   <- as.factor("PHA")
summary(lyford.stack)

pha.rwl <- rbind(harvard.stack, lyford.stack)

# PHO
howland.stack <- stack(howland.detrend)[,c(2,1)]
names(howland.stack) <- c("CoreID", "RWI")
howland.stack $RW     <- stack(rings.howland)[,1]
howland.stack$Year   <- as.numeric(row.names(howland.detrend))
howland.stack$TreeID <- as.factor(substr(howland.stack$CoreID, 1, 7))
howland.stack$PlotID <- as.factor(substr(howland.stack$CoreID, 1, 4))
howland.stack$Site   <- as.factor("PHO")
summary(howland.stack)

me029.stack <- stack(me029.detrend)[,c(2,1)]
names(me029.stack) <- c("CoreID", "RWI")
me029.stack$RW     <- stack(me029)[,1]
me029.stack$Year   <- as.numeric(row.names(me029.detrend))
me029.stack$TreeID <- as.factor(substr(me029.stack$CoreID, 1, 5))
me029.stack$PlotID <- as.factor("ME029")
me029.stack$Site   <- as.factor("PHO")
summary(me029.stack)

pho.rwl <- rbind(howland.stack, me029.stack)

# PMB

# PUN
mi006.stack <- stack(mi006.detrend)[,c(2,1)]
names(mi006.stack) <- c("CoreID", "RWI")
mi006.stack$RW     <- stack(mi006)[,1]
mi006.stack$Year   <- as.numeric(row.names(mi006.detrend))
mi006.stack$TreeID <- as.factor(substr(mi006.stack$CoreID, 1, 5))
mi006.stack$PlotID <- as.factor("MI006")
mi006.stack$Site   <- as.factor("PUN")
summary(mi006.stack)

wi002.stack <- stack(wi002.detrend)[,c(2,1)]
names(wi002.stack) <- c("CoreID", "RWI")
wi002.stack$RW     <- stack(wi002)[,1]
wi002.stack$Year   <- as.numeric(row.names(wi002.detrend))
wi002.stack$TreeID <- as.factor(substr(wi002.stack$CoreID, 1, 5))
wi002.stack$PlotID <- as.factor("WI002")
wi002.stack$Site   <- as.factor("PUN")
summary(wi002.stack)

pun.rwl <- rbind(mi006.stack, wi002.stack)

# Putting all the plot-level ring widths together
tree.rings <- rbind(pdl.rwl, pha.rwl, pho.rwl, pun.rwl)
summary(tree.rings)

# Convert CoreID field to all uppercase
tree.rings$CoreID <- as.factor(toupper(tree.rings$CoreID))
summary(tree.rings)
# -----------------

# -----------------
# Finding the best-guess age for each core
# For ITRDB samples, making the (probably bad) assumption that all cores go to pith
# -----------------

# Fusing a few cores that look like they got split apart
# Note: this was done based on crashes in the DBH Reconstruction
dim(tree.rings)
tree.rings <- tree.rings[!tree.rings$CoreID=="HOW3014I" | (tree.rings$CoreID=="HOW3014I" & !is.na(tree.rings$RW)),]
tree.rings[tree.rings$CoreID=="HOW3014I", "CoreID"] <- "HOW3014E"
dim(tree.rings)


# Adding in some core metadata
harvard.meta <- read.csv(file.path(dir.harvard, "TowerAllPlotsLogsheet.csv"), na.strings="")
lyford.meta <- read.csv(file.path(dir.lyford, "LyfordPlots_LogsheetAll_FINAL.csv"), na.strings="")
names(harvard.meta)[13] <- c("False")
summary(harvard.meta)
summary(lyford.meta)

# Merging core metadata into 1 object
core.meta <- rbind(harvard.meta, lyford.meta)
names(core.meta)[1:4] <- c("PlotID", "Tree", "Core", "Species")
core.meta$Tree    <- as.factor(core.meta$Tree)
core.meta$Tree    <- as.factor(ifelse(nchar(paste(core.meta$Tree))==1, paste0("00", core.meta$Tree), ifelse(nchar(paste(core.meta$Tree))==2, paste0("0", core.meta$Tree), core.meta$Tree)))
core.meta$TreeID  <- as.factor(paste0(core.meta$PlotID, core.meta$Tree))
core.meta$CoreID  <- as.factor(paste0(core.meta$TreeID, core.meta$Core))
core.meta$Pith.Yr <- core.meta$INNER_YEAR - core.meta$PITH
summary(core.meta)

# Finding the age for each tree 
summary(tree.rings)
length(unique(tree.rings$CoreID))

for(i in unique(tree.rings$CoreID)){
	if(i %in% core.meta$CoreID){
		pith <- core.meta[core.meta$CoreID==i, "Pith.Yr"]
	} else {
		pith <- min(tree.rings[tree.rings$CoreID==i & !is.na(tree.rings$RW), "Year"])
	}
	tree.rings[tree.rings$CoreID==i & !is.na(tree.rings$RW), "Age"] <-  tree.rings[tree.rings$CoreID==i & !is.na(tree.rings$RW), "Year"] - pith
}
summary(tree.rings)
# -----------------

# -----------------
# Coming up with a best-guess DBH for ITRDB trees: mean of sum of all rings widths
# -----------------
summary(tree.rings)
summary(tree.data); 
dim(tree.data)


for(t in unique(tree.rings$TreeID)[!(unique(tree.rings$TreeID) %in% tree.data$TreeID)]){
	# Making a dummy dataframe to insert info into 
	df.temp <- tree.data[1,]; df.temp[,] <- NA
	df.temp$TreeID <- t

	# Getting DBH estimates from all cores
	dbh <- vector()
	for(c in unique(tree.rings[tree.rings$TreeID==t, "CoreID"])){
		dbh <- c(dbh, sum(tree.rings[tree.rings$CoreID==c, "RW"], na.rm=T)*0.1*2)
	}
	# Taking the mean of DBH estimates
	df.temp$DBH <- mean(dbh, na.rm=T)

	# adding our dummy dataframe to the real data
	tree.data <- rbind(tree.data, df.temp)
}
summary(tree.data); 
dim(tree.data)
# -----------------


# -----------------
# Aggregating to the tree level
# -----------------
# Finding the mean RW and RWI
tree.rings2 <- aggregate(tree.rings[,c("RWI", "RW", "Age")], by=list(tree.rings$Site, tree.rings$PlotID, tree.rings$TreeID, tree.rings$Year), FUN=mean, na.rm=F)
names(tree.rings2)[1:4] <- c("Site", "PlotID", "TreeID", "Year")

# Using the Age off of the oldest core
tree.rings2$Age <- aggregate(tree.rings[,"Age"], by=list(tree.rings$Site, tree.rings$PlotID, tree.rings$TreeID, tree.rings$Year), FUN=max, na.rm=F)[,"x"]
# tree.rings2[tree.rings2$Age==-Inf, "Age"] <- NA
summary(tree.rings2)
# -----------------

# -----------------
# Reconstructing DBH, BA, and BAI
# -----------------
for(t in unique(tree.rings2$TreeID)){
	yr.start <- min(tree.rings2[tree.rings2$TreeID==t & !is.na(tree.rings2$RW), "Year"])
	yr.end   <- max(tree.rings2[tree.rings2$TreeID==t & !is.na(tree.rings2$RW), "Year"])
	dbh      <- tree.data[tree.data$TreeID==t, "DBH"]
	tree.rings2[tree.rings2$TreeID==t & tree.rings2$Year==yr.end, "DBH"] <- dbh
	for(y in (yr.end-1):yr.start){
		dbh  <- tree.rings2[tree.rings2$TreeID==t & tree.rings2$Year==(y+1), "DBH"]
		rw   <- tree.rings2[tree.rings2$TreeID==t & tree.rings2$Year==(y+1), "RW" ]
		dbh2 <- dbh - rw*2*.1
		bai  <- pi*((dbh/2)^2) - pi*((dbh2/2)^2)
		
		tree.rings2[tree.rings2$TreeID==t & tree.rings2$Year==(y+1), "BAI"] <- bai
		tree.rings2[tree.rings2$TreeID==t & tree.rings2$Year==(y)  , "DBH"] <- dbh2
	}
}

# Get rid of negative DBH and BAI reconstructions
tree.rings2$DBH <- ifelse(tree.rings2$DBH<0, NA, tree.rings2$DBH)
tree.rings2$BAI <- ifelse(tree.rings2$BAI<0, NA, tree.rings2$BAI)
summary(tree.rings2)
# -----------------

# -----------------
# merging in the met
# -----------------
summary(tree.data)
summary(tree.rings2)
summary(plots.clim)

tree.rings3 <- merge(tree.rings2, tree.data[,c("PlotID", "TreeID", "Species", "Canopy", "Status")], all.x=T, all.y=F)
summary(tree.rings3)

# Adding in Species/PFT
evergreen <- c("PIRE", "PIST", "PIAB", "PIRU", "THOC", "ABBA", "TSCA")
tree.rings3$Species <- as.factor(ifelse(tree.rings3$Site=="PDL","PIRE", paste(tree.rings3$Species)))
tree.rings3$Species <- as.factor(ifelse(tree.rings3$PlotID=="WI002","PIRE", paste(tree.rings3$Species)))
tree.rings3$Species <- as.factor(ifelse(tree.rings3$PlotID=="MI006","TSCA", paste(tree.rings3$Species)))
tree.rings3$Species <- as.factor(ifelse(tree.rings3$PlotID=="ME029","FRNI", paste(tree.rings3$Species)))
tree.rings3$PFT     <- as.factor(ifelse(tree.rings3$Species %in% evergreen, "Evergreen", "Deciduous"))
summary(tree.rings3)

ring.data <- merge(tree.rings3, plots.clim, all.x=T, all.y=F)
summary(ring.data)
# -----------------

# -----------------
# Calculating DBH, BA, & BAI through time
# For ITRDB, make (bad) assumption that DBH is equal to sum of all ring widths
# -----------------
# -----------------


# -----------------
# Doing some temporal smoothing to look at temporal resolution and sensitivity
# -----------------
# A lazy way of adding a Scale factor (this probably could be done more simply with a merge, but this works too)
ring.data.010 <- ring.data
ring.data$Resolution     <- as.factor("t.001")
ring.data.010$Resolution <- as.factor("t.010")
ring.data <- rbind(ring.data, ring.data.010)
summary(ring.data)

vars.all <- c("RW", "RWI", "BAI", "Age", "ppt.ann", "tmean.ann", "CO2")
# Doing the scale by model & site
for(t in unique(ring.data$TreeID)){
	for(v in vars.all){
		temp <- ring.data[ring.data$TreeID==t & ring.data$Resolution=="t.001", v]

		ring.data[ring.data$TreeID==t & ring.data$Resolution=="t.010", v] <- rollapply(temp, FUN=mean, width=10, align="center", fill=NA)
	}
}
summary(ring.data)
# -----------------



write.csv(ring.data, file.path(dat.base, "TreeRing_RingWidths.csv"), row.names=F)
# -----------------
# -------------------------------------------------

