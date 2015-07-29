# ----------------------------------------
# Temporal Scaling Analyses -- Model Comparison with Empirical Data
# Non-linear driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 28 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: How do models responses to environmental drivers compare with empirical climate relationships
# -------------------------
#
# -------------------------
# Data Inputs:
# -------------------------
# 1) Tree Ring Data from PalEON (possibly expand to ITRDB later?)
#     -- PHA (Harvard Forest): Lyford, EMS tower? (Neil)
#     -- PHO (Howland Forest): PalEON (Neil/Alex)
#     -- PUN(?) (UNDERC)     : Willow Creek? (Ross)
# -------------------------
#
# -------------------------
# Processing workflow
# -------------------------
# 1) Read in tree-level data
#    -- Raw Ring Widths (.rwl)
#    -- Tree metadata
# 2) Create plot-level RWI chronology
#	 -- Based on IV-weight individual tree chronologies
#	 -- RWI should be similar/easy to translate to the % change in NPP of the Models
# 3) Structure data for use in GAMMs
#    -- Fields Needed: Year, Site, Plot(?), "NPP", met drivers, 
#    -- Question: Average plots to get 1 site-level chronology or make a special empirical data gamm that 
#                 has nested Plot|Site effects?
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(dplR)
library(ggplot2); library(grid)
library(car)
# ----------------------------------------

# ----------------------------------------
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
# setwd("~/Desktop/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
in.base="raw_inputs"
dat.base="Data/TreeRings_processed"
fig.base="Figures/TreeRings_processed"

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- dat.base
dat.dir <- fig.base

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
#  A) Tree Data: Harvard Tower, Lyford, Howland
#      -- read in & format column names
#      -- recode variables for consistent naming conventions
#      -- add new factors for GAMMs
#  B) Tree-Ring Data: Harvard Tower, Lyford, Howland
# ----------------------------------------
dir.lyford  <- file.path(in.base, "Lyford_Data_13m_May2015")
dir.harvard <- file.path(in.base, "HarvardForest")
dir.howland <- file.path(in.base, "HowlandTreeRings")

# ----------------------------
# A) Tree Data
#  Note: lack of consistency among sites mean we have to do this piece-meal rather than smoothly all together;
#        this includes just going ahead and manually naming columns in the read.csv stage to make life a lot easier
# ----------------------------

# ------------
# A.1) Lyford Plot
# ------------
trees.lyford <- read.csv(file.path(dir.lyford, "LyfordAllPlots.csv"), skip=3, col.names=c("PlotID", "Tree", "Species", "Canopy", "Status", "DBH", "Distance", "Azimuth", "Plot.Lyford", "Tag.Harvard"))
summary(trees.lyford)

# Doing some recoding to make things consistent across data sets
trees.lyford$Site        <- as.factor("PHA")
trees.lyford$Site2       <- as.factor(substr(trees.lyford$PlotID,1,2))
trees.lyford$Tree        <- as.factor(ifelse(nchar(trees.lyford$Tree)==1, paste0("0", trees.lyford$Tree), trees.lyford$Tree))
trees.lyford$TreeID      <- as.factor(paste0(trees.lyford$PlotID, trees.lyford$Tree))
trees.lyford$Canopy      <- recode(trees.lyford$Canopy, "'codominant'='C'; 'dominant'='D'; 'intermediate'='I'; 'supressed'='S'; 'suppressed'='S'")
trees.lyford$Status       <- recode(trees.lyford$Status, "'alive'='Li'; 'dead'='Sn'")
summary(trees.lyford)

# Calculate individual density so we can do weighted IV

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

# Calculate individual density so we can do weighted IV

# ------------


# ------------
# A.3) Howland
# ------------
trees.howland <- read.csv(file.path(dir.howland, "HOWall_FieldData.csv"), col.names=c("PlotID", "Tree", "Species", "DBH", "Distance", "Azimuth", "Canopy", "Status", "Decay.Class", "Lat", "Lon"), na.strings="")
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
trees.howland$Decay.Class <- as.factor(trees.howland$Decay.Class)
summary(trees.howland)

# Calculate individual density so we can do weighted IV

# ------------

# ------------
# Merging Datasets together
# ------------
tree.data <- merge(trees.lyford, trees.harvard, all.x=T, all.y=T)
tree.data <- merge(tree.data, trees.howland, all.x=T, all.y=T)
summary(tree.data)
# ------------
# ----------------------------

# ----------------------------
# Ring Width Data
# ----------------------------
rings.harvard <- read.rwl(file.path(dir.harvard, "RW", "TP_All.rwl"))
rings.lyford  <- read.rwl(file.path(dir.lyford, "RW", "Lyford_All.rwl"))
rings.howland <- read.rwl(file.path(dir.howland, "RawRW_Data_by_sp", "Howland_All.rwl"))
# ----------------------------
# ----------------------------------------


# -------------------------------------------------
# Processing the Tree-Ring Data
# -------------------------------------------------
# ----------------------------
# Tree-level detrending
#  -- NOTE: first pass we're just doing a spline; in the future may want to try something else
# ----------------------------
harvard.detrend <- detrend(rings.harvard, method="Spline")
lyford.detrend  <- detrend(rings.lyford, method="Spline")
howland.detrend <- detrend(rings.howland, method="Spline")
# ----------------------------

# ----------------------------
# Aggregate to the plot level (plot-level chronology)
#   -- NOTE: first pass is just biweights mean to make a chronology;
#            in the future, want to weight by tree importance value in the plot to get a weighted mean
# ----------------------------

# ------------
# Creating indexes of which trees are in each plot; doing this the brute-force lazy way for now
# ------------
# Harvard Forest Tower
unique(substr(names(harvard.detrend),1,3))
harvard1 <- which(substr(names(harvard.detrend),1,3)=="TP1")
harvard2 <- which(substr(names(harvard.detrend),1,3)=="TP2")

# Lyford Plot (Harvard Forest)
unique(substr(names(lyford.detrend),1,3))
lyford1 <- which(substr(names(lyford.detrend),1,3)=="LF1")
lyford2 <- which(substr(names(lyford.detrend),1,3)=="LF2")
lyford3 <- which(substr(names(lyford.detrend),1,3)=="LF3")

# Howland Forest
unique(substr(names(howland.detrend),1,4))
howland1 <- which(substr(names(howland.detrend),1,4)=="HOW1")
howland2 <- which(substr(names(howland.detrend),1,4)=="HOW2")
howland3 <- which(substr(names(howland.detrend),1,4)=="HOW3")
# ------------

# ------------
# Create plot-level chronologies
# ------------
# Harvard Forest Tower
crn.harvard1 <- chron(harvard.detrend[,harvard1], biweght=T)
crn.harvard2 <- chron(harvard.detrend[,harvard2], biweght=T)

# Lyford Plot (Harvard Forest)
crn.lyford1  <- chron(lyford.detrend[,lyford1], biweight=T)
crn.lyford2  <- chron(lyford.detrend[,lyford2], biweight=T)
crn.lyford3  <- chron(lyford.detrend[,lyford3], biweight=T)

# Howland Forest Tower
crn.howland1 <- chron(howland.detrend[,howland1], biweight=T)
crn.howland2 <- chron(howland.detrend[,howland2], biweight=T)
crn.howland3 <- chron(howland.detrend[,howland3], biweight=T)
# ------------
# ----------------------------

# ----------------------------
# Format for GAMM models
# ----------------------------
# ------------
# Bind chronologies together
# ------------
# add in plot IDs
crn.harvard1$PlotID <- "TP1"
crn.harvard2$PlotID <- "TP2"
crn.lyford1 $PlotID <- "LF1"
crn.lyford2 $PlotID <- "LF2"
crn.lyford3 $PlotID <- "LF3"
crn.howland1$PlotID <- "HO1"
crn.howland2$PlotID <- "HO2"
crn.howland3$PlotID <- "HO3"


# add in Years
crn.harvard1$Year <- as.numeric(row.names(crn.harvard1))
crn.harvard2$Year <- as.numeric(row.names(crn.harvard2))
crn.lyford1 $Year <- as.numeric(row.names(crn.lyford1))
crn.lyford2 $Year <- as.numeric(row.names(crn.lyford2))
crn.lyford3 $Year <- as.numeric(row.names(crn.lyford3))
crn.howland1$Year <- as.numeric(row.names(crn.howland1))
crn.howland2$Year <- as.numeric(row.names(crn.howland2))
crn.howland3$Year <- as.numeric(row.names(crn.howland3))

# bind the chronologies together
crn.all <- rbind(crn.harvard1, crn.harvard2, crn.lyford1, crn.lyford2, crn.lyford3, crn.howland1, crn.howland2, crn.howland3)
crn.all$PlotID <- as.factor(crn.all$PlotID)
crn.all$Site   <- as.factor(ifelse(substr(crn.all$PlotID,1,2)=="HO", "PHO", "PHA"))
crn.all$Site2  <- as.factor(substr(crn.all$PlotID, 1, 2))
summary(crn.all)
# ------------

# ------------
# Merge in Met data
# ------------
sec2yr <- 1*60*60*24*365.25 # 1 sec * 60 sec/min * 60 min/hr * 24 hr/day * 365.25 days/yr

met.yr <- read.csv(file.path("Data", "analysis_drivers", "Drivers_Year_GrowingSeason.csv"))
met.yr$CO2 <- met.yr$CO2.yr 
met.yr[,substr(names(met.yr),1,6)=="precip"] <- met.yr[,substr(names(met.yr),1,6)=="precip"]*sec2yr

summary(met.yr)

dat.tr <- merge(crn.all, met.yr, all.x=T, all.y=F)
dat.tr <- dat.tr[complete.cases(dat.tr),]
summary(dat.tr)

write.csv(dat.tr, file.path(dat.base, "TreeRing_Chronologies.csv"), row.names=F)
# ------------
# ----------------------------
# -------------------------------------------------


# # -------------------------------------------------
# # Just running a couple test gams to see if I cam make things work
# # -------------------------------------------------
# predictors <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")

# dat.mod <- dat.tr[,c("Site", "Site2", "PlotID", "Year", paste0(predictors, ".gs"))]
# names(dat.mod)[(ncol(dat.mod)-length(predictors)+1):ncol(dat.mod)] <- predictors
# dat.mod$NPP <- dat.tr$xxxstd
# summary(dat.mod)

# library(mgcv)
# k=4
# gam1 <- gam(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=dat.mod, correlation=corARMA(form=~Year|PlotID, p=1))

# gam2 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, random=list(PlotID=~1, Site2=~1, Site=~1), data=dat.mod, correlation=corARMA(form=~Year|PlotID, p=1))

# gam3 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=dat.mod, correlation=corARMA(form=~Year|PlotID*Site2*Site, p=1))

# gam4 <- gam(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=dat.mod, correlation=corARMA(form=~Year|PlotID|Site2|Site, p=1))
# gam5 <- gam(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=dat.mod)
# gam6 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=dat.mod, correlation=corARMA(form=~Year|PlotID, p=1))


# summary(gam1)
# plot(gam1)

# # -------------------------------------------------
