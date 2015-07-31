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
fig.dir <- fig.base
dat.dir <- dat.base

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
#  B) Tree-Ring Data: 
#     Plot-Level: Harvard Tower, Lyford, Howland
#     ITRDB     : PDL (x5), PHO (x1), PUN (x2)
# ----------------------------------------
dir.lyford  <- file.path(in.base, "Lyford_Data_13m_May2015")
dir.harvard <- file.path(in.base, "HarvardForest")
dir.howland <- file.path(in.base, "HowlandTreeRings")
dir.itrdb   <- file.path(in.base, "ITRDB")
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

# ----------------------------
# Ring Width Data
# ----------------------------
# ------------------
# Plot-Level Tree Rings
# ------------------
# PHA
test <- read.rwl("THIS/IS/A/FOLDER/MINE.rwl")
rings.harvard <- read.rwl(file.path(dir.harvard, "RW", "TP_All.rwl"))
rings.lyford  <- read.rwl(file.path(dir.lyford, "RW", "Lyford_All.rwl"))

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
# ----------------------------
# ----------------------------------------


# -------------------------------------------------
# Processing the Tree-Ring Data
# -------------------------------------------------
# ----------------------------
# Tree-level detrending
#  -- NOTE: Right now we're just doing a 2/3 spline; in the future may want to try something else; but
#           Neil thought this sounded fine for now
# ----------------------------
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
# ----------------------------

# ----------------------------
# Merging plot-level data together to do relative importance weighting
# ----------------------------
harvard.stack <- stack(harvard.detrend)[c(2,1)]
names(harvard.stack) <- c("CoreID", "RWI")
harvard.stack$Year   <- as.numeric(row.names(harvard.detrend))
harvard.stack$TreeID <- as.factor(substr(harvard.stack$CoreID, 1, 6))
harvard.stack$PlotID <- as.factor(substr(harvard.stack$CoreID, 1, 3))
harvard.stack$Site   <- as.factor("PHA")
summary(harvard.stack)

lyford.stack <- stack(lyford.detrend)[c(2,1)]
names(lyford.stack) <- c("CoreID", "RWI")
lyford.stack$Year   <- as.numeric(row.names(lyford.detrend))
lyford.stack$TreeID <- as.factor(substr(lyford.stack$CoreID, 1, 5))
lyford.stack$PlotID <- as.factor(substr(lyford.stack$CoreID, 1, 3))
lyford.stack$Site   <- as.factor("PHA")
summary(lyford.stack)

howland.stack <- stack(howland.detrend)[c(2,1)]
names(howland.stack) <- c("CoreID", "RWI")
howland.stack$Year   <- as.numeric(row.names(howland.detrend))
howland.stack$TreeID <- as.factor(substr(howland.stack$CoreID, 1, 7))
howland.stack$PlotID <- as.factor(substr(howland.stack$CoreID, 1, 4))
howland.stack$Site   <- as.factor("PHO")
summary(howland.stack)

# Putting all the plot-level ring widths together
rwi.stack <- rbind(harvard.stack, lyford.stack, howland.stack)
summary(rwi.stack)

# Aggregate to the tree level
rwi.stack2 <- aggregate(rwi.stack$RWI, by=list(rwi.stack$Site, rwi.stack$PlotID, rwi.stack$TreeID, rwi.stack$Year), FUN=mean, na.rm=T)
names(rwi.stack2) <- c("Site", "PlotID", "TreeID", "Year", "RWI")
rwi.stack2 <- rwi.stack2[complete.cases(rwi.stack2),] # Get rid of missing values because we don't really care
summary(rwi.stack2)


# Merge in the tree-level data
summary(tree.data[,])

rwi.data <- merge(rwi.stack2, tree.data[,c("TreeID", "PlotID", "Site", "Site2", "Species", "Canopy", "Status", "DBH", "BA.tree")], all.x=T, all.y=F)
rwi.data <- rwi.data[complete.cases(rwi.data$DBH),] # Get rid of antying without a DBH/Basal Area 
summary(rwi.data)

# calculate the relative importance of each stem through time
for(p in unique(rwi.data$PlotID)){
	years <- unique(rwi.data[rwi.data$PlotID==p, "Year"])
	for(y in unique(years)){
		ba.plot <- sum(rwi.data[rwi.data$PlotID==p & rwi.data$Year==y, "BA.tree"])
		rwi.data[rwi.data$PlotID==p & rwi.data$Year==y, "RI"] <- rwi.data[rwi.data$PlotID==p & rwi.data$Year==y, "BA.tree"]/ba.plot
	}
}

summary(rwi.data)
# ----------------------------



# ----------------------------
# Create plot-level "chronologies"
# ----------------------------
# ------------
# Plot-Based Data: do a relative importance-weighted mean of RWI
# ------------
rwi.data$RWI.weight <- rwi.data$RWI * rwi.data$RI
summary(rwi.data)

rwi.plot <- aggregate(rwi.data$RWI.weight, by=list(rwi.data$Site, rwi.data$Site2, rwi.data$PlotID, rwi.data$Year), FUN=sum)
names(rwi.plot) <- c("Site", "Site2", "PlotID", "Year", "RWI")
summary(rwi.plot)
# ------------


# ------------
# ITRDB Data: Just make a standard biweighted chronology because we don't know anything about the trees themselves
# ------------
# PDL
chron.mn008 <- chron(mn008.detrend, biweight=T)
chron.mn008$Year <- as.numeric(row.names(chron.mn008))
chron.mn008$PlotID <- as.factor("mn008")

chron.mn009 <- chron(mn009.detrend, biweight=T)
chron.mn009$Year <- as.numeric(row.names(chron.mn009))
chron.mn009$PlotID <- as.factor("mn009")

chron.mn010 <- chron(mn010.detrend, biweight=T)
chron.mn010$Year <- as.numeric(row.names(chron.mn010))
chron.mn010$PlotID <- as.factor("mn010")

chron.mn017 <- chron(mn017.detrend, biweight=T)
chron.mn017$Year <- as.numeric(row.names(chron.mn017))
chron.mn017$PlotID <- as.factor("mn017")

chron.mn019 <- chron(mn019.detrend, biweight=T)
chron.mn019$Year <- as.numeric(row.names(chron.mn019))
chron.mn019$PlotID <- as.factor("mn019")

chron.pdl <- rbind(chron.mn008, chron.mn009, chron.mn010, chron.mn017, chron.mn019)
names(chron.pdl) <- c("RWI", "n.samples", "Year", "PlotID")
chron.pdl$Site <- as.factor("PDL")
summary(chron.pdl)

# PHO
chron.me029 <- chron(me029.detrend, biweight=T)
chron.me029$Year <- as.numeric(row.names(chron.me029))
chron.me029$PlotID <- as.factor("mn019")

chron.pho <- chron.me029
names(chron.pho) <- c("RWI", "n.samples", "Year", "PlotID")
chron.pho$Site <- as.factor("PHO")
summary(chron.pho)

# PUN
chron.mi006 <- chron(mi006.detrend, biweight=T)
chron.mi006$Year <- as.numeric(row.names(chron.mi006))
chron.mi006$PlotID <- as.factor("mi006")

chron.wi002 <- chron(wi002.detrend, biweight=T)
chron.wi002$Year <- as.numeric(row.names(chron.wi002))
chron.wi002$PlotID <- as.factor("wi002")

chron.pun <- rbind(chron.mi006, chron.wi002)
names(chron.pun) <- c("RWI", "n.samples", "Year", "PlotID")
chron.pun$Site <- as.factor("PUN")
summary(chron.pun)

# Merge all site-level chronologies together
chron.itrdb <- rbind(chron.pdl, chron.pho, chron.pun)
# ------------
# ----------------------------

# ----------------------------
# Format for GAMM models
# ----------------------------
# ------------
# merge chronologies together
crn.all <- rbind(chron.itrdb[,c("Site", "PlotID", "Year", "RWI")], rwi.plot[,c("Site", "PlotID", "Year", "RWI")])
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

# Note: only the met from 1901 & greater is observed (CRUNCEP), so lets not merge in false data
dat.tr <- merge(crn.all, met.yr[met.yr$Year>=1901,], all.x=T, all.y=F)
# dat.tr <- dat.tr[complete.cases(dat.tr),]
summary(dat.tr)

# Save the data
write.csv(dat.tr, file.path(dat.base, "TreeRing_Chronologies.csv"), row.names=F)

# Lets just look at the chronologies by site
pdf(file.path(fig.dir, "TreeRingData.pdf"))
ggplot(dat.tr) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=RWI, color=PlotID)) +
	scale_x_continuous(limits=c(1900,2000)) +
	theme_bw()
dev.off()
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
