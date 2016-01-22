# ----------------------------------------
# Objective: Format tree ring data so that it can be compared to models
# Christy Rollinson, crollinson@gmail.com
# Date Created: 21 January 2016
# ----------------------------------------
#
# -------------------------
# Workflow
# -------------------------
# 1. Load and format NPP data
#    a. Load tree-level aNPP data (Kg/tree)
#    b. Convert tree biomass to per area
#    c. Aggregate to species-plot (MgC/Ha)
# 2. Load plot metadata & merge with NPP
#    -- Save: TreeRing_NPP_PlotSpecies.csv
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
# ----------------------------------------

# ----------------------------------------
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/"
fig.base="Figures/"

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Tree Ring NPP Directories:
in.base = "raw_inputs/NPP_TreeRings"
# ----------------------------------------


# -------------------------------------------------
# Note: if you haven't run step 1 from script 3a that extracts met data, do 
#       so now & come back.  I was lazy and didn't want redundant code
# -------------------------------------------------

# -------------------------------------------------
# 1. Load and format NPP data
# -------------------------------------------------
# NOTE: Current data from Alex (as of 20 Jan 2016) is only at the site or tree
#       level.  Plot & Species break downs or scripts to do so shoudl be coming
#       soon.  
#
#       In the meanwhile, I'm going to go ahead and do the conversion from trees 
#       to plots here based on my understanding of how things were sampled so that
#       I don't have to wait.  Once I get data from them, I'll use it instead so 
#       that all of PalEON is using the same numbers.
#
# NOTE: Trees were sampled according to the following nested plot scheme
#       10-20 cm DBH: plot radius = 13 m 
#       20+   cm DBH: plot radius = 20 m
# -------------------------------------------------
# -------------------------
# 1.a. Load tree-level aNPP files
# -------------------------
lyford  <- read.csv(file.path(in.base, "Harvard_Lyford", "Lyf.ab.all.csv"))
tower   <- read.csv(file.path(in.base, "Harvard_Tower" , "tow.ab.all.csv"))
howland <- read.csv(file.path(in.base, "Howland"       , "how.ab.all.csv"))

summary(lyford)
summary(tower)
summary(howland)

# Merging the sites together
# Note the columns in the howland file don't line up with those from Harvard Forest
cols.use <- names(howland)[names(howland) %in% names(lyford)]
tree.npp <- rbind(lyford[,cols.use], tower[,cols.use], howland[,cols.use])
summary(tree.npp)

# Adding some indexing that makes more sense to me
# What they call "Site" is really a plot nested within sites; 
tree.npp$PlotID <- tree.npp$Site 
tree.npp$Site   <- as.factor(ifelse(nchar(paste(tree.npp$Site))>=4, substr(tree.npp$PlotID, 1, 3), substr(tree.npp$PlotID, 1, 2)))
summary(tree.npp)

# ####################
# NOTE -- we have trees listed that are missing biomass starting in 1964 -- need to ask about why
# ####################
# -------------------------

# -------------------------
# 1.b. Convert biomass per tree (kg/tree?) to per area (MgC/HA)
#
# Note: This will be called ABI.area (aboveground biomass increment per area)
# Note: because all sites/plots had the same sampling scheme, we can do this with a single step
#
# NOTE: Trees were sampled according to the following nested plot scheme
#       10-20 cm DBH: plot radius = 13 m 
#       20+   cm DBH: plot radius = 20 m
# -------------------------
# To convert from MG/tree to MGC/HA: biomass increment (mg/tree) * tree/HA * 0.5 (MG C)/(MG biomass)
# associate a stem density with each tree (just so there's a record)
# tree/HA = 1/(pi * radius^2) tree/m2 * 1e5 m2/HA
tree.npp$tree.HA <- ifelse(tree.npp$DBH<20, 1/(pi*(13^2)), 1/(pi*(20^2)))*1e4 # trees per hectare

# Convert tree biomass to per area basis
# biomass/tree * trees/HA * MgC/MgBiomass
tree.npp$ABI.area <- tree.npp$annAB * tree.npp$tree.HA * 0.5
tree.npp$AB.area <- tree.npp$AB * tree.npp$tree.HA * 0.5

summary(tree.npp)

# Note: it looks like we have some major outliers here (aNPP of a single tree is ~15% of model NPP)
#       -- This looks to be the work of a few very large QURUs (>60 cm) at Lyford; most is pre-1980
hist(tree.npp$ABI.area)
summary(tree.npp[tree.npp$ABI.area>3,])
summary(tree.npp[tree.npp$ABI.area>2,])

hist(tree.npp[tree.npp$Year>=1990,"ABI.area"])
# -------------------------

# -------------------------
# 1.c. Aggregate to species biomass per plot (MgC/HA)
# -------------------------
spp.npp <- aggregate(tree.npp[,c("BAI", "tree.HA", "AB.area", "ABI.area")], 
                     by=tree.npp[,c("Site", "PlotID", "Species", "Year")], 
                     FUN=sum, na.rm=T)
summary(spp.npp)

# quick QA/QC plots
pdf(file.path(fig.base, "TreeRings_NPP_Species.pdf"))
ggplot(data=spp.npp) + facet_wrap(~PlotID) +
	geom_line(aes(x=Year, y=ABI.area, color=Species), size=1) +
	scale_y_continuous(name=expression(bold(paste("Aboveground Biomass Increment (Mg C ha"^"-1","yr"^"-1",")")))) +
	scale_x_continuous(limits=c(1900,2010)) +
	theme_bw()
dev.off()

pdf(file.path(fig.base, "TreeRings_AGB_Species.pdf"))
ggplot(data=spp.npp) + facet_wrap(~PlotID) +
	geom_line(aes(x=Year, y=AB.area, color=Species), size=1) +
	scale_y_continuous(name=expression(bold(paste("Aboveground Biomass (Mg C ha"^"-1",")")))) +
	scale_x_continuous(limits=c(1900,2010)) +
	theme_bw()
dev.off()


# -------------------------

# -------------------------------------------------

# -------------------------------------------------
# 2. Load plot metadata & merge with NPP
# -------------------------------------------------
# -------------------------
# 2.a. Climate (from PRISM)
# -------------------------
climate <- read.csv(file.path("Data", "TreeRing_PRISM_Climate.csv"))
summary(climate)

# merge the data sets together, not taking any climate we don't actually need
spp.npp <- merge(spp.npp, climate, all.x=T, all.y=F)
summary(spp.npp)
# -------------------------

# -------------------------
# 2.b. Plant Functional Type Classifications
# -------------------------
unique(spp.npp$Species)
evergreen <- c("PIRE", "PIST", "PIAB", "PIRU", "THOC", "ABBA", "TSCA", "PCRU")

spp.npp$PFT <- as.factor(ifelse(spp.npp$Species %in% evergreen, "Evergreen", "Deciduous"))
summary(spp.npp)
# -------------------------

# -------------------------
# 2.c. Calculating Fcomp (by biomass) for each species
# -------------------------
for(i in unique(spp.npp$PlotID)){
	for(y in unique(spp.npp[spp.npp$PlotID==i, "Year"])){
		bm.tot <- sum(spp.npp[spp.npp$PlotID==i & spp.npp$Year==y, "AB.area"], na.rm=T)
		spp.npp[spp.npp$PlotID==i & spp.npp$Year==y,"Fcomp"] <- spp.npp[spp.npp$PlotID==i & spp.npp$Year==y, "AB.area"]/bm.tot
	}
}
summary(spp.npp)
# -------------------------

write.csv(spp.npp, file.path("Data", "TreeRing_NPP_PlotSpecies.csv"), row.names=F)
# -------------------------------------------------

