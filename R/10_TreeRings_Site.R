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
library(parallel)
library(mgcv)
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
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/gamms"
fig.base="Figures/gamms"

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- file.path(fig.base, "Sensitivity_TreeRings_Site")
dat.dir <- file.path(dat.base, "Sensitivity_TreeRings_Site")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData.Rdata"))
summary(ecosys)
model.colors

tree.rings <- read.csv("Data/TreeRing_RingWidths.csv")
tree.rings <- tree.rings[complete.cases(tree.rings[,c("tmean.ann", "ppt.ann", "CO2")]),]
tree.rings$Spp.Site <- as.factor(paste(tree.rings$Species, tree.rings$Site, sep=".")) # May want to add back in site 2 and use that instead
summary(tree.rings)

# source('R/0_calculate.sensitivity_TPC_Site.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)

# Read in model color scheme
model.colors
# ----------------------------------------

# -------------------------------------------------
# Settings for the rest of this script
# -------------------------------------------------
# Get rid of Models that have issues
ecosys <- ecosys[!ecosys$Model=="linkages",]

# get rid of site/year combos not in the tree rings
ecosys2 <- ecosys[ecosys$Site %in% unique(tree.rings$Site),]
dim(ecosys); dim(ecosys2)

for(s in unique(ecosys2$Site)){
	s.yrs <- unique(tree.rings[tree.rings$Site==s, "Year"])
	ecosys2 <- ecosys2[!ecosys2$Site==s | (ecosys2$Site==s & ecosys2$Year %in% s.yrs), ]
}
dim(ecosys2)
summary(ecosys2)

# Setting up a loop for 1 m.name, 1 temporal scale
resolutions <- c("t.010") # Note: Big models can't handle t.100 at the site level because there aren't enough data points
extents <- data.frame(Start=c(1895), End=c(2010)) 
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")
k=4
e=1	
# -------------------------------------------------

# ------------------------------------------------------------------------------------
# Tree Rings
# ------------------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Format & Process gamms -- 
# --------------------------------------------------------------------------
response <- c("RW")
tree.rings$Y                <- tree.rings[,response]
tree.rings[,predictors.all] <- tree.rings[,paste0(predictors.all, predictor.suffix)]
tree.rings$Model            <- as.factor("TreeRings")
tree.rings$Extent           <- as.factor("1895-2010")
# tree.rings$Resolution <- as.factor("t.001")
summary(tree.rings)

# tree.rings2 <- tree.rings[tree.rings$Resolution=="t.001" & (tree.rings$PlotID=="MN008" | tree.rings$PlotID=="ME029"), ]
tree.rings2 <- tree.rings[tree.rings$Resolution=="t.001" & complete.cases(tree.rings[,c(predictors.all, response)]), ]
summary(tree.rings2)
# summary(tree.rings2[tree.rings2$Site=="PHO",])

# detrend.only <- 	gam(Y ~ s(Age, by=Spp.Site), data= tree.rings2)
# summary(detrend.only)
# plot(detrend.only, ylim=c(-10,10))
# save(detrend.only, "gamm_TreeRings_DetrendOnly.Rdata")

# ----------------------------------------------------------------------------
# Site Intercepts (Baseline Model)
# ----------------------------------------------------------------------------

source('R/0_calculate.sensitivity_TPC_TreeRings.R', chdir = TRUE)

mod.tr <- paleon.gams.models(data=tree.rings2, response=response, k=k, predictors.all=predictors.all, site.effects=T)

# Adding Tree ID to help index the response prediction
mod.tr$ci.response$TreeID <- tree.rings2$TreeID

mod.tr$ci.terms$x <- as.numeric(mod.tr$ci.terms$x)
summary(mod.tr$ci.terms)
save(mod.tr, file=file.path(dat.dir, "gamm_TreeRings_RW_SiteIntercept.Rdata"))
# -------------------------------------------------

# -------------------------------------------------
# Diagnostic Graphs
# -------------------------------------------------
# m.order <- unique(mod.out$data$Model.Order)
# col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]
m.order <- "Tree Rings"
col.model="darkgreen"

pdf(file.path(fig.dir, "GAMM_TreeRingFit_RW_SiteIntercept.pdf"))
print(
ggplot(data=mod.tr$ci.response[,]) + facet_grid(Site~., scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[,], aes(x=Year, y=Y, color=TreeID), color="gray", alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill= TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color= TreeID), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.tr$ci.response[substr(mod.tr$ci.response$TreeID,1,3)=="125",]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[substr(mod.tr$data$TreeID,1,3)=="125",], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	scale_x_continuous(limits=c(1925,1975)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, "GAMM_TreeRing_DriverSensitivity_RW_SiteIntercept.pdf"))
# for(e in unique(mod.tr$ci.terms$Effect)[1:3]){
print(
# ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +	
ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap( ~ Effect, scales="free_x") + theme_bw() +	
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=1) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0(e, " Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
# }
dev.off()
# -------------------------------------------------

# ----------------------------------------------------------------------------
# Site Curves
# ----------------------------------------------------------------------------
source('R/0_calculate.sensitivity_TPC_TreeRings_Site.R', chdir = TRUE)

mod.tr <- paleon.gams.models(data=tree.rings2, response=response, k=k, predictors.all=predictors.all, comp.effects=F)

# Adding Tree ID to help index the response prediction
mod.tr$ci.response$TreeID <- tree.rings2$TreeID

mod.tr$ci.terms$x <- as.numeric(mod.tr$ci.terms$x)
summary(mod.tr$ci.terms)
save(mod.tr, file=file.path(dat.dir, "gamm_TreeRings_RW_SiteCurves.Rdata"))
# -------------------------------------------------

# -------------------------------------------------
# Diagnostic Graphs
# -------------------------------------------------
# m.order <- unique(mod.out$data$Model.Order)
# col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]
m.order <- "Tree Rings"
col.model="darkgreen"

pdf(file.path(fig.dir, "GAMM_TreeRingFit_RW_SiteCurves.pdf"))
print(
ggplot(data=mod.tr$ci.response[,]) + facet_grid(Site~., scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[,], aes(x=Year, y=Y, color=TreeID), color="gray", alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill= TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color= TreeID), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.tr$ci.response[substr(mod.tr$ci.response$TreeID,1,3)=="125",]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[substr(mod.tr$data$TreeID,1,3)=="125",], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	scale_x_continuous(limits=c(1925,1975)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Site Curves", response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, "GAMM_TreeRing_DriverSensitivity_RW_SiteCurves.pdf"))
# for(e in unique(mod.tr$ci.terms$Effect)[1:3]){
print(
# ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +	
ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap( ~ Effect, scales="free_x") + theme_bw() +	
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Site), size=1) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
# }
dev.off()
# -------------------------------------------------
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Site Curves + Composition Intercepts
# ----------------------------------------------------------------------------

source('R/0_calculate.sensitivity_TPC_TreeRings_Site.R', chdir = TRUE)

mod.tr <- paleon.gams.models(data=tree.rings2, response=response, k=k, predictors.all=predictors.all, comp.effects=T)

# Adding Tree ID to help index the response prediction
mod.tr$ci.response$TreeID <- tree.rings2$TreeID

mod.tr$ci.terms$x <- as.numeric(mod.tr$ci.terms$x)
summary(mod.tr$ci.terms)
save(mod.tr, file=file.path(dat.dir, "gamm_TreeRings_RW_SiteCurves_Comp.Rdata"))
# -------------------------------------------------

# -------------------------------------------------
# Diagnostic Graphs
# -------------------------------------------------
# m.order <- unique(mod.out$data$Model.Order)
# col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]
m.order <- "Tree Rings"
col.model="darkgreen"

pdf(file.path(fig.dir, "GAMM_TreeRingFit_RW_SiteCurves_Comp.pdf"))
print(
ggplot(data=mod.tr$ci.response[,]) + facet_grid(Site~., scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[,], aes(x=Year, y=Y, color=TreeID), color="gray", alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill= TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color= TreeID), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.tr$ci.response[substr(mod.tr$ci.response$TreeID,1,3)=="125",]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[substr(mod.tr$data$TreeID,1,3)=="125",], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	scale_x_continuous(limits=c(1925,1975)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, "GAMM_TreeRing_DriverSensitivity_RW_SiteCurves_Comp.pdf"))
# for(e in unique(mod.tr$ci.terms$Effect)[1:3]){
print(
# ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +	
ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap( ~ Effect, scales="free_x") + theme_bw() +	
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Site), size=1) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0(e, " Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
# }
dev.off()
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Composition Curves + Site Intercepts
# ----------------------------------------------------------------------------

source('R/0_calculate.sensitivity_TPC_TreeRings_Comp.R', chdir = TRUE)

summary(tree.rings2)

# modifying the PFT parameter
# tree.rings2$PFT1 <- as.factor(tree.rings2$PFT)
# tree.rings2$PFT  <- as.factor(paste(tree.rings2$PFT, tree.rings2$Site, sep="."))
summary(tree.rings2)

mod.tr <- paleon.gams.models(data=tree.rings2, response=response, k=k, predictors.all=predictors.all, site.effects=T)

# Adding Tree ID to help index the response prediction
mod.tr$ci.response$TreeID <- tree.rings2$TreeID

# making x numeric again
mod.tr$ci.terms$x <- as.numeric(mod.tr$ci.terms$x)
summary(mod.tr$ci.terms)
save(mod.tr, file=file.path(dat.dir, "gamm_TreeRings_RW_CompCurves.Rdata"))


pdf(file.path(fig.dir, "GAMM_TreeRingFit_RW_CompCurves.pdf"))
print(
ggplot(data=mod.tr$ci.response[,]) + facet_grid(Site~., scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[,], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.tr$ci.response[substr(mod.tr$ci.response$TreeID,1,3)=="125",]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[substr(mod.tr$data$TreeID,1,3)=="125",], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	scale_x_continuous(limits=c(1925,1975)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, "GAMM_TreeRing_DriverSensitivity_RW_CompCurves.pdf"))
# for(e in unique(mod.tr$ci.terms$Effect)[1:3]){
print(
# ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +	
ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap( ~ Effect, scales="free_x") + theme_bw() +	
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=PFT), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=PFT), size=1) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0(e, " Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
# }
dev.off()

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Species Curves (instead of PFT)
# ----------------------------------------------------------------------------

source('R/0_calculate.sensitivity_TPC_TreeRings_Comp.R', chdir = TRUE)

summary(tree.rings2)

# modifying the PFT parameter
tree.rings2$PFT1 <- as.factor(tree.rings2$PFT)
tree.rings2$PFT  <- tree.rings2$Species
summary(tree.rings2)

mod.tr <- paleon.gams.models(data=tree.rings2, response=response, k=k, predictors.all=predictors.all, site.effects=T)

# Adding Tree ID to help index the response prediction
mod.tr$ci.response$TreeID <- tree.rings2$TreeID

# making x numeric again
mod.tr$ci.terms$x <- as.numeric(mod.tr$ci.terms$x)
summary(mod.tr$ci.terms)
save(mod.tr, file=file.path(dat.dir, "gamm_TreeRings_RW_CompCurves_Species.Rdata"))


pdf(file.path(fig.dir, "GAMM_TreeRingFit_RW_CompCurves_Species.pdf"))
print(
ggplot(data=mod.tr$ci.response[,]) + facet_grid(Site~., scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[,], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.tr$ci.response[substr(mod.tr$ci.response$TreeID,1,3)=="125",]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.tr$data[substr(mod.tr$data$TreeID,1,3)=="125",], aes(x=Year, y=Y, color=TreeID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=TreeID), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=TreeID), size=0.35) +
	scale_x_continuous(limits=c(1925,1975)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, "GAMM_TreeRing_DriverSensitivity_RW_CompCurves_Species.pdf"))
# for(e in unique(mod.tr$ci.terms$Effect)[1:3]){
print(
# ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +	
ggplot(data=mod.tr$ci.terms[mod.tr$ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap( ~ Effect, scales="free_x") + theme_bw() +	
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=PFT), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=PFT), size=1) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0(e, " Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
# }
dev.off()

# -------------------------------------------------
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# -------------------------------------------------
# Setting up the data and putting it in a list to run the gamms in parallel
# -------------------------------------------------
source('R/0_calculate.sensitivity_TPC_Site.R', chdir = TRUE)

response="NPP"
ecosys <- ecosys[complete.cases(ecosys[,c(response, predictors.all)]) & ecosys$Year %in% tree.rings2$Year,]
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)

paleon.models <- list()
for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]

	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))

	# Note: Here we're renaming things that had the suffix to just be generalized tair, etc 
	dat.mod <- ecosys[ecosys$Model==m.name, c("Model", "Updated", "Model.Order", "Site", "Year", response, paste0(predictors.all, predictor.suffix), "Evergreen", "Grass")]
	names(dat.mod)[7:(7+length(predictors.all)-1)] <- predictors.all
	
	if(!max(dat.mod[,response], na.rm=T)>0) next # If a variable is missing, just skip over this model for now

# for(r in 1:length(resolutions)){ # Resolution loop

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.gs & temperature distributions
	run.end <- ifelse(substr(m.name,1,3)=="jul", 2009, 2010) # Note: Jules missing 2010, so 
	run.start <- 850
	inc <- 1 # making sure we're always dealign with whole numbers
	yrs <- seq(from=run.end, to=run.start, by=-1)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "Year", response, predictors.all, "Evergreen", "Grass")]

	# Making a note of the extent & resolution
	ext <- as.factor(paste(1850, 2010, sep="-"))
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions)

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,c(response)]),]

	data.temp$Y <- data.temp[,response]

	paleon.models[[paste(m.name)]] <- data.temp

} # End Model Loop
# --------------------------------


# --------------------------------------------------------------------------
# Run & Process gamms -- No Composition Effect
# --------------------------------------------------------------------------
cores.use <- min(12, length(paleon.models))
# cores.use <- length(paleon.models)

models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=cores.use, response=response, k=k, predictors.all=c(predictors.all, "Evergreen", "Grass"), comp.effects=F)
# -------------------------------------------------

# -------------------------------------------------
# Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_models_NPP_Site.Rdata"))
# -------------------------------------------------


# -------------------------------------------------
# Diagnostic Graphs
# -------------------------------------------------
m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

pdf(file.path(fig.dir, "GAMM_Models_ModelFit_NPP_Site.pdf"))
print(
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Model, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
)
# print(	
# ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	# geom_line(data= mod.out$data[,], aes(x=Year, y=Y), alpha=0.5) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	# scale_x_continuous(limits=c(1850,2010)) +
	# # scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	# labs(title=paste("Composition Effects", response, sep=" - "), x="Year", y=response)
# )
dev.off()


pdf(file.path(fig.dir, "GAMM_Models_DriverSensitivity_NPP_Site.pdf"))
for(e in unique(mod.out$ci.terms$Effect)[1:3]){
print(
ggplot(data=mod.out$ci.terms[mod.out$ci.terms$Effect==e,]) + facet_wrap( ~ Site, scales="fixed") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0(e, " Sensitivity (Non-Relativized)"), y=paste0("RW Contribution")) # +
)
}
dev.off()
# -------------------------------------------------
# --------------------------------------------------------------------------

# Clear the memory!
rm(mod.out, models.base)

# ------------------------------------------------------------------------------------
