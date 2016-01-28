# ----------------------------------------
# Objective: Create "baseline" sensitivity curves for each model/data stream 
#            that are based on all sites, all time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 28 July 2015
# ----------------------------------------
#
# -------------------------
# Workflow
# -------------------------
# 1. Set up Data
#    a. Ecosystem model output
#    b. Tree Ring NPP products
#    c. Raw Tree Ring widths
# 2. Run the gamms (with site intercept)
# 3. Bind Models into single list
# 4. Diagnostic Graphs
# -------------------------
# ----------------------------------------

rm(list=ls())

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
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")
resolutions <- "t.001"
k=4
# ----------------------------------------

# ----------------------------------------
# Set Directories & file paths
# ----------------------------------------
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/gamms"
fig.base="Figures/gamms"

# Source the gamm file
source('R/0_calculate.sensitivity_TPC.R', chdir = TRUE)

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- file.path(fig.base, "Sensitivity_Baseline")
dat.dir <- file.path(dat.base, "Sensitivity_Baseline")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# -------------------------------------------------------------------------------
# 1. Set up Data 
# -------------------------------------------------------------------------------
{
paleon.models <- list()
# ----------------------------------------
# 1.a. Load & set up Ecosystem Model Output first
# ----------------------------------------
{
# Define what our response variable will be
response <- "NPP"

# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData.Rdata"))

# Get rid of LINKAGES because it's weird & hasn't been updated
ecosys <- ecosys[!ecosys$Model=="linkages",]
summary(ecosys)

for(m in unique(ecosys$Model)){

	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m, " ------"))

	# Taking the subsets of data we want in a single gam
	dat.subsets <- ecosys$Resolution == "t.001" & 
		           ecosys$Model      == m

	# What will our spatio-temporal explanatory factor ("Time") be?
	if(!is.na(mean(ecosys[dat.subsets,"AGB"]))) time.mod="AGB" else time.mod="LAI"

	data.temp                  <- ecosys[dat.subsets, c("Model", "Model.Order", "Site", "Year")]
	data.temp$PlotID           <- ecosys[dat.subsets,"Site" ]
	data.temp$TreeID           <- as.factor(NA)
	data.temp$Y                <- ecosys[dat.subsets,response]
	data.temp$Time             <- ecosys[dat.subsets,time.mod]
	data.temp[,predictors.all] <- ecosys[dat.subsets, paste0(predictors.all, predictor.suffix)]
	data.temp$Resolution       <- ecosys[dat.subsets,"Resolution"]

	# Getting rid of NAs in predictors
	data.temp <- data.temp[complete.cases(data.temp[,c(predictors.all, "Y", "Time")]),]

	# Copy the response variable & some other things for the model
	paleon.models[[paste(m)]] <- data.temp

} # End Model Loop
} # End Model setup
# ----------------------------------------

# ----------------------------------------
# 1.b. Load & set up tree ring NPP
# ----------------------------------------
{
# Define what our response & time variables will be
response <- "ABI.area"
time.mod <- "AB.area"

# Load Tree ring NPP data
spp.npp <- read.csv(file.path("Data", "TreeRing_NPP_PlotSpecies.csv"))
summary(spp.npp)

# aggregate to total plot NPP (ABI.area)
plot.npp <- aggregate(spp.npp[,c(response, time.mod)], by=spp.npp[,c("Site", "Site2", "PlotID", "Year")], FUN=sum)
plot.npp[,c(paste0(predictors.all, predictor.suffix))] <- aggregate(spp.npp[,c(paste0(predictors.all, predictor.suffix))], by=spp.npp[,c("Site", "Site2", "PlotID", "Year")], FUN=mean)[,c(paste0(predictors.all, predictor.suffix))]
summary(plot.npp)

# Subset a period where we're not worried about 
plot.npp <- plot.npp[complete.cases(plot.npp) & plot.npp$Year>=(2010-30),]
summary(plot.npp)

# Add the data to paleon.models
paleon.models[["TreeRingNPP"]]             <- plot.npp[,c("Site", "PlotID", "Year")]
paleon.models$TreeRingNPP$Model            <- as.factor("TreeRingNPP")
paleon.models$TreeRingNPP$Model.Order      <- as.factor("Tree Ring NPP")
paleon.models$TreeRingNPP$TreeID           <- as.factor(NA)
paleon.models$TreeRingNPP$Y                <- plot.npp[,response]
paleon.models$TreeRingNPP$Time             <- plot.npp[,time.mod]
paleon.models$TreeRingNPP[,predictors.all] <- plot.npp[, paste0(predictors.all, predictor.suffix)]
paleon.models$TreeRingNPP$Resolution       <- as.factor("t.001")

# Make sure everything is complete cases
paleon.models$TreeRingNPP <- paleon.models$TreeRingNPP[complete.cases(paleon.models$TreeRingNPP[,c(predictors.all, "Y", "Time")]),]

# Order everything the same way to make life easier
paleon.models$TreeRingNPP <- paleon.models$TreeRingNPP[,names(paleon.models[[1]])]
summary(paleon.models$TreeRingNPP)
} # End Tree Ring NPP setup
# ----------------------------------------

# ----------------------------------------
# 1.c. Load & set up raw ring widths
# ----------------------------------------
{
response <- "RW"
time.mod <- "DBH"

tree.rings <- read.csv("Data/TreeRing_RingWidths.csv")
summary(tree.rings)

# subset only complete cases where we have met data & ring width measurements
tree.rings <- tree.rings[complete.cases(tree.rings[,c(response, paste0(predictors.all, predictor.suffix))]) & tree.rings$Resolution=="t.001",]
summary(tree.rings)

# -----------
# Adding in a few other factors that could be good predictors
# -----------
plot.age <- aggregate(tree.rings[,c("Age", "DBH")], by=tree.rings[,c("PlotID", "Year")], FUN=mean, na.rm=T)
plot.age$BA.sum <- aggregate(tree.rings[,c("DBH")], by=tree.rings[,c("PlotID", "Year")], FUN=function(x){sum(pi*((x/2)^2), na.rm=T)})[,3]
plot.age$BA.mean <- aggregate(tree.rings[,c("DBH")], by=tree.rings[,c("PlotID", "Year")], FUN=function(x){mean(pi*((x/2)^2), na.rm=T)})[,3]
names(plot.age)[3:ncol(plot.age)] <- paste0("plot.", names(plot.age)[3:ncol(plot.age)])
summary(plot.age)

tree.rings <- merge(tree.rings, plot.age, all.x=T, all.y=T)

tree.rings$BA <- pi*(tree.rings$DBH/2)^2
summary(tree.rings)
# -----------

# Add the data to paleon.models
paleon.models[["TreeRingRW"]]             <- tree.rings[,c("Site", "PlotID", "TreeID", "Year")]
paleon.models$TreeRingRW$Model            <- as.factor("TreeRingRW")
paleon.models$TreeRingRW$Model.Order      <- as.factor("Tree Ring RW")
paleon.models$TreeRingRW[,predictors.all] <- tree.rings[,paste0(predictors.all, predictor.suffix)]
paleon.models$TreeRingRW$Y                <- tree.rings[,response]
paleon.models$TreeRingRW$Time             <- tree.rings[,time.mod]
paleon.models$TreeRingRW$Resolution       <- tree.rings[,"Resolution"]

# Make sure everything is complete cases
paleon.models$TreeRingRW <- paleon.models$TreeRingRW[complete.cases(paleon.models$TreeRingRW[,c(predictors.all, "Y", "Time")]),]

# Order everything the same way to make life easier
paleon.models$TreeRingRW <- paleon.models$TreeRingRW[,names(paleon.models[[1]])]
summary(paleon.models$TreeRingRW)
} # End Ring Width setups
# ----------------------------------------
}
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# 2. Run the gamms 
# -------------------------------------------------------------------------------
cores.use <- min(12, length(paleon.models))
# cores.use <- length(paleon.models)

models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=cores.use, k=k, predictors.all=predictors.all, PFT=F)
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# 3. Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------------------------------------
{
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], "baseline", sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], "baseline", sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_baseline.Rdata"))
}
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# 4. Diagnostic Graphs
# -------------------------------------------------------------------------------
{
m.order <- unique(mod.out$data$Model.Order)
col.model <- c(paste(model.colors[model.colors$Model.Order %in% m.order,"color"]), "black", "gray30")

pdf(file.path(fig.dir, "GAMM_ModelFit_Baseline.pdf"))
{
print(	
ggplot(data=mod.out$ci.response[!substr(mod.out$ci.response$Model, 1, 8)=="TreeRing",]) + facet_grid(PlotID~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[!substr(mod.out$data$Model, 1, 8)=="TreeRing",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Baseline"), x="Year", y="NPP")
)
print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Model=="TreeRingNPP",]) + facet_wrap(~ PlotID, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Model=="TreeRingNPP",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Baseline"), x="Year", y="NPP")
)
print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Model=="TreeRingRW" & mod.out$ci.response$PlotID=="ME029",]) + facet_wrap(~ TreeID, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Model=="TreeRingRW" & mod.out$data$PlotID=="ME029",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Baseline"), x="Year", y="RW")
)
}
dev.off()


pdf(file.path(fig.dir, "GAMM_DriverSensitivity_Baseline.pdf"))
print(
ggplot(data=mod.out$ci.terms[!mod.out$ci.terms$Effect %in% c("Time", "Year"),]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (not Relativized)"), y=paste0("NPP Contribution")) # +
)
dev.off()
}
## -------------------------------------------------------------------------------


