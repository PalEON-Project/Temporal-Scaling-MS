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
# 2) Model NPP
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# Notes:
#  -- these relationships are truncated in time & space relative to the full MIP, so the other 
#     runs are needed to put the results into context
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(parallel)
library(mgcv)
# library(ncdf4)
# library(lme4)
# library(R2jags)
library(ggplot2); library(grid)
library(car)
# library(zoo)
# library(mvtnorm)
# library(MCMCpack)
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
dat.base="Data/gamms"
fig.base="Figures/gamms"

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- file.path(fig.base, "AllDrivers_GS_TreeRings")
dat.dir <- file.path(dat.base, "AllDrivers_GS_TreeRings")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
source('R/0_gamm.calculate2.R', chdir = TRUE)

# Read in Tree Ring data
dat.tr <- read.csv("Data/TreeRings_processed/TreeRing_Chronologies.csv")
summary(dat.tr)

# Subsetting tree rings to only be post-1901 because then we're running with CRUNCEP data with
#   observed variability rather than a model representation
dat.tr <- dat.tr[dat.tr$Year>=1901,]
dat.tr <- dat.tr[complete.cases(dat.tr),]
summary(dat.tr)

# Read ecosys file
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData_Raw.Rdata"))
summary(ecosys)
model.colors

# Read in model color scheme
model.colors

# Subset Ecosys to match the spatial & temporal extent of the tree ring data
for(s in unique(dat.tr$Site)){
	yr.min <- min(dat.tr[dat.tr$Site==s,"Year"])
	yr.max <- max(dat.tr[dat.tr$Site==s,"Year"])
	
	if(s==unique(dat.tr$Site)[1]){
		ecosys2 <- ecosys[ecosys$Site==s & ecosys$Year>=yr.min & ecosys$Year<=yr.max, ]
	} else {
		ecosys2 <- rbind(ecosys2, ecosys[ecosys$Site==s & ecosys$Year>=yr.min & ecosys$Year<=yr.max, ])
	}
}
summary(ecosys2)
# ----------------------------------------


# -------------------------------------------------
# Subsetting the datasets for Analysis
# -------------------------------------------------
resolutions <- c("t.001", "t.010")
response <- "NPP"
predictors.all <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
predictor.suffix <- c(".gs")
k=3
# r=1	
# -------------------------------------------------

# -------------------------------------------------
# Run the Tree ring model
# NOTE: Running it separately here because it requires a slightly different model
#       structure than the other models
# -------------------------------------------------
source('R/0_process.gamm.R', chdir = TRUE)
source('R/0_gamm_Plots.R', chdir = TRUE)

predictors <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")


# Making up a model & Model.Order
tr.type <- c("All", "ITRDB", "NPP")
NPP.prefix <- c("HO", "LF", "TP")

for(t in tr.type){

fig.dir <- file.path(fig.base, "AllDrivers_GS_TreeRings", paste0("TreeRings_",t))
dat.dir <- file.path(dat.base, "AllDrivers_GS_TreeRings", paste0("TreeRings_",t))

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)


if(t == "All") 	 rows.use <- 1:nrow(dat.tr)
if(t == "ITRDB") rows.use <- which(!(substr(dat.tr$PlotID,1,2) %in% NPP.prefix))
if(t == "NPP")   rows.use <- which(substr(dat.tr$PlotID,1,2) %in% NPP.prefix)


m.name="tree.rings"
m.order="Tree_Rings"

dat.mod <- dat.tr[rows.use,c("Site", "PlotID", "Year", paste0(predictors, ".gs"))]
names(dat.mod)[(ncol(dat.mod)-length(predictors)+1):ncol(dat.mod)] <- predictors
dat.mod$NPP <- dat.tr[rows.use, "RWI"]
dat.mod <- dat.mod[complete.cases(dat.mod),]
dat.mod$Model       <- as.factor(m.name)
dat.mod$Model.Order <- as.factor(m.order)
dat.mod$Updated      <- as.factor("Yes")
# summary(dat.mod)

dat.tr2 <- dat.mod
# Get rid of CLM-BGC because its actual drivers are messed up
ecosys3 <- ecosys2[!ecosys2$Model=="clm.bgc" & (ecosys2$Site %in% unique(dat.tr2$Site)) & ecosys2$Year>=min(dat.tr2$Year) & ecosys2$Year<=max(dat.mod$Year),]


for(r in 1:length(resolutions)){ # Resolution loop

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.yr & temperature distributions
	run.end <- max(dat.mod$Year)
	run.start <- min(dat.mod$Year)
	inc <- round(as.numeric(substr(resolutions[r],3,5)),0) # making sure we're always dealign with whole numbers
	yrs <- seq(from=run.end-round(inc/2,0), to=run.start+round(inc/2,0), by=-inc)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "PlotID", "Year")]

	# Making a note of the extent & resolution
	ext <- as.factor(paste(run.start, run.end, sep="-"))
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions[r])

	# Making place-holders for the response & predictors so the loop works correctly
	data.temp[,c(response, predictors.all)] <- NA

	plotIDs <- unique(data.temp$PlotID)
	# Calculating the mean for each wind.yrow for resolution
	# Note: because we're now only analyzing single points rathern than the full running mean, 
	#    we're now making the year in the middle of the resolution
	if(inc==1){ # if we're working at coarser than annual scale, we need to find the mean for each bin
		data.temp[,c(response, predictors.all)] <- dat.mod[dat.mod$Year %in% yrs,c(response, predictors.all)]
	} else {
		for(s in plotIDs){
			for(y in yrs){
				data.temp[data.temp$PlotID==s & data.temp$Year==y,c(response, predictors.all)] <- apply(dat.mod[dat.mod$PlotID==s & dat.mod$Year>=round(y-inc/2, 0) & dat.mod$Year<=round(y+inc/2, 0),c(response, predictors.all)], 2, FUN=mean, na.rm=F)
			}
		}
	}

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,response]),]

	model.name <- unique(data.temp$Model)
	ext.name   <- unique(data.temp$Extent)
	ext.index  <- regexpr("-", ext.name)[1] # Find the index so we can split apart extent into 2 numbers
	extent     <- c(as.numeric(substr(ext.name,1,ext.index-1)), as.numeric(substr(ext.name, ext.index+1,nchar(paste(ext.name)))))
	resolution <- unique(data.temp$Resolution)

	gam1 <- gam(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, data=data.temp, correlation=corARMA(form=~Year|PlotID, p=1))

	print(summary(gam1))	

	# get rid of values for predictors not used in the models for clarity later on
	data.temp[,predictors.all[!(predictors.all %in% predictors)]] <- NA

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1, newdata=data.temp)
	# ----------------------------------------

	# ----------------------------------------
	# Run all of the post-processing (calculate CIs, etc)
	# ----------------------------------------
	mod.temp <- process.gamm(gamm.model=gam1, data=data.temp, model.name=model.name, extent=extent, resolution=resolution, response=response, vars=predictors, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	# ----------------------------------------

	# putting stuff together
	if(r==1) {
		mod.out <- list()
		mod.out$data         <- mod.temp$data
		mod.out$weights      <- mod.temp$weights
		mod.out$ci.response  <- mod.temp$ci.response
		mod.out$sim.response <- mod.temp$sim.response
		mod.out$ci.terms     <- mod.temp$ci.terms
		mod.out$sim.terms    <- mod.temp$sim.terms
		mod.out[[paste("gamm", ext, substr(resolution,3,5), sep=".")]] <- mod.temp$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         mod.temp$data)
		mod.out$weights      <- rbind(mod.out$weights,      mod.temp$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  mod.temp$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, mod.temp$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     mod.temp$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    mod.temp$sim.terms)
		mod.out[[paste("gamm", ext, substr(resolution,3,5), sep=".")]] <- mod.temp$gamm
	}
} # End Resolution Loop



m.order <- unique(mod.out$data$Model.Order)
col.model <- "black"

save(mod.out, file=file.path(dat.dir, paste0("gamm_AllDrivers_Yr_", m.name, "_", response, ".Rdata")))

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_GS_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP, color=PlotID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP, color=PlotID), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	scale_x_continuous(limits=c(1950,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_AllDrivers_GS_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.terms[,]) + facet_wrap(~ Effect, scales="free_x") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Resolution), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Resolution), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_color_manual(values=c("red2", "blue", "green3")) +
	# scale_fill_manual(values=c("red2", "blue", "green3")) +
	labs(title=paste0("Driver Effects: ",m.order), y="Effect Size") # +
)
dev.off()
# } # End Tree-ring data type
# -------------------------------------------------


# -------------------------------------------------
# Run the Models!
# -------------------------------------------------
# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys3$Site)
model.name  <- unique(ecosys3$Model)
model.order <- unique(ecosys3$Model.Order)
response = "NPP"
for(m in 1:length(model.name)){
	paleon.models <- list()
	m.name  <- model.name[m]
	m.order <- model.order[m]

	print("-------------------------------------")
	print("-------------------------------------")
	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))

	# Note: Here we're renaming things that had the suffix to just be generalized tair, etc 
	dat.mod <- ecosys3[ecosys3$Model==m.name, c("Model", "Updated", "Model.Order", "Site", "Year", response, paste0(predictors.all, predictor.suffix))]
	names(dat.mod)[(ncol(dat.mod)-length(predictors.all)+1):ncol(dat.mod)] <- predictors.all

for(r in 1:length(resolutions)){ # Resolution loop

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.yr & temperature distributions
	run.end <- max(dat.mod$Year)
	run.start <- min(dat.mod$Year)
	inc <- round(as.numeric(substr(resolutions[r],3,5)),0) # making sure we're always dealign with whole numbers
	yrs <- seq(from=run.end-round(inc/2,0), to=run.start+round(inc/2,0), by=-inc)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "Year")]

	# Making a note of the extent & resolution
	ext <- as.factor(paste(run.start, run.end, sep="-"))
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions[r])


	# Making place-holders for the response & predictors so the loop works correctly
	data.temp[,c(response, predictors.all)] <- NA

	# Calculating the mean for each wind.yrow for resolution
	# Note: because we're now only analyzing single points rathern than the full running mean, 
	#    we're now making the year in the middle of the resolution
	if(inc==1){ # if we're working at coarser than annual scale, we need to find the mean for each bin
		data.temp[,c(response, predictors.all)] <- dat.mod[dat.mod$Year %in% yrs,c(response, predictors.all)]
	} else {
		for(s in sites){
			for(y in yrs){
				data.temp[data.temp$Site==s & data.temp$Year==y,c(response, predictors.all)] <- apply(dat.mod[dat.mod$Site==s & dat.mod$Year>=round(y-inc/2, 0) & dat.mod$Year<=round(y+inc/2, 0),c(response, predictors.all)], 2, FUN=mean)
			}
		}
	}

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,response]),]

	# Make a new variable called Y with the response variable so it can be generalized
	data.temp$Y <- data.temp[,response]

	paleon.models[[paste(resolutions[r])]] <- data.temp
} # End Resolution Loop
# --------------------------------

# -------------------------------------------------
# Run the gamms
# -------------------------------------------------
models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=length(paleon.models), response=response, k=k, predictors.all=predictors.all, site.effects=T)
# -------------------------------------------------


# -------------------------------------------------
# Bind Resolutions together to make them easier to work with
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
		mod.out[[paste("gamm", ext, substr(names(models.base)[i],3,nchar(paste(names(models.base)))), sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", ext, substr(names(models.base)[i],3,nchar(paste(names(models.base)))), sep=".")]] <- models.base[[i]]$gamm
	}
}

m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

save(mod.out, file=file.path(dat.dir, paste0("gamm_AllDrivers_Yr_", m.name, "_", response, ".Rdata")))

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_GS_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
print(	
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	scale_x_continuous(limits=c(1850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_AllDrivers_GS_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.terms[,]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Resolution), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Resolution), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_color_manual(values=c("red2", "blue", "green3")) +
	# scale_fill_manual(values=c("red2", "blue", "green3")) +
	labs(title=paste0("Driver Effects: ",m.order), y="Effect Size") # +
)
dev.off()
# -------------------------------------------------
} # End by Model Loop
} # End Tree Ring Subset