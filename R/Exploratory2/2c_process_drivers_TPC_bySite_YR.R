# ----------------------------------------
# Temporal Scaling Analyses -- Create baseline growing season model
# Non-linear driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 10 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: What is the relative control of different drivers within each model?
# Rationale: Not all models use all inputs, and many drivers are correlated, so we need to 
#            see if the temperature pattern is really a radiation pattern, etc. 
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# (Fit GAMM per site per m.name)
# 1) Temporal Grain (Resolution)
#    -- Fit GAMM over constant wind.gsow with different degrees of smoothing (1 yr - 250 yr)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# 1) Space-Time Comparison
#    -- Hypothesis: Driver responses across sites within a m.name converge at coarser temporal grains 
#       and larger extents because the models have time to adjust and seek equilibrium.
#
#    -- Analysis: Use the posterior CIs for each smoothing term to see if the driver curves for sites
#                 within a m.name are statstically different at different sites at different scales (or 
#                 alternatively if there are statistical differences in [parts of] the curves between
#                 scales)
#
#
# 2) Multi-Model Driver Comparison
#    -- Hypothesis: Because models were built & parameterized to perform well in the modern era, there 
#       will be greater agreement of the direction & primary driver of change in the more recent time 
#       periods than over the full extent of the PalEON runs.
#    -- Hypothesis about temporal grain?
#
#    -- Analysis: For each given scale, compare the m.name response curves for each driver.  Determine 
#                 which drivers are most similar/variable among models and at what scales?  Are there
#                 particular ranges of each driver response where models responses are most similar/different?
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(parallel)
library(ncdf4)
library(lme4)
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
fig.dir <- file.path(fig.base, "TPC_YR_byResolution_Site")
dat.dir <- file.path(dat.base, "TPC_YR_byResolution_Site")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData_Raw.Rdata"))
summary(ecosys)
model.colors

source('R/0_calculate.sensitivity_TPC.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)

# Read in model color scheme
model.colors
# ----------------------------------------


# -------------------------------------------------
# Settings for the rest of this script
# -------------------------------------------------
# Get rid of CLM-BGC because its actual drivers are messed up
ecosys <- ecosys[!substr(ecosys$Model, 1, 3)=="clm" & !ecosys$Model=="linkages",]

# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)
#resolutions <- c("t.001", "t.010", "t.050", "t.100")
resolutions <- c("t.001", "t.010", "t.050") # Note: Big models can't handle t.100 at the site level because there aren't enough data points
extents <- data.frame(Start=c(850, 1850, 1990), End=c(2010, 2010, 2010)) 
# response.all <- c("NPP", "NEE", "AGB.diff")
response.all <- c("NPP")
# predictors.all <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".yr")
k=4
e=1	
# -------------------------------------------------

# -------------------------------------------------
# Set up the appropriate data for each model into a list
# -------------------------------------------------
for(y in 1:length(response.all)){
	response <- response.all[y]
	print("-------------------------------------")
	print("-------------------------------------")
	print(paste0("------ Processing Var: ", response, " ------"))
	# M1 <- 1:length(model.name)
	# M2 <- which(!model.name=="jules.stat")
	#if(# response=="AGB.diff") models <- which(!model.name=="jules.stat") else models <- 1:length(model.name)

for(m in 1:length(model.name)){
	paleon.models <- list()
	m.name  <- model.name[m]
	m.order <- model.order[m]

	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))

	# Note: Here we're renaming things that had the suffix to just be generalized tair, etc 
	dat.mod <- ecosys[ecosys$Model==m.name, c("Model", "Updated", "Model.Order", "Site", "Year", response, paste0(predictors.all, predictor.suffix))]
	names(dat.mod)[(ncol(dat.mod)-length(predictors.all)+1):ncol(dat.mod)] <- predictors.all
	
	if(!max(dat.mod[,response], na.rm=T)>0) next # If a variable is missing, just skip over this model for now

for(r in 1:length(resolutions)){ # Resolution loop

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.gs & temperature distributions
	run.end <- ifelse(substr(m.name,1,3)=="jul", as.numeric(extents[e,2])-1, as.numeric(extents[e,2])) # Note: Jules missing 2010, so 
	run.start <- as.numeric(extents[e,1])
	inc <- round(as.numeric(substr(resolutions[r],3,5)),0) # making sure we're always dealign with whole numbers
	yrs <- seq(from=run.end-round(inc/2,0), to=run.start+round(inc/2,0), by=-1)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "Year")]

	# Making a note of the extent & resolution
	ext <- as.factor(paste(extents[e,1], extents[e,2], sep="-"))
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions[r])

	# Making place-holders for the response & predictors so the loop works correctly
	data.temp[,c(response, predictors.all)] <- NA

	# Calculating the mean for each wind.gsow for resolution
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

	data.temp$Y <- data.temp[,response]

	for(s in 1:length(sites)){
		data.temp2 <- data.temp[data.temp$Site==sites[s],]

		paleon.models[[paste0(sites[s], ".", resolutions[r])]] <- data.temp2
	}


} # End Resolution Loop
# --------------------------------

# -------------------------------------------------
# Run the gamms
# -------------------------------------------------
# cores.use <- min(12, length(paleon.models))
cores.use <- length(paleon.models)

models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=cores.use, response=response, k=k, predictors.all=predictors.all, site.effects=F)
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
		mod.out[[paste("gamm", ext, substr(names(models.base)[i],1,3), substr(names(models.base)[i],7,nchar(paste(names(models.base)))), sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", ext, substr(names(models.base)[i],1,3), substr(names(models.base)[i],7,nchar(paste(names(models.base)))), sep=".")]] <- models.base[[i]]$gamm
	}
}


# -------------------------------------------------
# Bind Resolutions together to make them easier to work with
# -------------------------------------------------

m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

save(mod.out, file=file.path(dat.dir, paste0("gamm_TPC_YR_", m.name,"_", response, ".Rdata")))

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_TPC_YR_Site_", m.order, "_", response, ".pdf")))
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


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_TPC_YR_Site_", m.order, "_", response, ".pdf")))
for(r in unique(mod.out$ci.terms$Resolution)){
print(
ggplot(data=mod.out$ci.terms[mod.out$ci.terms$Resolution==r,]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Site), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_color_manual(values=c("red2", "blue", "green3")) +
	# scale_fill_manual(values=c("red2", "blue", "green3")) +
	labs(title=paste0("Driver Effects: ", r, ": ", m.order), y=paste0(r, " Effect Size")) # +
)
}
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_Time_Site_", m.order, "_", response, "_0850-2010.pdf")))
print(plot.weights.time(df=mod.out$weights, xmin=851, xmax=2010, breaks=c(1000, 1250, 1500, 1750, 2000), plot.labs=labs(x="Year", title=paste0(m.order, " Driver Weights through Time: 0850 - 2010"))) )
print(plot.weights.time(df=mod.out$weights, xmin=1800, xmax=1900, breaks=c(1825, 1850, 1875), plot.labs=labs(x="Year", title=paste0(m.order, " Driver Weights through Time: 1800 - 1900"))) )
print( plot.weights.time(df=mod.out$weights, xmin=1900, xmax=2010, breaks=c(1925, 1950, 1975), plot.labs=labs(x="Year", title=paste0(m.order, " Driver Weights through Time: 1900 - 2010"))) )
dev.off()



# -------------------------------------------------
} # End by Model Loop
} # End Response loop
