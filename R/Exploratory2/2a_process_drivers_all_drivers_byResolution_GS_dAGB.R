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
fig.dir <- file.path(fig.base, "Big4_AGB_GS_lag1_byResolution")
dat.dir <- file.path(dat.base, "Big4_AGB_GS_lag1_byResolution")

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

source('R/0_calculate.sensitivity_4drivers_lag1.R', chdir = TRUE)

# Read in model color scheme
model.colors
# ----------------------------------------


# -------------------------------------------------
# Settings for the rest of this script
# -------------------------------------------------
# Get rid of CLM-BGC because its actual drivers are messed up
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)
resolutions <- c("t.001", "t.010", "t.050", "t.100")
response.all <- c("NPP", "AGB.diff", "NEE")
# predictors.all <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
predictors.all <- c("tair", "precipf", "swdown", "CO2")
predictor.suffix <- c(".gs")
k=4
y.lag=5
# r=1	
# -------------------------------------------------

# -------------------------------------------------
# Set up the appropriate data for each model into a list
# -------------------------------------------------
# for(y in 1:length(response.all)){
y=2
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

	for(s in sites){
		for(y in (min(dat.mod$Year)+y.lag):max(dat.mod$Year)){
			dat.mod[dat.mod$Site==s & dat.mod$Year==y,"Y.lag"] <- mean(dat.mod[dat.mod$Site==s & dat.mod$Year>=(y-y.lag) & dat.mod$Year<y,response], na.rm=T)
		}
	 }

for(r in 1:length(resolutions)){ # Loop through different resolutions

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.yr & temperature distributions
	run.end <- ifelse(substr(m.name,1,3)=="jul", max(ecosys$Year)-1, max(ecosys$Year)) # Note: Jules missing 2010, so 
	run.start <- 850
	inc <- round(as.numeric(substr(resolutions[r],3,5)),0) # making sure we're always dealing with whole numbers
	yrs <- seq(from=run.end-round(inc/2,0), to=run.start+round(inc/2,0), by=-inc)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "Year")]

	# Making a note of the extent & resolution
	ext <- as.factor("850-2010")
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions[r])

	# Making place-holders for the response & predictors so the loop works correctly
	data.temp[,c(response, predictors.all, "Y.lag")] <- NA

	# Calculating the mean for each wind.yrow for resolution
	# Note: because we're now only analyzing single points rathern than the full running mean, 
	#    we're now making the year in the middle of the resolution
	if(inc==1){ # if we're working at coarser than annual scale, we need to find the mean for each bin
		data.temp[,c(response, predictors.all, "Y.lag")] <- dat.mod[,c(response, predictors.all, "Y.lag")]
	} else {
		for(s in sites){
			for(y in yrs){
				data.temp[data.temp$Site==s & data.temp$Year==y,c(response, predictors.all, "Y.lag")] <- apply(dat.mod[dat.mod$Site==s & dat.mod$Year>=round(y-inc/2, 0) & dat.mod$Year<=round(y+inc/2, 0),c(response, predictors.all, "Y.lag")], 2, FUN=mean)
			}
		}
	}


	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp   <- data.temp[complete.cases(data.temp[,c(response, "Y.lag")]),]

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

save(mod.out, file=file.path(dat.dir, paste0("gamm_AllDrivers_Yr_", m.name, "_", response, "_Lag5.Rdata")))

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_GS_", m.order, "_", response, "_Lag5.pdf")))
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
	scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
dev.off()

# print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Resolution=="t.010",]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Resolution=="t.010",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	# scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, "Decadal", sep=" - "), x="Year", y=response)
# )

pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_AllDrivers_GS_", m.order, "_", response, "_Lag10.pdf")))
print(
ggplot(data=mod.out$ci.terms[!mod.out$ci.terms$Effect=="Y.lag",]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		
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
# } # End Response Loop