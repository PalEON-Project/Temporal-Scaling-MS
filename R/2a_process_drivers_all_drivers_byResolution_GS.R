# ----------------------------------------
# Temporal Scaling Analyses -- Annual Met Drivers
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
fig.dir <- file.path(fig.base, "AllDrivers_GS_byResolution")
dat.dir <- file.path(dat.base, "AllDrivers_GS_byResolution")

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

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2.gs
# 	predict.gamm.model.site.R = function to run a single 1 site - m.name combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal m.name responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal m.name fitting with random site effects very slow
source('R/0_process.gamm.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)


# Read in model color scheme
model.colors
# ----------------------------------------


# # -----------------------
# # Some exploratory Graphing
# # -----------------------
# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # ggtitle("Annual & Centennial AGB") +
  # theme_bw()

# # quantile(ecosys$AGB.diff, c(0.025, 0.33, 0.67, 0.975), na.rm=T)
# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # # geom_line(aes(x=Year, y=AGB.diff, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=AGB.diff.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # scale_y_continuous(limits=quantile(ecosys$AGB.diff.100, c(0.025, 0.975), na.rm=T)) +
  # ggtitle("Annual & Centennial dAGB") +
  # theme_bw()

# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=NPP.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # # scale_y_continuous(limits=quantile(ecosys$NPP.100, c(0.025, 0.975), na.rm=T)) +
  # ggtitle("Annual & Centennial NPP") +
  # theme_bw()
# # -----------------------

# ----------------------------------------



# ------------------------
# MegaLoop -- Looping through all models by Variable, by Extent
# ------------------------
# Get rid of CLM-BGC because its actual drivers are messed up
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

# -----------------
# Matrix of Models and Drivers
# -----------------
# Var    ED2  ED2-LU  CLM-BGC  CLM-CN  LPJ-WSL  LPJ-GUESS  JULES-STAT  JULES-TRIFFID  SIBCASA  LINKAGES
# tair.gs    X     X        X       X        X         X          X             X           X        X
# precipf.gs X     X        X       X        X         X          X             X           X        X
# swdown.gs  X     X        X       X        X         X          X             X           X
# lwdown.gs  X     X                         X                    X             X           X
# wind.gs    X     X        X       X                             X             X           X
# psurf.gs   X     X        X       X                             X             X           X
# qair.gs    X     X        X       X                             X             X           X
# CO2.gs     X     X        X       X        X         X          X             X           X
# Ndep                   ?       ?

# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)
resolutions <- c("t.001", "t.010", "t.050", "t.100")
response <- "NPP"
predictors.all <- c("tair.gs", "precipf.gs", "swdown.gs", "lwdown.gs", "psurf.gs", "qair.gs", "wind.gs", "CO2.gs")
k=4
	
for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]
	print("-------------------------------------")
	print("-------------------------------------")
	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))

	# if(m.name=="lpj.guess") resolutions = resolutions.all[1:3] else resolutions = resolutions.all
	dat.mod <- ecosys[ecosys$Model==m.name, c("Model", "Updated", "Model.Order", "Site", "Year", response, predictors.all)]
for(r in 1:length(resolutions)){ # Resolution loop
	print(       "-------------------------------------")
	print(paste0("------ Processing Resolution: ", resolutions[r], " ------"))

    # -----------
	# NOTE: To really get the resolution analysis right, we don't want running averages; we want
	#    single values that are the mean of a given window
    # -----------

	# Figure out which years to take: 
	# Note: working backwards to help make sure we get modern end of the CO2.gs & temperature distributions
	run.end <- ifelse(substr(m.name,1,3)=="jul", max(ecosys$Year)-1, max(ecosys$Year)) # Note: Jules missing 2010, so 
	run.start <- 850
	inc <- round(as.numeric(substr(resolutions[r],3,5)),0) # making sure we're always dealign with whole numbers
	yrs <- seq(from=run.end-round(inc/2,0), to=run.start+round(inc/2,0), by=-inc)

	data.temp <- dat.mod[(dat.mod$Year %in% yrs), c("Model", "Updated", "Model.Order", "Site", "Year")]

	# Making a note of the extent & resolution
	ext <- as.factor("850-2010")
	data.temp$Extent <- as.factor(ext)
	data.temp$Resolution <- as.factor(resolutions[r])

	# Making place-holders for the response & predictors so the loop works correctly
	data.temp[,c(response, predictors.all)] <- NA

	# Calculating the mean for each wind.gsow for resolution
	# Note: because we're now only analyzing single points rathern than the full running mean, 
	#    we're now making the year in the middle of the resolution
	if(inc>1){ # if we're working at coarser than annual scale, we need to find the mean for each bin
		for(s in sites){
			for(y in yrs){
				data.temp[data.temp$Site==s & data.temp$Year==y,c(response, predictors.all)] <- apply(dat.mod[dat.mod$Site==s & dat.mod$Year>=round(y-inc/2, 0) & dat.mod$Year<=round(y+inc/2, 0),c(response, predictors.all)], 2, FUN=mean)
			}
		}
	} else {
		data.temp[,c(response, predictors.all)] <- dat.mod[,c(response, predictors.all)]
	}

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,response]),]

    # -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	# This is different from model.site.gam in that it has an added random site slope.
	#	This random effect lets us gauge the overall model.name response to our fixed effects 
	#   regardless of the site.  
	#   Pros: Generalized and helps characterize the basic model responses
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)

	# Select which set of predictors based on which model it is
	# Each of the models is having different stability issues
	if(substr(m.name,1,2)=="ed"){
		predictors <- c("tair.gs", "precipf.gs", "swdown.gs", "lwdown.gs", "psurf.gs", "qair.gs", "wind.gs", "CO2.gs")
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(lwdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
	if(substr(m.name,1,3)=="clm") {
		predictors <- c("tair.gs", "precipf.gs", "swdown.gs", "psurf.gs", "qair.gs", "wind.gs", "CO2.gs")
		# if(m.name=="clm.bgc"){
			gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(niterEM=0, sing.tol=1e-20, opt="optim"))
		# } else {
			# gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))

		# }
	}
	if(substr(m.name,1,3)=="lpj") {
		predictors <- c("tair.gs", "precipf.gs", "swdown.gs", "CO2.gs")
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(niterEM=0, sing.tol=1e-20, method="optim")) 
	# , control=list(niterEM=0, sing.tol=1e-20, method="optim")
	}
	if(substr(m.name,1,3)=="jul") {
		predictors <- c("tair.gs", "precipf.gs", "swdown.gs", "lwdown.gs", "psurf.gs", "qair.gs", "wind.gs", "CO2.gs")
		# if(m.name=="jules.triffid"){
		# gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(lwdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(method="optim"))
		# } else {
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(lwdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(niterEM=0, sing.tol=1e-20, method="optim"))			
		# }
		# , control=list(niterEM=0, sing.tol=1e-20, method="optim")
 	}
	if(substr(m.name,1,3)=="sib") {
		predictors <- c("tair.gs", "precipf.gs", "swdown.gs", "lwdown.gs", "psurf.gs", "qair.gs", "wind.gs", "CO2.gs")
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + s(swdown.gs, k=k) + s(lwdown.gs, k=k) + s(qair.gs, k=k) + s(psurf.gs, k=k) + s(wind.gs, k=k) + s(CO2.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(niterEM=0, sing.tol=1e-20, opt="optim"))
	}
	if(substr(m.name,1,3)=="lin") {
		predictors <- c("tair.gs", "precipf.gs")
		if(resolutions[r]=="t.100"){
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(niterEM=0, sing.tol=1e-20, opt="optim"))
		} else {
		gam1 <- gamm(NPP ~ s(tair.gs, k=k) + s(precipf.gs, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
		}
	}    
	print(summary(gam1$gam))	

	# get rid of values for predictors not used in the models for clarity later on
	data.temp[,predictors.all[!(predictors.all %in% predictors)]] <- NA

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata=data.temp)

	mod.temp <- process.gamm(gamm.model=gam1, data=data.temp, model.name=m.name, extent=ext, resolution=resolutions[r], response=response, vars=predictors, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
	if(r==1) {
		mod.out <- list()
		mod.out$data         <- mod.temp$data
		mod.out$weights      <- mod.temp$weights
		mod.out$ci.response  <- mod.temp$ci.response
		mod.out$sim.response <- mod.temp$sim.response
		mod.out$ci.terms     <- mod.temp$ci.terms
		mod.out$sim.terms    <- mod.temp$sim.terms
		mod.out[[paste("gamm", ext, substr(resolutions[r],3,nchar(paste(resolutions[r]))), sep=".")]] <- mod.temp$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         mod.temp$data)
		mod.out$weights      <- rbind(mod.out$weights,      mod.temp$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  mod.temp$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, mod.temp$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     mod.temp$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    mod.temp$sim.terms)
		mod.out[[paste("gamm", ext, substr(resolutions[r],3,nchar(paste(resolutions[r]))), sep=".")]] <- mod.temp$gamm
	}
	
} # end resolutions
save(mod.out, file=file.path(dat.dir, paste0("gamm_AllDrivers_GS_", m.name,"_", response, ".Rdata")))

m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_GS_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Resolution, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP), alpha=0.5) +
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
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP), alpha=0.5) +
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

} # end model

# save(mod.out.AGB.diff, file=file.path(out.dir, "mod.out.dAGB.Rdata"))
# save(mod.out.NPP, file=file.path(out.dir, "mod.out.NPP.Rdata"))

