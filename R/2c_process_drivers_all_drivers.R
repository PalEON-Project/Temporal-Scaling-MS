# ----------------------------------------
# Temporal Scaling Analyses
# Non-constant driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 10 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: What is the relative control of different drivers within each model?
# Rationale: Not all models use all inputs, and many drivers are correlated, so we need to see if the temperature pattern is really a radiation pattern, etc. 
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# (Fit GAMM per site per m.name)
# 1) Temporal Grain (Resolution)
#    -- Fit GAMM over constant window with different degrees of smoothing (1 yr - 250 yr)
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
library(ncdf4)
library(lme4)
# library(R2jags)
library(ggplot2); library(grid)
library(car)
library(zoo)
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
dat.base="Data/gamms_byModel"
fig.base="Figures/gamms_byModel"

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData.Rdata"))

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2
# 	predict.gamm.model.site.R = function to run a single 1 site - m.name combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal m.name responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal m.name fitting with random site effects very slow
source('R/0_process.gamm.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)


# Read in model color scheme
# model.colors <- read.csv("raw_inputs/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
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


# ----------------------------------------
# Model approach: AGB ~ 3 non-interactive temporal smoothers: AGB, Temp, Precip
# ----------------------------------------
library(mgcv)

# ------------------------------------------------
# All Sites: (for 1 site, see m.name selection script)
# ------------------------------------------------

# ------------------------
# MegaLoop -- Looping through all models by Variable, by Extent
# ------------------------
# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)
scales      <- unique(ecosys$Scale)

# -----------------
# Matrix of Models and Drivers
# -----------------
# Var    ED2  ED2-LU  CLM-BGC  CLM-CN  LPJ-WSL  LPJ-GUESS  JULES-STAT  JULES-TRIFFID  SIBCASA  LINKAGES
# tair    X     X        X       X        X         X          X             X           X        X
# precip  X     X        X       X        X         X          X             X           X        X
# swdown  X     X        X       X        X         X          X             X           X
# lwdown  X     X                         X                    X             X           X
# wind    X     X        X       X                             X             X           X
# psurf   X     X        X       X                             X             X           X
# qair    X     X        X       X                             X             X           X
# co2     X     X        X       X        X         X          X             X           X
# Ndep                   ?       ?

k=4
response <- "NPP"
predictors.all <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
	
for(m in 1:length(model.name)){
	print("-------------------------------------")
	print("-------------------------------------")
	print("-------------------------------------")
	print(paste0("------ Processing Model: ", model.order[m], " ------"))
	m.name  <- model.name[m]
	m.order <- model.order[m]

	# Make sure folders for each model exist
	if(!dir.exists(file.path(dat.base, m.order))) dir.create(file.path(dat.base, m.order))
	if(!dir.exists(file.path(fig.base, m.order))) dir.create(file.path(fig, m.order))

	fig.dir <- file.path(fig.base, m.order, "AllDrivers")
	dat.dir <- file.path(dat.base, m.order, "AllDrivers")

	# Make sure the appropriate file paths are in place
	if(!dir.exists(dat.dir)) dir.create(dat.dir)
	if(!dir.exists(fig.dir)) dir.create(fig.dir)

for(t in 1:length(scales)){
	print(       "-------------------------------------")
	print(paste0("------ Processing Scale: ", scales[t], " ------"))

	data.temp <- ecosys[ecosys$Model==m.name & ecosys$Scale==scales[t], c("Model", "Updated", "Model.Order", "Site", "Year", "Scale", response, predictors.all)]
	# Making a note of the extent
	ext <- as.factor(paste(min(data.temp$Year), max(data.temp$Year), sep="-"))
	data.temp$Extent <- as.factor(ext)

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
		predictors <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
	if(substr(m.name,1,3)=="clm") {
		predictors <- c("tair", "precipf", "swdown", "psurf", "qair", "wind", "CO2")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(sing.tol=1e-20, opt="optim"))
	}
	if(substr(m.name,1,3)=="lpj") {
		predictors <- c("tair", "precipf", "swdown", "CO2")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
	if(substr(m.name,1,3)=="jul") {
		predictors <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
	if(substr(m.name,1,3)=="sib") {
		predictors <- c("tair", "precipf", "swdown", "lwdown", "psurf", "qair", "wind", "CO2")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(lwdown, k=k) + s(qair, k=k) + s(psurf, k=k) + s(wind, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
	if(substr(m.name,1,3)=="lin") {
		predictors <- c("tair", "precipf")
		gam1 <- gamm(NPP ~ s(tair, k=k) + s(precipf, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	}
    print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata=data.temp)

	mod.temp <- process.gamm(gamm.model=gam1, data=data.temp, model.name=m.name, extent=ext, scale=scales[t], response=response, vars=predictors, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
	if(t==1) {
		mod.out <- list()
		mod.out$data         <- mod.temp$data
		mod.out$weights      <- mod.temp$weights
		mod.out$ci.response  <- mod.temp$ci.response
		mod.out$sim.response <- mod.temp$sim.response
		mod.out$ci.terms     <- mod.temp$ci.terms
		mod.out$sim.terms    <- mod.temp$sim.terms
		mod.out[[paste("gamm", ext, substr(scales[t],3,nchar(paste(scales[t]))), sep=".")]] <- mod.temp$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         mod.temp$data)
		mod.out$weights      <- rbind(mod.out$weights,      mod.temp$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  mod.temp$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, mod.temp$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     mod.temp$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    mod.temp$sim.terms)
		mod.out[[paste("gamm", ext, substr(scales[t],3,nchar(paste(scales[t]))), sep=".")]] <- mod.temp$gamm
	}
	
} # end scales
save(mod.out, file=file.path(dat.dir, paste("gamm", m.name, response, "Rdata", sep=".")))

m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_", m.order, "_", response, "_0850-2010", ".pdf")))
print(
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Scale, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	# scale_x_continuous(limits=c(850,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_AllDrivers_", m.order, "_", response, "_1900-2010", ".pdf")))
print(	
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Scale, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
	scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=col.model) +
	# scale_color_manual(values=col.model) +		
	labs(title=paste(m.order, response, sep=" - "), x="Year", y=response)
)
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_AllDrivers_", m.order, "_", response, ".pdf")))
print(
ggplot(data=mod.out$ci.terms[,]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Scale), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Scale), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_color_manual(values=c("red2", "blue", "green3")) +
	# scale_fill_manual(values=c("red2", "blue", "green3")) +
	labs(title=paste0("Driver Effects: ",m.order), y="Effect Size") # +
)
dev.off()

} # end model

# save(mod.out.AGB.diff, file=file.path(out.dir, "mod.out.dAGB.Rdata"))
# save(mod.out.NPP, file=file.path(out.dir, "mod.out.NPP.Rdata"))

