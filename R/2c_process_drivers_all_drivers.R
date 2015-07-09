# ----------------------------------------
# Temporal Scaling Analyses
# Non-constant driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
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
setwd("~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
path.data <- "~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
fig.dir <- "~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path(path.data, "EcosysData.Rdata"))

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2
# 	predict.gamm.model.site.R = function to run a single 1 site - m.name combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal m.name responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal m.name fitting with random site effects very slow
source('~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_predict.gamm.model.site.R', chdir = TRUE)
source('~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_predict.gamm.model.R', chdir = TRUE)
source('~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_GAMM_Plots.R', chdir = TRUE)


# Read in model color scheme
model.colors <- read.csv("~/Dropbox/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
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
data.base="~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms_byModel"
fig.base="~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/gamms_byModel"

# ------------------------
# MegaLoop -- Looping through all models by Variable, by Extent
# ------------------------
# Setting up a loop for 1 m.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- unique(ecosys$Model)
model.order   <- unique(ecosys$Model.Order)

#### crashed on jules.triffid at one point; need to re-run that one alone
# model.name <- model.name[5:length(model.name)]
# model.order <- model.order[5:length(model.order)]


# var <- c("NPP", "AGB.diff")
var <- "NPP"
# scale    <- ""
scales <- c("", ".10", ".50", ".100", ".250")
# scales <- c(".100")
t.scales <- ifelse(scales=="", "t.001", paste0("t", scales))
extents <- data.frame(Start=c(850, 1900, 1990), End=c(2010, 2010, 2010)) 
extent <- c(850, 2010)
# -----------------
# Models are having varying stability issues in the gamm, so lets explicitly state 
# which set of gamm/lme controls to use
# Note: because each model uses different sets of drivers, we're jsut going to run them completely separately so we can explicitly state what's in it
# -----------------
control0 <- list()
control1 <- list(sing.tol=1e-20, opt="optim")
control2 <- list(niterEM=0,sing.tol=1e-20)
control3 <- list(niterEM=0, sing.tol=1e-20, opt="optim")

# NPP Settings
control.settings <- data.frame(Model=model.name, Control=c(rep("control3", 4), "control1", rep("control3", 2), rep("control1", 3)))
# ed2 = control3


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
# Ndep                           X




# ----------------------------------------
# ED2
# ----------------------------------------
# creating a working data.temp frame with just the data.temp we want
# creating a working data.temp frame with just the data.temp we want
 k=4
for(v in var){
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Variable: ",v, " ------"))
	ecosys2 <- ecosys[ecosys$Model=="ed2" & ecosys$Year >= 1900, ]
	ecosys2$Extent <- as.factor(paste(min(ecosys$Year), max(ecosys$Year), sep="-"))
	
for(t in 1:length(scales)){
	print(       "-------------------------------------")
	print(paste0("------ Processing Scale: ", t.scales[t], " ------"))

	t.scale <- ifelse(scale=="", "t.001", paste0("t", scale))

	data.temp          <- ecosys2[, c("Model", "Updated", "Model.Order", "Site", "Extent", "Year")]
	data.temp$Scale    <- as.factor(t.scale)
	data.temp$response <- ecosys2[, paste0(v, scales[t])]
	
	# These are the drivers used in each model
	# predictors <- c("tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair", "co2")
 	predictors <- c("Temp", "Precip", "CO2")
	for(p in predictors){
		data.temp[,p]  <- ecosys2[,paste0(p, scales[t])]
	}
	
	# Each of the models is having different stability issues
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(sing.tol=1e-20, opt="optim"))
  print(summary(gam1$gam))	


	mod.temp <- process.gamm(gamm.model=gam1, data=data.temp, model.name=m.name, extent=extent, scale=scales[t], response=v, vars=predictors, k=4, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T, control=control1)
	
	if(t==1 & e==1) {
		mod.out <- list()
		mod.out$data         <- mod.temp$data
		mod.out$weights      <- mod.temp$weights
		mod.out$ci.response  <- mod.temp$ci.response
		mod.out$sim.response <- mod.temp$sim.response
		mod.out$ci.terms     <- mod.temp$ci.terms
		mod.out$sim.terms    <- mod.temp$sim.terms
		mod.out[[paste("gamm", paste0(extent[1], "-", extent[2]), t.scales[t], sep=".")]] <- mod.temp$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         mod.temp$data)
		mod.out$weights      <- rbind(mod.out$weights,      mod.temp$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  mod.temp$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, mod.temp$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     mod.temp$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    mod.temp$sim.terms)
		mod.out[[paste("gamm", paste0(extent[1], "-", extent[2]), t.scales[t], sep=".")]] <- mod.temp$gamm
	}
	
	} # end scales
} # end extent
	save(mod.out, file=file.path(out.dir, paste("gamm", m.name, v, "Rdata", sep=".")))
	# assign(paste("gamm", m.name, v, sep="."), mod.out)
col.model <- model.colors[model.colors$Model.Order %in% unique(mod.out$data$Model.Order),"color"]

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Extent_", m.order, "_", v, "_0850-2010", ".pdf")))
print(
ggplot(data=mod.out$ci.response[mod.out$ci.response$Scale=="t.001",]) + facet_grid(Site~Extent, scales="free") + theme_bw() +
 		geom_line(data= mod.out$data[mod.out$data$Scale=="t.001",], aes(x=Year, y=response), alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		# scale_x_continuous(limits=c(850,2010)) +
		# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.order), x="Year", y=var)
)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Extent_", m.order, "_", v, "_1990-2010", ".pdf")))
print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Scale=="t.001",]) + facet_wrap(~Site, scales="free") + theme_bw() +
 		geom_line(data=mod.out$data[mod.out$data$Scale=="t.001",], aes(x=Year, y=response), size=1.5, alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
		geom_line(aes(x=Year, y=mean, color=Extent), size=1) +
		scale_x_continuous(limits=c(1990,2010)) +
		# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.order), x="Year", y=var)
)
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_Extent_", m.order, "_", v, ".pdf")))
print(
ggplot(data=mod.out$ci.terms[mod.out$ci.terms$Scale=="t.001",]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Extent), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title=paste0("Driver Effects: ",m.order), y="Effect Size") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
)
dev.off()


} # end var
} # end model

# save(mod.out.AGB.diff, file=file.path(out.dir, "mod.out.dAGB.Rdata"))
# save(mod.out.NPP, file=file.path(out.dir, "mod.out.NPP.Rdata"))

