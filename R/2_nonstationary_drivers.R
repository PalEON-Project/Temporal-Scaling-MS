# ----------------------------------------
# Temporal Scaling Analyses
# Non-constant driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: How does temporal scale affect our inference of Temp, Precip, & CO2 effects on ecosystems?
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# (Fit GAMM per site per model)
# 1) Temporal Grain (Resolution)
#    -- Fit GAMM over constant window with different degrees of smoothing (1 yr - 250 yr)
# 2) Temporal Extent 
#	 -- Fit GAMM to different windows (30 yrs 1990-2010, 100 yrs 1910-2010, full window)
# ** Response variables of interest: NPP, possibly dAGB (AGB 1st difference)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# 1) Space-Time Comparison
#    -- Hypothesis: Driver responses across sites within a model converge at coarser temporal grains 
#       and larger extents because the models have time to adjust and seek equilibrium.
#
#    -- Analysis: Use the posterior CIs for each smoothing term to see if the driver curves for sites
#                 within a model are statstically different at different sites at different scales (or 
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
#    -- Analysis: For each given scale, compare the model response curves for each driver.  Determine 
#                 which drivers are most similar/variable among models and at what scales?  Are there
#                 particular ranges of each driver response where models responses are most similar/different?
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ncdf4)
library(lme4)
library(R2jags)
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
setwd("~/Dropbox/PalEON CR/paleon_mip_site")
inputs <- "phase1a_output_variables"
fig.dir <- "~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
path.data <- "~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed model outputs
#	generated with 1_generate_ecosys.R
load(file.path(path.data, "EcosysData.Rdata"))

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2
# 	predict.gamm.model.site.R = function to run a single 1 site - model combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal model responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal model fitting with random site effects very slow
source('~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/predict.gamm.model.site.R', chdir = TRUE)
source('~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/predict.gamm.model.R', chdir = TRUE)

# ----------------------------------------


# -----------------------
# Some exploratory Graphing
# -----------------------
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
  geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
  scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  ggtitle("Annual & Centennial AGB") +
  theme_bw()
# -----------------------

# ----------------------------------------


# ----------------------------------------
# Model approach: AGB ~ 3 non-interactive temporal smoothers: AGB, Temp, Precip
# ----------------------------------------
library(mgcv)

# ------------------------------------------------
# All Sites: (for 1 site, see model selection script)
# ------------------------------------------------
outdir="~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"

# ------------------------
# LPJ-GUESS
# ------------------------
gam.lpj.guess <- model.gam(data=ecosys, model="lpj.guess", scale="", response="AGB", k=4, outdir=outdir)
summary(gam.lpj.guess)
summary(gam.lpj.guess$data)
summary(gam.lpj.guess$gam)
summary(gam.lpj.guess$weights)




pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.lpj.guess$gam)
dev.off()

lm.lpj.g <- lm(fit.gam ~ response, data=gam.lpj.guess$data)
summary(lm.lpj.g)

ggplot(data=gam.lpj.guess$data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$weights) + facet_wrap(~Site) +
	geom_line(data=gam.lpj.guess$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$data[gam.lpj.guess$data$Site=="PHA" & gam.lpj.guess$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$data[gam.lpj.guess$data$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# ------------------------

# ------------------------
# LPJ-WSL
# ------------------------
gam.lpj.wsl <- model.gam(data=ecosys, model="lpj.wsl", response="AGB", k=4, outdir)
summary(gam.lpj.wsl)
summary(gam.lpj.wsl$data)
summary(gam.lpj.wsl$gam)
summary(gam.lpj.wsl$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.lpj.wsl$gam)
dev.off()

lm.lpj.wsl <- lm(fit.gam ~ response, data=gam.lpj.wsl$data)
summary(lm.lpj.wsl)

ggplot(data=gam.lpj.wsl$data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$weights) + facet_wrap(~Site) +
	geom_line(data=gam.lpj.wsl$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$data[gam.lpj.wsl$data$Site=="PHA" & gam.lpj.wsl$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$data[gam.lpj.wsl$data$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------


# ------------------------
# LINKAGES
# ------------------------
# note: linkages didn't return temp or precip, but was updated, so we're going to assume it's the same as LPJ-Guess
data.linkages <- ecosys[ecosys$Model=="linkages",]
data.linkages$Temp <- ecosys[ecosys$Model=="lpj.guess", "Temp"]
data.linkages$Precip <- ecosys[ecosys$Model=="lpj.guess", "Precip"]
summary(data.linkages)

gam.linkages <- model.gam(data=data.linkages, model="linkages", response="AGB", k=3, outdir)
summary(gam.linkages)
summary(gam.linkages$data)
summary(gam.linkages$gam)
summary(gam.linkages$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.linkages$gam)
dev.off()

lm.linkages <- lm(fit.gam ~ response, data=gam.linkages$data)
summary(lm.linkages)

ggplot(data= gam.linkages $data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, Linkages, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$weights) + facet_wrap(~Site) +
	geom_line(data=gam.linkages$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$data[gam.linkages$data$Site=="PHA" & gam.linkages$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$data[gam.linkages$data$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------


# ------------------------
# ED2
# ------------------------
gam.ed2 <- model.gam(data=ecosys, model="ed2", response="AGB", k=3, outdir)
summary(gam.ed2)
summary(gam.ed2$data)
summary(gam.ed2$gam)
summary(gam.ed2$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.ed2$gam)
dev.off()

lm.ed2 <- lm(fit.gam ~ response, data=gam.ed2$data)
summary(lm.ed2)

ggplot(data= gam.ed2 $data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$weights) + facet_wrap(~Site) +
	geom_line(data=gam.ed2$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$data[gam.ed2$data$Site=="PHA" & gam.ed2$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$data[gam.ed2$data$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------