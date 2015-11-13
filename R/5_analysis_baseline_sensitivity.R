# ----------------------------------------
# Temporal Scaling Analyses
# Looking at whether temporal resolution or extent has greater impact on model responses to drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Which is the more important consideration in understanding ecosystem variation in response to drivers
#           of change: temporal resolution, or temporal scale?
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Individual model response curves to driver
#    -- generated with 2c_process_drivers_all_drivers.R & 2e_process_drivers_all_drivers_byExtent.R
#    -- NOTE: Taking growing season analysis because it had a slightly higher explanatory power
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) Determine & compare relative weight of individual drivers across model
#    -- pay special attention to the role of Tair, Precip & CO2
# B) Determine the change in driver response with change in scale
#    -- relativize to the annual res, 850-2010 extent curves
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(car)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
# setwd("..")
path.data <- "Data"
in.base <- "Data/gamms/Sensitivity_Models_Baseline"
out.dir <- "Data/analyses/analysis_baseline_sensitivity"
fig.dir <- "Figures/analyses/analysis_baseline_sensitivity"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]


load(file.path(in.base, "gamm_models_baseline_NPP_Site.Rdata"))

ci.terms   <- mod.out$ci.terms
dat.ecosys <- cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")])
wt.terms   <- mod.out$weights
sim.terms  <- mod.out$sim.terms

sim.terms$Effect  <- as.factor(sim.terms$Effect)
ci.terms$Extent   <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
wt.terms$Extent   <- as.factor(ifelse(wt.terms$Extent=="850-2010", "0850-2010", paste(wt.terms$Extent)))

summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)
summary(sim.terms[,1:10])
# ----------------------------------------



# ----------------------------------------
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){
		# for(s in unique(ci.terms[ci.terms$Model==m & ci.terms$Resolution==r, "Site"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# Using 1800:1850 mean (Settlement-era) as reference point
			# -----------------------
			# Find the start year for the extent
			# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year %in% 1800:1850 & dat.ecosys$Resolution==r, "NPP"], na.rm=T)			# -----------------------
			
			# -----------------------
			# Relativizing everythign in dat.ecosys to make it comparable to tree rings
			# -----------------------
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr"])/npp
			# -----------------------

			
			# -----------------------
			# Finding the percent change in NPP relative to the mean for that particular scale
			# -----------------------
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr"])/npp

			sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,7:ncol(sim.terms)] <- (sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,7:ncol(sim.terms)])/npp
			# -----------------------

			# -----------------------
			# Relativizing the factor fits through times and weights as well
			# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
			# -----------------------
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.tair.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.tair"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.precipf.rel"] <- 1+ (wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.precipf"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.CO2.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.CO2"])/npp
			# -----------------------

		# }
	}
}
summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])
# ----------------------------------------


# ----------------------------------------
# Comparing baseline sensitivity curves by Model
# ----------------------------------------
ggplot(data=ci.terms) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
# ----------------------------------------

# ----------------------------------------
# Aggregating to look at ensemble sensitivity
# ----------------------------------------
sim.terms[sim.terms$Effect=="CO2", "x"] <- round(sim.terms[sim.terms$Effect=="CO2", "x"], -1)
summary(sim.terms[sim.terms$Effect=="CO2", "x"])


x.CO2    <- unique(sim.terms[sim.terms$Effect=="CO2", "x"])
x.tair   <- unique(sim.terms[sim.terms$Effect=="tair", "x"])
x.precipf <- unique(sim.terms[sim.terms$Effect=="precipf", "x"])

terms.ensem <- data.frame(Effect=c(rep("tair", length(x.tair)), rep("precipf", length(x.precipf)), rep("CO2", length(x.CO2))),
                          x=c(x.tair, x.precipf, x.CO2))

for(e in unique(terms.ensem$Effect)){
	for(p in unique(terms.ensem[terms.ensem$Effect==e,"x"])){
		rows.use <- which(sim.terms$Effect==e & sim.terms$x==p)
		terms.ensem[terms.ensem$Effect==e & terms.ensem$x==p, "mean"] <-     mean(as.matrix(sim.terms[rows.use,7:ncol(sim.terms)]), na.rm=T)
		terms.ensem[terms.ensem$Effect==e & terms.ensem$x==p, "lwr"]  <- quantile(as.matrix(sim.terms[rows.use,7:ncol(sim.terms)]), 0.025, na.rm=T)
		terms.ensem[terms.ensem$Effect==e & terms.ensem$x==p, "upr"]  <- quantile(as.matrix(sim.terms[rows.use,7:ncol(sim.terms)]), 0.975, na.rm=T)
	}
}
summary(terms.ensem)

ggplot(data= terms.ensem) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Effect), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Effect), size=2) +
	scale_fill_manual(values=c("green3", "blue", "red")) +
	scale_color_manual(values=c("green3", "blue", "red")) +
	theme_bw()
	
# ----------------------------------------
