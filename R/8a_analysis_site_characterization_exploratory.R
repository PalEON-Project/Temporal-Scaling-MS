# ----------------------------------------
# Sensitivity & Scaling Analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 16 November 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
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
in.base <- "Data/gamms/"
out.dir <- "Data/analyses/analysis_site_characteristics"
fig.dir <- "Figures/analyses/analysis_site_characteristics"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]
summary(ecosys)

load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP_Site.Rdata"))
summary(mod.out)

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
# Relativizing to mean NPP
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# Using mean model NPP across sites since the GAMM response curves are for 
			#    the whole model & not site-specific are parameterized
			# -----------------------
			# Find the start year for the extent
			# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			npp <- mean(dat.ecosys[dat.ecosys$Model==m  & dat.ecosys$Resolution==r, "NPP"], na.rm=T)			
			# -----------------------
			
			# -----------------------
			# Relativizing everything in dat.ecosys to make it comparable to tree rings
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

			# Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
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

# ----------------------------------------

# ----------------------------------------
# Binning Factors to plot by model and site
# ----------------------------------------
resolutions     <- c("t.010") 
extents         <- data.frame(Start=c(850), End=c(2010)) 
response        <- c("NPP", "SoilMoist")
predictors      <- c("tair.gs", "precipf.gs", "CO2.gs", "Deciduous", "Evergreen", "Grass")
factors         <- c("Model", "Model.Order", "Site", "Year")

summary(dat.ecosys)
dat.ecosys2 <- merge(dat.ecosys, ecosys[,c("Model", "Site", "Year", "Resolution", "Evergreen", "Deciduous", "Grass", "SoilMoist")])
dat.ecosys2$tair.round   <- round(dat.ecosys2$tair   , 1)
dat.ecosys2$precip.round <- round(dat.ecosys2$precipf*1.5, -1)/1.5
summary(dat.ecosys2)
length(unique(dat.ecosys3$tair.round)); length(unique(dat.ecosys3$precip.round))

m.order <- unique(dat.ecosys2$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]
# ----------------------------------------


# ----------------------------------------
# Plotting Climate Space of Sites
# ----------------------------------------
ecosys.climate <- aggregate(dat.ecosys2[dat.ecosys2$Model=="ed2", "NPP"], by=list(dat.ecosys2[dat.ecosys2$Model=="ed2", "Site"], dat.ecosys2[dat.ecosys2$Model=="ed2", "tair.round"], dat.ecosys2[dat.ecosys2$Model=="ed2", "precip.round"]), FUN=length)
names(ecosys.climate) <- c("Site", "tair", "precip", "frequency")
ecosys.climate$tair   <- ecosys.climate$tair-273.15
summary(ecosys.climate)


pdf(file.path(fig.dir, "ClimateSpace_Sites.pdf"))
ggplot(data=ecosys.climate) + theme_bw() + 
	geom_tile(aes(x=precip, y=tair, fill=Site, alpha=frequency), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)") +
	scale_y_continuous(name="Mean Annual Temperature (C)") +
	scale_alpha_continuous(range=c(0.5,1)) +
	guides(alpha=F)
dev.off()
# ----------------------------------------

# ----------------------------------------
# Plotting Climate Space & NPP of Sites
# ----------------------------------------

ecosys.npp <- aggregate(dat.ecosys2[, c("NPP", "NPP.rel", "Evergreen", "Deciduous", "Grass", "SoilMoist")], by=list(dat.ecosys2$Model, dat.ecosys2$Model.Order, dat.ecosys2$Site, dat.ecosys2$tair.round, dat.ecosys2$precip.round), FUN=mean)
names(ecosys.npp)[1:5] <- c("Model", "Model.Order", "Site", "tair", "precip")
ecosys.npp$tair   <- ecosys.npp$tair-273.15
summary(ecosys.npp)


pdf(file.path(fig.dir, "ClimateSpace_NPP.pdf"))
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=NPP), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=NPP.rel), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
for(m in unique(ecosys.npp$Model)){
print(
ggplot(data=ecosys.npp[ecosys.npp$Model==m,]) + theme_bw() + facet_wrap(~Model) +
	geom_tile(aes(x=precip, y=tair, fill=NPP), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)") +
	scale_y_continuous(name="Mean Annual Temperature (C)")
)
}
dev.off()

pdf(file.path(fig.dir, "ClimateSpace_SoilMoist.pdf"))
print(
ggplot(data=dat.ecosys2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=SoilMoist, color=Model.Order)) +
	scale_color_manual(values=paste(col.model)) + 
	theme_bw()
)
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=SoilMoist), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
for(s in unique(ecosys.npp$Site)){
print(
ggplot(data=ecosys.npp[ecosys.npp$Site==s,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=SoilMoist), size=5) +
	ggtitle(s) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
}
dev.off()
# ---------------------
pdf(file.path(fig.dir, "ClimateSpace_Fcomp.pdf"))
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Evergreen), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Deciduous), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
print(
ggplot(data=ecosys.npp[,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Grass), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
dev.off()

pdf(file.path(fig.dir, "ClimateSpace_Evergreen_bySite.pdf"))
for(s in unique(ecosys.npp$Site)){
print(
ggplot(data=ecosys.npp[ecosys.npp$Site==s,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Evergreen), size=5) +
	ggtitle(s) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
}
dev.off()

pdf(file.path(fig.dir, "ClimateSpace_Deciduous_bySite.pdf"))
for(s in unique(ecosys.npp$Site)){
print(
ggplot(data=ecosys.npp[ecosys.npp$Site==s,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Deciduous), size=5) +
	ggtitle(s) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
}
dev.off()

pdf(file.path(fig.dir, "ClimateSpace_Grass_bySite.pdf"))
for(s in unique(ecosys.npp$Site)){
print(
ggplot(data=ecosys.npp[ecosys.npp$Site==s,]) + theme_bw() + facet_wrap(~Model, scales="free") +
	geom_tile(aes(x=precip, y=tair, fill=Grass), size=5) +
	ggtitle(s) +
	scale_x_continuous(name="Mean Annual Precip (mm)", expand=c(0,0)) +
	scale_y_continuous(name="Mean Annual Temperature (C)", expand=c(0,0))
)
}
dev.off()
# ----------------------------------------


# ----------------------------------------
# Comparing Model composition of sites
# ----------------------------------------
ecosys.comp <- aggregate(dat.ecosys2$decid.round, by=list(dat.ecosys2$Model, dat.ecosys2$Model.Order, dat.ecosys2$Site, dat.ecosys2$decid.round, dat.ecosys2$evg.round), FUN=length)
names(ecosys.comp) <- c("Model", "Model.Order", "Site", "Deciduous", "Evergreen", "frequency")
summary(ecosys.comp)

pdf(file.path(fig.dir, "PFTSpace_Sites.pdf"))
ggplot(data=ecosys.comp) + facet_wrap(~Site) + theme_bw() + 
	geom_tile(aes(x=Deciduous, y=Evergreen, fill=Model.Order, alpha=frequency), size=5) +
	scale_x_continuous(name="Fraction Deciduous") +
	scale_y_continuous(name="Fraction Evergreen") +
	scale_fill_manual(values=paste(col.model)) +
	scale_alpha_continuous(range=c(0.5,1)) +
	guides(alpha=F)
dev.off()
# ----------------------------------------
