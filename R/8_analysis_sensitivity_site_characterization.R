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


# # load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP_Site.Rdata"))
# mod.baseline <- mod.out; rm(mod.out)
# # names(mod.baseline)[7:length(mod.baseline)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
# summary(mod.baseline)

# load(file.path(in.base, "Sensitivity_Models_Site/gamm_models_NPP_Site.Rdata"))
# mod.site <- mod.out; rm(mod.out)
# summary(mod.site)

# ----------------------------------------

# ----------------------------------------
# Binning Factors to plot by model and site
# ----------------------------------------
resolutions     <- c("t.010") 
extents         <- data.frame(Start=c(850), End=c(2010)) 
response        <- c("NPP", "SoilMoist")
predictors      <- c("tair.gs", "precipf.gs", "CO2.gs", "Deciduous", "Evergreen", "Grass")
factors         <- c("Model", "Model.Order", "Site", "Year")

summary(ecosys)
ecosys2 <- ecosys[ecosys$Resolution==resolutions,c(factors, "Resolution", response, predictors)]
ecosys2$tair.round   <- round(ecosys2$tair.gs*3   , 1)/3
ecosys2$precip.round <- round(ecosys2$precipf.gs*3, -1)/3
ecosys2$decid.round  <- round(ecosys2$Deciduous*4, 1)/4
ecosys2$evg.round    <- round(ecosys2$Evergreen*4, 1)/4
summary(ecosys2)
length(unique(ecosys2$tair.round)); length(unique(ecosys2$precip.round))
length(unique(ecosys2$decid.round)); length(unique(ecosys2$evg.round))
# ----------------------------------------


# ----------------------------------------
# Plotting Climate Space of Sites
# ----------------------------------------
ecosys.climate <- aggregate(ecosys2[ecosys2$Model=="ed2", "tair.round"], by=list(ecosys2[ecosys2$Model=="ed2", "Site"], ecosys2[ecosys2$Model=="ed2", "tair.round"], ecosys2[ecosys2$Model=="ed2", "precip.round"]), FUN=length)
names(ecosys.climate) <- c("Site", "tair", "precip", "frequency")
ecosys.climate$tair   <- ecosys.climate$tair-273.15
summary(ecosys.climate)


pdf(file.path(fig.dir, "ClimateSpace_Sites.pdf"))
ggplot(data=ecosys.climate) + theme_bw() + 
	geom_tile(aes(x=precip, y=tair, fill=Site, alpha=frequency), size=5) +
	scale_x_continuous(name="Mean Annual Precip (mm)") +
	scale_y_continuous(name="Mean Annual Temperature (C)") +
	scale_alpha_continuous(range=c(0.7,1)) +
	guides(alpha=F)
dev.off()
# ----------------------------------------

# ----------------------------------------
# Comparing Model composition of sites
# ----------------------------------------
ecosys.comp <- aggregate(ecosys2$decid.round, by=list(ecosys2$Model, ecosys2$Model.Order, ecosys2$Site, ecosys2$decid.round, ecosys2$evg.round), FUN=length)
names(ecosys.comp) <- c("Model", "Model.Order", "Site", "Deciduous", "Evergreen", "frequency")
summary(ecosys.comp)

pdf(file.path(fig.dir, "PFTSpace_Sites.pdf"))
ggplot(data=ecosys.comp) + facet_wrap(~Site) + theme_bw() + 
	geom_tile(aes(x=Deciduous, y=Evergreen, fill=Model.Order, alpha=frequency), size=5) +
	scale_x_continuous(name="Fraction Deciduous") +
	scale_y_continuous(name="Fraction Evergreen") +
	scale_alpha_continuous(range=c(0.7,1)) +
	guides(alpha=F)
dev.off()
# ----------------------------------------
