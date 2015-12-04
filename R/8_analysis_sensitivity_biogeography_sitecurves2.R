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
out.dir <- "Data/analyses/analysis_biogeography"
fig.dir <- "Figures/analyses/analysis_biogeography"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]


load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP_Site.Rdata"))
mod.baseline <- mod.out; rm(mod.out)
# names(mod.baseline)[7:length(mod.baseline)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
summary(mod.baseline)

load(file.path(in.base, "Sensitivity_Models_Site/gamm_models_NPP_Site.Rdata"))
mod.site <- mod.out; rm(mod.out)
summary(mod.site)

# ----------------------------------------


# ----------------------------------------
# Binning things to plot distibutions of climate & fcomp
# ----------------------------------------
resolutions <- c("t.001") # Note: Big models can't handle t.100 at the site level because there aren't enough data points
extents <- data.frame(Start=c(850), End=c(2010)) 
response <- c("NPP")
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")


summary(ecosys)
ecosys2 <- ecosys[]

# ----------------------------------------
