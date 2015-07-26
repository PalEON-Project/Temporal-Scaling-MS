# ----------------------------------------
# Temporal Scaling Analyses -- Year vs. Growing Season
# Looking at whether temporal resolution or extent has greater impact on model responses to drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Does annual or growing season make a more accurate prediction of NPP ()
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Base-level gams (850-2010, annual resolution)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) Compare R2 value extracted from gamm
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
# setwd("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
path.data <- "Data"
in.base <- "Data/gamms"
in.yr  <- "AllDrivers_Yr_byResolution"
in.gs  <- "AllDrivers_GS_byResolution"
out.dir <- "Data/analysis_response_season"
fig.dir <- "Figures/analysis_response_season"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData_Raw.Rdata"))
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

models <- unique(ecosys$Model)
f.yr <- dir(file.path(in.base, in.yr))
f.gs <- dir(file.path(in.base, in.gs))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")

gam.summary <- data.frame(Model=models, Yr.R2=NA, GS.R2=NA, Yr.Dev=NA, GS.Dev=NA)
for(i in 1:length(models)){
	#-----------------
	# doing year first
	#-----------------
	fmod <- grep(models2[i], f.yr)
	load(file.path(in.base, in.yr, f.yr[fmod]))

	gam.summary[i,"Yr.R2"] <- summary(mod.out[["gamm.850-2010.001"]])$r.sq
	gam.summary[i,"Yr.Dev"] <- summary(mod.out[["gamm.850-2010.001"]])$dev
	#-----------------

	#-----------------
	# next growing season
	#-----------------
	fmod <- grep(models2[i], f.gs)
	load(file.path(in.base, in.gs, f.gs[fmod]))

	gam.summary[i,"GS.R2"] <- summary(mod.out[["gamm.850-2010.001"]])$r.sq
	gam.summary[i,"GS.Dev"] <- summary(mod.out[["gamm.850-2010.001"]])$dev
	#-----------------
	# Clear the mod.out to save space
	rm(mod.out)
}
gam.summary

# Summary of change in R2
mean(gam.summary$GS.R2-gam.summary$Yr.R2); sd(gam.summary$GS.R2-gam.summary$Yr.R2)

# Summary of change in Deviance
mean(gam.summary$GS.Dev-gam.summary$Yr.Dev); sd(gam.summary$GS.Dev-gam.summary$Yr.Dev)

write.csv(gam.summary, file.path(out.dir, "NPP_Yr_vs_GS.csv"), row.names=F)
# ----------------------------------------
