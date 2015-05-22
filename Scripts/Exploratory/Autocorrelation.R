# ---------------------------------------------------------------------------------------
# PalEON Project Paper: Temporal Scales
# Assessing Temporal Autocorrelation
# Christy Rollinson, Spring 2015, crollinson@gmail.com
# ---------------------------------------------------------------------------------------
#
# Objective: We want to know what environmental drivers (temp, precip etc) driver variation in 
# 	ecosystem properties (AGB, NPP) at different three temporal scales: annual, decadal, and centennial. 
#	To accomplish this, we need to deal with temporal autocorrelation, particularly at the annual time scales.
#
# We need to describe:
#	1) The temporal autocorrelation structure of each variable in each model (e.g. AGB probably has should have larger autocorrelation than NPP)
#	2) Quantify the effects of each driver in the context of temporal autocorrelation
#
# An example equation may be: NPP ~ Model*precip*temp*dlwf - 1, random=list(Site=~1), correlation=ARMA(??? form =1|Year)????
#
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# Set directories, load packages, read in data
# ---------------------------------------------------------------------------------------
setwd("../..")

library(ggplot2)
library(grid)

ecosys <- read.csv("phase1a_output_variables/PalEON_MIP_AGB_LAI_NPP_Yearly_2.csv")
summary(ecosys)

met.year <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Year.csv")
summary(met.year)

ecosys <- merge(ecosys, met.year, all.x=T, all.y=T)
summary(ecosys)

par(mfrow=c(6,1))
acf(ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA", "NPP"], 500, main="NPP, ED2")
acf(ecosys[ecosys$Model=="clm45" & ecosys$Site=="PHA", "NPP"], 500, main="NPP, CLM")
acf(ecosys[ecosys$Model=="jules.stat" & ecosys$Site=="PHA" & !is.na(ecosys$NPP), "NPP"], 500, main="NPP, Jules")
acf(ecosys[ecosys$Model=="lpj.guess" & ecosys$Site=="PHA", "NPP"], 500, main="NPP, LPJ-GUESS")
acf(ecosys[ecosys$Model=="lpj.wsl" & ecosys$Site=="PHA", "NPP"], 500, main="NPP, LPJ-WSL")
acf(ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA", "tair"], 500, main="Tair")

ccf(ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA", "NPP"],ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA", "tair"], 500, main="ED NPP vs Tair")

# ---------------------------------------------------------------------------------------
