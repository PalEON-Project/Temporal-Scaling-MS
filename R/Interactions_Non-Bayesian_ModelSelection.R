# ----------------------------------------
# Temporal Scaling Analyses
# Changes in Strength of Interactions with Temporal Scale
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 29 May 2015 
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
# Set Directories
# ----------------------------------------
setwd("~/Dropbox/PalEON CR/paleon_mip_site")
inputs <- "phase1a_output_variables"
# ----------------------------------------

# ----------------------------------------
# Load Data Sets
# ----------------------------------------
# Ecosystem Model Outputs
ecosys <- read.csv(file.path(inputs, "MIP_Data_Ann_2015.csv"))
ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
summary(ecosys)

# CO2 Record
nc.co2 <- nc_open("~/Desktop/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)

# Merging CO2 into Model Outputs
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("~/Desktop/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors

model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors
# ----------------------------------------

# ----------------------------------------
# Calculate Deviations from desired reference point
# Reference Point: 0850-0869 (spinup climate)
# ----------------------------------------
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp")
vars.climate <- c("Temp", "Precip", "CO2")

# vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

ref.window <- 850:869

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		# -----------------------
		# Model Variabiles -- Relative Change
		# Deviation = percent above or below the mean for the reference window
		#             (observed-ref.mean)/ref.mean 
		# -----------------------
		for(v in unique(vars)){
			ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".dev")] <- (ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean)/ref.mean

		}
		# -----------------------

		# -----------------------
		# Climate Drivers -- Absolute, not relative change
		# Deviation = absolute deviation from reference window
		#             observed - ref.mean
		# -----------------------
		for(v in unique(vars.climate)){
			ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".abs.dev")] <- ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean
		}
		# -----------------------
	}
}

summary(ecosys)
# ----------------------------------------


# ----------------------------------------
# Perform Temporal Smoothing on Data
# Note: Smoothing is performed over the PREVIOUS 100 years becuase ecosystems 
#       cannot respond to what they have not yet experienced
# ----------------------------------------
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp", "Temp", "Precip", "CO2")
vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		# -----------------------
		# 10-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 50-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollmean(temp, k=50, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollmean(temp, k=50, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 100-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 250-Year Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollmean(temp, k=250, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollmean(temp, k=250, align="right", fill=NA)
		}
		# -----------------------

	}
}
summary(ecosys)

# -----------------------
# Some exploratory Graphing
# -----------------------
colors.use <- 
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
	geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model %in% unique(ecosys$Model.Order),"color"])) +
	theme_bw()
# -----------------------

# -----------------------
# Subsetting individual models & sites for model prototyping
# -----------------------
ed2.pha <- ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA",]
summary(ed2.pha)

lpj.g.pha <- ecosys[ecosys$Model=="lpj.guess" & ecosys$Site=="PHA",]
summary(lpj.g.pha)

lpj.g <- ecosys[ecosys$Model=="lpj.guess",]
summary(lpj.g)

# -----------------------
# ----------------------------------------

# ----------------------------------------
# Testing things out with just LPJ-GUESS by itself
# ----------------------------------------

# -----------------------
# Test with AGB
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Site=="PHA",]

lm.lpj.g.1 <- lm(AGB ~ CO2*Temp*Precip, data=lpj.g2)
summary(lm.lpj.g.1)
lpj.g2$predict.1 <- fitted(lm.lpj.g.1)
# summary(lpj.g2)

lm.lpj.g.50 <- lm(AGB.50 ~ CO2.50*Temp.50*Precip.50, data=lpj.g2)
summary(lm.lpj.g.50)
lpj.g2$predict.50 <- fitted(lm.lpj.g.50)
# summary(lpj.g2)

lm.lpj.g.100 <- lm(AGB.100 ~ CO2.100*Temp.100*Precip.100, data=lpj.g2)
summary(lm.lpj.g.100)
lpj.g2$predict.100 <- fitted(lm.lpj.g.100)
# summary(lpj.g2)


lm.lpj.g.250 <- lm(AGB.250 ~ CO2.250*Temp.250*Precip.250, data=lpj.g2)
summary(lm.lpj.g.250)
lpj.g2$predict.250 <- fitted(lm.lpj.g.250)
# summary(lpj.g2)

AIC(lm.lpj.g.1)
AIC(lm.lpj.g.50)
AIC(lm.lpj.g.100)
AIC(lm.lpj.g.250)

# par(mfrow=c(4,1))
# plot(AGB ~ Year, data=lpj.g2, main="AGB, 1-yr")
	# lines(predict.1 ~ Year, data=lpj.g2, col="red")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB), size=3) +
	geom_line(aes(x=Year, y=predict.1), color="red", size=2) +
	ggtitle("AGB, 1 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.50), size=3) +
	geom_line(aes(x=Year, y=predict.50), color="red", size=2) +
	ggtitle("AGB, 50 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.100), size=3) +
	geom_line(aes(x=Year, y=predict.100), color="red", size=2) +
	ggtitle("AGB, 100 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.250), size=3) +
	geom_line(aes(x=Year, y=predict.250), color="red", size=2) +
	ggtitle("AGB, 250 yr")

# Extracting the coefficients
coef.1 <- data.frame(var=names(coef(lm.lpj.g.1)), coef.001=coef(lm.lpj.g.1), coef.050 =coef(lm.lpj.g.50), coef.100 =coef(lm.lpj.g.100), coef.250 =coef(lm.lpj.g.250)); row.names(coef.1) <- NULL
# coef.1$var2 <- c("CO2", "Temp", "Precip", rep("Site", 6), "CO2:Temp", "CO2:Precip", "Temp:Precip" )

# coef.1$level.int <- c(rep(1, 9),)
summary(coef.1)

# Stacking the coefficients to make them easier to graph
coef.2 <- stack(coef.1[,2:5])
coef.2$var <- coef.1$var
summary(coef.2)

ggplot(data=coef.2) + facet_grid(var~., scales="free")+
	geom_point(aes(x=ind, y=values))
# -----------------------

# -----------------------
# Test with AGB, randome site effect
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g),]
library(nlme)

lm.lpj.g.1 <- lme(AGB ~ CO2*Temp*Precip, random=list(Site=~1), data=lpj.g2)
summary(lm.lpj.g.1)
lpj.g2$predict.1 <- fitted(lm.lpj.g.1)
# summary(lpj.g2)

lm.lpj.g.50 <- lme(AGB.50 ~ CO2.50*Temp.50*Precip.50, random=list(Site=~1), data=lpj.g2)
summary(lm.lpj.g.50)
lpj.g2$predict.50 <- fitted(lm.lpj.g.50)
# summary(lpj.g2)

lm.lpj.g.100 <- lme(AGB.100 ~ CO2.100*Temp.100*Precip.100, random=list(Site=~1), data=lpj.g2)
summary(lm.lpj.g.100)
lpj.g2$predict.100 <- fitted(lm.lpj.g.100)
# summary(lpj.g2)


lm.lpj.g.250 <- lme(AGB.250 ~ CO2.250*Temp.250*Precip.250, random=list(Site=~1), data=lpj.g2)
summary(lm.lpj.g.250)
lpj.g2$predict.250 <- fitted(lm.lpj.g.250)
# summary(lpj.g2)

AIC(lm.lpj.g.1)
AIC(lm.lpj.g.50)
AIC(lm.lpj.g.100)
AIC(lm.lpj.g.250)

# par(mfrow=c(4,1))
# plot(AGB ~ Year, data=lpj.g2, main="AGB, 1-yr")
	# lines(predict.1 ~ Year, data=lpj.g2, col="red")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB), size=3) +
	geom_line(aes(x=Year, y=predict.1), color="red", size=2) +
	ggtitle("AGB, 1 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.50), size=3) +
	geom_line(aes(x=Year, y=predict.50), color="red", size=2) +
	ggtitle("AGB, 50 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.100), size=3) +
	geom_line(aes(x=Year, y=predict.100), color="red", size=2) +
	ggtitle("AGB, 100 yr")

ggplot(data=lpj.g2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.250), size=3) +
	geom_line(aes(x=Year, y=predict.250), color="red", size=2) +
	ggtitle("AGB, 250 yr")

# Extracting the coefficients
coef.1 <- data.frame(var=names(coef(lm.lpj.g.1)), coef.001=coef(lm.lpj.g.1), coef.050 =coef(lm.lpj.g.50), coef.100 =coef(lm.lpj.g.100), coef.250 =coef(lm.lpj.g.250)); row.names(coef.1) <- NULL
# coef.1$var2 <- c("CO2", "Temp", "Precip", rep("Site", 6), "CO2:Temp", "CO2:Precip", "Temp:Precip" )

coef.1$level.int <- c(rep(1, 4), rep(2, 3), 3)
summary(coef.1)

# Stacking the coefficients to make them easier to graph
coef.2 <- stack(coef.1[,2:5])
coef.2$var <- coef.1$var
coef.2$level.int <- coef.1$level.int
coef.2

ggplot(data=coef.2) + facet_grid(var~level.int, scales="free")+
	geom_point(aes(x=ind, y=values))
# -----------------------


# ----------------------------------------


