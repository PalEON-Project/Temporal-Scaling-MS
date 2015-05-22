# ----------------------------------------
# Temporal Scaling Analyses
# Non-Stationary Drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 7 May 2015 
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
# ----------------------------------------

# ----------------------------------------
# Load Data Sets
# ----------------------------------------
# Ecosystem Model Outputs
ecosys <- read.csv(file.path(inputs, "MIP_Data_Ann_2015.csv"))
ecosys$Model.Order <- recode(ecosys$Model, "'ed2'='1'; 'ed2.lu'='2'; 'clm45'='3'; 'lpj.wsl'='4'; 'lpj.guess'='5'; 'jules.stat'='6'; 'SiB'='7'; 'linkages'='8'")
levels(ecosys$Model.Order) <- c("ED2", "ED2-LU", "CLM4.5", "LPJ-WSL", "LPJ-GUESS", "JULES", "SiBCASA", "LINKAGES")
summary(ecosys)

# CO2 Record
nc.co2 <- nc_open("~/Dropbox/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)

# Merging CO2 into Model Outputs
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("~/Dropbox/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors$Model.Order <- recode(model.colors$Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'; 'LINKAGES'='8'")
levels(model.colors$Model.Order)[1:8] <- c("ED2", "ED2-LU", "CLM4.5", "LPJ-WSL", "LPJ-GUESS", "JULES", "SiBCASA", "LINKAGES")
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
		# Decadal Smoothing
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
		# Centennial Smoothing
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
# -----------------------
# ----------------------------------------



# ----------------------------------------
# Model approach: AGB ~ 3 non-interactive temporal smoothers: AGB, Temp, Precip
# ----------------------------------------
library(mgcv)

summary(lpj.g.pha)

# ------------------------------------------------
# 1 site:
# ------------------------------------------------
gam.lpj.g.pha <- gam(AGB ~ s(Year, by=CO2) + s(Year, by=Temp) + s(Year, by=Precip), data=lpj.g.pha)
summary(gam.lpj.g.pha)

par(mfrow=c(4,1))
#plot(gam.lpj.g.pha$gam)
plot(gam.lpj.g.pha)
acf(gam.lpj.g.pha$resid)

summary(gam.lpj.g.pha)

# -----------
# Trying to get the non-linear parameter estimates & variance
# Working from examples in predict.gam (?predict.gam)
# -----------
# Create the prediction matrix
Xp <- predict(gam.lpj.g.pha, newdata=lpj.g.pha, type="lpmatrix")
summary(Xp)

# just a quick example of how Xp times the coeff gets the predicted values 
fit5 <- Xp %*% coef(gam.lpj.g.pha)
summary(fit5)

plot(lpj.g.pha$AGB, type="l", lwd=2)
	lines(fit5, col="red")


d=10 # num. terms per effect
fit.int <- Xp[,1] * coef(gam.lpj.g)[1]
fit.co2 <- Xp[,(1+1*d+1-d):(1*d+1)] %*% coef(gam.lpj.g.pha)[(1+1*d+1-d):(1*d+1)]
fit.temp <- Xp[,(1+2*d+1-d):(2*d+1)] %*% coef(gam.lpj.g.pha)[(1+2*d+1-d):(2*d+1)]
fit.precip <- Xp[,(1+3*d+1-d):(3*d+1)] %*% coef(gam.lpj.g.pha)[(1+3*d+1-d):(3*d+1)]

fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
fit.spline <- fit.co2 + fit.temp + fit.precip
summary(fit.sum)
summary(fit5 - fit.sum)


fit.df <- data.frame(fit=fit.sum, intercept=fit.int, fit.spline=fit.spline, co2=fit.co2, temp=fit.temp, precip=fit.precip)
par(mfrow=c(3,1))
plot(gam.lpj.g.pha)

par(mfrow=c(1,1))
plot(fit.df$fit, type="l", lwd=2, ylim=range(fit.df))
	lines(fit.df$fit.int, col="gray50")
	lines(fit.df$co2, col="green3")
	lines(fit.df$temp, col="red2")
	lines(fit.df$precip, col="blue")

fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

factor.weights <- data.frame(Year = lpj.g.pha$Year, fit=fit.sum, intercept=fit.int/fit.spline2, fit.spline=fit.spline, co2=fit.co2/fit.spline2, temp=fit.temp/fit.spline2, precip=fit.precip/fit.spline2)
summary(factor.weights)


for(i in 1:nrow(factor.weights)){
	fweight <- abs(factor.weights[i,c("intercept", "co2", "temp", "precip")])
	factor.weights[i,"max"] <- max(fweight)
	factor.weights[i,"factor.max"] <- names(factor.weights[,3:6])[which(fweight==max(fweight))]
}
factor.weights$factor.max <- as.factor(factor.weights$factor.max)
summary(factor.weights)


ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), color=rgb(abs(factor.weights$temp),abs(factor.weights$co2),abs(factor.weights$precip)), size=4) +
	labs(y="Year", x="AGB kg/m2", title="Non-Stationary Drivers of AGB") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.75)), axis.text.y=element_text(color="black", size=rel(1.75)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.title=element_text(face="bold", size=rel(3)))
# ------------------------------------------------


# ------------------------------------------------
# All Sites:
# ------------------------------------------------

lm.lpj.g <- lm(AGB ~ CO2 + Temp + Precip + Site -1, data=lpj.g)
summary(lm.lpj.g)
lpj.g$predict <- fitted(lm.lpj.g)
summary(lpj.g)

ggplot(data=lpj.g) +
	geom_point(aes(x=AGB, y=predict, color=Site))

ggplot(data=lpj.g) + facet_wrap(~Site) +
	geom_point(aes(x=AGB, y=predict))

ggplot(data=lpj.g) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB), size=2) +
	geom_line(aes(x=Year, y=predict), color="red")


gam.lpj.g <- gam(AGB ~ s(Year, by=CO2) + s(Year, by=Temp) + s(Year, by=Precip) + Site -1, data=lpj.g)
summary(gam.lpj.g)
lpj.g$fit.gam <- fitted(gam.lpj.g)
summary(lpj.g)

pdf(file.path(fig.dir, "Non-StationaryDrivers_Drivers_Annual.pdf"), width=11, height=8.5)
ggplot(data=lpj.g) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Temp), size=1, color="red") +
	geom_line(aes(x=Year, y=Precip*sec2yr*.01+265), color="blue") +
	geom_line(aes(x=Year, y=CO2*.1+255), color="green", size=2)	 +
	ggtitle("Model Drivers - Annual") +
	theme_bw()
dev.off()

ggplot(data=lpj.g) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB), size=2) +
	geom_line(aes(x=Year, y=predict), color="orange") +
	geom_line(aes(x=Year, y=fit.gam), color="purple") +
	ggtitle("Model Responses")


# par(mfrow=c(4,1))
# #plot(gam.lpj.g.pha$gam)
# plot(gam.lpj.g)
# acf(gam.lpj.g$resid)

# par(mfrow=c(4,1))
# #plot(gam.lpj.g.pha$gam)
# plot(gam.lpj.g, ylim=c(-0.1,0.1))
# acf(gam.lpj.g$resid)

summary(gam.lpj.g)

# -----------
# Trying to get the non-linear parameter estimates & variance
# Working from examples in predict.gam (?predict.gam)
# -----------
# Create the prediction matrix
Xp <- predict(gam.lpj.g, newdata=lpj.g, type="lpmatrix")
summary(Xp)

# just a quick example of how Xp times the coeff gets the predicted values 
fit5 <- Xp %*% coef(gam.lpj.g)
summary(fit5)

d=10 # num. terms per effect
fit.int <- Xp[,1:6] %*% coef(gam.lpj.g)[1:6]
fit.co2 <- Xp[,(6+1*d+1-d):(1*d+6)] %*% coef(gam.lpj.g)[(6+1*d+1-d):(1*d+6)]
fit.temp <- Xp[,(6+2*d+1-d):(2*d+6)] %*% coef(gam.lpj.g)[(6+2*d+1-d):(2*d+6)]
fit.precip <- Xp[,(6+3*d+1-d):(3*d+6)] %*% coef(gam.lpj.g)[(6+3*d+1-d):(3*d+6)]

fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
fit.spline <- fit.co2 + fit.temp + fit.precip
summary(fit.sum)
summary(fit.spline)
summary(fit5 - fit.sum)


fit.df <- data.frame(Site=lpj.g$Site, Year=lpj.g$Year, fit=fit.sum, intercept=fit.int, fit.spline=fit.spline, co2=fit.co2, temp=fit.temp, precip=fit.precip)
# par(mfrow=c(3,1))
# plot(gam.lpj.g)
summary(fit.df)

ggplot(data=fit.df[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=intercept), color="gray50", size=1.5) +
	geom_line(aes(x=Year, y=co2), color="green4", size=1.0) +
	geom_line(aes(x=Year, y=temp), color="red2", size=1.0) +
	geom_line(aes(x=Year, y=precip), color="blue", size=1.0)

fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

factor.weights <- data.frame(Year = lpj.g$Year, Site=lpj.g$Site, fit=fit.sum, intercept=fit.int/fit.spline2, fit.spline=fit.spline, co2=fit.co2/fit.spline2, temp=fit.temp/fit.spline2, precip=fit.precip/fit.spline2)
summary(factor.weights)


for(i in 1:nrow(factor.weights)){
	fweight <- abs(factor.weights[i,c("intercept", "co2", "temp", "precip")])
	factor.weights[i,"max"] <- max(fweight)
	factor.weights[i,"factor.max"] <- names(factor.weights[,3:6])[which(fweight==max(fweight))]
}
factor.weights$factor.max <- as.factor(factor.weights$factor.max)
summary(factor.weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_Annual.pdf"), width=11, height=8.5)
ggplot(data=factor.weights) + facet_wrap(~Site) +
	geom_line(data=lpj.g, aes(x=Year, y=AGB), color="gray50", size=2) +
	geom_line(data=factor.weights[factor.weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHA","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHO","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PUN","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PBL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PDL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PMB","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# ------------------------------------------------
# All Sites: 10-yrs
# ------------------------------------------------
lpj.g.2 <- lpj.g[!is.na(lpj.g$AGB.10),]

lm.lpj.g.2 <- lm(AGB.10 ~ CO2.10 + Temp.10 + Precip.10 + Site -1, data=lpj.g.2)
summary(lm.lpj.g.2)

lpj.g.2$fit.lm <- fitted(lm.lpj.g.2)
summary(lpj.g.2)

ggplot(data=lpj.g.2) +
	geom_point(aes(x=AGB.10, y= fit.lm, color=Site))

ggplot(data=lpj.g.2) + facet_wrap(~Site) +
	geom_point(aes(x=AGB.10, y=fit.lm))

ggplot(data=lpj.g.2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.10), size=2) +
	geom_line(aes(x=Year, y=fit.lm), color="orange")


gam.lpj.g.2 <- gam(AGB ~ s(Year, by=CO2.10) + s(Year, by=Temp.10) + s(Year, by=Precip.10) + Site -1, data=lpj.g.2)
summary(gam.lpj.g.2)
lpj.g.2$fit.gam <- fitted(gam.lpj.g.2)
summary(lpj.g.2)


pdf(file.path(fig.dir, "Non-StationaryDrivers_Drivers_Decadal.pdf"), width=11, height=8.5)
ggplot(data=lpj.g.2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Temp.10), color="red", size=2) +
	geom_line(aes(x=Year, y=Precip.10*sec2yr*.01+265), color="blue", size=2) +
	geom_line(aes(x=Year, y=CO2.10*.1+255), color="green", size=2)	 +
	ggtitle("Model Drivers - Decadal") +
	theme_bw()
dev.off()

ggplot(data=lpj.g.2) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.10), size=3) +
	geom_line(aes(x=Year, y=fit.lm), color="orange", size=2) +
	geom_line(aes(x=Year, y=fit.gam), color="purple", size=2) +
	ggtitle("Model Responses")


par(mfrow=c(3,1))
plot(gam.lpj.g.2)
# acf(gam.lpj.g.2$resid)

# par(mfrow=c(4,1))
# #plot(gam.lpj.g.2.pha$gam)
# plot(gam.lpj.g.2, ylim=c(-0.1,0.1))
# acf(gam.lpj.g.2$resid)

summary(gam.lpj.g.2)

# -----------
# Trying to get the non-linear parameter estimates & variance
# Working from examples in predict.gam (?predict.gam)
# -----------
# Create the prediction matrix
Xp <- predict(gam.lpj.g.2, newdata=lpj.g.2, type="lpmatrix")
summary(Xp)

# just a quick example of how Xp times the coeff gets the predicted values 
fit5 <- Xp %*% coef(gam.lpj.g.2)
summary(fit5)

d=10 # num. terms per effect
fit.int <- Xp[,1:6] %*% coef(gam.lpj.g.2)[1:6]
fit.co2 <- Xp[,(6+1*d+1-d):(1*d+6)] %*% coef(gam.lpj.g.2)[(6+1*d+1-d):(1*d+6)]
fit.temp <- Xp[,(6+2*d+1-d):(2*d+6)] %*% coef(gam.lpj.g.2)[(6+2*d+1-d):(2*d+6)]
fit.precip <- Xp[,(6+3*d+1-d):(3*d+6)] %*% coef(gam.lpj.g.2)[(6+3*d+1-d):(3*d+6)]

fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
fit.spline <- fit.co2 + fit.temp + fit.precip
summary(fit.sum)
summary(fit.spline)
summary(fit5 - fit.sum)


fit.df <- data.frame(Site=lpj.g.2$Site, Year=lpj.g.2$Year, fit=fit.sum, intercept=fit.int, fit.spline=fit.spline, co2=fit.co2, temp=fit.temp, precip=fit.precip)
# par(mfrow=c(3,1))
# plot(gam.lpj.g.2)
summary(fit.df)

ggplot(data=fit.df[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=intercept), color="gray50", size=1.5) +
	geom_line(aes(x=Year, y=co2), color="green4", size=1.0) +
	geom_line(aes(x=Year, y=temp), color="red2", size=1.0) +
	geom_line(aes(x=Year, y=precip), color="blue", size=1.0)

fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

factor.weights <- data.frame(Year = lpj.g.2$Year, Site=lpj.g.2$Site, fit=fit.sum, intercept=fit.int/fit.spline2, fit.spline=fit.spline, co2=fit.co2/fit.spline2, temp=fit.temp/fit.spline2, precip=fit.precip/fit.spline2)
summary(factor.weights)


for(i in 1:nrow(factor.weights)){
	fweight <- abs(factor.weights[i,c("intercept", "co2", "temp", "precip")])
	factor.weights[i,"max"] <- max(fweight)
	factor.weights[i,"factor.max"] <- names(factor.weights[,3:6])[which(fweight==max(fweight))]
}
factor.weights$factor.max <- as.factor(factor.weights$factor.max)
summary(factor.weights)


pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_Decadal.pdf"), width=11, height=8.5)
ggplot(data=factor.weights) + facet_wrap(~Site) +
	geom_line(data=lpj.g.2, aes(x=Year, y=AGB.10), color="gray50", size=2) +
	geom_line(data=factor.weights[factor.weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHA","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHO","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PUN","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PBL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PDL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PMB","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Decadal") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------------------------------

# ------------------------------------------------
# All Sites: 100-yrs
# ------------------------------------------------
lpj.g.3 <- lpj.g[!is.na(lpj.g$AGB.100),]

lm.lpj.g.3 <- lm(AGB.100 ~ CO2.100 + Temp.100 + Precip.100 + Site -1, data=lpj.g.3)
summary(lm.lpj.g.3)

lpj.g.3$fit.lm <- fitted(lm.lpj.g.3)
summary(lpj.g.3)

ggplot(data=lpj.g.3) +
	geom_point(aes(x=AGB.100, y= fit.lm, color=Site))

ggplot(data=lpj.g.3) + facet_wrap(~Site) +
	geom_point(aes(x=AGB.100, y=fit.lm))

ggplot(data=lpj.g.3) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.100), size=2) +
	geom_line(aes(x=Year, y=fit.lm), color="orange")


gam.lpj.g.3 <- gam(AGB.100 ~ s(Year, by=CO2.100) + s(Year, by=Temp.100) + s(Year, by=Precip.100) + Site -1, data=lpj.g.3)
summary(gam.lpj.g.3)
lpj.g.3$fit.gam <- fitted(gam.lpj.g.3)
summary(lpj.g.3)


pdf(file.path(fig.dir, "Non-StationaryDrivers_Drivers_Centennial.pdf"), width=11, height=8.5)
ggplot(data=lpj.g.3) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Temp.100), color="red", size=2) +
	geom_line(aes(x=Year, y=Precip.100*sec2yr*.01+265), color="blue", size=2) +
	geom_line(aes(x=Year, y=CO2.100*.1+255), color="green", size=2)	 +
	ggtitle("Model Drivers - Centennial") +
	theme_bw()
dev.off()

ggplot(data=lpj.g.3) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB.100), size=3) +
	geom_line(aes(x=Year, y=fit.lm), color="orange", size=2) +
	geom_line(aes(x=Year, y=fit.gam), color="purple", size=2) +
	ggtitle("Model Responses")


par(mfrow=c(3,1))
plot(gam.lpj.g.3)
# acf(gam.lpj.g.3$resid)

# par(mfrow=c(4,1))
# #plot(gam.lpj.g.3.pha$gam)
# plot(gam.lpj.g.3, ylim=c(-0.1,0.1))
# acf(gam.lpj.g.3$resid)

summary(gam.lpj.g.3)

# -----------
# Trying to get the non-linear parameter estimates & variance
# Working from examples in predict.gam (?predict.gam)
# -----------
# Create the prediction matrix
Xp <- predict(gam.lpj.g.3, newdata=lpj.g.3, type="lpmatrix")
summary(Xp)

# just a quick example of how Xp times the coeff gets the predicted values 
fit5 <- Xp %*% coef(gam.lpj.g.3)
summary(fit5)

d=10 # num. terms per effect
fit.int <- Xp[,1:6] %*% coef(gam.lpj.g.3)[1:6]
fit.co2 <- Xp[,(6+1*d+1-d):(1*d+6)] %*% coef(gam.lpj.g.3)[(6+1*d+1-d):(1*d+6)]
fit.temp <- Xp[,(6+2*d+1-d):(2*d+6)] %*% coef(gam.lpj.g.3)[(6+2*d+1-d):(2*d+6)]
fit.precip <- Xp[,(6+3*d+1-d):(3*d+6)] %*% coef(gam.lpj.g.3)[(6+3*d+1-d):(3*d+6)]

fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
fit.spline <- fit.co2 + fit.temp + fit.precip
summary(fit.sum)
summary(fit.spline)
summary(fit5 - fit.sum)


fit.df <- data.frame(Site=lpj.g.3$Site, Year=lpj.g.3$Year, fit=fit.sum, intercept=fit.int, fit.spline=fit.spline, co2=fit.co2, temp=fit.temp, precip=fit.precip)
# par(mfrow=c(3,1))
# plot(gam.lpj.g.3)
summary(fit.df)

ggplot(data=fit.df[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=intercept), color="gray50", size=1.5) +
	geom_line(aes(x=Year, y=co2), color="green4", size=1.0) +
	geom_line(aes(x=Year, y=temp), color="red2", size=1.0) +
	geom_line(aes(x=Year, y=precip), color="blue", size=1.0)

fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

factor.weights <- data.frame(Year = lpj.g.3$Year, Site=lpj.g.3$Site, fit=fit.sum, intercept=fit.int/fit.spline2, fit.spline=fit.spline, co2=fit.co2/fit.spline2, temp=fit.temp/fit.spline2, precip=fit.precip/fit.spline2)
summary(factor.weights)


for(i in 1:nrow(factor.weights)){
	fweight <- abs(factor.weights[i,c("intercept", "co2", "temp", "precip")])
	factor.weights[i,"max"] <- max(fweight)
	factor.weights[i,"factor.max"] <- names(factor.weights[,3:6])[which(fweight==max(fweight))]
}
factor.weights$factor.max <- as.factor(factor.weights$factor.max)
summary(factor.weights)


pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_Centennial.pdf"), width=11, height=8.5)
ggplot(data=factor.weights) + facet_wrap(~Site) +
	geom_line(data=lpj.g.3, aes(x=Year, y=AGB.100), color="gray50", size=2) +
	geom_line(data=factor.weights[factor.weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHA","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHO","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PUN","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PBL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PDL","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PMB","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Centennial") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------------------------------
