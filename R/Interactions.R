# ----------------------------------------
# Temporal Scaling Analyses
# Changes in Strength of Interactions with Temporal Scale
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
# Set Directories
# ----------------------------------------
setwd("~/Desktop/PalEON CR/paleon_mip_site")
inputs <- "phase1a_mput_variables"
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
# -----------------------
# ----------------------------------------



# ----------------------------------------
# Model Specifications
# ----------------------------------------
# Model:  y[i] ~  dnorm(mu[i], sigma) # assumes constant variance 
#        mu[i] <- beta[1] + 
#				  beta[2]*TEMP[i]*PRECIP[i]*CO2 + 
#				  beta[3]*TEMP[i]*PRECIP[i] + beta[4]*TEMP[i]*PRECIP[i] + beta[5]*PRECIP[i]*CO2[i]
# 
# Priors: beta[1:5], sigma
#	      sigma ~ dnorm(0,0.001)
#		  for(i in 1:5) beta[i] ~ dnorm(0,0.001)
# ----------------------------------------

interactions <- function(){
	# Priors
	# for(i in 1:5) beta[i] ~ dunif(0,100)
	beta0 ~ dnorm(0, 1.0E-6)
	beta1 ~ dnorm(0, 1.0E-6)
	beta2 ~ dnorm(0, 1.0E-6)
	beta3 ~ dnorm(0, 1.0E-6)
	beta4 ~ dnorm(0, 1.0E-6)
	beta5 ~ dnorm(0, 1.0E-6)
	beta6 ~ dnorm(0, 1.0E-6)
	beta7 ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(.001,.001)
    sigma <- 1/tau
	# Model
	for(i in 1:n){
		# mu[i] <- beta*TEMP[i]
		mu[i] <- beta0 + 
				 beta1*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2*TEMP[i]*PRECIP[i] + beta3*TEMP[i]*CO2[i] + beta4*PRECIP[i]*CO2[i] + 
				 beta5*TEMP[i] + beta6*PRECIP[i] + beta7*CO2[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
}
# ----------------------------------------

# ----------------------------------------
# Run Jags -- Scales Test on LPJ-GUESS at PHA
# ----------------------------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]

# In a non-bayesian way
lm.lpj.g.pha <- lm(AGB ~ Temp*Precip*CO2, data=lpj.g.pha2)
agb.pred <- fitted(lm.lpj.g.pha)
summary(lm.lpj.g.pha)

lm.lpj.g.pha.250 <- lm(AGB.250 ~ Temp.250*Precip.250*CO2.250, data=lpj.g.pha2)
agb.pred.250 <- fitted(lm.lpj.g.pha.250)
summary(lm.lpj.g.pha.250)

par(mfrow=c(2,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2)
	lines(agb.pred ~ lpj.g.pha2$Year, col="red")
plot(AGB.250 ~ Year, data=lpj.g.pha2, type="l", lwd=2)
	lines(agb.pred.250 ~ lpj.g.pha2$Year, col="red")

# -----------------------
# 1-Yr Model, AGB
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp
PRECIP <- lpj.g.pha2$Precip
CO2    <- lpj.g.pha2$CO2
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "tau", "sigma", "mu")
# params <- c("beta", "sigma")


m1 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=5000, n.burnin=1000, model.file=interactions, DIC=F)

m1b <- as.mcmc(m1)
summary(m1b[,which(!substr(dimnames(m1b[[1]])[[2]], 1, 2)=="mu")])

plot(m1b[,which(!substr(dimnames(m1b[[1]])[[2]], 1, 2)=="mu")])


pulls <- 5000
y.predict1 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(m1b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(m1b[[c]]), 1, replace=T) # randomly pick an iteration 

	y.predict1[,r] <- m1b[[c]][r,"beta0"] + 
					  m1b[[c]][r,"beta1"]*TEMP*PRECIP*CO2+
					  m1b[[c]][r,"beta2"]*TEMP*PRECIP + 
					  m1b[[c]][r,"beta3"]*TEMP*CO2 + 
					  m1b[[c]][r,"beta4"]*PRECIP*CO2 + 
					  m1b[[c]][r,"beta5"]*TEMP + 
					  m1b[[c]][r,"beta6"]*PRECIP + 
					  m1b[[c]][r,"beta7"]*CO2
}


y.pred1 <- apply(y.predict1, 1, mean, na.rm=T)
y.lb1 <- apply(y.predict1, 1, quantile, 0.025, na.rm=T)
y.ub1 <- apply(y.predict1, 1, quantile, 0.975, na.rm=T)

summary(y.pred1)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


lm1 <- lm(y.pred1 ~ lpj.g.pha2$AGB )
summary(lm1)
ci1 <- data.frame(predict(lm1, newdata=lpj.g.pha2, interval="confidence"))
summary(ci1)

plot( y.pred1 ~ lpj.g.pha2$AGB)
	lines(ci1$fit ~ lpj.g.pha2$AGB, col="red")
	lines(ci1$lwr ~ lpj.g.pha2$AGB, lty="dashed", col="blue")
	lines(ci1$upr ~ lpj.g.pha2$AGB, lty="dashed", col="blue")
# -----------------------

# -----------------------
# 50-Yr Model, AGB
# -----------------------
# lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB.50
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.50
PRECIP <- lpj.g.pha2$Precip.50
CO2    <- lpj.g.pha2$CO2.50
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "tau", "sigma", "mu")
# params <- c("beta", "sigma")


m50 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=5000, n.burnin=1000, model.file=interactions, DIC=F)

m50b <- as.mcmc(m50)
summary(m50b[,which(!substr(dimnames(m50b[[1]])[[2]], 1, 2)=="mu")])

plot(m50b[,which(!substr(dimnames(m50b[[1]])[[2]], 1, 2)=="mu")])

pulls <- 5000
y.predict50 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(m50b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(m50b[[c]]), 1, replace=T) # randomly pick an iteration 

	y.predict50[,r] <- m50b[[c]][r,"beta0"] + 
					  m50b[[c]][r,"beta1"]*TEMP*PRECIP*CO2+
					  m50b[[c]][r,"beta2"]*TEMP*PRECIP + 
					  m50b[[c]][r,"beta3"]*TEMP*CO2 + 
					  m50b[[c]][r,"beta4"]*PRECIP*CO2 + 
					  m50b[[c]][r,"beta5"]*TEMP + 
					  m50b[[c]][r,"beta6"]*PRECIP + 
					  m50b[[c]][r,"beta7"]*CO2
}


y.pred50 <- apply(y.predict50, 1, mean, na.rm=T)
y.lb50 <- apply(y.predict50, 1, quantile, 0.025, na.rm=T)
y.ub50 <- apply(y.predict50, 1, quantile, 0.975, na.rm=T)

summary(y.pred50)
par(mfrow=c(1,1))
plot(AGB.50 ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred50 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb50 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub50 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


lm50 <- lm(y.pred50 ~ lpj.g.pha2$AGB.50)
summary(lm50)
ci50 <- data.frame(predict(lm50, newdata=lpj.g.pha2, interval="confidence"))
summary(ci50)

plot( y.pred50 ~ lpj.g.pha2$AGB.50)
	lines(ci50$fit ~ lpj.g.pha2$AGB.50, col="red")
	lines(ci50$lwr ~ lpj.g.pha2$AGB.50, lty="dashed", col="blue")
	lines(ci50$upr ~ lpj.g.pha2$AGB.50, lty="dashed", col="blue")
# -----------------------


# -----------------------
# 100-Yr Model, AGB
# -----------------------
# lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB.100
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.100
PRECIP <- lpj.g.pha2$Precip.100
CO2    <- lpj.g.pha2$CO2.100
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "tau", "sigma", "mu")
# params <- c("beta", "sigma")


m100 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=5000, n.burnin=1000, model.file=interactions, DIC=F)

m100b <- as.mcmc(m100)
summary(m100b[,which(!substr(dimnames(m100b[[1]])[[2]], 1, 2)=="mu")])

plot(m100b[,which(!substr(dimnames(m100b[[1]])[[2]], 1, 2)=="mu")])

pulls <- 10000
y.predict100 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(m100b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(m100b[[c]]), 1, replace=T) # randomly pick an iteration 

	y.predict100[,r] <- m100b[[c]][r,"beta0"] + 
					  m100b[[c]][r,"beta1"]*TEMP*PRECIP*CO2+
					  m100b[[c]][r,"beta2"]*TEMP*PRECIP + 
					  m100b[[c]][r,"beta3"]*TEMP*CO2 + 
					  m100b[[c]][r,"beta4"]*PRECIP*CO2 + 
					  m100b[[c]][r,"beta5"]*TEMP + 
					  m100b[[c]][r,"beta6"]*PRECIP + 
					  m100b[[c]][r,"beta7"]*CO2
}


y.pred100 <- apply(y.predict100, 1, mean, na.rm=T)
y.lb100 <- apply(y.predict100, 1, quantile, 0.025, na.rm=T)
y.ub100 <- apply(y.predict100, 1, quantile, 0.975, na.rm=T)

summary(y.pred100)
par(mfrow=c(1,1))
plot(AGB.100 ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred100 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb100 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub100 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


lm100 <- lm(y.pred100 ~ lpj.g.pha2$AGB.100)
summary(lm100)
ci100 <- data.frame(predict(lm100, newdata=lpj.g.pha2, interval="confidence"))
summary(ci100)

plot(y.pred100 ~ lpj.g.pha2$AGB.100)
	lines(ci100$fit ~ lpj.g.pha2$AGB.100, col="red")
	lines(ci100$lwr ~ lpj.g.pha2$AGB.100, lty="dashed", col="blue")
	lines(ci100$upr ~ lpj.g.pha2$AGB.100, lty="dashed", col="blue")
# -----------------------

# -----------------------
# 250-Yr Model, AGB
# -----------------------
# lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB.250
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.250
PRECIP <- lpj.g.pha2$Precip.250
CO2    <- lpj.g.pha2$CO2.250
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "tau", "sigma", "mu")
# params <- c("beta", "sigma")


m250 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=5000, n.burnin=1000, model.file=interactions, DIC=F)

m250b <- as.mcmc(m250)
summary(m250b[,which(!substr(dimnames(m250b[[1]])[[2]], 1, 2)=="mu")])

plot(m250b[,which(!substr(dimnames(m250b[[1]])[[2]], 1, 2)=="mu")])

pulls <- 25000
y.predict250 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(m250b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(m250b[[c]]), 1, replace=T) # randomly pick an iteration 

	y.predict250[,r] <- m250b[[c]][r,"beta0"] + 
					  m250b[[c]][r,"beta1"]*TEMP*PRECIP*CO2+
					  m250b[[c]][r,"beta2"]*TEMP*PRECIP + 
					  m250b[[c]][r,"beta3"]*TEMP*CO2 + 
					  m250b[[c]][r,"beta4"]*PRECIP*CO2 + 
					  m250b[[c]][r,"beta5"]*TEMP + 
					  m250b[[c]][r,"beta6"]*PRECIP + 
					  m250b[[c]][r,"beta7"]*CO2
}


y.pred250 <- apply(y.predict250, 1, mean, na.rm=T)
y.lb250 <- apply(y.predict250, 1, quantile, 0.025, na.rm=T)
y.ub250 <- apply(y.predict250, 1, quantile, 0.975, na.rm=T)

summary(y.pred250)
par(mfrow=c(1,1))
plot(AGB.250 ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred250 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb250 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub250 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")

lm250a <- lm(AGB.250 ~ Temp.250*Precip.250*CO2.250, data=lpj.g.pha2)
summary(lm250a)

lm250 <- lm(y.pred250 ~ lpj.g.pha2$AGB.250)
summary(lm250)
ci250 <- data.frame(predict(lm250, newdata=lpj.g.pha2, interval="confidence"))
summary(ci250)

plot(y.pred250 ~ lpj.g.pha2$AGB.250)
	lines(ci250$fit ~ lpj.g.pha2$AGB.250, col="red")
	lines(ci250$lwr ~ lpj.g.pha2$AGB.250, lty="dashed", col="blue")
	lines(ci250$upr ~ lpj.g.pha2$AGB.250, lty="dashed", col="blue")
# -----------------------
