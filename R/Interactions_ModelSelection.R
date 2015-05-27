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
setwd("~/Dropbox/PalEON CR/paleon_mip_site")
inputs <- "phase1a_output_variables"
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
	beta0 ~ dunif(0,100)
	beta1 ~ dunif(0,100)
	beta2 ~ dunif(0,100)
	beta3 ~ dunif(0,100)
	beta4 ~ dunif(0,100)
	sigma ~ dunif(0,100)
	# Model
	for(i in 1:n){
		# mu[i] <- beta*TEMP[i]
		mu[i] <- beta0 + beta1*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2*TEMP[i]*PRECIP[i] + beta3*TEMP[i]*CO2[i] + beta4*PRECIP[i]*CO2[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
}

interactions2 <- function(){
	# Priors
	# for(i in 1:5) beta[i] ~ dunif(0,100)
	beta0 ~ dunif(0,100)
	beta1 ~ dunif(0,100)
	beta2 ~ dunif(0,100)
	beta3 ~ dunif(0,100)
	beta4 ~ dunif(0,100)
	sigma ~ dunif(0,100)

	# Model
	for(i in 1:n){
		# mu[i] <- beta*TEMP[i]
		mu[i] <- lag[i] * (beta1*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2*TEMP[i]*PRECIP[i] + beta3*TEMP[i]*CO2[i] + beta4*PRECIP[i]*CO2[i])
		y[i] ~  dnorm(mu[i], sigma)
	}
}


interactions3 <- function(){
	# Priors
	# for(i in 1:5) beta[i] ~ dunif(0,100)
	beta0 ~ dunif(0,100)
	beta1 ~ dunif(0,100)
	beta2 ~ dunif(0,100)
	beta3 ~ dunif(0,100)
	beta4 ~ dunif(0,100)
	sigma ~ dunif(0,100)

	# Model
	for(i in 1:n){
		# mu[i] <- beta*TEMP[i]
		mu[i] <- beta0 * (beta1*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2*TEMP[i]*PRECIP[i] + beta3*TEMP[i]*CO2[i] + beta4*PRECIP[i]*CO2[i])
		y[i] ~  dnorm(mu[i], sigma)
	}
}

interactions4 <- function(){
	# Priors
	# for(i in 1:5) beta[i] ~ dunif(0,100)
	beta0 ~ dunif(0,100)
	beta1 ~ dunif(0,100)
	beta2 ~ dunif(0,100)
	beta3 ~ dunif(0,100)
	beta4 ~ dunif(0,100)
	beta5 ~ dunif(0,100)
	beta6 ~ dunif(0,100)
	beta7 ~ dunif(0,100)
	sigma ~ dunif(0,100)
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
# Run Jags -- Annual Test
# ----------------------------------------

# -----------------------
# Centennial Model 1: AGB in Year y = climate from past century + intercept (mean)
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.100
PRECIP <- lpj.g.pha2$Precip.100
CO2    <- lpj.g.pha2$CO2.100
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "sigma")
# params <- c("beta", "sigma")

out1 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=5000, model.file=interactions, DIC=F)

out1b <- as.mcmc(out1)

par(mfrow=c(length(params)/2, 2))
traceplot(out1b)

y.predict1 <- array(dim=c(n, nrow(out1b[[1]])))
for(i in 1:nrow(out1b[[1]])){
	y.predict1[,i] <- out1b[[1]][i,"beta0"] + out1b[[1]][i,"beta1"]*TEMP*PRECIP*CO2+
					  out1b[[1]][i,"beta2"]*TEMP*PRECIP + 
					  out1b[[1]][i,"beta3"]*TEMP*CO2 + 
					  out1b[[1]][i,"beta4"]*PRECIP*CO2 }

y.pred1 <- apply(y.predict1, 1, mean)
y.lb1 <- apply(y.predict1, 1, quantile, 0.025)
y.ub11 <- apply(y.predict1, 1, quantile, 0.975)

summary(y.pred)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g.pha, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g.pha, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g.pha, type="l", col="green3") 

colNames <- dimnames(out1b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out1b[[1]][,"beta1"]), out1b[[2]][,"beta1"],out1b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------

# -----------------------
# Centennial Model 2: centennial smoothed AGB in Year y = climate from past century + intercept (mean)
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB.100
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.100
PRECIP <- lpj.g.pha2$Precip.100
CO2    <- lpj.g.pha2$CO2.100
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "sigma")
# params <- c("beta", "sigma")

out2 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=5000, model.file=interactions, DIC=F)

out2b <- as.mcmc(out2)

par(mfrow=c(length(params)/2, 2))
traceplot(out2b)

y.predict <- array(dim=c(n, nrow(out2b[[1]])))
for(i in 1:nrow(out2b[[1]])){
	y.predict[,i] <- c((out2b[[1]][i,"beta0"]+out2b[[1]][i,"beta1"]*TEMP*PRECIP*CO2+out2b[[1]][i,"beta2"]*TEMP*PRECIP+out2b[[1]][i,"beta3"]*TEMP*CO2+out2b[[1]][i,"beta4"]*PRECIP*CO2))
}

y.pred <- apply(y.predict, 1, mean)
y.lb <- apply(y.predict, 1, quantile, 0.025)
y.ub <- apply(y.predict, 1, quantile, 0.975)

summary(y.pred)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g.pha, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g.pha, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g.pha, type="l", col="green3") 

colNames <- dimnames(out2b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out2b[[1]][,"beta1"]), out2b[[2]][,"beta1"],out2b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------

# -----------------------
# Centennial Model 3: AGB in Year y = (climate from past century) + (AGB y-1) 
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
lag    <- vector()
  for(i in 1:n) lag[i] <- ifelse(i==1, y[1], y[i-1])

TEMP   <- lpj.g.pha2$Temp.100
PRECIP <- lpj.g.pha2$Precip.100
CO2    <- lpj.g.pha2$CO2.100
params <- c("beta1", "beta2", "beta3", "beta4", "sigma")
# params <- c("beta", "sigma")

out3 <- jags(data=list(y=y, n=n, lag=lag, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=5000, model.file=interactions2, DIC=F)

out3b <- as.mcmc(out3)
summary(out3b)
par(mfrow=c(round((length(params)+.5)/2, 0), 2))
traceplot(out3b)

y.predict <- array(dim=c(n, nrow(out3b[[1]])))
for(i in 1:nrow(out3b[[1]])){
	y.predict[,i] <- lag*(out3b[[1]][i,"beta1"]*TEMP*PRECIP*CO2+out3b[[1]][i,"beta2"]*TEMP*PRECIP+out3b[[1]][i,"beta3"]*TEMP*CO2+out3b[[1]][i,"beta4"]*PRECIP*CO2)
}

y.pred <- apply(y.predict, 1, mean)
y.lb <- apply(y.predict, 1, quantile, 0.025)
y.ub <- apply(y.predict, 1, quantile, 0.975)


summary(y.pred)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


plot(AGB ~ Year, data=lpj.g.pha2[1:100,], type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="red")
	lines(y.lb[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g.pha, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g.pha, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g.pha, type="l", col="green3") 

colNames <- dimnames(out3b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out3b[[1]][,"beta1"]), out3b[[2]][,"beta1"],out3b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------

# -----------------------
# Centennial Model 4: AGB in Year y = (climate from past century) * Constant 
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp.100
PRECIP <- lpj.g.pha2$Precip.100
CO2    <- lpj.g.pha2$CO2.100
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "sigma")
# params <- c("beta", "sigma")

out4 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=5000, model.file=interactions3, DIC=F)

out4b <- as.mcmc(out4)
summary(out4b)
par(mfrow=c(round((length(params)+.5)/2, 0), 2))
traceplot(out4b)

y.predict <- array(dim=c(n, nrow(out4b[[1]])))
for(i in 1:nrow(out4b[[1]])){
	y.predict[,i] <- out4b[[1]][i,"beta0"]*(out4b[[1]][i,"beta1"]*TEMP*PRECIP*CO2+out4b[[1]][i,"beta2"]*TEMP*PRECIP+out4b[[1]][i,"beta3"]*TEMP*CO2+out4b[[1]][i,"beta4"]*PRECIP*CO2)
}

y.pred <- apply(y.predict, 1, mean)
y.lb <- apply(y.predict, 1, quantile, 0.025)
y.ub <- apply(y.predict, 1, quantile, 0.975)


summary(y.pred)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


plot(AGB ~ Year, data=lpj.g.pha2[1:100,], type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="red")
	lines(y.lb[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g.pha, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g.pha, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g.pha, type="l", col="green3") 

colNames <- dimnames(out4b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out4b[[1]][,"beta1"]), out4b[[2]][,"beta1"],out4b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------

# -----------------------
# Centennial Model 5: AGB in Year y = (climate from past century) + Constant 
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp
PRECIP <- lpj.g.pha2$Precip
CO2    <- lpj.g.pha2$CO2
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out5 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=5000, model.file=interactions4, DIC=F)

out5b <- as.mcmc(out5)
summary(out5b)
par(mfrow=c(round((length(params)+.5)/2, 0), 2))
traceplot(out5b)

pulls <- 5000

y.predict5 <- array(dim=c(n, pulls))
for(i in 1:pulls){
	c <- sample(1:length(out5b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(out5b[[c]]), 1, replace=T) # randomly pick an iteration 

	y.predict5[,r] <- out5b[[c]][r,"beta0"] + 
					  out5b[[c]][r,"beta1"]*TEMP*PRECIP*CO2+
					  out5b[[c]][r,"beta2"]*TEMP*PRECIP + 
					  out5b[[c]][r,"beta3"]*TEMP*CO2 + 
					  out5b[[c]][r,"beta4"]*PRECIP*CO2 + 
					  out5b[[c]][r,"beta5"]*TEMP + 
					  out5b[[c]][r,"beta6"]*PRECIP + 
					  out5b[[c]][r,"beta7"]*CO2
}

y.pred5 <- apply(y.predict5, 1, mean, na.rm=T)
y.lb5 <- apply(y.predict5, 1, quantile, 0.025, na.rm=T)
y.ub5 <- apply(y.predict5, 1, quantile, 0.975, na.rm=T)


summary(y.pred5)
par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")

plot(AGB ~ Year, data=lpj.g.pha2[1:100,], type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="red")
	lines(y.lb5[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")



lm5 <- lm(lpj.g.pha2$AGB ~ y.pred5)
summary(lm5)
ci5 <- data.frame(predict(lm5, newdata=lpj.g.pha2, interval="confidence"))
summary(ci5)
plot(lpj.g.pha2$AGB ~ y.pred5)
	lines(ci5$fit, col="red")
	lines(ci5$lwr, lty="dashed", col="blue")
	lines(ci5$upr, lty="dashed", col="blue")

plot(AGB ~ Year, data=lpj.g.pha2[1:100,], type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="red")
	lines(y.lb[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub[1:100] ~ lpj.g.pha2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g.pha, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g.pha, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g.pha, type="l", col="green3") 

colNames <- dimnames(out5b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out5b[[1]][,"beta1"]), out5b[[2]][,"beta1"],out5b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------
# ----------------------------------------
