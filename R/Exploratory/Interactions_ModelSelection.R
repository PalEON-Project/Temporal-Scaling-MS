
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
# ecosys <- ecosys[ecosys$Site=="PHA",]
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


ecosys.pha <- ecosys[ecosys$Site=="PHA",]
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

# Interactions only
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

# Adding a lag
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


# Adding individual effects 
interactions4 <- function(){
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
	sigma ~ dnorm(0, 1.0E-6)
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

# ADDING A SITE Factor
# Site should be a RANDOM effect
# Then will add Model as a FIXED effect
interactions5 <- function(){
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
	sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau

	# General structure:
	# for(m in model){ # Fixed
	# for(s in site){ # Random
	# for(y in year){ # Individual
		# mu[m,s,y] <- beta0[m] + 
				 # beta1[m]*TEMP[s]*PRECIP[s]*CO2[i] + 
				 # beta2[m]*TEMP[s]*PRECIP[s] + beta3[m]*TEMP[i]*CO2[i] + beta4[m]*PRECIP[i]*CO2[i] + 
				 # beta5[m]*TEMP[s] + beta6[m]*PRECIP[i] + beta7[m]*CO2[i]

	# }
	# }
	# }
	tau.s ~ dgamma(0.001, 0.001)
	for(s in 1:ns){	
		alpha.0[s] ~ dnorm(0, tau.s)
		# alpha.1[s] ~ dnorm(0, tau.s)
		# alpha.2[s] ~ dnorm(0, tau.s)
		# alpha.3[s] ~ dnorm(0, tau.s)
		# alpha.4[s] ~ dnorm(0, tau.s)
		# alpha.5[s] ~ dnorm(0, tau.s)
		# alpha.6[s] ~ dnorm(0, tau.s)
		# alpha.7[s] ~ dnorm(0, tau.s)
		}


	# Model
	for(s in 1:ns){
		# b0[s] <- beta0*alpha0[s]
		# b1[s] <- beta1*alpha1[s]
		# b2[s] <- beta2*alpha2[s]
		# b3[s] <- beta3*alpha2[s]
		# b4[s] <- beta4*alpha4[s]
		# b5[s] <- beta5*alpha5[s]
		# b6[s] <- beta6*alpha6[s]
		# b7[s] <- beta7*alpha7[s]
	for(i in 1:n){
	# Adding site random effect
		mu[s,i] <- beta0 + 
				 beta1*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2*TEMP[i]*PRECIP[i] + beta3*TEMP[i]*CO2[s,i] + beta4*PRECIP[s,i]*CO2[s,i] + 
				 beta5*TEMP[i] + beta6*PRECIP[s,i] + beta7*CO2[s,i]

		# mu[i] <- beta*TEMP[i]
		y[s,i] ~  dnorm(mu[s,i], sigma)
	}
	}
}


# Trying for a fixed model effect (instead of random site)
interactions6 <- function(){
	# Priors
	for(m in 1:nm){	
		beta0[m] ~ dnorm(0, 1.0E-6)
		beta1[m] ~ dnorm(0, 1.0E-6)
		beta2[m] ~ dnorm(0, 1.0E-6)
		beta3[m] ~ dnorm(0, 1.0E-6)
		beta4[m] ~ dnorm(0, 1.0E-6)
		beta5[m] ~ dnorm(0, 1.0E-6)
		beta6[m] ~ dnorm(0, 1.0E-6)
		beta7[m] ~ dnorm(0, 1.0E-6)
	}
		# beta0 ~ dnorm(0, 1.0E-6)
		# beta1 ~ dnorm(0, 1.0E-6)
		# beta2 ~ dnorm(0, 1.0E-6)
		# beta3 ~ dnorm(0, 1.0E-6)
		# beta4 ~ dnorm(0, 1.0E-6)
		# beta5 ~ dnorm(0, 1.0E-6)
		# beta6 ~ dnorm(0, 1.0E-6)
		# beta7 ~ dnorm(0, 1.0E-6)
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	# for(m in 1:nm){		
		# alpha[m] ~ dnorm(beta0[m], tau2)  # + 
	# }
	for(i in 1:n){
		mu[i] <- beta0[MODELS[i]] + beta7[MODELS[i]]*CO2[i] +
				 beta1[MODELS[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[MODELS[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[MODELS[i]]*TEMP[i]*CO2[i] + 
				 beta4[MODELS[i]]*PRECIP[i]*CO2[i] + 
				 beta5[MODELS[i]]*TEMP[i] + 
				 beta6[MODELS[i]]*PRECIP[i] + 
				 beta7[MODELS[i]]*CO2[i]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	
}

# Adding a random SITE Interceipt in
interactions7 <- function(){
	# Priors
	for(m in 1:nm){	
		beta0[m] ~ dnorm(0, 1.0E-3)
		beta1[m] ~ dnorm(0, 1.0E-3)
		beta2[m] ~ dnorm(0, 1.0E-3)
		beta3[m] ~ dnorm(0, 1.0E-3)
		beta4[m] ~ dnorm(0, 1.0E-3)
		beta5[m] ~ dnorm(0, 1.0E-3)
		beta6[m] ~ dnorm(0, 1.0E-3)
		beta7[m] ~ dnorm(0, 1.0E-3)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- beta0[MODELS[i]] + beta7[MODELS[i]]*CO2[i] +
				 beta1[MODELS[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[MODELS[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[MODELS[i]]*TEMP[i]*CO2[i] + 
				 beta4[MODELS[i]]*PRECIP[i]*CO2[i] + 
				 beta5[MODELS[i]]*TEMP[i] + 
				 beta6[MODELS[i]]*PRECIP[i] + 
				 beta7[MODELS[i]]*CO2[i] + alpha[SITE[i]]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# Random site slope on whole model
interactions8 <- function(){
	# Priors
	for(m in 1:nm){	
		beta0[m] ~ dnorm(0, 1.0E-3)
		beta1[m] ~ dnorm(0, 1.0E-3)
		beta2[m] ~ dnorm(0, 1.0E-3)
		beta3[m] ~ dnorm(0, 1.0E-3)
		beta4[m] ~ dnorm(0, 1.0E-3)
		beta5[m] ~ dnorm(0, 1.0E-3)
		beta6[m] ~ dnorm(0, 1.0E-3)
		beta7[m] ~ dnorm(0, 1.0E-3)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- (beta0[MODELS[i]] + beta7[MODELS[i]]*CO2[i] +
				 beta1[MODELS[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[MODELS[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[MODELS[i]]*TEMP[i]*CO2[i] + 
				 beta4[MODELS[i]]*PRECIP[i]*CO2[i] + 
				 beta5[MODELS[i]]*TEMP[i] + 
				 beta6[MODELS[i]]*PRECIP[i] + 
				 beta7[MODELS[i]]*CO2[i])*alpha[SITE[i]]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	
}
# ----------------------------------------




summary(ecosys)
library(nlme)
lm1 <- lme(AGB ~ Temp*Precip*CO2*Model -1, random=list(Site=~1), data=ecosys)


# -----------------------
# Model 3: AGB in Year y = (climate) + Constant 
# -----------------------
lpj.g.pha2 <- lpj.g.pha[complete.cases(lpj.g.pha$AGB),]
y      <- lpj.g.pha2$AGB
n      <- nrow(lpj.g.pha2)
TEMP   <- lpj.g.pha2$Temp
PRECIP <- lpj.g.pha2$Precip
CO2    <- lpj.g.pha2$CO2
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out5 <- jags(data=list(y=y, n=n, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions4, DIC=F)

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

lm5 <- lm(y.pred5 ~ lpj.g.pha2$AGB)
summary(lm5)

par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")

plot(AGB ~ Year, data=lpj.g.pha2, xlim=c(1850, 2010), type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	lines(y.lb5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")



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


# -----------------------
# Model 3: AGB in Year y = (climate) + Constant 
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g$AGB),]
y      <- lpj.g2$AGB
ny     <- length(unique(lpj.g2$Year))
n      <- nrow(lpj.g2)
TEMP   <- lpj.g2$Temp
PRECIP <- lpj.g2$Precip
CO2    <- lpj.g2$CO2
ns     <- length(unique(lpj.g2$Site))
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out5 <- jags(data=list(y=y, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ny=ny, ns=ns), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions5, DIC=F)

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

lm5 <- lm(y.pred5 ~ lpj.g2$AGB)
summary(lm5)

par(mfrow=c(1,1))
plot(AGB ~ Year, data=lpj.g2, type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5 ~ lpj.g2$Year, type="l", lwd=1, col="red")
	lines(y.lb5 ~ lpj.g2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5 ~ lpj.g2$Year, type="l", lwd=1, col="blue", lty="dashed")

plot(AGB ~ Year, data=lpj.g2, xlim=c(1850, 2010), type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred5 ~ lpj.g2$Year, type="l", lwd=1, col="red")
	lines(y.lb5 ~ lpj.g2$Year, type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub5 ~ lpj.g2$Year, type="l", lwd=1, col="blue", lty="dashed")



lm5 <- lm(lpj.g2$AGB ~ y.pred5)
summary(lm5)
ci5 <- data.frame(predict(lm5, newdata=lpj.g2, interval="confidence"))
summary(ci5)
plot(lpj.g2$AGB ~ y.pred5)
	lines(ci5$fit, col="red")
	lines(ci5$lwr, lty="dashed", col="blue")
	lines(ci5$upr, lty="dashed", col="blue")

plot(AGB ~ Year, data=lpj.g2[1:100,], type="l", lwd=2, ylab="AGB", xlab="Year")
	lines(y.pred[1:100] ~ lpj.g2$Year[1:100], type="l", lwd=1, col="red")
	lines(y.lb[1:100] ~ lpj.g2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")
	lines(y.ub[1:100] ~ lpj.g2$Year[1:100], type="l", lwd=1, col="blue", lty="dashed")

sec2yr <- 1*60*60*24*365
plot(AGB.dev.100 ~ Year, data=lpj.g, type="l", lwd=2, ylab="AGB", xlab="Year", ylim=c(-2,1)) 
	lines(Temp.abs.dev.100+1 ~ Year, data=lpj.g, type="l", col="red") 
	lines(Precip.abs.dev.100*sec2yr*.01+1 ~ Year, data=ed2.pha, type="l", col="blue") 
	lines(CO2.abs.dev.100*.01 ~ Year, data= lpj.g, type="l", col="green3") 

colNames <- dimnames(out5b[[2]])[[2]]
col.beta1 <- grep("beta1", colNames)
beta1 <- c(as.vector(out5b[[1]][,"beta1"]), out5b[[2]][,"beta1"],out5b[[3]][,"beta1"])

plot(y ~ ed2.pha$Year, type="l", col="red", ylab="AGB", xlab="Year")
# -----------------------
# ----------------------------------------

# -----------------------
# Model 6: Fixed model effect
# -----------------------
ecosys.pha2 <- ecosys.pha[complete.cases(ecosys.pha),]
ecosys.pha2$Model <- droplevels(ecosys.pha2$Model)
# ecosys.pha2 <- ecosys.pha2[ecosys.pha2$Model=="lpj.guess",]
y      <- ecosys.pha2$AGB
ny     <- length(unique(ecosys.pha2$Year))
n      <- nrow(ecosys.pha2)
TEMP   <- ecosys.pha2$Temp
PRECIP <- ecosys.pha2$Precip
CO2    <- ecosys.pha2$CO2
ns     <- length(unique(ecosys.pha2$Site))
MODELS <- as.numeric(ecosys.pha2$Model)
nm	   <- length(unique(MODELS))
params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out6 <- jags(data=list(y=y,n=n, nm=nm, MODELS=MODELS, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions6, DIC=F)


# -----------------------
# Model 7: Fixed model effect, Random Site Intercept
# -----------------------
ecosys2 <- ecosys[complete.cases(ecosys) & ecosys$Year>=1900,]
ecosys2$Model <- droplevels(ecosys2$Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
y      <- ecosys2$AGB
ny     <- length(unique(ecosys2$Year))
n      <- nrow(ecosys2)
TEMP   <- ecosys2$Temp
PRECIP <- ecosys2$Precip
CO2    <- ecosys2$CO2
ns     <- length(unique(ecosys2$Site))
MODELS <- as.numeric(ecosys2$Model)
nm	   <- length(unique(MODELS))
SITE   <- as.numeric(ecosys$Site)

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out7 <- jags(data=list(y=y,n=n, nm=nm, MODELS=MODELS, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions7, DIC=T)
out7$BUGSoutput$DIC	
out7b <- as.mcmc(out7)

# -----------------------
# Model 7: Fixed model effect, Random Site Slope
# -----------------------
ecosys2 <- ecosys[complete.cases(ecosys) & ecosys$Year>=1900,]
ecosys2$Model <- droplevels(ecosys2$Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
y      <- ecosys2$AGB
ny     <- length(unique(ecosys2$Year))
n      <- nrow(ecosys2)
TEMP   <- ecosys2$Temp
PRECIP <- ecosys2$Precip
CO2    <- ecosys2$CO2
ns     <- length(unique(ecosys2$Site))
MODELS <- as.numeric(ecosys2$Model)
nm	   <- length(unique(MODELS))
SITE   <- as.numeric(ecosys$Site)

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out8 <- jags(data=list(y=y,n=n, nm=nm, MODELS=MODELS, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions8, DIC=T)
out8$BUGSoutput$DIC	
out8b <- as.mcmc(out8)
