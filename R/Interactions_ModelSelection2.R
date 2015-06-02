
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
# # Adding a random SITE Interceipt in
# # Original model by MODEL
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

# Retooling a bit so that we run the analysis for a single model at a time, but across all scales 
interactions8 <- function(){
	# Priors
	for(t in 1:nt){	
		beta0[t] ~ dnorm(0, 1.0E-6)
		beta1[t] ~ dnorm(0, 1.0E-6)
		beta2[t] ~ dnorm(0, 1.0E-6)
		beta3[t] ~ dnorm(0, 1.0E-6)
		beta4[t] ~ dnorm(0, 1.0E-6)
		beta5[t] ~ dnorm(0, 1.0E-6)
		beta6[t] ~ dnorm(0, 1.0E-6)
		beta7[t] ~ dnorm(0, 1.0E-6)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- beta0[T.SCALE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[T.SCALE[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[T.SCALE[i]]*TEMP[i]*CO2[i] + 
				 beta4[T.SCALE[i]]*PRECIP[i]*CO2[i] + 
				 beta5[T.SCALE[i]]*TEMP[i] + 
				 beta6[T.SCALE[i]]*PRECIP[i] + 
				 beta7[T.SCALE[i]]*CO2[i] + 
				 alpha[SITE[i]]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# Retooling a bit so that we run the analysis for a single model at a time, but across all scales 
interactions9 <- function(){
	# Priors
	for(t in 1:nt){	
		beta0[t] ~ dnorm(0, 1.0E-6)
		beta1[t] ~ dnorm(0, 1.0E-6)
		beta2[t] ~ dnorm(0, 1.0E-6)
		beta3[t] ~ dnorm(0, 1.0E-6)
		beta4[t] ~ dnorm(0, 1.0E-6)
		beta5[t] ~ dnorm(0, 1.0E-6)
		beta6[t] ~ dnorm(0, 1.0E-6)
		beta7[t] ~ dnorm(0, 1.0E-6)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- (beta0[T.SCALE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[T.SCALE[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[T.SCALE[i]]*TEMP[i]*CO2[i] + 
				 beta4[T.SCALE[i]]*PRECIP[i]*CO2[i] + 
				 beta5[T.SCALE[i]]*TEMP[i] + 
				 beta6[T.SCALE[i]]*PRECIP[i] + 
				 beta7[T.SCALE[i]]*CO2[i])* 
				 alpha[SITE[i]]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# Trying random slopes on all parameters
## NOTE NOTE: NOT VERY STABLE... AT ALL!!
interactions10 <- function(){
	# Priors
	for(t in 1:nt){	
		beta0[t] ~ dnorm(0, 1.0E-2)
		beta1[t] ~ dnorm(0, 1.0E-2)
		beta2[t] ~ dnorm(0, 1.0E-2)
		beta3[t] ~ dnorm(0, 1.0E-2)
		beta4[t] ~ dnorm(0, 1.0E-2)
		beta5[t] ~ dnorm(0, 1.0E-2)
		beta6[t] ~ dnorm(0, 1.0E-2)
		beta7[t] ~ dnorm(0, 1.0E-6)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha0[s] ~ dnorm(0, tau2)  # + 
		alpha1[s] ~ dnorm(0, tau2)  # + 
		alpha2[s] ~ dnorm(0, tau2)  # + 
		alpha3[s] ~ dnorm(0, tau2)  # + 
		alpha4[s] ~ dnorm(0, tau2)  # + 
		alpha5[s] ~ dnorm(0, tau2)  # + 
		alpha6[s] ~ dnorm(0, tau2)  # + 
		alpha7[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- beta0[T.SCALE[i]]*alpha0[SITE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i]*alpha1[SITE[i]]  + 
				 beta2[T.SCALE[i]]*TEMP[i]*PRECIP[i]*alpha2[SITE[i]]  + 
				 beta3[T.SCALE[i]]*TEMP[i]*CO2[i]*alpha3[SITE[i]]  + 
				 beta4[T.SCALE[i]]*PRECIP[i]*CO2[i]*alpha4[SITE[i]]  + 
				 beta5[T.SCALE[i]]*TEMP[i]*alpha5[SITE[i]]  + 
				 beta6[T.SCALE[i]]*PRECIP[i]*alpha6[SITE[i]]  + 
				 beta7[T.SCALE[i]]*CO2[i]*alpha7[SITE[i]]

		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# Trying random slopes on just the intercept
interactions11 <- function(){
	# Priors
	for(t in 1:nt){	
		beta0[t] ~ dnorm(0, 1.0E-6)
		beta1[t] ~ dnorm(0, 1.0E-6)
		beta2[t] ~ dnorm(0, 1.0E-6)
		beta3[t] ~ dnorm(0, 1.0E-6)
		beta4[t] ~ dnorm(0, 1.0E-6)
		beta5[t] ~ dnorm(0, 1.0E-6)
		beta6[t] ~ dnorm(0, 1.0E-6)
		beta7[t] ~ dnorm(0, 1.0E-6)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha0[s] ~ dnorm(0, tau2)  # + 
		# alpha1[s] ~ dnorm(0, tau2)  # + 
		# alpha2[s] ~ dnorm(0, tau2)  # + 
		# alpha3[s] ~ dnorm(0, tau2)  # + 
		# alpha4[s] ~ dnorm(0, tau2)  # + 
		# alpha5[s] ~ dnorm(0, tau2)  # + 
		# alpha6[s] ~ dnorm(0, tau2)  # + 
		# alpha7[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- beta0[T.SCALE[i]]*alpha0[SITE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[T.SCALE[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[T.SCALE[i]]*TEMP[i]*CO2[i] + 
				 beta4[T.SCALE[i]]*PRECIP[i]*CO2[i] + 
				 beta5[T.SCALE[i]]*TEMP[i] + 
				 beta6[T.SCALE[i]]*PRECIP[i] + 
				 beta7[T.SCALE[i]]*CO2[i]
		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# Trying random slopes on just the intercept
interactions12 <- function(){
	# Priors
	for(t in 1:nt){	
		beta0[t] ~ dnorm(0, 1.0E-2)
		beta1[t] ~ dnorm(0, 1.0E-2)
		beta2[t] ~ dnorm(0, 1.0E-2)
		beta3[t] ~ dnorm(0, 1.0E-2)
		beta4[t] ~ dnorm(0, 1.0E-2)
		beta5[t] ~ dnorm(0, 1.0E-2)
		beta6[t] ~ dnorm(0, 1.0E-2)
		beta7[t] ~ dnorm(0, 1.0E-2)
	}
	# sigma ~ dnorm(0, 1.0E-6)
    tau ~ dgamma(0.001,0.001)
    sigma <- 1/tau
	tau2 <- 1/tau
	for(s in 1:ns){		
		alpha0[s] ~ dnorm(0, tau2)  # + 
		alpha1[s] ~ dnorm(0, tau2)  # + 
		# alpha2[s] ~ dnorm(0, tau2)  # + 
		# alpha3[s] ~ dnorm(0, tau2)  # + 
		# alpha4[s] ~ dnorm(0, tau2)  # + 
		# alpha5[s] ~ dnorm(0, tau2)  # + 
		# alpha6[s] ~ dnorm(0, tau2)  # + 
		# alpha7[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
		mu[i] <- beta0[T.SCALE[i]]*alpha0[SITE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 beta2[T.SCALE[i]]*TEMP[i]*PRECIP[i] + 
				 beta3[T.SCALE[i]]*TEMP[i]*CO2[i] + 
				 beta4[T.SCALE[i]]*PRECIP[i]*CO2[i] + 
				 beta5[T.SCALE[i]]*TEMP[i] + 
				 beta6[T.SCALE[i]]*PRECIP[i] + 
				 beta7[T.SCALE[i]]*CO2[i] + 
				 alpha1[SITE[i]] 
		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# ----------------------------------------




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
# Model 8: Fixed TEMPORAL SCALE slopes, Random Site Interceipt
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Year>=1900,]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out8 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=interactions8, DIC=T)
DIC.8 <- out8$BUGSoutput$DIC; DIC.8
out8b <- as.mcmc(out8)
summary(out8b)
plot(out8b)

# -----------------------
# Model 9: Fixed TEMPORAL SCALE slopes, Random Site Slope
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Year>=1900,]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out9 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=10000, n.burnin=5000, model.file=interactions9, DIC=T)
DIC.9 <- out9$BUGSoutput$DIC; DIC.9
out9b <- as.mcmc(out9)
summary(out9b)
plot(out9b)

# -----------------------
# Model 10: Fixed TEMPORAL SCALE slopes, Random Site Slopes on all terms
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Year>=1900,]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out10 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=10000, n.burnin=5000, model.file=interactions10, DIC=T)
DIC.10 <- out10$BUGSoutput$DIC; DIC.10
out10b <- as.mcmc(out10)
summary(out10b)
plot(out10b)

# -----------------------
# Model 11: Fixed TEMPORAL SCALE slopes, Random Site Slopes on just the intercept
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Year>=1900,]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out11 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=5000, n.burnin=2000, model.file=interactions11, DIC=T)
DIC.11 <- out11$BUGSoutput$DIC; DIC.11
out11b <- as.mcmc(out11)
summary(out11b)
plot(out11b)


# -----------------------
# Model 12: Fixed TEMPORAL SCALE slopes, Random Site Slope on the intercept plus a separate random site intercept
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g) & lpj.g$Year>=1900,]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma")
# params <- c("beta", "sigma")

out12 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=10000, n.burnin=5000, model.file=interactions12, DIC=T)
DIC.12 <- out12$BUGSoutput$DIC; DIC.12
out12b <- as.mcmc(out12)
summary(out12b)
plot(out12b)



















# --------------------------------------------------
# Model 12: with Random intercept and random slope on the Fixed-Effect Intercept best, so lets play with that at the moment
# -----------------------
# Model 12: Fixed TEMPORAL SCALE slopes, Random Site Slope on the intercept plus a separate random site intercept
# -----------------------
lpj.g2 <- lpj.g[complete.cases(lpj.g),]
lpj.g2$Model <- droplevels(lpj.g2 $Model)
# ecosys2 <- ecosys2[ecosys2$Model=="lpj.guess",]
ny      <- length(unique(lpj.g2$Year))
ns      <- length(unique(lpj.g2$Site))

y       <- c(lpj.g2$AGB,    lpj.g2$AGB.50, lpj.g2$AGB.100, lpj.g2$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))

TEMP    <- c(lpj.g2$Temp, lpj.g2$Temp.50, lpj.g2$Temp.100, lpj.g2$Temp.250)
PRECIP  <- c(lpj.g2$Precip, lpj.g2$Precip.50, lpj.g2$Precip.100, lpj.g2$Precip.250)
CO2     <- c(lpj.g2$CO2, lpj.g2$CO2.50, lpj.g2$CO2.100, lpj.g2$CO2.250)
SITE    <- rep(as.numeric(lpj.g2$Site), nt)
MODELS  <- rep(as.numeric(lpj.g2$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

out12 <- jags(data=list(y=y,n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE), parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=2500, model.file=interactions12, DIC=T)
DIC.12 <- out12$BUGSoutput$DIC; DIC.12
out12b <- as.mcmc(out12)
summary(out12b[,which(!substr(dimnames(out12b[[1]])[[2]], 1, 2)=="mu")])
plot(out12b[,which(!substr(dimnames(out12b[[1]])[[2]], 1, 2)=="mu")])

pulls <- 500
y.predict1 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(out12b), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(out12b[[c]]), 1, replace=T) # randomly pick an iteration 

	for(t in 1:nt){
	rows.t   <- which(T.SCALE==t)
	temp.t   <- TEMP[rows.t]
	precip.t <- PRECIP[rows.t]
	co2.t    <- CO2[rows.t]

	y.predict1[rows.t,i] <- out12b[[c]][r,paste0("beta0[",t,"]")] + 
					        out12b[[c]][r,paste0("beta1[",t,"]")]*temp.t*precip.t*co2.t +
					        out12b[[c]][r,paste0("beta2[",t,"]")]*temp.t* precip.t + 
					        out12b[[c]][r,paste0("beta3[",t,"]")]*temp.t*co2.t + 
					        out12b[[c]][r,paste0("beta4[",t,"]")]*precip.t*co2.t + 
					        out12b[[c]][r,paste0("beta5[",t,"]")]*temp.t + 
					        out12b[[c]][r,paste0("beta6[",t,"]")]*precip.t + 
					        out12b[[c]][r,paste0("beta7[",t,"]")]*co2.t
	}
}

lpj.g.analy <- data.frame(Scale=recode(T.SCALE, "'1'='t.001'; '2'='t.050'; '3'='t.100'; '4'='t.250'"), Site=rep(lpj.g2$Site, nt), Year=rep(lpj.g2$Year, nt), AGB=y, CO2=CO2, Temp=TEMP, Precip=PRECIP, Pred.y=apply(y.predict1, 1, mean, na.rm=T), Pred.LB=apply(y.predict1, 1, quantile, 0.025, na.rm=T), Pred.UB=apply(y.predict1, 1, quantile, 0.975, na.rm=T))

summary(lpj.g.analy)

# y.pred1 <- apply(y.predict1, 1, mean, na.rm=T)
# y.lb1 <- apply(y.predict1, 1, quantile, 0.025, na.rm=T)
# y.ub1 <- apply(y.predict1, 1, quantile, 0.975, na.rm=T)

summary(y.pred1)
ggplot(data=lpj.g.analy) + facet_grid(Site ~ Scale) +
	geom_line(aes(x=Year, y=AGB), siz=2)+
	geom_ribbon(aes(x=Year, ymin=Pred.LB, ymax=Pred.UB), alpha=0.5, fill="red") +
	geom_line(aes(x=Year, y=Pred.y), color="red")

# par(mfrow=c(1,1))
# plot(AGB ~ Year, data=lpj.g.pha2, type="l", lwd=2, ylab="AGB", xlab="Year")
	# lines(y.pred1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="red")
	# lines(y.lb1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")
	# lines(y.ub1 ~ lpj.g.pha2$Year, type="l", lwd=1, col="blue", lty="dashed")


