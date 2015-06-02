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
inputs <- "phase1a_output_variables"
path.data <- "~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
# ----------------------------------------

# # Note: Commented out because saved as EcosysData.RData 1 June 2015
# #       (with an increasing number of models, running this every time became cumbersome)
# # ----------------------------------------
# # Load Data Sets
# # ----------------------------------------
# # Ecosystem Model Outputs
# ecosys <- read.csv(file.path(inputs, "MIP_Data_Ann_2015.csv"))
# ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
# levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
# summary(ecosys)

# # CO2 Record
# nc.co2 <- nc_open("~/Desktop/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
# co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
# nc_close(nc.co2)

# # Merging CO2 into Model Outputs
# ecosys <- merge(ecosys, co2.ann)
# summary(ecosys)

# # Colors used for graphing
# model.colors <- read.csv("~/Desktop/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
# model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
# levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
# model.colors

# model.colors <- model.colors[order(model.colors$Model.Order),]
# model.colors
# # ----------------------------------------

# # ----------------------------------------
# # Calculate Deviations from desired reference point
# # Reference Point: 0850-0869 (spinup climate)
# # ----------------------------------------
# vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp")
# vars.climate <- c("Temp", "Precip", "CO2")

# # vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

# ref.window <- 850:869

# for(s in unique(ecosys$Site)){
	# for(m in unique(ecosys$Model)){
		# # -----------------------
		# # Model Variabiles -- Relative Change
		# # Deviation = percent above or below the mean for the reference window
		# #             (observed-ref.mean)/ref.mean 
		# # -----------------------
		# for(v in unique(vars)){
			# ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			# ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".dev")] <- (ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean)/ref.mean

		# }
		# # -----------------------

		# # -----------------------
		# # Climate Drivers -- Absolute, not relative change
		# # Deviation = absolute deviation from reference window
		# #             observed - ref.mean
		# # -----------------------
		# for(v in unique(vars.climate)){
			# ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			# ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".abs.dev")] <- ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean
		# }
		# # -----------------------
	# }
# }

# summary(ecosys)
# # ----------------------------------------


# # ----------------------------------------
# # Perform Temporal Smoothing on Data
# # Note: Smoothing is performed over the PREVIOUS 100 years becuase ecosystems 
# #       cannot respond to what they have not yet experienced
# # ----------------------------------------
# vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp", "Temp", "Precip", "CO2")
# vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

# for(s in unique(ecosys$Site)){
	# for(m in unique(ecosys$Model)){
		# # -----------------------
		# # 10-yr Smoothing
		# # -----------------------
		# ## Non-standardized
		# for(v in vars){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		# }

		# ## Non-standardized
		# for(v in vars.dev){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		# }
		# # -----------------------

		# # -----------------------
		# # 50-yr Smoothing
		# # -----------------------
		# ## Non-standardized
		# for(v in vars){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollmean(temp, k=50, align="right", fill=NA)
		# }

		# ## Non-standardized
		# for(v in vars.dev){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollmean(temp, k=50, align="right", fill=NA)
		# }
		# # -----------------------

		# # -----------------------
		# # 100-yr Smoothing
		# # -----------------------
		# ## Non-standardized
		# for(v in vars){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		# }

		# ## Non-standardized
		# for(v in vars.dev){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		# }
		# # -----------------------

		# # -----------------------
		# # 250-Year Smoothing
		# # -----------------------
		# ## Non-standardized
		# for(v in vars){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollmean(temp, k=250, align="right", fill=NA)
		# }

		# ## Non-standardized
		# for(v in vars.dev){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollmean(temp, k=250, align="right", fill=NA)
		# }
		# # -----------------------

	# }
# }
# summary(ecosys)
# save(ecosys, model.colors, file=file.path(path.data, "EcosysData.Rdata"))

# ----------------------------------------
# Load processing from previous step
load(file.path(path.data, "EcosysData.Rdata"))
# ----------------------------------------


# -----------------------
# Some exploratory Graphing
# -----------------------
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
	geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	theme_bw()
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



# ----------------------------------------
# Running the Bayesian Interaction Models 
# ----------------------------------------
# Note: The models are currently set up to be run individually by model
# ----------------------------------------

# -----------------------
# LPJ-GUESS
# -----------------------
dat.lpj.g <- ecosys[ecosys$Model=="lpj.guess",]
dat.lpj.g <- dat.lpj.g[complete.cases(dat.lpj.g$AGB.250),]
summary(dat.lpj.g)

ny      <- length(unique(dat.lpj.g$Year))
ns      <- length(unique(dat.lpj.g$Site))
y       <- c(dat.lpj.g$AGB, dat.lpj.g$AGB.50, dat.lpj.g$AGB.100, dat.lpj.g$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.lpj.g$Temp, dat.lpj.g$Temp.50, dat.lpj.g$Temp.100, dat.lpj.g$Temp.250)
PRECIP  <- c(dat.lpj.g$Precip, dat.lpj.g$Precip.50, dat.lpj.g$Precip.100, dat.lpj.g$Precip.250)
CO2     <- c(dat.lpj.g$CO2, dat.lpj.g$CO2.50, dat.lpj.g$CO2.100, dat.lpj.g$CO2.250)
SITE    <- rep(as.numeric(dat.lpj.g$Site), nt)
MODELS  <- rep(as.numeric(dat.lpj.g$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.lpj.g.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.lpj.g <- jags(data=dat.lpj.g.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.lpj.g.jags, m.lpj.g, file=filepath(path.data, "Interactions_LPJ-GUESS.RData"))
# -----------------------


# -----------------------
# LPJ-WSL
# -----------------------
dat.lpj.w <- ecosys[ecosys$Model=="lpj.wsl",]
dat.lpj.w <- dat.lpj.w[complete.cases(dat.lpj.w$AGB.250),]
summary(dat.lpj.w)

ny      <- length(unique(dat.lpj.w$Year))
ns      <- length(unique(dat.lpj.w$Site))
y       <- c(dat.lpj.w$AGB, dat.lpj.w$AGB.50, dat.lpj.w$AGB.100, dat.lpj.w$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.lpj.w$Temp, dat.lpj.w$Temp.50, dat.lpj.w$Temp.100, dat.lpj.w$Temp.250)
PRECIP  <- c(dat.lpj.w$Precip, dat.lpj.w$Precip.50, dat.lpj.w$Precip.100, dat.lpj.w$Precip.250)
CO2     <- c(dat.lpj.w$CO2, dat.lpj.w$CO2.50, dat.lpj.w$CO2.100, dat.lpj.w$CO2.250)
SITE    <- rep(as.numeric(dat.lpj.w$Site), nt)
MODELS  <- rep(as.numeric(dat.lpj.w$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.lpj.w.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.lpj.w <- jags(data=dat.lpj.w.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.lpj.w.jags, m.lpj.w, file=filepath(path.data, "Interactions_LPJ-WSL.RData"))
# -----------------------


# -----------------------
# JULES-STATIC
# -----------------------
dat.jules.s <- ecosys[ecosys$Model=="jules.stat",]
dat.jules.s <- dat.jules.s[complete.cases(dat.jules.s$AGB.250),]
summary(dat.jules.s)

ny      <- length(unique(dat.jules.s$Year))
ns      <- length(unique(dat.jules.s$Site))
y       <- c(dat.jules.s$AGB, dat.jules.s$AGB.50, dat.jules.s$AGB.100, dat.jules.s$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.jules.s$Temp, dat.jules.s$Temp.50, dat.jules.s$Temp.100, dat.jules.s$Temp.250)
PRECIP  <- c(dat.jules.s$Precip, dat.jules.s$Precip.50, dat.jules.s$Precip.100, dat.jules.s$Precip.250)
CO2     <- c(dat.jules.s$CO2, dat.jules.s$CO2.50, dat.jules.s$CO2.100, dat.jules.s$CO2.250)
SITE    <- rep(as.numeric(dat.jules.s$Site), nt)
MODELS  <- rep(as.numeric(dat.jules.s$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.jules.s.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.jules.s <- jags(data=dat.jules.s.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.jules.s.jags, m.jules.s, file=filepath(path.data, "Interactions_JULES-STATIC.RData"))
# -----------------------


# -----------------------
# JULES-TRIFFID
# -----------------------
dat.jules.triff <- ecosys[ecosys$Model=="jules.triffid",]
dat.jules.triff <- dat.jules.triff[complete.cases(dat.jules.triff$AGB.250),]
summary(dat.jules.triff)

ny      <- length(unique(dat.jules.triff$Year))
ns      <- length(unique(dat.jules.triff$Site))
y       <- c(dat.jules.triff$AGB, dat.jules.triff$AGB.50, dat.jules.triff$AGB.100, dat.jules.triff$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.jules.triff$Temp, dat.jules.triff$Temp.50, dat.jules.triff$Temp.100, dat.jules.triff$Temp.250)
PRECIP  <- c(dat.jules.triff$Precip, dat.jules.triff$Precip.50, dat.jules.triff$Precip.100, dat.jules.triff$Precip.250)
CO2     <- c(dat.jules.triff$CO2, dat.jules.triff$CO2.50, dat.jules.triff$CO2.100, dat.jules.triff$CO2.250)
SITE    <- rep(as.numeric(dat.jules.triff$Site), nt)
MODELS  <- rep(as.numeric(dat.jules.triff$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.jules.triff.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.jules.triff <- jags(data=dat.jules.triff.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.jules.triff.jags, m.jules.triff, file=filepath(path.data, "Interactions_JULES-TRIFFID.RData"))
# -----------------------


# -----------------------
# LINKAGES
# -----------------------
dat.linkages <- ecosys[ecosys$Model=="linkages",]
dat.linkages <- dat.linkages[complete.cases(dat.linkages$AGB.250),]
summary(dat.linkages)

ny      <- length(unique(dat.linkages$Year))
ns      <- length(unique(dat.linkages$Site))
y       <- c(dat.linkages$AGB, dat.linkages$AGB.50, dat.linkages$AGB.100, dat.linkages$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.linkages$Temp, dat.linkages$Temp.50, dat.linkages$Temp.100, dat.linkages$Temp.250)
PRECIP  <- c(dat.linkages$Precip, dat.linkages$Precip.50, dat.linkages$Precip.100, dat.linkages$Precip.250)
CO2     <- c(dat.linkages$CO2, dat.linkages$CO2.50, dat.linkages$CO2.100, dat.linkages$CO2.250)
SITE    <- rep(as.numeric(dat.linkages$Site), nt)
MODELS  <- rep(as.numeric(dat.linkages$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.linkages.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.linkages <- jags(data=dat.linkages.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.linkages.jags, m.linkages, file=filepath(path.data, "Interactions_LINKAGES.RData"))
# -----------------------


# -----------------------
# SiBCASA
# -----------------------
dat.sibcasa <- ecosys[ecosys$Model=="sibcasa",]
dat.sibcasa <- dat.sibcasa[complete.cases(dat.sibcasa$AGB.250),]
summary(dat.sibcasa)

ny      <- length(unique(dat.sibcasa$Year))
ns      <- length(unique(dat.sibcasa$Site))
y       <- c(dat.sibcasa$AGB, dat.sibcasa$AGB.50, dat.sibcasa$AGB.100, dat.sibcasa$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.sibcasa$Temp, dat.sibcasa$Temp.50, dat.sibcasa$Temp.100, dat.sibcasa$Temp.250)
PRECIP  <- c(dat.sibcasa$Precip, dat.sibcasa$Precip.50, dat.sibcasa$Precip.100, dat.sibcasa$Precip.250)
CO2     <- c(dat.sibcasa$CO2, dat.sibcasa$CO2.50, dat.sibcasa$CO2.100, dat.sibcasa$CO2.250)
SITE    <- rep(as.numeric(dat.sibcasa$Site), nt)
MODELS  <- rep(as.numeric(dat.sibcasa$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.sibcasa.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.sibcasa <- jags(data=dat.sibcasa.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.sibcasa.jags, m.sibcasa, file=filepath(path.data, "Interactions_SiBCASA.RData"))
# -----------------------


# -----------------------
# ED
# -----------------------
dat.ed <- ecosys[ecosys$Model=="ed2",]
dat.ed <- dat.ed[complete.cases(dat.ed$AGB.250),]
summary(dat.ed)

ny      <- length(unique(dat.ed$Year))
ns      <- length(unique(dat.ed$Site))
y       <- c(dat.ed$AGB, dat.ed$AGB.50, dat.ed$AGB.100, dat.ed$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.ed$Temp, dat.ed$Temp.50, dat.ed$Temp.100, dat.ed$Temp.250)
PRECIP  <- c(dat.ed$Precip, dat.ed$Precip.50, dat.ed$Precip.100, dat.ed$Precip.250)
CO2     <- c(dat.ed$CO2, dat.ed$CO2.50, dat.ed$CO2.100, dat.ed$CO2.250)
SITE    <- rep(as.numeric(dat.ed$Site), nt)
MODELS  <- rep(as.numeric(dat.ed$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.ed.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.ed <- jags(data=dat.ed.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.ed.jags, m.ed, file=filepath(path.data, "Interactions_ED.RData"))
# -----------------------

# -----------------------
# ED-LU
# -----------------------
dat.ed.lu <- ecosys[ecosys$Model=="ed2.lu",]
dat.ed.lu <- dat.ed.lu[complete.cases(dat.ed.lu$AGB.250),]
summary(dat.ed.lu)

ny      <- length(unique(dat.ed.lu$Year))
ns      <- length(unique(dat.ed.lu$Site))
y       <- c(dat.ed.lu$AGB, dat.ed.lu$AGB.50, dat.ed.lu$AGB.100, dat.ed.lu$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.ed.lu$Temp, dat.ed.lu$Temp.50, dat.ed.lu$Temp.100, dat.ed.lu$Temp.250)
PRECIP  <- c(dat.ed.lu$Precip, dat.ed.lu$Precip.50, dat.ed.lu$Precip.100, dat.ed.lu$Precip.250)
CO2     <- c(dat.ed.lu$CO2, dat.ed.lu$CO2.50, dat.ed.lu$CO2.100, dat.ed.lu$CO2.250)
SITE    <- rep(as.numeric(dat.ed.lu$Site), nt)
MODELS  <- rep(as.numeric(dat.ed.lu$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.ed.lu.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.ed.lu <- jags(data=dat.ed.lu.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.ed.lu.jags, m.ed.lu, file=filepath(path.data, "Interactions_ED-LU.RData"))
# -----------------------


# -----------------------
# CLM-CN
# -----------------------
dat.clm.cn <- ecosys[ecosys$Model=="clm.cn",]
dat.clm.cn <- dat.clm.cn[complete.cases(dat.clm.cn$AGB.250),]
summary(dat.clm.cn)

ny      <- length(unique(dat.clm.cn$Year))
ns      <- length(unique(dat.clm.cn$Site))
y       <- c(dat.clm.cn$AGB, dat.clm.cn$AGB.50, dat.clm.cn$AGB.100, dat.clm.cn$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.clm.cn$Temp, dat.clm.cn$Temp.50, dat.clm.cn$Temp.100, dat.clm.cn$Temp.250)
PRECIP  <- c(dat.clm.cn$Precip, dat.clm.cn$Precip.50, dat.clm.cn$Precip.100, dat.clm.cn$Precip.250)
CO2     <- c(dat.clm.cn$CO2, dat.clm.cn$CO2.50, dat.clm.cn$CO2.100, dat.clm.cn$CO2.250)
SITE    <- rep(as.numeric(dat.clm.cn$Site), nt)
MODELS  <- rep(as.numeric(dat.clm.cn$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.clm.cn.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.clm.cn <- jags(data=dat.clm.cn.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.clm.cn.jags, m.clm.cn, file=filepath(path.data, "Interactions_CLM-CN.RData"))
# -----------------------


# -----------------------
# CLM-BGC
# -----------------------
dat.clm.bgc <- ecosys[ecosys$Model=="clm.bgc",]
dat.clm.bgc <- dat.clm.bgc[complete.cases(dat.clm.bgc$AGB.250),]
summary(dat.clm.bgc)

ny      <- length(unique(dat.clm.bgc$Year))
ns      <- length(unique(dat.clm.bgc$Site))
y       <- c(dat.clm.bgc$AGB, dat.clm.bgc$AGB.50, dat.clm.bgc$AGB.100, dat.clm.bgc$AGB.250)
T.SCALE <- c(rep(1, ny*ns), rep(2, ny*ns),  rep(3, ny*ns),  rep(4, ny*ns))
n       <- length(y)
nt		<- length(unique(T.SCALE))
TEMP    <- c(dat.clm.bgc$Temp, dat.clm.bgc$Temp.50, dat.clm.bgc$Temp.100, dat.clm.bgc$Temp.250)
PRECIP  <- c(dat.clm.bgc$Precip, dat.clm.bgc$Precip.50, dat.clm.bgc$Precip.100, dat.clm.bgc$Precip.250)
CO2     <- c(dat.clm.bgc$CO2, dat.clm.bgc$CO2.50, dat.clm.bgc$CO2.100, dat.clm.bgc$CO2.250)
SITE    <- rep(as.numeric(dat.clm.bgc$Site), nt)
MODELS  <- rep(as.numeric(dat.clm.bgc$Model), nt)
nm	    <- length(unique(MODELS))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "sigma", "mu")
# params <- c("beta", "sigma")

dat.clm.bgc.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE),

m.clm.bgc <- jags(data=dat.clm.bgc.jags, parameters.to.save=params, n.chains=3, n.iter=50000, n.burnin=25000, model.file=interactions, DIC=F)

save(dat.clm.bgc.jags, m.clm.bgc, file=filepath(path.data, "Interactions_CLM-BGC.RData"))
# -----------------------

# ----------------------------------------


# ----------------------------------------
# Looking at the Output
# ----------------------------------------

# -----------------------
# LPJ-GUESS
# -----------------------
m.lpj.g2 <- as.mcmc(m.lpj.g)
summary(m.lpj.g2[,which(!substr(dimnames(m.lpj.g2[[1]])[[2]], 1, 2)=="mu")])

plot(m.lpj.g2[,which(!substr(dimnames(m.lpj.g2[[1]])[[2]], 1, 2)=="mu")])

pulls <- 5000
y.predict1 <- array(dim=c(n, pulls))

for(i in 1:pulls){
	c <- sample(1:length(m.lpj.g2), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(m.lpj.g2[[c]]), 1, replace=T) # randomly pick an iteration 

	for(t in 1:nt){
	rows.t   <- which(T.SCALE==t)
	temp.t   <- TEMP[rows.t]
	precip.t <- PRECIP[rows.t]
	co2.t    <- CO2[rows.t]

	y.predict1[rows.t,i] <- m.lpj.g2[[c]][r,paste0("beta0[",t,"]")] + 
					        m.lpj.g2[[c]][r,paste0("beta1[",t,"]")]*temp.t*precip.t*co2.t +
					        m.lpj.g2[[c]][r,paste0("beta2[",t,"]")]*temp.t* precip.t + 
					        m.lpj.g2[[c]][r,paste0("beta3[",t,"]")]*temp.t*co2.t + 
					        m.lpj.g2[[c]][r,paste0("beta4[",t,"]")]*precip.t*co2.t + 
					        m.lpj.g2[[c]][r,paste0("beta5[",t,"]")]*temp.t + 
					        m.lpj.g2[[c]][r,paste0("beta6[",t,"]")]*precip.t + 
					        m.lpj.g2[[c]][r,paste0("beta7[",t,"]")]*co2.t
	}
}


lpj.g.analy <- data.frame(Scale=recode(T.SCALE, "'1'='t.001'; '2'='t.050'; '3'='t.100'; '4'='t.250'"), Site=rep(lpj.g2$Site, nt), Year=rep(lpj.g2$Year, nt), AGB=y, CO2=CO2, Temp=TEMP, Precip=PRECIP, Pred.y=apply(y.predict1, 1, mean, na.rm=T), Pred.LB=apply(y.predict1, 1, quantile, 0.025, na.rm=T), Pred.UB=apply(y.predict1, 1, quantile, 0.975, na.rm=T))

summary(lpj.g.analy)

ggplot(data=lpj.g.analy) + facet_grid(Site ~ Scale) +
	geom_line(aes(x=Year, y=AGB), siz=2)+
	geom_ribbon(aes(x=Year, ymin=Pred.LB, ymax=Pred.UB), alpha=0.5, fill="red") +
	geom_line(aes(x=Year, y=Pred.y), color="red")

