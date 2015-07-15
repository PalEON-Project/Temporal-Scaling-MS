# ----------------------------------------
# Temporal Scaling Analyses
# Changes in Strength of Interactions with Temporal Scale
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: How important are interactions between Temp, Precip, & CO2 and how do they 
#                    vary with temporal scale?
# Rationale: The GAMMs get us non-linear response functions that allow temoral variation in what's driving
#			 patterns in models. However, getting driver interactions in these models is difficult at best.
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# (Fit Bayesian Linear Mixed Model per model)
# 1) Temporal Grain (Resolution)
#    -- Fit LMM over constant window with different degrees of smoothing (1 yr - 250 yr)
# 2) Temporal Extent 
#	 -- Fit LMM to different windows (30 yrs 1990-2010, 100 yrs 1910-2010, full window)
# ** Response variables of interest: NPP, possibly dAGB (AGB 1st difference)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# 1) Multi-Model Scale Comparison
#    -- Hypothesis: Interactions among drivers stronger at very fine resolutions (sub-annual, which we can't test)
#                   and large resolutions where combination of factors determines a relative equilibrium point.  
#                   At finer & intermediate grains, the interannual variability is too great for the interactions
#                   to be clear.
#
#    -- Analysis: Run mixed-model anova on posterior distributions of interaction coefficients to determine if 
#                 there are statstically significant trends in temporal scale and differences among levels of
#                 interactions (i.e. 3-way interaction shows stronger shifts than 2-way, etc)
# -------------------------
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
# setwd("~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
path.data <- "Data"
fig.dir <- "Figures"
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed model outputs
#	generated with 1_generate_ecosys.R
load(file.path(path.data, "EcosysData.Rdata"))

# Read in model color scheme
model.colors <- read.csv("raw_inputs/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors
# ----------------------------------------


# # -----------------------
# # Some exploratory Graphing
# # -----------------------
# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	# geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
	# geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
	# scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	# theme_bw()
# # -----------------------
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
		# alpha0[s] ~ dnorm(0, tau2)  # + 
		alpha1[s] ~ dnorm(0, tau2)  # + 
		# alpha2[s] ~ dnorm(0, tau2)  # + 
		# alpha3[s] ~ dnorm(0, tau2)  # + 
		# alpha4[s] ~ dnorm(0, tau2)  # + 
		# alpha5[s] ~ dnorm(0, tau2)  # + 
		# alpha6[s] ~ dnorm(0, tau2)  # + 
		# alpha7[s] ~ dnorm(0, tau2)  # + 
	}
	for(i in 1:n){
#		mu[i] <- beta0[T.SCALE[i]]*alpha0[SITE[i]] + 
		mu[i] <- beta0[T.SCALE[i]] + 
				 beta1[T.SCALE[i]]*TEMP[i] + 
				 beta2[T.SCALE[i]]*PRECIP[i] + 
				 beta3[T.SCALE[i]]*CO2[i] + 
				 beta4[T.SCALE[i]]*TEMP[i]*PRECIP[i] + 
				 beta5[T.SCALE[i]]*TEMP[i]*CO2[i] + 
				 beta6[T.SCALE[i]]*PRECIP[i]*CO2[i] + 
				 beta7[T.SCALE[i]]*TEMP[i]*PRECIP[i]*CO2[i] + 
				 alpha1[SITE[i]] 
		# mu[i] <- beta*TEMP[i]
		y[i] ~  dnorm(mu[i], sigma)
	}
	}

# ----------------------------------------



# ----------------------------------------
# Running the Bayesian Interaction Models 
# ----------------------------------------
# Note: Setting up a loop to go through each model
# ----------------------------------------
data.base="Data/interactions"
fig.base="Figures/interactions"

# Making sure the appropriate file paths exist
if(!dir.exists(data.base)) dir.create(data.base)
if(!dir.exists(fig.base)) dir.create(fig.base)


# Setting up a loop for 1 m.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- unique(ecosys$Model)
model.order   <- unique(ecosys$Model.Order)
# var <- c("NPP", "AGB.diff")
var <- "NPP"
# scale    <- ""
# scales <- c("", ".10", ".50", ".100", ".250")
# # scales <- c(".100")
# t.scales <- ifelse(scales=="", "t.001", paste0("t", scales))
extents <- data.frame(Start=c(850, 1900, 1990), End=c(2010, 2010, 2010)) 

for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]
	out.dir   <- file.path(data.base)
	fig.dir  <- file.path(fig.base)

	# Make sure the appropriate file paths are in place
	if(!dir.exists(out.dir)) dir.create(out.dir)
	if(!dir.exists(fig.dir)) dir.create(fig.dir)

	print(" ")
	print(" ")
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Model: ",m.name, " ------"))

    if(m.name=="jules.stat") var <- "NPP" else var <- c("NPP", "AGB.diff")

for(v in var){
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Variable: ",v, " ------"))

# -----------------------
# Organize the jags inputs
# -----------------------
dat <- ecosys[ecosys$Model==m.name,]
dat <- dat[complete.cases(dat[,v]),]
# summary(dat)
dim(dat)
# ny      <- length(unique(dat$Year))
# ns      <- length(unique(dat$Site))

# y <- T.SCALE <- TEMP <- PRECIP <- CO2 <- vector()
# for(i in 1:length(scales)){
	# T.SCALE <- c(T.SCALE, rep(i, ny*ns)) 
	# y       <- c(y,      dat[,v])	
	# TEMP    <- c(TEMP,   dat[,"tair")])
	# PRECIP  <- c(PRECIP, dat[,"precipf"])
	# CO2     <- c(CO2,    dat[,"CO2"])
# }

# n       <- length(y)
# nt		<- length(unique(T.SCALE))
# SITE    <- rep(as.numeric(dat$Site), nt)
# # MODELS  <- rep(as.numeric(dat$Model), nt)
# nm	    <- length(unique(MODELS))
y       <- dat[,v]
TEMP    <- dat$tair
PRECIP  <- dat$precipf
CO2     <- dat$CO2
SWDOWN  <- dat$swdown
T.SCALE <- as.numeric(dat$Scale)
SITE    <- as.numeric(dat$Site)

n       <- length(y)
nt		<- length(unique(T.SCALE))
ns		<- length(unique(SITE))

params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7","alpha1", "sigma")
# params <- c("beta", "sigma")
# -----------------------

# -----------------------
# Run and save the output
# -----------------------
dat.jags <- list(y=y, n=n, nt=nt, T.SCALE=T.SCALE, TEMP=TEMP, PRECIP=PRECIP, CO2=CO2, ns=ns, SITE=SITE)

jags.out <- jags(data=dat.jags, parameters.to.save=params, n.chains=3, n.iter=100, n.burnin=20, model.file=interactions, DIC=F)

out <- list(Data=dat.jags, jags.out=jags.out)

out.mcmc <- as.mcmc(jags.out)
print(summary(out.mcmc))
# class(summary(out.mcmc[,]))
# names(summary(out.mcmc))
# row.names(summary(out.mcmc)$quantiles)
# # summary(out.mcmc[,which(substr(dimnames(out.mcmc[[1]])[[2]], 1, 2)=="mu")])
# summary(out.mcmc)$statistics
pdf(file.path(fig.dir, paste0("Interactions_TracePlots_", m.order, "_", v, ".pdf")))
plot(out.mcmc)
dev.off()

# Write a data frame with the coefficient 95% ci from the full jags output
coefs.out <- data.frame(Model=m.name, var=row.names(summary(out.mcmc)$quantiles), mean=summary(out.mcmc)$statistics[,"Mean"], CI.025=summary(out.mcmc)$quantiles[,"2.5%"], CI.975=summary(out.mcmc)$quantiles[,"97.5%"])
coefs.out$Scale <- as.factor(ifelse(substr(coefs.out$var,1,4)=="beta", substr(coefs.out$var,7,7), NA))
coefs.out$Scale <- recode(coefs.out$Scale, "'1'='t.001'; '2'='t.010'; '3'='t.050'; '4'='t.100'; '5'='t.250'")
coefs.out$var2 <- as.factor(substr(coefs.out$var, 1, 5))
# summary(coefs.out)

out[["ci.coeff"]] <- coefs.out
# -----------------------

# -----------------------
# Pull from the MCMC iterations
# -----------------------
pulls <- 250
y.predict1 <- array(dim=c(n, pulls))
betas.samp <- data.frame(array(dim=c(n,0)))

for(i in 1:pulls){
	c <- sample(1:length(out.mcmc), 1, replace=T)    # randomly pick a chain
	r <- sample(1:nrow(out.mcmc[[c]]), 1, replace=T) # randomly pick an iteration 

	for(t in 1:nt){
	rows.t   <- which(T.SCALE==t)
	temp.t   <- TEMP[rows.t]
	precip.t <- PRECIP[rows.t]
	co2.t    <- CO2[rows.t]
	alpha    <- out.mcmc[[c]][r,paste0("alpha1[",SITE,"]")[rows.t]] # selecting the appropriate alpha for the vector

	# predicting y from the coefficients (rather than saving mu, 
	# which is very memory & time intensive)
	betas.samp[i,paste0("beta0[",t,"]")] <- out.mcmc[[c]][r,paste0("beta0[",t,"]")]
	betas.samp[i,paste0("beta1[",t,"]")] <- out.mcmc[[c]][r,paste0("beta1[",t,"]")]
	betas.samp[i,paste0("beta2[",t,"]")] <- out.mcmc[[c]][r,paste0("beta2[",t,"]")]
	betas.samp[i,paste0("beta3[",t,"]")] <- out.mcmc[[c]][r,paste0("beta3[",t,"]")]
	betas.samp[i,paste0("beta4[",t,"]")] <- out.mcmc[[c]][r,paste0("beta4[",t,"]")]
	betas.samp[i,paste0("beta5[",t,"]")] <- out.mcmc[[c]][r,paste0("beta5[",t,"]")]
	betas.samp[i,paste0("beta6[",t,"]")] <- out.mcmc[[c]][r,paste0("beta6[",t,"]")]
	betas.samp[i,paste0("beta7[",t,"]")] <- out.mcmc[[c]][r,paste0("beta7[",t,"]")]

	y.predict1[rows.t,i] <- out.mcmc[[c]][r,paste0("beta0[",t,"]")] + 
					        out.mcmc[[c]][r,paste0("beta1[",t,"]")]*temp.t*precip.t*co2.t +
					        out.mcmc[[c]][r,paste0("beta2[",t,"]")]*temp.t* precip.t + 
					        out.mcmc[[c]][r,paste0("beta3[",t,"]")]*temp.t*co2.t + 
					        out.mcmc[[c]][r,paste0("beta4[",t,"]")]*precip.t*co2.t + 
					        out.mcmc[[c]][r,paste0("beta5[",t,"]")]*temp.t + 
					        out.mcmc[[c]][r,paste0("beta6[",t,"]")]*precip.t + 
					        out.mcmc[[c]][r,paste0("beta7[",t,"]")]*co2.t +
					        alpha
	}
}

out.analy <- data.frame(Model=m.name, Scale=recode(T.SCALE, "'1'='t.001'; '2'='t.010'; '3'='t.050'; '4'='t.100'; '5'='t.250'"), Site=rep(dat$Site, nt), Year=rep(dat$Year, nt), Response=y, CO2=CO2, Temp=TEMP, Precip=PRECIP, Pred.y=apply(y.predict1, 1, mean, na.rm=T), Pred.LB=apply(y.predict1, 1, quantile, 0.025, na.rm=T), Pred.UB=apply(y.predict1, 1, quantile, 0.975, na.rm=T))

betas.samp2 <- stack(betas.samp)[,c(2,1)]
names(betas.samp2) <- c("beta.name", "value")
betas.samp2$Model <- as.factor(m.name)
betas.samp2$Beta <- as.factor(substr(betas.samp2$beta.name, 1, 5))
betas.samp2$Scale <- as.factor(recode(substr(betas.samp2$beta.name,7,7), "'1'='t.001'; '2'='t.010'; '3'='t.050'; '4'='t.100'; '5'='t.250'"))
betas.samp2$Interaction <- betas.samp2$Beta
levels(betas.samp2$Interaction) <- recode(levels(betas.samp2$Interaction), "'beta0'='Intercept'; 'beta1'='TEMP'; 'beta2'='PRECIP'; 'beta3'='CO2'; 'beta4'='TEMP x PRECIP'; 'beta5'='TEMP x CO2'; 'beta6'='PRECIP x CO2'; 'beta7'='TEMP x PRECIP x CO2'")
summary(betas.samp2)

out[["betas.sample"]] <- betas.samp2
out[["predicted"]] <- out.analy
# summary(out.analy)

col.model <- model.colors[model.colors$Model.Order %in% unique(dat$Model.Order),"color"]

pdf(file.path(fig.dir, paste0("Interactions_TimeFit_", m.order, "_", v, ".pdf")))
print(
ggplot(data=out.analy) + facet_grid(Site ~ Scale) +
	geom_line(aes(x=Year, y=Response), siz=2)+
	geom_ribbon(aes(x=Year, ymin=Pred.LB, ymax=Pred.UB), alpha=0.5, fill=col.model) +
	geom_line(aes(x=Year, y=Pred.y), color=col.model) + 
	theme_bw()
)
dev.off()

pdf(file.path(fig.dir, paste0("Interactions_BetasScale_", m.order, "_", v, ".pdf")))
print(
ggplot(betas.samp2[,]) + facet_grid(Interaction~., scales="free") +
	geom_hline(aes(yintercept=0), color="black") +
	geom_violin(aes(x=Scale, y=value, fill=Scale, color=Scale), alpha=0.3) +
	theme_bw() +
	theme(strip.text=element_text(size=rel(0.5)))
)
dev.off()

# summary(eval(parse(text=(paste("mcmc", m.name, v, sep=".")))))
# assign(paste("mcmc", m.name, v, sep="."), out)
save(out, file=file.path(data.base, paste0("Interactions_", m.order, "_", v, ".RData")))
# save(dat.jags, m.lpj.g, file=file.path(path.data, "Interactions_LPJ-GUESS.RData"))
# -----------------------

} # end var
} # end model




