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
# Gradually building up to an ARIMA model
# This is following the book "Time Series Analysis and Its Applications (With R Examples)
#    by Shumway & Stoffer (2011)
# ----------------------------------------
summary(lpj.g.pha)

# -----------------------
# Chapter 2
# -----------------------
# Just looking at the linear trend through time
lm1 <- lm(AGB ~ time(AGB), data=lpj.g.pha)
summary(lm1)
plot(lpj.g.pha$AGB, type="o")
	abline(lm1)

# Looking at drivers
pairs(lpj.g.pha[,c("AGB", "Temp", "Precip", "CO2")])

# doing an AOV
lm2 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=lpj.g.pha)
summary(lm2)
aov(lm2)
acf(lpj.g.pha$AGB, 500, main="LPJ-GUESS AGB Autocorr")
acf(ed2.pha$AGB, 500, main="ED2 AGB Autocorr")
acf(resid(lm2), 500, main="LPJ-Guess Resids")
acf(diff(lpj.g.pha$AGB), 500, main="LPJ-Guess First Difference")

# Looking at how some differencing influences things
pairs(lpj.g.pha[,c("AGB", "Temp", "Precip", "CO2")], main="Raw")
pairs(cbind(AGB=diff(lpj.g.pha$AGB), Temp=diff(lpj.g.pha$Temp), Precip=diff(lpj.g.pha$Precip), CO2=diff(lpj.g.pha$CO2)), main="First Differences, 1-year")
pairs(cbind(AGB=diff(lpj.g.pha$AGB, lag=10), Temp=diff(lpj.g.pha$Temp, lag=10), Precip=diff(lpj.g.pha$Precip, lag=10), CO2=diff(lpj.g.pha$CO2, lag=10)), main="First Differences, 10-year")
pairs(cbind(AGB=diff(lpj.g.pha$AGB, lag=100), Temp=diff(lpj.g.pha$Temp, lag=100), Precip=diff(lpj.g.pha$Precip, lag=100), CO2=diff(lpj.g.pha$CO2, lag=100)), main="First Differences, 100-year")
pairs(cbind(AGB=diff(lpj.g.pha$AGB, lag=250), Temp=diff(lpj.g.pha$Temp, lag=250), Precip=diff(lpj.g.pha$Precip, lag=250), CO2=diff(lpj.g.pha$CO2, lag=250)), main="First Differences, 250-year")
pairs(cbind(AGB=diff(lpj.g.pha$AGB, lag=500), Temp=diff(lpj.g.pha$Temp, lag=500), Precip=diff(lpj.g.pha$Precip, lag=500), CO2=diff(lpj.g.pha$CO2, lag=500)), main="First Differences, 500-year")

df.diff1 <- data.frame(AGB=diff(lpj.g.pha$AGB), Temp=diff(lpj.g.pha$Temp), Precip=diff(lpj.g.pha$Precip), CO2=diff(lpj.g.pha$CO2)) 
df.diff2 <- data.frame(AGB=diff(lpj.g.pha$AGB, lag=10), Temp=diff(lpj.g.pha$Temp, lag=10), Precip=diff(lpj.g.pha$Precip, lag=10), CO2=diff(lpj.g.pha$CO2, lag=10))
df.diff3 <- data.frame(AGB=diff(lpj.g.pha$AGB, lag=50), Temp=diff(lpj.g.pha$Temp, lag=50), Precip=diff(lpj.g.pha$Precip, lag=50), CO2=diff(lpj.g.pha$CO2, lag=50))
df.diff4 <- data.frame(AGB=diff(lpj.g.pha$AGB, lag=100), Temp=diff(lpj.g.pha$Temp, lag=100), Precip=diff(lpj.g.pha$Precip, lag=100), CO2=diff(lpj.g.pha$CO2, lag=100))
df.diff5 <- data.frame(AGB=diff(lpj.g.pha$AGB, lag=250), Temp=diff(lpj.g.pha$Temp, lag=250), Precip=diff(lpj.g.pha$Precip, lag=250), CO2=diff(lpj.g.pha$CO2, lag=250))
df.diff6 <- data.frame(AGB=diff(lpj.g.pha$AGB, lag=500), Temp=diff(lpj.g.pha$Temp, lag=500), Precip=diff(lpj.g.pha$Precip, lag=500), CO2=diff(lpj.g.pha$CO2, lag=500))

lm3.1 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff1) #   1-yr diffs
summary(lm3.1)

lm3.2 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff2) #  10-yr diffs
summary(lm3.2)
lm3.3 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff3) #  50-yr diffs
summary(lm3.3)
lm3.4 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff4) # 100-yr diffs 
summary(lm3.4)
lm3.5 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff5) # 250-yr diffs
summary(lm3.5)
lm3.6 <- lm(AGB ~ Temp + Precip + CO2 + I(CO2^2) + I(Temp^2), data=df.diff6) # 500-yr diffs
summary(lm3.6)


# Looking at the data with a loess smothing spline
plot(AGB ~ CO2, data=lpj.g.pha, type="o")
	lines(lowess(lpj.g.pha$CO2, lpj.g.pha$AGB), col="red", lwd=2)

plot(AGB ~ Year, data=lpj.g.pha, type="o")
# -----------------------

# -----------------------
# Trying to use interacting Time x Effect smoothers to get non-stationary climate effects
# -----------------------
library(mgcv)

gamm1 <- gamm(AGB ~ s(Year), data=lpj.g.pha)
plot(gamm1$gam)
summary(gamm1$gam)
summary(gamm1$lme)


gamm2 <- gamm(AGB ~ s(CO2), data=lpj.g.pha)
plot(gamm2$gam)
summary(gamm2$gam)
summary(gamm2$lme)

gamm3 <- gamm(AGB ~ s(Temp) + s(Precip) + s(CO2), data=lpj.g.pha)
par(mfrow=c(3,1))
plot(gamm3$gam)
summary(gamm3$gam)
summary(gamm3$lme)

gamm4 <- gamm(AGB ~ s(Temp) + s(Precip) + s(CO2) + s(Year), data=lpj.g.pha)
par(mfrow=c(4,1))
plot(gamm4$gam)
summary(gamm4$gam)
summary(gamm4$lme)

gamm5 <- gam(AGB ~ s(Year, by=CO2) + s(Year, by=Temp) + s(Year, by=Precip), data=lpj.g.pha)
summary(gamm5)

par(mfrow=c(4,1))
#plot(gamm5$gam)
plot(gamm5)
acf(gamm5$resid)

summary(gamm5$gam)

test1 <- predict(gamm5, newdata=lpj.g.pha, se.fit=T)
summary(test1)
summary(test1$fit)
summary(test1$se.fit)

test2 <- data.frame(predict(gamm5, newdata=lpj.g.pha, type="terms"))
summary(test2)
dim(test2)

test2$Pred.full <- test2[,1]*lpj.g.pha$Year + test2[,2]*lpj.g.pha$Year + test2[,3]*lpj.g.pha$Year
summary(test2)

test3 <- data.frame(predict(gamm5, newdata=lpj.g.pha, type="iterms"))
summary(test3)

test4 <- predict(gamm5, newdata=lpj.g.pha, type="lpmatrix")
summary(test4)


# -----------
# Trying to get the non-linear parameter estimates & variance
# Working from examples in predict.gam (?predict.gam)
# -----------
# Create the prediction matrix
Xp <- predict(gamm5, newdata=lpj.g.pha, type="lpmatrix")
summary(Xp)

# just a quick example of how Xp times the coeff gets the predicted values 
fit5 <- Xp %*% coef(gamm5)
summary(fit5)

d=10 # num. terms per effect
fit.int <- Xp[,1] * coef(gamm5)[1]
fit.co2 <- Xp[,(1+1*d+1-d):(1*d+1)] %*% coef(gamm5)[(1+1*d+1-d):(1*d+1)]
fit.temp <- Xp[,(1+2*d+1-d):(2*d+1)] %*% coef(gamm5)[(1+2*d+1-d):(2*d+1)]
fit.precip <- Xp[,(1+3*d+1-d):(3*d+1)] %*% coef(gamm5)[(1+3*d+1-d):(3*d+1)]

fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
fit.spline <- fit.co2 + fit.temp + fit.precip
summary(fit.sum)
summary(fit5 - fit.sum)


fit.df <- data.frame(fit=fit.sum, intercept=fit.int, fit.spline=fit.spline, co2=fit.co2, temp=fit.temp, precip=fit.precip)
par(mfrow=c(3,1))
plot(gamm5)

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


plot(fit.sum, type="l", lwd=2, ylim=c(2,max(fit.sum)))
	lines(abs(factor.weights$intercept)+3, col="gray50")
	lines(abs(factor.weights$co2)+3, col="green3")
	lines(abs(factor.weights$temp)+3, col="red2")
	lines(abs(factor.weights$precip)+3, col="blue")

# Now trying to get the color ramp
ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$temp)/factor.weights$max, color="red") +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$co2)/factor.weights$max*.5, color="green3") +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$intercept)/factor.weights$max, color="black") +	
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$precip), color="blue") +
	theme_bw()


ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$co2), color="green4", size=2) +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$temp)*.33, color="red", size=2) +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$intercept), color="black", size=2) +	
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$precip), color="blue", size=2) +
	theme_bw()

ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$temp), color="red", size=2) +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$co2)*.33, color="green4", size=2) +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$intercept), color="black", size=2) +	
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=abs(factor.weights$precip), color="blue", size=2) +
	theme_bw()

ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), color=rgb(abs(factor.weights$temp),abs(factor.weights$co2),abs(factor.weights$precip)), size=4) +
	labs(y="Year", x="AGB kg/m2", title="Non-Stationary Drivers of AGB") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.75)), axis.text.y=element_text(color="black", size=rel(1.75)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.title=element_text(face="bold", size=rel(3)))

	
ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=0.5, color="red") +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=0.5, color="green3")
ggplot() +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=0.5, color="green3") +
	geom_line(data=factor.weights, aes(x=Year, y=fit), alpha=0.5, color="red")


plot(fit.sum, type="l", lwd=2, ylim=range(factor.weights+fit.sum))
	lines(fit.sum + factor.weights$intercept, col="gray50")
	lines(fit.sum + factor.weights$co2, col="green3")
	lines(fit.sum + factor.weights$temp, col="red2")
	lines(fit.sum + factor.weights$precip, col="blue")


test <- coef(gamm5) %*% Xp

# Calculate variance of sum of predictions
a <- rep(1, 1161)
Xs <- t(a) %*% Xp
summary(Xs)
var.sum <- Xs %*% gamm5$Vp %*% t(Xs)
var.sum # This is the variance of sum of predictions

# -----------

# Getting variance on non-linear funciton of prediction by simulating the posterior distrib of params
# Function from ?predict.gamm
rmvn <- function(n, mu, sig) {
	L <- mroot(sig); m <- ncol(L)
	t(mu + L %*% matrix(rnorm(m*n),m,n))
}

# Replicate n param. vectors
n=1000
br <- rmvn(n, coef(gamm5), gamm5$Vp)
dim(br)
summary(br) 

res <- rep(0,n)
res <- colSums(log(abs(Xp %*% t(br))))
# # This is an expanded version of what is done above
# for(i in 1:n) {
	# pr <- Xp %*% br[i,] # replicate predictions
	# res[i] <- sum(log(abs(pr))) # example non-linear function	
# }
summary(res)
mean(res); var(res)

# -----------
# Differentiating smooths in the model (with CIs for derivates)
summary(gamm5)

# Evaluate derivatives of the smooths & associated SE by finite differencing
X0 <- predict(gamm5, lpj.g.pha, type="lpmatrix") # predict the lp matrix
eps <- 1e-7 # finite difference inteveral

# Shifting the data by eps
newdat <- lpj.g.pha[,c("Year", "CO2", "Temp", "Precip")]
newdat$CO2 <- lpj.g.pha$CO2 + eps
newdat$Temp<- lpj.g.pha$Temp + eps
newdat$Precip <- lpj.g.pha$Precip + eps
X1 <- predict(gamm5, newdat, type="lpmatrix")
summary(X1)

Xp <- (X1-X0)/eps
summary(Xp)

# Xi %*% coef(b) == smooth deriv i
# plot dervis & corresponding CIs
d <- 10 # d = number of smooths
par(mfrow=c(3,1))
for(i in 1:3) {
	Xi <- Xp*0
	Xi[,(i-1)*d+1:d+1] <- Xp[,(i-1)*d+1:d+1] # Xi %*% coef(b) == smooth deriv i
	df <- Xi %*% coef(gamm5) # ith smooth derivative
	df.sd <- rowSums(Xi %*% gamm5$Vp*Xi)^.5 # simple way to do diag(Xi %*% b$Vp %*% t(Xi))^0.5
	plot(newdat$Year, df, type="l", ylim=range(c(df+2*df.sd, df-2*df.sd)), ylab=names(newdat)[i+1])
	lines(newdat$Year, df+2*df.sd, lty=2); lines(newdat$Year, df-2*df.sd, lty=2)
}

# ----------------------------------------

