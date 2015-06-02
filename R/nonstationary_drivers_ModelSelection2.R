# ----------------------------------------
# Temporal Scaling Analyses
# Non-Stationary Drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 2 June 2015 
# Latest Updates: Tryign to add temporal correlation into the site by model form
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
setwd("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling")
inputs <- "phase1a_output_variables"
fig.dir <- "~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
path.data <- "~/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"

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
# nc.co2 <- nc_open("~/Dropbox/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
# co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
# nc_close(nc.co2)

# # Merging CO2 into Model Outputs
# ecosys <- merge(ecosys, co2.ann)
# summary(ecosys)

# # Colors used for graphing
# model.colors <- read.csv("~/Dropbox/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
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
	geom_line(aes(x=Year, y=AGB.250, color=Model.Order), size=1.5) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	theme_bw()

ggplot(data=ecosys[!ecosys$Model=="linkages",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.7, alpha=0.6) +
	geom_line(aes(x=Year, y=NPP.50, color=Model.Order), size=1.5) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
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
ecosys.pha <- ecosys[ecosys$Site=="PHA" | ecosys$Site=="PHO",]
# -----------------------
# ----------------------------------------



# ----------------------------------------
# Model approach: 
# ----------------------------------------
library(mgcv)

ecosys.pha <- ecosys[ecosys$Site=="PHA" | ecosys$Site=="PHO",]
# Copied from 
data=ecosys.pha; model="lpj.guess"; scale=""; response="NPP"; k=4; outdir=outdir
	# data     = data frame with data.temp in it
	# model    = which model to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata.temp file
	library(mgcv)
	
	# creating a working data.temp frame with just the data.temp we want
	data.temp          <- data[data$Model==model,c("Model", "Updated", "Model.Order", "Site", "Year")]
	data.temp$Scale    <- as.factor(paste0("t", scale))
	data.temp$response <- data[data$Model==model,paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model,paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model,paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model,paste0("CO2", scale)]	
	
data.temp2 <- data.temp[data.temp$Site=="PHA",]	
# data.temp <- data.temp2
    # -----------
	# Running the basic model
	# -----------
	# Running the gamm; note this has no temporal autocorrelation
	gam1 <- gam(response ~ s(Temp, by=Site, k=k) + s(Precip,by=Site, k=k) + s(CO2, by=Site, k=k) , data=data.temp)
	summary(gam1)
	acf(gam1$resid)
	pacf(gam1$resid)
	plot(gam1, residuals=T, pages=1)

	gam1b <- gam(response ~ s(Temp, by=Site, k=k) + s(Precip,by=Site, k=k) + s(CO2, by=Site, k=k) , data=data.temp[data.temp$Site=="PHA",])
	summary(gam1b)
	acf(gam1b$resid)
	pacf(gam1b$resid)
	plot(gam1b, residuals=T, pages=1)

	gam2 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k) , data=data.temp)
	summary(gam2$gam)
	acf(resid(gam2$lme, type="normalized"))
	pacf(resid(gam2$lme, type="normalized"))
	plot(gam2$gam, residuals=T, pages=1)

	gam3 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k), random=list(Site=~1), data=data.temp)
	summary(gam3$gam)
	acf(resid(gam3$lme, type="normalized"))
	pacf(resid(gam3$lme, type="normalized"))
	plot(gam3$gam, residuals=T, pages=1)

	gam4 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k), random=list(Site=~1), data=data.temp, correlation=corARMA(form=~Year, p=1))
	summary(gam4$gam)
	acf(resid(gam4$lme, type="normalized"))
	pacf(resid(gam4$lme, type="normalized"))
	plot(gam4$gam, residuals=T, pages=1)


	gam4b <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), data=data.temp[data.temp$Site=="PHA",], correlation=corARMA(form=~Year, p=1))
	summary(gam4b$gam)
	acf(resid(gam4b$lme, type="normalized"))
	pacf(resid(gam4b$lme, type="normalized"))
	plot(gam4b$gam, residuals=T, pages=1)

	gam5 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k),  random=list(Site=~1), data=data.temp, correlation=corARMA(form=~Year, p=1))
	summary(gam5$gam)
	acf(resid(gam5$lme, type="normalized"))
	pacf(resid(gam5$lme, type="normalized"))
	plot(gam5$gam, residuals=T, pages=1)

	# gam5 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k),  random=list(Site=~1), data=data.temp, correlation=corARMA(form=~Year, p=2))
	# summary(gam5$gam)
	# acf(resid(gam5$lme, type="normalized"))
	# pacf(resid(gam5$lme, type="normalized"))
	# plot(gam5$gam, residuals=T, pages=1)

	gam6 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k) + Site-1,  random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	summary(gam6$gam)
	acf(resid(gam6$lme, type="normalized"))
	pacf(resid(gam6$lme, type="normalized"))
	plot(gam6$gam, residuals=T, pages=1)

	gam7 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k) + Site -1,  random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	summary(gam7$gam)
	acf(resid(gam7$lme, type="normalized"))
	pacf(resid(gam7$lme, type="normalized"))
	plot(gam7$gam, residuals=T, pages=1)

	anova(gam2$lme, gam3$lme, gam4$lme, gam5$lme)
	anova(gam5$lme, gam6$lme)
	anova(gam6$lme, gam7$lme)

	gam8 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k) + s(Year, by=Site) + Site -1, data=data.temp)
	summary(gam8$gam)
	acf(resid(gam8$lme, type="normalized"))
	pacf(resid(gam8$lme, type="normalized"))
	plot(gam8$gam, residuals=T, pages=1)

	gam9 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k) + s(Year, by=Site) + Site -1,  random=list(Site=~Site), data=data.temp)
	summary(gam9$gam)
	acf(resid(gam9$lme, type="normalized"))
	pacf(resid(gam9$lme, type="normalized"))
	plot(gam9$gam, residuals=T, pages=1)

	anova(gam8$lme, gam9$lme)

	gam10 <- gamm(response ~ s(Temp, by=Site, k=k) + s(Precip, by=Site, k=k) + s(CO2, by=Site, k=k) + s(Year, by=Site) + Site -1,  random=list(Site=~1), data=data.temp, correlation=corARMA(form=~Year, p=1))
	summary(gam10$gam)
	acf(resid(gam10$lme, type="normalized"))
	pacf(resid(gam10$lme, type="normalized"))
	plot(gam10$gam, residuals=T, pages=1)

	anova(gam8$lme, gam10$lme, gam9$lme)
	# return(assign(paste0("gam.", model), gam1)) # Saving the gam a name based on the model

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1, newdata= data.temp)
	
	# return(assign(paste0("data.temp.", model), data.temp) )# saving the temporary data.temp frame for use in post-hoc QA/QC
	# -----------

	# -----------
	# calculating the weights for each of the factors
	# -----------
	# Create the prediction matrix
	Xp <- predict(gam7$gam, newdata=data.temp, type="lpmatrix")

	fit <- Xp %*% coef(gam7$gam) # The full predicted values; used for model QA/QC
	coef.gam <- coef(gam7$gam) # the gam coefficients
	
	# Some handy column indices
	cols.site <- which(substr(names(coef.gam),1,4)=="Site")
	cols.temp <- which(substr(names(coef.gam),1,7)=="s(Temp)")
	cols.precip <- which(substr(names(coef.gam),1,9)=="s(Precip)")
	cols.co2 <- which(substr(names(coef.gam),1,6)=="s(CO2)")

	# calculating the smoother for each effect
	fit.int <- Xp[,cols.site] %*% coef.gam[cols.site]
	fit.temp <- Xp[,cols.temp] %*% coef.gam[cols.temp]
	fit.precip <- Xp[,cols.precip] %*% coef.gam[cols.precip]
	fit.co2 <- Xp[,cols.co2] %*% coef.gam[cols.co2]

	# Summing the fixed effects to do QA/QC and throwing a warning if it's not very close to the predicted values
	fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
	fit.spline <- fit.co2 + fit.temp + fit.precip
	if(max(abs(fit - fit.sum))>1e-4) print("***** WARNING: sum of fixed effects not equal to predicted value *****")

	# summing the absolute values to get the weights for each fixed effect
	fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
	fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

	# Factor weights are determined by the relative strength of Temp, Precip, & CO2
	factor.weights <- data.frame(Year = data.temp$Year, Site=data.temp$Site, fit=fit, intercept=fit.int/fit.spline2, fit.spline=fit.spline, co2=fit.co2/fit.spline2, temp=fit.temp/fit.spline2, precip=fit.precip/fit.spline2)

	# doing a little bit of handy-dandy calculation to give a flag as to which factor is given the greatest weight in a given year
	for(i in 1:nrow(factor.weights)){
		fweight <- abs(factor.weights[i,c("co2", "temp", "precip")])
		factor.weights[i,"max"] <- max(fweight)
		factor.weights[i,"factor.max"] <- names(factor.weights[,c("co2", "temp", "precip")])[which(fweight==max(fweight))]
	}
	factor.weights$factor.max <- as.factor(factor.weights$factor.max)
	# summary(factor.weights)
	
ggplot(data=factor.weights) + facet_wrap(~Site) +
	geom_line(data=data.temp, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data= factor.weights[factor.weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHA","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=factor.weights[factor.weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(factor.weights[factor.weights$Site=="PHO","temp"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","co2"]),
			  			abs(factor.weights[factor.weights$Site=="PHO","precip"])), size=4) +
	theme_bw()
	# -----------
