# ----------------------------------------
# Making a funciton out of the gam prediciton script so that it can be farmed out using mclapply
# Christy Rollinson, crollinson@gmail.com
# Date Created: 10 July 2015
# ----------------------------------------

paleon.gams.models <- function(data, k, predictors.all, PFT=F){
	# data       = data frame with data for 1 model, 1 extent & 1 resolution
	# model.name = name of model (goes into data tables to make life easier)
	# Extent     = the temporal extent of the data (e.g. 850-2010, 1990-2010)
	# Resolution = the temporal resolution fo the data (t.001 = yearly, etc)
	# k          = number of knots in gam splines
	# predictors.all   = the full set of predictors that should be in data
	# predictor.suffix = any suffix attached to the predictors; 
	# 						 = yearly mean from drivers; 
	# 						.gs = growing season mean from drivers; 
	#						no suffix = yearly mean from model output


	# ----------------------------------------
	# Load Libaries
	# ----------------------------------------
	require(parallel)
	require(mgcv)
	
	source('R/0_process.gamm.R', chdir = TRUE)
	source('R/0_gamm_Plots.R', chdir = TRUE)
	# ----------------------------------------

	# extracting some info previously specified
	model.name <- unique(data$Model)
	ext.name   <- paste0(min(data$Year), "-", 2010)
    ext.index  <- regexpr("-", ext.name)[1] # Find the index so we can split apart extent into 2 numbers
	extent     <- c(as.numeric(substr(ext.name,1,ext.index-1)), as.numeric(substr(ext.name, ext.index+1,nchar(paste(ext.name)))))
	resolution <- unique(data$Resolution)

	data$Extent <- as.factor(ext.name)

	# ----------------------------------------
	# Select which set of predictors based on which model it is
	# Each of the models is having different stability issues & has different predictors
	# ----------------------------------------
	# Note: different model structure based on whether or not we have random sites
	# ----------------------------------------
	if(PFT==T){
		predictors=c("tair", "precipf", "CO2", "PFT", "Time", "Year")
		gam1 <- gam(Y ~  s(Year, k=3, bs="cr") + s(Time, k=4) + s(tair, k=k, by=PFT) + s(precipf, k=k, by=PFT) + s(CO2, k=k, by=PFT), data=data, correlation=corARMA(form=~Year|PlotID, p=1))
	# ----------------------------------------
	} else {
	# ----------------------------------------
		predictors=c("tair", "precipf", "CO2", "Time", "Year")
		gam1 <- gam(Y ~ s(Year, k=3, bs="cr") + s(Time, k=4) + s(tair, k=k) + s(precipf, k=k) + s(CO2, k=k), data=data)
	}
	# ----------------------------------------

	# ----------------------------------------
	# Storing the predicted values from the gam
	# ----------------------------------------
	data$fit.gam <- predict(gam1, newdata=data)
	# ----------------------------------------

	# ----------------------------------------
	# Run all of the post-processing (calculate CIs, etc)
	# ----------------------------------------
	mod.out <- process.gamm(gamm.model=gam1, data=data, model.name=model.name, extent=extent, resolution=resolution, vars=predictors, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	# ----------------------------------------
	
	return(mod.out)
}

