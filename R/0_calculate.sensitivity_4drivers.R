# ----------------------------------------
# Making a funciton out of the gam prediciton script so that it can be farmed out using mclapply
# Christy Rollinson, crollinson@gmail.com
# Date Created: 10 July 2015
# ----------------------------------------

# -----------------
# Matrix of Models and Drivers
# -----------------
# Var    ED2  ED2-LU  CLM-BGC  CLM-CN  LPJ-WSL  LPJ-GUESS  JULES-STAT  JULES-TRIFFID  SIBCASA  LINKAGES
# tair.gs    X     X        X       X        X         X          X             X           X        X
# precipf.gs X     X        X       X        X         X          X             X           X        X
# swdown.gs  X     X        X       X        X         X          X             X           X
# lwdown.gs  X     X                         X                    X             X           X
# wind.gs    X     X        X       X                             X             X           X
# psurf.gs   X     X        X       X                             X             X           X
# qair.gs    X     X        X       X                             X             X           X
# CO2.gs     X     X        X       X        X         X          X             X           X
# Ndep                   ?       ?


paleon.gams.models <- function(data, response, k, predictors.all, site.effects){
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
	ext.name   <- unique(data$Extent)
    ext.index  <- regexpr("-", ext.name)[1] # Find the index so we can split apart extent into 2 numbers
	extent     <- c(as.numeric(substr(ext.name,1,ext.index-1)), as.numeric(substr(ext.name, ext.index+1,nchar(paste(ext.name)))))
	resolution <- unique(data$Resolution)

	# ----------------------------------------
	# Running the gam; note this now has AR1 temporal autocorrelation
	# ----------------------------------------
	# This is different from model.site.gam in that it has an added random site slope.
	#	This random effect lets us gauge the overall model.name response to our fixed effects 
	#   regardless of the site.  
	#   Pros: Generalized and helps characterize the basic model responses
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)
	#------------------------------

	# ----------------------------------------
	# Select which set of predictors based on which model it is
	# Each of the models is having different stability issues & has different predictors
	# ----------------------------------------
	# Note: different model structure based on whether or not we have random sites
	# ----------------------------------------
	if(site.effects==T){
	if(substr(model.name,1,3)=="lin") {
		predictors <- c("tair", "precipf")
		gam1 <- gam(Y ~ s(tair, k=k) + s(precipf, k=k) + Site -1, data=data, correlation=corARMA(form=~Year|Site, p=1))
	} else {
		predictors <- c("tair", "precipf", "swdown", "CO2")
		gam1 <- gam(Y ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(CO2, k=k) + Site -1, data=data, correlation=corARMA(form=~Year|Site, p=1))
	}

	# ----------------------------------------
	} else {
	# ----------------------------------------
	if(substr(model.name,1,3)=="lin") {
		predictors <- c("tair", "precipf")
		gam1 <- gam(Y ~ s(tair, k=k) + s(precipf, k=k), data=data, correlation=corARMA(form=~Year, p=1))
	} else {
		predictors <- c("tair", "precipf", "swdown", "CO2")
		gam1 <- gam(Y ~ s(tair, k=k) + s(precipf, k=k) + s(swdown, k=k) + s(CO2, k=k) , data=data, correlation=corARMA(form=~Year, p=1))
	}

	}
	# ----------------------------------------

	# ----------------------------------------
	print("-------------------------------------")
	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))
	print(summary(gam1))	

	# get rid of values for predictors not used in the models for clarity later on
	data[,predictors.all[!(predictors.all %in% predictors)]] <- NA

	# Storing the predicted values from the gam
	data$fit.gam <- predict(gam1, newdata=data)
	# ----------------------------------------

	# ----------------------------------------
	# Run all of the post-processing (calculate CIs, etc)
	# ----------------------------------------
	mod.out <- process.gamm(gamm.model=gam1, data=data, model.name=model.name, extent=extent, resolution=resolution, response=response, vars=predictors, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	# ----------------------------------------
	
	return(mod.out)
}

