model.site.gam <- function(	data, model, site, response, scale="", k=4, outdir, 
						  	fweights=T, ci.model=T, ci.termps=T){
	# data     = data frame with data.temp in it
	# model    = which model to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata.temp file
	# fweights = do we compute the weights for each of the smoothing terms through time? (logical) 
	# ci.model = do we compute a CI for the model prediction of the response variable? (logical) 
	# ci.terms = do we compute a CI smoothing terms over the range of our observations? (logical) 

	library(mgcv)
	source("R/Calculate_GAMM_Weights.R")
	source("R/Calculate_GAMM_Posteriors.R")

	# creating a working data.temp frame with just the data.temp we want
	data.temp          <- data[data$Model==model & data$Site==site, c("Model", "Updated", "Model.Order", "Site", "Year")]
	data.temp$Scale    <- as.factor(paste0("t", scale))
	data.temp$response <- data[data$Model==model & data$Site==site, paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model & data$Site==site, paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model & data$Site==site, paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model & data$Site==site, paste0("CO2", scale)]	
	
	
    # -----------
	# Running the basic model
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), data=data.temp, correlation=corARMA(form=~Year, p=1))
	print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata=data.temp)

	out <- list(data=data.temp, gam=gam1)
	# -----------
	
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights=T){	
		f.weights <- factor.weights(model=gam1, newdata=data.temp); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model=T){
		ci.response <- factor.weights(model=gam1, newdata=data.temp, terms=F)
		out[["ci.response"]] <- ci.response 
	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms=T){
		n.out = 100
		new.dat <- data.frame(	Temp  =seq(min(data.temp$Temp),   max(data.temp$Temp),   length.out=n.out),
								Precip=seq(min(data.temp$Precip), max(data.temp$Precip), length.out=n.out),
								CO2   =seq(min(data.temp$CO2),    max(data.temp$CO2),    length.out=n.out))
		ci.terms.pred <- factor.weights(model=gam1, newdata=new.dat, terms=T)

		out[["ci.terms"]] <- ci.terms.pred 
	}	
	# -----------
	
	
	save(out, file=file.path(outdir, paste0("gam.", model, response, ".Rdata.temp")))
	return(out)
	# -----------

}