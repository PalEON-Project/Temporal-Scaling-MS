model.gam <- function(data, model.name, response, scale="", k=4, outdir, 
						  	fweights=T, ci.model=T, ci.terms=T){
	# data     = data frame with data.temp in it
	# model.name    = which model.name to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata file
	library(mgcv)
	source("R/Calculate_GAMM_Weights.R")
	source("R/Calculate_GAMM_Posteriors.R")

	# creating a working data.temp frame with just the data.temp we want
	data.temp          <- data[data$Model==model.name, c("Model", "Updated", "Model.Order", "Site", "Year")]
	data.temp$Scale    <- as.factor(paste0("t", scale))
	data.temp$response <- data[data$Model==model.name, paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model.name, paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model.name, paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model.name, paste0("CO2", scale)]	
	
	
    # -----------
	# Running the basic model.name
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	# This is different from model.site.gam in that it has an added random site slope.
	#	This random effect lets us gauge the overall model.name response to our fixed effects 
	#   regardless of the site.  
	#   Pros: Generalized and helps characterize the basic model.name
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata= data.temp)
	
	out <- list(data=data.temp, gamm=gam1)
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gam1, newdata=data.temp, sites=T); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gam1, newdata=data.temp, terms=F, sites=T)
		out[["ci.response"]] <- ci.response 
	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = 100
		new.dat <- data.frame(	Temp  =seq(min(data.temp$Temp),   max(data.temp$Temp),   length.out=n.out),
								Precip=seq(min(data.temp$Precip), max(data.temp$Precip), length.out=n.out),
								CO2   =seq(min(data.temp$CO2),    max(data.temp$CO2),    length.out=n.out))
		ci.terms.pred <- post.distns(model.gam=gam1, newdata=new.dat, terms=T, sites=T)

		out[["ci.terms"]] <- ci.terms.pred 
	}	
	# -----------
		
	save(out, file=file.path(outdir, paste("gamm", model.name, response, ifelse(scale=="","001", scale), "AllSites", "Rdata", sep=".")))
	return(out)
	# -----------

}