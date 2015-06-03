model.site.gam <- function(	data, model.name, site, response, scale="", k=4, outdir, 
						  	fweights=T, ci.model=T, ci.terms=T){
	# data     = data frame with data.temp in it
	# model.name    = which model.name to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata file
	# fweights = do we compute the weights for each of the smoothing terms through time? (logical) 
	# ci.model = do we compute a CI for the model.name prediction of the response variable? (logical) 
	# ci.terms = do we compute a CI smoothing terms over the range of our observations? (logical) 

	library(mgcv)
	source("R/Calculate_GAMM_Weights.R")
	source("R/Calculate_GAMM_Posteriors.R")

	# creating a working data.temp frame with just the data.temp we want
	data.temp          <- data[data$Model==model.name & data$Site==site, c("Model", "Updated", "Model.Order", "Site", "Year")]
	data.temp$Scale    <- as.factor(paste0("t", scale))
	data.temp$response <- data[data$Model==model.name & data$Site==site, paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model.name & data$Site==site, paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model.name & data$Site==site, paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model.name & data$Site==site, paste0("CO2", scale)]	
	
	
    # -----------
	# Running the basic model.name
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), data=data.temp, correlation=corARMA(form=~Year, p=1))
	print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata=data.temp)

	out <- list(data=data.temp, gamm=gam1)
	# -----------
	
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gam1, newdata=data.temp, sites=F); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gam1, newdata=data.temp, terms=F, sites=F)
		out[["ci.response"]] <- ci.response 
	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = 200
		new.dat <- data.frame(	Temp  =seq(min(data.temp$Temp),   max(data.temp$Temp),   length.out=n.out),
								Precip=seq(min(data.temp$Precip), max(data.temp$Precip), length.out=n.out),
								CO2   =seq(min(data.temp$CO2),    max(data.temp$CO2),    length.out=n.out))
		ci.terms.pred <- post.distns(model.gam=gam1, newdata=new.dat, terms=T, sites=F)

		out[["ci.terms"]] <- ci.terms.pred 
	}	
	# -----------
	
	
	save(out, file=file.path(outdir, paste("gamm", model.name, response, ifelse(scale=="","001", scale), site, "Rdata", sep=".")))
	return(out)
	# -----------

}