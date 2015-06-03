model.gam <- function(data, model, response, scale="", k=4, outdir){
	# data     = data frame with data.temp in it
	# model    = which model to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata file
	library(mgcv)
	
	# creating a working data.temp frame with just the data.temp we want
	data.temp          <- data[data$Model==model, c("Model", "Updated", "Model.Order", "Site", "Year")]
	data.temp$Scale    <- as.factor(paste0("t", scale))
	data.temp$response <- data[data$Model==model, paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model, paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model, paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model, paste0("CO2", scale)]	
	
	
    # -----------
	# Running the basic model
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	# This is different from model.site.gam in that it has an added random site slope.
	#	This random effect lets us gauge the overall model response to our fixed effects 
	#   regardless of the site.  
	#   Pros: Generalized and helps characterize the basic model
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1))
	print(summary(gam1$gam))	
	# return(assign(paste0("gam.", model), gam1)) # Saving the gam a name based on the model

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata= data.temp)
	
	# return(assign(paste0("data.temp.", model), data.temp) )# saving the temporary data.temp frame for use in post-hoc QA/QC
	# -----------

	# -----------
	# calculating the weights for each of the factors
	# -----------
	# Create the prediction matrix
	Xp <- predict(gam1$gam, newdata=data.temp, type="lpmatrix")

	fit <- Xp %*% coef(gam1$gam) # The full predicted values; used for model QA/QC
	coef.gam <- coef(gam1$gam) # the gam coefficients
	
	# Some handy column indices
	cols.site <- which(substr(names(coef.gam),1,4)=="Site")
	# cols.site    <- 1 # Running on 1 site so the "site" intercept is always 1
	cols.temp   <- which(substr(names(coef.gam),1,7)=="s(Temp)")
	cols.precip <- which(substr(names(coef.gam),1,9)=="s(Precip)")
	cols.co2    <- which(substr(names(coef.gam),1,6)=="s(CO2)")

	# calculating the smoother for each effect
	fit.int <- Xp[,cols.site] * coef.gam[cols.site]
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
	
	out <- list(data=data.temp, gam=gam1, weights=factor.weights)
	save(out, file=file.path(outdir, paste0("gam.", model, response, ".Rdata")))
	return(out)
	# -----------

}