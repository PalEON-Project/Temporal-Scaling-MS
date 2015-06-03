factor.weights <- function(model.gam, newdata, sites=F){
	# If the model.gam is a mixed model.gam (gamm) rather than a normal gam, extract just the gam portion
	if(class(model.gam)[[1]]=="gamm") model.gam <- model.gam$gam
	# -----------
	# calculating the weights for each of the factors
	# -----------
	# Create the prediction matrix
	Xp <- predict(model.gam, newdata=newdata, type="lpmatrix")

	fit <- Xp %*% coef(model.gam) # The full predicted values; used for model.gam QA/QC
	coef.gam <- coef(model.gam) # the gam coefficients
	
	# Some handy column indices
	cols.site <- if(sites==T) which(substr(names(coef.gam),1,4)=="Site") else 1
	cols.temp   <- which(substr(names(coef.gam),1,7)=="s(Temp)")
	cols.precip <- which(substr(names(coef.gam),1,9)=="s(Precip)")
	cols.co2    <- which(substr(names(coef.gam),1,6)=="s(CO2)")

	# calculating the smoother for each effect
	if(sites==T) {
		fit.int<- Xp[,cols.site]   %*% coef.gam[cols.site] # Note: no matrix multiplication because it's 1 x 1

	} else {
		fit.int<- Xp[,cols.site]    *  coef.gam[cols.site] # Note: no matrix multiplication because it's 1 x 1
	}
	fit.temp   <- Xp[,cols.temp]   %*% coef.gam[cols.temp]
	fit.precip <- Xp[,cols.precip] %*% coef.gam[cols.precip]
	fit.co2    <- Xp[,cols.co2]    %*% coef.gam[cols.co2]

	# Calculated the SD around each smoother
	if(sites==T){
		sd.int<- rowSums(Xp[,cols.site]   %*% model.gam$Vp[cols.site, cols.site]    *Xp[,cols.site]    )^0.5
	} else {
		sd.int<-    sum (Xp[,cols.site]    *  model.gam$Vp[cols.site, cols.site]    *Xp[,cols.site]    )^0.5
	}
	sd.temp   <- rowSums(Xp[,cols.temp]   %*% model.gam$Vp[cols.temp, cols.temp]    *Xp[,cols.temp]    )^0.5
	sd.precip <- rowSums(Xp[,cols.precip] %*% model.gam$Vp[cols.precip, cols.precip]*Xp[, cols.precip] )^0.5
	sd.co2    <- rowSums(Xp[,cols.co2]    %*% model.gam$Vp[cols.co2, cols.co2]      *Xp[, cols.co2]    )^0.5

	# Summing the fixed effects to do QA/QC and throwing a warning if it's not very close to the predicted values
	fit.sum <- fit.int + fit.co2 + fit.temp + fit.precip
	fit.spline <- fit.co2 + fit.temp + fit.precip
	if(max(abs(fit - fit.sum))>1e-4) print("***** WARNING: sum of fixed effects not equal to predicted value *****")

	# summing the absolute values to get the weights for each fixed effect
	fit.sum2 <- abs(fit.int) + abs(fit.co2) + abs(fit.temp) + abs(fit.precip)
	fit.spline2 <- abs(fit.co2) + abs(fit.temp) + abs(fit.precip)

	# Factor weights are determined by the relative strength of Temp, Precip, & CO2
	df.weights <- data.frame(Site=newdata$Site, Scale=newdata$Scale, Year=newdata$Year, Temp=newdata$Temp, Precip=newdata$Precip, CO2=newdata$CO2, fit=fit, fit.intercept=fit.int, fit.temp=fit.temp, fit.precip=fit.precip, fit.co2=fit.co2, sd.temp=sd.temp, sd.precip=sd.precip, sd.co2=sd.co2,  fit.spline=fit.spline, weight.co2=fit.co2/fit.spline2, weight.temp=fit.temp/fit.spline2, weight.precip=fit.precip/fit.spline2)
	# # Add in a couple factors that are useful if predicting on the old data
	# if(!is.null(newdata$Year)) df.weights$Year <- newdata$Year 
	# if(!is.null(newdata$Site)) df.weights$Year <- newdata$Site 

	# doing a little bit of handy-dandy calculation to give a flag as to which factor is given the greatest weight in a given year
	for(i in 1:nrow(df.weights)){
		fweight <- abs(df.weights[i,c("weight.co2", "weight.temp", "weight.precip")])
		df.weights[i,"max"] <- max(fweight)
		df.weights[i,"factor.max"] <- c("co2", "temp", "precip")[which(fweight==max(fweight))]
	}
	df.weights$factor.max <- as.factor(df.weights$factor.max)
	# summary(df.weights)
	return(df.weights)
}