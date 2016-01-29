factor.weights <- function(model.gam, model.name, newdata, extent, vars){
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
	# Note: In this script "Site" refers to any intercept (it can be PlotID, PFT, whatever)
	cols.list <- list()
	intercepts = which(!substr(names(coef.gam),1,2)=="s(")
	if(length(intercepts)>0){
		cols.list[["Site"]] <- intercepts
	}
	for(v in vars[!vars=="Y.lag"]){
		cols.list[[v]] <- which(substr(names(coef.gam),1,(nchar(v)+3))==paste0("s(",v,")"))
	}
	# calculating the smoother for each effect; 
	# now storing everything in a data frame
	gam.fits <- data.frame(intercept=rep(0, nrow(newdata)))

	if(length(cols.list[["Site"]])>1) {
		gam.fits[,"intercept"] <- as.vector(Xp[,cols.list[["Site"]]]   %*% coef.gam[cols.list[["Site"]]] )

	} else if(length(cols.list[["Site"]])==1) {
		gam.fits[["intercept"]] <- as.vector(t(Xp[,cols.list[["Site"]]]    *  coef.gam[cols.list[["Site"]]])) # Note: no matrix multiplication because it's 1 x 1
	}

	for(v in vars){
		gam.fits[,v] <- as.vector(Xp[,cols.list[[v]]]   %*% coef.gam[cols.list[[v]]])
	}

	# Calculated the SD around each smoother
	gam.sd <- data.frame(intercept=vector(length=nrow(newdata)))
	if(length(cols.list[["Site"]])>1){
		gam.sd[,"intercept"] <- rowSums(Xp[,cols.list[["Site"]]] %*% model.gam$Vp[cols.list[["Site"]], cols.list[["Site"]]] * Xp[,cols.list[["Site"]]] )^0.5
	} else if(length(cols.list[["Site"]])==1){
		gam.sd[,"intercept"] <-    sum (Xp[,cols.list[["Site"]]]  *  model.gam$Vp[cols.list[["Site"]], cols.list[["Site"]]] * Xp[,cols.list[["Site"]]] )^0.5
	}

	for(v in vars){
		gam.sd[,v] <- as.vector(rowSums(Xp[,cols.list[[v]]] %*% model.gam$Vp[cols.list[[v]], cols.list[[v]]] * Xp[,cols.list[[v]]] )^0.5)
	}

	# Summing the fixed effects to do QA/QC and throwing a warning if it's not very close to the predicted values
	fit.sum <- rowSums(gam.fits, na.rm=T)
	fit.spline <- rowSums(gam.fits[,2:ncol(gam.fits)], na.rm=T)
	if(max(abs(fit - fit.sum),na.rm=T)>1e-4) print("***** WARNING: sum of fixed effects not equal to predicted value *****")

	# summing the absolute values to get the weights for each fixed effect
	fit.sum2 <- rowSums(abs(gam.fits[,]), na.rm=T)
	fit.spline2 <- rowSums(abs(gam.fits[,2:ncol(gam.fits)]), na.rm=T)

	# Factor weights are determined by the relative strength of Temp, Precip, & CO2
	df.weights <- data.frame(Model=model.name, 
	                         Site=newdata$Site, 
	                         Extent=newdata$Extent, 
	                         Resolution=newdata$Resolution, 
	                         Year=newdata$Year, 
	                         fit.full=fit)

 	if("PlotID" %in% names(newdata)) df.weights$PlotID <- newdata$PlotID
	if("TreeID" %in% names(newdata)) df.weights$TreeID <- newdata$TreeID
	if("PFT"    %in% names(newdata)) df.weights$PFT    <- newdata$PFT

	for(v in vars){
		df.weights[,paste("fit", v, sep=".")   ] <- gam.fits[,v]
		df.weights[,paste( "sd", v, sep=".")   ] <- gam.sd[,v]
		df.weights[,paste("weight", v, sep=".")] <- gam.fits[,v]/fit.spline2
	}
	
	# doing a little bit of handy-dandy calculation to give a flag as to which factor is given the greatest weight in a given year
	# rows.na <- which(is.na(df.weights$weight.co2))	
	# rows.seq <- 1:nrow(df.weights) 
	# for(i in rows.seq[!(rows.seq %in% rows.na)]){
	cols.weights <- which(substr(names(df.weights),1,6)=="weight")
	for(i in 1:nrow(df.weights)){
		fweight <- abs(df.weights[i,cols.weights])
		df.weights[i,"max"] <- max(fweight, na.rm=T)
		df.weights[i,"factor.max"] <- vars[which(fweight==max(fweight))]
	}
	df.weights$factor.max <- as.factor(df.weights$factor.max)
	# summary(df.weights)
	return(df.weights)
}