model.gam <- function(data, model.name, response, scale="", extent=c(850,2010), k=4, 
						  	fweights=T, ci.model=T, ci.terms=T, n=250,
						  	write.out=T, outdir, control=list()){
	# data     = data frame with data.temp in it
	# model.name    = which model.name to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# n        = number of simulations for generating confidence intervals
	# outdir   = where to save the .Rdata file
	library(mgcv)
	source("R/0_Calculate_GAMM_Weights.R")
	source("R/0_Calculate_GAMM_Posteriors.R")

	# creating a working data.temp frame with just the data.temp we want
	# creating a working data.temp frame with just the data.temp we want
	t.scale <- ifelse(scale=="", "t.001", paste0("t", scale))

	data.temp          <- data[data$Model==model.name, c("Model", "Updated", "Model.Order", "Site", "Extent", "Year")]
	data.temp$Scale    <- as.factor(t.scale)
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
	#   Pros: Generalized and helps characterize the basic model responses
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)

	# Each of the models is having different stability issues
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k) + Site -1, random=list(Site=~Site), data=data.temp, correlation=corARMA(form=~Year, p=1), control=control)
  print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata= data.temp)
	
	out <- list(data=data.temp, gamm=gam1)
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gam1, model.name=model.name, newdata=data.temp, extent=extent, vars=c("Temp", "Precip", "CO2"), sites=T); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gam1, model.name=model.name, n=n, newdata=data.temp, vars=c("Temp", "Precip", "CO2"), terms=F, sites=T)
		out[["ci.response"]] <- ci.response$ci 
		out[["sim.response"]] <- ci.response$sims 

	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = n
		sites.dat <- unique(data.temp$Site)
		ns        <- length(sites.dat)
		for(i in 1:ns){
			if(i == 1){ 
				site.vec <- paste(rep(sites.dat[i], n.out) )
			} else {
				site.vec <- c(site.vec, paste(rep(sites.dat[i], n.out)))
			}
		}
	
		new.dat <- data.frame(	Site=site.vec,
								Extent=as.factor(paste(extent[1], extent[2], sep="-")), 
								Scale=rep(t.scale, n.out*ns),
							    Temp  =rep(seq(min(data.temp$Temp,   na.rm=T), max(data.temp$Temp,   na.rm=T), length.out=n.out), ns),
								Precip=rep(seq(min(data.temp$Precip, na.rm=T), max(data.temp$Precip, na.rm=T), length.out=n.out), ns),
								CO2   =rep(seq(min(data.temp$CO2,    na.rm=T), max(data.temp$CO2,    na.rm=T), length.out=n.out),ns))
		ci.terms.pred <- post.distns(model.gam=gam1, model.name=model.name, n=n, newdata=new.dat, vars=c("Temp", "Precip", "CO2"), terms=T, sites=T)

		out[["ci.terms"]] <- ci.terms.pred$ci 
		out[["sim.terms"]] <- ci.terms.pred$sims 
	}	
	# -----------
	if(write.out==T) save(out, file=file.path(outdir, paste("gamm", model.name, response, ifelse(scale=="",".001", scale), "AllSites", "Rdata", sep=".")))
	return(out)
	# -----------

}