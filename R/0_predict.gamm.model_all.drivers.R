process.gamm <- function(gamm.model, data, model.name, response, vars, scale="", extent=c(850,2010), k=4, 
						  	fweights=T, ci.model=T, ci.terms=T, n=250,
						  	write.out=T, outdir, control=list()){
	# data     = data frame with data in it
	# model.name    = which model.name to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# n        = number of simulations for generating confidence intervals
	# outdir   = where to save the .Rdata file
	library(mgcv)
	source("R/0_Calculate_GAMM_Weights.R")
	source("R/0_Calculate_GAMM_Posteriors.R")
    # -----------
	# Running the basic model.name
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	# This is different from model.site.gam in that it has an added random site slope.
	#	This random effect lets us gauge the overall model.name response to our fixed effects 
	#   regardless of the site.  
	#   Pros: Generalized and helps characterize the basic model responses
	#   Cons: Sloooooooooooow! (or at least slow with the PalEON data set)


	# Storing the predicted values from the gam
	data$fit.gam <- predict(gamm.model$gam, newdata= data)
	
	out <- list(data=data, gamm=gamm.model)
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gamm.model, model.name=model.name, newdata=data, extent=extent, vars=vars, sites=T); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=data, vars=vars, terms=F, sites=T)
		out[["ci.response"]] <- ci.response$ci 
		out[["sim.response"]] <- ci.response$sims 

	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = n
		sites.dat <- unique(data$Site)
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
								Scale=rep(t.scale, n.out*ns))
		for(v in vars){
			new.dat[,v] <- rep(seq(min(data[,v],   na.rm=T), max(data$[,v],   na.rm=T), length.out=n.out), ns)
		}								
								
		ci.terms.pred <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=new.dat, vars=vars, terms=T, sites=T)

		out[["ci.terms"]] <- ci.terms.pred$ci 
		out[["sim.terms"]] <- ci.terms.pred$sims 
	}	
	# -----------
	if(write.out==T) save(out, file=file.path(outdir, paste("gamm", model.name, response, ifelse(scale=="",".001", scale), "AllSites", "Rdata", sep=".")))
	return(out)
	# -----------

}