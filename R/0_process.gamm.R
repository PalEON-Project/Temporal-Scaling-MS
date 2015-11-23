process.gamm <- function(gamm.model, data, model.name, response, vars, resolution="t.001", extent=c(850,2010), 
						  	fweights=T, sites=T, ci.model=T, ci.terms=T, PFT=F, n=250,
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

	out <- list(data=data, gamm=gamm.model)
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gamm.model, model.name=model.name, newdata=data, extent=extent, vars=vars); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=data, vars=vars, terms=F)
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
		if(PFT==F){	
		  new.dat <- data.frame(	Site=site.vec,
								Extent=as.factor(paste(extent[1], extent[2], sep="-")), 
								Resolution=rep(resolution, n.out*ns))
		} else {
		  pfts <- unique(data$PFT)
		  new.dat <- data.frame(Site       = rep(site.vec, length(pfts)),
								Extent     = as.factor(paste(extent[1], extent[2], sep="-")), 
								Resolution = rep(rep(resolution, n.out*ns), length(pfts))
							    )
		  for(i in 1:length(pfts)){
		  	new.dat[(i*n.out*ns - n.out*ns + 1):(i*n.out*ns),"PFT"] <- pfts[i]
		  }
		}
		for(v in vars){
			if(is.factor(data[,v])){
				new.dat[,v] <- as.factor(unique(data[,v])[1])
			} else {
				new.dat[,v] <- rep(seq(min(data[,v],   na.rm=T), max(data[,v],   na.rm=T), length.out=n.out), ns)
			}
		}								
								
		ci.terms.pred <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=new.dat, vars=vars, terms=T, PFT=PFT)

		out[["ci.terms"]] <- ci.terms.pred$ci 
		out[["sim.terms"]] <- ci.terms.pred$sims 
	}	
	# -----------
	if(write.out==T) save(out, file=file.path(outdir, paste("gamm", model.name, response, resolution, "AllSites", "Rdata", sep=".")))
	return(out)
	# -----------

}