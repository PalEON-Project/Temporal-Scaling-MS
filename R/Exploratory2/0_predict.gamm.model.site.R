model.site.gam <- function(	data, model.name, site, response, scale="", extent=c(850,2010), k=4, 
						  	fweights=T, ci.model=T, ci.terms=T, n=250,
						  	write.out=T, outdir=""){
	# data     = data frame with data.temp in it
	# model.name    = which model.name to subset
	# response = which variable to use as response in the gam
	# k        = number of knots in the spline
	# outdir   = where to save the .Rdata file
	# fweights = do we compute the weights for each of the smoothing terms through time? (logical) 
	# ci.model = do we compute a CI for the model.name prediction of the response variable? (logical) 
	# ci.terms = do we compute a CI smoothing terms over the range of our observations? (logical) 

	library(mgcv)
	source("R/0_Calculate_GAMM_Weights.R")
	source("R/0_Calculate_GAMM_Posteriors.R")

	# creating a working data.temp frame with just the data.temp we want
	t.scale <- ifelse(scale=="", "t.001", paste0("t", scale))

	data.temp          <- data[data$Model==model.name & data$Site==site, c("Model", "Updated", "Model.Order", "Site", "Extent", "Year")]
	data.temp$Scale    <- as.factor(t.scale)
	data.temp$response <- data[data$Model==model.name & data$Site==site, paste0(response, scale)]
	data.temp$Temp     <- data[data$Model==model.name & data$Site==site, paste0("Temp", scale)]	
	data.temp$Precip   <- data[data$Model==model.name & data$Site==site, paste0("Precip", scale)]	
	data.temp$CO2      <- data[data$Model==model.name & data$Site==site, paste0("CO2", scale)]	
	
	
    # -----------
	# Running the basic model.name
	# -----------
	# Running the gamm; note this now has AR1 temporal autocorrelation
	gam1 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k) , data=data.temp, correlation=corARMA(form=~Year, p=1), control=list(sing.tol=1e-20, opt="optim"))
	# control=list(niterEM=0,sing.tol=1e-20)
	# gam2 <- gamm(response ~ s(Temp, k=k) + s(Precip, k=k) + s(CO2, k=k), data=data.temp)
	print(summary(gam1$gam))	

	# Storing the predicted values from the gam
	data.temp$fit.gam <- predict(gam1$gam, newdata=data.temp)

	out <- list(data=data.temp, gamm=gam1)
	# -----------
	
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gam1, model.name=model.name, newdata=data.temp, extent=extent, vars=c("Temp", "Precip", "CO2"), sites=F); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gam1, model.name=model.name, n=n, newdata=data.temp, terms=F, sites=F)
		out[["ci.response"]] <- ci.response$ci 
		out[["sim.response"]] <- ci.response$sims 
	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = n
		
		new.dat <- data.frame(	Site=rep(site, n.out), 
								Extent=as.factor(paste(extent[1], extent[2], sep="-")),
								Scale=rep(t.scale, n.out),
							    Temp  =seq(min(data.temp$Temp,   na.rm=T), max(data.temp$Temp,   na.rm=T), length.out=n.out),
								Precip=seq(min(data.temp$Precip, na.rm=T), max(data.temp$Precip, na.rm=T), length.out=n.out),
								CO2   =seq(min(data.temp$CO2,    na.rm=T), max(data.temp$CO2,    na.rm=T), length.out=n.out))
		ci.terms.pred <- post.distns(model.gam=gam1, model.name=model.name, n=n, newdata=new.dat, terms=T, sites=F)

		out[["ci.terms"]] <- ci.terms.pred$ci 
		out[["sim.terms"]] <- ci.terms.pred$sims 
	}	
	# -----------
	
	
	if(write.out==T) save(out, file=file.path(outdir, paste("gamm", model.name, response, ifelse(scale=="","001", scale), site, "Rdata", sep=".")))
	return(out)
	# -----------

}