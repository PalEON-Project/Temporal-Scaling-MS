post.distns <- function(model.gam, model.name, newdata, vars, n, terms=T, lwr=0.025, upr=0.975){
	# Note: this function can be used to generate a 95% CI on the full model.gam OR terms

	# -----------
	# Simulating a posterior distribution of Betas to get variance on non-linear functions
	# This is following Gavin Simpson's post here: 
	# http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/
	# His handy-dandy functions can be found here: https://github.com/gavinsimpson/random_code/
	#      Including the derivative funcition that will probably come in handy later
	# -----------
	library(MASS)
	set.seed(321)

	# If the model.gam is a mixed model.gam (gamm) rather than a normal gam, extract just the gam portion
	if(class(model.gam)[[1]]=="gamm") model.gam <- model.gam$gam


	coef.gam <- coef(model.gam)
	# Generate a random distribution of betas using the covariance matrix
	Rbeta <- mvrnorm(n=n, coef(model.gam), vcov(model.gam))

	# Create the prediction matrix
	Xp <- predict(model.gam, newdata=newdata, type="lpmatrix")

	# Some handy column indices
	cols.list <- list(Site = which(substr(names(coef.gam),1,4)=="Site" | substr(names(coef.gam),1,11)=="(Intercept)"))
	for(v in vars){
		cols.list[[v]] <- which(substr(names(coef.gam),1,(nchar(v)+3))==paste0("s(",v,")"))
	}

	# sim.list <- list()
	if(terms==T){
		for(v in vars){
			sim.tmp <- data.frame(Xp[,cols.list[[v]]] %*% t(Rbeta[,cols.list[[v]]]) )

			# Saving the quantiles into a data frame
			df.tmp <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Resolution=newdata$Resolution, Effect=v, x=newdata[,v],
							   mean=apply(sim.tmp, 1, mean), 
							   lwr=apply(sim.tmp, 1, quantile, lwr, na.rm=T), 
							   upr=apply(sim.tmp, 1, quantile, upr, na.rm=T))

			if(v == vars[1]) df.out <- df.tmp else df.out <- rbind(df.out, df.tmp)

  			# Creating a data frame storing all the simulations for more robust analyses
			sim.tmp$Model  <- model.name
			sim.tmp$Site   <- newdata$Site
			sim.tmp$Extent <- newdata$Extent
			sim.tmp$Resolution  <- newdata$Resolution
			sim.tmp$Effect <- v
			sim.tmp$x      <- newdata[,v]
			sim.tmp <- sim.tmp[,c((n+1):ncol(sim.tmp), 1:n)]
 
			if(v == vars[1]) df.sim <- sim.tmp else df.sim <- rbind(df.sim, sim.tmp)

		}

	} else {
		sim1 <- Xp %*% t(Rbeta) # simulates n predictions of the response variable in the model.gam
		
		df.out <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Resolution=newdata$Resolution, Year=newdata$Year, mean=apply(sim1, 1, mean, na.rm=T), lwr=apply(sim1, 1, quantile, lwr, na.rm=T), upr=apply(sim1, 1, quantile, upr, na.rm=T))
		df.sim <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Resolution=newdata$Resolution, Year=newdata$Year)
		for(v in vars){
			df.out[,v] <- newdata[,v]
			df.sim[,v] <- newdata[,v]
		}
		df.sim <- cbind(df.sim, sim1)

	}

	
	out <- list()
	out[["ci"]]	 <- df.out
	out[["sims"]] <- df.sim
			
	return(out)
}