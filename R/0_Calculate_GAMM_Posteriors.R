post.distns <- function(model.gam, model.name, newdata, n, terms=T, sites=F, lwr=0.025, upr=0.975){
	# Note: this function can be used to generate a 95% CI on the full model.gam OR terms

	# -----------
	# Simulating a posterior distribution of Betas to get variance on non-linear functions
	# This is following Gavin Simpson's post here: http://www.fromthebottomoftheheap.net/2011/06/12/additive-model.gamling-and-the-hadcrut3v-global-mean-temperature-series/
	# -----------
	library(MASS)
	set.seed(321)

	# If the model.gam is a mixed model.gam (gamm) rather than a normal gam, extract just the gam portion
	if(class(model.gam)[[1]]=="gamm") model.gam <- model.gam$gam


	coef.gam <- coef(model.gam)
	# Generate a random distribution of betas
	Rbeta <- mvrnorm(n=n, coef(model.gam), vcov(model.gam))

	# Create the prediction matrix
	Xp <- predict(model.gam, newdata=newdata, type="lpmatrix")

	# Some handy column indices
	cols.site   <- which(substr(names(coef.gam),1,4)=="Site" | substr(names(coef.gam),1,11)=="(Intercept)")
	cols.temp   <- which(substr(names(coef.gam),1,7)=="s(Temp)")
	cols.precip <- which(substr(names(coef.gam),1,9)=="s(Precip)")
	cols.co2    <- which(substr(names(coef.gam),1,6)=="s(CO2)")

	if(terms==T){
		sim.temp   <- Xp[,cols.temp]   %*% t(Rbeta[,cols.temp]) 
		sim.precip <- Xp[,cols.precip] %*% t(Rbeta[,cols.precip]) 
		sim.co2    <- Xp[,cols.co2]    %*% t(Rbeta[,cols.co2]) 
		
		out.temp <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Effect="Temp", x=newdata$Temp,
							   mean=apply(sim.temp, 1, mean), 
							   lwr=apply(sim.temp, 1, quantile, lwr, na.rm=T), 
							   upr=apply(sim.temp, 1, quantile, upr, na.rm=T))
		out.precip <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Effect="Precip", x=newdata$Precip, 
		                         mean=apply(sim.precip, 1, mean, na.rm=T), 
		                         lwr=apply(sim.precip, 1, quantile, lwr, na.rm=T), 
		                         upr=apply(sim.precip, 1, quantile, upr, na.rm=T))
		out.co2 <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Effect="CO2", x=newdata$CO2, 
		                      mean=apply(sim.co2, 1, mean, na.rm=T), 
		                      lwr=apply(sim.co2, 1, quantile, lwr, na.rm=T), 
		                      upr=apply(sim.co2, 1, quantile, upr, na.rm=T))

		df.out <- rbind(out.temp, out.precip, out.co2)

		# Creating a data frame storing all the simulations for more robust analyses
		sim.bind <- rbind(sim.temp, sim.precip, sim.co2)

		df.sim <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Effect=c(rep("Temp", length(newdata$Temp)), rep("Precip", length(newdata$Precip)), rep("CO2", length(newdata$CO2))), x=c(newdata$Temp, newdata$Precip, newdata$CO2))
		df.sim <- cbind(df.sim, sim.bind)

	} else {
		sim1 <- Xp %*% t(Rbeta) # simulates n predictions of the response variable in the model.gam
		
		df.out <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Year=newdata$Year, Temp=newdata$Temp, Precip=newdata$Precip, CO2=newdata$CO2, mean=apply(sim1, 1, mean, na.rm=T), lwr=apply(sim1, 1, quantile, lwr, na.rm=T), upr=apply(sim1, 1, quantile, upr, na.rm=T))
		
		df.sim <- data.frame(Model=model.name, Site=newdata$Site, Extent=newdata$Extent, Scale=newdata$Scale, Year=newdata$Year, Temp=newdata$Temp, Precip=newdata$Precip, CO2=newdata$CO2)
		df.sim <- cbind(df.sim, sim1)

	}

	
	out <- list()
	out[["ci"]]	 <- df.out
	out[["sims"]] <- df.sim
			
	return(out)
}