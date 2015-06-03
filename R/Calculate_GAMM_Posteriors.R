post.distns <- function(model.gam, newdata, n=1000, n2=25, terms=T, sites=F, lwr=0.025, upr=0.975){
	# Note: this function can be used to generate a 95% CI on the full model.gam OR terms

	# -----------
	# Simulating a posterior distribution of Betas to get variance on non-linear functions
	# This is following Gavin Simpson's post here: http://www.fromthebottomoftheheap.net/2011/06/12/additive-model.gamling-and-the-hadcrut3v-global-mean-temperature-series/
	# -----------
	library(MASS)
	set.seed(321)

	# If the model.gam is a mixed model.gam (gamm) rather than a normal gam, extract just the gam portion
	if(class(model.gam)[[1]]=="gamm") model.gam <- model.gam$gam



	# Generate a random distribution of betas
	Rbeta <- mvrnorm(n=n, coef(model.gam), vcov(model.gam))

	# Create the prediction matrix
	Xp <- predict(model.gam, newdata=newdata, type="lpmatrix")

	# Some handy column indices
	cols.site <- if(sites==T) which(substr(names(coef.gam),1,4)=="Site") else 1
	cols.temp   <- which(substr(names(coef.gam),1,7)=="s(Temp)")
	cols.precip <- which(substr(names(coef.gam),1,9)=="s(Precip)")
	cols.co2    <- which(substr(names(coef.gam),1,6)=="s(CO2)")

	want <- sample(n, n2) # sample n2 number of ranges from sim1
	
	if(terms==T){
		sim.temp   <- Xp[,cols.temp]   %*% t(Rbeta[,cols.temp]) 
		sim.precip <- Xp[,cols.precip] %*% t(Rbeta[,cols.precip]) 
		sim.co2    <- Xp[,cols.co2]    %*% t(Rbeta[,cols.co2]) 
		
		df.out <- list()
		df.out[["Temp"]] <- data.frame(Temp=newdata$Temp, mean=apply(sim.temp, 1, mean), lwr=apply(sim.temp, 1, quantile, lwr), upr=apply(sim.temp, 1, quantile, upr))
		df.out[["Precip"]] <- data.frame(Precip=newdata$Precip, mean=apply(sim.precip, 1, mean), lwr=apply(sim.precip, 1, quantile, lwr), upr=apply(sim.precip, 1, quantile, upr))
		df.out[["CO2"]] <- data.frame(CO2=newdata$CO2, mean=apply(sim.co2, 1, mean), lwr=apply(sim.co2, 1, quantile, lwr), upr=apply(sim.co2, 1, quantile, upr))

	} else {
		sim1 <- Xp %*% t(Rbeta) # simulates n predictions of the response variable in the model.gam
		
		df.out <- data.frame(Year=newdata$Year, Temp=newdata$Temp, Precip=newdata$Precip, CO2=newdata$CO2, mean=apply(sim1, 1, mean), lwr=apply(sim1, 1, quantile, lwr), upr=apply(sim1, 1, quantile, upr))
	}

	return(df.out)
}