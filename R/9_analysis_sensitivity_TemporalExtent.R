# ----------------------------------------
# Objective: Compare model/tree ring sensitivites by Temporal Extent
# Christy Rollinson, crollinson@gmail.com
# Date Created: 16 November 2015
# ----------------------------------------
#
# -------------------------
# Hypotheses: 
# -------------------------
# 1. Increased climate variability & range of conditions with longer 
#    temporal extents reduces climate sensitivity to the conditions 
#    in the modern era
# 2. Ecosystem adaptation such as changes in composition through time 
#    mediate long-term climate sensitivity 
# -------------------------
#
# -------------------------
# Workflow
# -------------------------
# 1. Set Directories
# 2. Load data files & function scripts
# 3. Refit all gamms to the range of conditions used in the full runs 
#    -- (we will shade out these "extra" regions in the graph, but I 
#        want to show what happens when you extrapolate)
# 4. Standardize driver responses to the mean model NPP to aid 
#    comparison (loop)
#     a. Find the NPP to relativize each set off of
#     b. Relativizing everything in dat.ecosys to make it comparable to 
#        tree rings
#     c. Finding the percent change in NPP relative to the mean for that 
#        particular scale
# 5. Graphing & Analyzing Sensitivity
#     a. Graphing
#     b. Quantitative Analysis
# 6. Graphing & Analyzing Emulated Change in NPP
#     a. Graphing
#     b. Quantitative Analysis
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid); library(scales); 
library(gridExtra)
library(car); library(zoo)
# ----------------------------------------

# ----------------------------------------
# 1. Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
# setwd("..")
path.data <- "Data"
in.base <- "Data/gamms/Sensitivity_TempExtent"
out.dir <- "Data/analyses/analysis_TempExtent"
fig.dir <- "Figures/analyses/analysis_TempExtent"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# 2. Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]

# Colors for the rest of the script
models.use <- unique(ecosys[,"Model.Order"])
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray40"))

# Load the statistical model results
load(file.path(in.base, "gamm_TempExtent.Rdata"))

dat.ecosys <- cbind(mod.out$data[, ], mod.out$ci.response[,c("mean", "lwr", "upr")])
# wt.terms   <- mod.out$weights[,]
# ci.terms   <- mod.out$ci.terms[,]
# sim.terms  <- mod.out$sim.terms
# sim.terms$Effect <- as.factor(sim.terms$Effect)

# Grouping the kind and source of the data
dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRW", "RW", "NPP"))
# wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRW", "RW", "NPP"))
# ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRW", "RW", "NPP"))
# sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRW", "RW", "NPP"))

dat.ecosys$data.type <- as.factor(ifelse(substr(dat.ecosys$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# ci.terms  $data.type <- as.factor(ifelse(substr(ci.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))

summary(dat.ecosys)
# summary(wt.terms)
# summary(ci.terms)
# summary(sim.terms[,1:10])
# ----------------------------------------



# ----------------------------------------
# 3. Refit all gamms to the range of conditions used in the full runs
# 3.a. Redoing sensitivity
# 3.b. "Emulating" change in NPP based on sensitivity curves
# ----------------------------------------
source("R/0_Calculate_GAMM_Posteriors.R")
source("R/0_Calculate_GAMM_Weights.R")

# set up and re-run the terms posterior distribution calclations
{
extents <- paste(unique(dat.ecosys$Extent))

out.new <- list()

n=250
var.ci <- data.frame(tair    = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "tair"   ], na.rm=T),
                                   max(dat.ecosys[dat.ecosys$data.type=="Model", "tair"   ], na.rm=T),
                                   length.out=n),
                     precipf = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "precipf"], na.rm=T),
                                   max(dat.ecosys[dat.ecosys$data.type=="Model", "precipf"], na.rm=T),
                                   length.out=n),
                     CO2     = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "CO2"    ], na.rm=T),
                                   max(dat.ecosys[dat.ecosys$data.type=="Model", "CO2"    ], na.rm=T),
                                   length.out=n)
                     )
summary(var.ci)

# Making a climate series so we can emulate the ecosystem dynamics using a base series
#  Note: this will assume the same AGB trajectories; just changing the NPP
clim.subset <- (dat.ecosys$Model=="ed2" & dat.ecosys$Extent=="850-2010" & dat.ecosys$Resolution=="t.001")
var.ts      <- data.frame(Year    = dat.ecosys[clim.subset,"Year"   ], 
                          tair    = dat.ecosys[clim.subset,"tair"   ],
                          precipf = dat.ecosys[clim.subset,"precipf"],
                          CO2     = dat.ecosys[clim.subset,"CO2"    ]
                          )

# basically going through process.gamm here
for(m in unique(dat.ecosys$Model)){
	for(e in unique(dat.ecosys[dat.ecosys$Model==m, "Extent"])){
		# The year index for the gams
		yr <- strsplit(paste(e), "-")[[1]][1]
		if(nchar(yr)<4) yr = paste0(0, yr)

   		# Setting up some stuff
   		dat.subset <- dat.ecosys$Model==m & dat.ecosys$Extent==e & dat.ecosys$Resolution=="t.001"
		gam.now    <- mod.out[[paste0("gamm.", m, "_", yr, ".TempExtent")]]

		# A new blank list for everything to go in
		out.new[[paste(m, yr, sep="_")]] <- list(data=dat.ecosys[dat.subset, ], gamm=gam.now)

		# -----------------------
		# 3.a. Extrapolated Model Sensitivity 
		# -----------------------
		# Setting up a data frame for this model-extent combo
		ci.dat <- data.frame(Model      = as.factor(m), 
		                     Extent     = as.factor(e), 
		                     Resolution = as.factor("t.001"),
		                     Site       = as.factor(unique(dat.ecosys[dat.subset, "Site"  ])[1]),
		                     PlotID     = as.factor(unique(dat.ecosys[dat.subset, "PlotID"])[1]),
		                     TreeID     = as.factor(unique(dat.ecosys[dat.subset, "TreeID"])[1]),
		                     Time       = seq(min(dat.ecosys[dat.ecosys$Model==m, "Time"]), 
		                                      max(dat.ecosys[dat.ecosys$Model==m, "Time"]),
		                                      length.out=n
		                                      )
		                      )

		# Adding in the met data that has the full range of values
		ci.dat <- cbind(ci.dat, var.ci)

		# doing the sensitivity calculations
		terms.out   <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ci.dat, 
		                           vars=c("tair", "precipf", "CO2"), terms=T)

		 # Adding the new sims to the list
		out.new[[paste(m, yr, sep="_")]][["ci.terms"    ]] <- terms.out$ci
		out.new[[paste(m, yr, sep="_")]][["sim.terms"   ]] <- terms.out$sims
		# -----------------------

		# -----------------------
		# 3.a. "Emulated" NPP 
		# -----------------------
		# Only run this for models!!
		if(!unique(dat.ecosys[dat.ecosys$Model==m, "data.type"])=="Model") next
		
		# Running the gamm through our full time series "emulator"
		ts.dat  <- data.frame(Model      = as.factor(m), 
		                      Extent     = as.factor(e), 
		                      Resolution = as.factor("t.001"),
		                      Site       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Site"  ],
		                      PlotID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","PlotID"],
		                      TreeID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","TreeID"],
		                      Time       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Time"  ]
		                      )

		 #NOTE: Jules is missing 2010, so we're just going to pretend it's exactly the same as 2009
		 if(substr(m, 1, 5)=="jules"){
		 	ts.dat <- rbind(ts.dat, ts.dat[nrow(ts.dat),])
		 }

		 # Adding in the met data that has the full range of values
		 ts.dat <- cbind(ts.dat, var.ts)	
		 
		 # Doing the emulation of NPP & what drives it
		 pred.out    <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ts.dat, 
		                            vars=c("tair", "precipf", "CO2"), terms=F) 
		 weight.out  <- factor.weights(model.gam=gam.now, model.name=m, newdata=ts.dat,extent=e, 
		                          vars=c("tair", "precipf", "CO2")) 
		 
		 # Adding the new sims to the list
		 out.new[[paste(m, yr, sep="_")]][["weights"     ]] <- weight.out
		 out.new[[paste(m, yr, sep="_")]][["ci.response" ]] <- pred.out $ci
		 out.new[[paste(m, yr, sep="_")]][["sim.response"]] <- pred.out $sims

		# -----------------------
		
	}
}
} # End sensitivity & NPP emulation

# Binding things together into a single data frame
summary(out.new)
{
for(i in 1:length(out.new)){
	if(i == 1){
		dat.ecosys2 <- out.new[[i]]$ci.response
		ci.terms    <- out.new[[i]]$ci.terms
		wt.terms    <- out.new[[i]]$weights
	} else {
		dat.ecosys2 <- rbind(dat.ecosys2, out.new[[i]]$ci.response)
		ci.terms    <- rbind(ci.terms   , out.new[[i]]$ci.terms   )
		wt.terms    <- rbind(wt.terms   , out.new[[i]]$weights    )
	}
}
} # end binding

# Adding i key factors from the actual 850-2010 runs 
dim(dat.ecosys2)
ecosys.vars <- c("Model", "Model.Order", "Y.type", "data.type", "Site", "PlotID", "TreeID", "Year", "Y", "tair", "precipf", "CO2")
data.bind <- rbind(dat.ecosys[dat.ecosys$data.type=="Model"   & dat.ecosys$Extent== "850-2010", ecosys.vars],
				   dat.ecosys[dat.ecosys$Model=="TreeRingRW"  & dat.ecosys$Extent=="1901-2010", ecosys.vars],
				   dat.ecosys[dat.ecosys$Model=="TreeRingNPP" & dat.ecosys$Extent=="1985-2010", ecosys.vars])

dat.ecosys2 <- merge(dat.ecosys2, data.bind, all.x=T, all.y=T)

summary(dat.ecosys[dat.ecosys$Model==m,])
summary(dat.ecosys2[dat.ecosys2$Model==m,])
# ----------------------------------------

# summary(wt.terms)
# summary(ci.terms)
# summary(sim.terms[,1:10])

# ----------------------------------------
# 4. Standardize driver responses to the mean model NPP to facilitate comparisons
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
summary(dat.ecosys2)

# Make sure all data sets are ordered by year, then treeID, then plotID, then Model
sort.order <- c("Model", "PlotID", "TreeID", "Year")
dat.ecosys2 <- dat.ecosys2[order(dat.ecosys2$Model, dat.ecosys2$PlotID, dat.ecosys2$TreeID, dat.ecosys2$Year),]
wt.terms <- wt.terms[order(wt.terms$Model, wt.terms$PlotID, wt.terms$TreeID, wt.terms$Year),]

# Double Check to make sure things are sorted by year so rollapply works
dat.ecosys2[which(dat.ecosys2$Model=="TreeRingRW")[1:20],]
wt.terms  [which(wt.terms  $Model=="TreeRingRW")[1:20],]

{
for(m in unique(ci.terms$Model)){

		# -----------------------
		# 4.a. Find the NPP to relativize each set off of
		# Using mean model NPP across sites since the GAMM response curves are for 
		#    the whole model & not site-specific are parameterized
		# Note: We're doing this based on the whole-model mean, from the full-length
		#       runs
		# -----------------------
		# Find the appropriate reference extent for the model type

		npp <- mean(dat.ecosys2[dat.ecosys2$Model==m, "Y"], na.rm=T)			
		# -----------------------
		
		# -----------------------
		# 4.b Relativizing everything in dat.ecosys2 to make it comparable to tree rings
		# -----------------------
		{		
		# Which factors to relativize
		y.rel <- c("Y", "mean", "lwr", "upr")

		# for some reason, I can create multiple new columns at once
		# Solution: use a loop to create blank columns and then fill them
		for(y in y.rel){
			dat.ecosys2[dat.ecosys2$Model==m,paste0(y, ".rel"       )] <- NA	
			dat.ecosys2[dat.ecosys2$Model==m,paste0(y, ".10"        )] <- NA	
			dat.ecosys2[dat.ecosys2$Model==m,paste0(y, ".rel", ".10")] <- NA	
		}
		dat.ecosys2[dat.ecosys2$Model==m,paste0(y.rel, ".rel")] <- dat.ecosys2[dat.ecosys2$Model==m, y.rel]/npp
		
		# Getting 10-year running means to make clearer figures
		for(s in unique(dat.ecosys2[dat.ecosys2$Model==m, "Site"])){
			# Note: If we're working with tree ring data, we need to go by plot for NPP products 
			#       & by Tree for individual-level tree rings products
			if(m=="TreeRingNPP"){
				for(p in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s, "PlotID"])){

					# Raw NPP (to add dark line over faded annual wiggles)
					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# Relativized NPP (to have generalized patterns for figures)
					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s  & dat.ecosys2$PlotID==p, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				}
			} else if(m=="TreeRingRW") {
				for(t in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s, "TreeID"])){
					# If we have too few data points, we need to skip that tree 
					if(length(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t, y.rel[1]]) < 10) next

					# Raw NPP (to add dark line over faded annual wiggles)
					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# Relativized NPP (to have generalized patterns for figures)
					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s  & dat.ecosys2$TreeID==t, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				}
			} else {
				# Raw NPP (to add dark line over faded annual wiggles)
				dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

				# Relativized NPP (to have generalized patterns for figures)
				dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Site==s, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
			}
		}
		}		
		# -----------------------

		
		# -----------------------
		# 4.c. Finding the percent change in NPP relative to the mean for that particular scale
		# -----------------------
		{
		y.rel <- c("mean", "lwr", "upr")
		for(y in y.rel){
			ci.terms[ci.terms$Model==m,paste0(y, ".rel"       )] <- NA	
		}		

		ci.terms[ci.terms$Model==m,paste0(y.rel,".rel")] <- ci.terms[ci.terms $Model==m, y.rel]/npp
		
		# # Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
		# cols.sim <- which(substr(names(sim.terms),1,1)=="X")
		# sim.terms[sim.terms$Model==m,cols.sim] <- sim.terms[sim.terms$Model==m,cols.sim]/npp
		}
		# -----------------------

		# -----------------------
		# 4.d. Relativizing the factor fits through times and weights as well
		# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
		# -----------------------
		{
		y.rel <- c("fit.full", "fit.tair", "fit.precipf", "fit.CO2")
		# for some reason, I can create multiple new columns at once
		# Solution: use a loop to create blank columns and then fill them
		for(y in y.rel){
			wt.terms[wt.terms$Model==m,paste0(y, ".rel"       )] <- NA	
			wt.terms[wt.terms$Model==m,paste0(y, ".rel", ".10")] <- NA	
		}
		
		wt.terms[wt.terms$Model==m,paste0(y.rel, ".rel")] <- 1+(wt.terms[wt.terms$Model==m,y.rel])/npp
		
		# We only really care about smoothing the relativized weights
		y.rel2 <- c(paste0(y.rel, ".rel"), "weight.tair", "weight.precipf", "weight.CO2")
		for(y in y.rel2){
			wt.terms[wt.terms$Model==m,paste0(y, ".10")] <- NA	
		}

		# Getting 10-year running means to make clearer figures
		for(s in unique(dat.ecosys2[dat.ecosys2$Model==m, "Site"])){
			if(m=="TreeRingNPP"){
				for(p in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "PlotID"])){
					# Relativized NPP (to have generalized patterns for figures)
					wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
				}
			} else if(m=="TreeRingRW"){
				for(t in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "TreeID"])){

					# If we have too few data points, we need to skip that tree 
					if(length(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2[1]]) < 10) next

					# Relativized NPP (to have generalized patterns for figures)
					wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
				}				
			} else {
				# Relativized NPP (to have generalized patterns for figures)
				wt.terms[wt.terms$Model==m & wt.terms$Site==s,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)
			}
		}
		}
		# -----------------------

}
} # End section block


summary(dat.ecosys2)
summary(ci.terms)
summary(wt.terms)
# ----------------------------------------

# ----------------------------------------
# 5. Graphing & Analyzing Sensitivity
# ----------------------------------------
# -----------------------
# 5.a. Graphing
# -----------------------
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
ci.terms.graph[which(ci.terms.graph$mean.rel>0.75),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(0.75),"lwr.rel"] <- 0.75 
ci.terms.graph[ci.terms.graph$upr.rel>(0.75),"upr.rel"] <- 0.75 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15

ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
summary(ci.terms.graph)

# Creating a mask for values outside of the model drivers for that time period
for(e in unique(ci.terms.graph$Extent)){
	yr <- as.numeric(strsplit(paste(e), "-")[[1]][1])

	tair    <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"tair"   ], na.rm=T) - 273.15
	precipf <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"precipf"], na.rm=T)
	co2     <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"CO2"    ], na.rm=T)

	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="tair"   , "line.min"] <- tair   [1]
	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="precipf", "line.min"] <- precipf[1]
	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="CO2"    , "line.min"] <- co2    [1]

	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="tair"   , "line.max"] <- tair   [2]
	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="precipf", "line.max"] <- precipf[2]
	ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="CO2"    , "line.max"] <- co2    [2]
}
ci.terms.graph$x.min <- ifelse(ci.terms.graph$x<ci.terms.graph$line.min, ci.terms.graph$x, ci.terms.graph$line.min)
ci.terms.graph$x.max <- ifelse(ci.terms.graph$x>ci.terms.graph$line.max, ci.terms.graph$x, ci.terms.graph$line.max)
ci.terms.graph$mask.min <- min(ci.terms.graph$lwr.rel, na.rm=T)
ci.terms.graph$mask.max <- max(ci.terms.graph$upr.rel, na.rm=T)


pdf(file.path(fig.dir, "Fig4_Sensitivity_Rel_Extent.pdf"), height=11, width=8,5)
{
print(
ggplot(data=ci.terms.graph[!ci.terms.graph$Effect=="Time",]) + facet_grid(Extent~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
	# Lower Shaded Region
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	# Upper Shaded Region
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme_bw()
)
}
dev.off()
# -----------------------

# -----------------------
# 5.b. Analysis: 
# -----------------------
# -----------------------
# ----------------------------------------

# ----------------------------------------
# 6. Graphing & Analyzing Emulated Change in NPP
# ----------------------------------------
# -----------------------
# 6.a. Graphing
# -----------------------
# Fixing 2010 Jules
dat.ecosys2[substr(dat.ecosys2$Model,1,5)=="jules" & dat.ecosys2$Year==2010, c("Model.Order", "Y.type", "data.type", "Y", "Y.rel", "Y.10", "Y.rel.10")] <- dat.ecosys2[substr(dat.ecosys2$Model,1,5)=="jules" & dat.ecosys2$Year==2009, c("Model.Order", "Y.type", "data.type", "Y", "Y.rel", "Y.10", "Y.rel.10")]

wt.terms2 <- merge(wt.terms, dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year", "Y", "Y.rel", "Y.10", "Y.rel.10")], all.x=T, all.y=F)
summary(wt.terms2)


wt.terms2[is.na(wt.terms2$weight.tair.10), c("weight.tair.10", "weight.precipf.10", "weight.CO2.10")] <- 0
summary(wt.terms2)
# -----------
# Aggregating to the ensemble-level
# -----------
{
factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")

ensemble.wts1 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

ensemble.wts.lo <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], FUN=quantile, 0.025, na.rm=T)
names(ensemble.wts.lo)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[5:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], FUN=quantile, 0.975, na.rm=T)
names(ensemble.wts.hi)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[5:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

ensemble.wts.final <- cbind(ensemble.wts1, ensemble.wts.lo[5:ncol(ensemble.wts.lo)], ensemble.wts.hi[5:ncol(ensemble.wts.hi)])
summary(ensemble.wts.final)
}
# -----------


# -----------
# Re-normalizing factor weights
# -----------
{
wts.sum.10 <- abs(ensemble.wts.final$weight.tair.10) + abs(ensemble.wts.final$weight.precipf.10) + abs(ensemble.wts.final$weight.CO2.10)
ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")] <- ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")]/wts.sum.10
ensemble.wts.final[is.na(ensemble.wts.final$weight.tair.10   ),"weight.tair.10"   ] <- 0
ensemble.wts.final[is.na(ensemble.wts.final$weight.precipf.10),"weight.precipf.10"] <- 0
ensemble.wts.final[is.na(ensemble.wts.final$weight.CO2.10    ),"weight.CO2.10"    ] <- 0


wts.sum <- abs(ensemble.wts.final$weight.tair) + abs(ensemble.wts.final$weight.precipf) + abs(ensemble.wts.final$weight.CO2)
ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")] <- ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")]/wts.sum
summary(ensemble.wts.final)
}
# -----------

summary(ensemble.wts.final)


# -----------
# The figure
# -----------
pdf(file.path(fig.dir, "Fig5_Ensemble_Drivers_Time_PHO_1900-2010_Annual.pdf"), width=11, height=8.5)
{
print(
ggplot(ensemble.wts.final) + facet_grid(Extent~., scales="free_y") +
 	geom_ribbon(data= ensemble.wts.final[,], aes(x=Year, ymin=fit.full.rel.lo*100, ymax=fit.full.rel.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, y=fit.full.rel*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1985-2010","weight.tair"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1985-2010","weight.CO2"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1985-2010","weight.precipf"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf"])), size=3) +
 	geom_hline(y=200, linetype="dashed") +
	scale_x_continuous(limits=c(1900,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts.final$Year), -2), round(max(ensemble.wts.final$Year), -2), by=100)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(size=rel(1), color="black"),
		  axis.text.y=element_text(size=rel(1), color="black"), 
		  axis.title.x=element_text(size=rel(1), face="bold"),  
		  axis.title.y=element_text(size=rel(1), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
dev.off()
# -----------

# -----------------------

# -----------------------
# 6.b. Analysis: 
# -----------------------
# -----------------------
# ----------------------------------------
