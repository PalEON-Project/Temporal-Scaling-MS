# ----------------------------------------
# Objective: Compare model/tree ring sensitivites by Temporal Extent
# Christy Rollinson, crollinson@gmail.com
# Date Created: January 2016
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
library(nlme)
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
{
load(file.path(path.data, "EcosysData.Rdata"))
# ecosys <- ecosys[!ecosys$Model=="linkages",]

# Colors for the rest of the script
models.use <- unique(ecosys[,"Model.Order"])
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray40"))

# Load the statistical model results
load(file.path(in.base, "gamm_TempExtent.Rdata"))

# dat.ecosys <- cbind(mod.out$data[, ], mod.out$ci.response[,c("mean", "lwr", "upr")])
# # wt.terms   <- mod.out$weights[,]
# # ci.terms   <- mod.out$ci.terms[,]
# # sim.terms  <- mod.out$sim.terms
# # sim.terms$Effect <- as.factor(sim.terms$Effect)
# 
# # Grouping the kind and source of the data
# dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRW", "RW", "NPP"))
# # wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRW", "RW", "NPP"))
# # ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRW", "RW", "NPP"))
# # sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRW", "RW", "NPP"))
# 
# dat.ecosys$data.type <- as.factor(ifelse(substr(dat.ecosys$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# # wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# # ci.terms  $data.type <- as.factor(ifelse(substr(ci.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# # sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# 
# summary(dat.ecosys)
# # summary(wt.terms)
# # summary(ci.terms)
# # summary(sim.terms[,1:10])
}
# ----------------------------------------



# ----------------------------------------
# 3. Refit all gamms to the range of conditions used in the full runs
# 3.a. Redoing sensitivity
# 3.b. "Emulating" change in NPP based on sensitivity curves
# ----------------------------------------
{
# source("R/0_Calculate_GAMM_Posteriors.R")
# source("R/0_Calculate_GAMM_Weights.R")
# 
# # set up and re-run the terms posterior distribution calclations
# {
# extents <- paste(unique(dat.ecosys$Extent))
# 
# out.new <- list()
# 
# n=250
# var.ci <- data.frame(tair    = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "tair"   ], na.rm=T),
#                                    max(dat.ecosys[dat.ecosys$data.type=="Model", "tair"   ], na.rm=T),
#                                    length.out=n),
#                      precipf = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "precipf"], na.rm=T),
#                                    max(dat.ecosys[dat.ecosys$data.type=="Model", "precipf"], na.rm=T),
#                                    length.out=n),
#                      CO2     = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "CO2"    ], na.rm=T),
#                                    max(dat.ecosys[dat.ecosys$data.type=="Model", "CO2"    ], na.rm=T),
#                                    length.out=n),
#                      Year    = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "Year"   ], na.rm=T),
#                                    max(dat.ecosys[dat.ecosys$data.type=="Model", "Year"   ], na.rm=T),
#                                    length.out=n)
#                      # Time    = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "Year"   ], na.rm=T),
#                                    # max(dat.ecosys[dat.ecosys$data.type=="Model", "Year"   ], na.rm=T),
#                                    # length.out=n)
#                      # Biomass = seq(min(dat.ecosys[dat.ecosys$data.type=="Model", "AGB"   ], na.rm=T),
#                                    # max(dat.ecosys[dat.ecosys$data.type=="Model", "AGB"   ], na.rm=T),
#                                    # length.out=n)
#                      )
# summary(var.ci)
# 
# # Making a climate series so we can emulate the ecosystem dynamics using a base series
# #  Note: this will assume the same AGB trajectories; just changing the NPP
# clim.subset <- (dat.ecosys$Model=="ed2" & dat.ecosys$Extent=="850-2010" & dat.ecosys$Resolution=="t.001")
# var.ts      <- data.frame(Site    = dat.ecosys[clim.subset,"Site"   ], 
#                           PlotID  = dat.ecosys[clim.subset,"PlotID" ], 
#                           Year    = dat.ecosys[clim.subset,"Year"   ], 
#                           # Time    = dat.ecosys[clim.subset,"Year"   ], 
#                           # Biomass = dat.ecosys[clim.subset,"AGB"    ], 
#                           tair    = dat.ecosys[clim.subset,"tair"   ],
#                           precipf = dat.ecosys[clim.subset,"precipf"],
#                           CO2     = dat.ecosys[clim.subset,"CO2"    ]
#                           )
# 
# # basically going through process.gamm here
# for(m in unique(dat.ecosys$Model)){
# 	for(e in unique(dat.ecosys[dat.ecosys$Model==m, "Extent"])){
# 		# The year index for the gams
# 		yr <- strsplit(paste(e), "-")[[1]][1]
# 		if(nchar(yr)<4) yr = paste0(0, yr)
# 
#    	# Setting up some stuff
#    	dat.subset <- dat.ecosys$Model==m & dat.ecosys$Extent==e & dat.ecosys$Resolution=="t.001"
# 		gam.now    <- mod.out[[paste0("gamm.", m, "_", yr, ".TempExtent")]]
# 
# 		# A new blank list for everything to go in
# 		out.new[[paste(m, yr, sep="_")]] <- list(data=dat.ecosys[dat.subset, ], gamm=gam.now)
# 
# 		# -----------------------
# 		# 3.a. Extrapolated Model Sensitivity 
# 		# -----------------------
# 		# Setting up a data frame for this model-extent combo
# 		ci.dat <- data.frame(Model      = as.factor(m), 
# 		                     Extent     = as.factor(e), 
# 		                     Resolution = as.factor("t.001"),
# 		                     Site       = as.factor(unique(dat.ecosys[dat.subset, "Site"  ])[1]),
# 		                     PlotID     = as.factor(unique(dat.ecosys[dat.subset, "PlotID"])[1]),
# 		                     TreeID     = as.factor(unique(dat.ecosys[dat.subset, "TreeID"])[1]),
# 		                     Time       = mean(dat.ecosys[dat.ecosys$Model==m, "Time"]), 
# 		                     Biomass    = seq(min(dat.ecosys[dat.ecosys$Model==m, "Biomass"], na.rm=T), 
# 		                                      max(dat.ecosys[dat.ecosys$Model==m, "Biomass"], na.rm=T), 
# 		                                      length.out=n)
# 		                      )
# 
# 		# Adding in the met data that has the full range of values
# 		ci.dat <- cbind(ci.dat, var.ci)
# 
# 		# doing the sensitivity calculations
# 		terms.out   <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ci.dat, 
# 		                           vars=c("tair", "precipf", "CO2", "Biomass"), terms=T)
# 
# 		# Adding the new sims to the list
# 		out.new[[paste(m, yr, sep="_")]][["ci.terms"    ]] <- terms.out$ci
# 		out.new[[paste(m, yr, sep="_")]][["sim.terms"   ]] <- terms.out$sims
# 		# -----------------------
# 
# 		# -----------------------
# 		# 3.a. "Emulated" NPP 
# 		# -----------------------
# 		# Only run this for models!!
# 		if(!unique(dat.ecosys[dat.ecosys$Model==m, "data.type"])=="Model") next
# 		
# 		# Running the gamm through our full time series "emulator"
# 		ts.dat  <- data.frame(Model      = as.factor(m), 
# 		                      Extent     = as.factor(e), 
# 		                      Resolution = as.factor("t.001"),
# 		                      Site       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Site"   ],
# 		                      PlotID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","PlotID" ],
# 		                      TreeID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","TreeID" ],
# 		                      Time       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Time"   ],
# 		                      Biomass    = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Biomass"],
# 		                      Year       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Year"   ],
# 		                      Y          = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Y"      ]
# 		                      )
# 		 
# 		 # Adding in the met data that has the full range of values
# 		 ts.dat <- merge(ts.dat, var.ts, all.x=T, all.y=T)	
# 
# 		 #NOTE: Jules is missing 2010, so we're just going to pretend it's exactly the same as 2009
# 		 if(substr(m, 1, 5) %in% c("jules", "linka") ){
# 		 	ts.dat[ts.dat$Year==2010,c("Y", "Time", "Biomass")] <- ts.dat[ts.dat$Year==2009,c("Y", "Time", "Biomass")]
# 		 }
# 
# 
# 		 # Doing the emulation of NPP & what drives it
# 		 pred.out    <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ts.dat, 
# 		                            vars=c("tair", "precipf", "CO2", "Year", "Time", "Biomass"), terms=F) 
# 		 weight.out  <- factor.weights(model.gam=gam.now, model.name=m, newdata=ts.dat,extent=e, 
# 		                          vars=c("tair", "precipf", "CO2", "Year", "Time", "Biomass")) 
# 		 
# 		 # Adding the new sims to the list
# 		 out.new[[paste(m, yr, sep="_")]][["weights"     ]] <- weight.out
# 		 out.new[[paste(m, yr, sep="_")]][["ci.response" ]] <- pred.out $ci
# 		 out.new[[paste(m, yr, sep="_")]][["sim.response"]] <- pred.out $sims
# 
# 		# -----------------------
# 		
# 	}
# }
# } # End sensitivity & NPP emulation
# 
# # Binding things together into a single data frame
# summary(out.new)
# {
# for(i in 1:length(out.new)){
# 	if(i == 1){
# 		dat.ecosys2 <- out.new[[i]]$data
# 		mod.resp    <- out.new[[i]]$ci.response
# 		ci.terms    <- out.new[[i]]$ci.terms
# 		sim.terms   <- out.new[[i]]$sim.terms
# 		wt.terms    <- out.new[[i]]$weights
# 	} else {
# 		dat.ecosys2 <- rbind(dat.ecosys2, out.new[[i]]$data       )
# 		mod.resp    <- rbind(mod.resp   , out.new[[i]]$ci.response)
# 		ci.terms    <- rbind(ci.terms   , out.new[[i]]$ci.terms   )
# 		sim.terms   <- rbind(sim.terms  , out.new[[i]]$sim.terms  )
# 		wt.terms    <- rbind(wt.terms   , out.new[[i]]$weights    )
# 	}
# }
# } # end binding
# 
# summary(dat.ecosys2)
# summary(mod.resp)
# 
# # combining the model emulation output (mod.resp) witht he inputs that didn't get transferred (dat.ecosys2) 
# dim(dat.ecosys2)
# ecosys.vars <- c("Model", "Model.Order", "data.type", "Y.type", "Site", "PlotID", "TreeID", "Y", "Time", "Year")
# 
# 
# dat.ecosys2a <- merge(mod.resp, dat.ecosys2[dat.ecosys2$data.type=="Model" & dat.ecosys2$Extent=="850-2010", ecosys.vars], all.x=F, all.y=F)
# dat.ecosys2b <- merge(dat.ecosys2a, dat.ecosys2[dat.ecosys2$Model=="TreeRingRW" & dat.ecosys2$Extent=="1901-2010", c(ecosys.vars, "tair", "precipf", "CO2", "Extent", "Resolution")], all.x=T, all.y=T)
# dat.ecosys2c <- merge(dat.ecosys2b, dat.ecosys2[dat.ecosys2$Model=="TreeRingNPP" & dat.ecosys2$Extent=="1980-2010", c(ecosys.vars, "tair", "precipf", "CO2", "Extent", "Resolution")], all.x=T, all.y=T)
# summary(dat.ecosys2c)
# 
# dat.ecosys2 <- dat.ecosys2c
# 
# summary(dat.ecosys[dat.ecosys$Model==m,])
# summary(dat.ecosys2[dat.ecosys2$Model==m,])
}
# ----------------------------------------

# summary(wt.terms)
# summary(ci.terms)
# summary(sim.terms[,1:10])

# ----------------------------------------
# 4. Standardize driver responses to the mean model NPP to facilitate comparisons
# ----------------------------------------
{
# # Across all scales (resolution) finding the mean NPP
# # NOTE: we ARE relativizing per site here since the response curves were site-specific
# summary(dat.ecosys2)
# 
# # Make sure all data sets are ordered by year, then treeID, then plotID, then Model
# sort.order <- c("Model", "PlotID", "TreeID", "Year")
# dat.ecosys2 <- dat.ecosys2[order(dat.ecosys2$Model, dat.ecosys2$PlotID, dat.ecosys2$TreeID, dat.ecosys2$Year),]
# wt.terms <- wt.terms[order(wt.terms$Model, wt.terms$PlotID, wt.terms$TreeID, wt.terms$Year),]
# 
# # Double Check to make sure things are sorted by year so rollapply works
# dat.ecosys2[which(dat.ecosys2$Model=="TreeRingRW")[1:20],]
# # wt.terms  [which(wt.terms  $Model=="TreeRingRW")[1:20],]
# 
# wt.terms <- wt.terms[!is.na(wt.terms$Resolution),]
# 
# {
# for(m in unique(ci.terms$Model)){
# 		for(e in unique(ci.terms[ci.terms$Model==m, "Extent"])){
# 		# -----------------------
# 		# 4.a. Find the NPP to relativize each set off of
# 		# Using mean model NPP across sites since the GAMM response curves are for 
# 		#    the whole model & not site-specific are parameterized
# 		# Note: We're doing this based on the whole-model mean, from the full-length
# 		#       runs
# 		# -----------------------
# 		# Find the appropriate reference extent for the model type
# 		if(m=="TreeRingRW") ext="1901-2010" else if(m=="TreeRingNPP") ext="1980-2010" else ext="850-2010"
# 		# ext="1980-2010"
# 		yr <- as.numeric(strsplit(paste(e), "-")[[1]][1])
# 
# 		npp <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==ext & dat.ecosys2$Year>=yr, "Y"], na.rm=T)			
# 		# npp.fit <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==ext, "mean"], na.rm=T)			
# 		# -----------------------
# 		
# 		# -----------------------
# 		# 4.b Relativizing everything in dat.ecosys2 to make it comparable to tree rings
# 		# -----------------------
# 		{		
# 		# Which factors to relativize
# 		y.rel <- c("Y", "mean", "lwr", "upr")
# 
# 		# for some reason, I can create multiple new columns at once
# 		# Solution: use a loop to create blank columns and then fill them
# 		for(y in y.rel){
# 			dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e,paste0(y, ".rel"       )] <- NA	
# 			dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e,paste0(y, ".10"        )] <- NA	
# 			dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e,paste0(y, ".rel", ".10")] <- NA	
# 		}
# 		# dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e,paste0("Y", ".rel")] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e, "Y"]/npp
# 		dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e,paste0(y.rel, ".rel")] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e, y.rel]/npp
# 		
# 		# Getting 10-year running means to make clearer figures
# 		for(s in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e, "Site"])){
# 			# Note: If we're working with tree ring data, we need to go by plot for NPP products 
# 			#       & by Tree for individual-level tree rings products
# 			if(m=="TreeRingNPP"){
# 				for(p in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s, "PlotID"])){
# 
# 					# Raw NPP (to add dark line over faded annual wiggles)
# 					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$PlotID==p,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s  & dat.ecosys2$PlotID==p, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 				}
# 			} else if(m=="TreeRingRW") {
# 				for(t in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s, "TreeID"])){
# 					# If we have too few data points, we need to skip that tree 
# 					if(length(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t, y.rel[1]]) < 10) next
# 
# 					# Raw NPP (to add dark line over faded annual wiggles)
# 					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s & dat.ecosys2$TreeID==t,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s  & dat.ecosys2$TreeID==t, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 				}
# 			} else {
# 				# Raw NPP (to add dark line over faded annual wiggles)
# 				dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 				# Relativized NPP (to have generalized patterns for figures)
# 				dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Site==s, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 			}
# 		}
# 		}		
# 		# -----------------------
# 
# 		
# 		# -----------------------
# 		# 4.c. Finding the percent change in NPP relative to the mean for that particular scale
# 		# -----------------------
# 		{
# 		y.rel <- c("mean", "lwr", "upr")
# 		for(y in y.rel){
# 			ci.terms[ci.terms$Model==m & ci.terms$Extent==e,paste0(y, ".rel"       )] <- NA	
# 		}		
# 
# 		ci.terms[ci.terms$Model==m & ci.terms$Extent==e,paste0(y.rel,".rel")] <- ci.terms[ci.terms$Model==m  & ci.terms$Extent==e, y.rel]/npp
# 		
# 		# Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
# 		cols.sim <- which(substr(names(sim.terms),1,1)=="X")
# 		sim.terms[sim.terms$Model==m,cols.sim] <- sim.terms[sim.terms$Model==m,cols.sim]/npp
# 		}
# 		# -----------------------
# 
# 		# -----------------------
# 		# 4.d. Relativizing the factor fits through times and weights as well
# 		# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
# 		# -----------------------
# 		{
# 		# Doing NPP calculation again for wt.terms because it does wonky things
# 		npp <- mean(wt.terms[wt.terms$Model==m & wt.terms$Extent==ext, "fit.full"], na.rm=T)			
# 
# 		y.rel <- c("fit.full", "fit.tair", "fit.precipf", "fit.CO2")
# 		# for some reason, I can create multiple new columns at once
# 		# Solution: use a loop to create blank columns and then fill them
# 		for(y in y.rel){
# 			wt.terms[wt.terms$Model==m & wt.terms$Extent==e,paste0(y, ".rel"       )] <- NA	
# 			wt.terms[wt.terms$Model==m & wt.terms$Extent==e,paste0(y, ".rel", ".10")] <- NA	
# 		}
# 		wt.terms[wt.terms$Model==m & wt.terms$Extent==e,paste0("fit.full", ".rel")] <- (wt.terms[wt.terms$Model==m & wt.terms$Extent==e,"fit.full"])/npp
# 		
# 		wt.terms[wt.terms$Model==m & wt.terms$Extent==e,paste0(y.rel[2:length(y.rel)], ".rel")] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Extent==e,y.rel[2:length(y.rel)]])/npp
# 		
# 		# We only really care about smoothing the relativized weights
# 		y.rel2 <- c(y.rel, paste0(y.rel, ".rel"), "weight.tair", "weight.precipf", "weight.CO2")
# 		for(y in y.rel2){
# 			wt.terms[wt.terms$Model==m & wt.terms$Extent==e,paste0(y, ".10")] <- NA	
# 		}
# 
# 		# Getting 10-year running means to make clearer figures
# 		for(s in unique(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e, "Site"])){
# 			if(m=="TreeRingNPP"){
# 				for(p in unique(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s, "PlotID"])){
# 					# Relativized NPP (to have generalized patterns for figures)
# 					wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s & wt.terms$PlotID==p,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s & wt.terms$PlotID==p, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
# 				}
# 			} else if(m=="TreeRingRW"){
# 				for(t in unique(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s, "TreeID"])){
# 
# 					# If we have too few data points, we need to skip that tree 
# 					if(length(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2[1]]) < 10) next
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s & wt.terms$TreeID==t,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
# 				}				
# 			} else {
# 				# Relativized NPP (to have generalized patterns for figures)
# 				wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Extent==e & wt.terms$Site==s, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 			}
# 		}
# 		}
# 		# -----------------------
# 
# 		} # End Extent Loop
# 	} # End Model Loop
# } # End section block
# 
# 
# summary(dat.ecosys2)
# summary(ci.terms)
# summary(wt.terms)
# 
# save(dat.ecosys, dat.ecosys2, ci.terms, sim.terms, wt.terms, file=file.path(out.dir, "post-process_TempExtent.RData"))
}
# ----------------------------------------

load(file.path(out.dir, "post-process_TempExtent.RData"))
summary(dat.ecosys2)
summary(ci.terms)
summary(sim.terms[,1:15])
summary(wt.terms)

# ----------------------------------------
# 5. Graphing & Analyzing Sensitivity
# ----------------------------------------
{
# -----------------------
# 5.a. Graphing
# -----------------------
{
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[!is.na(ci.terms.graph$mean.rel) & ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
ci.terms.graph[!is.na(ci.terms.graph$lwr.rel) & ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.00),"mean.rel"] <- NA 
ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$lwr.rel>(1.00),"lwr.rel"] <- 1.00
ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$upr.rel>(1.00),"upr.rel"] <- 1.00 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15

ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
summary(ci.terms.graph)

# Grouping the kind and source of the data
ci.terms.graph$Y.type <- as.factor(ifelse(ci.terms.graph$Model=="TreeRingRW", "RW", "NPP"))
ci.terms.graph$data.type <- as.factor(ifelse(substr(ci.terms.graph$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
summary(ci.terms.graph)
# summary(ci.terms.graph[ci.terms.graph$Model=="linkages",])

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

# Playing with the extent labels a bit so that "850-2010" is "all data" and "1901-2010" is left alone
ci.terms.graph$Extent2 <- as.factor(ifelse(ci.terms.graph$Extent=="850-2010" | 
                                             (ci.terms.graph$Extent=="1980-2010" & ci.terms.graph$Model=="TreeRingNPP") |
                                             (ci.terms.graph$Extent=="1901-2010" & ci.terms.graph$Model=="TreeRingRW")
                                           , "All Data", paste(ci.terms.graph$Extent)))

ci.terms.graph[ci.terms.graph$Extent2=="All Data", c("x.min", "x.max", "mask.min", "mask.max", "line.min", "line.max")] <- NA
summary(ci.terms.graph)


tree.rings.1901 <- ci.terms.graph[ci.terms.graph$Model=="TreeRingRW" & ci.terms.graph$Extent=="1901-2010",]
tree.rings.1901$Extent2 <- as.factor("1901-2010")
summary(tree.rings.1901)

tree.rings.npp <- ci.terms.graph[ci.terms.graph$Model=="TreeRingNPP" & ci.terms.graph$Extent=="1980-2010",]
tree.rings.npp$Extent2 <- as.factor("1980-2010")
summary(tree.rings.npp)


ci.terms.graph <- rbind(ci.terms.graph, tree.rings.1901, tree.rings.npp)
ci.terms.graph$Extent3 <- recode(ci.terms.graph$Extent2, "'All Data'='0'; '1901-2010'='1'; '1980-2010'='2'")
levels(ci.terms.graph$Extent3) <- c("0850-2010", "1901-2010", "1980-2010")
summary(ci.terms.graph)

ci.terms.graph <- ci.terms.graph[!(ci.terms.graph$Extent3=="0850-2010" & ci.terms.graph$data.type=="Tree Rings"),]
levels(ci.terms.graph$Effect) <- c("Temperature", "Precipitation", "CO2", "Biomass")
summary(ci.terms.graph)

pdf(file.path(fig.dir, "Sensitivity_Rel_Extent_AllExtents.pdf"), height=8.5, width=11)
{
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("Temperature", "Precipitation", "CO2") ,]) + facet_grid(Extent3~Effect, scales="free_x") +
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

pdf(file.path(fig.dir, "Sensitivity_Rel_Extent_Biomass.pdf"), height=11, width=8.5)
{
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect=="Biomass" ,]) + facet_wrap(~Model, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Extent), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Extent, linetype=Extent), size=1) +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	# scale_fill_manual(values=colors.use) +
	# scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme_bw()
)
}
dev.off()

fig3.tair <- {
  ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "Temperature",]) + 
  facet_grid(Extent3~Effect, scales="free_x") +
  geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
  geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
  # Lower Shaded Region
  geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
  geom_vline(aes(xintercept=line.min), linetype="dashed") +
  # Upper Shaded Region
  geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
  geom_vline(aes(xintercept=line.max), linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name=expression(bold(paste("Temperature ("^"o", "C)"))), breaks=c(10, 12.5, 15.0, 17.5)) +
  scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
  guides(fill=F, color=F, linetype=F) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
  theme(strip.text.x=element_text(size=12, face="bold"),
        strip.text.y=element_blank()) + 
  theme(axis.line=element_line(color="black", size=0.5), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_rect(fill=NA, color="black", size=0.5), 
        panel.background=element_blank(),
        panel.margin.x=unit(0, "lines"),
        panel.margin.y=unit(0, "lines"))  +
  theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
        axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
        axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
        axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.65,0,0,0), "lines")),
        axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"))
  
}
fig3.precip <- {
  ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "Precipitation",]) + 
    facet_grid(Extent3~Effect, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
    geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
    # Lower Shaded Region
    geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
    geom_vline(aes(xintercept=line.min), linetype="dashed") +
    # Upper Shaded Region
    geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
    geom_vline(aes(xintercept=line.max), linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name=expression(bold(paste("Precipitation (mm yr"^"-1", ")")))) +
    scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
    guides(fill=F, color=F, linetype=F) +
    scale_fill_manual(values=colors.use) +
    scale_color_manual(values=colors.use) +
    scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
    theme(strip.text.x=element_text(size=12, face="bold"),
          strip.text.y=element_blank()) + 
    theme(axis.line=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(fill=NA, color="black", size=0.5), 
          panel.background=element_blank(),
          panel.margin.x=unit(0, "lines"),
          panel.margin.y=unit(0, "lines"))  +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
          axis.title.y=element_blank(),  
          axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
          axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0.0,0.5,0.5), "lines"))
  
}
fig3.co2 <- {
  ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "CO2",]) + 
    facet_grid(Extent3~Effect, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
    geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
    # Lower Shaded Region
    geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
    geom_vline(aes(xintercept=line.min), linetype="dashed") +
    # Upper Shaded Region
    geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
    geom_vline(aes(xintercept=line.max), linetype="dashed") +
    scale_x_continuous(expand=c(0,0), name=expression(bold(paste("CO" ["2"], " (ppm)")))) +
    scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
    guides(fill=guide_legend(title="Model"), 
           color=guide_legend(title="Model"), 
           linetype=guide_legend(title="Model")) +
    scale_fill_manual(values=colors.use) +
    scale_color_manual(values=colors.use) +
    scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
    theme(legend.title=element_text(size=12, face="bold"),
          legend.text=element_text(size=10),
          legend.key=element_blank(),
          legend.key.size=unit(1.5, "lines"),
          legend.background=element_blank()) +
    theme(strip.text.x=element_text(size=12, face="bold"),
          strip.text.y=element_text(size=12, face="bold")) + 
    theme(axis.line=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(fill=NA, color="black", size=0.5), 
          panel.background=element_blank(),
          panel.margin.x=unit(0, "lines"),
          panel.margin.y=unit(0, "lines"))  +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
          axis.title.y=element_blank(),  
          axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.79,0,0,0), "lines")),
          axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"))
  
}

png(file.path(fig.dir, "Fig3_Sensitivity_Rel_extent.png"), width=8, height=6, units="in", res=180)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3, widths=c(1.3,1,2))))
print(fig3.tair  , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(fig3.precip, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
print(fig3.co2   , vp = viewport(layout.pos.row = 1, layout.pos.col=3))
dev.off()
}
# -----------------------

# -----------------------
# 5.b. Analysis: 
# Table 3 has the Effects of key components across scales
# -----------------------
# general summary stat for change in sensitivity across cales and effect

table3 <- data.frame()
{
  source("R/0_Calculate_GAMM_Posteriors.R")
  library(mgcv)
# --------
# 5.b.0. Adding some model-level stats to compare the relative sensitivities
# --------
summary(ci.terms)
{
# 1. Quantify with-in model sensitivity shifts due to change in temporal extent
{
# i. classify met by which extent it's observed in 
for(v in c("tair", "precipf", "CO2")){ # go through met (same for all models/data)
  range.1980 <- range(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Model", v], na.rm=T)
  range.1901 <- range(dat.ecosys2[dat.ecosys2$Year>=1901 & dat.ecosys2$data.type=="Model", v], na.rm=T)
    
  met.fill <-  as.factor(ifelse(ci.terms[ci.terms$Effect==v,"x"]>=range.1980[1] & 
                                  ci.terms[ci.terms$Effect==v,"x"]<=range.1980[2],
                                "1980-2010", 
                                ifelse(ci.terms[ci.terms$Effect==v,"x"]>=range.1901[1] & 
                                         ci.terms[ci.terms$Effect==v,"x"]<=range.1901[2],
                                       "1901-2010",
                                       "850-2010"
                                ))
                            )  
  ci.terms[ci.terms$Effect==v,"Met"] <- as.factor(met.fill)  
}
summary(ci.terms)
  
# ii.  Adjusting the sensitivity so that everywhere has the same x-intercept as 1980-2010
for(v in c("tair", "precipf", "CO2")){
  if(v=="CO2") r=0 else if(v=="tair") r=1 else r=-1
  
  for(m in unique(ci.terms$Model)){
    var.mid <- mean(ci.terms[ci.terms$Effect==v & ci.terms$Met=="1980-2010", "x"], na.rm=T)
    offset.1980 <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="1980-2010","mean.rel"])
    offset.1901 <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="1901-2010","mean.rel"])
    offset.850  <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="850-2010","mean.rel"])
  
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","mean.rel"] - offset.1980
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","mean.rel"] - offset.1901
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"mean.rel"] - offset.850

    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","upr.rel"] - offset.1980
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","upr.rel"] - offset.1901
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"upr.rel"] - offset.850
    
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","lwr.rel"] - offset.1980
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","lwr.rel"] - offset.1901
    ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"lwr.rel"] - offset.850
    
  }  
}


# Creating a mask for values outside of the model drivers for that time period
for(e in unique(ci.terms$Extent)){
  yr <- as.numeric(strsplit(paste(e), "-")[[1]][1])
  
  tair    <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"tair"   ], na.rm=T) 
  precipf <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"precipf"], na.rm=T)
  co2     <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"CO2"    ], na.rm=T)
  
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="tair"   , "line.min"] <- tair   [1]
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="precipf", "line.min"] <- precipf[1]
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="CO2"    , "line.min"] <- co2    [1]
  
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="tair"   , "line.max"] <- tair   [2]
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="precipf", "line.max"] <- precipf[2]
  ci.terms[ci.terms$Extent==e & ci.terms$Effect=="CO2"    , "line.max"] <- co2    [2]
}
ci.terms$x.min <- ifelse(ci.terms$x<ci.terms$line.min, ci.terms$x, ci.terms$line.min)
ci.terms$x.max <- ifelse(ci.terms$x>ci.terms$line.max, ci.terms$x, ci.terms$line.max)
ci.terms$mask.min <- min(ci.terms$lwr.rel, na.rm=T)
ci.terms$mask.max <- max(ci.terms$upr.rel, na.rm=T)
# ci.terms$mask.min.cent <- min(ci.terms$mean.rel.cent, na.rm=T)
# ci.terms$mask.max.cent <- max(ci.terms$mean.rel.cent, na.rm=T)
ci.terms$mask.min.cent <- -0.75
ci.terms$mask.max.cent <- 1

summary(ci.terms)


pdf(file.path(fig.dir, "Sensitivitity_Extents_Centered_1980.pdf"), height=11, width=8.5)
ggplot(data=ci.terms[ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
  facet_grid(Extent ~ Effect, scales="free_x") +
  geom_line(aes(x=x, y=mean.rel.cent, color=Model), size=2) +
  # Lower Shaded Region
  geom_ribbon(aes(x=x.min, ymin=mask.min.cent, ymax=mask.max.cent), alpha=0.3) +
  geom_vline(aes(xintercept=line.min), linetype="dashed") +
  # Upper Shaded Region
  geom_ribbon(aes(x=x.max, ymin=mask.min.cent, ymax=mask.max.cent), alpha=0.3) +
  geom_vline(aes(xintercept=line.max), linetype="dashed") +
  scale_x_continuous(expand=c(0,0), name="") +
  scale_y_continuous(limits=c(-0.75,1), expand=c(0,0)) +
  scale_color_manual(values=colors.use) +
  theme_bw()
dev.off()
}

# 2. Condensing model variability across space and time to get general model characteristics
{
  vars.agg <- c("Y", "Y.rel", "Biomass")
  mod.agg                           <- aggregate(dat.ecosys2[,vars.agg], 
                                                 by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=mean)
  mod.agg[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys2[,vars.agg], 
                                                 by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=sd)[,vars.agg]
  mod.agg
  
  
  # Finding the change in key variables in the modern era
  for(v in vars.agg){
    mod.agg[mod.agg$Model==m,paste0("dModern.", v)] <- NA  
  }
  
  for(m in unique(mod.agg$Model)){
    mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=1990 & dat.ecosys2$Year<=2010, vars.agg], na.rm=T) -
      colMeans(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=1830 & dat.ecosys2$Year<=1850, vars.agg], na.rm=T)
  }
  mod.agg
  
  # Add in the vegetation scheme
  mod.agg$veg.scheme <- as.factor(ifelse(mod.agg$Model %in% c("clm.bgc", "clm.cn", "sibcasa", "jules.stat"), "Static", "Dynamic"))
  
  # Finding mean fire return
  {
  output.all <- read.csv("../../phase1a_output_variables/PalEON_MIP_Yearly.csv")
  summary(output.all)
  
  agg.fire <- aggregate(output.all[,c("Fire")], by=output.all[,c("Model", "Updated")], FUN=mean, na.rm=T)
  names(agg.fire)[3] <- "Fire"
  agg.fire[is.na(agg.fire$Fire),"Fire"] <- 0
  agg.fire$Fire <- agg.fire$Fire*(60*60*24*365)*(100*100)*1e-3
  agg.fire
  
  mod.agg <- merge(mod.agg, agg.fire[,c("Model", "Fire")], all.x=T, all.y=F)
  mod.agg$Fire <- mod.agg$Fire # changing fire from KgC/m2/s to MgC/HA/yr (same units as NPP)
  mod.agg$fire.scheme <- as.factor(ifelse(mod.agg$Fire>0, "Yes", "No"))
  
  summary(mod.agg)
  }
}

ci.terms <- merge(ci.terms, mod.agg, all.x=T, all.y=F)
summary(ci.terms)


# 3. Getting the variability of slow processes at each time scale to use as explanatory factors
#    Slow Processes = Composition, Biomass
# NOTE: This overwrites a lot of what we already calcualted in mod.ag
{
# --------
# Need to aggregate to relativize Biomass & Evergreen so that we can get the sd across 
#   space & time that isn't confounded by baseline differences in model biomass, etc
# --------
# First need to get the composition variables in there
dat.ecosys <- merge(dat.ecosys, ecosys[ecosys$Resolution=="t.001",c("Model", "Site", "Year", "Evergreen", "Deciduous", "Grass")], all.x=T, all.y=F)
dat.ecosys2 <- merge(dat.ecosys2, ecosys[ecosys$Resolution=="t.001",c("Model", "Site", "Year", "Evergreen", "Deciduous", "Grass")], all.x=T, all.y=F)
summary(dat.ecosys)
summary(dat.ecosys2)

for(m in unique(dat.ecosys2[dat.ecosys2$data.type=="Model", "Model"])){
  if(m == "TreeRingNPP") ext.full <- "1980-2010" else if(m=="TreeRingRW") ext.full <- "1901-2010" else ext.full <- "850-2010"
#   biomass.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Biomass"], na.rm=T)
#   dat.ecosys2[dat.ecosys2$Model==m, "Biomass.rel"] <- dat.ecosys2[dat.ecosys2$Model==m, "Biomass"]/biomass.mean
  for(p in unique(dat.ecosys$PlotID)){
    biomass.mean <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass"], na.rm=T)
    dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass.rel"] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass"]/biomass.mean  
  }
  
  if(is.na(mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Evergreen"], na.rm=T))) next # Skip things that we don't have composition shift handy for
#   evergreen.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Evergreen"], na.rm=T)
#   dat.ecosys2[dat.ecosys2$Model==m, "Evergreen.rel"] <- dat.ecosys2[dat.ecosys2$Model==m, "Evergreen"]/evergreen.mean  
  for(p in unique(dat.ecosys$PlotID)){
    evergreen.mean <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen"], na.rm=T)
    dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen.rel"] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen"]/evergreen.mean  
  } 
}
summary(dat.ecosys2)

for(m in unique(ci.terms[ci.terms$data.type=="Model","Model"])){
  for(e in unique(ci.terms$Extent)){
    yr.min <- as.numeric(strsplit(e, "-")[[1]][1])
    
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Y.sd"]             <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Y"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Y.rel.sd"]         <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Y.rel"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Biomass.mean"]     <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Biomass"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Biomass.rel.sd"]   <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Biomass.rel"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.mean"]   <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Deciduous.mean"]   <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Deciduous"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Grass.mean"]       <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Grass"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.sd"]     <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen"], na.rm=T)
    ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.rel.sd"] <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen.rel"], na.rm=T)
    
  }
}
}
# --------
summary(ci.terms)


# 4. Condensing the full sensitivity curves to the values at the 25, 50, and 75% for 
#    analysis of continuous characteristics of models (i.e. NPP, modern change, etc)
{
  summary(ci.terms)
  summary(dat.ecosys2)

  factors.analy <- c("Y.sd", "Y.rel.sd","Biomass.sd", "Evergreen.sd")
  factors.agg <- c("Effect", "Extent", "veg.scheme", "fire.scheme", "Model", "x", "Quantile")
  
  co2.driver     <- dat.ecosys2[dat.ecosys2$Model=="ed2","CO2"    ]
  tair.driver    <- dat.ecosys2[dat.ecosys2$Model=="ed2","tair"   ]
  precipf.driver <- dat.ecosys2[dat.ecosys2$Model=="ed2","precipf"]
  
  ci.terms[,"Quantile"] <- "Other"
  for(e in c("tair", "precipf", "CO2")){
    # Get the distirbution of met drivers for each driver & specify to what precision we want to round
    effect.driver <- dat.ecosys2[dat.ecosys2$Model=="ed2", e]
    if(e=="CO2") r=0 else if(e=="tair") r=1 else r=-1
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.05, na.rm=T),r), "Quantile"] <- "q05"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.25, na.rm=T),r), "Quantile"] <- "q25"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.50, na.rm=T),r), "Quantile"] <- "q50"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.75, na.rm=T),r), "Quantile"] <- "q75"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.95, na.rm=T),r), "Quantile"] <- "q95"
  }
  ci.terms$Quantile <- as.factor(ci.terms$Quantile)
  summary(ci.terms)
  
  factors.agg <- c("x", factors.analy, "Fire", "mean.rel", "lwr.rel", "upr.rel")
  ci.terms.agg <- aggregate(ci.terms[,factors.agg], 
                            by=ci.terms[,c("Effect", "Extent", "data.type", "Y.type", "Model", "Model.Order", "veg.scheme", "fire.scheme", "Quantile")],
                            FUN=mean, na.rm=T)
  summary(ci.terms.agg)
}

# 5. Find the model-relative shift in sensitivity by comparing the slopes (1st derivatives) relative to a base extent
#    -- Base Extent == "1980-2010
{
  pdf(file.path(fig.dir, "Sensitivity_byModel_byExtent.pdf"))
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[ci.terms$Effect==v,]) +
        facet_wrap(~Model) +
        geom_ribbon(aes(x=x, ymin=lwr.rel.cent, ymax=upr.rel.cent, fill=Extent),alpha=0.5) +
        geom_line(aes(x=x, y=mean.rel.cent, color=Extent, linetype=Extent), size=2) +
        ggtitle(v) +
        scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()
    )
  }
  dev.off()

  png(file.path(fig.dir, "SuppFig1_Sensitivity_byModel_tair.png"), height=8, width=8, units="in", res=120)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="tair",]) +
      facet_wrap(~Model.Order, scales="free_x") +
      geom_ribbon(aes(x=x-273.15, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
      geom_line(aes(x=x-273.15, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
      scale_x_continuous(name=expression(bold(paste("Temperature, May-Sep ("^"o","C)"))), expand=c(0,0), breaks=c(10, 12.5, 15, 17.5)) +
      scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
      scale_fill_manual(values=c("blue3", "red3", "green3")) +
      scale_color_manual(values=c("blue3", "red3", "green3")) +
      theme(legend.title=element_text(size=10, face="bold"),
            legend.text=element_text(size=9),
            legend.key=element_blank(),
            legend.key.size=unit(1.5, "lines"),
            legend.background=element_blank(), 
            legend.position=c(0.83, 0.2)) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_text(size=12, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0.25, "lines"))  +
      theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
    }
  dev.off()

  png(file.path(fig.dir, "SuppFig2_Sensitivity_byModel_precipf.png"), height=8, width=8, units="in", res=120)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="precipf",]) +
      facet_wrap(~Model.Order, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
      geom_line(aes(x=x, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
      scale_x_continuous(name=expression(bold(paste("Precipitation, May-Sep (mm yr"^"-1",")"))), expand=c(0,0)) +
      scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
      scale_fill_manual(values=c("blue3", "red3", "green3")) +
      scale_color_manual(values=c("blue3", "red3", "green3")) +
      theme(legend.title=element_text(size=10, face="bold"),
            legend.text=element_text(size=9),
            legend.key=element_blank(),
            legend.key.size=unit(1.5, "lines"),
            legend.background=element_blank(), 
            legend.position=c(0.92, 0.09)) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_text(size=12, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0.25, "lines"))  +
      theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
  dev.off()

  png(file.path(fig.dir, "SuppFig3_Sensitivity_byModel_CO2.png"), height=8, width=8, units="in", res=120)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="CO2",]) +
      facet_wrap(~Model.Order, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
      geom_line(aes(x=x, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
#       geom_hline(yintercept=0, linetype="dashed", size=0.5) +
      scale_x_continuous(name=expression(bold(paste("CO"["2"], " (ppm)"))), expand=c(0,0)) +
      scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
      scale_fill_manual(values=c("blue3", "red3", "green3")) +
      scale_color_manual(values=c("blue3", "red3", "green3")) +
      theme(legend.title=element_text(size=12, face="bold"),
            legend.text=element_text(size=10),
            legend.key=element_blank(),
            legend.key.size=unit(1.5, "lines"),
            legend.background=element_blank(), 
            legend.position=c(0.85, 0.1)) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_text(size=12, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0.25, "lines"))  +
      theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
  dev.off()

png(file.path(fig.dir, "SuppFig3_Sensitivity_byModel_CO2_1980-2010.png"), height=8, width=8, units="in", res=120)
{
    ggplot(data=ci.terms[ci.terms$Effect=="CO2" & ci.terms$Extent=="1980-2010",]) +
      facet_wrap(~Model.Order, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
      geom_line(aes(x=x, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
      geom_hline(yintercept=0, linetype="dashed", size=0.5) +
      scale_x_continuous(name=expression(bold(paste("CO"["2"], " (ppm)"))), expand=c(0,0)) +
      scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
      scale_fill_manual(values=c("blue3", "red3", "green3")) +
      scale_color_manual(values=c("blue3", "red3", "green3")) +
      theme(legend.title=element_text(size=12, face="bold"),
            legend.text=element_text(size=10),
            legend.key=element_blank(),
            legend.key.size=unit(1.5, "lines"),
            legend.background=element_blank(), 
            legend.position=c(0.85, 0.1)) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_text(size=12, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0.25, "lines"))  +
      theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
dev.off()

  # Trying to relativize the sensitivity to look at the change from the base extent
  # Using absolute value so that positive indicates more extreme even though this 
  #   leads to some very weird shapes
  ext.base <- "1980-2010"
  for(m in unique(ci.terms[!is.na(ci.terms$mean.rel.cent),"Model"])){
    for(v in c("tair", "precipf", "CO2")){
      mod.ext <- unique(ci.terms[ci.terms$Model==m & !is.na(ci.terms$mean.rel.cent),"Extent"]) 
      for(e in mod.ext){
        dif.x   <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "x"            ], lag=1))
        dif.y   <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"], lag=1))
        # dif.lwr <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "lwr.rel.cent" ], lag=1))
        # dif.upr <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "upr.rel.cent" ], lag=1))
        
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"]   <- dif.y
        # ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "lwr.cent.deriv" ] <- dif.lwr/dif.x
        # ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "upr.cent.deriv" ] <- dif.upr/dif.x
      }
      
      # Looking at the change in slope relative to a baseline -- looking at % change r
      if(!length(mod.ext)>1) next # Skip models/datasets where we only have one extent
      deriv.base <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==ext.base, "mean.cent.deriv"]
      cent.base  <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==ext.base, "mean.rel.cent"]
      for(e in mod.ext[which(!mod.ext==ext.base)]){
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev"]     <-     ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"] -      deriv.base
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev.abs"] <- abs(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"]) - abs(deriv.base)
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev.per"] <- (ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"] - deriv.base)/deriv.base
        
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev"]     <-     ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"] -      cent.base
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev.per"] <- (ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"] - cent.base)/cent.base
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev.abs"] <-  abs(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"]) - abs(cent.base)
        
      }
    }
  }
  
  summary(ci.terms)
  summary(ci.terms[ci.terms$Effect=="precipf",])
  
  pdf(file.path(fig.dir, "Sensitivity_ExtentComparison_Slope.pdf"))
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[ci.terms$Effect==v,]) +
        facet_wrap(~Model)+
        # geom_ribbon(aes(x=x, ymin=lwr.cent.deriv, ymax=upr.cent.deriv, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean.cent.deriv, color=Extent), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        ggtitle(paste0("First Derivative: ",v)) +
        # scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()
    )  
  }
  dev.off()
  
  pdf(file.path(fig.dir, "Sensitivity_ExtentComparison_SlopeChange.pdf"))
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[!ci.terms$Extent==ext.base & ci.terms$Effect==v,]) +
        facet_wrap(~Extent)+
        # geom_ribbon(aes(x=x, ymin=lwr.cent.deriv, ymax=upr.cent.deriv, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean.deriv.dev, color=Model), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        ggtitle(paste0("deviation from 1980: ", v)) +
        # scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=colors.use) +
        theme_bw()
    )  
  }
  dev.off()
  
  pdf(file.path(fig.dir, "Sensitivity_ExtentComparison_SlopeChangeAbs.pdf"))
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[ci.terms$Effect==v,]) +
        facet_wrap(~Model)+
        # geom_ribbon(aes(x=x, ymin=lwr.cent.deriv, ymax=upr.cent.deriv, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean.deriv.dev.abs, color=Extent), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        ggtitle(paste0("deviation from abs(1980): ",v)) +
        # scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()
    )  
  }
  dev.off()
  
  pdf(file.path(fig.dir, "Sensitivity_ExtentComparison_RelSens.pdf"))
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[ci.terms$Effect==v,]) +
        facet_wrap(~Model)+
        # geom_ribbon(aes(x=x, ymin=lwr.cent.deriv, ymax=upr.cent.deriv, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean.cent.dev, color=Extent), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        ggtitle(paste0("deviation from 1980: ",v)) +
        # scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()
    )  
  }
  for(v in c("tair", "precipf", "CO2")){
    print(
      ggplot(data=ci.terms[ci.terms$Effect==v,]) +
        facet_wrap(~Model)+
        # geom_ribbon(aes(x=x, ymin=lwr.cent.deriv, ymax=upr.cent.deriv, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean.cent.dev.abs, color=Extent), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        ggtitle(paste0("absolute deviation from abs(1980): ",v)) +
        # scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()
    )  
  }
  dev.off()
  
  
  # Getting a few quick summary stats
  summary(ci.terms)
  mean(ci.terms$mean.cent.dev.abs, na.rm=T); sd(ci.terms$mean.cent.dev.abs, na.rm=T)
  mean(ci.terms[ci.terms$Effect=="tair", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Effect=="tair", "mean.cent.dev.abs"], na.rm=T)
  mean(ci.terms[ci.terms$Effect=="precipf", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Effect=="precipf", "mean.cent.dev.abs"], na.rm=T)
  mean(ci.terms[ci.terms$Effect=="CO2", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Effect=="CO2", "mean.cent.dev.abs"], na.rm=T)
  
  # Difference in climate effect on Ring Width
  mean(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair", "mean.cent.dev.abs"], na.rm=T)
  mean(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf", "mean.cent.dev.abs"], na.rm=T)
  mean(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2", "mean.cent.dev.abs"], na.rm=T); sd(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2", "mean.cent.dev.abs"], na.rm=T)
  
  tair.rw <- lm(mean.cent.dev.abs ~ 1 ,data=ci.terms[!is.na(ci.terms$mean.cent.dev.abs) & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair",])
  precipf.rw <- lm(mean.cent.dev.abs ~ 1 ,data=ci.terms[!is.na(ci.terms$mean.cent.dev.abs) & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf",])
  CO2.rw <- lm(mean.cent.dev.abs ~ 1 ,data=ci.terms[!is.na(ci.terms$mean.cent.dev.abs) & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2",])
  summary(tair.rw)
  summary(precipf.rw)
  summary(CO2.rw)
  
  hist(ci.terms$mean.cent.dev.abs)
}

# 6. Making the data frame so we can graph & compare the curves quantitatively
{
  factors.analy <- c("Y.sd", "Y.rel.sd","Biomass.sd", "Evergreen.sd")
  factors.agg <- c("Effect", "Extent", "veg.scheme", "fire.scheme", "Model", "x", "Quantile")
  df.co2 <- aggregate(ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",factors.analy], 
                      by=ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",factors.agg],
                      FUN=mean, na.rm=T)
  summary(df.co2)
  
  
  df.tair <- aggregate(ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",factors.analy], 
                       by=ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",factors.agg],
                       FUN=mean, na.rm=T)
  summary(df.tair)
  
  df.precipf <- aggregate(ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",factors.analy], 
                          by=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",factors.agg],
                          FUN=mean, na.rm=T)
  summary(df.precipf)
}

# 7. setting up some null models for the full sensitivity curves
{
  co2.null     <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.null    <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.null <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  co2.null2     <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="CO2"    ,])
  tair.null2    <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="tair"   ,])
  precipf.null2 <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="precipf",])
  
  summary(co2.null)
  summary(tair.null)
  summary(precipf.null)
  
  # Plot the different sensitivity curves by characteristic
  co2.null.post <- post.distns(model.gam=co2.null, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.null.post$null.scheme <- df.co2$null.scheme
  co2.null.post$Effect <- as.factor("CO2")
  summary(co2.null.post)
  
  tair.null.post <- post.distns(model.gam=tair.null, model.name="Tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.null.post$null.scheme <- df.tair$null.scheme
  tair.null.post$Effect <- as.factor("tair")
  summary(tair.null.post)
  
  precipf.null.post <- post.distns(model.gam=precipf.null, model.name="Precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.null.post$null.scheme <- df.precipf$null.scheme
  precipf.null.post$Effect <- as.factor("precipf")
  summary(precipf.null.post)
  
  null.post <- rbind(tair.null.post, precipf.null.post, co2.null.post)
  summary(null.post)
}

# General differences in sensitivty across Models & effects & extents
mean(ci.terms[!ci.terms$Effect=="Biomass" & !ci.terms$Extent=="1980-2010" & ci.terms$data.type=="Model","mean.cent.dev.abs"])
sd(ci.terms[!ci.terms$Effect=="Biomass" & !ci.terms$Extent=="1980-2010" & ci.terms$data.type=="Model","mean.cent.dev.abs"])

mean(ci.terms[!ci.terms$Effect=="Biomass" & !ci.terms$Extent=="1980-2010" & ci.terms$data.type=="Model","mean.cent.dev"])
mean(ci.terms[!ci.terms$Effect=="Biomass" & !ci.terms$Extent=="1980-2010" & ci.terms$data.type=="Model","mean.rel.cent"])

}
# --------


# --------
# 5.b.1 Ensemble-level Cross-Scale shifts in sensitivity
#   Hypothesis: Model sensitivities are more similar at short temporal scales
#               because it's the differences in long-term feedbacks that cuase
#               them to diverge
# --------
{
  # Calculate ensemeble mean sensitivitiy & deviation of individual models from 
  #   the ensemble at each temporal scale
{  
  co2.ext     <- gam(mean.rel.cent ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.ext    <- gam(mean.rel.cent ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.ext <- gam(mean.rel.cent ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])

  # Plot the different sensitivity curves by characteristic
  co2.ext.post <- post.distns(model.gam=co2.ext, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.ext.post[,factors.agg] <- df.co2[,factors.agg]
  summary(co2.ext.post)
  
  tair.ext.post <- post.distns(model.gam=tair.ext, model.name="tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.ext.post[,factors.agg] <- df.tair[,factors.agg]
  summary(tair.ext.post)

  precipf.ext.post <- post.distns(model.gam=precipf.ext, model.name="precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.ext.post[,factors.agg] <- df.precipf[,factors.agg]  
  summary(precipf.ext.post)
  
  ext.post <- rbind(tair.ext.post, precipf.ext.post, co2.ext.post)
  summary(ext.post)
  
  ext.post2 <- aggregate(ext.post[,c("mean", "lwr", "upr")], by=ext.post[,c("Effect", "Extent", "x", "Quantile")], FUN=mean)
  summary(ext.post2)
  
  # Plotting everything out
  pdf(file.path(fig.dir, "Sensitivity_Ensemble_Extent.pdf"))
  {
  print(
  ggplot(data=ext.post2[,]) +
    facet_wrap(~Effect, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
    geom_line(aes(x=x, y=mean, color=Extent, linetype=Extent), size=2) +
    ggtitle("Models Ensemble by Temporal Extent") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()  
  )
  print(
  ggplot(data=ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel.cent, ymax=upr.rel.cent, fill=Extent), alpha=0.5) +
    geom_line(aes(x=x, y=mean.rel.cent, color=Extent, linetype=Extent), size=2) +
    ggtitle("Tree Rings by Temporal Extent") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()  
  )
  }  
  dev.off()
  
  
  # quantitatively comparing the sensitivities
  # Models
  {
    for(e in unique(ext.post$Extent)){
      ext.post[ext.post$Extent==e,"mean.dev"] <- ext.post[ext.post$Extent==e,"mean"]-ext.post[ext.post$Extent=="1980-2010","mean"]
    }
    summary(ext.post)
    
    # tair 
    tair.mod <- lm(abs(mean.dev) ~ Extent-1, data=ext.post[ext.post$Effect=="tair" & !ext.post$Quantile=="Other",] )
    summary(tair.mod)

    tair.mod2 <- lme(abs(mean.dev) ~ Extent-1, random=list(Quantile=~1), data=ext.post[ext.post$Effect=="tair" & !ext.post$Quantile=="Other",] )
    summary(tair.mod2)
    
    mean(abs(ext.post[ext.post$Effect=="tair" & ext.post$Extent=="1901-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="tair" & ext.post$Extent=="1901-2010","mean.dev"]))
    mean(abs(ext.post[ext.post$Effect=="tair" & ext.post$Extent=="850-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="tair" & ext.post$Extent=="850-2010","mean.dev"]))
    
    # precipf
    precipf.mod <- lm(abs(mean.dev) ~ Extent-1, data=ext.post[ext.post$Effect=="precipf" & !ext.post$Quantile=="Other",] )
    summary(precipf.mod)
  
    mean(abs(ext.post[ext.post$Effect=="precipf" & ext.post$Extent=="1901-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="precipf" & ext.post$Extent=="1901-2010","mean.dev"]))
    mean(abs(ext.post[ext.post$Effect=="precipf" & ext.post$Extent=="850-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="precipf" & ext.post$Extent=="850-2010","mean.dev"]))

    # CO2
    CO2.mod <- lm(abs(mean.dev) ~ Extent-1, data=ext.post[ext.post$Effect=="CO2" & !ext.post$Quantile=="Other",] )
    summary(CO2.mod)
    
    mean(abs(ext.post[ext.post$Effect=="CO2" & ext.post$Extent=="1901-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="CO2" & ext.post$Extent=="1901-2010","mean.dev"]))
    mean(abs(ext.post[ext.post$Effect=="CO2" & ext.post$Extent=="850-2010","mean.dev"])); sd(abs(ext.post[ext.post$Effect=="CO2" & ext.post$Extent=="850-2010","mean.dev"]))
  }
  
  # Tree Rings
  {
    # tair
    tair.rw <- lm(abs(mean.cent.dev) ~ 1, data=ci.terms[ci.terms$Extent=="1901-2010" & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair" & !ci.terms$Quantile=="Other",])
    summary(tair.rw)
    mean(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair" & ci.terms$Extent=="1901-2010","mean.cent.dev"])); sd(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="tair" & ci.terms$Extent=="1901-2010","mean.cent.dev"]))

    # precipf
    precipf.rw <- lm(abs(mean.cent.dev) ~ 1, data=ci.terms[ci.terms$Extent=="1901-2010" & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf" & !ci.terms$Quantile=="Other",])
    summary(precipf.rw)
    mean(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf" & ci.terms$Extent=="1901-2010","mean.cent.dev"])); sd(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="precipf" & ci.terms$Extent=="1901-2010","mean.cent.dev"]))
    
    # CO2
    CO2.rw <- lm(abs(mean.cent.dev) ~ 1, data=ci.terms[ci.terms$Extent=="1901-2010" & ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2" & !ci.terms$Quantile=="Other",])
    summary(CO2.rw)
    mean(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2" & ci.terms$Extent=="1901-2010","mean.cent.dev"])); sd(abs(ci.terms[ci.terms$Model=="TreeRingRW" & ci.terms$Effect=="CO2" & ci.terms$Extent=="1901-2010","mean.cent.dev"]))
    
  }
}


# Calculate model deviation from the ensemble   
{
  summary(ci.terms)
  for(m in unique(ci.terms[ci.terms$data.type=="Model", "Model"])){
    for(v in c("tair", "precipf", "CO2")){
      for(e in unique(ci.terms[ci.terms$Model==m, "Extent"])){
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"dev.ensemble"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"mean.rel.cent"] - ext.post[ext.post$Extent==e & ext.post$Effect==v & ext.post$Model==m , "mean"]
      }
    }
  }
 
# Re-graphing everything now with the relative 
{
  # Setting up the graphing data frame
  {
  models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
  colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))
  
  # Creating a cheat data frame that lets values go off the graph
  ci.terms.graph <- ci.terms
  ci.terms.graph[ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
  ci.terms.graph[!is.na(ci.terms.graph$mean.rel.cent) & ci.terms.graph$mean.rel.cent< -0.65,"mean.rel.cent"] <- NA 
  ci.terms.graph[ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
  ci.terms.graph[ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
  ci.terms.graph[which(ci.terms.graph$mean.rel>1.00),"mean.rel"] <- NA 
  ci.terms.graph[which(ci.terms.graph$mean.rel.cent>1.00),"mean.rel.cent"] <- NA 
  ci.terms.graph[ci.terms.graph$lwr.rel>(1.00),"lwr.rel"] <- 1.00
  ci.terms.graph[ci.terms.graph$upr.rel>(1.00),"upr.rel"] <- 1.00 
  ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
  
  ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
  summary(ci.terms.graph)
  
  # Grouping the kind and source of the data
  ci.terms.graph$Y.type <- as.factor(ifelse(ci.terms.graph$Model=="TreeRingRW", "RW", "NPP"))
  ci.terms.graph$data.type <- as.factor(ifelse(substr(ci.terms.graph$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
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
  
  # Playing with the extent labels a bit so that "850-2010" is "all data" and "1901-2010" is left alone
  ci.terms.graph$Extent2 <- as.factor(ifelse(ci.terms.graph$Extent=="850-2010" | 
                                               (ci.terms.graph$Extent=="1980-2010" & ci.terms.graph$Model=="TreeRingNPP") |
                                               (ci.terms.graph$Extent=="1901-2010" & ci.terms.graph$Model=="TreeRingRW")
                                             , "All Data", paste(ci.terms.graph$Extent)))
  
  ci.terms.graph[ci.terms.graph$Extent2=="All Data", c("x.min", "x.max", "mask.min", "mask.max", "line.min", "line.max")] <- NA
  summary(ci.terms.graph)
  
  
  tree.rings.1901 <- ci.terms.graph[ci.terms.graph$Model=="TreeRingRW" & ci.terms.graph$Extent=="1901-2010",]
  tree.rings.1901$Extent2 <- as.factor("1901-2010")
  summary(tree.rings.1901)
  
  tree.rings.npp <- ci.terms.graph[ci.terms.graph$Model=="TreeRingNPP" & ci.terms.graph$Extent=="1980-2010",]
  tree.rings.npp$Extent2 <- as.factor("1980-2010")
  summary(tree.rings.npp)
  
  
  ci.terms.graph <- rbind(ci.terms.graph, tree.rings.1901, tree.rings.npp)
  ci.terms.graph$Extent3 <- recode(ci.terms.graph$Extent2, "'All Data'='0'; '1901-2010'='1'; '1980-2010'='2'")
  levels(ci.terms.graph$Extent3) <- c("All Data", "1901-2010", "1980-2010")
  summary(ci.terms.graph)
  }
  ext.post$Extent3 <- as.factor(ifelse(ext.post$Extent=="850-2010", "All Data", paste(ext.post$Extent)))
  ext.post[ext.post$Effect=="tair", "x"] <- ext.post[ext.post$Effect=="tair", "x"]-273.15
  
  
  pdf(file.path(fig.dir, "Sensitivity_Ensemble_Deviation_Models.pdf"), height=8.5, width=11)
  {
    models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
    colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))
    print(
      ggplot(data=ci.terms.graph[ci.terms.graph$data.type=="Model" & ci.terms.graph$Effect %in% c("tair", "precipf", "CO2") ,]) + 
        facet_grid(Extent3~Effect, scales="free_x") +
        geom_line(aes(x=x, y=mean.rel.cent*100, color=Model.Order, linetype=Model.Order), size=1) +
        geom_ribbon(data=ext.post[ext.post$Model=="ed2",], aes(x=x, ymin=lwr*100, ymax=upr*100), size=3, alpha=0.7, fill="black") +
        geom_line(data=ext.post[ext.post$Model=="ed2",], aes(x=x, y=mean*100), size=3, alpha=0.7, color="black") +
        # Lower Shaded Region
        geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
        geom_vline(aes(xintercept=line.min), linetype="dashed") +
        # Upper Shaded Region
        geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
        geom_vline(aes(xintercept=line.max), linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="") +
        scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=c(-65, 100)) +
        scale_fill_manual(values=colors.use) +
        scale_color_manual(values=colors.use) +
        scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
        theme_bw()
    )
    print(
      ggplot(data=ci.terms.graph[ci.terms.graph$data.type=="Model" & ci.terms.graph$Effect %in% c("tair", "precipf", "CO2") ,]) + 
        facet_grid(Extent3~Effect, scales="free_x") +
        geom_line(aes(x=x, y=dev.ensemble*100, color=Model.Order, linetype=Model.Order), size=1) +
        geom_hline(yintercept=0, linetype="dashed") +
        # Lower Shaded Region
        geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
        geom_vline(aes(xintercept=line.min), linetype="dashed") +
        # Upper Shaded Region
        geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
        geom_vline(aes(xintercept=line.max), linetype="dashed") +
        scale_x_continuous(expand=c(0,0), name="") +
        scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=c(-65, 100)) +
        scale_fill_manual(values=colors.use) +
        scale_color_manual(values=colors.use) +
        scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
        theme_bw()
    )
  }
  dev.off()
}

}

# Statistical test of model agreement at different scales (Using the quantile comparisons)
{
  pdf(file.path(fig.dir, "Ensemble_Sensitivity_Model_Agreement.pdf"))
  {
  print(
  ggplot(data=ci.terms[!is.na(ci.terms$Quantile) & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + 
    facet_grid(~Effect) +
    geom_boxplot(aes(x=Extent, y=abs(dev.ensemble), fill=Extent)) +
    ggtitle("Full Met Range, quantiles") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()
  )
  print(
  ggplot(data=ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) + 
    facet_wrap(~Effect) +
    geom_violin(aes(x=Extent, y=abs(dev.ensemble), fill=Extent), adjust=2, scale="width") +
    ggtitle("1980-2010 Met Range, all values") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()
  )
  }
  dev.off()
  
  # Comparing the variation at the quantiles
  ensemble.co2     <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[!ci.terms$Quantile=="Other" & ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  ensemble.tair    <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  ensemble.precipf <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  summary(ensemble.co2)  
  summary(ensemble.tair)  
  summary(ensemble.precipf)  
  
  # Restricting our comparison to just the range of met observed in 1980-2010
  ensemble.co2.2     <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  ensemble.tair.2    <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  ensemble.precipf.2 <- lm(abs(dev.ensemble) ~ Extent, data=ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  summary(ensemble.co2.2)  
  summary(ensemble.tair.2)  
  summary(ensemble.precipf.2)  
  
}

# Quantifying the percent increase in met range relative to 1980-2010
# tair
tair.base.max <- max(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="tair", "x"], na.rm=T)
tair.base.min <- min(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="tair", "x"], na.rm=T)
tair.1901 <- max(max(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="tair", "x"], na.rm=T), tair.base.max) - min(min(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="tair", "x"], na.rm=T), tair.base.min)
tair.850  <- max(max(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="tair", "x"], na.rm=T), tair.base.max) - min(min(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="tair", "x"], na.rm=T), tair.base.min)

tair.1901/(tair.base.max - tair.base.min)
tair.850/(tair.base.max - tair.base.min)

# precipf
precipf.base.max <- max(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="precipf", "x"], na.rm=T)
precipf.base.min <- min(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="precipf", "x"], na.rm=T)
precipf.1901 <- max(max(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="precipf", "x"], na.rm=T), precipf.base.max) - min(min(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="precipf", "x"], na.rm=T), precipf.base.min)
precipf.850  <- max(max(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="precipf", "x"], na.rm=T), precipf.base.max) - min(min(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="precipf", "x"], na.rm=T), precipf.base.min)

precipf.1901/(precipf.base.max - precipf.base.min)
precipf.850/(precipf.base.max - precipf.base.min)

# CO2
CO2.base.max <- max(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="CO2", "x"], na.rm=T)
CO2.base.min <- min(ci.terms[ci.terms$Met=="1980-2010" & ci.terms$Effect=="CO2", "x"], na.rm=T)
CO2.1901 <- max(max(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="CO2", "x"], na.rm=T), CO2.base.max) - min(min(ci.terms[ci.terms$Met=="1901-2010" & ci.terms$Effect=="CO2", "x"], na.rm=T), CO2.base.min)
CO2.850  <- max(max(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="CO2", "x"], na.rm=T), CO2.base.max) - min(min(ci.terms[ci.terms$Met=="850-2010"  & ci.terms$Effect=="CO2", "x"], na.rm=T), CO2.base.min)

CO2.1901/(CO2.base.max - CO2.base.min)
CO2.850/(CO2.base.max - CO2.base.min)
}
# --------

# --------
# 5.b.2 Finding patterns in Model Shifts: dynamic or static veg
#  Hypothesis: Climate sensitivities of models with static veg are more consistent
#              across temporal scales because they less capacity for ecosystem-level
#              adaptation to climate variaibility & directional change
#  NOTE:  Here we're looking specifically at the within-model shifts rather than finding 
#         ensemble-level patterns
# --------
{

  # Making a combo site & Extent factor
  ci.terms$veg.extent <- as.factor(paste(ci.terms$veg.scheme, ci.terms$Extent, sep="."))
  summary(ci.terms)
  
  df.co2    $veg.extent <- as.factor(paste(df.co2    $veg.scheme, df.co2$Extent, sep="."))
  df.tair   $veg.extent <- as.factor(paste(df.tair   $veg.scheme, df.co2$Extent, sep="."))
  df.precipf$veg.extent <- as.factor(paste(df.precipf$veg.scheme, df.co2$Extent, sep="."))
# ------------  
# Veg scheme & Sensitivity across scales
# ------------
{

  co2.veg.lm1     <- lm(abs(mean.rel.cent) ~ Extent*veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.veg.lm1    <- lm(abs(mean.rel.cent) ~ veg.scheme*Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.veg.lm1 <- lm(abs(mean.rel.cent) ~ Extent*veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  anova(co2.veg.lm1)
  anova(tair.veg.lm1)
  anova(precipf.veg.lm1) # Precipf = no veg  or veg x extent interaction 
  
  tair.veg.t.veg    <- t.test(abs(mean.rel.cent) ~ veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  tair.veg.lm.ext    <- lm(mean.cent.dev.abs ~ Extent-1, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  tair.veg.t.veg
  summary(tair.veg.lm.ext)
  
  co2.veg.lm2     <- lm(abs(mean.rel.cent) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.veg.lm2    <- lm(abs(mean.rel.cent) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.veg.lm2 <- lm(abs(mean.rel.cent) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.veg.lm2)
  summary(tair.veg.lm2)
  summary(precipf.veg.lm2) 
  
  veg.precip <- data.frame(char="veg.scheme (static)", 
                           effect="precip", 
                           extent=c("1980-2010","1901-2010", "850-2010"),
                           estimate=summary(precipf.veg.lm2)$coefficients[4:6,1],
                           std.err = summary(precipf.veg.lm2)$coefficients[4:6,2],
                           t.stat  = summary(precipf.veg.lm2)$coefficients[4:6,3],
                           p.value = summary(precipf.veg.lm2)$coefficients[4:6,4])
  veg.tair <- data.frame(char="veg.scheme (static)", 
                           effect="tair", 
                           extent=c("1980-2010","1901-2010", "850-2010"),
                           estimate=summary(tair.veg.lm2)$coefficients[4:6,1],
                           std.err = summary(tair.veg.lm2)$coefficients[4:6,2],
                           t.stat  = summary(tair.veg.lm2)$coefficients[4:6,3],
                           p.value = summary(tair.veg.lm2)$coefficients[4:6,4])
  veg.co2 <- data.frame(char="veg.scheme (static)", 
                         effect="co2", 
                         extent=c("1980-2010","1901-2010", "850-2010"),
                         estimate=summary(co2.veg.lm2)$coefficients[4:6,1],
                         std.err = summary(co2.veg.lm2)$coefficients[4:6,2],
                         t.stat  = summary(co2.veg.lm2)$coefficients[4:6,3],
                         p.value = summary(co2.veg.lm2)$coefficients[4:6,4])
  
  veg.stats <- rbind(veg.tair, veg.precip, veg.co2)
  
  
  
  co2.veg.lme     <- lme(abs(mean.rel.cent) ~ veg.scheme, random=list(Extent=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.veg.lme    <- lme(abs(mean.rel.cent) ~ veg.scheme, random=list(Extent=~1, x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.veg.lme <- lme(abs(mean.rel.cent) ~ veg.scheme, random=list(Extent=~1, x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.veg.lme)
  summary(tair.veg.lme)
  summary(precipf.veg.lme)
  
  
  
  
  co2.veg.ext     <- gam(mean.rel.cent ~ s(x, by=veg.extent), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.veg.ext    <- gam(mean.rel.cent ~ s(x, by=veg.extent), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.veg.ext <- gam(mean.rel.cent ~ s(x, by=veg.extent), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  # Plot the different sensitivity curves by characteristic
  co2.veg.ext.post <- post.distns(model.gam=co2.veg.ext, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.veg.ext.post[,factors.agg] <- df.co2[,factors.agg]
  summary(co2.veg.ext.post)
  
  tair.veg.ext.post <- post.distns(model.gam=tair.veg.ext, model.name="tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.veg.ext.post[,factors.agg] <- df.tair[,factors.agg]
  summary(tair.veg.ext.post)
  
  precipf.veg.ext.post <- post.distns(model.gam=precipf.veg.ext, model.name="precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.veg.ext.post[,factors.agg] <- df.precipf[,factors.agg]  
  summary(precipf.veg.ext.post)
  
  veg.ext.post <- rbind(tair.veg.ext.post, precipf.veg.ext.post, co2.veg.ext.post)
  summary(veg.ext.post)
  
  veg.ext.post2 <- aggregate(veg.ext.post[,c("mean", "lwr", "upr")], by=veg.ext.post[,c("veg.scheme", "Effect", "Extent", "x", "Quantile")], FUN=mean)
  summary(veg.ext.post2)
  
  # Plotting everything out
  pdf(file.path(fig.dir, "Sensitivity_VegScheme_Extent.pdf"))
  {
    print(
      ggplot(data=veg.ext.post2[,]) +
        facet_grid(Extent~Effect, scales="free_x") +
        geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=veg.scheme), alpha=0.5) +
        geom_line(aes(x=x, y=mean, color=veg.scheme, linetype=veg.scheme), size=2) +
        ggtitle("Models Ensemble by Temporal Extent") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()  
    )
    print(
      ggplot(data=veg.ext.post2[,]) +
        facet_grid(veg.scheme~Effect, scales="free_x") +
        geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean, color=Extent, linetype=Extent), size=2) +
        ggtitle("Models Ensemble by Temporal Extent") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()  
    )
  }
  dev.off()
}
# ------------    
  
  
    
# ------------
# Veg scheme & scalability (change in sensitivity) across scales
# ------------
{    
  pdf(file.path(fig.dir, "Sensitivity_Change_VegScheme.pdf"))
  {
  print(
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other" & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_grid(Effect~Extent, scales="free") +
    geom_violin(aes(x=veg.scheme, y=mean.cent.deriv, fill=veg.scheme), adjust=2) +
    scale_fill_manual(values=c("blue3", "red3")) +
    theme_bw()
  )
  print(
  ggplot(data=ci.terms[ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_grid(Effect~Extent, scales="free") +
    geom_violin(aes(x=veg.scheme, y=mean.cent.deriv, fill=veg.scheme), adjust=2) +
    scale_fill_manual(values=c("blue3", "red3")) +
    theme_bw()
  )
  print(
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_grid(Effect~Extent, scales="free") +
    geom_violin(aes(x=veg.scheme, y=mean.deriv.dev, fill=veg.scheme), adjust=3) +
    scale_fill_manual(values=c("blue3", "red3")) +
    theme_bw()
  )
  print(
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_grid(Effect~Extent, scales="free") +
    geom_violin(aes(x=veg.scheme, y=mean.deriv.dev.abs, fill=veg.scheme), adjust=3) +
    scale_fill_manual(values=c("blue3", "red3")) +
    theme_bw()
  )
  }
  dev.off()
  

  co2.veg.lm     <- lm(mean.cent.dev.abs ~ veg.scheme*Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$Effect=="CO2",])
  tair.veg.lm    <- lm(mean.cent.dev.abs ~ veg.scheme*Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$Effect=="tair",])
  precipf.veg.lm <- lm(mean.cent.dev.abs ~ veg.scheme*Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$Effect=="precipf",])
  anova(co2.veg.lm)
  anova(tair.veg.lm)
  anova(precipf.veg.lm) # Precipf = no veg  or veg x extent interaction 

  precipf.veg.aov <- aov(mean.cent.dev.abs ~ veg.scheme*Extent, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$Effect=="precipf",])
  
  TukeyHSD(precipf.veg.aov)
  
  co2.veg.lme     <- lme(mean.cent.dev.abs ~ veg.scheme-1, random=list(Extent=~1,x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.veg.lme    <- lme(mean.cent.dev.abs ~ veg.scheme-1, random=list(Extent=~1,x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.veg.lme <- lme(mean.cent.dev.abs ~ veg.scheme-1, random=list(Extent=~1,x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.veg.lme)
  summary(tair.veg.lme)
  summary(precipf.veg.lme)
  
  
  precipf.veg.lme2 <- lme(mean.cent.dev.abs ~ Extent*(veg.scheme-1), random=list(x=~1), data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(precipf.veg.lme2)
  anova(precipf.veg.lme2)
  
  anova(co2.veg.lme)
  anova(tair.veg.lme)
  anova(precipf.veg.lme)
  
  # Looking at veg x extent interactions
  {
    co2.veg.gam     <- gam(mean.cent.dev.abs ~ s(x,by=veg.extent)-1, data=ci.terms[!ci.terms$Extent==ext.base & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
    tair.veg.gam    <- gam(mean.cent.dev.abs ~ s(x,by=veg.extent)-1, data=ci.terms[!ci.terms$Extent==ext.base & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])  
    precipf.veg.gam <- gam(mean.cent.dev.abs ~ s(x,by=veg.extent)-1, data=ci.terms[!ci.terms$Extent==ext.base & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
    summary(co2.veg.gam)
    summary(tair.veg.gam)
    summary(precipf.veg.gam)
  
    # doing the posterior predictions
    df.co2    $veg.extent <- as.factor(paste(df.co2    $veg.scheme, df.co2    $Extent, sep="."))
    df.tair   $veg.extent <- as.factor(paste(df.tair   $veg.scheme, df.tair   $Extent, sep="."))
    df.precipf$veg.extent <- as.factor(paste(df.precipf$veg.scheme, df.precipf$Extent, sep="."))

    co2.veg.post <- post.distns(model.gam=co2.veg.gam, model.name="CO2", n=50, newdata=df.co2[!df.co2$Extent==ext.base, ], vars="x", terms=F)$ci
    co2.veg.post[,factors.agg] <- df.co2[!df.co2$Extent==ext.base, factors.agg]
    co2.veg.post$veg.extent <- df.co2[!df.co2$Extent==ext.base, "veg.extent"]
    co2.veg.post$Effect <- as.factor("CO2")
    summary(co2.veg.post)
  
    tair.veg.post <- post.distns(model.gam=tair.veg.gam, model.name="tair", n=50, newdata=df.tair[!df.tair$Extent==ext.base, ], vars="x", terms=F)$ci
    tair.veg.post[,factors.agg] <- df.tair[!df.tair$Extent==ext.base, factors.agg]
    tair.veg.post$veg.extent <- df.tair[!df.tair$Extent==ext.base, "veg.extent"]
    tair.veg.post$Effect <- as.factor("tair")
    summary(tair.veg.post)

    precipf.veg.post <- post.distns(model.gam=precipf.veg.gam, model.name="precipf", n=50, newdata=df.precipf[!df.precipf$Extent==ext.base, ], vars="x", terms=F)$ci
    precipf.veg.post[,factors.agg] <- df.precipf[!df.precipf$Extent==ext.base, factors.agg]
    precipf.veg.post$veg.extent <- df.precipf[!df.precipf$Extent==ext.base, "veg.extent"]
    precipf.veg.post$Effect <- as.factor("precipf")
    summary(precipf.veg.post)
  
    veg.post <- rbind(tair.veg.post, precipf.veg.post, co2.veg.post)
    veg.post <- aggregate(veg.post[,c("mean", "lwr", "upr")], 
                          by=veg.post[,c("Effect", "veg.scheme", "Extent", "veg.extent", "x", "Quantile")],
                          FUN=mean)
    summary(veg.post)
    
    pdf(file.path(fig.dir, "Sensitivity_ExtentDeviation_VegExtent.pdf"))
    {
    print(
    ggplot(data=veg.post[!veg.post$Extent==ext.base, ]) +
      facet_grid(Extent~Effect,scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=veg.scheme), alpha=0.5)+
      geom_line(aes(x=x, y=mean, color=veg.scheme, linetype=veg.scheme), size=2) +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_y_continuous(name="mean.cent.dev.abs") +
      scale_fill_manual(values=c("red3", "green3")) +
      scale_color_manual(values=c("red3", "green3")) +
      theme_bw()
      )
    }
    dev.off()
  }

  # Looking at overall behavior of veg scheme while taking into account the extent
  co2.veg.gam2     <- gam(mean.deriv.dev ~ s(x,by=veg.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.veg.gam2    <- gam(mean.deriv.dev ~ s(x,by=veg.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.veg.gam2 <- gam(mean.deriv.dev ~ s(x,by=veg.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
#  summary(co2.veg.gam1$lme)  
  summary(co2.veg.gam2); anova(co2.veg.gam2)
  summary(tair.veg.gam2); anova(tair.veg.gam2)
  summary(precipf.veg.gam2); anova(tair.veg.gam2)

  co2.veg.post2 <- post.distns(model.gam=co2.veg.gam2, model.name="CO2", n=50, newdata=df.co2[!df.co2$Extent==ext.base, ], vars="x", terms=T)$ci
  co2.veg.post2$veg.scheme <- df.co2[!df.co2$Extent==ext.base, "veg.scheme"]
  co2.veg.post2$veg.extent <- df.co2[!df.co2$Extent==ext.base, "veg.extent"]
  co2.veg.post2$Model      <- df.co2[!df.co2$Extent==ext.base, "Model"     ]  
  co2.veg.post2$Effect <- as.factor("CO2")
  summary(co2.veg.post2)

  tair.veg.post2 <- post.distns(model.gam=tair.veg.gam2, model.name="tair", n=50, newdata=df.tair[!df.tair$Extent==ext.base, ], vars="x", terms=T)$ci
  tair.veg.post2$veg.scheme <- df.tair[!df.tair$Extent==ext.base, "veg.scheme"]
  tair.veg.post2$veg.extent <- df.tair[!df.tair$Extent==ext.base, "veg.extent"]
  tair.veg.post2$Model      <- df.tair[!df.tair$Extent==ext.base, "Model"     ]  
  tair.veg.post2$Effect <- as.factor("tair")
  summary(tair.veg.post2)

  precipf.veg.post2 <- post.distns(model.gam=precipf.veg.gam2, model.name="precipf", n=50, newdata=df.precipf[!df.precipf$Extent==ext.base, ], vars="x", terms=T)$ci
  precipf.veg.post2$veg.scheme <- df.precipf[!df.precipf$Extent==ext.base, "veg.scheme"]
  precipf.veg.post2$veg.extent <- df.precipf[!df.precipf$Extent==ext.base, "veg.extent"]
  precipf.veg.post2$Model      <- df.precipf[!df.precipf$Extent==ext.base, "Model"     ]  
  precipf.veg.post2$Effect <- as.factor("precipf")
  summary(precipf.veg.post2)

  veg.post2 <- rbind(tair.veg.post2, precipf.veg.post2, co2.veg.post2)
  summary(veg.post2)

  pdf(file.path(fig.dir, "Sensitivity_ExtentDeviation_VegScheme.pdf"))
  {
    print(
      ggplot(data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & !ci.terms$Extent==ext.base,]) +
        facet_grid(Extent~Effect,scales="free_x") +
#         geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=veg.scheme), alpha=0.5)+
        geom_line(aes(x=x, y=mean.deriv.dev, color=Model), size=2) +
        geom_hline(yintercept=0, linetype="dashed") +
        scale_fill_manual(values=colors.use) +
        scale_color_manual(values=colors.use) +
        theme_bw()
    )
    print(
    ggplot(data=veg.post2[veg.post2$Model%in% c("ed2", "clm.bgc") & !veg.post2$Extent==ext.base, ]) +
      facet_wrap(~Effect,scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=veg.scheme), alpha=0.5)+
      geom_line(aes(x=x, y=mean, color=veg.scheme, linetype=veg.scheme), size=2) +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_fill_manual(values=c("blue3", "red3")) +
      scale_color_manual(values=c("blue3", "red3")) +
      theme_bw()
    )
  }
  dev.off()
}
# --------

}
# --------

# --------
# 5.b.4 Slow processes & sensitivity change: Composition (Fraction Evergreen)
#  Hypothesis: Like Biomass, shifts in composition is a slow process that can mediate
#              climate sensitivity at more long time scales but has limited influence
#              at short time scales
# --------
{
  
  # ------------
  # Correlation of Composition with sensitivity across scales
  # ------------
{  
  # -----
  # Evergreen!
  # -----
{
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_vs_Evergreen.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Evergreen.mean, y=abs(mean.rel.cent), color=Extent), size=5) +
    stat_smooth(aes(x=Evergreen.mean, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.evg.lm     <- lm(abs(mean.rel.cent) ~ Evergreen.mean*(Extent-1) - Evergreen.mean , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm    <- lm(abs(mean.rel.cent) ~ Evergreen.mean*(Extent-1) - Evergreen.mean , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm <- lm(abs(mean.rel.cent) ~ Evergreen.mean*(Extent-1) - Evergreen.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lm)
  summary(tair.evg.lm)
  summary(precipf.evg.lm)
}
# -----

# -----
# Deciduous!
# -----
{
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_vs_Evergreen.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Deciduous.mean, y=abs(mean.rel.cent), color=Extent), size=5) +
    stat_smooth(aes(x=Deciduous.mean, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.dec.lm     <- lm(abs(mean.rel.cent) ~ Deciduous.mean*(Extent-1) - Deciduous.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.dec.lm    <- lm(abs(mean.rel.cent) ~ Deciduous.mean*(Extent-1) - Deciduous.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.dec.lm <- lm(abs(mean.rel.cent) ~ Deciduous.mean*(Extent-1) - Deciduous.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.dec.lm)
  summary(tair.dec.lm)
  summary(precipf.dec.lm)
}
# -----

# -----
# Grass!
# -----
{
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_vs_Grass.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Grass.mean, y=abs(mean.rel.cent), color=Extent), size=5) +
    stat_smooth(aes(x=Grass.mean, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.grass.lm     <- lm(abs(mean.rel.cent) ~ Grass.mean*(Extent-1) - Grass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.grass.lm    <- lm(abs(mean.rel.cent) ~ Grass.mean*(Extent-1) - Grass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.grass.lm <- lm(abs(mean.rel.cent) ~ Grass.mean*(Extent-1) - Grass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.grass.lm)
  summary(tair.grass.lm)
  summary(precipf.grass.lm)
}
# -----

summary(co2.evg.lm)
summary(co2.dec.lm)
summary(co2.grass.lm)

summary(tair.evg.lm)
summary(tair.dec.lm)
summary(tair.grass.lm)

summary(precipf.evg.lm)
summary(precipf.dec.lm)
summary(precipf.grass.lm)

}
# ------------


# ------------
# Correlation of Evergreen *Variability* with sensitivity across scales
# ------------
{
  # Correlation of Evergreen with relative sensitivity across all scales
  pdf(file.path(fig.dir, "Sensitivity_Slope_vs_EvergreenVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Evergreen.sd, y=abs(mean.cent.deriv), color=Extent), size=5) +
    stat_smooth(aes(x=Evergreen.sd, y=abs(mean.cent.deriv), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()
  dev.off()
  

  co2.evg.lm0     <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent) , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm0    <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent) , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm0 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  anova(co2.evg.lm0)
  anova(tair.evg.lm0)
  anova(precipf.evg.lm0)
  
  
  co2.evg.lm1     <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm1    <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm1 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lm1)
  summary(tair.evg.lm1)
  summary(precipf.evg.lm1)
  
  

  
  summary(ci.terms[ci.terms$Model=="ed2",])
  summary(ci.terms[ci.terms$Model=="ed2",])
  
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_vs_EvergreenVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Evergreen.sd, y=abs(mean.rel.cent), color=Extent), size=5) +
    stat_smooth(aes(x=Evergreen.sd, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.evg.lm2     <- lm(abs(mean.rel.cent) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm2    <- lm(abs(mean.rel.cent) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm2 <- lm(abs(mean.rel.cent) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lm2)
  summary(tair.evg.lm2)
  summary(precipf.evg.lm2)
  
  
  evg.precip <- data.frame(char="evg.sd", 
                           effect="precip", 
                           extent=c("1980-2010","1901-2010", "850-2010"),
                           estimate=summary(precipf.evg.lm2)$coefficients[4:6,1],
                           std.err = summary(precipf.evg.lm2)$coefficients[4:6,2],
                           t.stat  = summary(precipf.evg.lm2)$coefficients[4:6,3],
                           p.value = summary(precipf.evg.lm2)$coefficients[4:6,4])
  evg.tair <- data.frame(char="evg.sd", 
                         effect="tair", 
                         extent=c("1980-2010","1901-2010", "850-2010"),
                         estimate=summary(tair.evg.lm2)$coefficients[4:6,1],
                         std.err = summary(tair.evg.lm2)$coefficients[4:6,2],
                         t.stat  = summary(tair.evg.lm2)$coefficients[4:6,3],
                         p.value = summary(tair.evg.lm2)$coefficients[4:6,4])
  evg.co2 <- data.frame(char="evg.sd", 
                        effect="co2", 
                        extent=c("1980-2010","1901-2010", "850-2010"),
                        estimate=summary(co2.evg.lm2)$coefficients[4:6,1],
                        std.err = summary(co2.evg.lm2)$coefficients[4:6,2],
                        t.stat  = summary(co2.evg.lm2)$coefficients[4:6,3],
                        p.value = summary(co2.evg.lm2)$coefficients[4:6,4])
  
  evg.stats <- rbind(evg.tair, evg.precip, evg.co2)
  
  #   co2.evg.lme     <- lme(abs(mean.rel.cent) ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  #   tair.evg.lme    <- lme(abs(mean.rel.cent) ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  #   precipf.evg.lme <- lme(abs(mean.rel.cent) ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  #   summary(co2.evg.lme)
  #   summary(tair.evg.lme)
  #   summary(precipf.evg.lme)
}
# ------------

# ------------
# Correlation of Evergreen with Scalability of short-term sensitivity
# ------------
{
  # Correlation of Evergreen with relative sensitivity across all scales
  pdf(file.path(fig.dir, "Sensitivity_Change_Slope_vs_EvergreenVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & !ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Evergreen.sd, y=mean.deriv.dev.abs, color=Extent), size=5) +
    stat_smooth(aes(x=Evergreen.sd, y=mean.deriv.dev.abs, color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("red3", "green3")) +
    scale_fill_manual(values=c("red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.evg.lm1     <- lm(abs(mean.deriv.dev.abs) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm1    <- lm(abs(mean.deriv.dev.abs) ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm1 <- lm(abs(mean.deriv.dev.abs) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lm1)
  summary(tair.evg.lm1)
  summary(precipf.evg.lm1)
  
  summary(ci.terms[ci.terms$Model=="ed2",])
  summary(ci.terms[ci.terms$Model=="jules.stat",])
  summary(ci.terms[ci.terms$Evergreen.sd>=0.55,])
  summary(ci.terms[ci.terms$Evergreen.sd<=0.35,])
  summary(ci.terms$Evergreen.sd)
  
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_Change_vs_EvergreenVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & !ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Evergreen.sd, y=mean.cent.dev.abs, color=Extent), size=5) +
    stat_smooth(aes(x=Evergreen.sd, y=mean.cent.dev.abs, color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.evg.lm2     <- lm(mean.cent.dev.abs ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lm2    <- lm(mean.cent.dev.abs ~ Evergreen.sd*(Extent-1) - Evergreen.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lm2 <- lm(mean.cent.dev.abs ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lm2)
  summary(tair.evg.lm2)
  summary(precipf.evg.lm2)
  
  
  co2.evg.lme     <- lme(mean.rel.cent ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & !is.na(ci.terms$Evergreen.sd) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.evg.lme    <- lme(mean.rel.cent ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & !is.na(ci.terms$Evergreen.sd) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.evg.lme <- lme(mean.rel.cent ~ Evergreen.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & !is.na(ci.terms$Evergreen.sd) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.evg.lme)
  summary(tair.evg.lme)
  summary(precipf.evg.lme)
}
# ------------


}
# --------

# --------
# 5.b.3 Finding patterns in Model Shifts: presence of fire
# --------
{
  
  # Making a combo site & Extent factor
  ci.terms$fire.extent <- as.factor(paste(ci.terms$fire.scheme, ci.terms$Extent, sep="."))
  summary(ci.terms)
  
  df.co2    $fire.extent <- as.factor(paste(df.co2    $fire.scheme, df.co2$Extent, sep="."))
  df.tair   $fire.extent <- as.factor(paste(df.tair   $fire.scheme, df.co2$Extent, sep="."))
  df.precipf$fire.extent <- as.factor(paste(df.precipf$fire.scheme, df.co2$Extent, sep="."))
  # ------------  
  # fire scheme & Sensitivity across scales
  # ------------
{
  co2.fire.lm1     <- lm(abs(mean.rel.cent) ~ Extent*fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.fire.lm1    <- lm(abs(mean.rel.cent) ~ Extent*fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.fire.lm1 <- lm(abs(mean.rel.cent) ~ Extent*fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  anova(co2.fire.lm1)
  anova(tair.fire.lm1)
  anova(precipf.fire.lm1) # Precipf = no veg  or veg x extent interaction 
  
  
  co2.fire.lm2     <- lm(abs(mean.rel.cent) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.fire.lm2    <- lm(abs(mean.rel.cent) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.fire.lm2 <- lm(abs(mean.rel.cent) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.fire.lm2)
  summary(tair.fire.lm2)
  summary(precipf.fire.lm2) 
  
  fire.precip <- data.frame(char="fire (yes)", 
                           effect="precip", 
                           extent=c("1980-2010","1901-2010", "850-2010"),
                           estimate=summary(precipf.fire.lm2)$coefficients[4:6,1],
                           std.err = summary(precipf.fire.lm2)$coefficients[4:6,2],
                           t.stat  = summary(precipf.fire.lm2)$coefficients[4:6,3],
                           p.value = summary(precipf.fire.lm2)$coefficients[4:6,4])
  fire.tair <- data.frame(char="fire (yes)", 
                         effect="tair", 
                         extent=c("1980-2010","1901-2010", "850-2010"),
                         estimate=summary(tair.fire.lm2)$coefficients[4:6,1],
                         std.err = summary(tair.fire.lm2)$coefficients[4:6,2],
                         t.stat  = summary(tair.fire.lm2)$coefficients[4:6,3],
                         p.value = summary(tair.fire.lm2)$coefficients[4:6,4])
  fire.co2 <- data.frame(char="fire (yes)", 
                        effect="co2", 
                        extent=c("1980-2010","1901-2010", "850-2010"),
                        estimate=summary(co2.fire.lm2)$coefficients[4:6,1],
                        std.err = summary(co2.fire.lm2)$coefficients[4:6,2],
                        t.stat  = summary(co2.fire.lm2)$coefficients[4:6,3],
                        p.value = summary(co2.fire.lm2)$coefficients[4:6,4])
  
  fire.stats <- rbind(fire.tair, fire.precip, fire.co2)
  
  
  
  co2.fire.ext     <- gam(mean.rel.cent ~ s(x, by=fire.extent), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.fire.ext    <- gam(mean.rel.cent ~ s(x, by=fire.extent), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.fire.ext <- gam(mean.rel.cent ~ s(x, by=fire.extent), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  # Plot the different sensitivity curves by characteristic
  co2.fire.ext.post <- post.distns(model.gam=co2.fire.ext, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.fire.ext.post[,factors.agg] <- df.co2[,factors.agg]
  summary(co2.fire.ext.post)
  
  tair.fire.ext.post <- post.distns(model.gam=tair.fire.ext, model.name="tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.fire.ext.post[,factors.agg] <- df.tair[,factors.agg]
  summary(tair.fire.ext.post)
  
  precipf.fire.ext.post <- post.distns(model.gam=precipf.fire.ext, model.name="precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.fire.ext.post[,factors.agg] <- df.precipf[,factors.agg]  
  summary(precipf.fire.ext.post)
  
  fire.ext.post <- rbind(tair.fire.ext.post, precipf.fire.ext.post, co2.fire.ext.post)
  summary(fire.ext.post)
  
  fire.ext.post2 <- aggregate(fire.ext.post[,c("mean", "lwr", "upr")], by=fire.ext.post[,c("fire.scheme", "Effect", "Extent", "x", "Quantile")], FUN=mean)
  summary(fire.ext.post2)
  
  # Plotting everything out
  pdf(file.path(fig.dir, "Sensitivity_FireScheme_Extent.pdf"))
  {
    print(
      ggplot(data=fire.ext.post2[,]) +
        facet_grid(Extent~Effect, scales="free_x") +
        geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=fire.scheme), alpha=0.5) +
        geom_line(aes(x=x, y=mean, color=fire.scheme, linetype=fire.scheme), size=2) +
        ggtitle("Models Ensemble by Temporal Extent") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()  
    )
    print(
      ggplot(data=fire.ext.post2[,]) +
        facet_grid(fire.scheme~Effect, scales="free_x") +
        geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
        geom_line(aes(x=x, y=mean, color=Extent, linetype=Extent), size=2) +
        ggtitle("Models Ensemble by Temporal Extent") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("blue3", "red3", "green3")) +
        scale_color_manual(values=c("blue3", "red3", "green3")) +
        theme_bw()  
    )
  }
  dev.off()
}
# ------------    

# ------------    
# Fire presence & frequencey effect on scalability of observations
# ------------    
{
co2.fire.lm     <- lm(mean.cent.dev.abs ~ fire.scheme*(Extent), data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
tair.fire.lm    <- lm(mean.cent.dev.abs ~ fire.scheme*(Extent), data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
precipf.fire.lm <- lm(mean.cent.dev.abs ~ fire.scheme*(Extent), data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model"& !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model"  & ci.terms$Effect=="precipf",])
anova(co2.fire.lm)
anova(tair.fire.lm)
anova(precipf.fire.lm) # Precipf = no veg  or veg x extent interaction 

summary(co2.fire.lm)


# Looking at overall behavior of fire scheme while taking into account the extent
co2.fire.gam2     <- gam(mean.cent.dev.abs ~ s(x,by=fire.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
tair.fire.gam2    <- gam(mean.cent.dev.abs ~ s(x,by=fire.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
precipf.fire.gam2 <- gam(mean.cent.dev.abs ~ s(x,by=fire.scheme) + Extent-1, data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
#  summary(co2.fire.gam1$lme)  
summary(co2.fire.gam2); anova(co2.fire.gam2)
summary(tair.fire.gam2); anova(tair.fire.gam2)
summary(precipf.fire.gam2); anova(tair.fire.gam2)

co2.fire.post2 <- post.distns(model.gam=co2.fire.gam2, model.name="CO2", n=50, newdata=df.co2[!df.co2$Extent==ext.base, ], vars="x", terms=T)$ci
co2.fire.post2$fire.scheme <- df.co2[!df.co2$Extent==ext.base, "fire.scheme"]
co2.fire.post2$fire.extent <- df.co2[!df.co2$Extent==ext.base, "fire.extent"]
co2.fire.post2$Model      <- df.co2[!df.co2$Extent==ext.base, "Model"     ]  
co2.fire.post2$Effect <- as.factor("CO2")
summary(co2.fire.post2)

tair.fire.post2 <- post.distns(model.gam=tair.fire.gam2, model.name="tair", n=50, newdata=df.tair[!df.tair$Extent==ext.base, ], vars="x", terms=T)$ci
tair.fire.post2$fire.scheme <- df.tair[!df.tair$Extent==ext.base, "fire.scheme"]
tair.fire.post2$fire.extent <- df.tair[!df.tair$Extent==ext.base, "fire.extent"]
tair.fire.post2$Model      <- df.tair[!df.tair$Extent==ext.base, "Model"     ]  
tair.fire.post2$Effect <- as.factor("tair")
summary(tair.fire.post2)

precipf.fire.post2 <- post.distns(model.gam=precipf.fire.gam2, model.name="precipf", n=50, newdata=df.precipf[!df.precipf$Extent==ext.base, ], vars="x", terms=T)$ci
precipf.fire.post2$fire.scheme <- df.precipf[!df.precipf$Extent==ext.base, "fire.scheme"]
precipf.fire.post2$fire.extent <- df.precipf[!df.precipf$Extent==ext.base, "fire.extent"]
precipf.fire.post2$Model      <- df.precipf[!df.precipf$Extent==ext.base, "Model"     ]  
precipf.fire.post2$Effect <- as.factor("precipf")
summary(precipf.fire.post2)

fire.post2 <- rbind(tair.fire.post2, precipf.fire.post2, co2.fire.post2)
summary(fire.post2)

pdf(file.path(fig.dir, "Sensitivity_ExtentDeviation_FireScheme.pdf"))
{
  print(
    ggplot(data=ci.terms[!is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & !ci.terms$Extent==ext.base,]) +
      facet_grid(Extent~Effect,scales="free_x") +
      #         geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=fire.scheme), alpha=0.5)+
      geom_line(aes(x=x, y=mean.cent.dev.abs, color=Model), size=2) +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_fill_manual(values=colors.use) +
      scale_color_manual(values=colors.use) +
      theme_bw()
  )
  print(
    ggplot(data=fire.post2[fire.post2$Model%in% c("ed2", "clm.bgc") & !fire.post2$Extent==ext.base, ]) +
      facet_wrap(~Effect,scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=fire.scheme), alpha=0.5)+
      geom_line(aes(x=x, y=mean, color=fire.scheme, linetype=fire.scheme), size=2) +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_fill_manual(values=c("blue3", "red3")) +
      scale_color_manual(values=c("blue3", "red3")) +
      theme_bw()
  )
}
dev.off()
}
# ------------    


}
# --------


# --------
# 5.b.3 Slow processes & sensitivity change: Biomass
#  Hypothesis: Models with greater changes in biomass have greater differences
#              among scales because change in biomass is a slow stabilizing process 
#              that can't happen at short time scales
# --------
{

  # ------------
  # Correlation of Mean Biomass with sensitivity across scales
  # NOTE: Can't do this because "Biomass" means different things for differnet models!!
  # ------------
{
  # Correlation of Mean Biomass with relative sensitivity across all scales
  pdf(file.path(fig.dir, "Sensitivity_Slope_vs_Biomass.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Biomass.mean, y=abs(mean.cent.deriv), color=Extent), size=5) +
    stat_smooth(aes(x=Biomass.mean, y=abs(mean.cent.deriv), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.bm.lm1     <- lm(abs(mean.cent.deriv) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lm1    <- lm(abs(mean.cent.deriv) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lm1 <- lm(abs(mean.cent.deriv) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lm1)
  summary(tair.bm.lm1)
  summary(precipf.bm.lm1)
  
  summary(ci.terms[ci.terms$Model=="ed2",])
  summary(ci.terms[ci.terms$Model=="ed2",])
  
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_vs_Biomass.pdf"))
  ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Biomass.mean, y=abs(mean.rel.cent), color=Extent), size=5) +
    stat_smooth(aes(x=Biomass.mean, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.bm.lm2     <- lm(abs(mean.rel.cent) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lm2    <- lm(abs(mean.rel.cent) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lm2 <- lm(abs(mean.rel.cent) ~ Biomass.mean*(Extent-1) - Biomass.mean, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lm2)
  summary(tair.bm.lm2)
  summary(precipf.bm.lm2)
  
  
  co2.bm.lme     <- lme(mean.rel.cent ~ Biomass.mean, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lme    <- lme(mean.rel.cent ~ Biomass.mean, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lme <- lme(mean.rel.cent ~ Biomass.mean, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lme)
  summary(tair.bm.lme)
  summary(precipf.bm.lme)

  # ------------

}
# ------------

  
# ------------
# Correlation of Biomass *Variability* (% change through time) with sensitivity across scales
# ------------
{
# Correlation of Biomass with relative sensitivity across all scales
pdf(file.path(fig.dir, "Sensitivity_Slope_vs_BiomassVar.pdf"))
ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
  facet_wrap(~Effect, scales="fixed") +
  geom_point(aes(x=Biomass.rel.sd, y=abs(mean.rel.cent), color=Extent), size=5) +
  stat_smooth(aes(x=Biomass.rel.sd, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=c("blue3", "red3", "green3")) +
  scale_fill_manual(values=c("blue3", "red3", "green3")) +
  theme_bw()
dev.off()

co2.bm.lm1     <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
tair.bm.lm1    <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
precipf.bm.lm1 <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
summary(co2.bm.lm1)
summary(tair.bm.lm1)
summary(precipf.bm.lm1)

summary(ci.terms[ci.terms$Model=="ed2",])
summary(ci.terms[ci.terms$Model=="ed2",])

# Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
pdf(file.path(fig.dir, "Sensitivity_vs_BiomassVar.pdf"))
ggplot(data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
  facet_wrap(~Effect, scales="fixed") +
  geom_point(aes(x=Biomass.rel.sd, y=abs(mean.rel.cent), color=Extent), size=5) +
  stat_smooth(aes(x=Biomass.rel.sd, y=abs(mean.rel.cent), color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
  scale_color_manual(values=c("blue", "red3", "green3")) +
  scale_fill_manual(values=c("blue", "red3", "green3")) +
  theme_bw()
dev.off()

co2.bm.lm2     <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
tair.bm.lm2    <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
precipf.bm.lm2 <- lm(abs(mean.rel.cent) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
summary(co2.bm.lm2)
summary(tair.bm.lm2)
summary(precipf.bm.lm2)

bm.precip <- data.frame(char="bm.sd", 
                          effect="precip", 
                          extent=c("1980-2010","1901-2010", "850-2010"),
                          estimate=summary(precipf.bm.lm2)$coefficients[4:6,1],
                          std.err = summary(precipf.bm.lm2)$coefficients[4:6,2],
                          t.stat  = summary(precipf.bm.lm2)$coefficients[4:6,3],
                          p.value = summary(precipf.bm.lm2)$coefficients[4:6,4])
bm.tair <- data.frame(char="bm.sd", 
                        effect="tair", 
                        extent=c("1980-2010","1901-2010", "850-2010"),
                        estimate=summary(tair.bm.lm2)$coefficients[4:6,1],
                        std.err = summary(tair.bm.lm2)$coefficients[4:6,2],
                        t.stat  = summary(tair.bm.lm2)$coefficients[4:6,3],
                        p.value = summary(tair.bm.lm2)$coefficients[4:6,4])
bm.co2 <- data.frame(char="bm.sd", 
                       effect="co2", 
                       extent=c("1980-2010","1901-2010", "850-2010"),
                       estimate=summary(co2.bm.lm2)$coefficients[4:6,1],
                       std.err = summary(co2.bm.lm2)$coefficients[4:6,2],
                       t.stat  = summary(co2.bm.lm2)$coefficients[4:6,3],
                       p.value = summary(co2.bm.lm2)$coefficients[4:6,4])

bm.stats <- rbind(bm.tair, bm.precip, bm.co2)



co2.bm.lme     <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
tair.bm.lme    <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
precipf.bm.lme <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
summary(co2.bm.lme)
summary(tair.bm.lme)
summary(precipf.bm.lme)
}
# ------------

# ------------
# Correlation of Biomass *Variability* (% change through time) with Scalability of short-term sensitivity
# ------------
{
  # Correlation of Biomass with relative sensitivity across all scales
  pdf(file.path(fig.dir, "Sensitivity_Change_Slope_vs_BiomassVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & !ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Biomass.rel.sd, y=mean.deriv.dev.abs, color=Extent), size=5) +
    stat_smooth(aes(x=Biomass.rel.sd, y=mean.deriv.dev.abs, color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("red3", "green3")) +
    scale_fill_manual(values=c("red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.bm.lm1     <- lm(abs(mean.deriv.dev.abs) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lm1    <- lm(abs(mean.deriv.dev.abs) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lm1 <- lm(abs(mean.deriv.dev.abs) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.cent.deriv) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lm1)
  summary(tair.bm.lm1)
  summary(precipf.bm.lm1)
  
  summary(ci.terms[ci.terms$Model=="ed2",])
  summary(ci.terms[ci.terms$Model=="jules.stat",])
  summary(ci.terms[ci.terms$Biomass.rel.sd>=0.55,])
  summary(ci.terms[ci.terms$Biomass.rel.sd<=0.35,])
  summary(ci.terms$Biomass.rel.sd)
  
  # Looking at the asbolute value of the sensitivity (closer to 0 = less sensitivite)
  pdf(file.path(fig.dir, "Sensitivity_Change_vs_BiomassVar.pdf"))
  ggplot(data=ci.terms[!ci.terms$Extent==ext.base & !ci.terms$Quantile=="Other"  & ci.terms$Effect %in% c("tair", "precipf", "CO2"),]) +
    facet_wrap(~Effect, scales="fixed") +
    geom_point(aes(x=Biomass.rel.sd, y=mean.cent.dev.abs, color=Extent), size=5) +
    stat_smooth(aes(x=Biomass.rel.sd, y=mean.cent.dev.abs, color=Extent, fill=Extent, linetype=Extent), method="lm", size=2) +
    scale_color_manual(values=c("blue", "red3", "green3")) +
    scale_fill_manual(values=c("blue", "red3", "green3")) +
    theme_bw()
  dev.off()
  
  co2.bm.lm2     <- lm(mean.cent.dev.abs ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lm2    <- lm(mean.cent.dev.abs ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd , data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lm2 <- lm(mean.cent.dev.abs ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lm2)
  summary(tair.bm.lm2)
  summary(precipf.bm.lm2)
  
  
  co2.bm.lme     <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  tair.bm.lme    <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  precipf.bm.lme <- lme(mean.rel.cent ~ Biomass.sd, random=list(Extent=~1, x=~1), data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.deriv.dev.abs) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
  summary(co2.bm.lme)
  summary(tair.bm.lme)
  summary(precipf.bm.lme)
}
# ------------


}
# --------


}
# -----------------------
}
# ----------------------------------------

veg.stats
evg.stats
fire.stats
bm.stats

model.stats <- rbind(veg.stats, evg.stats, fire.stats, bm.stats)
model.stats$effect.ext <- as.factor(paste(model.stats$effect, model.stats$extent, sep="."))
model.stats[,c("estimate", "std.err", "t.stat", "p.value")] <- round(model.stats[,c("estimate", "std.err", "t.stat", "p.value")], 2)
model.stats$std.err <- round(model.stats[,c("std.err")], 2)
model.stats$sig     <- ifelse(model.stats$p.value<0.05, "*", "")
model.stats$est.err <- paste(model.stats$estimate, "+/-", model.stats$std.err, " ", model.stats$sig)
model.stats
write.csv(model.stats, file.path(out.dir, "ModelCharStats_full.csv"), row.names=F)


library(reshape2)
table3 <- merge(unique(model.stats$char), unique(model.stats$effect))
names(table3) <- c("Character", "Effect")
table3 <- table3[,c(2,1)]
table3

for(i in unique(table3$Character)){
  for(j in unique(table3$Effect)){
    for(e in unique(model.stats$extent)){
      table3[table3$Character==i & table3$Effect==j, e] <- model.stats[model.stats$char==i & model.stats$effect==j & model.stats$extent==e, "est.err"]
    }
  }
}
table3

write.csv(table3, file.path(out.dir, "Table3_ModelCharStats.csv"), row.names=F)
