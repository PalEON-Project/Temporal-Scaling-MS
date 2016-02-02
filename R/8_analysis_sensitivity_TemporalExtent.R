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
{
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]

# Colors for the rest of the script
models.use <- unique(ecosys[,"Model.Order"])
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray40"))

# # Load the statistical model results
# load(file.path(in.base, "gamm_TempExtent.Rdata"))
# 
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
#                      )
# summary(var.ci)
# 
# # Making a climate series so we can emulate the ecosystem dynamics using a base series
# #  Note: this will assume the same AGB trajectories; just changing the NPP
# clim.subset <- (dat.ecosys$Model=="ed2" & dat.ecosys$Extent=="850-2010" & dat.ecosys$Resolution=="t.001")
# var.ts      <- data.frame(Site    = dat.ecosys[clim.subset,"Site"   ], 
#                           PlotID  = dat.ecosys[clim.subset,"PlotID" ], 
#                           Year    = dat.ecosys[clim.subset,"Year"   ], 
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
#    		# Setting up some stuff
#    		dat.subset <- dat.ecosys$Model==m & dat.ecosys$Extent==e & dat.ecosys$Resolution=="t.001"
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
# 		                     Time       = mean(dat.ecosys[dat.ecosys$Model==m, "Time"]) 
# 		                      )
# 
# 		# Adding in the met data that has the full range of values
# 		ci.dat <- cbind(ci.dat, var.ci)
# 
# 		# doing the sensitivity calculations
# 		terms.out   <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ci.dat, 
# 		                           vars=c("tair", "precipf", "CO2"), terms=T)
# 
# 		 # Adding the new sims to the list
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
# 		                      Site       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Site"  ],
# 		                      PlotID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","PlotID"],
# 		                      TreeID     = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","TreeID"],
# 		                      Time       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Time"  ],
# 		                      Year       = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Year"  ],
# 		                      Y          = dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent=="850-2010","Y"     ]
# 		                      )
# 		 
# 		 # Adding in the met data that has the full range of values
# 		 ts.dat <- merge(ts.dat, var.ts, all.x=T, all.y=T)	
# 
# 		 #NOTE: Jules is missing 2010, so we're just going to pretend it's exactly the same as 2009
# 		 if(substr(m, 1, 5)=="jules"){
# 		 	ts.dat[ts.dat$Year==2010,c("Y", "Time")] <- ts.dat[ts.dat$Year==2009,c("Y", "Time")]
# 		 }
# 
# 
# 		 # Doing the emulation of NPP & what drives it
# 		 pred.out    <- post.distns(model.gam=gam.now, model.name=m, n=n, newdata=ts.dat, 
# 		                            vars=c("tair", "precipf", "CO2", "Year", "Time"), terms=F) 
# 		 weight.out  <- factor.weights(model.gam=gam.now, model.name=m, newdata=ts.dat,extent=e, 
# 		                          vars=c("tair", "precipf", "CO2", "Year", "Time")) 
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
# dat.ecosys2c <- merge(dat.ecosys2b, dat.ecosys2[dat.ecosys2$Model=="TreeRingNPP" & dat.ecosys2$Extent=="1985-2010", c(ecosys.vars, "tair", "precipf", "CO2", "Extent", "Resolution")], all.x=T, all.y=T)
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

----------------------------------------
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
# 		if(m=="TreeRingRW") ext="1901-2010" else if(m=="TreeRingNPP") ext="1985-2010" else ext="850-2010"
# 		# ext="1985-2010"
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
ci.terms.graph[ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.00),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(1.00),"lwr.rel"] <- 1.00
ci.terms.graph[ci.terms.graph$upr.rel>(1.00),"upr.rel"] <- 1.00 
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


pdf(file.path(fig.dir, "Fig4_Sensitivity_Rel_Extent.pdf"), height=8.6, width=11)
{
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "gray30"))

print(
ggplot(data=ci.terms.graph[!ci.terms.graph$Effect=="Time" & !ci.terms.graph$Extent=="1985-2010",]) + facet_grid(Extent~Effect, scales="free_x") +
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
}
# -----------------------

# -----------------------
# 5.b. Analysis: 
# -----------------------
# --------
# 5.b.0. Getting the first derivative (first diference) of each line so we can take the mean slope
# --------
{
# First make sure the effects are sorted by x to make this easier
ci.terms <- ci.terms[order(ci.terms$Model,ci.terms$Extent, ci.terms$Effect, ci.terms$x),]
sim.terms <- sim.terms[order(sim.terms$Model, sim.terms$PFT, sim.terms$Effect, sim.terms$x),]
ci.terms[1:20,1:11]
dim(ci.terms)
summary(ci.terms)
ci.terms$Y.type     <- as.factor(ifelse(ci.terms $Model=="TreeRingRW", "RW", "NPP"))
sim.terms$Y.type    <- as.factor(ifelse(sim.terms$Model=="TreeRingRW", "RW", "NPP"))
ci.terms$data.type  <- as.factor(ifelse(substr(ci.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
sim.terms$data.type <- as.factor(ifelse(substr(sim.terms$Model,1,8)=="TreeRing", "Tree Rings", "Model"))

# Looking at the range of the 95% CI
ci.terms$rel.ci.range  <- ci.terms$upr.rel - ci.terms$lwr.rel

# Making a new dataframe dedicated to the derivatives
cols.sims <- which(substr(names(sim.terms),1,1)=="X")
sim.deriv <- sim.terms[,]
sim.deriv[,cols.sims] <- NA

for(e in unique(ci.terms$Effect)){
  for(m in unique(ci.terms$Model)){
    for(p in unique(ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "Extent"])){
      x.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "x"], lag=1), NA)
      y.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "mean.rel"], lag=1), NA)
      ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "deriv"] <- y.dif/x.dif
      
      # For the full simiulation for robust analysis
      y.dif2 <- rbind(apply(sim.terms[sim.terms$Model==m & sim.terms$Effect==e & sim.terms$Extent==p, cols.sims], 2, FUN=diff), NA)
      x.dif <- c(diff(sim.terms[sim.terms$Model==m & sim.terms$Effect==e  & sim.terms$Extent==p, "x"], lag=1), NA)
      sim.deriv[sim.deriv$Model==m & sim.deriv$Effect==e & sim.terms$Extent==p, cols.sims] <- apply(y.dif2, 2, FUN=function(y){y/x.dif})
    }
  } 
}
summary(ci.terms)

# Stacking and aggregating the simulations
deriv.stack <- stack(sim.deriv[,cols.sims])
names(deriv.stack)  <- c("deriv", "sim")
deriv.stack[,names(sim.deriv)[which(!(1:ncol(sim.deriv)) %in% cols.sims)]] <- sim.deriv[,which(!(1:ncol(sim.deriv)) %in% cols.sims)]
deriv.stack$Model  <- as.factor(deriv.stack$Model)
deriv.stack$Effect <- as.factor(deriv.stack$Effect)
summary(deriv.stack)

sim.stack <- stack(sim.terms[,cols.sims])
names(sim.stack)  <- c("Y", "sim")
sim.stack[,names(sim.terms)[which(!(1:ncol(sim.terms)) %in% cols.sims)]] <- sim.terms[,which(!(1:ncol(sim.terms)) %in% cols.sims)]
sim.stack$Model  <- as.factor(sim.stack$Model)
sim.stack$Effect <- as.factor(sim.stack$Effect)
summary(sim.stack)

# Getting the mean derivative for each model & effet & Extent
# -- Need to do this in two bins: overlap climate, extremes
for(e in unique(ci.terms$Effect)){
  range.modern <- range(dat.ecosys[dat.ecosys$Extent=="1901-2010", e], na.rm=T)
  ci.terms[ci.terms$Effect==e, "Climate"] <- as.factor(ifelse(ci.terms[ci.terms$Effect==e, "x"]>=range.modern[1] & 
                                                                ci.terms[ci.terms$Effect==e, "x"]<=range.modern[2],
                                                       "Overlap", "Extended"))
  sim.stack[sim.stack$Effect==e, "Climate"] <- as.factor(ifelse(sim.stack[sim.stack$Effect==e, "x"]>=range.modern[1] & 
                                                                  sim.stack[sim.stack$Effect==e, "x"]<=range.modern[2],
                                                              "Overlap", "Extended"))
  deriv.stack[deriv.stack$Effect==e, "Climate"] <- as.factor(ifelse(deriv.stack[deriv.stack$Effect==e, "x"]>=range.modern[1] & 
                                                                      deriv.stack[deriv.stack$Effect==e, "x"]<=range.modern[2],
                                                              "Overlap", "Extended"))
}
summary(ci.terms)
summary(sim.stack)
summary(deriv.stack)

vars.deriv <- c("deriv", "x")
deriv.agg                             <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                   FUN=mean, na.rm=T)
deriv.agg[,paste0(vars.deriv, ".sd")] <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                   FUN=sd, na.rm=T)[,vars.deriv]
summary(deriv.agg)

vars.sim <- c("Y", "x")
sim.agg                            <- aggregate(sim.stack[,vars.sim], 
                                                 by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                 FUN=mean, na.rm=T)
sim.agg[,paste0(vars.sim, ".sd")]  <- aggregate(sim.stack[,vars.sim], 
                                               by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                               FUN=sd, na.rm=T)[,vars.sim]
summary(sim.agg)


vars.agg <- c("x", "mean.rel", "lwr.rel", "upr.rel", "rel.ci.range", "deriv")
ci.terms.agg <- aggregate(ci.terms[,vars.agg],
                          by=ci.terms[,c("Model", "Y.type", "data.type", "Extent", "Effect", "Climate")],
                          FUN=mean, na.rm=T)

# Ignoring 1985-2010 for analyses & presentation
ci.terms.agg <- ci.terms.agg[!ci.terms.agg$Extent=="1985-2010",]
sim.agg      <- sim.agg[!sim.agg$Extent=="1985-2010",]
deriv.agg    <- deriv.agg[!deriv.agg$Extent=="1985-2010",]
summary(ci.terms.agg)
summary(sim.agg)
summary(deriv.agg)

}
# --------

# --------
# 5.b.1. Looking at changes in sensitivity based on extent
# --------
library(nlme)
{
summary(ci.terms.agg)
ci.terms <- ci.terms[!ci.terms$Extent=="1985-2010",]
sim.stack <- sim.stack[!sim.stack$Extent=="1985-2010",]

# Tair
tair.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
summary(tair.overlap); anova(tair.overlap)
tair.overlap2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
summary(tair.overlap2); anova(tair.overlap2)

tair.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
summary(tair.overlap3); anova(tair.overlap3)
tair.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
summary(tair.overlap4); anova(tair.overlap4)


tair.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Extended",])
summary(tair.extended); anova(tair.extended)
tair.extended2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Extended",])
summary(tair.extended2); anova(tair.extended2)

# Precipf
precipf.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
summary(precipf.overlap); anova(precipf.overlap)
precipf.overlap2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
summary(precipf.overlap2); anova(precipf.overlap2)

precipf.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
summary(precipf.overlap3); anova(precipf.overlap3)
precipf.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
summary(precipf.overlap4); anova(precipf.overlap4)

precipf.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Extended",])
summary(precipf.extended); anova(precipf.extended)

# CO2
CO2.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
summary(CO2.overlap); anova(CO2.overlap)

CO2.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
summary(CO2.overlap3); anova(CO2.overlap3)
CO2.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
summary(CO2.overlap4); anova(CO2.overlap4)

CO2.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Extended",])
summary(CO2.extended); anova(CO2.extended)
}
# --------

# --------
# 5.b.2. Looking explicitly at change in sensitivity among models
# --------
{
summary(deriv.stack)
deriv.stack2 <- deriv.stack[deriv.stack$Extent=="850-2010",]
summary(deriv.stack2)

for(m in unique(deriv.stack2$Model)){
  deriv.stack2[deriv.stack2$Model==m, "dDeriv"] <- deriv.stack[deriv.stack$Model==m & deriv.stack$Extent=="850-2010","deriv"] - deriv.stack[deriv.stack$Model==m & deriv.stack$Extent=="1901-2010","deriv"]
}
summary(deriv.stack2)

# Aggregating a bit; ignoring the periods of climate overlap
deriv.agg2 <- aggregate(deriv.stack2[,c("deriv", "dDeriv", "x")], 
                        by=deriv.stack2[,c("Model", "Y.type", "data.type", "Effect", "sim")], 
                        FUN=mean, na.rm=T)
mean(abs(deriv.agg2[deriv.agg2$Effect=="CO2","dDeriv"]), na.rm=T); sd(abs(deriv.agg2[deriv.agg3$Effect=="CO2","dDeriv"]), na.rm=T)
mean(deriv.agg2[deriv.agg2$Effect=="tair","dDeriv"], na.rm=T); sd(deriv.agg2[deriv.agg3$Effect=="tair","dDeriv"], na.rm=T)
mean(deriv.agg2[deriv.agg2$Effect=="precipf","dDeriv"], na.rm=T); sd(deriv.agg2[deriv.agg3$Effect=="precipf","dDeriv"], na.rm=T)



deriv.agg3 <- aggregate(deriv.agg2[,c("deriv", "dDeriv", "x")], 
                        by=deriv.agg2[,c("Model", "Y.type", "data.type", "Effect")], 
                        FUN=mean, na.rm=T)
deriv.agg3
mean(deriv.agg3[deriv.agg3$Effect=="CO2","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="CO2","dDeriv"], na.rm=T)
mean(deriv.agg3[deriv.agg3$Effect=="tair","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="tair","dDeriv"], na.rm=T)
mean(deriv.agg3[deriv.agg3$Effect=="precipf","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="precipf","dDeriv"], na.rm=T)

d.co2.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="CO2",])
d.tair.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="tair",])
d.precipf.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="precipf",])
summary(d.co2.mod)
summary(d.tair.mod)
summary(d.precipf.mod)

d.co2.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="CO2",])
d.tair.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="tair",])
d.precipf.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="precipf",])
summary(d.co2.mod2)
summary(d.tair.mod2)
summary(d.precipf.mod2)
}
# --------

# --------
# 5.b.3. Attributing change in sensitivity to change in biomass and composition
# --------
{
summary(ecosys)
vars.agg <- c("NPP", "AGB", "LAI", "Evergreen", "Deciduous", "Grass")
# Calculating precent change in time in the two reference periods
mod.agg1 <- aggregate(ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=850,vars.agg], 
                     by=ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=850,c("Model", "Model.Order")],
                     FUN=mean, na.rm=T)
mod.agg2 <- aggregate(ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=1901,vars.agg], 
                      by=ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=1901,c("Model", "Model.Order")],
                      FUN=mean, na.rm=T)
mod.agg1$Extent <- "850-2010"
mod.agg2$Extent <- "1901-2010"

mod.agg <- mod.agg1[,1:2]
mod.agg[,vars.agg] <- (mod.agg1[,vars.agg]-mod.agg2[,vars.agg])/mod.agg1[,vars.agg]

mod.agg[mod.agg$Model=="jules.stat", "AGB"] <- mod.agg[mod.agg$Model=="jules.stat", "LAI"] 
mod.agg[is.na(mod.agg)] <- 0
mod.agg

# Merging the percent change in time with 
deriv.agg3 <- merge(deriv.agg3, mod.agg, all.x=T, all.y=T)
deriv.agg3

effect.AGB <- lm(abs(dDeriv) ~ AGB*(Effect-1)-Effect-AGB, data=deriv.agg3[,])
summary(effect.AGB)

effect.evg <- lm(abs(dDeriv) ~ Evergreen*(Effect-1)-Effect-Evergreen, data=deriv.agg3[,])
summary(effect.evg)
}
# --------


}
# -----------------------
# ----------------------------------------

# ----------------------------------------
# 6. Graphing & Analyzing Emulated Change in NPP
# ----------------------------------------
# -----------------------
# 6.a. Graphing
# -----------------------
wt.terms2 <- merge(wt.terms, dat.ecosys2[dat.ecosys2$Extent=="850-2010",c("Model", "Model.Order", "Y.type", "Site", "data.type", "Year", "Y", "Y.rel", "Y.10", "Y.rel.10")], all.x=T, all.y=F)
{
# summary(wt.terms2)
dim(wt.terms); dim(wt.terms2)
# --------------------------
# Adjusting CO2 Effect
# --------------------------
# Note: because the gam makes the smoother cross 0 at the MEAN CO2 (which is in the 1800s), 
# it's saying the region is pretty CO2-limited at periods where that doesn't really make 
# sense, so we're going to off relativize it to whatever the starting point for the run is
# --------------------------
{
for(m in unique(wt.terms2$Model)){
	for(e in unique(wt.terms2[wt.terms2$Model==m, "Extent"])){
		yr1       <- min(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "Year"]) # find the minimum year
		yr2       <- min(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & !is.na(wt.terms2$weight.CO2), "Year"]) # find the minimum year
		co2.base <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & wt.terms2$Year<=(yr1+5), "fit.CO2"], na.rm=T) # mean of the 5 years around the starting 
		co2.base.10 <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & wt.terms2$Year<=(yr2+5),"fit.CO2.10"],na.rm=T) # mean of the 5 years around the starting point
		wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.adj"] <- wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2"] - co2.base
		wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.10.adj"] <- wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.10"] - co2.base.10
	}
}
summary(wt.terms2)

wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")] <- abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")])/rowSums(abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")]))
wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")] <- abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")])/rowSums(abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")]))

wt.terms2[is.na(wt.terms2$weight.tair.10), c("weight.tair.10", "weight.precipf.10", "weight.CO2.10")] <- 0
wt.terms2[is.na(wt.terms2$weight.tair.10.adj), c("weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")] <- 0

summary(rowSums(wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")]))

summary(wt.terms2)
}
# -----------
# Aggregating to the region ensemble-level
# -----------
{
# factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "Y.rel", "Y.rel.10", "weight.tair.adj", "weight.precipf.adj", "weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")

# aggregate across plots
ensemble.wts0 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=mean, na.rm=T)
ensemble.wts0[,paste0(factors.aggregate, ".lo")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=quantile, 0.025, na.rm=T)[,factors.aggregate[1:5]]
ensemble.wts0[,paste0(factors.aggregate, ".hi")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=quantile, 0.975, na.rm=T)[,factors.aggregate[1:5]]
summary(ensemble.wts0)

ensemble.wts1 <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

ensemble.wts.lo <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=quantile, 0.025, na.rm=T)
names(ensemble.wts.lo)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[5:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=quantile, 0.975, na.rm=T)
names(ensemble.wts.hi)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[5:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

ensemble.wts.final <- cbind(ensemble.wts1, ensemble.wts.lo[5:ncol(ensemble.wts.lo)], ensemble.wts.hi[5:ncol(ensemble.wts.hi)])
summary(ensemble.wts.final)
}
# -----------

# -----------
# Aggregating to the site ensemble-level
# -----------
{
# factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "Y.rel", "Y.rel.10", "weight.tair.adj", "weight.precipf.adj", "weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")

# aggregate across plots
ensemble.wts.site                                    <- aggregate(wt.terms2[,factors.aggregate], 
                                                                  by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                  FUN=mean, na.rm=T)
ensemble.wts.site[,paste0(factors.aggregate, ".lo")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], 
                                                                  by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                  FUN=quantile, 0.025, na.rm=T)[,factors.aggregate[1:5]]
ensemble.wts.site[,paste0(factors.aggregate, ".hi")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], 
                                                                  by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                  FUN=quantile, 0.975, na.rm=T)[,factors.aggregate[1:5]]
summary(ensemble.wts.site)

}
# -----------

# -----------
# Re-normalizing factor weights
# -----------
{
summary(rowSums(ensemble.wts.final[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")]))


# {
# wts.sum.10 <- abs(ensemble.wts.final$weight.tair.10) + abs(ensemble.wts.final$weight.precipf.10) + abs(ensemble.wts.final$weight.CO2.10)
# ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")] <- ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")]/wts.sum.10
# ensemble.wts.final[is.na(ensemble.wts.final$weight.tair.10   ),"weight.tair.10"   ] <- 0
# ensemble.wts.final[is.na(ensemble.wts.final$weight.precipf.10),"weight.precipf.10"] <- 0
# ensemble.wts.final[is.na(ensemble.wts.final$weight.CO2.10    ),"weight.CO2.10"    ] <- 0


# wts.sum <- abs(ensemble.wts.final$weight.tair) + abs(ensemble.wts.final$weight.precipf) + abs(ensemble.wts.final$weight.CO2)
# ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")] <- ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")]/wts.sum
# summary(ensemble.wts.final)
# }
# -----------

summary(ensemble.wts.final)

ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
# ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$ci.min)
# ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
}

# -----------
# Region figures
# -----------
{
pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_1500-2010_Annual_AllExtent.pdf"), width=11, height=8.5)
{
print(
ggplot(ensemble.wts.final) + facet_grid(Extent~., scales="fixed") +
	scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
 	geom_ribbon(data= ensemble.wts.final[,], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, y=fit.graph*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1985-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1985-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1985-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
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
		  # axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
dev.off()

pdf(file.path(fig.dir, "Fig5_Ensemble_Drivers_Time_Region_1500-2010_Annual_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0.5,0.5,ensemble.wts.final$fit.full.rel.10.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
{
print(
ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
	scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
 	geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
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
		  # axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
dev.off()

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_0850-2010_Decadal_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
{
print(
ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
	scale_x_continuous(limits=c(850,2010), expand=c(0,0)) +
 	geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
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
		  # axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(900,1300), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      ggtitle("Medieval Warm Period") + 
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
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1400,1800), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      ggtitle("Little Ice Age") + 
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
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
dev.off()

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_Annual_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.hi>2,2,ensemble.wts.final$fit.full.rel.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.lo<0,0,ensemble.wts.final$fit.full.rel.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel<0,NA,ensemble.wts.final$fit.full.rel)
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1850,2010), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
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
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1750,1850), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
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
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
dev.off()

}
# -----------

# -----------
# Site Figures
# -----------
pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Sites_1500-2010_Decadal.pdf"), width=11, height=8.5)
ensemble.wts.site$ci.max    <- ifelse(ensemble.wts.site$fit.full.rel.10.hi > 2  , 2  , ensemble.wts.site $fit.full.rel.10.hi)
ensemble.wts.site$ci.min    <- ifelse(ensemble.wts.site$fit.full.rel.10.lo < 0.25, 0.25, ensemble.wts.site $fit.full.rel.10.lo)
ensemble.wts.site$fit.graph <- ifelse(ensemble.wts.site$fit.full.rel.10   < 0.25   , NA , ensemble.wts.site $fit.full.rel.10)
for(s in unique(ensemble.wts.site$Site)){
print(
ggplot(ensemble.wts.site[ensemble.wts.site$Site==s & !ensemble.wts.site$Extent=="1985-2010",]) + facet_grid(Extent~., scales="fixed") +
	scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
 	geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Site==s & !ensemble.wts.site$Extent=="1985-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="1901-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="1901-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="1901-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
	ggtitle(s) + 
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
		  # axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
dev.off()


pdf(file.path(fig.dir, "Ensemble_Drivers_Time_AllSites_1500-2010_Decadal.pdf"), width=11, height=8.5)
{
print(
ggplot(ensemble.wts.site[ensemble.wts.site$Extent=="850-2010",]) + 
	facet_grid(Site~., scales="fixed") +
	scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
 	geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Extent=="850-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
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
		  # axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
}
dev.off()
# -----------

}

# -----------------------

# -----------------------
# 6.b. Analysis: 
# -----------------------
summary(ensemble.wts0)
summary(ensemble.wts.final)

weight.CO2     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010",])
weight.tair    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010",])
weight.precipf <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010",])
summary(weight.CO2)
summary(weight.tair)
summary(weight.precipf)


weight.CO22     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year>=1850,])
weight.tair2    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year>=1850,])
weight.precipf2 <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year>=1850,])
summary(weight.CO22)
summary(weight.tair2)
summary(weight.precipf2)

weight.CO22     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year<=1850,])
weight.tair2    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year<=1850,])
weight.precipf2 <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010" & ensemble.wts0$Year<=1850,])
summary(weight.CO22)
summary(weight.tair2)
summary(weight.precipf2)

summary(ensemble.wts0)
mean(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])

mean(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])

mean(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
quantile(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel.10"], c(0.025, 0.975))


mean(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
quantile(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel.10"], c(0.025, 0.975))

vars.compare <- c("Model", "Year", "fit.full.rel", "weight.CO2.adj", "weight.tair.adj", "weight.precipf.adj")
ensemble.ext1 <- ensemble.wts0[ensemble.wts0$Extent=="850-2010",vars.compare]
names(ensemble.ext1)[3:length(vars.compare)] <- paste0(vars.compare[3:length(vars.compare)],".850")
ensemble.ext2 <- ensemble.wts0[ensemble.wts0$Extent=="1901-2010",vars.compare]
names(ensemble.ext2)[3:length(vars.compare)] <- paste0(vars.compare[3:length(vars.compare)],".1901")

ensemble.comparison <- merge(ensemble.ext1, ensemble.ext2, all.x=T, all.y=T)
summary(ensemble.comparison)

ggplot(data=ensemble.comparison) + 
  facet_wrap(~Model, scales="fixed") +
  geom_point(aes(x= fit.full.rel.850, y= fit.full.rel.1901)) +
  coord_cartesian() +
  theme_bw()


plot(abs(weight.CO2.adj) ~ Extent, data=ensemble.wts0[!ensemble.wts0$Extent=="1985-2010",])
# -----------------------
# ----------------------------------------
