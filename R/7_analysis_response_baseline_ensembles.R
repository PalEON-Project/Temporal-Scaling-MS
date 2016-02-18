# ----------------------------------------
# Objective: Compare NPP deviations through time & driver of "wiggles" among 
#            models, tree ring NPP, and tree ring RW products
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
#
# -------------------------
# Hypotheses: 
# -------------------------
# 1. Models that have more similar NPP or NPP variability to the data will 
#    have greater agreement regarding climate sensitivity and what drives 
#    change through time
# 2. Models that have land use & realistic disturbance regimes will more 
#    closely resemble the empirical data
# -------------------------

# -------------------------
# Workflow
# -------------------------
# 1. Set Directories
# 2. Load data files & function scripts
# 3. Standardize driver responses to the mean model NPP to aid comparison (loop)
#    a. Find the NPP to relativize each set off of
#    b. Relativizing everything in dat.ecosys to make it comparable to tree rings
#    c. Finding the percent change in NPP relative to the mean for that particular scale
#    d. Relativizing the factor fits through times and weights as well
# 4. Graphing & Analyzing Raw NPP & RW
#    a. Graphing
#       1. All Sites 1800-2010 (Fig. 1)
#       2. Models Only, All Sites 0850-2010 (Supplemental Figure 1)
#       3. Just PHA 1900-2010
#       4. PHA, PHO (locations with NPP products, 1950-2010)
#    b. Summary Statistics
#       1. Mean NPP, range of NPP variability within & among models/data
# 5. Graphing & Analyzing Sensitivity
#    a. Graphing
#    b. Quantitative Analysis
# 6. Graphing & Analyzing ensemble means of Drivers of change and sensitivities through time
#    a. Graphing
#    b. Quantitative Analysis
# -------------------------
# ----------------------------------------
rm(list=ls())
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
path.data <- "Data"
in.base <- "Data/gamms/Sensitivity_Baseline"
out.dir <- "Data/analyses/analysis_baseline"
fig.dir <- "Figures/analyses/analysis_baseline"

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
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray30"))

# # Load the statistical model results
# load(file.path(in.base, "gamm_baseline.Rdata"))
# 
# dat.ecosys <- cbind(mod.out$data[, ], mod.out$ci.response[,c("mean", "lwr", "upr")])
# ci.terms   <- mod.out$ci.terms[,]
# wt.terms   <- mod.out$weights[,]
# sim.terms  <- mod.out$sim.terms
# sim.terms $Effect <- as.factor(sim.terms$Effect)
# 
# # Grouping the kind and source of the data
# dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRW", "RW", "NPP"))
# ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRW", "RW", "NPP"))
# wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRW", "RW", "NPP"))
# sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRW", "RW", "NPP"))
# 
# dat.ecosys$data.type <- as.factor(ifelse(substr(dat.ecosys$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# ci.terms  $data.type <- as.factor(ifelse(substr(ci.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# 
# summary(ci.terms)
# summary(dat.ecosys)
# summary(wt.terms)
# summary(sim.terms[,1:10])
}
# ----------------------------------------


# ----------------------------------------
# 3. Standardize driver responses to the mean model NPP to facilitate comparisons
#	 -- Two things for standardization and graphing:
# 		1. Relative by model-mean NPP
#       2. Decadal smoothing to help show generalized patterns
# ----------------------------------------
{
# # Across all scales (resolution) finding the mean NPP
# # NOTE: we ARE relativizing per site here since the response curves were site-specific
# 
# # Make sure all data sets are ordered by year, then treeID, then plotID, then Model
# # sort.order <- c("Model", "PlotID", "TreeID", "Year")
# dat.ecosys <- dat.ecosys[order(dat.ecosys$Model, dat.ecosys$PlotID, dat.ecosys$TreeID, dat.ecosys$Year),]
# wt.terms   <- wt.terms[order(wt.terms$Model, wt.terms$PlotID, wt.terms$TreeID, wt.terms$Year),]
# 
# # Double Check to make sure things are sorted by year so rollapply works
# dat.ecosys[which(dat.ecosys$Model=="TreeRingRW")[1:20],]
# wt.terms  [which(wt.terms  $Model=="TreeRingRW")[1:20],]
# 
# summary(dat.ecosys)
# summary(wt.terms)
# 
# {
# for(m in unique(ci.terms$Model)){
# 
# 		# -----------------------
# 		# 3.a. Find the NPP to relativize each set off of
# 		# Using mean model NPP across sites since the GAMM response curves are for 
# 		#    the whole model & not site-specific are parameterized
# 		# -----------------------
# 		# Find the start year for the extent
# 		# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))
# 
# 		npp <- mean(dat.ecosys[dat.ecosys$Model==m, "Y"], na.rm=T)			
# 		# -----------------------
# 		
# 		# -----------------------
# 		# 3.b Relativizing everything in dat.ecosys to make it comparable to tree rings
# 		# -----------------------
# 		{		
# 		# Which factors to relativize
# 		y.rel <- c("Y", "fit.gam", "mean", "lwr", "upr")
# 
# 		# for some reason, I can create multiple new columns at once
# 		# Solution: use a loop to create blank columns and then fill them
# 		for(y in y.rel){
# 			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel"       )] <- NA	
# 			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".10"        )] <- NA	
# 			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel", ".10")] <- NA	
# 		}
# 		dat.ecosys[dat.ecosys$Model==m,paste0(y.rel, ".rel")] <- dat.ecosys[dat.ecosys$Model==m, y.rel]/npp
# 		
# 		# Getting 10-year running means to make clearer figures
# 		for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
# 			# Note: If we're working with tree ring data, we need to go by plot for NPP products 
# 			#       & by Tree for individual-level tree rings products
# 			if(m=="TreeRingNPP"){
# 				for(p in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "PlotID"])){
# 
# 					# Raw NPP (to add dark line over faded annual wiggles)
# 					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$PlotID==p, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 				}
# 			} else if(m=="TreeRingRW") {
# 				for(t in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "TreeID"])){
# 					# If we have too few data points, we need to skip that tree 
# 					if(length(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel[1]]) < 10) next
# 
# 					# Raw NPP (to add dark line over faded annual wiggles)
# 					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$TreeID==t, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 				}
# 			} else {
# 				# Raw NPP (to add dark line over faded annual wiggles)
# 				dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 
# 				# Relativized NPP (to have generalized patterns for figures)
# 				dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 			}
# 		}
# 		}		
# 		# -----------------------
# 
# 		
# 		# -----------------------
# 		# 3.c. Finding the percent change in NPP relative to the mean for that particular scale
# 		# -----------------------
# 		{
# 		y.rel <- c("mean", "lwr", "upr")
# 		for(y in y.rel){
# 			ci.terms[ci.terms$Model==m,paste0(y, ".rel"       )] <- NA	
# 		}		
# 
# 		ci.terms[ci.terms$Model==m,paste0(y.rel,".rel")] <- ci.terms[ci.terms $Model==m, y.rel]/npp
# 		
# 		# Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
# 		cols.sim <- which(substr(names(sim.terms),1,1)=="X")
# 		sim.terms[sim.terms$Model==m,cols.sim] <- sim.terms[sim.terms$Model==m,cols.sim]/npp
# 		}
# 		# -----------------------
# 
# 		# -----------------------
# 		# 3.d. Relativizing the factor fits through times and weights as well
# 		# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
# 		# -----------------------
# 		{
# 		y.rel <- c("fit.tair", "fit.precipf", "fit.CO2")
# 		# for some reason, I can create multiple new columns at once
# 		# Solution: use a loop to create blank columns and then fill them
# 		for(y in y.rel){
# 			wt.terms[wt.terms$Model==m,paste0(y, ".rel"       )] <- NA	
# 			wt.terms[wt.terms$Model==m,paste0(y, ".rel", ".10")] <- NA	
# 		}
# 		
# 		wt.terms[wt.terms$Model==m,paste0(y.rel, ".rel")] <- 1+(wt.terms[wt.terms$Model==m,y.rel])/npp
# 		
# 		# We only really care about smoothing the relativized weights
# 		y.rel2 <- c(y.rel, paste0(y.rel, ".rel"), "weight.tair", "weight.precipf", "weight.CO2")
# 		for(y in y.rel2){
# 			wt.terms[wt.terms$Model==m,paste0(y, ".10")] <- NA	
# 		}
# 
# 		# Getting 10-year running means to make clearer figures
# 		for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
# 			if(m=="TreeRingNPP"){
# 				for(p in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "PlotID"])){
# 					# Relativized NPP (to have generalized patterns for figures)
# 					wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
# 				}
# 			} else if(m=="TreeRingRW"){
# 				for(t in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "TreeID"])){
# 
# 					# If we have too few data points, we need to skip that tree 
# 					if(length(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2[1]]) < 10) next
# 
# 					# Relativized NPP (to have generalized patterns for figures)
# 					wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
# 				}				
# 			} else {
# 				# Relativized NPP (to have generalized patterns for figures)
# 				wt.terms[wt.terms$Model==m & wt.terms$Site==s,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)
# 			}
# 		}
# 		}
# 		# -----------------------
# 
# }
# } # End section block
# 
# 
# summary(dat.ecosys)
# summary(ci.terms)
# summary(wt.terms)
# summary(sim.terms[,1:10])
# 
# save(dat.ecosys, ci.terms, wt.terms, sim.terms, file=file.path(out.dir, "post-process_baseline.RData"))
}
# ----------------------------------------

load(file.path(out.dir, "post-process_baseline.RData"))

# ----------------------------------------
# 4. Analyzing Region-Level NPP, change, & climate effects
# ----------------------------------------
summary(dat.ecosys)

# ---------------------
# 4.a. Aggregate to region level
#  -- Note: going to site first for tree ring products
# ---------------------
{
summary(dat.ecosys)
summary(wt.terms)

fac.ecosys <- c("Y", "Y.10", "Y.rel", "Y.rel.10", "tair", "precipf", "CO2", "Time")
dat.ecosys2 <- aggregate(dat.ecosys[,fac.ecosys],  
                         by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site", "Year")],
                         FUN=mean, na.rm=T)
summary(dat.ecosys2)

fac.wt <- c("fit.full", "fit.tair", "fit.precipf", "fit.CO2", "fit.tair.10", "fit.precipf.10", "fit.CO2.10", "fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel", "fit.tair.rel.10", "fit.precipf.rel.10", "fit.CO2.rel.10", "weight.tair", "weight.precipf", "weight.CO2", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
wt.terms2 <- aggregate(wt.terms[,fac.wt],  
                       by=wt.terms[,c("Model","Site", "Year")],
                         FUN=mean, na.rm=T)
summary(wt.terms2)


# ----------
# Adjusting CO2 Effect
# ----------
# Note: because the gam makes the smoother cross 0 at the MEAN CO2 (which is in the 1800s), 
# it's saying the region is pretty CO2-limited at periods where that doesn't really make 
# sense, so we're going to off relativize it to whatever the starting point for the run is
# ----------
{
for(m in unique(wt.terms2$Model)){
	yr1       <- min(wt.terms2[wt.terms2$Model==m, "Year"]) # find the minimum year
	yr2       <- min(wt.terms2[wt.terms2$Model==m & !is.na(wt.terms2$weight.CO2.10), "Year"]) # find the minimum year
	co2.base <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr1+5), "fit.CO2"], na.rm=T) # mean of the 5 years around the starting 
	co2.base.10 <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr2+5),"fit.CO2.10"],na.rm=T) # mean of the 5 years around the starting point

	co2.rel.base <- 1-mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr1+5), "fit.CO2.rel"], na.rm=T) # mean of the 5 years around the starting 
	co2.rel.base.10 <- 1-mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr2+5),"fit.CO2.rel.10"],na.rm=T) # mean of the 5 years around the starting point
	
  wt.terms2[wt.terms2$Model==m, "fit.CO2.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2"] - co2.base
	wt.terms2[wt.terms2$Model==m, "fit.CO2.10.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.10"] - co2.base.10
	wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.rel"] + co2.rel.base # The adjustment needs to be relative to 1
	wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.10.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.10"] + co2.rel.base.10 # The adjustment needs ot be relative to 1
}

wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")] <- abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")])/rowSums(abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")]))
wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")] <- abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")])/rowSums(abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")]))

summary(wt.terms2)

wt.check <- rowSums(wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")])
wt.check.10 <- rowSums(wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")])
summary(wt.check)
summary(wt.check.10)
}
# ----------

dat.ecosys2 <- merge(dat.ecosys2, wt.terms2, all.x=T, all.y=T)
summary(dat.ecosys2)

fac.wt <- c(fac.wt, "fit.CO2.rel.adj", "fit.CO2.rel.10.adj")
fac.wt.2 <- c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")
fac.wt.ci <- c("fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel.adj", "fit.tair.rel.10", "fit.precipf.rel.10", "fit.CO2.rel.10.adj")
factors.agg <- c(fac.ecosys, fac.wt.2, fac.wt.ci)

dat.region                             <- aggregate(dat.ecosys2[,factors.agg], 
                                                   by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                   FUN=mean, na.rm=T)
dat.region[,paste0(factors.agg, ".lo")] <- aggregate(dat.ecosys2[,factors.agg], 
                                                    by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                    FUN=quantile, 0.025, na.rm=T)[,factors.agg]
dat.region[,paste0(factors.agg, ".hi")] <- aggregate(dat.ecosys2[,factors.agg], 
                                                    by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                    FUN=quantile, 0.975, na.rm=T)[,factors.agg]
summary(dat.region)


# Ensemble-level stats
# Aggregate models, them merge in tree ring region level
dev.region <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg], 
                                   by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                   FUN=mean, na.rm=T)
dev.region[,paste0(factors.agg, ".lo")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg], 
                                                    by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                                    FUN=quantile, 0.025, na.rm=T)[,factors.agg]
dev.region[,paste0(factors.agg, ".hi")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg],
                                                    by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                                    FUN=quantile, 0.975, na.rm=T)[,factors.agg]


# Adding in the standard deviation for a few other variables we want to use to describe the models through time
factors.sd <- c("Y.rel", "fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel.adj")
dev.region[,paste0(factors.sd, ".sd")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.sd],
                                          by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                          FUN=sd, na.rm=T)[,factors.sd]
summary(dev.region)

dev.region <- merge(dev.region, dat.region[dat.region$data.type=="Tree Rings",], all.x=T, all.y=T)

dev.region[is.na(dev.region$weight.tair.10.adj   ),"weight.tair.10.adj"   ] <- 0
dev.region[is.na(dev.region$weight.precipf.10.adj),"weight.precipf.10.adj"] <- 0
dev.region[is.na(dev.region$weight.CO2.10.adj    ),"weight.CO2.10.adj"    ] <- 0
}
# ---------------------


# ---------------------
# 4.b. Plotting 3 Regional patterns
# 4.b.1. Raw NPP by model
# 4.b.2. NPP deviation by data stream
# 4.b.3. Factor weights
# ---------------------
{
# --------
# 4.b.1. Raw NPP by model
# --------
plot.npp <- {
	ggplot(data=dat.region) + 
	scale_x_continuous(limits=c(1700,2010), expand=c(0,0), name="Year") +
	scale_y_continuous(expand=c(0,0), name="NPP MgC/HA/yr") +
	facet_grid(Y.type~., scales="free_y", space="free") +
	geom_ribbon(aes(x=Year, ymin=Y.10.lo, ymax=Y.10.hi, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=Year, y=Y.10, color=Model.Order, linetype=Model.Order), size=1.5) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	guides(color=guide_legend(nrow=2),
	       fill =guide_legend(nrow=2)) +
	theme(legend.title=element_text(size=rel(1), face="bold"),
	      legend.text=element_text(size=rel(1)),
	      legend.position=c(0.5, 0.9),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(fill=NA, color="black", size=0.5), 
	      panel.background=element_blank(),
	      panel.margin.x=unit(0, "lines"),
	      panel.margin.y=unit(0, "lines"))  +
# 	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
# 		  axis.text.y=element_text(color="black", size=rel(1.5)), 
# 		  axis.title.x=element_text(size=rel(2), face="bold", vjust=-0.4),  
# 		  axis.title.y=element_text(size=rel(2), face="bold"),
# 		  axis.ticks.length=unit(-0.5, "lines"),
# 	    axis.ticks.margin=unit(1.0, "lines")) +
# 	theme(plot.margin=unit(c(1.5,1,1.15,1), "lines"))
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=rel(1.5)), 
        axis.title.x=element_blank(),  
        axis.title.y=element_text(size=rel(2), face="bold"),
        axis.ticks.length=unit(-0.5, "lines"),
        axis.ticks.margin=unit(1.0, "lines")) +
  theme(plot.margin=unit(c(1.5,1,1,1), "lines"))
}
plot.npp <- ggplotGrob(plot.npp )
plot.npp $heights[[5]] <- unit(5, "null")
plot(plot.npp)
# --------

# --------
# 4.b.2. NPP deviation by data stream
# --------
dat.plot.dev <- dev.region[dev.region$Year >=1700,]

plot.dev <- {
ggplot() + 
	scale_x_continuous(expand=c(0,0), name="Year") +
	scale_y_continuous(expand=c(0,0), name="frac. NPP") +
	facet_grid(Y.type~., scales="free_y", space="free") +
	geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Model",], aes(x=Year, ymin=Y.rel.10.lo, ymax=Y.rel.10.hi), alpha=0.5) +
	geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Model",], aes(x=Year, y=Y.rel.10), size=2,
	          color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.tair.10.adj"    ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.CO2.10.adj"     ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.precipf.10.adj" ])), 
                        size=3) +

	geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP",],
	            aes(x=Year, ymin=Y.rel.10.lo, ymax=Y.rel.10.hi), alpha=0.5) +
	geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings"  & dat.plot.dev$Y.type=="NPP",], aes(x=Year, y=Y.rel.10), size=2,
	          color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.tair.10.adj"    ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.CO2.10.adj"     ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.precipf.10.adj" ])), 
                        size=3) +

	geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW",],
	            aes(x=Year, ymin=Y.rel.10.lo, ymax=Y.rel.10.hi), alpha=0.5) +
	geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings"  & dat.plot.dev$Y.type=="RW",], aes(x=Year, y=Y.rel.10), size=2,
	          color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.tair.10.adj"    ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.CO2.10.adj"     ]),
                        abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.precipf.10.adj" ])), 
                        size=3) +
	geom_hline(yintercept=1, linetype="dashed") +
  scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme(legend.title=element_text(size=rel(1), face="bold"),
	      legend.text=element_text(size=rel(1)),
	      # legend.position=c(0.2, 0.18),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(fill=NA, color="black", size=0.5), 
	      panel.background=element_blank(),
	      panel.margin.x=unit(0, "lines"),
	      panel.margin.y=unit(0, "lines"))  +
# 	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
# 		  axis.text.y=element_text(color="black", size=rel(1.5)), 
# 		  axis.title.x=element_text(size=rel(2), face="bold", vjust=-0.4),  
# 		  axis.title.y=element_text(size=rel(2), face="bold"),
# 		  axis.ticks.length=unit(-0.5, "lines"),
# 	    axis.ticks.margin=unit(1.0, "lines")) +
# 	theme(plot.margin=unit(c(1.5,1,1.15,1), "lines"))
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=rel(1.5)), 
        axis.title.x=element_blank(),  
        axis.title.y=element_text(size=rel(2), face="bold"),
        axis.ticks.length=unit(-0.5, "lines"),
        axis.ticks.margin=unit(1.0, "lines")) +
  theme(plot.margin=unit(c(0,1,1,1), "lines"))
}
plot.dev <- ggplotGrob(plot.dev )
plot.dev $heights[[5]] <- unit(0.5, "null")
plot(plot.dev)
# --------

# --------
# 4.b.3. Factor weights
# --------
{
fit.stack <- stack(dev.region[,c("fit.CO2.rel.10.adj","fit.tair.rel.10","fit.precipf.rel.10")])
names(fit.stack) <- c("fit.mean", "Effect")
fit.stack$Effect <- as.factor(ifelse(substr(fit.stack$Effect,5,7)=="CO2", "CO2", ifelse(substr(fit.stack$Effect,5,8)=="tair", "tair", "precipf")))
fit.stack[,c("data.type", "Y.type", "Year")] <- dev.region[,c("data.type", "Y.type", "Year")]
fit.stack[,"ci.lo"] <- stack(dev.region[,c("fit.CO2.rel.10.adj.lo","fit.tair.rel.10.lo","fit.precipf.rel.10.lo")])[,1]
fit.stack[,"ci.hi"] <- stack(dev.region[,c("fit.CO2.rel.10.adj.hi","fit.tair.rel.10.hi","fit.precipf.rel.10.hi")])[,1]

fit.stack[!is.na(fit.stack$ci.lo) & fit.stack$ci.lo<0.44, "ci.lo"] <- 0.44
fit.stack[!is.na(fit.stack$ci.lo) & fit.stack$ci.hi<0.44, "ci.hi"] <- 0.44
fit.stack[!is.na(fit.stack$ci.hi) & fit.stack$ci.hi>1.95, "ci.hi"] <- 1.95
fit.stack[!is.na(fit.stack$fit.mean) & (fit.stack$fit.mean<0.44 | fit.stack$fit.mean>1.95),"fit.mean"] <- NA
summary(fit.stack)


plot.wts <- {
ggplot(fit.stack[!(fit.stack$data.type=="Tree Rings" & fit.stack$Y.type=="NPP"),]) + 
  scale_x_continuous(limits=c(1700,2010), expand=c(0,0), name="Year") +
	scale_y_continuous(expand=c(0,0), name="frac. NPP") +
	facet_grid(Y.type~., scales="free_y", space="free") +
	geom_ribbon(aes(x=Year, ymin=ci.lo, ymax=ci.hi, fill=Effect), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.mean, color=Effect), size=2) +
	geom_hline(yintercept=1, linetype="dashed") +
  scale_color_manual(values=c("green3", "blue", "red2")) +
	scale_fill_manual(values=c("green3", "blue", "red2")) +
	guides(color=guide_legend(nrow=1), 
         fill =guide_legend(nrow=1))+
  theme(legend.title=element_text(size=rel(1), face="bold"),
	      legend.text=element_text(size=rel(1)),
	      legend.position=c(0.2, 0.18),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(fill=NA, color="black", size=0.5), 
	      panel.background=element_blank(),
	      panel.margin.x=unit(0, "lines"),
	      panel.margin.y=unit(0, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
		  axis.text.y=element_text(color="black", size=rel(1.5)), 
		  axis.title.x=element_text(size=rel(2), face="bold", vjust=-0.4),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(0,1,1.15,1), "lines"))
}
plot.wts <- ggplotGrob(plot.wts )
plot.wts $heights[[5]] <- unit(0.5, "null")
plot(plot.wts)
}
# --------

# --------
# Putting NPP & change through time in context
# --------
pdf(file.path(fig.dir, "Fig1_NPP_Dev_1700-2010.pdf"), height=11, width=8.5)
grid.arrange(plot.npp, plot.dev, plot.wts, ncol=1)
dev.off()

# --------
}
# ---------------------

# ---------------------
# 4.c. Some summary statistics on model agreement through time
# ---------------------
{
summary(dat.region)

# --------
# Model agreement in terms of relative NPP
# --------
{
# pre-1900
summary(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"])
mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"])

# 1980-2010
summary(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"])
mean(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"]); sd(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"])

t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "Y.rel.sd"])
}
# --------

# --------
# Comparing Factor weights and the shift in their agreement 
# --------
{
# tair
t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.tair.rel.lo", "fit.tair.rel.hi", "fit.tair.rel")],2, mean)

# precipf
t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.precipf.rel.lo", "fit.precipf.rel.hi", "fit.precipf.rel")],2, mean)

# CO2
t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
summary(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.CO2.rel.adj.lo", "fit.CO2.rel.adj.hi", "fit.CO2.rel.adj")],2, mean)
}
# --------
}
# ---------------------
# ----------------------------------------



# ----------------------------------------
# 5. Graphing & Summary Statistics of Raw NPP
# ----------------------------------------
{
summary(dat.ecosys)

# Getting some quick site stats for Table 1:
{
# ecosys2 <- aggregate(dat.ecosys[dat.ecosys$Model=="ed2" & dat.ecosys$Resolution=="t.001", c("tair", "precipf", "CO2")],
#                      by=dat.ecosys[dat.ecosys$Model=="ed2" & dat.ecosys$Extent=="850-2010" & dat.ecosys$Resolution=="t.001",c("Site", "Extent")], 
#                      FUN=mean)
# ecosys2[,c("tair.sd", "precipf.sd", "CO2.sd")] <- aggregate(dat.ecosys[dat.ecosys$Model=="ed2" & dat.ecosys$Resolution=="t.001", c("tair", "precipf", "CO2")],
#                                                            by=dat.ecosys[dat.ecosys$Model=="ed2" & dat.ecosys$Extent=="850-2010" & dat.ecosys$Resolution=="t.001",c("Site", "Extent")], 
#                                                            FUN=sd)[,c("tair", "precipf", "CO2")]
# ecosys2[,c("tair")] <- ecosys2[,c("tair")]-273.15
# ecosys2
# 
# sites.table <- data.frame(Site  = c( "PHA" ,  "PHO" ,  "PUN" ,  "PBL" ,  "PDL" ,    "PMB"   ), 
#                           Lon   = c(-72.18 , -68.73 , -89.53 , -94.58 , -95.17 ,   -82.83   ), 
#                           Lat   = c( 42.54 ,  45.25 ,  62.22 ,  46.28 ,  47.17 ,    43.61   ),
#                           Biome = c("Decid", "Mixed", "Mixed", "Mixed", "Mixed", "Evergreen"))
# sites.table <- merge(sites.table, ecosys2[,c("Site", "tair", "tair.sd", "precipf", "precipf.sd")], all.x=T, all.y=T)
# sites.table[,c("tair", "tair.sd")] <- round(sites.table[,c("tair", "tair.sd")], 1)
# sites.table[,c("precipf", "precipf.sd")] <- round(sites.table[,c("precipf", "precipf.sd")], 0)
# 
# write.csv(sites.table, file=file.path(out.dir, "Table1_SiteDescriptions.csv"), row.names=F)
}
# models.use <- unique(dat.ecosys[,"Model.Order"])
# colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])

# ---------------------
# Getting a site mean & CI for the tree ring NPP product NPP
# ---------------------
tr.npp <- dat.ecosys[dat.ecosys$Model %in% c("TreeRingNPP", "TreeRingRW") ,]
summary(tr.npp)


tr.npp.site <- aggregate(tr.npp[,c("Y", "Y.10")], by= tr.npp[,c("Model", "Model.Order", "Y.type", "Site", "Year")], FUN=mean, na.rm=T)
tr.npp.site[,paste0(c("Y", "Y.10"), ".lwr")] <- aggregate(tr.npp[,c("Y", "Y.10")], by= tr.npp[,c("Model", "Model.Order", "Y.type", "Site", "Year")], FUN=quantile, 0.025, na.rm=T)[,c("Y", "Y.10")]
tr.npp.site[,paste0(c("Y", "Y.10"), ".upr")] <- aggregate(tr.npp[,c("Y", "Y.10")], by= tr.npp[,c("Model", "Model.Order", "Y.type", "Site", "Year")], FUN=quantile, 0.975, na.rm=T)[,c("Y", "Y.10")]
summary(tr.npp.site)
# ---------------------

# ---------------------
# 5.a. Figures (site-level)
# ---------------------
{
# --------
# 5.a.1. All Sites 1800-2010 (Fig. 1)
# --------

ggFig1 <- {
	ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRW", "TreeRingBAI", "TreeRingNPP"),])  + facet_grid(Y.type~Site, scales="free_y", space="free") +
	geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.1, alpha=0.3) + 
	geom_line(aes(x=Year, y=Y.10, color=Model.Order), size=0.75, alpha=1) + 
	geom_ribbon(data=tr.npp.site[,], aes(x=Year, ymin=Y.10.lwr, ymax=Y.10.upr, fill=Model.Order), size=0.5, alpha=0.5) +
	geom_line(data=tr.npp.site[,], aes(x=Year, y=Y.10, color=Model.Order), size=1, alpha=0.8) +
	scale_x_continuous(limits=c(1800, 2010), expand=c(0,0), breaks=seq(min(dat.ecosys$Year), max(dat.ecosys$Year), by=100)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_fill_manual(values=c("black", "gray50")) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2), fill=F) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1.25)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black"), 
	      axis.text.y=element_text(color="black"), 
	      axis.title.x=element_text(face="bold", vjust=-0.5),  
	      axis.title.y=element_text(face="bold", vjust=1))
}

Fig1 <- ggplotGrob(ggFig1 )
Fig1$heights[[7]] <- unit(10, "null")
plot(Fig1)

pdf(file.path(fig.dir, "NPP_Raw_AllSites_1800-2010_Simple.pdf"), height=8.5, width=11)
grid.newpage()
grid.draw(Fig1)
dev.off()

# --------

# --------
# 5.a.2. Models Only, All Sites 0850-2010 (Supplemental Figure 1)
# --------
pdf(file.path(fig.dir, "SuppFig1_NPP_Raw_AllSites_0850-2010_Simple_Models.pdf"), height=8.5, width=11)
{
print(
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRW", "TreeRingBAI", "TreeRingNPP"),])  + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.1, alpha=0.3) + 
	geom_line(aes(x=Year, y=Y.10, color=Model.Order), size=0.75, alpha=1) + 
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0), breaks=seq(min(dat.ecosys$Year), max(dat.ecosys$Year), by=250)) +
	# scale_y_continuous(limits=c(0,25), expand=c(0,0)) +
	scale_fill_manual(values=c("black", "gray50")) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2), fill=F) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1.25)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black"), 
	      axis.text.y=element_text(color="black"), 
	      axis.title.x=element_text(face="bold", vjust=-0.5),  
	      axis.title.y=element_text(face="bold", vjust=1))
)
}
dev.off()
# --------

# --------
# 5.a.3. Just PHA 1900-2010
# --------
pdf(file.path(fig.dir, "NPP_Raw_PHA_1900-2010_Simple.pdf"), height=8.5, width=11)
{
print(
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRW", "TreeRingBAI", "TreeRingNPP") & dat.ecosys$Site=="PHA",])  +
	facet_grid(Y.type~., space="free", scales="free_y") +
	geom_line(aes(x=Year, y=Y, color=Model.Order), size=1, alpha=0.3) + 
	geom_line(aes(x=Year, y=Y.10, color=Model.Order), size=2, alpha=1) + 
	geom_ribbon(data=tr.npp.site[tr.npp.site$Site=="PHA",], aes(x=Year, ymin=Y.lwr, ymax=Y.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=tr.npp.site[tr.npp.site$Site=="PHA",], aes(x=Year, y=Y, color=Model.Order), size=1.5, alpha=1) +
	# geom_line(data=tr.npp.site[tr.npp.site$Site=="PHA",], aes(x=Year, y=Y.10, color=Model.Order), size=2, alpha=1) +
	scale_x_continuous(limits=c(1900, 2010), expand=c(0,0)) +
	# scale_y_continuous(limits=c(0,20), expand=c(0,0)) +
	scale_fill_manual(values=colors.use[(length(colors.use)-1):length(colors.use)]) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2), fill=F) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1.25)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black"), 
	      axis.text.y=element_text(color="black"), 
	      axis.title.x=element_text(face="bold", vjust=-0.5),  
	      axis.title.y=element_text(face="bold", vjust=1))
)
}
dev.off()
# --------

# --------
# 5.a.4. PHA, PHO (locations with NPP products, 1950-2010)
# --------
pdf(file.path(fig.dir, "NPP_Raw_PHA_PHO_1900-2010_Simple.pdf"), height=8.5, width=11)
{
print(
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRW", "TreeRingBAI", "TreeRingNPP") & dat.ecosys$Site %in% c("PHA", "PHO"),])  + facet_grid(Y.type~Site, space="free", scales="free_y") +
	geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.5, alpha=0.3) + 
	geom_line(aes(x=Year, y=Y.10, color=Model.Order), size=1.5, alpha=0.8) + 
	geom_ribbon(data=tr.npp.site[tr.npp.site$Site %in% c("PHA", "PHO"),], aes(x=Year, ymin=Y.lwr, ymax=Y.upr, fill=Model.Order), size=0.5, alpha=0.5) +
	geom_line(data=tr.npp.site[tr.npp.site$Site %in% c("PHA", "PHO"),], aes(x=Year, y=Y, color=Model.Order), size=1, alpha=1) +
	scale_x_continuous(limits=c(1900, 2010), expand=c(0,0)) +
	# scale_y_continuous(limits=c(0,25), expand=c(0,0)) +
	scale_fill_manual(values=c("black", "gray50")) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2), fill=F) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1.25)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black"), 
	      axis.text.y=element_text(color="black"), 
	      axis.title.x=element_text(face="bold", vjust=-0.5),  
	      axis.title.y=element_text(face="bold", vjust=1))
)
}
dev.off()
# --------
}
# ---------------------

# ---------------------
# 5.b. Summary Statistics
# ---------------------
# --------
# 5.b.1 Description of Model variability across space & Time
# --------
{
# Aggregate to the site-level by model
models.sites                         <- aggregate(dat.ecosys[,c("Y", "Y.rel")], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                                  FUN=mean)
models.sites[,c("Y.sd", "Y.rel.sd")] <- aggregate(dat.ecosys[,c("Y", "Y.rel")], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                                  FUN=sd)[,c("Y", "Y.rel")]
summary(models.sites)

# Aggregate to the model level by time
models.time                         <- aggregate(dat.ecosys[,c("Y", "Y.rel")], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                  FUN=mean)
models.time[,c("Y.sd", "Y.rel.sd")] <- aggregate(dat.ecosys[,c("Y", "Y.rel")], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                  FUN=sd)[,c("Y", "Y.rel")]
summary(models.time)



# Compute some sumary statistics
models.stats                         <- aggregate(models.sites[,c("Y", "Y.rel")], 
                                                  by=models.sites[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                  FUN=mean)
models.stats[,c("Y.sd", "Y.rel.sd")] <- aggregate(models.sites[,c("Y", "Y.rel")], 
                                                  by=models.sites[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                  FUN=sd)[,c("Y", "Y.rel")]
models.stats

models.stats.time                         <- aggregate(models.time[,c("Y", "Y.rel")], 
                                                       by=models.time[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                       FUN=mean)
models.stats.time[,c("Y.sd", "Y.rel.sd")] <- aggregate(models.time[,c("Y", "Y.rel")], 
                                                       by=models.time[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                       FUN=sd)[,c("Y", "Y.rel")]
models.stats.time
}
# --------

# --------
# 5.b.2. Description of Site variability
# --------
{
site.stats      <- aggregate(models.sites[,c("Y", "Y.rel")], 
                             by=models.sites[,c("Site", "Y.type", "data.type")], 
                             FUN=mean)
site.stats[,c("Y.sd", "Y.rel.sd")]   <- aggregate(models.sites[,c("Y", "Y.rel")], 
                                                     by=models.sites[,c("Site", "Y.type", "data.type")], 
                                                     FUN=sd)[,c("Y", "Y.rel")]
site.stats
}
# --------


# ---------------------

}
# ----------------------------------------ˀ


# ----------------------------------------
# 6. Graphing & Analyzing Sensitivity
# ----------------------------------------
{
summary(ci.terms)
# summary(sim.terms[,1:10])
# -----------------------
# 6.a. Graphing
# -----------------------
{
# Trying out the basic plot to compare model responses to drivers
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))

colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-0.75),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.75),"lwr.rel"] <- -0.75 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.75),"upr.rel"] <- -0.75 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.0),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(1.0),"lwr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$upr.rel>(1.0),"upr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15

ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
summary(ci.terms.graph)


# Plot the relativized
pdf(file.path(fig.dir, "Sensitivity_Models_Rel_Baseline.pdf"), height=8.5, width=11)
{
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("tair", "precipf", "CO2"),]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme_bw()
)
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("Biomass"),]) + facet_wrap(~Model, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme_bw()
)
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("Time"),]) + facet_wrap(~Model, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
	theme_bw()
)
}
dev.off()

} # End graphing section
# -----------------------


# -----------------------
# 6.b. Quantitative Analysis
# -----------------------
source("R/0_Calculate_GAMM_Posteriors.R")
library(mgcv)

# --------
# 6.b.0. Adding some model-level stats to compare the relative sensitivities
# --------
{
# Condensing model variability across space and time to get general model characteristics
{
vars.agg <- c("Y", "Y.rel", "Biomass")
mod.site                           <- aggregate(dat.ecosys[,vars.agg], 
                                               by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                               FUN=mean)
mod.site[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys[,vars.agg], 
                                               by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                               FUN=sd)[,vars.agg]
mod.site$Biomass.sd.per <- mod.site$Biomass.sd/mod.site$Biomass
summary(mod.site)

mod.agg                           <- aggregate(mod.site[,c(vars.agg, paste0(vars.agg, ".sd"), "Biomass.sd.per")], 
                                               by=mod.site[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                               FUN=mean)
mod.agg



# Finding the change in key variables in the modern era
for(v in vars.agg){
  mod.agg[,paste0("dModern.", v)] <- NA  
}

for(m in unique(mod.agg$Model)){
  mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1990 & dat.ecosys$Year<=2010, vars.agg], na.rm=T) -
    colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1830 & dat.ecosys$Year<=1850, vars.agg], na.rm=T)
}
mod.agg

# Add in the vegetation scheme
mod.agg$veg.scheme <- as.factor(ifelse(mod.agg$Model %in% c("clm.bgc", "clm.cn", "sibcasa", "jules.stat"), "Static", "Dynamic"))

# Finding mean fire return
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

ci.terms <- merge(ci.terms, mod.agg, all.x=T, all.y=F)
summary(ci.terms)

# Making the data frame so we can graph & compare the curves quantitatively
{
factors.analy <- c("Y", "Y.rel", "Y.sd", "Y.rel.sd", "dModern.Y", "dModern.Y.rel", "dModern.Biomass", "Biomass", "Biomass.sd", "Biomass.sd.per")
  
df.co2 <- aggregate(ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",c(factors.analy)], 
                    by=ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                    FUN=mean, na.rm=T)
summary(df.co2)

df.tair <- aggregate(ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",c(factors.analy)], 
                     by=ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                     FUN=mean, na.rm=T)
summary(df.tair)

df.precipf <- aggregate(ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",c(factors.analy)], 
                        by=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                        FUN=mean, na.rm=T)
summary(df.precipf)
}

# setting up some null models for the full sensitivity curves
{
co2.null     <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
tair.null    <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
precipf.null <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])

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


# Some quick plots for sanity check
{
ggplot(data=ci.terms[ci.terms$Effect %in% c("tair", "precipf", "CO2") & ci.terms$data.type=="Model",])+
  facet_wrap(~Effect, scales="free_x") +
#   geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
#   geom_line(aes(x=x, y=mean.rel, color=Model), size=2) +
  geom_ribbon(data=null.post, aes(x=x, ymin=lwr, ymax=upr), fill="black", alpha=0.5) +
  geom_line(data=null.post, aes(x=x, y=mean, color=Effect), size=2, color="black") +
  scale_fill_manual(values=c(colors.use)) +
  scale_color_manual(values=colors.use) +
  theme_bw()

ggplot(data=ci.terms[ci.terms$Effect %in% c("tair", "precipf", "CO2") & ci.terms$data.type=="Model",])+
  facet_wrap(~Effect, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
    geom_line(aes(x=x, y=mean.rel, color=Model), size=2) +
  geom_ribbon(data=null.post, aes(x=x, ymin=lwr, ymax=upr), fill="black", alpha=0.5) +
  geom_line(data=null.post, aes(x=x, y=mean, color=Effect), size=3, color="black") +
#   geom_point(aes(x=x, y=mean.rel, color=Model)) +
#   stat_smooth(aes(x=x, y=mean.rel)) +
  scale_fill_manual(values=c(colors.use)) +
  scale_color_manual(values=colors.use) +
  theme_bw()
}
}

# Condensing the full sensitivity curves to the values at the 25, 50, and 75% for 
# analysis of continuous characteristics of models (i.e. NPP, modern change, etc)
{
summary(ci.terms)
summary(dat.ecosys)

co2.driver     <- dat.ecosys[dat.ecosys$Model=="ed2","CO2"    ]
tair.driver    <- dat.ecosys[dat.ecosys$Model=="ed2","tair"   ]
precipf.driver <- dat.ecosys[dat.ecosys$Model=="ed2","precipf"]

for(e in c("tair", "precipf", "CO2")){
  # Get the distirbution of met drivers for each driver & specify to what precision we want to round
  effect.driver <- dat.ecosys[dat.ecosys$Model=="ed2", e]
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
                          by=ci.terms[,c("Effect", "data.type", "Y.type", "Model", "Model.Order", "veg.scheme", "fire.scheme", "Quantile")],
                          FUN=mean, na.rm=T)
summary(ci.terms.agg)

# # Just a few quick graphs to make sure this looks reasonable
# ggplot(data=ci.terms.agg) +
#   facet_wrap(~Effect) +
#   geom_point(aes(x=Y, y=mean.rel, color=Model), size=8) +
#   scale_color_manual(values=colors.use) +
#   stat_smooth(aes(x=Y, y=mean.rel), method="lm") +
#   ggtitle("Aggregated CIs") +
#   theme_bw()
# 
# ggplot(data=ci.terms.agg) +
#   facet_wrap(~Effect, scales="free_x") +
#   geom_point(aes(x=x, y=mean.rel, color=Model), size=8) +
#   scale_color_manual(values=colors.use) +
# #   stat_smooth(aes(x=Y, y=mean.rel), method="lm") +
#   ggtitle("Aggregated CIs") +
#   theme_bw()

}

}
# --------


# --------
# 6.b.1 Analyzing Static vs. Dynamic Veg
# --------
{
# Doing the comparisons with the full sensitivity curve
co2.veg     <- gam(mean.rel ~ s(x, by=veg.scheme) + veg.scheme, data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
tair.veg    <- gam(mean.rel ~ s(x, by=veg.scheme) + veg.scheme, data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
precipf.veg <- gam(mean.rel ~ s(x, by=veg.scheme) + veg.scheme, data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
tair.veg2   <- gam(mean.rel ~ s(x, by=veg.scheme) + veg.scheme, data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model" & !ci.terms$Model=="jules.stat",])

summary(co2.veg)
summary(tair.veg)
summary(precipf.veg)

co2.veg.lm     <- lm(mean.rel ~ veg.scheme, data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
tair.veg.lm    <- lm(mean.rel ~ veg.scheme, data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
precipf.veg.lm <- lm(mean.rel ~ veg.scheme, data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
summary(co2.veg.lm)
summary(tair.veg.lm)
summary(precipf.veg.lm)

# Plot the different sensitivity curves by characteristic
co2.veg.post <- post.distns(model.gam=co2.veg, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
co2.veg.post$veg.scheme <- df.co2$veg.scheme
co2.veg.post$Effect <- as.factor("CO2")
summary(co2.veg.post)

tair.veg.post <- post.distns(model.gam=tair.veg, model.name="Tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
tair.veg.post$veg.scheme <- df.tair$veg.scheme
tair.veg.post$Effect <- as.factor("tair")
summary(tair.veg.post)

tair.veg.post2 <- post.distns(model.gam=tair.veg2, model.name="Tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
tair.veg.post2$veg.scheme <- df.tair$veg.scheme
tair.veg.post2$veg.scheme <- as.factor(paste0(tair.veg.post2$veg.scheme, "-No JULES-STATIC"))
tair.veg.post2$Effect <- as.factor("tair")
summary(tair.veg.post2)

precipf.veg.post <- post.distns(model.gam=precipf.veg, model.name="Precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
precipf.veg.post$veg.scheme <- df.precipf$veg.scheme
precipf.veg.post$Effect <- as.factor("precipf")
summary(precipf.veg.post)

veg.post <- rbind(tair.veg.post, precipf.veg.post, co2.veg.post, tair.veg.post2)
summary(veg.post)

# ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("tair", "precipf", "CO2") & ci.terms.graph$data.type=="Model",])+
#   facet_wrap(~Effect, scales="free_x") +
#   geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
#   geom_line(aes(x=x, y=mean.rel, color=Model), size=2) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_fill_manual(values=colors.use) +
#   scale_color_manual(values=colors.use) +
#   theme_bw()

pdf(file.path(fig.dir, "Sensitivity_VegScheme.pdf"))
print(
ggplot(data=veg.post[!veg.post$veg.scheme=="Dynamic-No JULES-STATIC",]) +
  facet_wrap(~Effect, scales="free_x") +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=veg.scheme), alpha=0.5) +
  geom_line(aes(x=x, y=mean, color=veg.scheme, linetype=veg.scheme), size=2) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("blue3", "red3", "green3")) +
  scale_color_manual(values=c("blue3", "red3", "green3")) +
  theme_bw()
)
dev.off()
}

# Comparing the relative change in sensitivity from dynamic vs. static veg
summary(veg.post)
veg.post2 <- veg.post[veg.post$veg.scheme %in% c("Dynamic", "Static"),]
summary(veg.post2)
veg.post2 <- aggregate(veg.post2[,c("mean", "lwr", "upr")], by=veg.post2[,c("Effect", "veg.scheme", "x")], FUN=mean)
summary(veg.post2)

for(e in unique(veg.post2$Effect)){
  dynamic <- veg.post2[veg.post2$Effect==e & veg.post2$veg.scheme=="Dynamic","mean"]
  static  <- veg.post2[veg.post2$Effect==e & veg.post2$veg.scheme=="Static" ,"mean"]
  
  dif.temp <- data.frame(Effect=e, x=veg.post2[veg.post2$Effect==e & veg.post2$veg.scheme=="Dynamic","x"], diffs=static-dynamic)
  if(e == unique(veg.post2$Effect)[1]){
    veg.diffs <- dif.temp  
  } else {
    veg.diffs <- rbind(veg.diffs, dif.temp)
  }
}
summary(veg.diffs)
mean(abs(veg.diffs[veg.diffs$Effect=="tair","diffs"])); sd(abs(veg.diffs[veg.diffs$Effect=="tair","diffs"]))
mean(abs(veg.diffs[veg.diffs$Effect=="precipf","diffs"])); sd(abs(veg.diffs[veg.diffs$Effect=="precipf","diffs"]))
mean(abs(veg.diffs[veg.diffs$Effect=="CO2","diffs"])); sd(abs(veg.diffs[veg.diffs$Effect=="CO2","diffs"]))

# --------

# --------
# 6.b.2 By Fire presence/absences
# --------
{
# Doing the comparisons with the full sensitivity curve
co2.fire     <- gam(mean.rel ~ s(x, by=fire.scheme) + fire.scheme, data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
tair.fire    <- gam(mean.rel ~ s(x, by=fire.scheme) + fire.scheme, data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
precipf.fire <- gam(mean.rel ~ s(x, by=fire.scheme) + fire.scheme, data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])

# Plot the different sensitivity curves by characteristic
co2.fire.post <- post.distns(model.gam=co2.fire, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
co2.fire.post$fire.scheme <- df.co2$fire.scheme
co2.fire.post$Effect <- as.factor("CO2")
summary(co2.fire.post)

co2.fire.post <- aggregate(co2.fire.post[,c("mean", "lwr", "upr")], by=co2.fire.post[,c("fire.scheme", "Effect", "x")], FUN=mean)

tair.fire.post <- post.distns(model.gam=tair.fire, model.name="Tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
tair.fire.post$fire.scheme <- df.tair$fire.scheme
tair.fire.post$Effect <- as.factor("tair")
summary(tair.fire.post)

precipf.fire.post <- post.distns(model.gam=precipf.fire, model.name="Precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
precipf.fire.post$fire.scheme <- df.precipf$fire.scheme
precipf.fire.post$Effect <- as.factor("precipf")
summary(precipf.fire.post)

fire.post <- rbind(tair.fire.post, precipf.fire.post, co2.fire.post)
summary(fire.post)

ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("tair", "precipf", "CO2") & ci.terms.graph$data.type=="Model",])+
  facet_wrap(~Effect, scales="free_x") +
  geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
  geom_line(aes(x=x, y=mean.rel, color=Model), size=2) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  theme_bw()

pdf(file.path(fig.dir, "Sensitivity_FireScheme.pdf"))
print(
ggplot(data=fire.post) +
  facet_wrap(~Effect, scales="free_x") +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=fire.scheme), alpha=0.5) +
  geom_line(aes(x=x, y=mean, color=fire.scheme, linetype=fire.scheme), size=2) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("blue3", "red3")) +
  scale_color_manual(values=c("blue3", "red3")) +
  theme_bw()
)
dev.off()
}

# Comparing the relative change in sensitivity from dynamic vs. static veg
summary(fire.post)
fire.post2 <- fire.post[,]
summary(fire.post2)
fire.post2 <- aggregate(fire.post[,c("mean", "lwr", "upr")], by=fire.post[,c("Effect", "fire.scheme", "x")], FUN=mean)
summary(fire.post2)

fire.post2 <- fire.post2[!(fire.post2$Effect=="CO2" & fire.post2$fire.scheme=="No" & !(fire.post2$x %in% fire.post2[fire.post2$Effect=="CO2" & fire.post2$fire.scheme=="Yes","x"])),]
summary(fire.post2)

for(e in unique(fire.post2$Effect)){
  dynamic <- fire.post2[fire.post2$Effect==e & fire.post2$fire.scheme=="Yes","mean"]
  static  <- fire.post2[fire.post2$Effect==e & fire.post2$fire.scheme=="No" ,"mean"]
  
  dif.temp <- data.frame(Effect=e, x=fire.post2[fire.post2$Effect==e & fire.post2$fire.scheme=="No","x"], diffs=static-dynamic)
  if(e == unique(fire.post2$Effect)[1]){
    fire.diffs <- dif.temp  
  } else {
    fire.diffs <- rbind(fire.diffs, dif.temp)
  }
}
summary(fire.diffs)
summary((fire.diffs[fire.diffs$Effect=="tair","diffs"]))
summary((fire.diffs[fire.diffs$Effect=="precipf","diffs"]))
summary((fire.diffs[fire.diffs$Effect=="CO2","diffs"]))

mean(abs(fire.diffs[fire.diffs$Effect=="tair","diffs"])); sd(abs(fire.diffs[fire.diffs$Effect=="tair","diffs"]))
mean(abs(fire.diffs[fire.diffs$Effect=="precipf","diffs"])); sd(abs(fire.diffs[fire.diffs$Effect=="precipf","diffs"]))
mean(abs(fire.diffs[fire.diffs$Effect=="CO2","diffs"])); sd(abs(fire.diffs[fire.diffs$Effect=="CO2","diffs"]))

# --------

# --------
# 6.b.3 By Fire intensity (only where present)
# --------
{
co2.fire.lm     <- lm(mean.rel ~ Fire, data=ci.terms.agg[ci.terms.agg$Effect=="CO2"     & 
                                                           ci.terms.agg$data.type=="Model" & 
                                                           ci.terms.agg$fire.scheme=="Yes",])
tair.fire.lm    <- lm(mean.rel ~ Fire, data=ci.terms.agg[ci.terms.agg$Effect=="tair"    & 
                                                           ci.terms.agg$data.type=="Model" & 
                                                           ci.terms.agg$fire.scheme=="Yes",])
precipf.fire.lm <- lm(mean.rel ~ Fire, data=ci.terms.agg[ci.terms.agg$Effect=="precipf" & 
                                                           ci.terms.agg$data.type=="Model" & 
                                                           ci.terms.agg$fire.scheme=="Yes",])
summary(co2.fire.lm    )
summary(tair.fire.lm   )
summary(precipf.fire.lm)

models.fire <- data.frame(Model=unique(ci.terms.agg[ci.terms.agg$fire.scheme=="Yes","Model"]), 
                          Model.Order=unique(ci.terms.agg[ci.terms.agg$fire.scheme=="Yes","Model.Order"]))
colors.fire <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.fire$Model.Order, "color"])))

pdf(file.path(fig.dir, "Sensitivity_FireIntensity.pdf"))
print(
ggplot(data=ci.terms.agg[ci.terms.agg$fire.scheme=="Yes",]) + 
  facet_wrap(~ Effect, scales="free_x") +
  geom_point(aes(x=Fire, y=mean.rel, color=Model), size=5) +
  stat_smooth(aes(x=Fire, y=mean.rel), method="lm", size=2) +
  scale_color_manual(values=colors.fire) +
  theme_bw()
)
dev.off()
}
# --------

# --------
# 6.b.4 By NPP
# --------
{
co2.npp.lm     <- lm(mean.rel ~ Y, data=ci.terms.agg[ci.terms.agg$Effect=="CO2"     & 
                                                           ci.terms.agg$data.type=="Model",])
tair.npp.lm    <- lm(mean.rel ~ Y, data=ci.terms.agg[ci.terms.agg$Effect=="tair"    & 
                                                           ci.terms.agg$data.type=="Model",])
precipf.npp.lm <- lm(mean.rel ~ Y, data=ci.terms.agg[ci.terms.agg$Effect=="precipf" & 
                                                           ci.terms.agg$data.type=="Model",])
summary(co2.npp.lm    )
summary(tair.npp.lm   )
summary(precipf.npp.lm)


pdf(file.path(fig.dir, "Sensitivity_NPP.pdf"))
print(
  ggplot(data=ci.terms.agg[,]) + 
    facet_wrap(~ Effect, scales="free_x") +
    geom_point(aes(x=Y, y=mean.rel, color=Model), size=5) +
    stat_smooth(aes(x=Y, y=mean.rel), method="lm", size=2) +
    scale_color_manual(values=colors.use) +
    theme_bw()
)
dev.off()
}
# --------

# --------
# 6.b.5 By Biomass variability
# --------
{
  co2.Biomass.sd.lm     <- lm(mean.rel ~ Biomass.sd.per, data=ci.terms[!is.na(ci.terms$Quantile) & ci.terms$Effect=="CO2"     & 
                                                             ci.terms$data.type=="Model",])
  tair.Biomass.sd.lm    <- lm(mean.rel ~ Biomass.sd.per, data=ci.terms[!is.na(ci.terms$Quantile) & ci.terms$Effect=="tair"    & 
                                                             ci.terms$data.type=="Model",])
  precipf.Biomass.sd.lm <- lm(mean.rel ~ Biomass.sd.per, data=ci.terms[!is.na(ci.terms$Quantile) & ci.terms$Effect=="precipf" & 
                                                             ci.terms$data.type=="Model",])
  summary(co2.Biomass.sd.lm    )
  summary(tair.Biomass.sd.lm   )
  summary(precipf.Biomass.sd.lm)
  
  
  pdf(file.path(fig.dir, "Sensitivity_Biomass_sd_per.pdf"))
  print(
    ggplot(data=ci.terms[!is.na(ci.terms$Quantile) & ci.terms$data.type=="Model",]) + 
      facet_wrap(~ Effect, scales="free_x") +
      geom_point(aes(x=Biomass.sd.per, y=mean.rel, color=Model), size=5) +
      stat_smooth(aes(x=Biomass.sd.per, y=mean.rel), method="lm", size=2) +
      scale_color_manual(values=colors.use) +
      theme_bw()
  )
  dev.off()
}
# --------


}
# -----------------------
# ----------------------------------------






