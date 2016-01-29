# ----------------------------------------
# Objective: Compare model/tree ring sensitivites by forest type/PFT
# Christy Rollinson, crollinson@gmail.com
# Date Created: 16 November 2015
# ----------------------------------------
#
# -------------------------
# Hypotheses: 
# -------------------------
# 1. Models with the climate sensitivity closes to the data at the 
#    regional scale have biogeography right
# 2. Comparing climate sensitivity of specific forest types will 
#    result in more agreement among models & data than comparing at 
#    the region scale.
# -------------------------
#
# -------------------------
# Workflow
# -------------------------
# 1. Set Directories
# 2. Load data files & function scripts
# 3. Standardize driver responses to the mean model NPP to aid comparison (loop)
#    a. Find the NPP to relativize each set off of
#    b. Relativizing everything in dat.ecosys to make it comparable to tree rings
#    c. Finding the percent change in NPP relative to the mean for that particular scale
# 4. Graphing & Analyzing Sensitivity
#    a. Graphing
#    b. Quantitative Analysis
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
in.base <- "Data/gamms/Sensitivity_PFT"
out.dir <- "Data/analyses/analysis_biogeography"
fig.dir <- "Figures/analyses/analysis_biogeography"

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
load(file.path(in.base, "gamm_PFT.Rdata"))

dat.ecosys <- cbind(mod.out$data[, ], mod.out$ci.response[,c("mean", "lwr", "upr")])
ci.terms   <- mod.out$ci.terms[,]
wt.terms   <- mod.out$weights[,]
sim.terms  <- mod.out$sim.terms
sim.terms$Effect <- as.factor(sim.terms$Effect)

# Force the terms ci & sims "x" variable to numeric
ci.terms $x <- as.numeric(paste(ci.terms $x))
sim.terms$x <- as.numeric(paste(sim.terms$x))

# Grouping the kind and source of the data
dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRW", "RW", "NPP"))
ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRW", "RW", "NPP"))
wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRW", "RW", "NPP"))
sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRW", "RW", "NPP"))

dat.ecosys$data.type <- as.factor(ifelse(substr(dat.ecosys$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
ci.terms  $data.type <- as.factor(ifelse(substr(ci.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree 
Rings", "Model"))
sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "Tree 
Rings", "Model"))

summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)
summary(sim.terms[,1:10])
# ----------------------------------------


# ----------------------------------------
# 3. Standardize driver responses to the mean model NPP to facilitate comparisons
#	 -- Two things for standardization and graphing:
# 		1. Relative by model-mean NPP
#       2. Decadal smoothing to help show generalized patterns
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
summary(dat.ecosys)

# Make sure all data sets are ordered by year, then treeID, then plotID, then Model
sort.order <- c("Model", "PlotID", "TreeID", "Year")
dat.ecosys <- dat.ecosys[order(dat.ecosys$Model, dat.ecosys$PlotID, dat.ecosys$TreeID, dat.ecosys$Year),]
wt.terms <- wt.terms[order(wt.terms$Model, wt.terms$PlotID, wt.terms$TreeID, wt.terms$Year),]

# Double Check to make sure things are sorted by year so rollapply works
dat.ecosys[which(dat.ecosys$Model=="TreeRingRW")[1:20],]
wt.terms  [which(wt.terms  $Model=="TreeRingRW")[1:20],]

{
for(m in unique(ci.terms$Model)){

		# -----------------------
		# 3.a. Find the NPP to relativize each set off of
		# Using mean model NPP across sites since the GAMM response curves are for 
		#    the whole model & not site-specific are parameterized
		# -----------------------
		# Find the start year for the extent
		# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

		npp <- mean(dat.ecosys[dat.ecosys$Model==m, "Y"], na.rm=T)			
		# -----------------------
		
		# -----------------------
		# 3.b Relativizing everything in dat.ecosys to make it comparable to tree rings
		# -----------------------
		{		
		# Which factors to relativize
		y.rel <- c("Y", "fit.gam", "mean", "lwr", "upr")

		# for some reason, I can create multiple new columns at once
		# Solution: use a loop to create blank columns and then fill them
		for(y in y.rel){
			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel"       )] <- NA	
			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".10"        )] <- NA	
			dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel", ".10")] <- NA	
		}
		dat.ecosys[dat.ecosys$Model==m,paste0(y.rel, ".rel")] <- dat.ecosys[dat.ecosys$Model==m, y.rel]/npp
		
		# Getting 10-year running means to make clearer figures
		for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
			# Note: If we're working with tree ring data, we need to go by plot for NPP products 
			#       & by Tree for individual-level tree rings products
			if(m=="TreeRingNPP"){
				for(p in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "PlotID"])){

					# Raw NPP (to add dark line over faded annual wiggles)
					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# Relativized NPP (to have generalized patterns for figures)
					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$PlotID==p, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				}
			} else if(m=="TreeRingRW") {
				for(t in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "TreeID"])){
					# If we have too few data points, we need to skip that tree 
					if(length(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel[1]]) < 10) next

					# Raw NPP (to add dark line over faded annual wiggles)
					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# Relativized NPP (to have generalized patterns for figures)
					dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$TreeID==t, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				}
			} else {
				# Raw NPP (to add dark line over faded annual wiggles)
				dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

				# Relativized NPP (to have generalized patterns for figures)
				dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
			}
		}
		}		
		# -----------------------

		
		# -----------------------
		# 3.c. Finding the percent change in NPP relative to the mean for that particular scale
		# -----------------------
		{
		y.rel <- c("mean", "lwr", "upr")
		for(y in y.rel){
			ci.terms[ci.terms$Model==m,paste0(y, ".rel"       )] <- NA	
		}		

		ci.terms[ci.terms$Model==m,paste0(y.rel,".rel")] <- ci.terms[ci.terms $Model==m, y.rel]/npp
		
		# Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
		cols.sim <- which(substr(names(sim.terms),1,1)=="X")
		sim.terms[sim.terms$Model==m,cols.sim] <- sim.terms[sim.terms$Model==m,cols.sim]/npp
		}
		# -----------------------

		# -----------------------
		# 3.d. Relativizing the factor fits through times and weights as well
		# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
		# -----------------------
		{
		y.rel <- c("fit.tair", "fit.precipf", "fit.CO2")
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
		for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
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


summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])
# ----------------------------------------


# ----------------------------------------
# 4. Graphing & Analyzing Sensitivity
# ----------------------------------------

# -----------------------
# 4.a. Graphing
# -----------------------
{
# Trying out the basic plot to compare model responses to drivers
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-2.0),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-2.0),"lwr.rel"] <- -2.0 
ci.terms.graph[ci.terms.graph$upr.rel<(-2.0),"upr.rel"] <- -2.0 
ci.terms.graph[which(ci.terms.graph$mean.rel>2.00),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(2.0),"lwr.rel"] <- 2.0 
ci.terms.graph[ci.terms.graph$upr.rel>(2.0),"upr.rel"] <- 2.0 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)

ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
summary(ci.terms.graph)

ci.terms <- merge(ci.terms, models.df, all.x=T, all.y=F)
summary(ci.terms)

# Plot the relativized
pdf(file.path(fig.dir, "Fig3_Sensitivity_Models_Rel_PFT.pdf"), height=8.5, width=11)
{
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("tair", "precipf", "CO2"),]) + facet_grid(PFT~Effect, scales="free_x") +
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

pdf(file.path(fig.dir, "Sensitivity_Models_Rel_PFT_breakdown.pdf"), height=8.5, width=11)
for(e in c("tair", "precipf", "CO2")){
# print(summary(ci.terms[ci.terms$Effect==e,]))
print(
ggplot(data= ci.terms.graph[ci.terms.graph$Effect==e,]) + facet_wrap(~PFT, scales="fixed") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
	scale_x_continuous(expand=c(0,0), name=e) +
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
# 4.b. Analysis: Looking at the change in Effect Sensitivity from added site effect
# Note: Looking at absolute value of the effect to really be able to get at if we're getting MORE (negative) or LESS (positive) sensitive
# -----------------------
# -----------------------
# ----------------------

# ----------------------
# Analyzing the differences
# ----------------------
# ----------------------
# ----------------------------------------
