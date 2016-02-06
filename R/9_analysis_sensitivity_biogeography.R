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
{
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
wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "TreeRings", "Model"))

summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)
summary(sim.terms[,1:10])
}
# ----------------------------------------


# ----------------------------------------
# 3. Standardize driver responses to the mean model NPP to facilitate comparisons
#	 -- Two things for standardization and graphing:
# 		1. Relative by model-mean NPP
#       2. Decadal smoothing to help show generalized patterns
# ----------------------------------------
{
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

save(dat.ecosys, ci.terms, wt.terms, sim.terms, file=file.path(out.dir, "post-process_composition.RData"))
}
# ----------------------------------------


# ----------------------------------------
# 4. Graphing & Analyzing Sensitivity
# ----------------------------------------
load(file.path(out.dir, "post-process_composition.RData"))
summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])


# -----------------------
# 4.a. Graphing
# -----------------------
{
# Trying out the basic plot to compare model responses to drivers
models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-1.0),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-1.0),"lwr.rel"] <- -1.0 
ci.terms.graph[ci.terms.graph$upr.rel<(-1.0),"upr.rel"] <- -1.0 
ci.terms.graph[which(ci.terms.graph$mean.rel>2.00),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(2.0),"lwr.rel"] <- 2.0 
ci.terms.graph[ci.terms.graph$upr.rel>(2.0),"upr.rel"] <- 2.0 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)

ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
summary(ci.terms.graph)

ci.terms <- merge(ci.terms, models.df, all.x=T, all.y=F)
summary(ci.terms)

ci.terms.graph$PFT2 <- recode(ci.terms.graph$PFT, "'Deciduous'='1'; 'Mixed'='2'; 'Evergreen'='3'; 'Grass/Savanna'='4'")
levels(ci.terms.graph$PFT2) <- c("Deciduous", "Mixed", "Evergreen", "Grass/Savanna")
  
# Plot the relativized
pdf(file.path(fig.dir, "Fig3_Sensitivity_Models_Rel_PFT.pdf"), height=8.5, width=11)
{
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect %in% c("tair", "precipf", "CO2"),]) + facet_grid(PFT2~Effect, scales="free_x") +
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
ggplot(data= ci.terms.graph[ci.terms.graph$Effect==e,]) + facet_wrap(~PFT2, scales="fixed") +
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
# --------
# 4.b.0. Getting the first derivative (first diference) of each line so we can take the mean slope
# --------
{
# First make sure the effects are sorted by x to make this easier
ci.terms <- ci.terms[order(ci.terms$Model, ci.terms$PFT, ci.terms$Effect, ci.terms$x),]
sim.terms <- sim.terms[order(sim.terms$Model, sim.terms$PFT, sim.terms$Effect, sim.terms$x),]
ci.terms[1:20,1:15]
dim(ci.terms)
summary(ci.terms)

# Making a new dataframe dedicated to the derivatives
cols.sims <- which(substr(names(sim.terms),1,1)=="X")
sim.deriv <- sim.terms[,]
sim.deriv[,cols.sims] <- NA

for(e in unique(ci.terms$Effect)){
  for(m in unique(ci.terms$Model)){
    for(p in unique(ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "PFT"])){
      x.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$PFT==p, "x"], lag=1), NA)
      y.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$PFT==p, "mean.rel"], lag=1), NA)
      ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$PFT==p, "deriv"] <- y.dif/x.dif

      # For the full simiulation for robust analysis
      y.dif2 <- rbind(apply(sim.terms[sim.terms$Model==m & sim.terms$Effect==e & sim.terms$PFT==p, cols.sims], 2, FUN=diff), NA)
      x.dif <- c(diff(sim.terms[sim.terms$Model==m & sim.terms$Effect==e  & sim.terms$PFT==p, "x"], lag=1), NA)
      sim.deriv[sim.deriv$Model==m & sim.deriv$Effect==e & sim.terms$PFT==p, cols.sims] <- apply(y.dif2, 2, FUN=function(y){y/x.dif})
    }
  } 
}
summary(ci.terms)
summary(sim.deriv[,1:10])

# Stacking and aggregating the simulations
deriv.stack <- stack(sim.deriv[,cols.sims])
names(deriv.stack)  <- c("deriv", "sim")
deriv.stack[,names(sim.deriv)[which(!(1:ncol(sim.deriv)) %in% cols.sims)]] <- sim.deriv[,which(!(1:ncol(sim.deriv)) %in% cols.sims)]
summary(deriv.stack)

sim.stack <- stack(sim.terms[,cols.sims])
names(sim.stack)  <- c("Y", "sim")
sim.stack[,names(sim.terms)[which(!(1:ncol(sim.terms)) %in% cols.sims)]] <- sim.terms[,which(!(1:ncol(sim.terms)) %in% cols.sims)]
summary(sim.stack)


vars.deriv <- c("deriv", "x")
deriv.agg                             <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "PFT", "Effect", "sim")], 
                                                   FUN=mean, na.rm=T)
deriv.agg[,paste0(vars.deriv, ".sd")] <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "PFT", "Effect", "sim")], 
                                                   FUN=sd, na.rm=T)[,vars.deriv]
summary(deriv.agg)

vars.sim <- c("Y", "x")
sim.agg                             <- aggregate(sim.stack[,vars.sim], 
                                                   by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "PFT", "Effect", "sim")], 
                                                   FUN=mean, na.rm=T)
sim.agg[,paste0(vars.sim, ".sd")] <- aggregate(sim.stack[,vars.sim], 
                                                   by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "PFT", "Effect", "sim")], 
                                                   FUN=sd, na.rm=T)[,vars.sim]
summary(sim.agg)

# Getting the mean derivative for each model & effet & PFT
vars.agg <- c("x", "mean.rel", "deriv")
ci.terms.agg <- aggregate(ci.terms[,vars.agg],
                          by=ci.terms[,c("Model", "Y.type", "data.type", "PFT", "Effect")],
                          FUN=mean, na.rm=T)
summary(ci.terms.agg)
}
# --------

# --------
# 4.b.1 Model agreement by PFT
# --------
vars.agg <- c("x", "mean.rel", "deriv")


pdf(file.path(fig.dir, "Sensitivity_versus_PFT_Model.pdf"), )
{
deriv.agg.graph <- deriv.agg
deriv.agg.graph[!is.na(deriv.agg.graph$deriv) & deriv.agg.graph$deriv>0.1,"deriv"] <- 0.1
deriv.agg.graph[!is.na(deriv.agg.graph$deriv) & deriv.agg.graph$deriv< -0.15,"deriv"] <- -0.15
deriv.agg.graph[deriv.agg.graph$Effect=="CO2" & !is.na(deriv.agg.graph$deriv) & deriv.agg.graph$deriv< -0.01,"deriv"] <- -0.01
deriv.agg.graph[deriv.agg.graph$Effect=="CO2" & !is.na(deriv.agg.graph$deriv) & deriv.agg.graph$deriv> 0.06,"deriv"] <- 0.06
print(
ggplot(data=deriv.agg.graph[deriv.agg.graph$data.type=="Model" & deriv.agg$Effect %in% c("tair", "precipf", "CO2"),]) +
  facet_grid(Effect ~ PFT, scales="free") +
  geom_violin(aes(x=Model, y=deriv, fill=Model), scale="width", width=3, position=position_dodge(width=0.5), alpha=0.5) +
  scale_fill_manual(values=colors.use) +
  theme_bw()
)}
dev.off()

# Differences in sensitivity among Tree RIng Width
co2.pft     <- lm(deriv ~ PFT, data=deriv.agg[deriv.agg$Model=="TreeRingRW" & deriv.agg$Effect=="CO2"    ,])
tair.pft    <- lm(deriv ~ PFT, data=deriv.agg[deriv.agg$Model=="TreeRingRW" & deriv.agg$Effect=="tair"   ,])
precipf.pft <- lm(deriv ~ PFT, data=deriv.agg[deriv.agg$Model=="TreeRingRW" & deriv.agg$Effect=="precipf",])
summary(co2.pft)
summary(tair.pft)
summary(precipf.pft)

# Aggregating to the effect level directly from deriv.terms.agg
# -- each model-PFT combo gets treated separately
effect.agg                           <- aggregate(ci.terms.agg[,vars.agg], 
                                                  by=ci.terms.agg[,c("Y.type", "data.type", "Effect")], 
                                                  FUN=mean)
effect.agg[,paste0(vars.agg, ".sd")] <- aggregate(ci.terms.agg[,vars.agg], 
                                                  by=ci.terms.agg[,c("Y.type", "data.type", "Effect")], 
                                                  FUN=sd)[,vars.agg]
effect.agg
effect.agg[effect.agg$data.type=="Model",]
effect.agg[effect.agg$data.type=="Tree Rings",]

clim.pft <- lm(deriv ~ Effect*(PFT-1)-Effect, data=ci.terms.agg[ci.terms.agg$data.type=="Model",])
summary(clim.pft); anova(clim.pft)

co2.pft     <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.pft    <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.pft <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])

summary(co2.pft)    ; anova(co2.pft)
summary(tair.pft)   ; anova(tair.pft)
summary(precipf.pft); anova(precipf.pft)

# Restricting comparisons to just Evergreens & Deciduous
co2.pft2     <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$PFT %in% c("Evergreen", "Deciduous") & ci.terms.agg$Effect=="CO2"    ,])
tair.pft2    <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$PFT %in% c("Evergreen", "Deciduous") & ci.terms.agg$Effect=="tair"   ,])
precipf.pft2 <- lm(deriv ~ PFT-1, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$PFT %in% c("Evergreen", "Deciduous") & ci.terms.agg$Effect=="precipf",])
anova(co2.pft2)
anova(tair.pft2)
anova(precipf.pft2)


ci.terms.agg.graph <- ci.terms.agg
ci.terms.agg.graph[!is.na(ci.terms.agg.graph$deriv) & ci.terms.agg.graph$Effect=="CO2" & ci.terms.agg.graph$deriv>0.1, "deriv"] <- 0.1
ci.terms.agg.graph[!is.na(ci.terms.agg.graph$deriv) & ci.terms.agg.graph$deri< -0.15, "deriv"] <- -0.15
summary(ci.terms.agg.graph)
ggplot(data=ci.terms.agg.graph[ci.terms.agg.graph$data.type=="Model" & ci.terms.agg.graph$Effect %in% c("tair", "precipf", "CO2"),]) +
  facet_wrap(~Effect) +
  geom_violin(aes(x=PFT, y=deriv), fill="gray30") +
  geom_point(aes(x=PFT, y=deriv, color=Model), position=position_jitter(width=0.1), size=5) +
  scale_color_manual(values=colors.use) +
#   scale_y_continuous(limits=c(-0.11, 0.11)) +
  theme_bw()

# --------

# --------
# 4.b.2 PFT agreement by Effect
# --------
library(nlme)
co2.mod     <- lme(deriv ~ PFT, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.mod    <- lme(deriv ~ PFT, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.mod <- lme(deriv ~ PFT, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
anova(co2.mod); summary(co2.mod)
anova(tair.mod); summary(tair.mod)
anova(precipf.mod); summary(precipf.mod)

co2.mod     <- lme(deriv ~ PFT, random=list(Model=~1, x=~1), data=deriv.stack[deriv.stack$data.type=="Model" & deriv.stack$Effect=="CO2"    ,])
tair.mod    <- lme(deriv ~ PFT, random=list(Model=~1, x=~1), data=deriv.stack[deriv.stack$data.type=="Model" & deriv.stack$Effect=="tair"   ,])
precipf.mod <- lme(deriv ~ PFT, random=list(Model=~1, x=~1), data=deriv.stack[deriv.stack$data.type=="Model" & deriv.stack$Effect=="precipf",])
anova(co2.mod); summary(co2.mod)
anova(tair.mod); summary(tair.mod)
anova(precipf.mod); summary(precipf.mod)



# Aggregating to effects & PFTs
vars.agg <- c("x", "mean.rel", "deriv")
pft.effect                              <- aggregate(ci.terms.agg[,vars.agg], 
                                                     by=ci.terms.agg[,c("Y.type", "data.type", "PFT", "Effect")], 
                                                     FUN=mean, na.rm=T)
pft.effect[,c(paste0(vars.agg, ".sd"))] <- aggregate(ci.terms.agg[,vars.agg], 
                                                     by=ci.terms.agg[,c("Y.type", "data.type", "PFT", "Effect")], 
                                                     FUN=sd, na.rm=T)[,vars.agg]

pft.effect <- pft.effect[order(pft.effect$data.type, pft.effect$Y.type, pft.effect$Effect, pft.effect$PFT),]
pft.effect[pft.effect$data.type=="Model"  & pft.effect$Effect %in% c("tair", "precipf", "CO2"),]


pft.effect <- pft.effect[order(pft.effect$data.type, pft.effect$Y.type, pft.effect$PFT, pft.effect$Effect),]
pft.effect[pft.effect$data.type=="Model" & pft.effect$Effect %in% c("tair", "precipf", "CO2"),]


# Getting the general PFT climate sensitivity stats
pft.agg                              <- aggregate(pft.effect[pft.effect$Effect %in% c("tair", "precipf", "CO2"), vars.agg], 
                                                  by=pft.effect[pft.effect$Effect %in% c("tair", "precipf", "CO2"), c("Y.type", "data.type", "PFT")], 
                                                  FUN=mean, na.rm=T)
pft.agg[,c(paste0(vars.agg, ".sd"))] <- aggregate(pft.effect[pft.effect$Effect %in% c("tair", "precipf", "CO2"), vars.agg], 
                                                  by=pft.effect[pft.effect$Effect %in% c("tair", "precipf", "CO2"), c("Y.type", "data.type", "PFT")], 
                                                  FUN=sd, na.rm=T)[,vars.agg]
pft.agg <- pft.agg[order(pft.agg$data.type, pft.agg$Y.type, pft.agg$PFT),]
pft.agg

# Comparing with tree ring data streams
ci.terms.agg <- ci.terms.agg[order(ci.terms.agg$Y.type),]
ci.terms.agg[ci.terms.agg$data.type=="Tree Rings" & ci.terms.agg$Effect %in% c("tair", "precipf", "CO2"),]


# --------

# --------


# ----------------------

# ----------------------
# Analyzing the differences
# ----------------------
# ----------------------
# ----------------------------------------
