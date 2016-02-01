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

# dat.ecosys <- cbind(mod.out$data[, ], mod.out$ci.response[,c("mean", "lwr", "upr")])
# ci.terms   <- mod.out$ci.terms[,]
# wt.terms   <- mod.out$weights[,]
# sim.terms  <- mod.out$sim.terms
# sim.terms $Effect <- as.factor(sim.terms$Effect)

# # Grouping the kind and source of the data
# dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRW", "RW", "NPP"))
# ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRW", "RW", "NPP"))
# wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRW", "RW", "NPP"))
# sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRW", "RW", "NPP"))

# dat.ecosys$data.type <- as.factor(ifelse(substr(dat.ecosys$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# ci.terms  $data.type <- as.factor(ifelse(substr(ci.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# wt.terms  $data.type <- as.factor(ifelse(substr(wt.terms  $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
# sim.terms $data.type <- as.factor(ifelse(substr(sim.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))

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

# # Make sure all data sets are ordered by year, then treeID, then plotID, then Model
# # sort.order <- c("Model", "PlotID", "TreeID", "Year")
# dat.ecosys <- dat.ecosys[order(dat.ecosys$Model, dat.ecosys$PlotID, dat.ecosys$TreeID, dat.ecosys$Year),]
# wt.terms   <- wt.terms[order(wt.terms$Model, wt.terms$PlotID, wt.terms$TreeID, wt.terms$Year),]

# # Double Check to make sure things are sorted by year so rollapply works
# dat.ecosys[which(dat.ecosys$Model=="TreeRingRW")[1:20],]
# wt.terms  [which(wt.terms  $Model=="TreeRingRW")[1:20],]

# summary(dat.ecosys)
# summary(wt.terms)

# {
# for(m in unique(ci.terms$Model)){

		# # -----------------------
		# # 3.a. Find the NPP to relativize each set off of
		# # Using mean model NPP across sites since the GAMM response curves are for 
		# #    the whole model & not site-specific are parameterized
		# # -----------------------
		# # Find the start year for the extent
		# # yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

		# npp <- mean(dat.ecosys[dat.ecosys$Model==m, "Y"], na.rm=T)			
		# # -----------------------
		
		# # -----------------------
		# # 3.b Relativizing everything in dat.ecosys to make it comparable to tree rings
		# # -----------------------
		# {		
		# # Which factors to relativize
		# y.rel <- c("Y", "fit.gam", "mean", "lwr", "upr")

		# # for some reason, I can create multiple new columns at once
		# # Solution: use a loop to create blank columns and then fill them
		# for(y in y.rel){
			# dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel"       )] <- NA	
			# dat.ecosys[dat.ecosys$Model==m,paste0(y, ".10"        )] <- NA	
			# dat.ecosys[dat.ecosys$Model==m,paste0(y, ".rel", ".10")] <- NA	
		# }
		# dat.ecosys[dat.ecosys$Model==m,paste0(y.rel, ".rel")] <- dat.ecosys[dat.ecosys$Model==m, y.rel]/npp
		
		# # Getting 10-year running means to make clearer figures
		# for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
			# # Note: If we're working with tree ring data, we need to go by plot for NPP products 
			# #       & by Tree for individual-level tree rings products
			# if(m=="TreeRingNPP"){
				# for(p in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "PlotID"])){

					# # Raw NPP (to add dark line over faded annual wiggles)
					# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# # Relativized NPP (to have generalized patterns for figures)
					# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$PlotID==p,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$PlotID==p, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				# }
			# } else if(m=="TreeRingRW") {
				# for(t in unique(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, "TreeID"])){
					# # If we have too few data points, we need to skip that tree 
					# if(length(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel[1]]) < 10) next

					# # Raw NPP (to add dark line over faded annual wiggles)
					# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

					# # Relativized NPP (to have generalized patterns for figures)
					# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s & dat.ecosys$TreeID==t,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s  & dat.ecosys$TreeID==t, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
				# }
			# } else {
				# # Raw NPP (to add dark line over faded annual wiggles)
				# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, y.rel], FUN=mean, width=10, align="center", fill=NA, by.column=T)

				# # Relativized NPP (to have generalized patterns for figures)
				# dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s,paste0(y.rel, ".rel", ".10" )] <- rollapply(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Site==s, paste0(y.rel, ".rel")], FUN=mean, width=10, align="center", fill=NA, by.column=T)
			# }
		# }
		# }		
		# # -----------------------

		
		# # -----------------------
		# # 3.c. Finding the percent change in NPP relative to the mean for that particular scale
		# # -----------------------
		# {
		# y.rel <- c("mean", "lwr", "upr")
		# for(y in y.rel){
			# ci.terms[ci.terms$Model==m,paste0(y, ".rel"       )] <- NA	
		# }		

		# ci.terms[ci.terms$Model==m,paste0(y.rel,".rel")] <- ci.terms[ci.terms $Model==m, y.rel]/npp
		
		# # Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
		# cols.sim <- which(substr(names(sim.terms),1,1)=="X")
		# sim.terms[sim.terms$Model==m,cols.sim] <- sim.terms[sim.terms$Model==m,cols.sim]/npp
		# }
		# # -----------------------

		# # -----------------------
		# # 3.d. Relativizing the factor fits through times and weights as well
		# # Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
		# # -----------------------
		# {
		# y.rel <- c("fit.tair", "fit.precipf", "fit.CO2")
		# # for some reason, I can create multiple new columns at once
		# # Solution: use a loop to create blank columns and then fill them
		# for(y in y.rel){
			# wt.terms[wt.terms$Model==m,paste0(y, ".rel"       )] <- NA	
			# wt.terms[wt.terms$Model==m,paste0(y, ".rel", ".10")] <- NA	
		# }
		
		# wt.terms[wt.terms$Model==m,paste0(y.rel, ".rel")] <- 1+(wt.terms[wt.terms$Model==m,y.rel])/npp
		
		# # We only really care about smoothing the relativized weights
		# y.rel2 <- c(y.rel, paste0(y.rel, ".rel"), "weight.tair", "weight.precipf", "weight.CO2")
		# for(y in y.rel2){
			# wt.terms[wt.terms$Model==m,paste0(y, ".10")] <- NA	
		# }

		# # Getting 10-year running means to make clearer figures
		# for(s in unique(dat.ecosys[dat.ecosys$Model==m, "Site"])){
			# if(m=="TreeRingNPP"){
				# for(p in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "PlotID"])){
					# # Relativized NPP (to have generalized patterns for figures)
					# wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$PlotID==p, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
				# }
			# } else if(m=="TreeRingRW"){
				# for(t in unique(wt.terms[wt.terms$Model==m & wt.terms$Site==s, "TreeID"])){

					# # If we have too few data points, we need to skip that tree 
					# if(length(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2[1]]) < 10) next

					# # Relativized NPP (to have generalized patterns for figures)
					# wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s & wt.terms$TreeID==t, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)			
				# }				
			# } else {
				# # Relativized NPP (to have generalized patterns for figures)
				# wt.terms[wt.terms$Model==m & wt.terms$Site==s,paste0(y.rel2, ".10" )] <- rollapply(wt.terms[wt.terms$Model==m & wt.terms$Site==s, y.rel2], FUN=mean, width=10, align="center", fill=NA, by.column=T)
			# }
		# }
		# }
		# # -----------------------

# }
# } # End section block


# summary(dat.ecosys)
# summary(ci.terms)
# summary(wt.terms)
# summary(sim.terms[,1:10])

# save(dat.ecosys, ci.terms, wt.terms, sim.terms, file=file.path(out.dir, "post-process_baseline.RData"))
}
# ----------------------------------------

load(file.path(out.dir, "post-process_baseline.RData"))

# ----------------------------------------
# 4. Graphing & Summary Statistics of Raw NPP
# ----------------------------------------
{
summary(dat.ecosys)

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
# 4.a. Figures
# ---------------------
{
# --------
# 4.a.1. All Sites 1800-2010 (Fig. 1)
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

pdf(file.path(fig.dir, "Fig1_NPP_Raw_AllSites_1800-2010_Simple.pdf"), height=8.5, width=11)
grid.newpage()
grid.draw(Fig1)
dev.off()

# --------

# --------
# 4.a.2. Models Only, All Sites 0850-2010 (Supplemental Figure 1)
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
# 4.a.3. Just PHA 1900-2010
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
# 4.a.4. PHA, PHO (locations with NPP products, 1950-2010)
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
# 4.b. Summary Statistics
# ---------------------
# --------
# 4.b.1 Description of Model variability across space & Time
# --------
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

# --------

# --------
# 4.b.2. Description of Site variability
# --------
site.stats      <- aggregate(models.sites[,c("Y", "Y.rel")], 
                             by=models.sites[,c("Site", "Y.type", "data.type")], 
                             FUN=mean)
site.stats[,c("Y.sd", "Y.rel.sd")]   <- aggregate(models.sites[,c("Y", "Y.rel")], 
                                                     by=models.sites[,c("Site", "Y.type", "data.type")], 
                                                     FUN=sd)[,c("Y", "Y.rel")]
site.stats
# --------


# ---------------------

}
# ----------------------------------------


# ----------------------------------------
# 5. Graphing & Analyzing Sensitivity
# ----------------------------------------
{
summary(ci.terms)
summary(sim.terms[,1:10])
# -----------------------
# 5.a. Graphing
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
pdf(file.path(fig.dir, "Fig2_Sensitivity_Models_Rel_Baseline.pdf"), height=8.5, width=11)
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
}
dev.off()

} # End graphing section
# -----------------------


# -----------------------
# 5.b. Quantitative Analysis
# -----------------------
# --------
# 5.b.0. Getting the first derivative (first diference) of each line so we can take the mean slope
# --------
{
# First make sure the effects are sorted by x to make this easier
ci.terms  <- ci.terms [order(ci.terms $Model, ci.terms $Effect, ci.terms $x),]
sim.terms <- sim.terms[order(sim.terms$Model, sim.terms$Effect, sim.terms$x),]
ci.terms[1:20,1:10]
dim(ci.terms)
summary(ci.terms)

# Making a new dataframe dedicated to the derivatives
cols.sims <- which(substr(names(sim.terms),1,1)=="X")
sim.deriv <- sim.terms[,]
sim.deriv[,cols.sims] <- NA

for(e in unique(ci.terms$Effect)){
 for(m in unique(ci.terms$Model)){
   x.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "x"], lag=1), NA)
   y.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "mean.rel"], lag=1), NA)
   ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "deriv"] <- y.dif/x.dif
 
   # For the full simiulation for robust analysis
   y.dif2 <- rbind(apply(sim.terms[sim.terms$Model==m & sim.terms$Effect==e, cols.sims], 2, FUN=diff), NA)
   x.dif <- c(diff(sim.terms[sim.terms$Model==m & sim.terms$Effect==e, "x"], lag=1), NA)
   sim.deriv[sim.deriv$Model==m & sim.deriv$Effect==e, cols.sims] <- apply(y.dif2, 2, FUN=function(y){y/x.dif})
 } 
}

# Stacking and aggregating the simulations
deriv.stack <- stack(sim.deriv[,cols.sims])
names(deriv.stack)  <- c("deriv", "sim")
deriv.stack[,names(sim.deriv)[which(!(1:ncol(sim.deriv)) %in% cols.sims)]] <- sim.deriv[,which(!(1:ncol(sim.deriv)) %in% cols.sims)]
summary(deriv.stack)

vars.deriv <- c("deriv", "x")
deriv.agg                             <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Effect", "sim")], 
                                                   FUN=mean, na.rm=T)
deriv.agg[,paste0(vars.deriv, ".sd")] <- aggregate(deriv.stack[,vars.deriv], 
                                                   by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Effect", "sim")], 
                                                   FUN=sd, na.rm=T)[,vars.deriv]
summary(deriv.agg)

# Condensing model variability across space and time (since that's what we did for calculating sensitivity)
vars.agg <- c("Y", "Y.rel", "Time")
mod.agg                           <- aggregate(dat.ecosys[,vars.agg], 
                                               by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                               FUN=mean)
mod.agg[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys[,vars.agg], 
                                               by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                               FUN=sd)[,vars.agg]
mod.agg


# Finding the change in key variables in the modern era
for(v in vars.agg){
  mod.agg[mod.agg$Model==m,paste0("dModern.", v)] <- NA  
}

for(m in unique(mod.agg$Model)){
  mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1990 & dat.ecosys$Year<=2010, vars.agg], na.rm=T) -
                                                             colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1830 & dat.ecosys$Year<=1850, vars.agg], na.rm=T)
}
mod.agg

# Combining mod.agg and the mean deriviative for each variable and model
summary(ci.terms)
ci.terms.agg <- aggregate(ci.terms[,c("mean.rel", "deriv")], by=ci.terms[,c("Model", "Effect")], FUN=mean, na.rm=T)
ci.terms.agg <- merge(mod.agg, ci.terms.agg, all.x=T, all.y=T)
summary(ci.terms.agg)


deriv.agg <- merge(deriv.agg, mod.agg, all.x=T, all.y=T)

# --------

# --------
# 5.b.1. Sensitivity comparisons as a function of baseline NPP & variability
#  -- Note: Here we're comparing characteristics of models, so we're going to go with 
#           the data aggregated so that there's one value per model.  (We will NOT use
#           the simulations)
# --------
{
# Sensitivity range vs. NPP
summary(ci.terms.agg)
summary(deriv.agg)

# Correlations with mean NPP
pdf(file.path(fig.dir, "Sensitivity_versus_NPP_Models.pdf"), height=8.5, width=11)
ggplot(data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect %in% c("tair", "precipf", "CO2"),]) + 
  facet_wrap(~Effect, scales="free") +
  geom_point(aes(x=Y, y=deriv, color=Model), size=10) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  stat_smooth(aes(x=Y, y=deriv), method="lm", size=2, color="black") +
  theme_bw()
dev.off()

co2.npp     <- lm(deriv ~ Y, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.npp    <- lm(deriv ~ Y, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.npp <- lm(deriv ~ Y, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
summary(co2.npp)
summary(tair.npp)
summary(precipf.npp)


# Correlations with NPP variability
summary(ci.terms.agg)
pdf(file.path(fig.dir, "Sensitivity_versus_Variability_Models.pdf"), height=8.5, width=11)
ggplot(data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect %in% c("tair", "precipf", "CO2"),]) + 
  facet_wrap(~Effect, scales="free") +
  geom_point(aes(x=Y.rel.sd, y=deriv, color=Model), size=10) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  stat_smooth(aes(x=Y.rel.sd, y=deriv), method="lm", size=2, color="black") +
  theme_bw()
dev.off()

co2.sd     <- lm(deriv ~ Y.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.sd    <- lm(deriv ~ Y.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.sd <- lm(deriv ~ Y.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
summary(co2.sd)
summary(tair.sd)
summary(precipf.sd)

# Correlations with modern change in NPP
summary(ci.terms.agg)
pdf(file.path(fig.dir, "Sensitivity_versus_Change_Models.pdf"), height=8.5, width=11)
ggplot(data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect %in% c("tair", "precipf", "CO2"),]) + 
  facet_wrap(~Effect, scales="free") +
  geom_point(aes(x=dModern.Y.rel, y=deriv, color=Model), size=10) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  stat_smooth(aes(x=dModern.Y.rel, y=deriv), method="lm", size=2, color="black") +
  theme_bw()
dev.off()

co2.mod     <- lm(deriv ~ dModern.Y.rel, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.mod    <- lm(deriv ~ dModern.Y.rel, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.mod <- lm(deriv ~ dModern.Y.rel, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
summary(co2.mod)
summary(tair.mod)
summary(precipf.mod)


# A couple additional stats on sensitivity agreement 
summary(ci.terms)
clim.agg <- c("mean.rel", "deriv")
clim.effects                           <- aggregate(ci.terms.agg[,clim.agg], 
                                                    by=ci.terms.agg[,c("Y.type", "data.type", "Effect")], 
                                                    FUN=mean, na.rm=T)
clim.effects[,paste0(clim.agg, ".sd")] <- aggregate(ci.terms.agg[,clim.agg], 
                                                    by=ci.terms.agg[,c("Y.type", "data.type", "Effect")], 
                                                    FUN=sd, na.rm=T)[,clim.agg]
clim.effects
clim.effects[clim.effects$data.type=="Model",]
}
# --------

# --------
# 5.b.2. Sensitivity as a function of disturbance
# --------
# Locating the fire file that didn't get put into ecosys
output.all <- read.csv("../../phase1a_output_variables/PalEON_MIP_Yearly.csv")
summary(output.all)

agg.fire <- aggregate(output.all[,c("Fire")], by=output.all[,c("Model", "Updated")], FUN=mean, na.rm=T)
names(agg.fire)[3] <- "Fire"
agg.fire[is.na(agg.fire$Fire),] <- 0
agg.fire

ci.terms.agg <- merge(ci.terms.agg, agg.fire[,c("Model", "Fire")], all.x=T, all.y=F)
ci.terms.agg

# Correlations with NPP variability
summary(ci.terms.agg)
pdf(file.path(fig.dir, "Sensitivity_versus_Fire_Models.pdf"), height=8.5, width=11)
ggplot(data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect %in% c("tair", "precipf", "CO2"),]) + 
  facet_wrap(~Effect, scales="free") +
  geom_point(aes(x=Fire, y=deriv, color=Model), size=10) +
  scale_fill_manual(values=colors.use) +
  scale_color_manual(values=colors.use) +
  stat_smooth(aes(x=Fire, y=deriv), method="lm", size=2, color="black") +
  theme_bw()
dev.off()

co2.fire     <- lm(deriv ~ Fire, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2"    ,])
tair.fire    <- lm(deriv ~ Fire, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair"   ,])
precipf.fire <- lm(deriv ~ Fire, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
summary(co2.fire)
summary(tair.fire)
summary(precipf.fire)

summary(co2.fire)
summary(tair.fire)
summary(precipf.fire)
# --------

# -----------------------

}
# ----------------------------------------



# ---------------------------------
# 6. Graphing & Analyzing ensemble means of Drivers of change and sensitivities through time
# ---------------------------------
{
# -----------------------
# 6.a. Graphing
# -----------------------
# combining wt.terms & dat.ecosys... they shouldn't need a merge because dimensions are the same
wt.terms2 <- cbind(wt.terms, dat.ecosys[,c("Model.Order","Y", "Y.rel", "Y.10", "Y.rel.10")])

# change the data type
wt.terms2$data.type <- as.factor(ifelse(substr(wt.terms2$Model,1,8)=="TreeRing", paste(wt.terms2$Model), "Model"))
summary(wt.terms2)

# --------------------------
# Adjusting CO2 Effect
# --------------------------
# Note: because the gam makes the smoother cross 0 at the MEAN CO2 (which is in the 1800s), 
# it's saying the region is pretty CO2-limited at periods where that doesn't really make 
# sense, so we're going to off relativize it to whatever the starting point for the run is
# --------------------------
{
for(m in unique(wt.terms2$Model)){
	yr1       <- min(wt.terms2[wt.terms2$Model==m, "Year"]) # find the minimum year
	yr2       <- min(wt.terms2[wt.terms2$Model==m & !is.na(wt.terms2$weight.CO2), "Year"]) # find the minimum year
	co2.base <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr1+5), "fit.CO2"], na.rm=T) # mean of the 5 years around the starting 
	co2.base.10 <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr2+5),"fit.CO2.10"],na.rm=T) # mean of the 5 years around the starting point
	wt.terms2[wt.terms2$Model==m, "fit.CO2.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2"] - co2.base
	wt.terms2[wt.terms2$Model==m, "fit.CO2.10.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.10"] - co2.base.10
}
summary(wt.terms2)


ggplot(data=wt.terms2[wt.terms2$Model=="lpj.guess",]) + 
	facet_wrap(~PlotID, scales="free") +
	scale_x_continuous(limits=c(1500,2010)) +
	geom_line(aes(x=Year, y=fit.CO2.10), color="green3", size=1, linetype="dashed") +
	geom_line(aes(x=Year, y=fit.CO2.10.adj), color="green3", size=2) +
	geom_line(aes(x=Year, y=fit.precipf.10), color="blue2") +
	geom_line(aes(x=Year, y=fit.tair.10), color="red2") +
	theme_bw()


wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")] <- abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")])/rowSums(abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")]))
wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")] <- abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")])/rowSums(abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")]))

ggplot(data=wt.terms2[wt.terms2$Model=="lpj.guess",]) + 
	facet_wrap(~PlotID, scales="free") +
	scale_x_continuous(limits=c(1500,2010)) +
	# geom_line(aes(x=Year, y=abs(weight.CO2.10)), color="green3", size=1, linetype="dashed") +
	# geom_line(aes(x=Year, y=abs(weight.precipf.10)), color="blue2", size=0.5, linetype="dashed") +
	# geom_line(aes(x=Year, y=abs(weight.tair.10)), color="red2", size=0.5, linetype="dashed") +
	geom_line(aes(x=Year, y=weight.CO2.10.adj), color="green3") +
	geom_line(aes(x=Year, y=weight.precipf.10.adj), color="blue2") +
	geom_line(aes(x=Year, y=weight.tair.10.adj), color="red2") +
	theme_bw()


wt.terms2[is.na(wt.terms2$weight.tair.10       ),"weight.tair.10"       ] <- 0
wt.terms2[is.na(wt.terms2$weight.precipf.10    ),"weight.precipf.10"    ] <- 0
wt.terms2[is.na(wt.terms2$weight.CO2.10        ),"weight.CO2.10"        ] <- 0
wt.terms2[is.na(wt.terms2$weight.tair.10.adj   ),"weight.tair.10.adj"   ] <- 0
wt.terms2[is.na(wt.terms2$weight.precipf.10.adj),"weight.precipf.10.adj"] <- 0
wt.terms2[is.na(wt.terms2$weight.CO2.10.adj    ),"weight.CO2.10.adj"    ] <- 0

summary(rowSums(wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")]))
}
# --------------------------

summary(wt.terms2)
summary(ci.terms)

# factors.aggregate <- c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10", "weight.CO2.adj",)
factors.aggregate <- c("fit.full", "Y.rel", "Y.rel.10", "weight.tair.adj", "weight.precipf.adj", "weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")
# factors.aggregate <- c("fit.full", "Y.rel", "Y.rel.10", "weight.tair", "weight.precipf", "weight.CO2", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")

# --------
# 6.a.1. Aggregate by site
# --------
{
summary(wt.terms2[wt.terms2$data.type=="TreeRingRW",])

ensemble.wts1 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

ensemble.wts.lo <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year")], FUN=quantile, 0.025, na.rm=T)
names(ensemble.wts.lo)[4:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[4:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "data.type", "Year")], FUN=quantile, 0.975, na.rm=T)
names(ensemble.wts.hi)[4:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[4:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

dim(ensemble.wts1); dim(ensemble.wts.lo); dim(ensemble.wts.hi)

ensemble.wts.site <- cbind(ensemble.wts1, ensemble.wts.lo[,4:ncol(ensemble.wts.lo)], ensemble.wts.hi[,4:ncol(ensemble.wts.hi)])
summary(ensemble.wts.site)

# Re-normalizing factor weights
summary(rowSums(abs(ensemble.wts.site[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")])))

# wts.sum.10 <- abs(ensemble.wts.site$weight.tair.10.adj) + abs(ensemble.wts.site$weight.precipf.10.adj) + abs(ensemble.wts.site$weight.CO2.10.adj)
# ensemble.wts.site[,c("weight.tair.10.adj","weight.precipf.10.adj", "weight.CO2.10.adj")] <- ensemble.wts.site[,c("weight.tair.10.adj","weight.precipf.10.adj", "weight.CO2.10.adj")]/wts.sum.10
# ensemble.wts.site[is.na(ensemble.wts.site$weight.tair.10.adj   ),"weight.tair.10.adj"   ] <- 0
# ensemble.wts.site[is.na(ensemble.wts.site$weight.precipf.10.adj),"weight.precipf.10.adj"] <- 0
# ensemble.wts.site[is.na(ensemble.wts.site$weight.CO2.10.adj    ),"weight.CO2.10.adj"    ] <- 0

# wts.sum <- abs(ensemble.wts.site$weight.tair) + abs(ensemble.wts.site$weight.precipf) + abs(ensemble.wts.site$weight.CO2)
# ensemble.wts.site[,c("weight.tair","weight.precipf", "weight.CO2")] <- ensemble.wts.site[,c("weight.tair","weight.precipf", "weight.CO2")]/wts.sum
# summary(ensemble.wts.site)


# summary(ensemble.wts.site)
} # End aggregation

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_AllSites_1850-2010_Decadal.pdf"), width=11, height=8.5)
 source("R/5a_graph_ensembles_rgb_AllSites.R")	# THis is in a separate script because it's so dang long
dev.off()

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_PHA_1700-2010_Decadal.pdf"), width=11, height=8.5)
{
print(
ggplot() + facet_grid(data.type~., scales="free_y") +
	scale_x_continuous(limits=c(1700,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts.site$Year), -2), round(max(ensemble.wts.site$Year), -2), by=100)) +
 	geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA",], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="Model","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="Model","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="Model","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingRW","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingRW","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingRW","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingNPP","weight.tair.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingNPP","weight.CO2.10.adj"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$data.type=="TreeRingNPP","weight.precipf.10.adj"])), size=3) +
 	geom_hline(y=100, linetype="dashed") +
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

# --------

# --------
# 6.a.s. Aggregate to model first across models to get a regional signal
#        Note: Going to aggregate to site first because of the tree rings so that they aren't quite so spatially skewed
# --------
{
# First get everything to the 1 record per site level
ensemble.wts1 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Model", "Site", "data.type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

# Now get to 1 record per Model (regional patterns)
ensemble.wts2 <- aggregate(ensemble.wts1[,factors.aggregate], by=ensemble.wts1[,c("Model", "data.type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts2)

# Now agregate across models
ensemble.wts3 <- rbind(
				aggregate(ensemble.wts2[ensemble.wts2$data.type=="Model", factors.aggregate], 
					   by=ensemble.wts2[ensemble.wts2$data.type=="Model",c("data.type", "Year")], 
					   FUN=mean, na.rm=T),
				aggregate(ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing", factors.aggregate], 
					   by=ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing",c("data.type", "Year")], 
					   FUN=mean, na.rm=T)
					      )
summary(ensemble.wts3)

ensemble.wts.lo <- rbind(
				aggregate(ensemble.wts2[ensemble.wts2$data.type=="Model", factors.aggregate], 
					   by=ensemble.wts2[ensemble.wts2$data.type=="Model",c("data.type", "Year")], 
					   FUN=quantile, 0.025, na.rm=T),
				aggregate(ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing", factors.aggregate], 
					   by=ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing",c("data.type", "Year")], 
					   FUN=quantile, 0.025, na.rm=T)
					      )
names(ensemble.wts.lo)[3:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[3:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- rbind(
				aggregate(ensemble.wts2[ensemble.wts2$data.type=="Model", factors.aggregate], 
					   by=ensemble.wts2[ensemble.wts2$data.type=="Model",c("data.type", "Year")], 
					   FUN=quantile, 0.975, na.rm=T),
				aggregate(ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing", factors.aggregate], 
					   by=ensemble.wts1[substr(ensemble.wts1$data.type,1,8)=="TreeRing",c("data.type", "Year")], 
					   FUN=quantile, 0.975, na.rm=T)
					      )
names(ensemble.wts.hi)[3:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[3:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

ensemble.wts <- cbind(ensemble.wts3, ensemble.wts.lo[,3:ncol(ensemble.wts.lo)], ensemble.wts.hi[,3:ncol(ensemble.wts.hi)])
summary(ensemble.wts)

# Fill weights in the smoothed data that are NA with 0
summary(rowSums(abs(ensemble.wts.site[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")])))

# ensemble.wts[is.na(ensemble.wts$weight.tair.10   ),"weight.tair.10"   ] <- 0
# ensemble.wts[is.na(ensemble.wts$weight.precipf.10),"weight.precipf.10"] <- 0
# ensemble.wts[is.na(ensemble.wts$weight.CO2.10    ),"weight.CO2.10"    ] <- 0
} # End aggregation
summary(ensemble.wts) 

# Re-normalizing factor weights
# wts.sum.10 <- abs(ensemble.wts$weight.tair.10) + abs(ensemble.wts$weight.precipf.10) + abs(ensemble.wts$weight.CO2.10)
# ensemble.wts[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")] <- ensemble.wts[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")]/wts.sum.10
# ensemble.wts[is.na(ensemble.wts$weight.tair.10   ),"weight.tair.10"   ] <- 0
# ensemble.wts[is.na(ensemble.wts$weight.precipf.10),"weight.precipf.10"] <- 0
# ensemble.wts[is.na(ensemble.wts$weight.CO2.10    ),"weight.CO2.10"    ] <- 0

# wts.sum <- abs(ensemble.wts$weight.tair) + abs(ensemble.wts$weight.precipf) + abs(ensemble.wts$weight.CO2)
# ensemble.wts[,c("weight.tair","weight.precipf", "weight.CO2")] <- ensemble.wts[,c("weight.tair","weight.precipf", "weight.CO2")]/wts.sum
summary(ensemble.wts)



pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_1500-2010_Decadal.pdf"), width=11, height=8.5)
{
ggplot() + facet_grid(data.type~., space="free", scales="free_y") +
	scale_x_continuous(limits=c(1500,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts$Year), -2), round(max(ensemble.wts$Year), -2), by=100)) +
 	geom_ribbon(data= ensemble.wts[,], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts[ensemble.wts$data.type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$data.type=="Model","weight.tair.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="Model","weight.CO2.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="Model","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts[ensemble.wts$data.type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$data.type=="TreeRingRW","weight.tair.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="TreeRingRW","weight.CO2.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="TreeRingRW","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts[ensemble.wts$data.type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$data.type=="TreeRingNPP","weight.tair.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="TreeRingNPP","weight.CO2.10.adj"]),
                        abs(ensemble.wts[ensemble.wts$data.type=="TreeRingNPP","weight.precipf.10.adj"])), size=3) +
 	geom_hline(y=100, linetype="dashed") +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0), breaks=c(50, 100, 150, 200)) +
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
}
dev.off()
# --------

### --------------------------------
### Note: See script 5a for graphing all sites (too big to put in here)
### --------------------------------

# -----------------------


# -----------------------
# 6.b. Quantitative Analysis
# -----------------------
# -----------------------

}
# ----------------------------------------






