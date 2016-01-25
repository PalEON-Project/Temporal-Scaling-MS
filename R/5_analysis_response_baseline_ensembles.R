# ----------------------------------------
# Objective: Compare NPP deviations through time & driver of "wiggles" among 
#            models, tree ring NPP, and tree ring RWI products
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
# 4. Graphing & Analyzing Raw NPP & RWI
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
in.base <- "Data/gamms/Sensitivity_Baseline"
out.dir <- "Data/analyses/response_baseline"
fig.dir <- "Figures/analyses/response_baseline"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# 2. Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]

load(file.path(in.base, "gamm_baseline_Models.Rdata"))
mod.models <- mod.out

load(file.path(in.base, "gamm_baseline_TreeRingNPP.Rdata"))
mod.npp <- mod.out

load(file.path(in.base, "gamm_baseline_TreeRings.Rdata"))
mod.tr <- mod.out

names.dat <- names(mod.models$data     )[names(mod.models$data     ) %in% names(mod.tr$data     )]
names.ci  <- names(mod.models$ci.terms )[names(mod.models$ci.terms ) %in% names(mod.tr$ci.terms )]
names.wt  <- names(mod.models$weights  )[names(mod.models$weights  ) %in% names(mod.tr$weights  )]
names.sim <- names(mod.models$sim.terms)[names(mod.models$sim.terms) %in% names(mod.tr$sim.terms)]

dat.ecosys <- rbind(cbind(mod.models$data[,  names.dat],  PlotID=NA, TreeID=NA, mod.models$ci.response[,c("mean", "lwr", "upr")]), 
                    cbind(mod.npp   $data[,c(names.dat , "PlotID")], TreeID=NA, mod.npp   $ci.response[,c("mean", "lwr", "upr")]), 
                    cbind(mod.tr    $data[,c(names.dat , "PlotID",  "TreeID")], mod.tr    $ci.response[,c("mean", "lwr", "upr")])
                    )
ci.terms   <- rbind(cbind(mod.models$ci.terms[,  names.ci],  PlotID=NA,  TreeID=NA), 
                    cbind(mod.npp   $ci.terms[,c(names.ci , "PlotID")],  TreeID=NA), 
                          mod.tr    $ci.terms[,c(names.ci , "PlotID"  , "TreeID")]
                    )
wt.terms   <- rbind(cbind(mod.models$weights[,  names.wt],  PlotID=NA,  TreeID=NA), 
                    cbind(mod.npp   $weights[,c(names.wt , "PlotID")],  TreeID=NA), 
                          mod.tr    $weights[,c(names.wt , "PlotID"  , "TreeID")]
                    )
sim.terms  <- rbind(cbind(mod.models$sim.terms[,  names.sim],  PlotID=NA,  TreeID=NA), 
                    cbind(mod.npp   $sim.terms[,c(names.sim , "PlotID")],  TreeID=NA), 
                          mod.tr    $sim.terms[,c(names.sim , "PlotID"  , "TreeID")])

sim.terms $Effect <- as.factor(sim.terms$Effect)
sim.terms $Extent <- as.factor(ifelse(sim.terms $Extent=="850-2010", "0850-2010", paste(ci.terms  $Extent)))
ci.terms  $Extent <- as.factor(ifelse(ci.terms  $Extent=="850-2010", "0850-2010", paste(ci.terms  $Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
wt.terms  $Extent <- as.factor(ifelse(wt.terms  $Extent=="850-2010", "0850-2010", paste(wt.terms  $Extent)))

# Making sure PlotID is showing up where it needs to
dat.ecosys$PlotID <- as.factor(dat.ecosys$PlotID)
wt.terms  $PlotID <- as.factor(wt.terms  $PlotID)
ci.terms  $PlotID <- as.factor(ci.terms  $PlotID)
sim.terms $PlotID <- as.factor(sim.terms $PlotID)

summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)
summary(sim.terms[,1:10])

# Use RWI as our tree ring measure
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="TreeRingBAI", ]
sim.terms  <- sim.terms [!sim.terms $Model=="TreeRingBAI", ]
ci.terms   <- ci.terms  [!ci.terms  $Model=="TreeRingBAI", ]
wt.terms   <- wt.terms  [!wt.terms  $Model=="TreeRingBAI", ]

# Double Check to make sure things are sorted by year
dat.ecosys[which(dat.ecosys$Model=="TreeRingRWI")[1:20],]
wt.terms  [which(wt.terms  $Model=="TreeRingRWI")[1:20],]


models.use <- unique(dat.ecosys[,"Model.Order"])
colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray30"))
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
			} else if(m=="TreeRingRWI") {
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
			} else if(m=="TreeRingRWI"){
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
dat.ecosys$Y.type <- as.factor(ifelse(dat.ecosys$Model=="TreeRingRWI", "RWI", "NPP"))
ci.terms  $Y.type <- as.factor(ifelse(ci.terms  $Model=="TreeRingRWI", "RWI", "NPP"))
wt.terms  $Y.type <- as.factor(ifelse(wt.terms  $Model=="TreeRingRWI", "RWI", "NPP"))
sim.terms $Y.type <- as.factor(ifelse(sim.terms $Model=="TreeRingRWI", "RWI", "NPP"))
} # End section block


summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])
# ----------------------------------------

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
tr.npp <- dat.ecosys[dat.ecosys$Model %in% c("TreeRingNPP", "TreeRingRWI") ,]
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
	ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRWI", "TreeRingBAI", "TreeRingNPP"),])  + facet_grid(Y.type~Site, scales="free_y", space="free") +
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
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRWI", "TreeRingBAI", "TreeRingNPP"),])  + facet_wrap(~Site) +
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
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRWI", "TreeRingBAI", "TreeRingNPP") & dat.ecosys$Site=="PHA",])  +
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
ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRWI", "TreeRingBAI", "TreeRingNPP") & dat.ecosys$Site %in% c("PHA", "PHO"),])  + facet_grid(Y.type~Site, space="free", scales="free_y") +
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
# 4.b.1 Range of NPP variability within & among models/data
# --------
# --------

# ---------------------

}
# ----------------------------------------


# ----------------------------------------
# 5. Graphing & Analyzing Sensitivity
# ----------------------------------------
{
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
pdf(file.path(fig.dir, "Fig2_NPP_Sensitivity_Models_Rel_Baseline.pdf"), height=8.5, width=11)
{
print(
ggplot(data=ci.terms.graph) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model.Order)) +
	scale_x_continuous(expand=c(0,0), name="") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
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
# 5.b.1. Range of change in NPP with each factor
# --------
# --------

# --------
# 5.b.2. Model closest to data for the ranges observed
# --------
# --------

# -----------------------

}
# ----------------------------------------



# ---------------------------------
# 6. Graphing & Analyzing ensemble means of Drivers of change and sensitivities through time
# ---------------------------------
# -----------------------
# 6.a. Graphing
# -----------------------
# combining wt.terms & dat.ecosys... they shouldn't need a merge because dimensions are the same
wt.terms2 <- cbind(wt.terms, dat.ecosys[,c("Model", "Model.Order", "Site", "Year", "Y", "Y.rel", "Y.10", "Y.rel.10")])
wt.terms2$Data.Type <- as.factor(ifelse(substr(wt.terms2$Model,1,8)=="TreeRing", paste(wt.terms2$Model), "Model"))
summary(wt.terms2)

indices.wt2 <- wt.terms2$Site=="PHA" & wt.terms2$Resolution=="t.001" & wt.terms2$Extent=="0850-2010" 

summary(wt.terms2)
summary(ci.terms)

factors.aggregate <- c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")

# --------
# 6.a.1. Aggregate by site
# --------
{
summary(wt.terms2[wt.terms2$Data.Type=="TreeRingRWI",])

ensemble.wts1 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "Data.Type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

ensemble.wts.lo <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "Data.Type", "Year")], FUN=quantile, 0.025, na.rm=T)
names(ensemble.wts.lo)[4:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[4:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Site", "Data.Type", "Year")], FUN=quantile, 0.975, na.rm=T)
names(ensemble.wts.hi)[4:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[4:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

dim(ensemble.wts1); dim(ensemble.wts.lo); dim(ensemble.wts.hi)

ensemble.wts.site <- cbind(ensemble.wts1, ensemble.wts.lo[,4:ncol(ensemble.wts.lo)], ensemble.wts.hi[,4:ncol(ensemble.wts.hi)])
summary(ensemble.wts.site)

# Fill weights in the smoothed data that are NA with 0
ensemble.wts.site[is.na(ensemble.wts.site$weight.tair.10   ),"weight.tair.10"   ] <- 0
ensemble.wts.site[is.na(ensemble.wts.site$weight.precipf.10),"weight.precipf.10"] <- 0
ensemble.wts.site[is.na(ensemble.wts.site$weight.CO2.10    ),"weight.CO2.10"    ] <- 0


summary(ensemble.wts.site)
}

pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_AllSites_1900-2010_Decadal.pdf"), width=11, height=8.5)
 source("R/5a_graph_ensembles_rgb_AllSites.R")	# THis is in a separate script because it's so dang long
dev.off()

pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_PHA_1850-2010_Decadal.pdf"), width=11, height=8.5)
{
print(
ggplot() + facet_grid(Data.Type~., space="free", scales="free_y") +
 	geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA",], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="Model","weight.tair"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="Model","weight.CO2"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="Model","weight.precipf"])), size=3) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingRWI",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingRWI","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingRWI","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingRWI","weight.precipf.10"])), size=3) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingNPP","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingNPP","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingNPP","weight.precipf.10"])), size=3) +
 	geom_hline(y=100, linetype="dashed") +
	scale_x_continuous(limits=c(1700,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts.site$Year), -2), round(max(ensemble.wts.site$Year), -2), by=100)) +
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
ensemble.wts1 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Model", "Site", "Data.Type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts1)

# Now get to 1 record per Model (regional patterns)
ensemble.wts2 <- aggregate(ensemble.wts1[,factors.aggregate], by=ensemble.wts1[,c("Model", "Data.Type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts2)

# Now agregate across models
ensemble.wts3 <- aggregate(ensemble.wts2[, factors.aggregate], by=ensemble.wts2[,c("Data.Type", "Year")], FUN=mean, na.rm=T)
summary(ensemble.wts3)

ensemble.wts.lo <- aggregate(ensemble.wts2[,factors.aggregate], by=ensemble.wts2[,c("Data.Type", "Year")], FUN=quantile, 0.025, na.rm=T)
names(ensemble.wts.lo)[3:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[3:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(ensemble.wts2[,factors.aggregate], by=ensemble.wts2[,c("Data.Type", "Year")], FUN=quantile, 0.975, na.rm=T)
names(ensemble.wts.hi)[3:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[3:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

ensemble.wts <- cbind(ensemble.wts3, ensemble.wts.lo[,3:ncol(ensemble.wts.lo)], ensemble.wts.hi[,3:ncol(ensemble.wts.hi)])
summary(ensemble.wts)

# Fill weights in the smoothed data that are NA with 0
ensemble.wts[is.na(ensemble.wts$weight.tair.10   ),"weight.tair.10"   ] <- 0
ensemble.wts[is.na(ensemble.wts$weight.precipf.10),"weight.precipf.10"] <- 0
ensemble.wts[is.na(ensemble.wts$weight.CO2.10    ),"weight.CO2.10"    ] <- 0
}
summary(ensemble.wts)

ggplot() + facet_grid(Data.Type~., space="free", scales="free_y") +
 	geom_ribbon(data= ensemble.wts[,], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts[ensemble.wts$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Data.Type=="Model","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="Model","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="Model","weight.precipf"])), size=3) +
	geom_line(data= ensemble.wts[ensemble.wts$Data.Type=="TreeRingRWI",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingRWI","weight.tair.10"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingRWI","weight.CO2.10"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingRWI","weight.precipf.10"])), size=3) +
	geom_line(data= ensemble.wts[ensemble.wts$Data.Type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingNPP","weight.tair.10"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingNPP","weight.CO2.10"]),
                        abs(ensemble.wts[ensemble.wts$Data.Type=="TreeRingNPP","weight.precipf.10"])), size=3) +
 	geom_hline(y=100, linetype="dashed") +
	scale_x_continuous(limits=c(1700,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts$Year), -2), round(max(ensemble.wts$Year), -2), by=100)) +
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
# --------

### --------------------------------
### Note: See script 5a for graphing all sites (too big to put in here)
### --------------------------------

# -----------------------


# -----------------------
# 6.b. Quantitative Analysis
# -----------------------
# -----------------------

# ----------------------------------------






