# ----------------------------------------
# Temporal Scaling Analyses
# Looking at whether temporal resolution or extent has greater impact on model responses to drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Which is the more important consideration in understanding ecosystem variation in response to drivers
#           of change: temporal resolution, or temporal scale?
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Individual model response curves to driver
#    -- generated with 2c_process_drivers_all_drivers.R & 2e_process_drivers_all_drivers_byExtent.R
#    -- NOTE: Taking growing season analysis because it had a slightly higher explanatory power
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) Determine & compare relative weight of individual drivers across model
#    -- pay special attention to the role of Tair, Precip & CO2
# B) Determine the change in driver response with change in scale
#    -- relativize to the annual res, 850-2010 extent curves
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(car)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
# setwd("..")
path.data <- "Data"
in.base <- "Data/gamms"
out.dir <- "Data/analyses/response_ensemble_baseline"
fig.dir <- "Figures/analyses/response_ensemble_baseline"

# Prescribing our baseline resolution
resolution="t.001"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]


load(file.path(in.base, "Sensitivity_All_TempRes/gamm_Models_NPP_Resolutions_ExtentFull.Rdata"))

ci.terms   <- mod.out$ci.terms
dat.ecosys <- cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")])
wt.terms   <- mod.out$weights
sim.terms  <- mod.out$sim.terms

sim.terms$Effect  <- as.factor(sim.terms$Effect)
ci.terms$Extent   <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
wt.terms$Extent   <- as.factor(ifelse(wt.terms$Extent=="850-2010", "0850-2010", paste(wt.terms$Extent)))

summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)
summary(sim.terms[,1:10])


models.use <- unique(dat.ecosys[,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])


pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Raw_Site.pdf"), height=8.5, width=11)
# png(file.path(fig.dir, "NPP_Sensitivity_Models_Raw_NoSite.png"))
print(
ggplot(data=ci.terms[ci.terms$Resolution==resolution,]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off()


pdf(file.path(fig.dir, "PHA_Evergreen_quick.pdf"))
print(
ggplot(data=ecosys[ecosys$Resolution=="t.050" & ecosys$Site=="PUN",])  +
	geom_line(aes(x=Year, y=Evergreen*100, color=Model.Order), size=3) +
	scale_y_continuous(name="Precent Evergreen", limits=c(25,100)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	guides(color=F) +
	theme(legend.title=element_text(size=rel(2), face="bold"),
	      legend.text=element_text(size=rel(1.5)),
	      # legend.position=c(0.2, 0.18),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text=element_text(size=rel(2.5), face="bold")) + 
	theme(plot.title=element_text(size=rel(3), face="bold", vjust=1.5)) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(fill=NA, color="black", size=0.5), 
	      panel.background=element_blank(),
	      panel.margin.x=unit(0, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_text(color="black", size=rel(3)), 
		  axis.title.x=element_text(size=rel(2.5), face="bold", vjust=-0.4),  
		  axis.title.y=element_text(size=rel(2.5), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))  +
	theme(plot.margin=unit(c(1.5,2,1.15,1), "lines"))
)
dev.off()
# ----------------------------------------

# ----------------------------------------
# Some Plotting etc of NPP
# ----------------------------------------
summary(dat.ecosys)

ecosys$Site <- recode(ecosys$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(ecosys$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
ecosys$Site2 <- ecosys$Site
levels(ecosys$Site2) <- c("Demming (MN)", "Billy's (MN)", "UNDERC (WI)", "Minden (MI)", "Harvard (MA)", "Howland (ME)")


# models.use <- unique(dat.ecosys[,"Model.Order"])
# colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])
pdf(file.path(fig.dir, "NPP_Scales_SitesAll_0850-2010_Simple.pdf"), height=8.5, width=11)
print(
ggplot(data=ecosys[,])  + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=NPP, color=Model.Order), size=0.3, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=NPP, color=Model.Order), size=1.5, alpha=1) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
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
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1))
)
dev.off()

pdf(file.path(fig.dir, "NPP_Scales_Sites3_0850-2010_Simple.pdf"), height=9, width=16)
print(
ggplot(data=ecosys[ecosys$Site %in% c("PHA", "PUN", "PDL"),])  + facet_wrap(~Site2) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site %in% c("PHA", "PUN", "PDL"),], aes(x=Year, y=NPP, color=Model.Order), size=0.3, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site %in% c("PHA", "PUN", "PDL"),], aes(x=Year, y=NPP, color=Model.Order), size=2.5, alpha=1) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
    theme(plot.margin=unit(c(1,2,1,2), "lines")) +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(strip.text=element_text(size=rel(2), face="bold")) + 
	theme(legend.text=element_text(size=rel(1.5)), 
	      legend.title=element_text(size=rel(1.75)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1))
)
dev.off()



pdf(file.path(fig.dir, "NPP_Scales_PHA_0850-2010_Simple.pdf"), height=8.5, width=11)
print(
ggplot(data=ecosys[ecosys$Site=="PHA",])  + facet_wrap(~Site2) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=0.3, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=2.5, alpha=1) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
    theme(plot.margin=unit(c(1,2,1,2), "lines")) +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(strip.text=element_text(size=rel(2), face="bold")) + 
	theme(legend.text=element_text(size=rel(1.5)), 
	      legend.title=element_text(size=rel(1.75)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1))
)
dev.off()


pdf(file.path(fig.dir, "NPP_Scales_SitesIndividual_0850-2010_Simple.pdf"), height=8.5, width=11)
for(s in unique(ecosys$Site2)){
print(
ggplot(data= ecosys[ecosys$Site2==s,]) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site2==s,], aes(x=Year, y=NPP, color=Model.Order), size=0.5, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site2==s,], aes(x=Year, y=NPP, color=Model.Order), size=3, alpha=1) +
	geom_text(aes(label=s, x=1250, y=20), size=rel(10), face="bold") +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(limits=range(dat.ecosys$NPP, na.rm=T), expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
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
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1))
)
}
dev.off()
# ----------------------------------------






# ----------------------------------------
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# Using mean model NPP across sites since the GAMM response curves are for 
			#    the whole model & not site-specific are parameterized
			# -----------------------
			# Find the start year for the extent
			# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			npp <- mean(dat.ecosys[dat.ecosys$Model==m  & dat.ecosys$Resolution==r, "NPP"], na.rm=T)			
			# -----------------------
			
			# -----------------------
			# Relativizing everything in dat.ecosys to make it comparable to tree rings
			# -----------------------
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr"])/npp
			# -----------------------

			
			# -----------------------
			# Finding the percent change in NPP relative to the mean for that particular scale
			# -----------------------
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr"])/npp

			# Tacking on the simulated distributions so we can do ensemble CIs or robust comparisons
			sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,7:ncol(sim.terms)] <- (sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,7:ncol(sim.terms)])/npp
			# -----------------------

			# -----------------------
			# Relativizing the factor fits through times and weights as well
			# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
			# -----------------------
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.tair.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.tair"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.precipf.rel"] <- 1+ (wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.precipf"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r,"fit.CO2.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r,"fit.CO2"])/npp
			# -----------------------

		# }
	}
}
summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
ci.terms.graph[which(ci.terms.graph$mean.rel>0.9),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(0.9),"lwr.rel"] <- 0.9
ci.terms.graph[ci.terms.graph$upr.rel>(0.9),"upr.rel"] <- 0.9 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)


# Plot the relativized
# pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Rel_Site.pdf"), height=8.5, width=11)
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Resolution==resolution, ]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
# dev.off()

levels(ci.terms.graph$Effect) <- c("Temperature", "Precipitation", "CO2")

plot.tair <- ggplot(data=ci.terms.graph[ci.terms.graph$Resolution==resolution & ci.terms.graph$Effect=="Temperature", ]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name=expression(bold(paste("Temp (May-Sep, "^"o","C)"))), breaks=seq(9,19, by=3)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text=element_text(size=rel(2.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_text(color="black", size=rel(3)), 
		  axis.title.x=element_text(size=rel(2.5), face="bold", vjust=-0.4),  
		  axis.title.y=element_text(size=rel(2.5), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(4,0,1.15,1), "lines"))

plot.precip <- ggplot(data=ci.terms.graph[ci.terms.graph$Resolution==resolution & ci.terms.graph$Effect=="Precipitation", ]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name="Precip (May-Sep, mm)") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text=element_text(size=rel(2.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold", vjust=-1),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(4,0,1.95,-1), "lines"))

plot.co2 <- ggplot(data=ci.terms.graph[ci.terms.graph$Resolution==resolution & ci.terms.graph$Effect=="CO2", ]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name=expression(bold(paste("CO"["2"]," (pmm)"))), breaks=seq(280, 400, by=30)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	ggtitle("Cross-Site Sensitivity") +
	theme(plot.title=element_text(size=rel(3), face="bold", hjust=18, vjust=2)) + 
	guides(fill=guide_legend(title.position="top", ncol=3), 
	       color=guide_legend(title.position="top", ncol=3)) +
	theme(legend.title=element_text(size=rel(2), face="bold"),
	      legend.text=element_text(size=rel(1.5)),
	      legend.position=c(0.2, 0.18),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text=element_text(size=rel(2.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      # panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold", vjust=-1),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1.7,1,1.45,-1.1), "lines"))


pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Rel_Site_Presentation.pdf"), height=9, width=16)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3, widths=c(1.1,1,1))))
print(plot.tair    , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.precip, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
print(plot.co2, vp = viewport(layout.pos.row = 1, layout.pos.col=3))
dev.off()

# ---------------------------------
# Looking at the ensemble mean of Drivers of change and sensitivities through time
# ---------------------------------

# Merging wt.terms & dat.ecosys to get the relativized NPP lined up
wt.terms2 <- merge(wt.terms, dat.ecosys[,c("Model", "Model.Order", "Site", "Year", "Resolution", "Extent", "NPP", "NPP.rel")], all.x=T, all.y=F)
summary(wt.terms2)

# indices.wt2 <- wt.terms2$Site=="PHA" & wt.terms2$Resolution=="t.001" & wt.terms2$Extent=="0850-2010" 
# indices.wt2b <- wt.terms2$Site=="PHA" & wt.terms2$Resolution=="t.010" & wt.terms2$Extent=="0850-2010" 


summary(wt.terms2)
ensemble.base.wts   <- wt.terms2[,]
ensemble.base.terms <- ci.terms[,]

summary(ensemble.base.wts)
summary(ensemble.base.terms)

ensemble.wts1 <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Resolution, ensemble.base.wts$Year), FUN=mean)
names(ensemble.wts1)[1:3] <- c("Site", "Resolution", "Year") 
summary(ensemble.wts1)


ensemble.wts.lo <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Resolution, ensemble.base.wts$Year), FUN=quantile, 0.025)
names(ensemble.wts.lo) <- c("Site", "Resolution", "Year", paste0(names(ensemble.wts.lo[4:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Resolution, ensemble.base.wts$Year), FUN=quantile, 0.975)
names(ensemble.wts.hi) <- c("Site", "Resolution", "Year", paste0(names(ensemble.wts.hi[4:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

dim(ensemble.wts1); dim(ensemble.wts.lo); dim(ensemble.wts.hi)

ensemble.wts <- cbind(ensemble.wts1, ensemble.wts.lo[,2:ncol(ensemble.wts.lo)], ensemble.wts.hi[,2:ncol(ensemble.wts.hi)])
summary(ensemble.wts)

ensemble.wts2 <- stack(ensemble.wts[,c("fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel")])[,c(2,1)]
names(ensemble.wts2) <- c("Driver", "fit.rel")
ensemble.wts2$Driver <- as.factor(substr(ensemble.wts2$Driver, 5, nchar(paste(ensemble.wts2$Driver))-4 ))
levels(ensemble.wts2$Driver) <- c("CO2", "Precip", "Temp")
ensemble.wts2[,c("Site", "Resolution", "Year")] <- ensemble.wts[,c("Site", "Resolution", "Year")]
ensemble.wts2$fit.rel.lo <- stack(ensemble.wts[,c("fit.tair.rel.lo", "fit.precipf.rel.lo", "fit.CO2.rel.lo")])[,1]
ensemble.wts2$fit.rel.hi <- stack(ensemble.wts[,c("fit.tair.rel.hi", "fit.precipf.rel.hi", "fit.CO2.rel.hi")])[,1]
summary(ensemble.wts2)


ensemble.wts$Site <- recode(ensemble.wts$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(ensemble.wts$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")

ensemble.wts2$Site <- recode(ensemble.wts$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(ensemble.wts2$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")


pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_SitesAll_1850-2010_t.001.pdf"), width=11, height=7.5)
# png(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_AllSites_1800-2010_Annual.png"), width=1100, height=750)
print( 
ggplot() + facet_grid(Site~., scales="free_y") +
 	geom_ribbon(data= ensemble.wts[ensemble.wts$Resolution=="t.001",], aes(x=Year, ymin=NPP.rel.lo*100, ymax=NPP.rel.hi*100), alpha=0.5) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHO",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHO","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHO","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHO","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHA",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHA","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHA","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PHA","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PMB",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PMB","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PMB","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PMB","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PUN",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PUN","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PUN","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PUN","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PBL",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PBL","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PBL","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PBL","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PDL",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PDL","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PDL","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site=="PDL","weight.precipf"])), size=4) +
	scale_x_continuous(limits=c(1800,2005), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("% Mean NPP"))), expand=c(0,0)) +
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
	      panel.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(2.5))) + 
	theme(strip.text=element_text(size=rel(0.75), face="bold")) +
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
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
dev.off()

pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_SitesIndividual_1800-2010_t.001.pdf"), width=15, height=8)
for(s in unique(ensemble.wts$Site)){
plot.ensemble.npp <- ggplot() +
 	geom_ribbon(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site==s,], aes(x=Year, ymin=NPP.rel.lo*100, ymax=NPP.rel.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site==s,], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site==s,"weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site==s,"weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.001" & ensemble.wts$Site==s,"weight.precipf"])), size=6) +
	geom_hline(y=100, linetype="dashed") +
	geom_text(aes(label=s, x=1815, y=250), size=rel(15), face="bold") +
	scale_x_continuous(limits=c(1800,2005), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), limits=range(ensemble.wts[ensemble.wts$Year>=1800,c("NPP.rel.lo", "NPP.rel.hi")], na.rm=T)*100, expand=c(0,0)) +
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
	theme(axis.text.x=element_blank(),
		  axis.text.y=element_text(size=rel(2), color="black"), 
		  axis.title.x=element_blank(),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,1,-0.9,0), "lines"))


plot.ensemble.weights <- ggplot() +
 	geom_ribbon(data= ensemble.wts2[ensemble.wts$Resolution=="t.001" & ensemble.wts2$Site==s,], aes(x=Year, ymin=fit.rel.lo*100, ymax= fit.rel.hi*100, fill=Driver), alpha=0.35) +
 	geom_line(data= ensemble.wts2[ensemble.wts$Resolution=="t.001" & ensemble.wts2$Site==s,], aes(x=Year, y=fit.rel*100, color=Driver, linetype=Driver), size=1.5) +
	scale_x_continuous(limits=c(1800,2005), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), limits=range(ensemble.wts2[ensemble.wts2$Year>=1800,c("fit.rel.lo", "fit.rel.hi")], na.rm=T)*100, expand=c(0,0), breaks=seq(50, 200, by=50)) +
	scale_fill_manual(values=c("green","blue", "red")) +
	scale_color_manual(values=c("green4", "blue2", "red3")) +
	scale_linetype_manual(values=c("dashed", "solid", "longdash")) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1.5)), 
	      legend.title=element_text(size=rel(1.5)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines"),
	      legend.position=c(0.25, 0.83),
	      legend.direction="horizontal") + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_text(color="black", size=rel(2)), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(0,1,0.5,0), "lines"))
	      
	      
# 2-panel Figure!!
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
print(plot.ensemble.npp    , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.ensemble.weights, vp = viewport(layout.pos.row = 2, layout.pos.col=1))
}
dev.off()

# --------





pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_SitesIndividual_1800-2010_t.010.pdf"), width=15, height=8)
for(s in unique(ensemble.wts$Site)){
plot.ensemble.npp <- ggplot() +
 	geom_ribbon(data= ensemble.wts[ensemble.wts$Resolution=="t.010" & ensemble.wts$Site==s,], aes(x=Year, ymin=NPP.rel.lo*100, ymax=NPP.rel.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts[ensemble.wts$Resolution=="t.010" & ensemble.wts$Site==s,], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Resolution=="t.010" & ensemble.wts$Site==s,"weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.010" & ensemble.wts$Site==s,"weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Resolution=="t.010" & ensemble.wts$Site==s,"weight.precipf"])), size=8) +
	geom_hline(y=100, linetype="dashed") +
	geom_text(aes(label=s, x=1815, y=250), size=rel(15), face="bold") +
	scale_x_continuous(limits=c(1800,2005), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), limits=range(ensemble.wts[ensemble.wts$Year>=1800,c("NPP.rel.lo", "NPP.rel.hi")], na.rm=T)*100, expand=c(0,0)) +
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
	theme(axis.text.x=element_blank(),
		  axis.text.y=element_text(size=rel(2), color="black"), 
		  axis.title.x=element_blank(),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,1,-0.9,0), "lines"))


plot.ensemble.weights <- ggplot() +
 	geom_ribbon(data= ensemble.wts2[ensemble.wts$Resolution=="t.010" & ensemble.wts2$Site==s,], aes(x=Year, ymin=fit.rel.lo*100, ymax= fit.rel.hi*100, fill=Driver), alpha=0.35) +
 	geom_line(data= ensemble.wts2[ensemble.wts$Resolution=="t.010" & ensemble.wts2$Site==s,], aes(x=Year, y=fit.rel*100, color=Driver, linetype=Driver), size=1.5) +
	scale_x_continuous(limits=c(1800,2005), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), limits=range(ensemble.wts2[ensemble.wts2$Year>=1800,c("fit.rel.lo", "fit.rel.hi")], na.rm=T)*100, expand=c(0,0), breaks=seq(50, 200, by=50)) +
	scale_fill_manual(values=c("green","blue", "red")) +
	scale_color_manual(values=c("green4", "blue2", "red3")) +
	scale_linetype_manual(values=c("dashed", "solid", "longdash")) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1.5)), 
	      legend.title=element_text(size=rel(1.5)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines"),
	      legend.position=c(0.25, 0.83),
	      legend.direction="horizontal") + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_text(color="black", size=rel(2)), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),   
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(0,1,0.5,0), "lines"))
	      
	      
# 2-panel Figure!!
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
print(plot.ensemble.npp    , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.ensemble.weights, vp = viewport(layout.pos.row = 2, layout.pos.col=1))
}
dev.off()

# ----------------------------------------






