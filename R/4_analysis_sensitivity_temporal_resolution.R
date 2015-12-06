# ----------------------------------------
# Sensitivity & Scaling Analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 16 November 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
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
in.base <- "Data/gamms/"
out.dir <- "Data/analyses/analysis_temporal_resolution"
fig.dir <- "Figures/analyses/analysis_temporal_resolution"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]


load(file.path(in.base, "Sensitivity_All_TempRes/gamm_Models_NPP_Resolutions_ExtentFull.Rdata"))
summary(mod.out)

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

# ----------------------------------------

# ----------------------------------------
# Some Plotting etc of NPP
# ----------------------------------------
summary(dat.ecosys)

# models.use <- unique(dat.ecosys[,"Model.Order"])
# colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])
pdf(file.path(fig.dir, "NPP_Scales_SitesAll_0850-2010_Simple.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001",])  + facet_wrap(~Site) +
	geom_line(data=dat.ecosys[dat.ecosys$Resolution=="t.001",], aes(x=Year, y=NPP, color=Model.Order), size=0.25, alpha=0.3) +
	geom_line(data=dat.ecosys[dat.ecosys$Resolution=="t.010",], aes(x=Year, y=NPP, color=Model.Order), size=0.7, alpha=1) +
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

pdf(file.path(fig.dir, "NPP_Scales_Sites3_0850-2010_Simple.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Site %in% c("PHA", "PUN", "PBL"),])  + facet_wrap(~Site) +
	geom_line(data=dat.ecosys[dat.ecosys$Site %in% c("PHA", "PUN", "PBL") & dat.ecosys$Resolution=="t.001",], aes(x=Year, y=NPP, color=Model.Order), size=0.25, alpha=0.3) +
	geom_line(data=dat.ecosys[dat.ecosys$Site %in% c("PHA", "PUN", "PBL") & dat.ecosys$Resolution=="t.010",], aes(x=Year, y=NPP, color=Model.Order), size=0.7, alpha=1) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
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
for(s in unique(dat.ecosys$Site)){
print(
ggplot(data=dat.ecosys[dat.ecosys$Site==s,]) +
	geom_line(data=dat.ecosys[dat.ecosys$Site==s & dat.ecosys$Resolution=="t.001",], aes(x=Year, y=NPP, color=Model.Order), size=0.5, alpha=0.3) +
	geom_line(data=dat.ecosys[dat.ecosys$Site==s & dat.ecosys$Resolution=="t.010",], aes(x=Year, y=NPP, color=Model.Order), size=1, alpha=1) +
	geom_text(aes(label=s, x=1000, y=20), size=rel(15), face="bold") +
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
# Extract gam fit information
# ----------------------------------------
summary.stats <- data.frame(Model=unique(ecosys$Model))
summary.stats

# summary(mod.out$gamm.clm.bgc)
# summary(mod.site)
# summary(mod.comp)

for(m in unique(summary.stats$Model)){
	summary.stats[summary.stats$Model==m, "R2.t.001" ] <- summary(mod.out[[paste0("gamm.", m, "_t.001")]])$r.sq
	summary.stats[summary.stats$Model==m, "R2.t.010" ] <- summary(mod.out[[paste0("gamm.", m, "_t.010")]])$r.sq
	summary.stats[summary.stats$Model==m, "Dev.t.001"] <- summary(mod.out[[paste0("gamm.", m, "_t.001")]])$dev.expl
	summary.stats[summary.stats$Model==m, "Dev.t.010"] <- summary(mod.out[[paste0("gamm.", m, "_t.010")]])$dev.expl
}
summary.stats$dR2  <- summary.stats$R2.t.010  - summary.stats$R2.t.001
summary.stats$dDev <- summary.stats$Dev.t.010 - summary.stats$Dev.t.001
summary.stats
# mean(summary.stats$dSite); sd(summary.stats$dSite)

write.csv(summary.stats, file.path(out.dir, "Summary_ModelFits_TemporalRes.csv"), row.names=F)
# ----------------------------------------

# ==========================================================================
# **************************************************************************
# ==========================================================================
# Results!
# ==========================================================================
# **************************************************************************
# ==========================================================================


# ----------------------------------------
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Standardize responses to mean model NPP
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
			sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,8:ncol(sim.terms)] <- (sim.terms[sim.terms$Model==m & sim.terms$Resolution==r,8:ncol(sim.terms)])/npp
			# -----------------------

		# }
	}
}
summary(dat.ecosys)
summary(ci.terms)

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-0.75),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.75),"lwr.rel"] <- -0.75 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.75),"upr.rel"] <- -0.75 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.0),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(1.0),"lwr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$upr.rel>(1.0),"upr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)

# Getting Rid of values outside the range of those observed in the decadal record
ci.terms.graph2 <- ci.terms.graph
for(e in unique(ci.terms.graph2$Effect)){
	ci.terms.graph2 <- 	ci.terms.graph2[
			!ci.terms.graph2$Effect==e | ci.terms.graph2$Resolution=="t.010" | 
			(ci.terms.graph2$Effect==e & ci.terms.graph2$Resolution=="t.001" & 
			ci.terms.graph2$x >= min(ci.terms.graph2[ci.terms.graph2$Effect==e & ci.terms.graph2$Resolution=="t.010","x"]) & 
			ci.terms.graph2$x <= max(ci.terms.graph2[ci.terms.graph2$Effect==e & ci.terms.graph2$Resolution=="t.010","x"])),]
}
summary(ci.terms.graph2)

# Plot the relativized Effects
pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Comparison_Resolution.pdf"), height=8.5, width=11)
print(
ggplot(data=ci.terms.graph[,]) + facet_grid(Resolution~Effect, scale="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	ggtitle("Full Ranges") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
print(
ggplot(data=ci.terms.graph2[,]) + facet_grid(Resolution~Effect, scale="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	# ggtitle("Decadal Ranges") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	guides(col=guide_legend(nrow=2), fill=guide_legend(nrow=2)) +
	theme(legend.title=element_text(size=rel(1.5), face="bold"),
	      legend.text=element_text(size=rel(1.25)),
	      legend.position=c(0.5, 0.5),
	      legend.direction="horizontal",
	      legend.key=element_blank(),
	      legend.key.size=unit(1.25, "lines")) +
	theme(strip.text=element_text(size=rel(1.25), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      # panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
		  axis.text.y=element_text(color="black", size=rel(1.5)), 
		  axis.title.x=element_text(size=rel(1.5), face="bold"),  
		  axis.title.y=element_text(size=rel(1.5), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,1,1,1), "lines"))
)
dev.off()


levels(ci.terms.graph2$Effect)     <- c("Temperature", "Precipitation", "CO2")
levels(ci.terms.graph2$Resolution) <- c("Annual", "Decadal")

plot.tair <- ggplot(data=ci.terms.graph2[ci.terms.graph2$Effect=="Temperature", ]) + facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name="Temp (May - Sep, C)", breaks=c(13:17)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph2[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text.x=element_text(size=rel(2), face="bold"),
	      strip.text.y=element_blank()) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
		  axis.text.y=element_text(color="black", size=rel(1.5)), 
		  axis.title.x=element_text(size=rel(1.5), face="bold"),  
		  axis.title.y=element_text(size=rel(1.5), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,-1,1,1), "lines"))

plot.precip <- ggplot(data=ci.terms.graph2[ci.terms.graph2$Effect=="Precipitation", ]) + facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name="Precip (May - Sep, mm)", breaks=c(350, 450, 550)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph2[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(col=guide_legend(ncol=2), fill=guide_legend(ncol=2)) +
	theme(legend.title=element_text(size=rel(1.), face="bold"),
	      legend.text=element_text(size=rel(1.)),
	      legend.position=c(0.0, 0.4),
	      legend.direction="horizontal",
	      legend.key=element_blank(),
	      legend.key.size=unit(1.25, "lines")) +
	theme(strip.text.x=element_text(size=rel(2), face="bold"),
	      strip.text.y=element_blank()) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(1.5), face="bold"),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,-1,1,-1), "lines"))

plot.co2 <- ggplot(data=ci.terms.graph2[ci.terms.graph2$Effect=="CO2", ]) + facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name="CO2 (pmm)", breaks=c(275, 300, 325, 350, 375)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph2[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text=element_text(size=rel(1.6), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      # panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(1.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(1.5), face="bold"),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,1,1,-1), "lines"))

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Rel_Resolution_Presentation.pdf"), height=8.5, width=11)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,3, widths=c(1.25,1,1.25))))
print(plot.tair    , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.precip, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
print(plot.co2, vp = viewport(layout.pos.row = 1, layout.pos.col=3))
dev.off()

# ----------------------------------------


# ----------------------------------------
# Looking at distribution of NPP or climate senstivity with
# ----------------------------------------
