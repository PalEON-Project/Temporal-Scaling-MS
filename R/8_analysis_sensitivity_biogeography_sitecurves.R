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
out.dir <- "Data/analyses/analysis_biogeography"
fig.dir <- "Figures/analyses/analysis_biogeography"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))
ecosys <- ecosys[!ecosys$Model=="linkages",]


load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP_Site.Rdata"))
mod.baseline <- mod.out; rm(mod.out)
# names(mod.baseline)[7:length(mod.baseline)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
summary(mod.baseline)

load(file.path(in.base, "Sensitivity_Models_Site/gamm_models_NPP_Site.Rdata"))
mod.site <- mod.out; rm(mod.out)
summary(mod.site)

# ----------------------------------------


# ----------------------------------------
# Extract gam fit information
# ----------------------------------------
summary.stats <- data.frame(Model=unique(ecosys$Model))
summary.stats

# summary(mod.baseline$gamm.clm.bgc)
# summary(mod.site)
# summary(mod.comp)

for(m in unique(summary.stats$Model)){
	summary.stats[summary.stats$Model==m, "R2.baseline" ] <- summary(mod.baseline [[paste0("gamm.", m, ".baseline")]])$r.sq
	summary.stats[summary.stats$Model==m, "R2.site"     ] <- summary(mod.site     [[paste0("gamm.", m)]])$r.sq
}
summary.stats$dSite <- summary.stats$R2.site - summary.stats$R2.baseline
summary.stats
mean(summary.stats$dSite); sd(summary.stats$dSite)

write.csv(summary.stats, file.path(out.dir, "GAMM_ModelFits_Site_Curves.csv"), row.names=F)
# ----------------------------------------

# ==========================================================================
# **************************************************************************
# ==========================================================================
# Results!
#    Site curves had highest R2, although differences were modest
#    --> Slightly higher change in R2 with site curves in dynamic veg models
#
#    In general, static veg models had the highest R2 values
#    --> makes sense because there are no veg shifts muddying the picture over
#        long time periods (hypothesis to test)
# ==========================================================================
# **************************************************************************
# ==========================================================================


# ----------------------------------------
# Looking at changes in climate sensitivity with addition of site intercept
# Note: These will all be based on the "baseline" and "Site" intercept models
#       because this is where we saw the greatest change
# ----------------------------------------

# Binding the input data together
dat.ecosys <- cbind(mod.baseline$data, mod.baseline$ci.response[,c("mean", "lwr", "upr")], GAM="Baseline")
dat.ecosys <- rbind(dat.ecosys, cbind(mod.site$data[, !(names(mod.site$data) %in% c("Evergreen", "Grass"))], mod.site$ci.response[,c("mean", "lwr", "upr")], GAM="Site"))
summary(dat.ecosys)

# Binding the term sensitivity together;
# NOTE: Only using PHA because the sensitivity is the same across sites (there's only a site intercept)
ci.terms <- rbind(cbind(GAM="Baseline", mod.baseline$ci.terms[, ]), cbind(GAM="Site", mod.site$ci.terms[,]))
sim.terms <- rbind(cbind(GAM="Baseline", mod.baseline$sim.terms[, ]), cbind(GAM="Site", mod.site$sim.terms[, ]))
sim.terms$Effect <- as.factor(sim.terms$Effect)

summary(ci.terms)
summary(sim.terms[,1:10])

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

# Trimming out the part of the site-specific curves that lies outside of the range of observatiosn
ci.terms.graph2 <- ci.terms.graph
                                
for(s in unique(ci.terms.graph$Site)){
	tair     <- range(dat.ecosys[dat.ecosys$Site==s, "tair"   ], na.rm=T) - 273
	precipf  <- range(dat.ecosys[dat.ecosys$Site==s, "precipf"], na.rm=T) 
	co2      <- range(dat.ecosys[dat.ecosys$Site==s, "CO2"    ], na.rm=T) 

	# ----------------------------------
	# Creating a mask for ci.terms.graph
	# ----------------------------------
	# Tair
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="tair","line.min"] <- tair[1]
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="tair","line.max"] <- tair[2]

	# Precipf
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="precipf","line.min"] <- precipf[1]
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="precipf","line.max"] <- precipf[2]

	# CO2
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="CO2","line.min"] <- co2[1]
	ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$Effect=="CO2","line.max"] <- co2[2]
	# ----------------------------------

	# ----------------------------------
	# Getting rid of observations outside of the range
	# ----------------------------------
	ci.terms.graph2 <- ci.terms.graph2[!(ci.terms.graph2$Site==s) |  (ci.terms.graph2$Site==s & (
	                (ci.terms.graph2$Effect=="tair" & ci.terms.graph2$x >= tair[1] & ci.terms.graph2$x <= tair[2] ) |
	                (ci.terms.graph2$Effect=="precipf" & ci.terms.graph2$x >= precipf[1] & ci.terms.graph2$x <= precipf[2]) |
	                (ci.terms.graph2$Effect=="CO2" & ci.terms.graph2$x >= co2[1] & ci.terms.graph2$x <= co2[2])
	                )),]
	# ----------------------------------

}
ci.terms.graph$x.min    <- ifelse(ci.terms.graph$x < ci.terms.graph$line.min, ci.terms.graph$x, ci.terms.graph$line.min)
ci.terms.graph$x.max    <- ifelse(ci.terms.graph$x > ci.terms.graph$line.max, ci.terms.graph$x, ci.terms.graph$line.max)
ci.terms.graph$mask.min <- min(ci.terms.graph$lwr.rel, na.rm=T)
ci.terms.graph$mask.max <- max(ci.terms.graph$upr.rel, na.rm=T)

summary(ci.terms.graph)
summary(ci.terms.graph2)

# Plot the relativized Effects
pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Comparison_SiteCurves_byEffect.pdf"), height=8.5, width=11)
for(e in unique(ci.terms.graph$Effect)){
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Effect==e,]) + facet_grid(GAM~Site) +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	ggtitle(e) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Comparison_SiteCurves_bySite.pdf"), height=8.5, width=11)
for(s in unique(ci.terms.graph$Site)){
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Site==s,]) + facet_grid(GAM~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	ggtitle(s) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
}
dev.off()

# Making the overall model a pseudosite to be able to compare sites & the overall fit
ci.terms.graph3 <- ci.terms.graph
ci.terms.graph3[,"Site"] <- as.factor(ifelse(ci.terms.graph3$GAM=="Baseline", "Baseline", paste(ci.terms.graph3$Site)))
summary(ci.terms.graph3)

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Comparison_SiteCurves_byModel.pdf"), height=8.5, width=11)
for(m in unique(ci.terms.graph2$Model)){
print(
ggplot(data=ci.terms.graph3[ci.terms.graph3$Model==m,]) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(data=ci.terms.graph3[ci.terms.graph3$Model==m & ci.terms.graph3$Site=="Baseline",], 
	            aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100), alpha=0.5, fill="black") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Site), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Site), size=2) +
	geom_line(data=ci.terms.graph3[ci.terms.graph3$Model==m & ci.terms.graph3$Site=="Baseline",], 
	          aes(x=x, y=mean.rel*100), size=4, color="black") +
	# geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.min), linetype="dashed") +
	# geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.max), linetype="dashed") +
	ggtitle(m) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=c("black", "red", "chocolate1", "darkgoldenrod1", "green3", "blue", "darkorchid")) +
	scale_color_manual(values=c("black", "red", "chocolate1", "darkgoldenrod1", "green3", "blue", "darkorchid")) +
	theme_bw()
)
}
dev.off()
# ----------------------------------------

# ----------------------------------------
# Looking at the change in Effect Sensitivity from added site effect
# Note: Looking at absolute value of the effect to really be able to get at if we're getting MORE (negative) or LESS (positive) sensitive
# ----------------------------------------
dif.terms <- cbind(sim.terms[sim.terms$GAM=="Baseline",2:7], 
                   mean.rel = apply(abs(sim.terms[sim.terms$GAM=="Site", which(substr(names(sim.terms),1,1)=="X")]) - 
                                    abs(sim.terms[sim.terms$GAM=="Baseline"    , which(substr(names(sim.terms),1,1)=="X")]),
                                    1, FUN=mean),
                   lwr.rel  = apply(abs(sim.terms[sim.terms$GAM=="Site", which(substr(names(sim.terms),1,1)=="X")]) - 
                                    abs(sim.terms[sim.terms$GAM=="Baseline"    , which(substr(names(sim.terms),1,1)=="X")]),
                                    1, FUN=quantile, 0.025),
                   upr.rel  = apply(abs(sim.terms[sim.terms$GAM=="Site", which(substr(names(sim.terms),1,1)=="X")]) - 
                                    abs(sim.terms[sim.terms$GAM=="Baseline"    , which(substr(names(sim.terms),1,1)=="X")]),
                                    1, FUN=quantile, 0.975)
                              )
dif.terms$Effect <- recode(dif.terms$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(dif.terms$Effect) <- c("tair", "precipf", "CO2")
summary(dif.terms)

write.csv(dif.terms, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteCurves.csv"), row.names=F)

# ----------------------
# Visualizing the differences
# ----------------------
dif.graph <- dif.terms
dif.graph[dif.graph$mean<(-0.5),"mean.rel"] <- NA 
dif.graph[dif.graph$lwr<(-0.5),"lwr.rel"] <- -0.5 
dif.graph[dif.graph$upr<(-0.5),"upr.rel"] <- -0.5 
dif.graph[which(dif.graph$mean>0.5),"mean.rel"] <- NA 
dif.graph[dif.graph$lwr>(0.5),"lwr.rel"] <- 0.5 
dif.graph[dif.graph$upr>(0.5),"upr.rel"] <- 0.5 
dif.graph[dif.graph$Effect=="tair", "x"] <- dif.graph[dif.graph $Effect=="tair", "x"]-273.15
dif.graph$GAM <- as.factor("Difference")
summary(dif.graph)


for(s in unique(dif.graph$Site)){
	tair     <- range(dat.ecosys[dat.ecosys$Site==s, "tair"   ], na.rm=T) - 273
	precipf  <- range(dat.ecosys[dat.ecosys$Site==s, "precipf"], na.rm=T) 
	co2      <- range(dat.ecosys[dat.ecosys$Site==s, "CO2"    ], na.rm=T) 

	# ----------------------------------
	# Creating a mask for dif.graph
	# ----------------------------------
	# Tair
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="tair","line.min"] <- tair[1]
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="tair","line.max"] <- tair[2]

	# Precipf
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="precipf","line.min"] <- precipf[1]
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="precipf","line.max"] <- precipf[2]

	# CO2
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="CO2","line.min"] <- co2[1]
	dif.graph[dif.graph$Site==s & dif.graph$Effect=="CO2","line.max"] <- co2[2]
	# ----------------------------------
}
dif.graph$x.min    <- ifelse(dif.graph$x < dif.graph$line.min, dif.graph$x, dif.graph$line.min)
dif.graph$x.max    <- ifelse(dif.graph$x > dif.graph$line.max, dif.graph$x, dif.graph$line.max)
dif.graph$mask.min <- min(dif.graph$lwr.rel, na.rm=T)
dif.graph$mask.max <- max(dif.graph$upr.rel, na.rm=T)

summary(dif.graph)

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Difference_SiteCurves.pdf"), width=8.5, height=11)
print(
ggplot(data= dif.graph) + facet_grid(Site~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	geom_hline(yintercept=0, linetype="dashed") +
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	ggtitle("Change in Sensitivity (+ = more sensitive; - = less)") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="Change in Sensitivity (difference %NPP, + = more, - = less)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off() 

# Graphing only significant changes
dif.graph2 <- dif.graph
dif.graph2[!((dif.graph2$upr.rel<0 & dif.graph2$lwr.rel<0) | (dif.graph2$upr.rel>0 & dif.graph2$lwr.rel>0) ),c("mean.rel", "upr.rel", "lwr.rel")] <- NA 
summary(dif.graph2)

ggplot(data= dif.graph2[,]) + facet_grid(Site~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	geom_hline(yintercept=0, linetype="dashed") +
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="Change in Sensitivity (difference %NPP, + = more, - = less)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()

# ----------------------

# ----------------------
# Analyzing the differences
# ----------------------
summary(dif.terms)


dif.terms2 <- data.frame(Model  = rep(rep(unique(sim.terms$Model),3), 6), 
                         Site   = c(rep("PBL", length(unique(sim.terms$Model))*3),
                                    rep("PDL", length(unique(sim.terms$Model))*3),
                                    rep("PHA", length(unique(sim.terms$Model))*3),
                                    rep("PHO", length(unique(sim.terms$Model))*3),
                                    rep("PMB", length(unique(sim.terms$Model))*3),
                                    rep("PUN", length(unique(sim.terms$Model))*3)
                                    ),
                         Effect = rep(c(rep("tair", length(unique(sim.terms$Model))), 
                                    rep("precipf", length(unique(sim.terms$Model))), 
                                    rep("CO2", length(unique(sim.terms$Model)))), 6)
                                    )
summary(dif.terms2)

dif.terms.site <- data.frame(Site   = rep(c("PBL", "PDL", "PHA", "PHO", "PMB", "PUN"), 3),
                             Effect = c(rep("tair", 6), rep("precipf", 6), rep("CO2", 6))
                            )
dif.terms.site

dif.terms.mod <- data.frame(Model  = rep(unique(sim.terms$Model),3), 
                            Effect = c(rep("tair", length(unique(sim.terms$Model))), 
                                       rep("precipf", length(unique(sim.terms$Model))), 
                                       rep("CO2", length(unique(sim.terms$Model))))
                                      )
dif.terms.mod


for(m in unique(dif.terms2$Model)){
	for(e in unique(dif.terms2[dif.terms2$Model==m, "Effect"])){
		# Averaging across sites for models
		sim.dif <- as.matrix(abs(sim.terms[sim.terms$GAM=="Site" & sim.terms$Model==m & sim.terms$Effect==e, which(substr(names(sim.terms),1,1)=="X")]) - 
		                     abs(sim.terms[sim.terms$GAM=="Baseline" & sim.terms$Model==m & sim.terms$Effect==e, which(substr(names(sim.terms),1,1)=="X")]))

		dif.terms.mod[dif.terms.mod$Model==m & dif.terms.mod$Effect==e, "mean" ] <- mean(sim.dif)
		dif.terms.mod[dif.terms.mod$Model==m & dif.terms.mod$Effect==e, "stdev"] <- sd(sim.dif)
		dif.terms.mod[dif.terms.mod$Model==m & dif.terms.mod$Effect==e, "lwr"  ] <- quantile(sim.dif, 0.025)
		dif.terms.mod[dif.terms.mod$Model==m & dif.terms.mod$Effect==e, "upr"  ] <- quantile(sim.dif, 0.975)
		rm(sim.dif)	

		for(s in unique(dif.terms2[dif.terms2$Model==m, "Site"])){
			# Averaging across the range of x per site per model
			sim.dif <- as.matrix(abs(sim.terms[sim.terms$GAM=="Site" & sim.terms$Model==m & sim.terms$Effect==e & sim.terms$Site==s, which(substr(names(sim.terms),1,1)=="X")]) - 
		                     abs(sim.terms[sim.terms$GAM=="Baseline" & sim.terms$Model==m & sim.terms$Effect==e & sim.terms$Site==s, which(substr(names(sim.terms),1,1)=="X")]))
			dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e & dif.terms2$Site==s, "mean" ] <- mean(sim.dif)
			dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e & dif.terms2$Site==s, "stdev"] <- sd(sim.dif)
			dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e & dif.terms2$Site==s, "lwr"  ] <- quantile(sim.dif, 0.025)
			dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e & dif.terms2$Site==s, "upr"  ] <- quantile(sim.dif, 0.975)
			rm(sim.dif)	

			# Averaging across models at each site
			sim.dif <- as.matrix(abs(sim.terms[sim.terms$GAM=="Site" & sim.terms$Effect==e & sim.terms$Site==s, which(substr(names(sim.terms),1,1)=="X")]) - 
		                     abs(sim.terms[sim.terms$GAM=="Baseline" & sim.terms$Effect==e & sim.terms$Site==s, which(substr(names(sim.terms),1,1)=="X")]))
			dif.terms.site[dif.terms.site$Effect==e & dif.terms.site$Site==s, "mean" ] <- mean(sim.dif)
			dif.terms.site[dif.terms.site$Effect==e & dif.terms.site$Site==s, "stdev"] <- sd(sim.dif)
			dif.terms.site[dif.terms.site$Effect==e & dif.terms.site$Site==s, "lwr"  ] <- quantile(sim.dif, 0.025)
			dif.terms.site[dif.terms.site$Effect==e & dif.terms.site$Site==s, "upr"  ] <- quantile(sim.dif, 0.975)
			rm(sim.dif)	
		}
	}
}
dif.terms2$Effect <- recode(dif.terms2$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(dif.terms2$Effect) <- c("tair", "precipf", "CO2")
summary(dif.terms2)

dif.terms.site$Effect <- recode(dif.terms.site$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(dif.terms.site$Effect) <- c("tair", "precipf", "CO2")
dif.terms.site

dif.terms.mod$Effect <- recode(dif.terms.mod$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(dif.terms.mod$Effect) <- c("tair", "precipf", "CO2")
dif.terms.mod


write.csv(dif.terms2, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteCurves_Summary.csv"), row.names=F)
write.csv(dif.terms.site, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteCurves_Summary_Sites.csv"), row.names=F)
write.csv(dif.terms.mod, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteCurves_Summary_Models.csv"), row.names=F)


pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Difference_SiteCurves_2.pdf"), height=8.5, width=11)
print(
ggplot(data=dif.terms2) + facet_grid(Site~Effect) +
	geom_point(aes(x=Model, y=mean, color=Model), size=5) +
	geom_errorbar(aes(x=Model, ymin=lwr, ymax=upr, color=Model), size=1.5, width=0) +
	geom_hline(aes(yinterceipt=0), linetype="dashed") +
	# geom_errorbar(aes(x=Model, ymin=mean-stdev, ymax=mean+stdev, color=Model), size=1.5, width=0) +
 	scale_color_manual(values=colors.use) +
	theme_bw()
)
print(
ggplot(data=dif.terms.mod) + facet_grid(.~Effect) +
	geom_point(aes(x=Model, y=mean, color=Model), size=5) +
	geom_errorbar(aes(x=Model, ymin=lwr, ymax=upr, color=Model), size=1.5, width=0) +
	geom_hline(aes(yinterceipt=0), linetype="dashed") +
	# geom_errorbar(aes(x=Model, ymin=mean-stdev, ymax=mean+stdev, color=Model), size=1.5, width=0) +
 	scale_color_manual(values=colors.use) +
	theme_bw()
)
print(
ggplot(data=dif.terms.site) + facet_grid(.~Effect) +
	geom_point(aes(x=Site, y=mean, color=Site), size=5) +
	geom_errorbar(aes(x=Site, ymin=lwr, ymax=upr, color=Site), size=1.5, width=0) +
	geom_hline(aes(yinterceipt=0), linetype="dashed") +
	# geom_errorbar(aes(x=Model, ymin=mean-stdev, ymax=mean+stdev, color=Model), size=1.5, width=0) +
 	# scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off() 
# ----------------------
# ----------------------------------------
