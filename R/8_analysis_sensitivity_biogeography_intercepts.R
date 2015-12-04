# ----------------------------------------
# Sensitivity & Scaling Analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 13 November 2015
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


load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP.Rdata"))
mod.baseline <- mod.out; rm(mod.out)
# names(mod.baseline)[7:length(mod.baseline)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
summary(mod.baseline)

load(file.path(in.base, "Sensitivity_Models_Baseline/gamm_models_baseline_NPP_Site.Rdata"))
mod.site <- mod.out; rm(mod.out)
# names(mod.site)[7:length(mod.site)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
summary(mod.site)

load(file.path(in.base, "Sensitivity_Models_Comp/gamm_models_NPP_Comp.Rdata"))
mod.comp <- mod.out; rm(mod.out)
summary(mod.comp)

load(file.path(in.base, "Sensitivity_Models_Comp/gamm_models_NPP_Comp_Site.Rdata"))
mod.comp.site <- mod.out; rm(mod.out)
summary(mod.comp.site)
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
	summary.stats[summary.stats$Model==m, "R2.site"     ] <- summary(mod.site     [[paste0("gamm.", m, ".baseline")]])$r.sq
  if(!is.null(mod.comp     [[paste0("gamm.", m)]])){
	summary.stats[summary.stats$Model==m, "R2.comp"     ] <- summary(mod.comp     [[paste0("gamm.", m)]])$r.sq
	summary.stats[summary.stats$Model==m, "R2.comp.site"] <- summary(mod.comp.site[[paste0("gamm.", m)]])$r.sq
  }
}
summary.stats$dSite <- summary.stats$R2.site - summary.stats$R2.baseline
summary.stats$dComp <- summary.stats$R2.comp - summary.stats$R2.baseline
summary.stats$dSite.Comp  <- summary.stats$R2.site - summary.stats$R2.comp
summary.stats$dSite.Comp2 <- summary.stats$R2.comp.site - summary.stats$R2.site
summary.stats

write.csv(summary.stats, file.path(out.dir, "GAMM_ModelFits_Composition_Intercepts.csv"), row.names=F)
# ----------------------------------------

# ==========================================================================
# **************************************************************************
# ==========================================================================
# Results!
#    Adding Comp doesn't always make the fit better by much
#     --> Least difference between site & comp for dynamic veg models (for 
#         static veg models, site-intercept R2 is considerably higher)
#    R2 improvements from adding composition on top of site intercept are
#    small
#
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
dat.ecosys <- rbind(dat.ecosys, cbind(mod.site$data, mod.site$ci.response[,c("mean", "lwr", "upr")], GAM="Site"))
summary(dat.ecosys)

# Binding the term sensitivity together;
# NOTE: Only using PHA because the sensitivity is the same across sites (there's only a site intercept)
ci.terms <- rbind(cbind(GAM="Baseline", mod.baseline$ci.terms[mod.baseline$ci.terms$Site=="PHA", ]), cbind(GAM="Site", mod.site$ci.terms[mod.site$ci.terms$Site=="PHA", ]))
sim.terms <- rbind(cbind(GAM="Baseline", mod.baseline$sim.terms[mod.baseline$sim.terms$Site=="PHA", ]), cbind(GAM="Site", mod.site$sim.terms[mod.site$sim.terms$Site=="PHA", ]))
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


# Plot the relativized
pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Comparison_SiteIntercept.pdf"), height=8.5, width=11)
print(
ggplot(data=ci.terms.graph) + facet_grid(GAM~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off()

# ----------------------------------------

# ----------------------------------------
# Looking at the change in Effect Sensitivity from added site effect
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

write.csv(dif.terms, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteIntercept.csv"), row.names=F)

# ----------------------
# Visualizing the differences
# ----------------------
dif.graph <- dif.terms
dif.graph[dif.graph$mean<(-0.75),"mean.rel"] <- NA 
dif.graph[dif.graph$lwr<(-0.75),"lwr.rel"] <- -0.75 
dif.graph[dif.graph$upr<(-0.75),"upr.rel"] <- -0.75 
dif.graph[which(dif.graph$mean>1.0),"mean.rel"] <- NA 
dif.graph[dif.graph$lwr>(1.0),"lwr.rel"] <- 1.0 
dif.graph[dif.graph$upr>(1.0),"upr.rel"] <- 1.0 
dif.graph[dif.graph$Effect=="tair", "x"] <- dif.graph[dif.graph $Effect=="tair", "x"]-273.15
dif.graph$GAM <- as.factor("Difference")
summary(dif.graph)

ggplot(data= dif.graph) + facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()

# # Graphing only significant changes
# dif.graph2 <- dif.graph
# dif.graph2[!((dif.graph2$upr.rel<0 & dif.graph2$lwr.rel<0) | (dif.graph2$upr.rel>0 & dif.graph2$lwr.rel>0) ),c("mean.rel", "upr.rel", "lwr.rel")] <- NA 
# summary(dif.graph2)

# ggplot(data= dif.graph2[,]) + facet_wrap(~Effect, scales="free_x") +
	# geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	# geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	# geom_hline(yintercept=0, linetype="dashed") +
	# scale_x_continuous(expand=c(0,0)) +
	# scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,1)) +
	# scale_fill_manual(values=colors.use) +
	# scale_color_manual(values=colors.use) +
	# theme_bw()

ci.terms.graph2 <- merge(ci.terms.graph, dif.graph, all.x=T, all.y=T)
summary(ci.terms.graph2)

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Difference_SiteIntercept.pdf"), height=8.5, width=11)
print(
ggplot(data=ci.terms.graph2) + facet_grid(GAM~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off() 
# ----------------------

# ----------------------
# Analyzing the differences
# ----------------------
summary(dif.terms)


dif.terms2 <- data.frame(Model  = rep(unique(sim.terms$Model),3), 
                         Effect = c(rep("tair", length(unique(sim.terms$Model))), 
                                    rep("precipf", length(unique(sim.terms$Model))), 
                                    rep("CO2", length(unique(sim.terms$Model))))
                                    )
summary(dif.terms2)

for(m in unique(dif.terms2$Model)){
	for(e in unique(dif.terms2[dif.terms2$Model==m, "Effect"])){
		sim.dif <- as.matrix(sim.terms[sim.terms$GAM=="Baseline" & sim.terms$Model==m & sim.terms$Effect==e, which(substr(names(sim.terms),1,1)=="X")] - 
		                     sim.terms[sim.terms$GAM=="Site" & sim.terms$Model==m & sim.terms$Effect==e, which(substr(names(sim.terms),1,1)=="X")])
		dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e, "mean" ] <- mean(sim.dif)
		dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e, "stdev"] <- sd(sim.dif)
		dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e, "lwr"  ] <- quantile(sim.dif, 0.025)
		dif.terms2[dif.terms2$Model==m & dif.terms2$Effect==e, "upr"  ] <- quantile(sim.dif, 0.975)
	}
}
dif.terms2$Effect <- recode(dif.terms2$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(dif.terms2$Effect) <- c("tair", "precipf", "CO2")
summary(dif.terms2)
dif.terms2

write.csv(dif.terms2, file.path(out.dir, "GAMM_Sensitivity_Difference_SiteIntercept_Summary.csv"), row.names=F)


pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Difference_SiteIntercept_2.pdf"), height=8.5, width=11)
print(
ggplot(data=ci.terms.graph2) + facet_grid(GAM~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
)
print(
ggplot(data=dif.terms2) + facet_grid(.~Effect) +
	geom_point(aes(x=Model, y=mean, color=Model), size=5) +
	geom_errorbar(aes(x=Model, ymin=lwr, ymax=upr, color=Model), size=1.5, width=0) +
	# geom_errorbar(aes(x=Model, ymin=mean-stdev, ymax=mean+stdev, color=Model), size=1.5, width=0) +
 	scale_color_manual(values=colors.use) +
	theme_bw()
)
dev.off() 
# ----------------------
# ----------------------------------------
