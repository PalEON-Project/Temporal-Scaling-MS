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
# library(R2jags)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
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


load(file.path(in.base, "Sensitivity_All_TempRes/gamm_Models_NPP_Resolutions_ExtentFull.Rdata"))
mod.baseline <- mod.out; rm(mod.out)
# names(mod.baseline)[7:length(mod.baseline)] <- paste("gamm", c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "lpj.guess", "lpj.wsl", "sibcasa"), sep=".")
summary(mod.baseline)

load(file.path(in.base, "Sensitivity_Models_Site/t.001/gamm_models_NPP_SiteCurves.Rdata"))
mod.site <- mod.out; rm(mod.out)
summary(mod.site)



# Binding the input data together
dat.ecosys <- cbind(mod.baseline$data[,!(names(mod.site$data) %in% c("Evergreen", "Grass"))], mod.baseline$ci.response[,c("mean", "lwr", "upr")], GAM="Baseline")
dat.ecosys <- rbind(dat.ecosys, cbind(mod.site$data[, !(names(mod.site$data) %in% c("Evergreen", "Grass"))], mod.site$ci.response[,c("mean", "lwr", "upr")], GAM="Site"))
dat.ecosys <- merge(dat.ecosys, ecosys[,c("Model", "Site", "Year", "Resolution", "Evergreen", "Deciduous", "Grass", "SoilMoist", "LAI", "AGB")])
summary(dat.ecosys)


# Binding the term sensitivity together;
ci.terms <- rbind(cbind(GAM="Baseline", mod.baseline$ci.terms[, ]), cbind(GAM="Site", mod.site$ci.terms[,]))
wt.terms   <- rbind(cbind(GAM="Baseline", mod.baseline$weights), cbind(GAM="Site", mod.baseline$weights))
sim.terms <- rbind(cbind(GAM="Baseline", mod.baseline$sim.terms[, ]), cbind(GAM="Site", mod.site$sim.terms[, ]))
sim.terms$Effect <- as.factor(sim.terms$Effect)

dat.ecosys$Site <- recode(dat.ecosys$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(dat.ecosys$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
dat.ecosys$Site2 <- dat.ecosys$Site
levels(dat.ecosys$Site2) <- c("Demming (MN)", "Billy's (MN)", "UNDERC (WI)", "Minden (MI)", "Harvard (MA)", "Howland (ME)")

ci.terms$Site <- recode(ci.terms$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(ci.terms$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
ci.terms$Site2 <- ci.terms$Site
levels(ci.terms$Site2) <- c("Demming (MN)", "Billy's (MN)", "UNDERC (WI)", "Minden (MI)", "Harvard (MA)", "Howland (ME)")

wt.terms$Site <- recode(wt.terms$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(wt.terms$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
wt.terms$Site2 <- wt.terms$Site
levels(wt.terms$Site2) <- c("Demming (MN)", "Billy's (MN)", "UNDERC (WI)", "Minden (MI)", "Harvard (MA)", "Howland (ME)")

sim.terms$Site <- recode(sim.terms$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(sim.terms$Site) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
sim.terms$Site2 <- sim.terms$Site
levels(sim.terms$Site2) <- c("Demming (MN)", "Billy's (MN)", "UNDERC (WI)", "Minden (MI)", "Harvard (MA)", "Howland (ME)")

dat.ecosys <- dat.ecosys[dat.ecosys$Resolution=="t.001", ]
ci.terms   <- ci.terms  [ci.terms  $Resolution=="t.001", ]
wt.terms   <- wt.terms  [wt.terms  $Resolution=="t.001", ]
sim.terms  <- sim.terms [sim.terms $Resolution=="t.001", ]

summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)
summary(sim.terms[,1:10])



models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])
# ----------------------------------------

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
summary(sim.terms[,1:10])
# ----------------------------------------

# ----------------------------------------
# Splitting Effects & Merging in possible predictors
# ----------------------------------------
# Creatinga  single sim.terms that stores the overall model sensitivity 
# (no site curves) at a single site (since they're all the same)
sim.terms.baseline <- sim.terms[sim.terms$Site=="PHA" & sim.terms$GAM=="Baseline",]
summary(sim.terms.baseline[,1:10])

# --------------------------
# Removing values outside those observed at each site to allow for valid comparisons
# --------------------------
summary(dat.ecosys)
sim.terms.site <- sim.terms
for(s in unique(sim.terms.site$Site)){
	print(paste0(" -------- ", s, " -------- "))
	for(e in unique(sim.terms.site$Effect)){
		print(paste0("     ", e))
		effect.min <- min(dat.ecosys[dat.ecosys$Site==s, e])
		effect.max <- max(dat.ecosys[dat.ecosys$Site==s, e])
		print(paste0("     -- min : ", effect.min))
		print(paste0("     -- max : ", effect.max))
		sim.terms.site <- sim.terms.site[!sim.terms.site$Site==s | !sim.terms.site$Effect==e | 
			(sim.terms.site$Site==s & sim.terms.site$Effect==e & 
			 sim.terms.site$x>=effect.min & sim.terms.site$x<=effect.max)
			,]
	}
}
summary(sim.terms.site[,1:10])
dim(sim.terms.site); dim(sim.terms)
# --------------------------


# --------------------------
# Aggregating to get a single statistic per simulation per effect to simplify analyses
#   statistics: mean, min/max, range
# --------------------------
# Mean
agg.temp <- aggregate(sim.terms.site[,which(substr(names(sim.terms.site),1,1)=="X")], 
                          by=list(sim.terms.site$GAM, 
                                  sim.terms.site$Model, 
                                  sim.terms.site$Site, 
                                  sim.terms.site$Effect),
                          FUN=mean)
names(agg.temp)[1:4] <- c("GAM", "Model", "Site", "Effect")
summary(agg.temp[,1:10])

sim.agg <- stack(agg.temp[,which(substr(names(agg.temp),1,1)=="X")])
names(sim.agg) <- c("Sens.mean", "Sim")
sim.agg[,c("GAM", "Model", "Site", "Effect")] <- agg.temp[,c("GAM", "Model", "Site", "Effect")]
summary(sim.agg)

# Range
agg.temp <- aggregate(sim.terms.site[,which(substr(names(sim.terms.site),1,1)=="X")], 
                          by=list(sim.terms.site$GAM, 
                                  sim.terms.site$Model, 
                                  sim.terms.site$Site, 
                                  sim.terms.site$Effect),
                          FUN=function(x) c(range=max(x)-min(x)))
summary(agg.temp[,1:10])
sim.agg$Sens.range <- stack(agg.temp[,which(substr(names(agg.temp),1,1)=="X")])[,1]

summary(sim.agg)
rm(agg.temp)


# Aggregating the Ecosystem properties to single values for comparison
factors <- c("NPP", "NPP.rel", "tair", "precipf", "CO2", "Evergreen", "Deciduous", "Grass", "SoilMoist", "AGB", "LAI")
ecosys.agg <- aggregate(dat.ecosys[,factors], 
					    by=list(dat.ecosys$Model, 
                                dat.ecosys$Site),
                        FUN=mean)
names(ecosys.agg) <- c("Model", "Site", paste0(factors, ".mean"))
summary(ecosys.agg)

ecosys.agg[, paste0(factors, ".range")] <- aggregate(dat.ecosys[,factors], by=list(dat.ecosys$Model, dat.ecosys$Site), FUN=function(x) c(range=(max(x) - min(x))))[,3:(length(factors)+2)]
ecosys.agg[, paste0(factors, ".ci.lo")] <- aggregate(dat.ecosys[,factors], by=list(dat.ecosys$Model, dat.ecosys$Site), FUN=quantile, 0.025, na.rm=T)[,3:(length(factors)+2)]
ecosys.agg[, paste0(factors, ".ci.hi")] <- aggregate(dat.ecosys[,factors], by=list(dat.ecosys$Model, dat.ecosys$Site), FUN=quantile, 0.975, na.rm=T)[,3:(length(factors)+2)]
summary(ecosys.agg)


pfts <- as.factor(c("Evergreen", "Deciduous", "Grass"))
for(i in 1:nrow(ecosys.agg)){
	if(!is.na(max(ecosys.agg[i,paste0(pfts, ".mean")]))){
	ecosys.agg[i,"PFT.dominant"] <- pfts[which(ecosys.agg[i,paste0(pfts, ".mean")]==max(ecosys.agg[i,paste0(pfts, ".mean")]))]
	}
}
summary(ecosys.agg)

# Merging ecosys.agg into sim.agg
summary(sim.agg)
dim(sim.agg)

sim.agg <- merge(sim.agg, ecosys.agg, all.x=T, all.y=F)
summary(sim.agg)
dim(sim.agg)


# Re-ordering so we can plot sites W to E
sim.agg$Site.Order <- recode(sim.agg$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
levels(sim.agg$Site.Order) <- c("PDL", "PBL", "PUN", "PMB", "PHA", "PHO")
sim.agg$Effect <- recode(sim.agg$Effect, "'tair'='1'; 'precipf'='2'; 'CO2'='3'")
levels(sim.agg$Effect) <- c("Temp", "Precip", "CO2")


veg.stat <- c("clm.bgc", "clm.cn", "jules.stat", "sibcasa")
P.Collatz <- c("clm.bgc", "clm.cn", "lpj.guess", "lpj.wsl", "jules.stat", "jules.triffid", "sibcasa")
P.Farquar <- c("ed2", "ed2.lu")
met.all <- c("ed2", "ed2.lu", "jules.stat", "jules.triffid", "sibcasa", "clm.cn", "clm.bgc")

N.lim <- c("lpj.guess", "clm.cn", "clm.bgc", "jules.stat", "jules.triffid")

sim.agg$VegScheme <- as.factor(ifelse(sim.agg$Model %in% veg.stat, "Prescribed", "Emergent"))
sim.agg$Nlim <- as.factor(ifelse(sim.agg$Model %in% N.lim, "Yes", "No"))
sim.agg$nMet <- as.factor(ifelse(sim.agg$Model %in% met.all, 8, ifelse(ci.terms2$Model=="lpj.wsl", 5, 4)))
summary(sim.agg)



summary(ecosys.agg)


# Aggregating the met
summary(ecosys.agg)
met.agg <- stack(ecosys.agg[,paste0(factors, ".mean")])[,c(2,1)]
names(met.agg) <- c("Factor", "Mean")
met.agg[,c("Model", "Site")] <- ecosys.agg[,c("Model", "Site")]
met.agg$Factor <- as.factor(substr(met.agg$Factor, 1, nchar(paste(met.agg$Factor))-5))
met.agg$ci.lo <- stack(ecosys.agg[,paste0(factors, ".ci.lo")])[,1]
met.agg$ci.hi <- stack(ecosys.agg[,paste0(factors, ".ci.hi")])[,1]
# met.agg$Factor <- recode(met.agg$Factor, "'tair'=''; 'precipf'='2'; 'CO2'='3'; 'SoilMoist'='4'")
# levels(met.agg$Factor) <- c("Temp", "Precip", "CO2", "SoilMoist")
summary(met.agg)


# --------------------------

# --------------------------
# Graphing some of the summary statistics
# --------------------------
summary(sim.agg)

pdf(file.path(fig.dir, "SiteMean_Met.pdf"))
ggplot(data= met.agg[met.agg$Model=="ed2" & (met.agg$Factor %in% c("tair", "precipf")), ]) + facet_wrap(~Factor, scales="free_y") +
	geom_errorbar(aes(x=Site, ymin=ci.lo, ymax=ci.hi, color=Site), width=0, size=2) +
	geom_point(aes(x=Site, y=Mean, color=Site), size=5) +
	ggtitle("Mean Temp & Precip") +
	theme_bw()
ggplot(data= met.agg[(met.agg$Factor == "SoilMoist"), ]) + facet_wrap(~Model) +
	geom_errorbar(aes(x=Site, ymin=ci.lo, ymax=ci.hi, color=Site), width=0, size=2) +
	geom_point(aes(x=Site, y=Mean, color=Site), size=5) +
	ggtitle("Mean Soil Moisture") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "SiteMean_PFT"))
ggplot(data= met.agg[!met.agg$Model=="sibcasa" & (met.agg$Factor %in% c("Evergreen", "Deciduous", "Grass")), ]) + facet_wrap(~Model) +
	geom_errorbar(aes(x=Site, ymin=ci.lo, ymax=ci.hi, color=Factor), width=0, size=0.8) +
	geom_point(aes(x=Site, y=Mean, color=Factor, shape=Factor), size=5) +
	scale_y_continuous(name="Fraction PFT") +
	scale_shape_manual(values=c(1, 8, 5)) +
	ggtitle("Mean Fraction Deciduous") +
	theme_bw()
dev.off()

# ggplot(data=sim.agg[sim.agg$GAM=="Site",]) + facet_grid(Model~Effect, scales="free_y") +
	# geom_violin(aes(x=Site, y=Sens.range, fill=Site, color=Site)) +
	# # scale_color_manual(values=paste(colors.use)) +
	# # scale_fill_manual(values=paste(colors.use)) +
	# theme_bw()

	
# for(e in unique(sim.agg$Effect)){
# print(
# ggplot(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e,]) + facet_wrap(~Site) +
	# geom_violin(aes(x=Model, y=Sens.range, fill=Model, color=Model)) +
	# scale_color_manual(values=paste(colors.use)) +
	# scale_fill_manual(values=paste(colors.use)) +
	# ggtitle(paste0(e, " Sensitivity - Range")) +
	# theme_bw()
# )
# }

# for(m in unique(sim.agg$Model)){
# print(
# ggplot(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Model==m,]) + facet_wrap(~Effect) +
	# geom_violin(aes(x=Site, y=Sens.range, fill=Site, color=Site)) +
	# # scale_color_manual(values=paste(colors.use)) +
	# # scale_fill_manual(values=paste(colors.use)) +
	# ggtitle(paste0(m, " Sensitivity - Range")) +
	# theme_bw()
# )
# }


# for(e in unique(sim.agg$Effect)){
# print(
# ggplot(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e,]) + facet_wrap(~Site) +
	# geom_violin(aes(x=Model, y=Mean, fill=Model, color=Model)) +
	# scale_color_manual(values=paste(colors.use)) +
	# scale_fill_manual(values=paste(colors.use)) +
	# ggtitle(paste0(e, " Sensitivity - Mean")) +
	# theme_bw()
# )
# }

# --------------------------


# --------------------------
# Doing some quick & dirty analyses to see how different factors influence sensitivity
# --------------------------
library(nlme)
library(MuMIn)

sens.co2 <- sim.agg[sim.agg$Effect=="CO2",]
summary(sens.co2)


for(e in unique(sim.agg$Effect)){
print(
ggplot(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e,]) + facet_wrap(~Site) +
	geom_violin(aes(x=Model, y=Sens.range, fill=Model, color=Model)) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	ggtitle(paste0(e, " Sensitivity - Range")) +
	theme_bw()
)
}

pdf(file.path(fig.dir, "EffectRanges_byModel_bySite.pdf"))
for(e in unique(sim.agg$Effect)){
print(
ggplot(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e,]) + facet_wrap(~Model) +
	geom_violin(aes(x=Site, y=Sens.range, fill=Site, color=Site)) +
	# scale_color_manual(values=paste(colors.use)) +
	# scale_fill_manual(values=paste(colors.use)) +
	ggtitle(paste0(e, " Sensitivity - Range")) +
	theme_bw()
)
}
dev.off()
# -------------------------------------------------------
# The Quick & Dirty Way -- LM on site-level sensitivity that ignores Model Effects 
# -------------------------------------------------------
mod.comp <- lm(Sens.range ~ (Evergreen.mean + Grass.mean + Deciduous.mean)*(Effect-1) , data=sim.agg[,])
summary(mod.comp)
anova(mod.comp)

mod.evg   <- lm(Sens.range ~ Evergreen.mean*(Effect-1) - Evergreen.mean, data=sim.agg[,])
mod.decid <- lm(Sens.range ~ Deciduous.mean*(Effect-1) - Deciduous.mean, data=sim.agg[,])
mod.grass <- lm(Sens.range ~ Grass.mean*(Effect-1) - Grass.mean, data=sim.agg[,])
summary(mod.evg)
summary(mod.decid)
summary(mod.grass)


mod.moist <- lm(Sens.range ~ (SoilMoist.mean)*(Effect-1), data=sim.agg[,])
summary(mod.moist)
anova(mod.moist)

tair.evg <- lm(Sens.range ~ Evergreen.mean, data=sim.agg[sim.agg$Effect=="Temp",])
tair.decid <- lm(Sens.range ~ Deciduous.mean, data=sim.agg[sim.agg$Effect=="Temp",])
tair.grass <- lm(Sens.range ~ Grass.mean, data=sim.agg[sim.agg$Effect=="Temp",])
summary(tair.evg)
summary(tair.decid)
summary(tair.grass)

precip.evg <- lm(Sens.range ~ Evergreen.mean, data=sim.agg[sim.agg$Effect=="Precip",])
precip.decid <- lm(Sens.range ~ Deciduous.mean, data=sim.agg[sim.agg$Effect=="Precip",])
precip.grass <- lm(Sens.range ~ Grass.mean, data=sim.agg[sim.agg$Effect=="Precip",])
summary(precip.evg)
summary(precip.decid)
summary(precip.grass)


co2.comp <- lm(Sens.range ~ Evergreen.mean + Grass.mean, data=sim.agg[sim.agg$Effect=="CO2",])
tair.comp <- lm(Sens.range ~ Evergreen.mean + Grass.mean, data=sim.agg[sim.agg$Effect=="Temp",])
precipf.comp <- lm(Sens.range ~ Evergreen.mean + Grass.mean, data=sim.agg[sim.agg$Effect=="Precip",])
summary(co2.comp)
summary(tair.comp)
summary(precipf.comp)
anova(co2.comp)


# ------------------------------------------------- #
# Results summary (Annual Resolution)               #
# ------------------------------------------------- #
# 1) composition correlates greater with rough      #
#    climate sensitivity than soil moisture (which  #
#    is affect by comp (R2 0.78 vs 0.70)            #
# 2) All PFTs have similar correlations across clim.#
#    effects (R2=0.74-0.77)                         #
#    -- Evergreen associated with INCREASED sens.   #
#       across the board                            #
#    -- Deciduous associated with DECREASED sens.   #
#       for all clim. drivers                       #
#    -- Grass associated with decreased sens. for   # 
#       precip & has no effect on CO2 sensitivity   #
# 3) raw composition has the about the same corr.   #
#    across factors (CO2 R2=0.13, precip=0.15,      #
#    tair0.10)                                      #
# ------------------------------------------------- #
# -------------------------------------------------------

# -------------------------------------------------------
# Getting a bit more complicated: Looking at composition effects given the overal Model Trend
# -------------------------------------------------------
# ----------------
# Explaining differences among sites within a model
# ----------------
sens.moist.mod <- lme(Sens.range ~ (Effect-1)*SoilMoist.mean, random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.mod); 
anova(sens.moist.mod)
r.squaredGLMM(sens.moist.mod)


sens.comp.mod1 <- lme(Sens.range ~ (Effect-1)*Evergreen.mean*Grass.mean*Deciduous.mean, random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.mod); 
anova(sens.comp.mod1)
r.squaredGLMM(sens.comp.mod1)

sens.comp.mod2 <- lme(Sens.range ~ (Effect-1)*(Evergreen.mean + Grass.mean + Deciduous.mean), random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.mod2); 
anova(sens.comp.mod2)
r.squaredGLMM(sens.comp.mod2)

sens.evg.mod <- lme(Sens.range ~ (Effect-1)*Evergreen.mean - Evergreen.mean, random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.evg.mod); 
anova(sens.evg.mod)
r.squaredGLMM(sens.evg.mod)

sens.decid.mod <- lme(Sens.range ~ (Effect-1)*Deciduous.mean - Deciduous.mean, random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.decid.mod); 
anova(sens.decid.mod)
r.squaredGLMM(sens.decid.mod)

sens.grass.mod <- lme(Sens.range ~ (Effect-1)*Grass.mean - Grass.mean, random=list(Model=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.grass.mod); 
anova(sens.grass.mod)
r.squaredGLMM(sens.grass.mod)
# -----------------

# ----------------
# Explaining differences among Models given site-based differences
# ----------------
sens.moist.site <- lme(Sens.range ~ (Effect-1)*SoilMoist.mean, random=list(Site=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.site); 
anova(sens.moist.site)
r.squaredGLMM(sens.moist.site)


sens.comp.site1 <- lme(Sens.range ~ (Effect-1)*Evergreen.mean*Grass.mean*Deciduous.mean, random=list(Site=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.site); 
anova(sens.comp.site1)
r.squaredGLMM(sens.comp.site1)

sens.comp.site2 <- lme(Sens.range ~ (Effect-1)*(Evergreen.mean + Grass.mean + Deciduous.mean), random=list(Site=~1), data=sim.agg[complete.cases(sim.agg),])
# summary(sens.comp.site2); 
anova(sens.comp.site2)
r.squaredGLMM(sens.comp.site2)

sens.evg.site <- lme(Sens.range ~ (Effect-1)*Evergreen.mean - Evergreen.mean, random=list(Site=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.evg.site); 
# anova(sens.evg.site)
r.squaredGLMM(sens.evg.site)

sens.decid.site <- lme(Sens.range ~ (Effect-1)*Deciduous.mean - Deciduous.mean, random=list(Site=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.decid.site); 
# anova(sens.decid.site)
r.squaredGLMM(sens.decid.site)

sens.grass.site <- lme(Sens.range ~ (Effect-1)*Grass.mean - Grass.mean, random=list(Evergreen=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.grass.site); 
# anova(sens.grass.site)
r.squaredGLMM(sens.grass.site)
# -----------------

# Looking at some overall model characteristics given differences in composition among sites & models
sens.veg.scheme <- lme(Sens.range ~ (Effect-1)*(VegScheme), random=list(Site=~1, Evergreen.mean=~1, Deciduous.mean=~1, Grass.mean=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.veg.scheme); 
anova(sens.veg.scheme)
r.squaredGLMM(sens.veg.scheme)

sens.nlim <- lme(Sens.range ~ (Effect-1)*(Nlim), random=list(Site=~1, Evergreen.mean=~1, Deciduous.mean=~1, Grass.mean=~1), data=sim.agg[complete.cases(sim.agg),])
summary(sens.nlim); 
anova(sens.nlim)
r.squaredGLMM(sens.nlim)

sens.nlim0 <- lm(Sens.range ~ (Effect-1)*(Nlim), data=sim.agg[complete.cases(sim.agg),])
summary(sens.nlim0)
# -------------------------------------------------- #
# Results summary                                    #
# -------------------------------------------------- #
# Effects of Composition/Biogeogrpahy                #
# -----------------------------------                #
# 1) Once you look at explanations within the model, #
#    composition has the greatest effect on precip   #
#    (marginal R2 = 0.06) and least effect on CO2    #
#    (marginal R2 = 0.02)                            #
# 2) Within a site, composition has greatest effect  #
#    on CO2 (marginal R2=0.225)  --> composition     #
#    (Evergreen & Grass separately) explains a       #
#    considerable portion of differences in model    #
#    CO2 sensitivity (moreso than tair & precip)     #
#                                                    #
#
# Effects of Soil Moisture (Site Characterization)   #
# -----------------------------------                #
# 1) Within a model, soil moisture explains much     #
#    more variation in tair & precip sensitivity     #
#    than composition (makes sense)                  #
# 2) Within a site, biogeography explains a lot more #
#    of the variation in temp & co2 sensitivity;     #
#    explanatory of precipf sensitivity is similar   #
# 3) Combining biogeography & soil moisture explains #
#    29% of the variation (R2=0.29) among site and   #
#    27% of the variation (R2=0.27) among models     #
#
#                                                    #
# Effects of Model Properties (Veg Scheme, N Lim)    # 
# -----------------------------------                #
# 1) Veg Scheme -- Static Veg has lower precip sens  #
#    & higher CO2 sensitivity                        #
# 2) Nlim significantly decreased sensitivity to     #
#    precip & CO2 ONCE you take into acount comp.    #
#    & site variation                                #
#    -- Nlim explains a more variation among models  #
#       than veg scheme (mar, R2 = 0.33 vs. 0.24)    #
# -------------------------------------------------- #

# Graphing some of the results
summary(sim.agg)
sim.simple1 <- aggregate(sim.agg[,c("Sens.range", "NPP.mean", "NPP.rel.mean", "tair.mean", "precipf.mean", "CO2.mean", "Evergreen.mean", "Deciduous.mean", "Grass.mean", "SoilMoist.mean", "AGB.mean")], by=list(sim.agg$GAM, sim.agg$Model, sim.agg$Site, sim.agg$Effect), FUN=mean)
names(sim.simple1)[1:4] <- c("GAM", "Model", "Site", "Effect")
sim.simple1$Sens.range.lb <- aggregate(sim.agg$Sens.range, by=list(sim.agg$GAM, sim.agg$Model, sim.agg$Site, sim.agg$Effect), FUN=quantile, 0.025)[,"x"]
sim.simple1$Sens.range.ub <- aggregate(sim.agg$Sens.range, by=list(sim.agg$GAM, sim.agg$Model, sim.agg$Site, sim.agg$Effect), FUN=quantile, 0.975)[,"x"]


pfts <- as.factor(c("Evergreen", "Deciduous", "Grass"))
for(i in 1:nrow(sim.simple1)){
	if(!is.na(max(sim.simple1[i,paste0(pfts, ".mean")]))){
	sim.simple1[i,"PFT.dominant"] <- pfts[which(sim.simple1[i,paste0(pfts, ".mean")]==max(sim.simple1[i,paste0(pfts, ".mean")]))]
	}
}
summary(sim.simple1)

dim(sim.agg)
dim(sim.simple1)




pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracEvergreen_byModel.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_line(stat="smooth", aes(x= Evergreen.mean, y=Sens.range, group=Model, color=Model, fill=Model), method="lm", level=0, size=2, alpha=1) +
	geom_errorbar(aes(x= Evergreen.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x=Evergreen.mean, y=Sens.range, color=Model), size=5, alpha=0.4) +
	scale_x_continuous(name="Fraction Evergreen") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracEvergreen_Overall.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_ribbon(stat="smooth", aes(x= Evergreen.mean, y=Sens.range), method="lm", size=2, alpha=0.3) +
	geom_line(stat="smooth", aes(x= Evergreen.mean, y=Sens.range), method="lm", size=2, alpha=1) +
	geom_errorbar(aes(x= Evergreen.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x=Evergreen.mean, y=Sens.range, color=Model), size=5, alpha=0.8) +
	scale_x_continuous(name="Fraction Evergreen") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()


pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracDeciduous_byModel.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_line(stat="smooth", aes(x= Deciduous.mean, y=Sens.range, group=Model, color=Model, fill=Model), method="lm", level=0, size=2, alpha=1) +
	geom_errorbar(aes(x= Deciduous.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x= Deciduous.mean, y=Sens.range, color=Model), size=5, alpha=0.4) +
	scale_x_continuous(name="Fraction Deciduous") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracDeciduous_Overall.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_ribbon(stat="smooth", aes(x= Deciduous.mean, y=Sens.range), method="lm", size=2, alpha=0.3) +
	geom_line(stat="smooth", aes(x= Deciduous.mean, y=Sens.range), method="lm", size=2, alpha=1) +
	geom_errorbar(aes(x= Deciduous.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x= Deciduous.mean, y=Sens.range, color=Model), size=5, alpha=0.8) +
	scale_x_continuous(name="Fraction Deciduous") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()


pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracGrass_byModel.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_line(stat="smooth", aes(x= Grass.mean, y=Sens.range, group=Model, color=Model, fill=Model), method="lm", level=0, size=2, alpha=1) +
	geom_errorbar(aes(x= Grass.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x=Grass.mean, y=Sens.range, color=Model), size=5, alpha=0.4) +
	scale_x_continuous(name="Fraction Grass") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_FracGrass_Overall.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	geom_ribbon(stat="smooth", aes(x= Grass.mean, y=Sens.range), method="lm", size=2, alpha=0.3) +
	geom_line(stat="smooth", aes(x= Grass.mean, y=Sens.range), method="lm", size=2, alpha=1) +
	geom_errorbar(aes(x= Grass.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x= Grass.mean, y=Sens.range, color=Model), size=5, alpha=0.8) +
	scale_x_continuous(name="Fraction Deciduous") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_SoilMoist.pdf"))
for(e in unique(sim.simple1$Effect)){
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$Effect==e,]) + facet_wrap(~Effect) +
	# geom_ribbon(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e, ], stat="smooth", aes(x=SoilMoist.mean, y=Sens.range, group=Model, fill=Model), method="lm", size=2, alpha=0.5) +
	geom_line(data=sim.agg[sim.agg$GAM=="Site" & sim.agg$Effect==e, ], stat="smooth", aes(x=SoilMoist.mean, y=Sens.range, group=Model, color=Model, fill=Model), method="lm", level=0, size=2, alpha=1) +
	geom_errorbar(aes(x= SoilMoist.mean, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x= SoilMoist.mean, y=Sens.range, color=Model), size=5, alpha=0.4) +
	scale_x_continuous(name="Soil Moisture") +
	scale_y_continuous(name=paste0(e, " Sensitivity Range")) +
	scale_color_manual(values=paste(colors.use)) +
	scale_fill_manual(values=paste(colors.use)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()

model.colors$Model2 <- c("clm.bgc", "clm.cn", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "linkages", "lpj.guess", "lpj.wsl", "sibcasa", "clm", "clm.ed")
# 
summary(sim.simple1[sim.simple1$Site=="PDL" & sim.simple1$GAM=="Site" & sim.simple1$Effect=="Precip",])

pdf(file.path(fig.dir, "Model_NPP_DriverSensitivity_byEffect_byDominantPFT.pdf"))
for(pft in levels(sim.simple1$PFT.dominant)){
models.use2 <- unique(sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$PFT.dominant==pft & !is.na(sim.simple1$PFT.dominant),"Model"])
colors.use2 <- as.vector(model.colors[model.colors$Model2 %in% models.use2, "color"])

sim.simple1$x.graph <- sim.simple1[,paste0(pft, ".mean")]
print(
ggplot(data=sim.simple1[sim.simple1$GAM=="Site" & sim.simple1$PFT.dominant==pft & !is.na(sim.simple1$PFT.dominant),]) + facet_grid(Effect~PFT.dominant) +
	geom_line(stat="smooth", aes(x= x.graph, y=Sens.range, group=Model, color=Model, fill=Model), method="lm", level=0, size=2, alpha=1) +
	geom_errorbar(aes(x= x.graph, ymin=Sens.range.lb, ymax=Sens.range.ub, color=Model), alpha=0.4) +
	geom_point(aes(x= x.graph, y=Sens.range, color=Model), size=5, alpha=0.4) +
	scale_x_continuous(name=paste0("Fraction ", pft)) +
	scale_y_continuous(name="Sensitivity Range") +
	scale_color_manual(values=paste(colors.use2)) +
	scale_fill_manual(values=paste(colors.use2)) +
	guides(fill=F) +
	theme_bw()
)
}
dev.off()
# -------------------------------------------------------

# -------------------------------------------------------
# Going by
# -------------------------------------------------------
summary(ci.terms)
summary(ecosys.agg)
veg.stat <- c("clm.bgc", "clm.cn", "jules.stat", "sibcasa")
P.Collatz <- c("clm.bgc", "clm.cn", "lpj.guess", "lpj.wsl", "jules.stat", "jules.triffid", "sibcasa")
P.Farquar <- c("ed2", "ed2.lu")
met.all <- c("ed2", "ed2.lu", "jules.stat", "jules.triffid", "sibcasa", "clm.cn", "clm.bgc")

N.lim <- c("lpj.guess", "clm.cn", "clm.bgc", "jules.stat", "jules.triffid")

ci.terms2 <- merge(ci.terms, ecosys.agg[,c("Model", "Site", "PFT.dominant")])
# ci.terms2 <- ci.terms2[!is.na(ci.terms2$PFT.dominant),]
ci.terms2$VegScheme <- as.factor(ifelse(ci.terms2$Model %in% veg.stat, "Prescribed", "Emergent"))
ci.terms2$Nlim <- as.factor(ifelse(ci.terms2$Model %in% N.lim, "Yes", "No"))
ci.terms2$nMet <- as.factor(ifelse(ci.terms2$Model %in% met.all, 8, ifelse(ci.terms2$Model=="lpj.wsl", 5, 4)))
summary(ci.terms2)

models.use2 <- unique(ci.terms2$Model)
colors.use2 <- as.vector(model.colors[model.colors$Model2 %in% models.use2, "color"])

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms2
ci.terms.graph[ci.terms.graph$mean.rel<(-0.75),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.75),"lwr.rel"] <- -0.75 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.75),"upr.rel"] <- -0.75 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.0),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(1.0),"lwr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$upr.rel>(1.0),"upr.rel"] <- 1.0 
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)

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

}
ci.terms.graph$x.min    <- ifelse(ci.terms.graph$x < ci.terms.graph$line.min, ci.terms.graph$x, ci.terms.graph$line.min)
ci.terms.graph$x.max    <- ifelse(ci.terms.graph$x > ci.terms.graph$line.max, ci.terms.graph$x, ci.terms.graph$line.max)
ci.terms.graph$mask.min <- min(ci.terms.graph$lwr.rel, na.rm=T)
ci.terms.graph$mask.max <- max(ci.terms.graph$upr.rel, na.rm=T)

summary(ci.terms.graph)

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_bySite_byDominantPFT.pdf"))
for(s in unique(ci.terms.graph$Site)){
print(
ggplot(data=ci.terms.graph[ci.terms.graph$Site==s & ci.terms.graph$GAM=="Site" & !is.na(ci.terms.graph$PFT.dominant),]) + facet_grid(PFT.dominant~Effect, scales="free_x") +
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

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_bySite_byDominantPFT_byEffect_PDL.pdf"), width=15.5, height=8.5)
for(e in unique(ci.terms.graph$Effect)){

# Specify the x-axis for each type
if(e=="tair") x.axis <- scale_x_continuous(name=expression(bold(paste("Temp (May-Sep, "^"o","C)"))), breaks=seq(9,19, by=3), expand=c(0,0))
if(e=="precipf") x.axis <- scale_x_continuous(name="Precip (May-Sep, mm)", breaks=seq(250, 800, by=250), expand=c(0,0))
if(e=="CO2") x.axis <- scale_x_continuous(name=expression(bold(paste("CO"["2"]," (ppm)"))), expand=c(0,0))

print(
ggplot(data=ci.terms.graph[ci.terms.graph$Site=="PDL" & ci.terms.graph$GAM=="Site" & ci.terms.graph$Effect==e & !is.na(ci.terms.graph$PFT.dominant),]) + facet_grid(.~PFT.dominant) +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.min), linetype="dashed") +
	geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	geom_vline(aes(xintercept=line.max), linetype="dashed") +
	ggtitle("Demming Lake (MN)") +
	x.axis +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
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
	      axis.ticks.margin=unit(1.0, "lines")) +
	theme(plot.margin=unit(c(1.5,1,1.15,1), "lines"))
)
}
dev.off()



pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Overall_byVegScheme.pdf"))

levels(ci.terms.graph$Effect) <- c("Temperature", "Precipitation", "CO2")

plot.tair <- ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline" & ci.terms.graph$Effect=="Temperature", ]) + facet_grid(VegScheme~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name=expression(bold(paste("Temp (May-Sep, "^"o","C)"))), breaks=seq(9,19, by=3)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text.x=element_text(size=rel(2.5), face="bold"),
	      strip.text.y=element_blank()) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.margin=unit(0, "lines"),
	      panel.background=element_blank()) +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_text(color="black", size=rel(3)), 
		  axis.title.x=element_text(size=rel(2.25), face="bold", vjust=-0.4),  
		  axis.title.y=element_text(size=rel(2.5), face="bold", vjust=1),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	theme(plot.margin=unit(c(4,-1,1.6,1), "lines"))

plot.precip <- ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline" & ci.terms.graph$Effect=="Precipitation", ]) + facet_grid(VegScheme~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name="Precip (May-Sep, mm)") +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	guides(fill=F, color=F) +
	theme(strip.text.x=element_text(size=rel(2.5), face="bold"),
	      strip.text.y=element_blank()) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.margin=unit(0, "lines"),
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.25), face="bold", vjust=-1),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(4,-1,2.33,-1), "lines"))

plot.co2 <- ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline" & ci.terms.graph$Effect=="CO2", ]) + facet_grid(VegScheme~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	scale_x_continuous(name=expression(bold(paste("CO"["2"]," (pmm)"))), breaks=seq(280, 400, by=30)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0), limits=range(ci.terms.graph[,c("lwr.rel", "upr.rel")])*100) +
	scale_fill_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	scale_color_manual(values=colors.use, labels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "SibCASA")) +
	ggtitle("Cross-Site Sensitivity") +
	theme(plot.title=element_text(size=rel(3), face="bold", hjust=3.5, vjust=2)) + 
	guides(fill=guide_legend(title.position="top", ncol=3), 
	       color=guide_legend(title.position="top", ncol=3)) +
	theme(legend.title=element_text(size=rel(1.5), face="bold"),
	      legend.text=element_text(size=rel(1.25)),
	      legend.position=c(0.2, 0.1),
	      legend.key=element_blank(),
	      legend.key.size=unit(1.5, "lines")) +
	theme(strip.text.x=element_text(size=rel(2.5), face="bold"),
          strip.text.y=element_text(size=rel(2.5), face="bold")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_rect(size=0.5, color="black", fill=NA), 
	      panel.margin=unit(0, "lines"),
	      panel.background=element_blank())  +
	theme(axis.text.x=element_text(color="black", size=rel(3)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold", vjust=-0.8),  
		  axis.title.y=element_blank(),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	theme(plot.margin=unit(c(1.7,1,1.7,-1.1), "lines"))


pdf(file.path(fig.dir, "NPP_Sensitivity_Models_VegScheme_Presentation.pdf"), height=9, width=14)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3, widths=c(1.1,0.85,0.9))))
print(plot.tair  , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.precip, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
print(plot.co2   , vp = viewport(layout.pos.row = 1, layout.pos.col=3))
dev.off()


# for(s in unique(ci.terms.graph$Site)){
# print(
ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline",]) + facet_grid(VegScheme~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	# geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.min), linetype="dashed") +
	# geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.max), linetype="dashed") +
	# ggtitle(s) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
# )
# }
dev.off()

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Overall_byNLim.pdf"))
# for(s in unique(ci.terms.graph$Site)){
# print(
ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline",]) + facet_grid(Nlim~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	# geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.min), linetype="dashed") +
	# geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.max), linetype="dashed") +
	# ggtitle(s) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
# )
# }
dev.off()

pdf(file.path(fig.dir, "NPP_Sensitivity_Models_Overall_nMet.pdf"))
# for(s in unique(ci.terms.graph$Site)){
# print(
ggplot(data=ci.terms.graph[ci.terms.graph$GAM=="Baseline",]) + facet_grid(nMet~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model)) +
	# geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.min), linetype="dashed") +
	# geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.5) +
	# geom_vline(aes(xintercept=line.max), linetype="dashed") +
	# ggtitle(s) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	theme_bw()
# )
# }
dev.off()


# # -------------------------------------------------------
# # the Bayesian way
# # -------------------------------------------------------
# # Interactions only
# comp.effects <- function(){
	# # Priors
	# beta0 ~ dunif(0,100)
	# beta1 ~ dunif(0,100)
	# beta2 ~ dunif(0,100)
	# beta3 ~ dunif(0,100)
	# sigma ~ dunif(0,100)
	# # Model
	# for(i in 1:n){
		# # mu[i] <- beta*TEMP[i]
		# mu[i] <- beta0 + beta1*EVERGREEN[i]*GRASS[i]+ 
				 # beta2*EVERGREEN[i] + beta3*GRASS[i]
		# y[i] ~  dnorm(mu[i], sigma)
	# }
# }

# SENS      <- sim.agg[sim.agg$Effect=="CO2" & complete.cases(sim.agg), "Sens.range"    ]
# n         <- length(SENS)
# EVERGREEN <- sim.agg[sim.agg$Effect=="CO2" & complete.cases(sim.agg), "Evergreen.mean"]
# GRASS     <- sim.agg[sim.agg$Effect=="CO2" & complete.cases(sim.agg), "Grass.mean"    ]
# params    <- c("beta0", "beta1", "beta2", "beta3", "sigma")

# co2.comp1 <- jags(data=list(y=SENS, n=n, EVERGREEN=EVERGREEN, GRASS=GRASS), parameters.to.save=params, n.chains=3, n.iter=2000, n.burnin=500, model.file=comp.effects, DIC=F)
# summary(co2.comp1$model)

# co2.comp1 <- as.mcmc(co2.comp1)
# summary(co2.comp1)
# par(mfrow=c(round((length(params)+.5)/2, 0), 2))
# traceplot(co2.comp1)

# par(mfrow=c(length(params), 2))
# plot(co2.comp1)

# # Generating the predictions
# pulls <- 500
# y.comp1 <- array(dim=c(n, pulls))
# for(i in 1:pulls){
	# c <- sample(1:length(co2.comp1), 1, replace=T)    # randomly pick a chain
	# r <- sample(1:nrow(co2.comp1[[c]]), 1, replace=T) # randomly pick an iteration 

	# y.comp1[,i] <- co2.comp1[[c]][r,"beta0"] + 
				   # co2.comp1[[c]][r,"beta1"]*EVERGREEN*GRASS +
				   # co2.comp1[[c]][r,"beta2"]*EVERGREEN + 
				   # co2.comp1[[c]][r,"beta3"]*GRASS 
# }

# y.pred1 <- apply(y.comp1, 1, mean, na.rm=T)
# y.lb1 <- apply(y.comp1, 1, quantile, 0.025, na.rm=T)
# y.ub1 <- apply(y.comp1, 1, quantile, 0.975, na.rm=T)

# lm5 <- lm(y.pred1 ~ SENS)
# summary(lm5)
# # -------------------------------------------------------


# ----------------------------------------

