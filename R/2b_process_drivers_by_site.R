# ----------------------------------------
# Temporal Scaling Analyses
# Non-constant driver effects through time
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Driving Questions: How does temporal scale affect our inference of Temp, Precip, & CO2 effects on ecosystems?
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# (Fit GAMM per site per model.name)
# 1) Temporal Grain (Resolution)
#    -- Fit GAMM over constant window with different degrees of smoothing (1 yr - 250 yr)
# 2) Temporal Extent 
#	 -- Fit GAMM to different windows (30 yrs 1990-2010, 100 yrs 1910-2010, full window)
# ** Response variables of interest: NPP, possibly dAGB (AGB 1st difference)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# 1) Space-Time Comparison
#    -- Hypothesis: Driver responses across sites within a model.name converge at coarser temporal grains 
#       and larger extents because the models have time to adjust and seek equilibrium.
#
#    -- Analysis: Use the posterior CIs for each smoothing term to see if the driver curves for sites
#                 within a model.name are statstically different at different sites at different scales (or 
#                 alternatively if there are statistical differences in [parts of] the curves between
#                 scales)
#
#
# 2) Multi-Model Driver Comparison
#    -- Hypothesis: Because models were built & parameterized to perform well in the modern era, there 
#       will be greater agreement of the direction & primary driver of change in the more recent time 
#       periods than over the full extent of the PalEON runs.
#    -- Hypothesis about temporal grain?
#
#    -- Analysis: For each given scale, compare the model.name response curves for each driver.  Determine 
#                 which drivers are most similar/variable among models and at what scales?  Are there
#                 particular ranges of each driver response where models responses are most similar/different?
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ncdf4)
library(lme4)
library(R2jags)
library(ggplot2); library(grid)
library(car)
library(zoo)
# library(mvtnorm)
# library(MCMCpack)
# ----------------------------------------

# ----------------------------------------
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
path.data <- "~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
fig.dir <- "~/Desktop/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
# ----------------------------------------


# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed model.name outputs
#	generated with 1_generate_ecosys.R
load(file.path(path.data, "EcosysData.Rdata"))

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2
# 	predict.gamm.model.site.R = function to run a single 1 site - model.name combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal model.name responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal model.name fitting with random site effects very slow
source('~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/predict.gamm.model.site.R', chdir = TRUE)
source('~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/predict.gamm.model.R', chdir = TRUE)
source('~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/GAMM_Plots.R', chdir = TRUE)
# ----------------------------------------


# # -----------------------
# # Some exploratory Graphing
# # -----------------------
# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # ggtitle("Annual & Centennial AGB") +
  # theme_bw()

# # quantile(ecosys$AGB.diff, c(0.025, 0.33, 0.67, 0.975), na.rm=T)
# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # # geom_line(aes(x=Year, y=AGB.diff, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=AGB.diff.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # scale_y_continuous(limits=quantile(ecosys$AGB.diff.100, c(0.025, 0.975), na.rm=T)) +
  # ggtitle("Annual & Centennial dAGB") +
  # theme_bw()

# ggplot(data=ecosys[,]) + facet_wrap(~Site) +
  # geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1, alpha=0.6) +
  # geom_line(aes(x=Year, y=NPP.100, color=Model.Order), size=1.5) +
  # scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
  # # scale_y_continuous(limits=quantile(ecosys$NPP.100, c(0.025, 0.975), na.rm=T)) +
  # ggtitle("Annual & Centennial NPP") +
  # theme_bw()

# # -----------------------

# ----------------------------------------


# ----------------------------------------
# Model approach: AGB ~ 3 non-interactive temporal smoothers: AGB, Temp, Precip
# ----------------------------------------
library(mgcv)

# ------------------------------------------------
# All Sites: (for 1 site, see model.name selection script)
# ------------------------------------------------
outdir="~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"

# ------------------------
# LPJ-GUESS -- Loop to cover all sites for 1 model.name, 1 temporal scale
# ------------------------
# Setting up a loop for 1 model.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- "lpj.guess"
var <- "NPP"
# scale    <- ""
scales <- c("", ".10", ".50", ".100", ".250")
t.scales <- ifelse(scales=="", "t.001", paste0("t", scales))
out.dir   <- "~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms/LPJ-GUESS/"
fig.dir  <- "~/Desktop/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures/gamms/LPJ-GUESS/"

for(t in 1:length(scales)){
	print(       "-------------------------------------")
	print(paste0("------ Processing Scale: ", t.scales[t], " ------"))
for(s in 1:length(sites)){
	print(       "      -------------------------      ")
	print(paste0("------ Processing Site: ", sites[s], " ------"))
	# object.name <- paste("gam", model.name, s, sep=".")
	mod.temp <- model.site.gam(data=ecosys, model.name=model.name, site=sites[s], scale=scales[t], response=var, k=4, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
	if(t==1 & s==1) {
		mod.out <- list()
		mod.out$data        <- mod.temp$data
		mod.out$weights     <- mod.temp$weights
		mod.out$ci.response <- mod.temp$ci.response
		mod.out$ci.terms    <- mod.temp$ci.terms
		mod.out[[paste("gamm", t.scales[t], sites[s], sep=".")]] <- mod.temp$gamm
		# summary(mod.out)
	} else {
		mod.out$data        <- rbind(mod.out$data,        mod.temp$data)
		mod.out$weights     <- rbind(mod.out$weights,     mod.temp$weights)
		mod.out$ci.response <- rbind(mod.out$ci.response, mod.temp$ci.response)
		mod.out$ci.terms    <- rbind(mod.out$ci.terms,    mod.temp$ci.terms)
		mod.out[[paste("gamm", t.scales[t], sites[s], sep=".")]] <- mod.temp$gamm
	}
	
	if(s==length(sites)) assign(paste("gamm", model.name, var, sep="."), mod.out)
}
}

# save(gamm.lpj.guess.AGB.diff, file=file.path(out.dir, "gamm.lpj.guess.dAGB.Rdata"))
save(gamm.lpj.guess.NPP, file=file.path(out.dir, "gamm.lpj.guess.NPP.Rdata"))

gamm.lpj.guess <- gamm.lpj.guess.NPP
summary(gamm.lpj.guess)
summary(gamm.lpj.guess$data)

pdf(file.path(fig.dir, "GAMM_ResponsePrediction_LPJ-GUESS_NPP_0850-2010.pdf"))
	col.model <- model.colors[model.colors$Model.Order %in% unique(gamm.lpj.guess$data$Model.Order),"color"]
	ggplot(data= gamm.lpj.guess$ci.response) + facet_grid(Scale~Site, scales="free") + theme_bw() +
		geom_line(data= gamm.lpj.guess$data, aes(x=Year, y=response), alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		scale_x_continuous(limits=c(850,2010)) +
		# scale_y_continuous(limits=quantile(gamm.lpj.guess$data$response, c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", model.name), x="Year", y=var)
dev.off()

pdf(file.path(fig.dir, "GAMM_ResponsePrediction_LPJ-GUESS_NPP_1900-2010.pdf"))
	col.model <- model.colors[model.colors$Model.Order %in% unique(df$data$Model.Order),"color"]
	ggplot(data= gamm.lpj.guess$ci.response) + facet_grid(Scale~Site, scales="free") + theme_bw() +
		geom_line(data= gamm.lpj.guess$data, aes(x=Year, y=response), size=1.5, alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=1, color= col.model) +
		scale_x_continuous(limits=c(1901,2010)) +
		# scale_y_continuous(limits=quantile(gamm.lpj.guess$data[gamm.lpj.guess$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", model.name), x="Year", y=var)
dev.off()


summary(gamm.lpj.guess$ci.terms)

pdf(file.path(fig.dir, "GAMM_DriverEffects_LPJ-GUESS_NPP.pdf"))
	# ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +
 	ggplot(data=gamm.lpj.guess$ci.terms) + facet_grid(Scale ~ Effect, scales="free") + theme_bw() +		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title="Driver Effects: LPJ-GUESS", y="Effect Size") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverEffects_LPJ-GUESS_NPP_CO2.pdf"))
	# ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +
 	ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +		
 		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title="Driver Effects: LPJ-GUESS", y="Effect Size", x="CO2") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverEffects_LPJ-GUESS_NPP_Temp.pdf"))
	# ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +
 	ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="Temp",]) + 
 		# facet_grid(Effect ~ Scale, scales="free") + theme_bw() +		
		facet_grid(Scale ~ Effect, scales="free") + theme_bw() +		
 		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title="Driver Effects: LPJ-GUESS", y="Effect Size", x="Temp (K)") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverEffects_LPJ-GUESS_NPP_Precip.pdf"))
	# ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +
 	ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="Precip",]) + 
 		# facet_grid(Effect ~ Scale, scales="free") + theme_bw() +		
		facet_grid(Scale ~ Effect, scales="free") + theme_bw() +		
 		geom_ribbon(aes(x=x, ymin=lwr*sec2yr, ymax=upr*sec2yr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean*sec2yr, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title="Driver Effects: LPJ-GUESS", y="Effect Size", x="Precip (mm/yr)") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
dev.off()

summary(gam.lpj.guess.pha$weights)
# gam.lpj.guess.pha$weights$Year <- gam.lpj.guess.pha$data$Year

site="PHA"
df=gamm.lpj.guess
pdf(file.path(fig.dir, "GAMM_DriverTime_LPJ-GUESS_AGB_0850-2010.pdf"), width=11, height=8.5)
ggplot(data= gamm.lpj.guess$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= gamm.lpj.guess$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y="NPP kgC/m2/s", title="Drivers through Time: LPJ-GUESS, 0850-2010") +
	scale_x_continuous(limits=c(850,2010)) +
	scale_y_continuous(limits=quantile(gamm.lpj.guess$data$response, c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverTime_LPJ-GUESS_AGB_1900-2010.pdf"), width=11, height=8.5)
ggplot(data= gamm.lpj.guess$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= gamm.lpj.guess$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y="NPP kgC/m2/s", title="Drivers through Time: LPJ-GUESS, 1900-2010") +
	scale_x_continuous(limits=c(1900,2010)) +
	scale_y_continuous(limits=quantile(gamm.lpj.guess$data[gamm.lpj.guess$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverTime_LPJ-GUESS_AGB_1800-1900.pdf"), width=11, height=8.5)
ggplot(data= gamm.lpj.guess$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= gamm.lpj.guess$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=gamm.lpj.guess, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y="NPP kgC/m2/s", title="Drivers through Time: LPJ-GUESS, 1800-1900") +
	scale_x_continuous(limits=c(1800,1900)) +
	scale_y_continuous(limits=quantile(gamm.lpj.guess$data[gamm.lpj.guess$data$Year>=1800 & gamm.lpj.guess$data$Year<=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# ------------------------




# ------------------------
# ED -- Loop to cover all sites for 1 model.name, 1 temporal scale
# ------------------------
# Setting up a loop for 1 model.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- "ed2"
var <- "AGB.diff"
scale    <- ""
	t.scale <- ifelse(scale=="", "t001", paste0("t", scale))
out.dir   <- "~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms/ED/dAGB.001.sites"
fig.dir  <- "~/Desktop/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures/gamms/ED/dAGB.001.sites"

for(s in 1:length(sites)){
	print(paste0("------ Processing Site: ", sites[s], " ------"))
	# object.name <- paste("gam", model.name, s, sep=".")
	mod.temp <- model.site.gam(data=ecosys, model.name=model.name, site=sites[s], scale="", response=var, k=4, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
	if(s==1) {
		mod.out <- list()
		mod.out$data        <- mod.temp$data
		mod.out$weights     <- mod.temp$weights
		mod.out$ci.response <- mod.temp$ci.response
		mod.out$ci.terms    <- mod.temp$ci.terms
		mod.out[[paste("gamm", sites[s], sep=".")]] <- mod.temp$gamm
		# summary(mod.out)
	} else {
		mod.out$data        <- rbind(mod.out$data,        mod.temp$data)
		mod.out$weights     <- rbind(mod.out$weights,     mod.temp$weights)
		mod.out$ci.response <- rbind(mod.out$ci.response, mod.temp$ci.response)
		mod.out$ci.terms    <- rbind(mod.out$ci.terms,    mod.temp$ci.terms)
		mod.out[[paste("gamm", sites[s], sep=".")]] <- mod.temp$gamm
	}
	
	if(s==length(sites)) assign(paste("gamm", model.name, var, t.scale, sep="."), mod.out)
}
save(gamm.ed2.AGB.diff.t.001, file=file.path(out.dir, "gamm.ed2.dAGB.001.Rdata"))
# gamm.ed2 <- gamm.ed2.AGB.diff.t.001
# summary(gamm.ed2)

pdf(file.path(fig.dir, "GAMM_ResponsePrediction_ED_dAGB_0850-2010.pdf"))
	col.model <- model.colors[model.colors$Model.Order %in% unique(df$data$Model.Order),"color"]
	ggplot(data= gamm.ed2$ci.response) + facet_wrap(~Site) + theme_bw() +
		geom_line(data= gamm.ed2$data, aes(x=Year, y=response), alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		scale_x_continuous(limits=c(850,2010)) +
		scale_y_continuous(limits=quantile(gamm.ed2$data$response, c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ", ", ifelse(scale=="", "  1", scale), "yr: ", model.name), x="Year", y=var)
dev.off()

pdf(file.path(fig.dir, "GAMM_ResponsePrediction_ED_dAGB_1900-2010.pdf"))
	col.model <- model.colors[model.colors$Model.Order %in% unique(df$data$Model.Order),"color"]
	ggplot(data= gamm.ed2$ci.response) + facet_wrap(~Site) + theme_bw() +
		geom_line(data= gamm.ed2$data, aes(x=Year, y=response), size=1.5, alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=1, color= col.model) +
		scale_x_continuous(limits=c(1900,2010)) +
		scale_y_continuous(limits=quantile(gamm.ed2$data[gamm.ed2$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ", ", ifelse(scale=="", "  1", scale), "yr: ", model.name), x="Year", y=var)
dev.off()


summary(gamm.ed2$ci.terms)

pdf(file.path(fig.dir, "GAMM_DriverEffects_ED_dAGB.pdf"))
	ggplot(data=gamm.ed2$ci.terms) + facet_wrap(~Effect, scales="free_x",ncol=2) + theme_bw() +
		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title="Driver Effects: ED", y="Effect Size") +
		theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2)), legend.key.size=unit(2, "line"))
dev.off()

summary(gam.ed2.pha$weights)
# gam.ed2.pha$weights$Year <- gam.ed2.pha$data$Year

site="PHA"
df=gamm.ed2
pdf(file.path(fig.dir, "GAMM_DriverTime_ED_dAGB_0850-2010.pdf"), width=11, height=8.5)
ggplot(data= gamm.ed2$weights) + facet_wrap(~Site) +
	geom_line(data= gamm.ed2$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	# geom_line(aes(x=Year, y=fit, color=rgb(abs(weight.temp), abs(weight.co2), abs(weight.precip))), size=4) + 
	rgb.line(df=gamm.ed2, site="PHA", size=2) +
	rgb.line(df=gamm.ed2, site="PHO", size=2) +
	rgb.line(df=gamm.ed2, site="PUN", size=2) +
	rgb.line(df=gamm.ed2, site="PBL", size=2) +
	rgb.line(df=gamm.ed2, site="PDL", size=2) +
	rgb.line(df=gamm.ed2, site="PMB", size=2) +
	labs(x="Year", y="dAGB kg/m2", title="Drivers through Time: ED, 0850-2010") +
	scale_x_continuous(limits=c(850,2010)) +
	scale_y_continuous(limits=quantile(gamm.ed2$data$response, c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverTime_ED_dAGB_1900-2010.pdf"), width=11, height=8.5)
ggplot(data= gamm.ed2$weights) + facet_wrap(~Site) +
	geom_line(data= gamm.ed2$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	# geom_line(aes(x=Year, y=fit, color=rgb(abs(weight.temp), abs(weight.co2), abs(weight.precip))), size=4) + 
	rgb.line(df=gamm.ed2, site="PHA", size=2) +
	rgb.line(df=gamm.ed2, site="PHO", size=2) +
	rgb.line(df=gamm.ed2, site="PUN", size=2) +
	rgb.line(df=gamm.ed2, site="PBL", size=2) +
	rgb.line(df=gamm.ed2, site="PDL", size=2) +
	rgb.line(df=gamm.ed2, site="PMB", size=2) +
	labs(x="Year", y="dAGB kg/m2", title="Drivers through Time: ED, 0850-2010") +
	scale_x_continuous(limits=c(1900,2010)) +
	scale_y_continuous(limits=quantile(gamm.ed2$data[gamm.ed2$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverTime_ED_dAGB_1800-1900.pdf"), width=11, height=8.5)
ggplot(data= gamm.ed2$weights) + facet_wrap(~Site) +
	geom_line(data= gamm.ed2$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	# geom_line(aes(x=Year, y=fit, color=rgb(abs(weight.temp), abs(weight.co2), abs(weight.precip))), size=4) + 
	rgb.line(df=gamm.ed2, site="PHA", size=2) +
	rgb.line(df=gamm.ed2, site="PHO", size=2) +
	rgb.line(df=gamm.ed2, site="PUN", size=2) +
	rgb.line(df=gamm.ed2, site="PBL", size=2) +
	rgb.line(df=gamm.ed2, site="PDL", size=2) +
	rgb.line(df=gamm.ed2, site="PMB", size=2) +
	labs(x="Year", y="dAGB kg/m2", title="Drivers through Time: ED, 0850-2010") +
	scale_x_continuous(limits=c(1800,1900)) +
	scale_y_continuous(limits=quantile(gamm.ed2$data[gamm.ed2$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# ------------------------

