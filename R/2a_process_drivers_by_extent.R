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
# (Fit GAMM per site per m.name)
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
#    -- Hypothesis: Driver responses across sites within a m.name converge at coarser temporal grains 
#       and larger extents because the models have time to adjust and seek equilibrium.
#
#    -- Analysis: Use the posterior CIs for each smoothing term to see if the driver curves for sites
#                 within a m.name are statstically different at different sites at different scales (or 
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
#    -- Analysis: For each given scale, compare the m.name response curves for each driver.  Determine 
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
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path(path.data, "EcosysData.Rdata"))

# Scripts to run the gamms to predict a response variable as a function of Temp, Precip, & CO2
# 	predict.gamm.model.site.R = function to run a single 1 site - m.name combo at a time (fit curves independently)
# 	predict.gamm.mode.R		= function to get overal m.name responses with random site effects 
# 	Note: these two functions were split because they now incorporate AR1 autocorrelation that can make the 
#		  overal m.name fitting with random site effects very slow
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
# All Sites: (for 1 site, see m.name selection script)
# ------------------------------------------------
data.base="~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms"
fig.base="~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/gamms"

# ------------------------
# MegaLoop -- Looping through all models by Variable, by Extent
# ------------------------
# Setting up a loop for 1 m.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- unique(ecosys$Model)
model.order   <- unique(ecosys$Model.Order)
var <- c("NPP", "AGB.diff")
# scale    <- ""
scales <- c("", ".10", ".50", ".100", ".250")
# scales <- c(".100")
t.scales <- ifelse(scales=="", "t.001", paste0("t", scales))
extents <- data.frame(Start=c(850, 1900, 1990), End=c(2010, 2010, 2010)) 

for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]
	out.dir   <- file.path(data.base, m.order, "ByExtent")
	fig.dir  <- file.path(fig.base, m.order, "ByExtent")

	print(" ")
	print(" ")
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Model: ",m.name, " ------"))

for(v in var){
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Variable: ",v, " ------"))
for(e in 1:nrow(extents)){
	print(       "      ----------------------      ")
	print(paste0("------ Processing Extent: ", extents[e,1], " - ", extents[e,2], " ------"))

	extent <- as.numeric(extents[e,])
	# ecosys2 <- ecosys[ecosys$Model==m.name & ecosys$Year>=extent[1] & ecosys$Year<=extent[2], ]
	ecosys2 <- ecosys[ecosys$Model==m.name & ecosys$Year>=extent[1] & ecosys$Year<=extent[2], ]
	ecosys2$Extent <- as.factor(paste(extent[1], extent[2], sep="-"))
	
for(t in 1:length(scales)){
	print(       "-------------------------------------")
	print(paste0("------ Processing Scale: ", t.scales[t], " ------"))

	mod.temp <- model.gam(data=ecosys2, model.name=m.name, extent=extent, scale=scales[t], response=v, k=4, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
	if(t==1 & e==1) {
		mod.out <- list()
		mod.out$data        <- mod.temp$data
		mod.out$weights     <- mod.temp$weights
		mod.out$ci.response <- mod.temp$ci.response
		mod.out$ci.terms    <- mod.temp$ci.terms
		mod.out[[paste("gamm", paste0(extent[1], "-", extent[2]), t.scales[t], sep=".")]] <- mod.temp$gamm
	} else {
		mod.out$data        <- rbind(mod.out$data,        mod.temp$data)
		mod.out$weights     <- rbind(mod.out$weights,     mod.temp$weights)
		mod.out$ci.response <- rbind(mod.out$ci.response, mod.temp$ci.response)
		mod.out$ci.terms    <- rbind(mod.out$ci.terms,    mod.temp$ci.terms)
		mod.out[[paste("gamm", paste0(extent[1], "-", extent[2]), t.scales[t], sep=".")]] <- mod.temp$gamm
	}
	
	} # end scales
} # end extent
	save(mod.out, file=file.path(out.dir, paste("gamm", m.name, v, "Rdata", sep=".")))
	# assign(paste("gamm", m.name, v, sep="."), mod.out)

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Extent_", m.order, "_", v, "_0850-2010", ".pdf")))
	col.model <- model.colors[model.colors$Model.Order %in% unique(gamm.lpj.guess$data$Model.Order),"color"]
	# ggplot(data= gamm.lpj.guess$ci.response) + facet_wrap(~Scale, scales="free") + theme_bw() +
 ggplot(data=gamm.lpj.guess$ci.response[gamm.lpj.guess$ci.response$Scale=="t.001",]) + facet_grid(Site~Extent, scales="free") + theme_bw() +
 		geom_line(data= gamm.lpj.guess$data[gamm.lpj.guess$data$Scale=="t.001",], aes(x=Year, y=response), alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		# scale_x_continuous(limits=c(850,2010)) +
		# scale_y_continuous(limits=quantile(gamm.lpj.guess$data$response, c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.name), x="Year", y=var)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Extent_", m.order, "_", v, "_1990-2010", ".pdf")))
	col.model <- model.colors[model.colors$Model.Order %in% unique(df$data$Model.Order),"color"]
	ggplot(data=gamm.lpj.guess$ci.response[gamm.lpj.guess$ci.response$Scale=="t.001",]) + facet_wrap(~Site, scales="free") + theme_bw() +
 		geom_line(data=gamm.lpj.guess$data[gamm.lpj.guess$data$Scale=="t.001",], aes(x=Year, y=response), size=1.5, alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
		geom_line(aes(x=Year, y=mean, color=Extent), size=1) +
		scale_x_continuous(limits=c(1990,2010)) +
		# scale_y_continuous(limits=quantile(gamm.lpj.guess$data[gamm.lpj.guess$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.name), x="Year", y=var)
dev.off()


pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_Extent_", m.order, "_", v, ".pdf")))
	# ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Effect=="CO2",]) + facet_grid(Effect ~ Scale, scales="free") + theme_bw() +
 	ggplot(data=gamm.lpj.guess$ci.terms[gamm.lpj.guess$ci.terms$Scale=="t.001",]) + facet_wrap(~ Effect, scales="free") + theme_bw() +		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Extent), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title=paste0("Driver Effects: ",m.order), y="Effect Size") # +
		# theme(legend.position=c(0.75,0.3), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)), legend.key.size=unit(1.5, "line"))
dev.off()


} # end var
} # end model

# save(gamm.lpj.guess.AGB.diff, file=file.path(out.dir, "gamm.lpj.guess.dAGB.Rdata"))
# save(gamm.lpj.guess.NPP, file=file.path(out.dir, "gamm.lpj.guess.NPP.Rdata"))

