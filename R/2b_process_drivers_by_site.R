# ----------------------------------------
# Temporal Scaling Analyses: Driver Effects Through Time
# Space-Time interactions by separating out model responses by site 
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
# 2) Temporal Extent -- not currently looked at at the site level -- model only
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
setwd("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
path.data <- "~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
fig.dir <- "~/Desktop/Research/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
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
source('~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_predict.gamm.model.site.R', chdir = TRUE)
source('~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_predict.gamm.model.R', chdir = TRUE)
source('~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/0_GAMM_Plots.R', chdir = TRUE)



# Read in model color scheme
model.colors <- read.csv("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors

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
data.base="~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms_byModel"
fig.base="~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Figures/gamms_byModel"

# ------------------------
# MegaLoop -- Looping through all models by Variable, by Site
# ------------------------
# Setting up a loop for 1 model.name, 1 temporal scale
sites    <- unique(ecosys$Site)
model.name    <- unique(ecosys$Model)
model.order   <- unique(ecosys$Model.Order)
var <- c("NPP", "AGB.diff")
# scale    <- ""
scales <- c("", ".10", ".50", ".100", ".250")
t.scales <- ifelse(scales=="", "t.001", paste0("t", scales))

extent <- c(850, 2010)

# out.dir   <- "~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data/gamms/LPJ-GUESS/"
# fig.dir  <- "~/Desktop/Research/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures/gamms/LPJ-GUESS/"

for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]
	out.dir   <- file.path(data.base, m.order, "BySite")
	fig.dir  <- file.path(fig.base, m.order, "BySite")

	print(" ")
	print(" ")
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Model: ",m.name, " ------"))
	print(paste0("------ Fig Dir: ",fig.dir, " ------"))
	print(paste0("------ Dat Dir: ",out.dir, " ------"))

for(v in var){
	print(" ")
	print(" ")
	print(       "      ----------------------      ")
	print(paste0("------ Processing Variable: ",v, " ------"))

	ecosys2 <- ecosys[ecosys$Model==m.name & ecosys$Year>=extent[1] & ecosys$Year<=extent[2], ]
	ecosys2$Extent <- as.factor(paste(extent[1], extent[2], sep="-"))

for(t in 1:length(scales)){
	print(       "-------------------------------------")
	print(paste0("------ Processing Scale: ", t.scales[t], " ------"))

for(s in 1:length(sites)){
	print(       "      -------------------------      ")
	print(paste0("------ Processing Site: ", sites[s], " ------"))
	# object.name <- paste("gam", model.name, s, sep=".")
	mod.temp <- model.site.gam(data=ecosys2, model.name=m.name, extent=extent, site=sites[s], scale=scales[t], response=v, k=4, write.out=F, outdir=out.dir, fweights=T, ci.model=T, ci.terms=T)
	
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
	
} # end sites
} # end scales
	save(mod.out, file=file.path(out.dir, paste("gamm", m.name, v, "Rdata", sep=".")))
	# assign(paste("gamm", model.name, var, sep="."), mod.out)


# summary(mod.out)
# summary(mod.out$data)
col.model <- model.colors[model.colors$Model.Order %in% unique(mod.out$data$Model.Order),"color"]

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Site_", m.order, "_", v, "_0850-2010", ".pdf")))
print(
	ggplot(data= mod.out$ci.response) + facet_grid(Scale~Site, scales="free") + theme_bw() +
		geom_line(data= mod.out$data, aes(x=Year, y=response), alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		scale_x_continuous(limits=c(850,2010), breaks=c(1300,1800)) +
		# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.order, " 0850-2010"), x="Year", y=var)
		)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_ResponsePrediction_Site_", m.order, "_", v, "_1900-2010", ".pdf")))
print(
	ggplot(data= mod.out$ci.response) + facet_grid(Scale~Site, scales="free") + theme_bw() +
		geom_line(data= mod.out$data, aes(x=Year, y=response), size=1.5, alpha=0.5) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=1, color= col.model) +
		scale_x_continuous(limits=c(1901,2010), breaks=c(1925,1975)) +
		# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ": ", m.order, " 1900-2010"), x="Year", y=var)
		)
dev.off()


# summary(mod.out$ci.terms)

pdf(file.path(fig.dir, paste0("GAMM_DriverEffects_Site_", m.order, "_", v, ".pdf")))
print(
 	ggplot(data=mod.out$ci.terms[,]) + facet_wrap(Scale ~ Effect, scales="free", ncol=3) + theme_bw() +		
 		geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Site), alpha=0.5) +
		geom_line(aes(x=x, y=mean, color=Site), size=2) +
		geom_hline(yintercept=0, linetype="dashed") +
		# scale_color_manual(values=c("red2", "blue", "green3")) +
		# scale_fill_manual(values=c("red2", "blue", "green3")) +
		labs(title=paste0("Driver Effects: ",m.order), y="Effect Size")  +
		theme(axis.text.x=element_text(angle=0, color="black", size=rel(0.5)), axis.text.y=element_text(color="black", size=rel(0.5)), axis.title.x=element_text(face="bold", size=rel(0.75)),  axis.title.y=element_text(face="bold", size=rel(0.75)))
)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_DriverTime_Site_", m.order, "_", v, "_0850-2010.pdf")))
print(
ggplot(data= mod.out$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= mod.out$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y=v, title=paste0("Driver Effects through Time: ",m.order, ", 0850-2010")) +
	scale_x_continuous(limits=c(850,2010), breaks=c(1300,1800)) +
	scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	theme_bw() # + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_DriverTime_Site_", m.order, "_", v, "_1900-2010.pdf")))
print(
ggplot(data= mod.out$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= mod.out$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y=v, title=paste0("Driver Effects through Time: ",m.order, ", 1900-2010")) +
	scale_x_continuous(limits=c(1900,2010), breaks=c(1925,1975)) +
	scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() # + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
)
dev.off()

pdf(file.path(fig.dir, paste0("GAMM_DriverTime_Site_", m.order, "_", v, "_1800-1900.pdf")))
print(
ggplot(data= mod.out$weights) + facet_grid(Scale~Site, scales="free") +
	geom_line(data= mod.out$data, aes(x=Year, y=response), color="gray50", alpha=0.5, size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHA", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PHO", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PHO", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PUN", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PUN", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PBL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PBL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PDL", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PDL", scale="t.250", size=2) +

	rgb.line2(df=mod.out, site="PMB", scale="t.001", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.10", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.50", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.100", size=2) +
	rgb.line2(df=mod.out, site="PMB", scale="t.250", size=2) +

	labs(x="Year", y=v, title=paste0("Driver Effects through Time: ",m.order, ", 1800-1900")) +
	scale_x_continuous(limits=c(1800,1900), breaks=c(1825,1875)) +
	scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1800 & mod.out$data$Year<=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	theme_bw() #+ theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
)
dev.off()

} # end var
} # end model

# save(mod.out.AGB.diff, file=file.path(out.dir, "mod.out.dAGB.Rdata"))
# save(mod.out.NPP, file=file.path(out.dir, "mod.out.NPP.Rdata"))


# ------------------------



