# ----------------------------------------
# Temporal Scaling Analyses
# Changes in Strength of Interactions with Temporal Scale
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Across models, does the temporal extent change the temp, precip, & CO2 response curves?
#
# -------------------------
#
# -------------------------
# Data/Results Generation:
# -------------------------
# Inputs from 2a_process_drivers_by_extent.R
#
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# 
#
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(nlme)
library(car)
library(zoo)
# library(mvtnorm)
# library(MCMCpack)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
# setwd("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
path.data <- "Data/"
# fig.dir <- "Figures"
# ----------------------------------------

----------------------------------------
Load data files & function scripts
----------------------------------------
#------------
# Putting all the models into 1 object
#------------
models <- dir(file.path(path.data, "gamms_byModel"))

# By Extent
for(m in 1:length(models)){
	path.temp <- file.path(path.data, "gamms_byModel", models[m], "byExtent")
	dir.temp <- dir(path.temp, "NPP")
	load(file.path(path.temp, dir.temp))

	model       <- unique(mod.out$data$Model)
	model.order <- unique(mod.out$data$Model.Order)

	mod.out$weights$Model           <- model
	mod.out$weights$Model.Order     <- model.order
	mod.out$ci.response$Model       <- model
	mod.out$ci.response$Model.Order <- model.order
	mod.out$ci.terms$Model          <- model
	mod.out$ci.terms$Model.Order    <- model.order

	if(m == 1){
		mod.extent <- mod.out[1:4]
	} else { 
		mod.extent$data           <- rbind(mod.extent$data, mod.out$data)
		mod.extent$weights        <- rbind(mod.extent$weights, mod.out$weights)
		mod.extent$ci.response    <- rbind(mod.extent$ci.response, mod.out$ci.response)
		mod.extent$ci.terms       <- rbind(mod.extent$ci.terms, mod.out$ci.terms)
	}
}
save(mod.extent, file=file.path(path.data, "GAMM_byExtent_AllModels.RData"))

# By Site
for(m in 1:length(models)){
	path.temp <- file.path(path.data, "gamms_byModel", models[m], "bySite")
	dir.temp <- dir(path.temp, "NPP")
	load(file.path(path.temp, dir.temp))

	model       <- unique(mod.out$data$Model)
	model.order <- unique(mod.out$data$Model.Order)

	mod.out$weights$Model           <- model
	mod.out$weights$Model.Order     <- model.order
	mod.out$ci.response$Model       <- model
	mod.out$ci.response$Model.Order <- model.order
	mod.out$ci.terms$Model          <- model
	mod.out$ci.terms$Model.Order    <- model.order

	if(m == 1){
		mod.site <- mod.out[1:4]
	} else { 
		mod.site$data           <- rbind(mod.site$data, mod.out$data)
		mod.site$weights        <- rbind(mod.site$weights, mod.out$weights)
		mod.site$ci.response    <- rbind(mod.site$ci.response, mod.out$ci.response)
		mod.site$ci.terms       <- rbind(mod.site$ci.terms, mod.out$ci.terms)
	}
}
save(mod.site, file=file.path(path.data, "GAMM_bySite_AllModels.RData"))

# #------------

# Load gamm by extent outputs
load(file.path(path.data, "GAMM_byExtent_AllModels.RData"))
load(file.path(path.data, "GAMM_bySite_AllModels.RData"))

# Read in model color scheme
model.colors <- read.csv("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors

# ----------------------------------------


# ----------------------------------------
# Exploratory graphing & analyses
# ----------------------------------------
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% unique(mod.extent$ci.terms$Model.Order) & !model.colors$Model.Order=="CLM-BGC" & !model.colors$Model.Order=="LINKAGES","color"])

# ------------------
# Comparing effect of extent, colors by model
# ------------------
ggplot(data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Scale=="t.001" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",]) + facet_grid(Extent~Effect, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use) +
	scale_fill_manual(values=colors.use) +
	ggtitle("NPP Responses by Temporal Extent, Grain = Annual") +
	theme_bw()


ggplot(data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Scale=="t.001" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",]) + facet_wrap(Extent~Effect, scales="free", ncol=3) +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use) +
	scale_fill_manual(values=colors.use) +
	ggtitle("NPP Responses by Temporal Extent, Grain = Annual") +
	theme_bw()

# Basic mixed linear effects model; may want to repeat on bootstrapped data used to generate the CI
lm.extent.grain <- lme(mean ~ Effect*Extent*Scale - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA",])
anova(lm.extent.grain)

lm.extent.t001 <- lme(mean ~ Effect*Extent - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Scale=="t.001" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.extent.t001)

lm.extent.t10 <- lme(mean ~ Effect*Extent - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Scale=="t.10" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.extent.t10)

lm.extent.t100 <- lme(mean ~ Effect*Extent - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Scale=="t.100" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.extent.t100)
# ------------------


# ------------------
# Comparing effect of resolution, color by model
# ------------------
ggplot(data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Extent=="850-2010" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",]) + facet_grid(Scale~Effect, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use) +
	scale_fill_manual(values=colors.use) +
	ggtitle("NPP Responses by Temporal Scale, Extent = 850-2010") +
	theme_bw()

ggplot(data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Extent=="850-2010" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",]) + facet_wrap(Scale~Effect, scales="free", ncol=3) +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use) +
	scale_fill_manual(values=colors.use) +
	ggtitle("NPP Responses by Temporal Scale, Extent = 850-2010") +
	theme_bw()

# Basic mixed linear effects model; may want to repeat on bootstrapped data used to generate the CI
lm.extent.grain <- lme(mean ~ Effect*Extent*Scale - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.extent.grain)

lm.res.850 <- lme(mean ~ Effect*Scale - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Extent=="850-2010" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.res.850)

lm.res.1900 <- lme(mean ~ Effect*Scale - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Extent=="1900-2010" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.res.1900)

lm.res.1990 <- lme(mean ~ Effect*Scale - 1, random=list(Model=~1), data=mod.extent$ci.terms[mod.extent$ci.terms$Site=="PHA" & mod.extent$ci.terms$Extent=="1990-2010" & !mod.extent$ci.terms$Model.Order=="CLM-BGC" & !mod.extent$ci.terms$Model.Order=="LINKAGES",])
anova(lm.res.1990)
# ------------------


# ------------------
# Comparing effect of site and resolution, color by model
# ------------------
colors.use2 <- as.vector(model.colors[model.colors$Model.Order %in% unique(mod.extent$ci.terms$Model.Order) & !model.colors$Model.Order=="CLM-BGC" & !model.colors$Model.Order=="LINKAGES" & !model.colors$Model.Order=="ED2-LU","color"])

ggplot(data= mod.site$ci.terms[mod.site$ci.terms$Scale=="t.001" & !mod.site$ci.terms$Model.Order=="CLM-BGC" & !mod.site$ci.terms$Model.Order=="LINKAGES" & !mod.site$ci.terms$Model.Order=="ED2-LU" ,]) + 
	facet_grid(Site~Effect, scales="free") +
	# facet_grid(Site~Scale, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order), size=0.75) +
	geom_hline(aes(yintercept=0), linetype="dashed") +
	scale_color_manual(values=colors.use2) +
	scale_fill_manual(values=colors.use2) +
	labs(y="Effect Size", x="") +
	ggtitle("NPP Responses by Effect, Grain = Annual") +
	theme_bw()


ggplot(data= mod.site$ci.terms[mod.site$ci.terms$Scale=="t.10" & !mod.site$ci.terms$Model.Order=="CLM-BGC" & !mod.site$ci.terms$Model.Order=="LINKAGES" & !mod.site$ci.terms$Model.Order=="ED2-LU" ,]) + 
	facet_grid(Site~Effect, scales="free") +
	# facet_grid(Site~Scale, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use2) +
	scale_fill_manual(values=colors.use2) +
	labs(y="Effect Size", x="") +
	ggtitle("NPP Responses by Effect, Grain = Decadal") +
	theme_bw()


ggplot(data= mod.site$ci.terms[mod.site$ci.terms$Scale=="t.100" & !mod.site$ci.terms$Model.Order=="CLM-BGC" & !mod.site$ci.terms$Model.Order=="LINKAGES" & !mod.site$ci.terms$Model.Order=="ED2-LU" ,]) + 
	facet_grid(Site~Effect, scales="free") +
	# facet_grid(Site~Scale, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use2) +
	scale_fill_manual(values=colors.use2) +
	labs(y="Effect Size", x="") +
	ggtitle("NPP Responses by Effect, Grain = Centennial") +
	theme_bw()


ggplot(data= mod.site$ci.terms[mod.site$ci.terms$Effect=="Temp" & ! mod.site$ci.terms$Model.Order=="CLM-BGC" & !mod.site$ci.terms$Model.Order=="LINKAGES" & !mod.site$ci.terms$Model.Order=="ED2-LU",]) + 
	facet_grid(Scale~Site, scales="free") +
	# facet_grid(Site~Scale, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use2) +
	scale_fill_manual(values=colors.use2) +
	ggtitle("NPP Responses to Temperature by Grain") +
	theme_bw()

ggplot(data= mod.site$ci.terms[mod.site$ci.terms$Effect=="Temp" & ! mod.site$ci.terms$Model.Order=="CLM-BGC" & !mod.site$ci.terms$Model.Order=="LINKAGES" & !mod.site$ci.terms$Model.Order=="ED2-LU",]) + 
	facet_wrap(Scale~Site, scales="free", ncol=6) +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.3) +
	geom_line(aes(x=x, y=mean, color=Model.Order)) +
	scale_color_manual(values=colors.use2) +
	scale_fill_manual(values=colors.use2) +
	ggtitle("NPP Responses by Temporal Extent, Grain = Annual") +
	theme_bw()

# ------------------
	
# ----------------------------------------
	
# ----------------------------------------
