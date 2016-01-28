# ----------------------------------------
# Objective: Compare climate sensitivities of models fit for different periods of data:
#             1. 1990-2010 (modern data, no fading record issues)
#             2. 1900-2010 (empirical record, modern met data, scale of succession)
#             3.  850-2010 (scale of migrations & multiple succession, climate change, etc)
# Christy Rollinson, crollinson@gmail.com
# Date Created: 28 July 2015
# ----------------------------------------
#
# -------------------------
# Workflow
# -------------------------
# 1. Set up Data
#    a. Ecosystem model output
#    b. Tree Ring NPP products
#    c. Raw Tree Ring widths
# 2. Run the gamms (with site intercept)
# 3. Bind Models into single list
# 4. Diagnostic Graphs
# -------------------------
# ----------------------------------------
rm(list=ls())
# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(parallel)
library(mgcv)
library(ggplot2); library(grid)
library(car)
# ----------------------------------------

# ----------------------------------------
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
start.yrs <- c(1985, 1901, 850)
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")
resolutions <- "t.001"
k=4

# wanted to do PFT, but we don't get even representation,
# so we're going to restrict it to PHO which has NPP data &
# the greatest composition consistency among models and data
sites <- c("PHO") 

# ----------------------------------------

# ----------------------------------------
# Set Directories & file paths
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/gamms"
fig.base="Figures/gamms"

# Source the gamm file
source('R/0_calculate.sensitivity_TPC.R', chdir = TRUE)

# Making sure the appropriate file paths exist
if(!dir.exists(dat.base)) dir.create(dat.base)
if(!dir.exists(fig.base)) dir.create(fig.base)

# Setting the data & figure directories
fig.dir <- file.path(fig.base, "Sensitivity_TempExtent")
dat.dir <- file.path(dat.base, "Sensitivity_TempExtent")

# Make sure the appropriate file paths are in place
if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# -------------------------------------------------------------------------------
# 1. Set up Data 
# -------------------------------------------------------------------------------
{
paleon.models <- list()
# ----------------------------------------
# 1.a. Load & set up Ecosystem Model Output first
# ----------------------------------------
{
# Define what our response variable will be
response <- "NPP"

# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData.Rdata"))

# Get rid of LINKAGES because it's weird & hasn't been updated
ecosys <- ecosys[!ecosys$Model=="linkages",]
summary(ecosys)


# ------------------
## Adding a biome classification so we can restrict analyses to 1 plot
#
#  # NOTE: This doesn't work with different temporal extents, so we have to go with
#          specific plots
# ------------------
ecosys$Fcomp_check <- rowSums(ecosys[,c("Evergreen", "Deciduous", "Grass")])
ecosys$PFT <- as.factor(
			  ifelse(ecosys$Evergreen/ecosys$Fcomp_check>=0.7, "Evergreen", 
              ifelse(ecosys$Deciduous/ecosys$Fcomp_check>=0.7, "Deciduous",
              ifelse(ecosys$Grass    /ecosys$Fcomp_check>=0.7, "Grass"    ,
              ifelse(rowSums(ecosys[,c("Evergreen", "Deciduous")])/ecosys$Fcomp_check>=0.7, "Mixed",
              ifelse(ecosys$Grass/ecosys$Fcomp_check>=0.3 & 
                     ecosys$Grass/ecosys$Fcomp_check<0.7, "Savanna",
              ifelse(!is.na(ecosys$Evergreen), "Other", NA
              )))))))
summary(ecosys[,c(1:18,(ncol(ecosys)-3):ncol(ecosys))])

summary(ecosys[!(ecosys$Model=="sibcasa") & ecosys$PFT=="Evergreen",c(1:18,(ncol(ecosys)-3):ncol(ecosys))])
summary(ecosys[!(ecosys$Model=="sibcasa") & ecosys$PFT=="Evergreen" & ecosys$Resolution=="t.001","PFT"])
# ------------------


for(m in unique(ecosys$Model)){

	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m, " ------"))

	for(y in start.yrs){
		# Taking the subsets of data we want in a single gam
		dat.subsets <- ecosys$Resolution == "t.001" & 
		          	   ecosys$Model      == m       & 
		               ecosys$Site      %in% sites  &
		               ecosys$Year       >= y

		# Skip models without the right data combo
		if(!length(which(dat.subsets)==TRUE)>0) next 

		# What will our spatio-temporal explanatory factor ("Time") be?
		if(!is.na(mean(ecosys[dat.subsets,"AGB"]))) time.mod="AGB" else time.mod="LAI"

		data.temp                  <- ecosys[dat.subsets, c("Model", "Model.Order", "Site", "Year")]
		data.temp$PlotID           <- ecosys[dat.subsets,"Site" ]
		data.temp$TreeID           <- as.factor(NA)
		data.temp$Y                <- ecosys[dat.subsets,response]
		data.temp$Time             <- ecosys[dat.subsets,time.mod]
		data.temp[,predictors.all] <- ecosys[dat.subsets, paste0(predictors.all, predictor.suffix)]
		data.temp$Resolution       <- ecosys[dat.subsets,"Resolution"]

		# Getting rid of NAs in predictors
		data.temp <- data.temp[complete.cases(data.temp[,c(predictors.all, "Y", "Time")]),]

		# Copy the response variable & some other things for the model
		paleon.models[[paste(m, ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp

	}
} # End Model Loop
} # End Model setup
# ----------------------------------------

# ----------------------------------------
# 1.b. Load & set up tree ring NPP
# ----------------------------------------
{
# Define what our response & time variables will be
response <- "ABI.area"
time.mod <- "AB.area"

# Load Tree ring NPP data
spp.npp <- read.csv(file.path("Data", "TreeRing_NPP_PlotSpecies.csv"))
summary(spp.npp)

# Getting PFT Fcomp per plot
pft.npp <- aggregate(spp.npp[,c("AB.area", "ABI.area", "tree.HA", "Fcomp")], by=spp.npp[,c("Site", "Site2", "PlotID", "PFT", "Year")], FUN=sum)
pft.npp[,c(paste0(predictors.all, predictor.suffix))] <- aggregate(spp.npp[,c(paste0(predictors.all, predictor.suffix))], by=spp.npp[,c("Site", "Site2", "PlotID", "PFT", "Year")], FUN=mean)[,c(paste0(predictors.all, predictor.suffix))]
summary(pft.npp)


# aggregate to total plot NPP (ABI.area)
plot.npp <- aggregate(spp.npp[,c(response, time.mod)], by=spp.npp[,c("Site", "Site2", "PlotID", "Year")], FUN=sum)
plot.npp[,c(paste0(predictors.all, predictor.suffix))] <- aggregate(spp.npp[,c(paste0(predictors.all, predictor.suffix))], by=spp.npp[,c("Site", "Site2", "PlotID", "Year")], FUN=mean)[,c(paste0(predictors.all, predictor.suffix))]
summary(plot.npp)

# Adding in the PFT factor (just for extra info at this point)
pft.vector <- vector()
for(i in 1:nrow(plot.npp)){
	p=plot.npp[i,"PlotID"]
	y=plot.npp[i,"Year"]

	f.decid <- pft.npp[pft.npp$PlotID==p & pft.npp$Year==y & pft.npp$PFT=="Deciduous", "Fcomp"]
	f.evg   <- pft.npp[pft.npp$PlotID==p & pft.npp$Year==y & pft.npp$PFT=="Evergreen", "Fcomp"]
		
	if(length(f.decid)==0) f.decid = 0
	if(length(f.evg  )==0) f.evg   = 0

	pft.type <- ifelse(f.evg>=0.7, "Evergreen", 
              ifelse(f.decid>=0.7, "Deciduous",
              ifelse((f.decid + f.evg)>=0.7, "Mixed", 
              ifelse(!is.na(ecosys$Evergreen), "Other", NA
              ))))

	pft.vector <- c(pft.vector, pft.type)
}
plot.npp$PFT <- as.factor(pft.vector)
summary(plot.npp)

# Subset a period where we're not worried about 
plot.npp <- plot.npp[complete.cases(plot.npp) & plot.npp$Year>=(2010-30),]
# summary(plot.npp)

for(y in start.yrs){
	# Taking the subsets of data we want in a single gam
	dat.subsets <- plot.npp$Site      %in% sites &
	               plot.npp$Year       >=  y

	# Skip data that doesn't meet the criteria for the model
	if(!length(which(dat.subsets)==TRUE)>0 | min(plot.npp$Year)>y) next 

	# If a variable is missing, just skip over this model for now
	if(!max(plot.npp[,response], na.rm=T)>0) next 

	# Add the data to paleon.models
	data.temp                  <- plot.npp[dat.subsets,c("Site", "PlotID", "Year", "PFT")]
	data.temp$Model            <- as.factor("TreeRingNPP")
	data.temp$Model.Order      <- as.factor("Tree Ring NPP")
	data.temp$TreeID           <- as.factor(NA)
	data.temp$Y                <- plot.npp[dat.subsets,response]
	data.temp$Time             <- plot.npp[dat.subsets,time.mod]
	data.temp[,predictors.all] <- plot.npp[dat.subsets, paste0(predictors.all, predictor.suffix)]
	data.temp$Resolution       <- as.factor("t.001")

	# Make sure everything is complete cases
	data.temp <- data.temp[complete.cases(data.temp[,c(predictors.all, "Y", "Time")]),]

	# Order everything the same way to make life easier
	data.temp <- data.temp[,names(paleon.models[[1]])]

	paleon.models[[paste("TreeRingNPP", ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp
} # End year loop

} # End Tree Ring NPP setup
# ----------------------------------------

# ----------------------------------------
# 1.c. Load & set up raw ring widths
# ----------------------------------------
{
response <- "RW"
time.mod <- "DBH"

tree.rings <- read.csv("Data/TreeRing_RingWidths.csv")
summary(tree.rings)

# subset only complete cases where we have met data & ring width measurements
tree.rings <- tree.rings[complete.cases(tree.rings[,c(response, paste0(predictors.all, predictor.suffix))]) & tree.rings$Resolution=="t.001",]
summary(tree.rings)


for(y in start.yrs){
	# Taking the subsets of data we want in a single gam
	dat.subsets <- tree.rings$Site      %in% sites &
	               tree.rings$Year       >=  y

	# Skip data that doesn't meet the criteria for the model
	if(!length(which(dat.subsets)==TRUE)>0 | min(tree.rings$Year)>y) next 

	# If a variable is missing, just skip over this model for now
	if(!max(tree.rings[,response], na.rm=T)>0) next 

	# Add the data to paleon.models
	data.temp                  <- tree.rings[dat.subsets,c("Site", "PlotID", "TreeID", "Year", "PFT")]
	data.temp$Model            <- as.factor("TreeRingRW")
	data.temp$Model.Order      <- as.factor("Tree Ring RW")
	data.temp$Y                <- tree.rings[dat.subsets,response]
	data.temp$Time             <- log(tree.rings[dat.subsets,time.mod])
	data.temp[,predictors.all] <- tree.rings[dat.subsets, paste0(predictors.all, predictor.suffix)]
	data.temp$Resolution       <- as.factor("t.001")

	# Make sure everything is complete cases
	data.temp <- data.temp[complete.cases(data.temp[,c(predictors.all, "Y", "Time")]),]

	# Order everything the same way to make life easier
	data.temp <- data.temp[,names(paleon.models[[1]])]

	paleon.models[[paste("TreeRingRW", ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp
} # End year loop

} # End Ring Width setups
# ----------------------------------------
}
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# 2. Run the gamms 
# -------------------------------------------------------------------------------
cores.use <- min(12, length(paleon.models))
# cores.use <- length(paleon.models)

models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=cores.use, k=k, predictors.all=predictors.all, PFT=F)
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# 3. Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------------------------------------
{
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], "TempExtent", sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], "TempExtent", sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_TempExtent.Rdata"))
}
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# 4. Diagnostic Graphs
# -------------------------------------------------------------------------------
{
m.order <- unique(mod.out$data$Model.Order)
col.model <- c(paste(model.colors[model.colors$Model.Order %in% m.order,"color"]), "black", "gray30")

pdf(file.path(fig.dir, "GAMM_ModelFit_TempExtent.pdf"))
{
for(e in unique(mod.out$data$Extent)){
print(	
ggplot(data=mod.out$ci.response[!substr(mod.out$ci.response$Model, 1, 8)=="TreeRing",]) + facet_wrap(~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Extent==e & !substr(mod.out$data$Model, 1, 8)=="TreeRing",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(data=mod.out$ci.response[mod.out$ci.response$Extent==e & !substr(mod.out$ci.response$Model, 1, 8)=="TreeRing",], aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(data=mod.out$ci.response[mod.out$ci.response$Extent==e & !substr(mod.out$ci.response$Model, 1, 8)=="TreeRing",], aes(x=Year, y=mean, color=Model), size=0.35) +
	# scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Extent", response, sep=" - "), x="Year", y=response)
)
}

print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Model=="TreeRingNPP",]) + facet_wrap(~ PlotID, scales="free_x") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Model=="TreeRingNPP",], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("TempExtent", response, sep=" - "), x="Year", y=response)
)

trees.subset <- c("30-01", "30-02", "30-03", "30-04", "30-21", "30-22")
print(	
ggplot(data=mod.out$ci.response[mod.out$ci.response$Model=="TreeRingRW" & mod.out$ci.response$PlotID=="ME029" & mod.out$ci.response$TreeID %in% trees.subset,]) + 
	facet_grid(TreeID~ Extent, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Model=="TreeRingRW" & mod.out$data$PlotID=="ME029" & mod.out$data$TreeID %in% trees.subset,], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	# scale_x_continuous(limits=c(1900,2010)) +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("TempExtent", response, sep=" - "), x="Year", y=response)
)
}
dev.off()

mod.out$ci.terms$x <- as.numeric(paste(mod.out$ci.terms$x))
summary(mod.out$ci.terms)

pdf(file.path(fig.dir, "GAMM_DriverSensitivity_TempExtent.pdf"))
{
m.order <- unique(mod.out$data[,"Model.Order"])
col.model <- c(paste(model.colors[model.colors$Model.Order %in% m.order,"color"]), "black", "gray30")
print(
ggplot(data=mod.out$ci.terms[,]) + facet_grid(Extent ~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (not Relativized)"), y=paste0("NPP Contribution")) )

}
dev.off()
}
## -------------------------------------------------------------------------------


