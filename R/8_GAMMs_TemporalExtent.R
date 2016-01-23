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
# Note: To really isolate data-model comparisons to temporal extent, we need to try to
#       minimize differences among models & sites.  I had originally tried to restrict 
#       things to a single PFT, but that doesn't pan out.  Instead we'll just restrict
#       the temporal extent analysis to the two sites with best data: Harvard & Howland 
#
# 1. Models
#    a. Load model data files & function scripts
#    b. Settings for the rest of this section
#    c. Setting up to run gamms in parallel
#    d. Run the gamms (with site intercept)
#    e. Bind Models into single list
#    f. Diagnostic Graphs
# 2. Tree Rings NPP 
#    a. Load model data files & function scripts
#    b. Settings for the rest of this section
#    c. Setting up to run gamms 
#    d. Run the gamms (with site & intercepts); save the data
#    e. Bind Models into single list
#    f. Diagnostic Graphs
# 3. Tree Rings: BAI & RWI
#    a. Load model data files & function scripts
#    b. Settings for the rest of this section
#    c. Setting up to run gamms in parallel
#    d. Run the gamms (with site & intercepts); save the data
#    e. Bind Models into single list
#    f. Diagnostic Graphs
# -------------------------
# ----------------------------------------

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
sec2yr    <- 1*60*60*24*365
start.yrs <- c(1985, 1901, 850)
# ----------------------------------------

# ----------------------------------------
# Set Directories & file paths
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
dat.base="Data/gamms"
fig.base="Figures/gamms"

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
# 1. Models
# -------------------------------------------------------------------------------
# ----------------------------------------
# 1.a. Load model data files & function scripts
# ----------------------------------------
# Ecosys file = organized, post-processed m.name outputs
#	generated with 1_generate_ecosys.R
load(file.path("Data", "EcosysData.Rdata"))
summary(ecosys)
model.colors

source('R/0_calculate.sensitivity_TPC.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)

# Read in model color scheme
model.colors
# ----------------------------------------


# -------------------------------------------------
# 1.b. Settings for the rest of this section
# -------------------------------------------------
# Get rid of LINKAGES because its acting kind of funny
ecosys <- ecosys[!ecosys$Model=="linkages",]

# Setting up a loop for 1 m.name, 1 temporal scale
sites       <- unique(ecosys$Site)
model.name  <- unique(ecosys$Model)
model.order <- unique(ecosys$Model.Order)
resolutions <- c("t.001") # Note: Big models can't handle t.100 at the site level because there aren't enough data points
response <- "NPP"
predictors.all <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")
k=4
e=1	
# -------------------------------------------------

# -------------------------------------------------
# 1.c. Setting up the data and putting it in a list to run the gamms in parallel
# ------------------------------------------------
paleon.models <- list()
for(m in 1:length(model.name)){
	m.name  <- model.name[m]
	m.order <- model.order[m]

	print("-------------------------------------")
	print(paste0("------ Processing Model: ", m.order, " ------"))
	
	for(y in start.yrs){
		# Taking the subsets of data we want in a single gam
		dat.subsets <- ecosys$Resolution ==  resolutions     & 
		               ecosys$Model      ==  m.name          &
		               ecosys$Site      %in% c("PHO", "PHA") &
		               ecosys$Year       >=  y

		if(!length(which(dat.subsets)==TRUE)>0) next 

		data.temp <- ecosys[dat.subsets, c("Model", "Model.Order", "Site", "Year", response, paste0(predictors.all, predictor.suffix))]

		# renaming the met var
		names(data.temp)[(ncol(data.temp)-length(predictors.all)+1):ncol(data.temp)] <- predictors.all 
	
		# If a variable is missing, just skip over this model for now
		if(!max(data.temp[,response], na.rm=T)>0) next 

		# Making a note of the resolution
		data.temp$Resolution <- as.factor(resolutions)

		# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
		data.temp <- data.temp[complete.cases(data.temp[,response]),]

		data.temp$Y <- data.temp[,response]

		paleon.models[[paste(m.name, ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp

	} # End extent loop
} # End Model Loop
# --------------------------------


# -------------------------------------------------
# 1.d. Run the gamms -- WITH site intercept
# -------------------------------------------------
cores.use <- min(12, length(paleon.models))
# cores.use <- length(paleon.models)

models.base <- mclapply(paleon.models, paleon.gams.models, mc.cores=cores.use, k=k, predictors.all=predictors.all, PFT=F)
# -------------------------------------------------

# -------------------------------------------------
# 1.e. Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], "TempExt", sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], "TempExt", sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_TempExt_Models.Rdata"))
# -------------------------------------------------


# -------------------------------------------------
# 1.f. Diagnostic Graphs
# -------------------------------------------------
m.order <- unique(mod.out$data$Model.Order)
col.model <- model.colors[model.colors$Model.Order %in% m.order,"color"]

pdf(file.path(fig.dir, "GAMM_ModelFit_TempExt_Models.pdf"))
# print(
# ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~Model, scales="free") + theme_bw() +
 	# geom_line(data= mod.out$data[,], aes(x=Year, y=Y), alpha=0.5) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(aes(x=Year, y=mean, color=Model), size=0.35) +
	# # scale_x_continuous(limits=c(850,2010)) +
	# # scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99),na.rm=T)) +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	# labs(title=paste("Baseline, No Site Effect", response, sep=" - "), x="Year", y=response)
# )
for(e in unique(mod.out$data$Extent)){
print(	
ggplot(data=mod.out$ci.response[,]) + facet_grid(Site~ Model, scales="free") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Extent==e,], aes(x=Year, y=Y), alpha=0.5) +
	geom_ribbon(data=mod.out$ci.response[mod.out$ci.response$Extent==e,], aes(x=Year, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(data=mod.out$ci.response[mod.out$ci.response$Extent==e,], aes(x=Year, y=mean, color=Model), size=0.35) +
	scale_x_continuous() +
	# scale_y_continuous(limits=quantile(mod.out$data[mod.out$data$Year>=1900,"response"], c(0.01, 0.99),na.rm=T)) +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste("Baseline, No Site Effect", response, sep=" - "), x="Year", y=response)
)
}
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverSensitivity_TempExt_Models.pdf"))
print(
ggplot(data=mod.out$ci.terms[,]) + facet_grid(Extent ~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	scale_fill_manual(values=paste(col.model)) +
	scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (not Relativized)"), y=paste0("NPP Contribution")) # +
)
dev.off()
# -------------------------------------------------

# Clear the memory!
rm(mod.out, models.base, ecosys, dat.mod)
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# 2. Tree Rings NPP
# -------------------------------------------------------------------------------
# ----------------------------------------
# 2.a. Load & format data files & function scripts
# ----------------------------------------
source('R/0_calculate.sensitivity_TPC_TreeRingNPP.R', chdir = TRUE)
source('R/0_GAMM_Plots.R', chdir = TRUE)

# What Climate predictors we're interested in
predictors.all   <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")

# Load Tree ring NPP data
spp.npp <- read.csv(file.path("Data", "TreeRing_NPP_PlotSpecies.csv"))
summary(spp.npp)

# Getting PFT Fcomp per plot
pft.npp <- aggregate(spp.npp[,c("AB.area", "ABI.area", "tree.HA", "Fcomp")], by=spp.npp[,c("Site", "PlotID", "PFT", "Year")], FUN=sum)
pft.npp[,c(paste0(predictors.all, predictor.suffix))] <- aggregate(spp.npp[,c(paste0(predictors.all, predictor.suffix))], by=spp.npp[,c("Site", "PlotID", "PFT", "Year")], FUN=mean)[,c(paste0(predictors.all, predictor.suffix))]
summary(pft.npp)

# aggregate to total plot NPP (ABI.area)
plot.npp <- aggregate(spp.npp[,c("AB.area", "ABI.area", "tree.HA", "Fcomp")], by=spp.npp[,c("Site", "PlotID", "Year")], FUN=sum)
plot.npp[,c(paste0(predictors.all, predictor.suffix))] <- aggregate(spp.npp[,c(paste0(predictors.all, predictor.suffix))], by=spp.npp[,c("Site", "PlotID", "Year")], FUN=mean)[,c(paste0(predictors.all, predictor.suffix))]
summary(plot.npp)

plot.npp$Model       <- as.factor("plotNPP")
plot.npp$Model.Order <- as.factor("Plot NPP")
summary(plot.npp)

# subset only complete cases where we have met data and data from the past 30 years
plot.npp <- plot.npp[complete.cases(plot.npp) & plot.npp$Year>=(2010-30),]
summary(plot.npp)
# -------------------------------------------------

# -------------------------------------------------
# 2.b. Settings for the rest of this section
# -------------------------------------------------
# Setting up a loop for 1 m.name, 1 temporal scale
response    <- "ABI.area"
resolutions <- "t.001" # Note: Big models can't handle t.100 at the site level because there aren't enough data points
m.name      <- unique(plot.npp$Model)
k=4
e=1	
# -------------------------------------------------

# -------------------------------------------------
# 2.c. Setting up the data for the model
# -------------------------------------------------
dat.mod <- list()
for(y in start.yrs){
	# Taking the subsets of data we want in a single gam
	dat.subsets <- plot.npp$Site      %in% c("PHO", "PHA") &
	               plot.npp$Year       >=  y

	if(!length(which(dat.subsets)==TRUE)>0 | min(plot.npp$Year)>y) next 

	data.temp <- plot.npp[dat.subsets, c("Model", "Model.Order", "Site", "PlotID", "Year", response, paste0(predictors.all, predictor.suffix))]

	# renaming the met var
	names(data.temp)[(ncol(data.temp)-length(predictors.all)+1):ncol(data.temp)] <- predictors.all 

	# If a variable is missing, just skip over this model for now
	if(!max(data.temp[,response], na.rm=T)>0) next 

	# Making a note of the resolution
	data.temp$Resolution <- as.factor(resolutions)

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,response]),]

	data.temp$Y <- data.temp[,response]

	dat.mod[[paste(m.name, ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp

} # End extent loop
summary(dat.mod)
# -------------------------------------------------

# -------------------------------------------------
# 2.d. Run the gamms (with site & intercepts); save the data
# -------------------------------------------------
cores.use <- min(12, length(dat.mod))

models.base <- mclapply(dat.mod, paleon.gams.models, mc.cores=cores.use, k=k, predictors.all=predictors.all, PFT=F)
# -------------------------------------------------


# -------------------------------------------------
# 1.e. Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_TempExt_TreeRingNPP.Rdata"))
# -------------------------------------------------

# -------------------------------------------------
# 2.f. Diagnostic Graphs
# -------------------------------------------------
pdf(file.path(fig.dir, "GAMM_ModelFit_TempExt_TreeRingNPP.pdf"))
print(
ggplot(data=mod.out$data) + facet_wrap(PlotID~Model, scales="fixed") + theme_bw() +
 	geom_line(data= mod.out$data[mod.out$data$Model=="plotNPP",], aes(x=Year, y=Y), alpha=0.5, size=1) +
	geom_ribbon(data=mod.out$ci.response[mod.out$ci.response$Model=="plotNPP",], aes(x=Year, ymin=lwr, ymax=upr, fill=PlotID), alpha=0.5) +
	geom_line(data=mod.out$ci.response[mod.out$ci.response$Model=="plotNPP",], aes(x=Year, y=mean, color=PlotID), size=0.35) +
	labs(title=paste("Baseline, No Site Effect", response, sep=" - "), x="Year", y=response)
)
dev.off()

pdf(file.path(fig.dir, "GAMM_DriverSensitivity_TempExt_TreeRingNPP.pdf"))
print(
ggplot(data=mod.out$ci.terms[!mod.out$ci.terms$Effect=="PFT", ]) + facet_grid(Extent~ Effect, scales="free_x") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (not Relativized)"), y=paste0("NPP Contribution")) # +
)
dev.off()
# -------------------------------------------------
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# 3. Tree Rings BAI
# -------------------------------------------------------------------------------
# -------------------------------------------------
# 3.a. Load model data files & function scripts
# -------------------------------------------------
source('R/0_calculate.sensitivity_TPC_TreeRings.R', chdir = TRUE)

# What Climate predictors we're interested in
predictors.all   <- c("tair", "precipf", "CO2")
predictor.suffix <- c(".gs")
response         <- "RWI"

# Load tree ring width data
tree.rings <- read.csv("Data/TreeRing_RingWidths.csv")
summary(tree.rings)

tree.rings$Model       <- as.factor("TreeRingRWI")
tree.rings$Model.Order <- as.factor("Tree Ring RWI")

# subset only complete cases where we have met data & ring width measurements
tree.rings <- tree.rings[complete.cases(tree.rings[,c("BAI", "RWI", paste0(predictors.all, predictor.suffix))]) & tree.rings$Resolution=="t.001",]
summary(tree.rings)
# -------------------------------------------------

# -------------------------------------------------
# 3.b. Settings for the rest of this section
# -------------------------------------------------
k=4
e=1	
# -------------------------------------------------

# -------------------------------------------------
# 3.c. Setting up the data and putting it in a list to run the gamms in parallel
# -------------------------------------------------
dat.mod <- list()
for(y in start.yrs){
	m.name=unique(tree.rings$Model)
	# Taking the subsets of data we want in a single gam
	dat.subsets <- tree.rings$Site      %in% c("PHO", "PHA") &
	               tree.rings$Year       >=  y

	if(!length(which(dat.subsets)==TRUE)>0 | min(tree.rings$Year)>y) next 

	data.temp <- tree.rings[dat.subsets, c("Model", "Model.Order", "Site", "PlotID", "TreeID", "Year", response, paste0(predictors.all, predictor.suffix))]

	# renaming the met var
	names(data.temp)[(ncol(data.temp)-length(predictors.all)+1):ncol(data.temp)] <- predictors.all 

	# If a variable is missing, just skip over this model for now
	if(!max(data.temp[,response], na.rm=T)>0) next 

	# Making a note of the resolution
	data.temp$Resolution  <- as.factor(resolutions)

	# Getting rid of NAs; note: this has to happen AFTER extent definition otherwise scale & extent are compounded
	data.temp <- data.temp[complete.cases(data.temp[,response]),]

	data.temp$Y <- data.temp[,response]

	dat.mod[[paste(m.name, ifelse(nchar(y)==3, paste0(0, y), y), sep="_")]] <- data.temp

} # End extent loop
summary(dat.mod)
# -------------------------------------------------

# -------------------------------------------------
# 3.d. Run the gamms -- WITH site intercept
# -------------------------------------------------
# # This isn't working well in parallel, so we'll run it 1 at a time
# cores.use <- min(12, length(dat.mod))
# models.base <- mclapply(dat.mod, paleon.gams.models, mc.cores=cores.use, k=k, predictors.all=predictors.all, PFT=F)

models.base <- list()
models.base[[names(dat.mod)[1]]] <- paleon.gams.models(data=dat.mod[[1]], k=k, predictors.all=predictors.all, PFT=F)
models.base[[names(dat.mod)[2]]] <- paleon.gams.models(data=dat.mod[[2]], k=k, predictors.all=predictors.all, PFT=F)

save(models.base, file=file.path(dat.dir, "gamm_TempExt_TreeRings.Rdata"))
# -------------------------------------------------
# -------------------------------------------------

# -------------------------------------------------
# 3.e. Bind Models together to put them in a single object to make them easier to work with
# -------------------------------------------------
for(i in 1:length(models.base)){
	if(i==1) {
		mod.out <- list()
		mod.out$data         <- models.base[[i]]$data
		mod.out$weights      <- models.base[[i]]$weights
		mod.out$ci.response  <- models.base[[i]]$ci.response
		mod.out$sim.response <- models.base[[i]]$sim.response
		mod.out$ci.terms     <- models.base[[i]]$ci.terms
		mod.out$sim.terms    <- models.base[[i]]$sim.terms
		mod.out[[paste("gamm", names(models.base)[i], "PFT", sep=".")]] <- models.base[[i]]$gamm
	} else {
		mod.out$data         <- rbind(mod.out$data,         models.base[[i]]$data)
		mod.out$weights      <- rbind(mod.out$weights,      models.base[[i]]$weights)
		mod.out$ci.response  <- rbind(mod.out$ci.response,  models.base[[i]]$ci.response)
		mod.out$sim.response <- rbind(mod.out$sim.response, models.base[[i]]$sim.response)
		mod.out$ci.terms     <- rbind(mod.out$ci.terms,     models.base[[i]]$ci.terms)
		mod.out$sim.terms    <- rbind(mod.out$sim.terms,    models.base[[i]]$sim.terms)
		mod.out[[paste("gamm", names(models.base)[i], "PFT", sep=".")]] <- models.base[[i]]$gamm
	}
}

save(mod.out, file=file.path(dat.dir, "gamm_TempExt_TreeRings.Rdata"))
# -------------------------------------------------

# -------------------------------------------------
# 3.f. Diagnostic Graphs
# -------------------------------------------------
pdf(file.path(fig.dir, "GAMM_DriverSensitivity_TempExt_TreeRings.pdf"))
print(
ggplot(data=mod.out$ci.terms[, ]) + facet_grid(.~ Effect, scales="free") + theme_bw() +		
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Extent), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Extent), size=2) +
	geom_hline(yintercept=0, linetype="dashed") +
	# scale_fill_manual(values=paste(col.model)) +
	# scale_color_manual(values=paste(col.model)) +		
	labs(title=paste0("Driver Sensitivity (not Relativized)"), y=paste0("NPP Contribution")) # +
)
dev.off()
# -------------------------------------------------
# -------------------------------------------------------------------------------
