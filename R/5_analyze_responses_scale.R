# ----------------------------------------
# Temporal Scaling Analyses
# Looking at whether temporal resolution or extent has greater impact on model responses to drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Which is the more important consideration in understanding ecosystem variation in response to drivers
#           of change: temporal resolution, or temporal scale?
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Individual model response curves to driver
#    -- generated with 2c_process_drivers_all_drivers.R & 2e_process_drivers_all_drivers_byExtent.R
#    -- NOTE: Taking growing season analysis because it had a slightly higher explanatory power
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) Determine & compare relative weight of individual drivers across model
#    -- pay special attention to the role of Tair, Precip & CO2
# B) Determine the change in driver response with change in scale
#    -- relativize to the annual res, 850-2010 extent curves
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
# setwd("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
path.data <- "Data"
in.base <- "Data/gamms"
in.res  <- "AllDrivers_GS_byResolution"
in.ext  <- "AllDrivers_GS_byExtent"
out.dir <- "Data/analysis_response_scale"
fig.dir <- "Figures/analysis_response_scale"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData_Raw.Rdata"))

# Figure out what models we have to work with
models <- dir(in.base)
models <- models[2:length(models)] # CLM-BGC was being buggy so we had to stop

for(i in 1:length(models)){
	# loop through by resolution
	fmod <- dir(file.path(in.base, models[i], in.res))
	load(file.path(in.base, models[i], in.res, fmod))
	
	if(i==1) {
		ci.terms  <- mod.out$ci.terms[mod.out$ci.terms$Extent=="850-2010-NA",]
		# sim.terms <- mod.out$sim.terms
	} else {
		ci.terms  <- rbind(ci.terms, mod.out$ci.terms )
		# sim.terms <- rbind(sim.terms, mod.out$sim.terms)
	}

	# loop through by extent
	fmod <- dir(file.path(in.base, models[i], in.ext))
	load(file.path(in.base, models[i], in.ext, fmod))
	# if(i==1) {
		# ci.terms  <- mod.out$ci.terms 
		# # sim.terms <- mod.out$sim.terms
	# } else {
		# Note: because we're lumping everything together, let's not mess with reiterating the base level
		ci.terms  <- rbind(ci.terms, mod.out$ci.terms[!(mod.out$ci.terms$Scale=="t.001" & substr(mod.out$ci.terms$Extent,1,3)=="850"),] )
		# # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
	# }

	# Clear the mod.out to save space
	rm(mod.out)
}

ci.terms$Extent <- as.factor(substr(ci.terms$Extent, 1, nchar(as.character(ci.terms$Extent))-3))
ci.terms$Extent <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
summary(ci.terms$Extent)


# Write the files to csv so I don't have to mess with loading large .Rdata files again if I don't have to
# write.csv(ci.terms.res,  file.path(out.dir, "Driver_Responses_CI_Resolution.csv"), row.names=F)
# write.csv(sim.terms.res, file.path(out.dir, "Driver_Responses_Sims_Resolution.csv"), row.names=F)
# write.csv(ci.terms.ext,  file.path(out.dir, "Driver_Responses_CI_Extent.csv"), row.names=F)
# write.csv(sim.terms.ext, file.path(out.dir, "Driver_Responses_Sims_Extent.csv"), row.names=F)
# ----------------------------------------

# ----------------------------------------
# Standardize driver responses to the mean model NPP
#
# NOTE: Right now, just going to work with the mean responses & ignore the CIs for the sake of simplicity
# ----------------------------------------
# Across all scales (res + extent) finding 
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Scale"])){
		for(e in unique(ci.terms[ci.terms$Model==m, "Extent"])){
			# Find the start year for the extent
			yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			# Find the NPP to relativize each set off of
			npp <- mean(ecosys[ecosys$Model==m & ecosys$Scale==r & ecosys$Year>=yr, "NPP"], na.rm=T)			
			
			# Finding the percent change in NPP relative to the mean for that particular scale
			ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Scale==r & ci.terms$Extent==e,"upr"])/npp
		}
	}
}
summary(ci.terms)

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(ecosys[ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])


pdf(file.path(fig.dir, "NPP_RelChange_Resolution.pdf"))
ggplot(data= ci.terms[ci.terms$Extent=="0850-2010",]) + 
	facet_grid(Scale~Effect, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_Extent.pdf"))
ggplot(data= ci.terms[ci.terms$Scale=="t.001",]) + 
	facet_grid(Extent~Effect, scales="free") +
	geom_ribbon(aes(x=x, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_ResAnn_Ext0850.pdf"))
ggplot(data=ci.terms[ci.terms$Scale=="t.001" & ci.terms$Extent=="0850-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()


# ----------------------------------------

# ----------------------------------------
# Standardize changes in scale to be deviation from the annual 850-2010 response
#
# Note: because in what we already ran, we only did the CI & sims over the range of the condition in a given model, 
#       we actually can't just subtract. (booo!)  This means that we need to go back through, load the appropriate gamms 
#       from the rdata file & re-predict the models over the same range of values.
# ----------------------------------------

# -----------------
# Loading libaries & functions
# -----------------
library(mgcv)
source("R/0_Calculate_GAMM_Posteriors.R")
# -----------------

# Some variables to get set
n.out <- n <- 250
sites.dat <- unique(ecosys$Site)
ns        <- length(sites.dat)

for(m in unique(ci.terms$Model)){
	# Select which variables to extract
	vars <- unique(ci.terms[ci.terms$Model==m, "Effect"])

	# re-running everything on the base level (t.001, 850-2010) because that's what everything will get compared to
	# -----------------
	# Create the baseline data frame
	# -----------------
	for(i in 1:ns){
		if(i == 1){ 
			site.vec <- paste(rep(sites.dat[i], n.out) )
		} else {
			site.vec <- c(site.vec, paste(rep(sites.dat[i], n.out)))
		}
	}

	new.dat <- data.frame(	Site=site.vec,
							Extent=as.factor("850-2010"), 
							Scale="t.001")
	for(v in vars){
		new.dat[,v] <- rep(seq(min(ecosys[ecosys$Scale=="t.001",v],   na.rm=T), max(ecosys[ecosys$Scale=="t.001",v],   na.rm=T), length.out=n.out), ns)
	}								
	# -----------------

	# -----------------
	# finding the right Rdata files to open
	# -----------------
	model.order <- unique(ecosys[ecosys$Model==m, "Model.Order"])
	fmod <- dir(file.path(in.base, model.order, in.res))
	load(file.path(in.base, model.order, in.res, fmod))
	# -----------------

	# -----------------
	# getting the baseline post.distns
	# -----------------
	mod1 <- mod.out[["gamm.850-2010.001"]]
	ci.terms.pred1 <- post.distns(model.gam=mod1, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
	
	# relativizing the effect to the mean npp
	npp <- mean(ecosys[ecosys$Model==m & ecosys$Scale=="t.001" & ecosys$Year>=850, "NPP"], na.rm=T)			
	ci.terms.pred1$ci[,c("mean", "lwr", "upr")] <- ci.terms.pred1$ci[,c("mean", "lwr", "upr")]/npp
	# -----------------

	for(r in unique(ci.terms[ci.terms$Model==m, "Scale"])){
		mod2 <- mod.out[[paste0("gamm.850-2010.", substr(r, 3, 5))]]
		# just changing the scale (resolution) in new.dat to the appropriate r
		new.dat$Scale <- as.factor(r)

		ci.terms.pred2 <- post.distns(model.gam=mod2, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
		ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)] <- (ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)]/npp)-ci.terms.pred1$ci$mean

		df.out <- data.frame(Model=as.factor(m), Site=ci.terms.pred2$sims$Site, Extent=ci.terms.pred2$sims$Extent, Scale=ci.terms.pred2$sims$Scale, Effect=as.factor(ci.terms.pred2$sims$Effect), x=ci.terms.pred2$sims$x,  mean.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, mean, na.rm=T), lwr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.025, na.rm=T), upr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.975, na.rm=T))

		if(m==unique(ci.terms$Model)[1] & r==unique(ci.terms[ci.terms$Model==m, "Scale"])[1]) {
			ci.terms.rel <- df.out
		} else {
			ci.terms.rel <- rbind(ci.terms.rel, df.out)
		}
	}
	
	# changing the resolution back to "t.001"
	new.dat$Scale <- as.factor("t.001")

	ext <- unique(ci.terms[ci.terms$Model==m, "Extent"])
	for(e in ext[2:length(ext)]){  # note, not doing the first one because it's essentially aready done
		# Need to switch and load the extent file
		fmod <- dir(file.path(in.base, model.order, in.ext))
		load(file.path(in.base, model.order, in.ext, fmod))

		r <- ifelse(e=="0850-2010", "001", NA)
		mod2 <- mod.out[[paste0("gamm.", e, ".", r)]]
		# just changing the scale (resolution) in new.dat to the appropriate r
		new.dat$Extent <- as.factor(e)

		ci.terms.pred2 <- post.distns(model.gam=mod2, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
		ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)] <- (ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)]/npp)-ci.terms.pred1$ci$mean

		df.out <- data.frame(Model=as.factor(m), Site=ci.terms.pred2$sims$Site, Extent=ci.terms.pred2$sims$Extent, Scale=ci.terms.pred2$sims$Scale, Effect=as.factor(ci.terms.pred2$sims$Effect), x=ci.terms.pred2$sims$x, mean.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, mean, na.rm=T), lwr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.025, na.rm=T), upr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.975, na.rm=T))

		ci.terms.rel <- rbind(ci.terms.rel, df.out)
	}
}
# Fixing Extent Labels
ci.terms.rel$Extent <- as.factor(ifelse(ci.terms.rel$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
summary(ci.terms.rel)
save(ci.terms.rel, file=file.path(out.dir, "Driver_Responses_DevAnn.Rdata"))

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Resolution.pdf"))
ggplot(data= ci.terms.rel[ci.terms.rel$Extent=="0850-2010",]) + 
	facet_grid(Scale~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Extent.pdf"))
ggplot(data= ci.terms.rel[ci.terms.rel$Scale=="t.001",]) + 
	facet_grid(Extent~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_ResAnn_Ext1990.pdf"))
ggplot(data=ci.terms.rel[ci.terms.rel$Scale=="t.001" & ci.terms.rel$Extent=="1990-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Change in Driver Sensitivity, Resolution: Annual, Extent: 1990-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Res100_Ext0850.pdf"))
ggplot(data=ci.terms.rel[ci.terms.rel$Scale=="t.100" & ci.terms.rel$Extent=="0850-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Change in Driver Sensitivity, Resolution: Centennial, Extent: 0850-2010") +
	theme_bw()
dev.off()
# ----------------------------------------

# ----------------------------------------
# Summary graphs & statistics
# Selecting 3 sets of curves
# 1) Base Effect (BE) = the percent change in NPP per unit change in the driver (res=001, ext=850-2010)
# 2) Res. Effect (RE) = the effect of change in resolution; deviation from the BE at res=250 (ext=850-2010)
# 3) Ext. Effect (EE) = the effect of change in extent; deviation from the BE at ext=1990-2010 (res=001)
# ----------------------------------------
# load(file.path(out.dir, "Driver_Responses_DevAnn.Rdata"))

# Try to figure out how to make it a 3-panel figure
data.final <- ci.terms[ci.terms$Scale=="t.001" & ci.terms$Extent=="0850-2010", c("Model", "Site", "Extent", "Scale", "Effect", "x", "mean.rel", "lwr.rel", "upr.rel")]
data.final$std <- as.factor("model.mean")
summary(data.final)

ci.terms.rel$std <- as.factor("t.001.0850-2010")
summary(ci.terms.rel)

data.final <- rbind(data.final, ci.terms.rel[ci.terms.rel$Scale=="t.100" | ci.terms.rel$Extent=="1990-2010",])
data.final$Panel <- as.factor(paste(data.final$std, data.final$Scale, data.final$Extent, sep="-"))
data.final$Panel <- recode(data.final$Panel, "'model.mean-t.001-0850-2010'='1'; 't.001.0850-2010-t.100-850-2010'='2'; 't.001.0850-2010-t.001-1990-2010'='3'")
levels(data.final$Panel) <- c("Base Effect", "Resolution Deviation", "Extent Deviation")
summary(data.final)

# Making a figure that illustrates each of the Effects
ecosys2 <- ecosys[!ecosys$Model=="clm.bgc" & !ecosys$Model=="linkages",]
col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys2$Model.Order),"color"])
yrs.100  <- seq(from=2009, to=min(ecosys$Year), by=-100)


pdf(file.path(fig.dir, "NPP_PHA_Scales.pdf"))
ggplot(data=ecosys2[ecosys2$Site=="PHA",]) +
	geom_line(data=ecosys2[ecosys2$Scale=="t.001" & ecosys2$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys2[ecosys2$Scale=="t.100" & (ecosys2$Year %in% yrs.100) & ecosys2$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	geom_vline(xintercept=1990, linetype="dashed") +
  scale_x_continuous(limits=c(850,2010)) +
  scale_color_manual(values=col.model) +
 	labs(y="NPP MgC ha-1 yr-1") +
 	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_SummaryFigure_byVar.pdf"))
for(v in unique(data.final$Effect)){
# Getting the right color set for each variable
models.use <- unique(ecosys[ecosys$Model %in% data.final[data.final$Effect==v,"Model"],"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])

print(
ggplot(data=data.final[data.final$Effect==v,]) + 
	facet_wrap( ~ Panel, scales="free_x", nrow=1) +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_y_continuous(limits=range(data.final[data.final$Effect==v,c("mean.rel", "lwr.rel", "upr.rel")], na.rm=T)) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title=paste0("Change in Driver Sensitivity: ", v)) +
	theme_bw()
)	
}
dev.off()
# ----------------------------------------
