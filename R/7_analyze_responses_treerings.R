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
in.tr  <- "AllDrivers_GS_TreeRings/TreeRings_All"
out.dir <- "Data/analysis_response_scale_TreeRings"
fig.dir <- "Figures/analysis_response_scale_TreeRings"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData_Raw.Rdata"))
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

# Figure out what models we have to work with
models <- c(paste(unique(ecosys$Model)))
f.tr <- dir(file.path(in.base, in.tr))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")

# Note Ignoring tree rings for the moment becaus they're a bit different
for(i in 1:length(models)){
	# loop through by resolution
	fmod <- grep(models2[i], f.tr)
	load(file.path(in.base, in.tr, f.tr[fmod]))
	
	if(i==1) {
		ci.terms   <- mod.out$ci.terms
		dat.ecosys <- cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")])
		# sim.terms <- mod.out$sim.terms
	} else {
		ci.terms  <- rbind(ci.terms, mod.out$ci.terms )
		dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]))
		# sim.terms <- rbind(sim.terms, mod.out$sim.terms)
	}

	# Clear the mod.out to save space
	rm(mod.out)
}
summary(ci.terms)
summary(dat.ecosys)

# Add in the tree rings
fmod <- grep("tree.rings", f.tr)
load(file.path(in.base, in.tr, f.tr[fmod]))
summary(mod.out$data)

ci.terms <- rbind(ci.terms, mod.out$ci.terms)
dat.ecosys <- merge(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]), all.x=T, all.y=T)
summary(dat.ecosys)


# Get rid of Linkages, because it's just weird
ci.terms   <- ci.terms[!ci.terms$Model=="linkages",]
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="linkages",]

# Write the files to csv so I don't have to mess with loading large .Rdata files again if I don't have to
# write.csv(ci.terms.res,  file.path(out.dir, "Driver_Responses_CI_Resolution.csv"), row.names=F)
# write.csv(sim.terms.res, file.path(out.dir, "Driver_Responses_Sims_Resolution.csv"), row.names=F)
# write.csv(ci.terms,  file.path(out.dir, "Driver_Responses_CI_Extent.csv"), row.names=F)
# write.csv(sim.terms, file.path(out.dir, "Driver_Responses_Sims_Extent.csv"), row.names=F)
# ----------------------------------------

# ----------------------------------------
# Some Plotting etc of NPP
# ----------------------------------------
summary(dat.ecosys)

models.use <- unique(dat.ecosys[,"Model.Order"])
colors.use <- c(as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"]))

# relativizing NPP
for(m in unique(dat.ecosys$Model)){
	for(r in unique(dat.ecosys[dat.ecosys$Model==m, "Resolution"])){
		npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r, "NPP"], na.rm=T)

		# Relativizing everythign in dat.ecosys to make it comparable to tree rings
		dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"NPP"])/npp
		dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
		dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"mean"])/npp
		dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"lwr"])/npp
		dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r,"upr"])/npp

		# Finding the percent change in NPP relative to the mean for that particular scale
		ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"mean"])/npp
		ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"lwr"])/npp
		ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r,"upr"])/npp

	}
}
summary(dat.ecosys)
dat.ecosys$Model.Order2 <- as.factor(ifelse(dat.ecosys$Model=="tree.rings", paste0("Tree Rings, ", dat.ecosys$PlotID), paste(dat.ecosys$Model.Order)))
levels(dat.ecosys$Model.Order2)

colors.use2 <- c(colors.use, rep("gray50", length(unique(dat.ecosys[dat.ecosys$Model=="tree.rings", "Model.Order2"]))))

# Graphing Individual Models with the Base Effect Predictions
pdf(file.path(fig.dir, "NPP_Annual_AllSites.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order2), size=0.25) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	# geom_line(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Model=="tree.rings",], 
			  # aes(x=Year, y=NPP), color="black", size=0.25) +
	labs(x="Year", y="NPP (Mg C ha-1 yr-1)", title="NPP") +
	scale_color_manual(values=colors.use2) +
	theme_bw()
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP.rel*100, color=Model.Order2), size=0.25) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	# geom_line(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Model=="tree.rings",], 
			  # aes(x=Year, y=NPP.rel*100, color=Model), color="black", size=0.25) +
	labs(x="Year", y="%dNPP", title="Relative NPP, full tree ring time") +
	scale_y_continuous(limits=c(0,200)) +
	scale_color_manual(values=colors.use2) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      # panel.border=element_blank(), 
	      panel.background=element_blank())
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP.rel*100, color=Model.Order2), size=0.8) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	# geom_line(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Model=="tree.rings",], 
			  # aes(x=Year, y=NPP.rel*100, color=PlotID), size=0.5) +
	labs(x="Year", y="%dNPP", title="Relative NPP, 1990-2010") +
	# scale_y_continuous(limits=c(0,200)) +
	scale_x_continuous(limits=c(1950,1970)) +
	scale_color_manual(values=colors.use2) +
	# guides(col=guide_legend(nrow=2)) +
	# theme(legend.position="top") +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      # panel.border=element_blank(), 
	      panel.background=element_blank())
)

for(s in unique(dat.ecosys$Site)){
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Site==s,])  +
	geom_line(aes(x=Year, y=NPP.rel, color=Model.Order2), size=0.5) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	scale_y_continuous(limits=c(0,2)) +
	scale_x_continuous(limits=c(1900,2010)) +
	scale_color_manual(values=c(colors.use, rep("black", 5))) +
	ggtitle(paste0("NPP - ", s)) +
	theme_bw()
)
}
dev.off()


extent.box <- data.frame(Year=850:1850, Min=0, Max=20)


pdf(file.path(fig.dir, "NPP.Rel_ModelFits_PHA.pdf"))
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order2) +
	geom_line(aes(x=Year, y=NPP.rel, color=Model.Order2), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr.rel, ymax=upr.rel), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam.rel), alpha=0.8, size=0.5, color="gray40") +
	scale_color_manual(values=c(colors.use, rep("black", 5))) +
	# scale_x_continuous(limits=c(1900,2010)) +
	ggtitle("NPP 0850-2010") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order2) +
	geom_line(aes(x=Year, y=NPP.rel, color=Model.Order2), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr.rel, ymax=upr.rel), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam.rel), alpha=0.8, size=0.5, color="gray40") +
	scale_color_manual(values=c(colors.use, rep("black", 5))) +
	scale_x_continuous(limits=c(1900,2010)) +
	ggtitle("NPP 1900-2010") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order2) +
	geom_line(aes(x=Year, y=NPP.rel, color=Model.Order2), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr.rel, ymax=upr.rel), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam.rel), alpha=0.8, size=0.5, color="gray40") +
	scale_color_manual(values=c(colors.use, rep("black", 5))) +
	scale_x_continuous(limits=c(1950,2010)) +
	ggtitle("NPP 1950-2010") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order2) +
	geom_line(aes(x=Year, y=NPP.rel, color=Model.Order2), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr.rel, ymax=upr.rel), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam.rel), alpha=0.8, size=0.5, color="gray40") +
	scale_color_manual(values=c(colors.use, rep("black", 5))) +
	scale_x_continuous(limits=c(1990,2010)) +
	ggtitle("NPP 1990-2010") +
	theme_bw() + guides(color=F)
)
dev.off()
# ----------------------------------------


# ----------------------------------------
# Standardize driver responses to the mean model NPP
# ----------------------------------------

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- c(as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black")

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-0.5),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-0.5),"lwr.rel"] <- -0.5 
ci.terms.graph[ci.terms.graph$upr.rel<(-0.5),"upr.rel"] <- -0.5 
ci.terms.graph[which(ci.terms.graph$mean.rel>0.5),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(0.5),"lwr.rel"] <- 0.5 
ci.terms.graph[ci.terms.graph$upr.rel>(0.5),"upr.rel"] <- 0.5 
summary(ci.terms.graph)

pdf(file.path(fig.dir, "NPP_RelChange_Resolution.pdf"))
ggplot(data= ci.terms.graph[,]) + 
	facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_BaseEffect.pdf"))
ggplot(data=ci.terms.graph[ci.terms.graph$Resolution=="t.001",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 1901-2010") +
	theme_bw()
dev.off()


big3 <- c("tair", "precipf", "CO2")
pdf(file.path(fig.dir, "NPP_RelChange_BaseEffect_Big3.pdf"))
ggplot(data=ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Effect %in% big3,]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 1901-2010") +
	theme_bw()
dev.off()


# ----------------------------------------

