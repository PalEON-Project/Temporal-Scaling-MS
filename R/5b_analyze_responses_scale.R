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
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

# Figure out what models we have to work with
#response.all <- c("NPP", "AGB.diff", "NEE")
response.all <- c("NPP")
models <- unique(ecosys$Model)
f.res <- dir(file.path(in.base, in.res))
f.ext <- dir(file.path(in.base, in.ext))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")

# Put all responses & models & scales into single data frames
# for(r in 1:length(response.all)){
  #response <- response.all[r]
  response="NPP"
  response.res <- grep(response, f.res)
  response.ext <- grep(response, f.ext)
  for(i in 1:length(models)){
    # First narrow to the models

	fmod <- response.res[grep(models2[i], f.res[response.res])]

    if(!length(fmod)>0) next

    load(file.path(in.base, in.res, f.res[fmod]))
    
    if(i==1) {
      ci.terms   <- mod.out$ci.terms
      dat.ecosys <- cbind(mod.out$data[,], mod.out$ci.response[,c("mean", "lwr", "upr")])
	  wt.terms   <- mod.out$weights
      # sim.terms <- mod.out$sim.terms
    } else {
      ci.terms   <- rbind(ci.terms, mod.out$ci.terms )
      dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]))
	  wt.terms   <- merge(wt.terms, mod.out$weights, all.x=T, all.y=T) # Note: need to use merge bc this was only done with the relevant drivers
      # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    }
    
    # loop through by extent
    fmod <- response.ext[grep(models2[i], f.ext[response.ext])]
    load(file.path(in.base, in.ext, f.ext[fmod]))
    
    # Note: because we're lumping everything together, let's not mess with reiterating the base level
    ci.terms  <- rbind(ci.terms, mod.out$ci.terms[!(mod.out$ci.terms$Resolution=="t.001" & substr(mod.out$ci.terms$Extent,1,3)=="850"),] )
    dat.ecosys <- rbind(dat.ecosys,
                        cbind(mod.out$data[!(mod.out$data$Resolution=="t.001" & substr(mod.out$data$Extent,1,3)=="850"),], 
                              mod.out$ci.response[!(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850"),c("mean", "lwr", "upr")]))
	  wt.terms   <- merge(wt.terms, mod.out$weights[!(mod.out$weights$Resolution=="t.001" & substr(mod.out$weights$Extent,1,3)=="850"),], all.x=T, all.y=T)
    # # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    # !(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850")
    # }
    
    # Clear the mod.out to save space
    rm(mod.out)
  }
  #ci.terms$Response <- response.all[r]
  
#}# Fix Extent Labels for consistency
ci.terms$Extent   <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
wt.terms$Extent   <- as.factor(ifelse(wt.terms$Extent=="850-2010", "0850-2010", paste(wt.terms$Extent)))
summary(ci.terms)
summary(dat.ecosys)
summary(wt.terms)

# Get rid of Linkages, because it's just weird
ci.terms   <- ci.terms[!ci.terms$Model=="linkages",]
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="linkages",]
wt.terms   <- wt.terms[!wt.terms$Model=="linkages",]
ecosys     <- ecosys[!ecosys$Model=="linkages",]


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
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])



# Graphing Individual Models with the Base Effect Predictions
pdf(file.path(fig.dir, "NPP_Annual_AllSites.pdf"))
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.1) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	scale_color_manual(values=colors.use) +
	theme_bw()
for(s in unique(dat.ecosys$Site)){
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site==s,])  +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.25) +
	# geom_line(aes(x=Year, y=fit.gam, color=Model.Order), alpha=0.8, size=0.25) +
	scale_color_manual(values=colors.use) +
	ggtitle(paste0("NPP - ", s)) +
	theme_bw()
)
}
dev.off()


extent.box <- data.frame(Year=850:1850, Min=0, Max=20)

pdf(file.path(fig.dir, "NPP_Scales_PHA_0850-2010.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",])  +
	# geom_ribbon(data=extent.box, aes(x=Year, ymin=Min, ymax=Max), alpha=0.2) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.5, alpha=0.35) +
	geom_line(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_point(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=5) +
	geom_vline(xintercept=1850, linetype="dashed", size=1.5) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(limits=c(0,20), expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
	labs(color="Model", x="Year (A.D.)", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1.75)), 
	      legend.title=element_text(size=rel(2)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1))
)
dev.off()


pdf(file.path(fig.dir, "NPP_ModelFits_PHA.pdf"))
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.5, color="gray40") +
	scale_color_manual(values=colors.use) +
	# scale_x_continuous(limits=c(1900,2010)) +
	ggtitle("PHA NPP 0850-2010") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	scale_color_manual(values=colors.use) +
	scale_x_continuous(limits=c(1900,2010)) +
	ggtitle("PHA NPP 1900-2010") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	scale_color_manual(values=colors.use) +
	scale_x_continuous(limits=c(1800,1820)) +
	ggtitle("PHA NPP 1800-1820") +
	theme_bw() + guides(color=F)
)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	scale_color_manual(values=colors.use) +
	scale_x_continuous(limits=c(1990,2010)) +
	ggtitle("PHA NPP 1990-2010") +
	theme_bw() + guides(color=F)
)
dev.off()
# ----------------------------------------






# ----------------------------------------
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Across all scales (res + extent) finding the mean NPP
# NOTE: not worrying about Site because the terms are across sites
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){
		for(e in unique(ci.terms[ci.terms$Model==m & ci.terms$Resolution==r, "Extent"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# -----------------------
			# Find the start year for the extent
			yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=yr & dat.ecosys$Resolution==r, "NPP"], na.rm=T)			# -----------------------
			
			# -----------------------
			# Relativizing everythign in dat.ecosys to make it comparable to tree rings
			# -----------------------
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"NPP"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"lwr"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Extent==e,"upr"])/npp
			# -----------------------

			
			# -----------------------
			# Finding the percent change in NPP relative to the mean for that particular scale
			# -----------------------
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Extent==e,"upr"])/npp
			# -----------------------
		}
	}
}
summary(dat.ecosys)
summary(ci.terms)

# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms$Model,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])

# Creating a cheat data frame that lets values go off the graph
ci.terms.graph <- ci.terms
ci.terms.graph[ci.terms.graph$mean.rel<(-1.25),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel<(-1.25),"lwr.rel"] <- -1.25 
ci.terms.graph[ci.terms.graph$upr.rel<(-1.25),"upr.rel"] <- -1.25 
ci.terms.graph[which(ci.terms.graph$mean.rel>1.25),"mean.rel"] <- NA 
ci.terms.graph[ci.terms.graph$lwr.rel>(1.25),"lwr.rel"] <- 1.25 
ci.terms.graph[ci.terms.graph$upr.rel>(1.25),"upr.rel"] <- 1.25 
summary(ci.terms.graph)

pdf(file.path(fig.dir, "NPP_RelChange_Resolution.pdf"))
ggplot(data= ci.terms.graph[ci.terms.graph$Extent=="0850-2010",]) + 
	facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_Extent.pdf"))
ggplot(data= ci.terms[ci.terms$Resolution=="t.001",]) + 
	facet_grid(Extent~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_BaseEffect.pdf"))
ggplot(data=ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="0850-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity to May-Sep Met, 850-2010") +
	guides(col=guide_legend(ncol=2), fill=guide_legend(ncol=2)) +
	theme(legend.position=c(0.85, 0.15)) +
	theme(plot.title=element_text(face="bold", size=rel(1))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()) 
dev.off()


big3 <- c("tair", "precipf", "CO2")
pdf(file.path(fig.dir, "NPP_RelChange_BaseEffect_Big3.pdf"))
ggplot(data=ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="0850-2010" & ci.terms$Effect %in% big3,]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

# A very crude way of showing the effects at 3 different scales
E1 <- ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="0850-2010",]
E1$Scale <- as.factor(1)
levels(E1$Scale) <- "Base Effect"
summary(E1)

E2 <- ci.terms[ci.terms$Resolution=="t.050" & ci.terms$Extent=="0850-2010",]
E2$Scale <- as.factor(2)
levels(E2$Scale) <- "Resolution: 50-yr"
summary(E2)

E3 <- ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="1850-2010",]
E3$Scale <- as.factor(3)
levels(E3$Scale) <- "Extent: 1850-2010"
summary(E3)

ci.terms2 <- rbind(E1, E2, E3)
summary(ci.terms2)

pdf(file.path(fig.dir, "NPP_RelChange_byDriver.pdf"))
for(E in unique(ci.terms$Effect)){
print(
ggplot(data=ci.terms2[ci.terms2$Effect==E,]) + 
	facet_grid(.~Scale) +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title=E) +
	theme_bw()
)
}
dev.off()

big3 <- c("tair", "precipf", "CO2")

pdf(file.path(fig.dir, "NPP_RelChange_Big3.pdf"))
print(
ggplot(data=ci.terms2[ci.terms2$Effect %in% big3,]) + 
	facet_grid(Scale~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP") +
	guides(col=guide_legend(nrow=2), fill=guide_legend(nrow=2)) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(
		  # legend.text=element_text(size=rel(1.5)), 
	      # legend.title=element_text(size=rel(1.5)),
	      legend.key=element_blank() #,
	      # legend.key.size=unit(2, "lines"),  
	      # legend.key.width=unit(2, "lines") 
	      ) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank() #, 
	      # axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      # axis.text.y=element_text(color="black", size=rel(2.5)), 
	      # axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      # axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1)
	      )
)
dev.off()
# ----------------------------------------


# ----------------------------------------
# Using the model-wide weights to do the drivers through time by site
# ----------------------------------------
source("R/0_GAMM_Plots.R")

models.all <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
models.use <- models.all[as.numeric(levels(wt.terms$Model.Order))]

wt.terms$Model.Order <- recode(wt.terms$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
levels(wt.terms$Model.Order) <- models.all[as.numeric(levels(wt.terms$Model.Order))]
summary(wt.terms)

indices.wt <- wt.terms$Site=="PHA" & wt.terms$Resolution=="t.001" & wt.terms$Extent=="0850-2010" 
indices.dat <- dat.ecosys$Site=="PHA" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Extent=="0850-2010" 

pdf(file.path(fig.dir, "NPP_Drivers_Time_PHA.pdf"))
print(
ggplot(data=wt.terms[indices.wt,]) + facet_wrap(~Model) +
 	geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="clm.cn",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2.lu",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.stat",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.triffid",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.guess",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.precipf"])), size=2) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.wsl",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.precipf"])), size=2)+
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="sibcasa",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.precipf"])), size=2) +
	scale_x_continuous(limits=c(850,2010), expand=c(0,0)) +
    theme_bw()
)
print(
ggplot(data=wt.terms[indices.wt,]) + facet_wrap(~Model) +
 	geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="clm.cn",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2.lu",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.stat",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.triffid",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.guess",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.wsl",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.precipf"])), size=3)+
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="sibcasa",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.precipf"])), size=3) +
	scale_x_continuous(limits=c(1800,1900), expand=c(0,0)) +
    theme_bw()
)
print(
ggplot(data=wt.terms[indices.wt,]) + facet_wrap(~Model) +
 	geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="clm.cn",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2.lu",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.stat",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.triffid",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.guess",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.wsl",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.precipf"])), size=3)+
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="sibcasa",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.precipf"])), size=3) +
	scale_x_continuous(limits=c(1901,2010), expand=c(0,0)) +
    theme_bw()
)
dev.off()

pdf(file.path(fig.dir, "NPP_Drivers_Time_PHA_1800-2010.pdf"), width=11, height=8.5)
print(
ggplot(data=wt.terms[indices.wt,]) + facet_wrap(~Model.Order) +
 	geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="clm.cn",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2.lu",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.stat",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.triffid",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.guess",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.precipf"])), size=3) +
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.wsl",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.precipf"])), size=3)+
	geom_line(data=wt.terms[indices.wt & wt.terms$Model=="sibcasa",], aes(x=Year, y=fit.full),
	          color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.precipf"])), size=3) +
	scale_x_continuous(limits=c(1700,2010), expand=c(0,0), breaks=c(1750, 1850, 1950)) +
	scale_y_continuous(name=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")"))), expand=c(0,0)) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(2.5)), 
	      axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1)) +
	theme(strip.text=element_text(size=rel(1.5)))
)
dev.off()
# ----------------------------------------








# ----------------------------------------
# Look at scale-based changes in driver responses
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
sites.dat <- unique(dat.ecosys$Site)
ns        <- length(sites.dat)

# Figure out what models we have to work with
models <- unique(dat.ecosys$Model)
f.res <- dir(file.path(in.base, in.res))
f.ext <- dir(file.path(in.base, in.ext))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")


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
							Resolution="t.001")
	for(v in vars){
		new.dat[,v] <- rep(seq(min(dat.ecosys[dat.ecosys$Resolution=="t.001",v], na.rm=T), max(dat.ecosys[dat.ecosys$Resolution=="t.001",v], na.rm=T), length.out=n.out), ns)
	}								
	# -----------------

	# -----------------
	# finding the right Rdata files to open
	# -----------------
	model.order <- unique(dat.ecosys[dat.ecosys $Model==m, "Model.Order"])
	
	m.name2 <- ifelse(m=="ed2", "ed2_", m)
	fmod <- grep(m.name2, f.res)
	load(file.path(in.base, in.res, f.res[fmod]))
	# -----------------

	# -----------------
	# getting the baseline post.distns
	# -----------------
	mod1 <- mod.out[["gamm.850-2010.001"]]
	ci.terms.pred1 <- post.distns(model.gam=mod1, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
	
	# relativizing the effect to the mean npp
	npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution=="t.001" & dat.ecosys$Year>=850, "NPP"], na.rm=T)			
	ci.terms.pred1$ci[,c("mean", "lwr", "upr")] <- ci.terms.pred1$ci[,c("mean", "lwr", "upr")]/npp
	# -----------------

	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){
		mod2 <- mod.out[[paste0("gamm.850-2010.", substr(r, 3, 5))]]
		# just changing the scale (resolution) in new.dat to the appropriate r
		new.dat$Resolution <- as.factor(r)

		ci.terms.pred2 <- post.distns(model.gam=mod2, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
		ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)] <- (ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)]/npp)-ci.terms.pred1$ci$mean

		df.out <- data.frame(Model=as.factor(m), Site=ci.terms.pred2$sims$Site, Extent=ci.terms.pred2$sims$Extent, Resolution=ci.terms.pred2$sims$Resolution, Effect=as.factor(ci.terms.pred2$sims$Effect), x=ci.terms.pred2$sims$x,  mean.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, mean, na.rm=T), lwr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.025, na.rm=T), upr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.975, na.rm=T))

		if(m==unique(ci.terms$Model)[1] & r==unique(ci.terms[ci.terms$Model==m, "Resolution"])[1]) {
			ci.terms.rel <- df.out
		} else {
			ci.terms.rel <- rbind(ci.terms.rel, df.out)
		}
	}
	
	# changing the resolution back to "t.001"
	new.dat$Resolution <- as.factor("t.001")

	ext <- unique(ci.terms[ci.terms$Model==m, "Extent"])
	# Need to switch and load the extent file
	fmod <- grep(m.name2, f.ext)
	load(file.path(in.base, in.ext, f.ext[fmod]))
	for(e in ext[2:length(ext)]){  # note, not doing the first one because it's essentially aready done

		# r <- ifelse(e=="850-2010", "001", NA)
		mod2 <- mod.out[[paste0("gamm.", e)]]
		# just changing the scale (resolution) in new.dat to the appropriate r
		new.dat$Extent <- as.factor(e)

		ci.terms.pred2 <- post.distns(model.gam=mod2, model.name=m, n=n, newdata=new.dat, vars=vars, terms=T)
		ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)] <- (ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)]/npp)-ci.terms.pred1$ci$mean

		df.out <- data.frame(Model=as.factor(m), Site=ci.terms.pred2$sims$Site, Extent=ci.terms.pred2$sims$Extent, Resolution=ci.terms.pred2$sims$Resolution, Effect=as.factor(ci.terms.pred2$sims$Effect), x=ci.terms.pred2$sims$x, mean.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, mean, na.rm=T), lwr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.025, na.rm=T), upr.rel=apply(ci.terms.pred2$sims[,7:ncol(ci.terms.pred2$sims)], 1, quantile, 0.975, na.rm=T))

		ci.terms.rel <- rbind(ci.terms.rel, df.out)
	}
}

# Fixing Extent Labels
ci.terms.rel$Extent <- as.factor(ifelse(ci.terms.rel$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
summary(ci.terms.rel)
save(ci.terms.rel, file=file.path(out.dir, "Driver_Responses_DevAnn.Rdata"))

# ------------------------------------
# Trimming out the x values outside the range seen in the actual data for graphing etc
# ------------------------------------
for(f in unique(ci.terms.rel$Effect)){
	for(e in unique(ci.terms.rel$Extent)){
		x.min <- min(ci.terms[ci.terms$Effect==f & ci.terms$Extent==e & ci.terms$Resolution=="t.001", "x" ]) 
		x.max <- max(ci.terms[ci.terms$Effect==f & ci.terms$Extent==e & ci.terms$Resolution=="t.001", "x" ]) 
		if(f==unique(ci.terms.rel$Effect)[1] & e==unique(ci.terms.rel$Extent)[1]){
			ci.terms.rel2 <- ci.terms.rel[ci.terms.rel$Effect==f & 
			                              ci.terms.rel$Extent==e &
			                              ci.terms.rel$Resolution=="t.001" & 
			                              ci.terms.rel$x >= x.min & 
			                              ci.terms.rel$x <= x.max,
			                              ]
		} else {
			ci.terms.rel2 <-rbind(ci.terms.rel2, ci.terms.rel[ci.terms.rel$Effect==f & 
			                              ci.terms.rel$Extent==e &
			                              ci.terms.rel$Resolution=="t.001" & 
			                              ci.terms.rel$x >= x.min & 
			                              ci.terms.rel$x <= x.max,
			                              ])
		}
		}
	for(r in unique(ci.terms.rel$Resolution)){
		x.min <- min(ci.terms[ci.terms$Effect==f & ci.terms$Extent=="0850-2010" & ci.terms$Resolution==r, "x" ]) 
		x.max <- max(ci.terms[ci.terms$Effect==f & ci.terms$Extent=="0850-2010" & ci.terms$Resolution==r, "x" ]) 
		ci.terms.rel2 <-rbind(ci.terms.rel2, ci.terms.rel[ci.terms.rel$Effect==f & 
			                              ci.terms.rel$Extent=="0850-2010" &
			                              ci.terms.rel$Resolution==r & 
			                              ci.terms.rel$x >= x.min & 
			                              ci.terms.rel$x <= x.max,
			                              ])
	}
}
# Also just going ahead and changing the limits to make pretty graphs as well
ci.terms.rel2[which(ci.terms.rel2$mean.rel<(-0.5)),"mean.rel"] <- NA 
ci.terms.rel2[ci.terms.rel2$lwr.rel<(-0.5),"lwr.rel"] <- -0.5
ci.terms.rel2[ci.terms.rel2$upr.rel<(-0.5),"upr.rel"] <- -0.5
ci.terms.rel2[which(ci.terms.rel2$mean.rel>0.5),"mean.rel"] <- NA 
ci.terms.rel2[ci.terms.rel2$lwr.rel>(0.5),"lwr.rel"] <- 0.5
ci.terms.rel2[ci.terms.rel2$upr.rel>(0.5),"upr.rel"] <- 0.5 
summary(ci.terms.rel2)
# ------------------------------------
# Trying out the basic plot to compare model responses to drivers
models.use <- unique(dat.ecosys[dat.ecosys$Model %in% ci.terms.rel2$Model,"Model.Order"])
models.use <- models.use[order(models.use)]
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])[order(models.use)]

models.df <- data.frame(Model=unique(dat.ecosys$Model), Model.Order=unique(dat.ecosys$Model.Order))
ci.terms.rel2 <- merge(ci.terms.rel2, models.df, all.x=T, all.y=F)

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Resolution.pdf"))
ggplot(data= ci.terms.rel2[ci.terms.rel2$Extent=="0850-2010",]) + 
	facet_grid(Resolution~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model.Order), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Extent.pdf"))
ggplot(data= ci.terms.rel2[ci.terms.rel2$Resolution=="t.001",]) + 
	facet_grid(Extent~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model.Order), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity, Resolution: Annual, Extent: 850-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Res001_Ext1990.pdf"))
ggplot(data=ci.terms.rel2[ci.terms.rel2$Resolution=="t.001" & ci.terms.rel2$Extent=="1990-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model.Order), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Change in Driver Sensitivity, Resolution: Annual, Extent: 1990-2010") +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_fromAnn_Res100_Ext0850.pdf"))
ggplot(data=ci.terms.rel2[ci.terms.rel2$Resolution=="t.100" & ci.terms.rel2$Extent=="0850-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model.Order), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Model.Order), size=0.75) +
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
data.final <- ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="0850-2010", c("Model", "Site", "Extent", "Resolution", "Effect", "x", "mean.rel", "lwr.rel", "upr.rel")]
data.final$std <- as.factor("model.mean")
summary(data.final)

ci.terms.rel$std <- as.factor("t.001.0850-2010")
summary(ci.terms.rel)

data.final <- rbind(data.final, ci.terms.rel[ci.terms.rel$Resolution=="t.100" | ci.terms.rel$Extent=="1990-2010",])
data.final$Panel <- as.factor(paste(data.final$std, data.final$Resolution, data.final$Extent, sep="-"))
data.final$Panel <- recode(data.final$Panel, "'model.mean-t.001-0850-2010'='1'; 't.001.0850-2010-t.100-850-2010'='2'; 't.001.0850-2010-t.001-1990-2010'='3'")
levels(data.final$Panel) <- c("Base Effect", "Resolution Deviation", "Extent Deviation")
summary(data.final)

# Making a figure that illustrates each of the Effects
ecosys2 <- ecosys[!ecosys$Model=="clm.bgc" & !ecosys$Model=="linkages",]
col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys2$Model.Order),"color"])
yrs.100  <- seq(from=2009, to=min(ecosys$Year), by=-100)


pdf(file.path(fig.dir, "NPP_PHA_Resolutions.pdf"))
ggplot(data=ecosys2[ecosys2$Site=="PHA",]) +
	geom_line(data=ecosys2[ecosys2$Resolution=="t.001" & ecosys2$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys2[ecosys2$Resolution=="t.100" & (ecosys2$Year %in% yrs.100) & ecosys2$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=1.5) +
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
