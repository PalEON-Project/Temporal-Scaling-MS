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
in.res  <- "Big4_GS_byResolution"
in.ext  <- "Big4_GS_byExtent"
out.dir <- "Data/analysis_response_scale_big4"
fig.dir <- "Figures/analysis_response_scale_big4"

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
      # ci.terms   <- mod.out$ci.terms
      dat.ecosys <- cbind(mod.out$data[,], mod.out$ci.response[,c("mean", "lwr", "upr")])
	  # wt.terms   <- mod.out$weights
      # sim.terms <- mod.out$sim.terms
    } else {
      # ci.terms   <- rbind(ci.terms, mod.out$ci.terms )
      dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]))
	  # wt.terms   <- merge(wt.terms, mod.out$weights, all.x=T, all.y=T) # Note: need to use merge bc this was only done with the relevant drivers
      # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    }
    
    # loop through by extent
    fmod <- response.ext[grep(models2[i], f.ext[response.ext])]
    load(file.path(in.base, in.ext, f.ext[fmod]))
    
    # Note: because we're lumping everything together, let's not mess with reiterating the base level
    # ci.terms  <- rbind(ci.terms, mod.out$ci.terms[!(mod.out$ci.terms$Resolution=="t.001" & substr(mod.out$ci.terms$Extent,1,3)=="850"),] )
    dat.ecosys <- rbind(dat.ecosys,
                        cbind(mod.out$data[!(mod.out$data$Resolution=="t.001" & substr(mod.out$data$Extent,1,3)=="850"),], 
                              mod.out$ci.response[!(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850"),c("mean", "lwr", "upr")]))
	  # wt.terms   <- merge(wt.terms, mod.out$weights[!(mod.out$weights$Resolution=="t.001" & substr(mod.out$weights$Extent,1,3)=="850"),], all.x=T, all.y=T)
    # # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    # !(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850")
    # }
    
    # Clear the mod.out to save space
    rm(mod.out)
  }
  #ci.terms$Response <- response.all[r]
  
#}# Fix Extent Labels for consistency
# ci.terms$Extent   <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
# wt.terms$Extent   <- as.factor(ifelse(wt.terms$Extent=="850-2010", "0850-2010", paste(wt.terms$Extent)))
# summary(ci.terms)
summary(dat.ecosys)
# summary(wt.terms)

# Get rid of Linkages, because it's just weird
# ci.terms   <- ci.terms[!ci.terms$Model=="linkages",]
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="linkages",]
# wt.terms   <- wt.terms[!wt.terms$Model=="linkages",]
ecosys     <- ecosys[!ecosys$Model=="linkages",]

# Read in tree-ring data
dat.tr <- read.csv("Data/analysis_response_scale_TreeRings/TreeRings_PHA.csv")
summary(dat.tr)

# Adding a random scalar to show uncertainty in magnitude of tree rings
plotIDs <- unique(dat.tr$PlotID)
for(i in 1:length(plotIDs)){
	dat.tr[dat.tr$PlotID==plotIDs[i],"NPP.Graph"] <- dat.tr[dat.tr$PlotID==plotIDs[i],"NPP"] + i*2
}
summary(dat.tr)

# ----------------------------------------
# Some Plotting etc of NPP
# ----------------------------------------
summary(dat.ecosys)

models.use <- unique(dat.ecosys[,"Model.Order"])
colors.use <- as.vector(model.colors[model.colors$Model.Order %in% models.use, "color"])


pdf(file.path(fig.dir, "NPP_Scales_PHA_0850-2010_Simple_TreeRings.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",])  +
	# geom_ribbon(data=extent.box, aes(x=Year, ymin=Min, ymax=Max), alpha=0.2) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.5, alpha=0.35) +
	geom_line(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=1) +
	geom_point(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=5) +
	geom_line(data=dat.tr, aes(x=Year, y=NPP.Graph, color=Model.Order2), size=1.5, alpha=0.95) +
	geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(limits=c(0,20), expand=c(0,0)) +
	scale_color_manual(values=c(colors.use, rep("chartreuse", length(plotIDs)))) +
	labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
	guides(col=guide_legend(nrow=2)) +
	theme(legend.position="top") +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1.25)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
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

