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
in.res  <- "AllDrivers_Yr_byResolution"
in.ext  <- "AllDrivers_Yr_byExtent"
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
      # sim.terms <- mod.out$sim.terms
    } else {
      ci.terms  <- rbind(ci.terms, mod.out$ci.terms )
      dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]))
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
    # # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    # !(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850")
    # }
    
    # Clear the mod.out to save space
    rm(mod.out)
  }
  #ci.terms$Response <- response.all[r]
  
#}# Fix Extent Labels for consistency
ci.terms$Extent <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
summary(ci.terms)
summary(dat.ecosys)

# Get rid of Linkages, because it's just weird
ci.terms   <- ci.terms[!ci.terms$Model=="linkages",]
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="linkages",]
ecosys     <- ecosys[!ecosys$Model=="linkages",]


# Write the files to csv so I don't have to mess with loading large .Rdata files again if I don't have to
# write.csv(ci.terms.res,  file.path(out.dir, "Driver_Responses_CI_Resolution.csv"), row.names=F)
# write.csv(sim.terms.res, file.path(out.dir, "Driver_Responses_Sims_Resolution.csv"), row.names=F)
# write.csv(ci.terms,  file.path(out.dir, "Driver_Responses_CI_Extent.csv"), row.names=F)
# write.csv(sim.terms, file.path(out.dir, "Driver_Responses_Sims_Extent.csv"), row.names=F)
# ----------------------------------------


# ----------------------------------------
# Standardize driver responses to the mean model NPP
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


pdf(file.path(fig.dir, "NPP_RelChange_BaseEffect_Yr.pdf"))
ggplot(data=ci.terms[ci.terms$Resolution=="t.001" & ci.terms$Extent=="0850-2010",]) + 
	facet_wrap(~Effect, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel*100, color=Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title="Driver Sensitivity to Annual Met, 850-2010") +
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

