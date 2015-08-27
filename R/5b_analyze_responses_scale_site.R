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
in.res  <- "Big4_GS_byResolution_Site"
# in.ext  <- "Big4_GS_byExtent"
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
# f.ext <- dir(file.path(in.base, in.ext))

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
    
    # # loop through by extent
    # fmod <- response.ext[grep(models2[i], f.ext[response.ext])]
    # load(file.path(in.base, in.ext, f.ext[fmod]))
    
    # # Note: because we're lumping everything together, let's not mess with reiterating the base level
    # ci.terms  <- rbind(ci.terms, mod.out$ci.terms[!(mod.out$ci.terms$Resolution=="t.001" & substr(mod.out$ci.terms$Extent,1,3)=="850"),] )
    # dat.ecosys <- rbind(dat.ecosys,
                        # cbind(mod.out$data[!(mod.out$data$Resolution=="t.001" & substr(mod.out$data$Extent,1,3)=="850"),], 
                              # mod.out$ci.response[!(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850"),c("mean", "lwr", "upr")]))
	  # wt.terms   <- merge(wt.terms, mod.out$weights[!(mod.out$weights$Resolution=="t.001" & substr(mod.out$weights$Extent,1,3)=="850"),], all.x=T, all.y=T)
    # # # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
    # # !(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850")
    # # }
    
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
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Across all scales (res + extent) finding the mean NPP
# NOTE: not worrying about Site because the terms are across sites
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){
		for(e in unique(ci.terms[ci.terms$Model==m & ci.terms$Resolution==r, "Site"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# -----------------------
			# Find the start year for the extent
			# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))
			yr <- 850

			npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=yr & dat.ecosys$Resolution==r & dat.ecosys$Site==e, "NPP"], na.rm=T)			# -----------------------
			
			# -----------------------
			# Relativizing everythign in dat.ecosys to make it comparable to tree rings
			# -----------------------
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"NPP"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"lwr"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==e,"upr"])/npp
			# -----------------------

			
			# -----------------------
			# Finding the percent change in NPP relative to the mean for that particular scale
			# -----------------------
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==e,"upr"])/npp
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
ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
summary(ci.terms.graph)

pdf(file.path(fig.dir, "NPP_RelChange_Models_by_Site.pdf"))
for(e in unique(ci.terms.graph$Effect)){
print(
ggplot(data= ci.terms.graph[ci.terms.graph$Effect==e & ci.terms.graph$Resolution=="t.001",]) + 
	facet_wrap(~Model, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Site), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color=Site), size=0.75) +
	# scale_fill_manual(values=colors.use) +
	# scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title=paste0(e, " sensitivity, Models by Site")) +
	theme_bw()
)
}
dev.off()

pdf(file.path(fig.dir, "NPP_RelChange_Sites_by_Model.pdf"))
for(e in unique(ci.terms.graph$Effect)){
print(
ggplot(data= ci.terms.graph[ci.terms.graph$Effect==e & ci.terms.graph$Resolution=="t.001",]) + 
	facet_wrap(~Site, scales="free_x") +
	geom_ribbon(aes(x=x, ymin=lwr.rel, ymax=upr.rel, fill=Model), alpha=0.5) +
	geom_line(aes(x=x, y=mean.rel, color= Model), size=0.75) +
	scale_fill_manual(values=colors.use) +
	scale_color_manual(values=colors.use) +
	labs(y="% Change NPP", title=paste0(e, " sensitivity, Sites by Model")) +
	theme_bw()
)
}
dev.off()
