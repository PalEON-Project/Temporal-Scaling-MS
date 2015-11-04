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
in.res  <- "TPC_YR_byResolution_Site"
# in.ext  <- "Big4_GS_byExtent"
out.dir <- "Data/analysis_response_ensemble_TPC"
fig.dir <- "Figures/analysis_response_ensemble_TPC"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData_Raw.Rdata"))
ecosys <- ecosys[!substr(ecosys$Model, 1, 3)=="clm" & !ecosys$Model=="linkages",]

# Figure out what models we have to work with
#response.all <- c("NPP", "AGB.diff", "NEE")
response.all <- c("NPP")
models <- unique(ecosys$Model)
f.res <- dir(file.path(in.base, in.res))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")

# Put all responses & models & scales into single data frames
# for(r in 1:length(response.all)){
  #response <- response.all[r]
  response="NPP"
  response.res <- grep(response, f.res)
  # response.ext <- grep(response, f.ext)
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

summary(wt.terms)

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


pdf(file.path(fig.dir, "NPP_Scales_PHA_0850-2010_Simple.pdf"), height=8.5, width=11)
print(
ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",])  +
	# geom_ribbon(data=extent.box, aes(x=Year, ymin=Min, ymax=Max), alpha=0.2) +
	geom_line(aes(x=Year, y=NPP, color=Model.Order), size=0.5, alpha=0.35) +
	geom_line(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=2) +
	# geom_point(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.050" & dat.ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model.Order), size=5) +
	# geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_x_continuous(limits=c(0850, 2010), expand=c(0,0)) +
	scale_y_continuous(limits=c(0,20), expand=c(0,0)) +
	scale_color_manual(values=colors.use) +
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


# pdf(file.path(fig.dir, "NPP_ModelFits_PHA.pdf"))
# print(
# ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	# geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	# geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.5, color="gray40") +
	# scale_color_manual(values=colors.use) +
	# # scale_x_continuous(limits=c(1900,2010)) +
	# ggtitle("PHA NPP 0850-2010") +
	# theme_bw() + guides(color=F)
# )
# print(
# ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	# geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	# geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	# scale_color_manual(values=colors.use) +
	# scale_x_continuous(limits=c(1900,2010)) +
	# ggtitle("PHA NPP 1900-2010") +
	# theme_bw() + guides(color=F)
# )
# print(
# ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	# geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	# geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	# scale_color_manual(values=colors.use) +
	# scale_x_continuous(limits=c(1800,1820)) +
	# ggtitle("PHA NPP 1800-1820") +
	# theme_bw() + guides(color=F)
# )
# print(
# ggplot(data=dat.ecosys[dat.ecosys$Extent=="0850-2010" & dat.ecosys$Resolution=="t.001" & dat.ecosys$Site=="PHA",]) + facet_wrap(~Model.Order) +
	# geom_line(aes(x=Year, y=NPP, color=Model.Order), size=1) +
	# geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5) +
	# geom_line(aes(x=Year, y=fit.gam), alpha=0.8, size=0.8, color="gray40") +
	# scale_color_manual(values=colors.use) +
	# scale_x_continuous(limits=c(1990,2010)) +
	# ggtitle("PHA NPP 1990-2010") +
	# theme_bw() + guides(color=F)
# )
# dev.off()
# ----------------------------------------






# ----------------------------------------
# Compare driver responses across models by standardizing driver responses to the mean model NPP
# ----------------------------------------
# Across all scales (resolution) finding the mean NPP
# NOTE: we ARE relativizing per site here since the response curves were site-specific
for(m in unique(ci.terms$Model)){
	for(r in unique(ci.terms[ci.terms$Model==m, "Resolution"])){
		for(s in unique(ci.terms[ci.terms$Model==m & ci.terms$Resolution==r, "Site"])){

			# -----------------------
			# Find the NPP to relativize each set off of
			# Using 1800:1850 mean (Settlement-era) as reference point
			# -----------------------
			# Find the start year for the extent
			# yr <- ifelse(nchar(as.character(e))==8, as.numeric(substr(e,1,3)), as.numeric(substr(e,1,4)))

			npp <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year %in% 1800:1850 & dat.ecosys$Resolution==r & dat.ecosys$Site==s, "NPP"], na.rm=T)			# -----------------------
			
			# -----------------------
			# Relativizing everythign in dat.ecosys to make it comparable to tree rings
			# -----------------------
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"NPP.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"NPP"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"fit.gam.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"mean.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"mean"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"lwr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"lwr"])/npp
			dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"upr.rel"] <- (dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Resolution==r & dat.ecosys$Site==s,"upr"])/npp
			# -----------------------

			
			# -----------------------
			# Finding the percent change in NPP relative to the mean for that particular scale
			# -----------------------
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"mean.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"mean"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"lwr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"lwr"])/npp
			ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"upr.rel"] <- (ci.terms[ci.terms$Model==m & ci.terms$Resolution==r & ci.terms$Site==s,"upr"])/npp
			# -----------------------

			# -----------------------
			# Relativizing the factor fits through times and weights as well
			# Note: because a fit of 0 means no change from the mean, we need to add 1 to all of these
			# -----------------------
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r & wt.terms$Site==s,"fit.tair.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r & wt.terms$Site==s,"fit.tair"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r & wt.terms$Site==s,"fit.precipf.rel"] <- 1+ (wt.terms[wt.terms$Model==m & wt.terms$Resolution==r & wt.terms$Site==s,"fit.precipf"])/npp
			wt.terms[wt.terms$Model==m & wt.terms $Resolution==r & wt.terms$Site==s,"fit.CO2.rel"] <- 1+(wt.terms[wt.terms$Model==m & wt.terms$Resolution==r & wt.terms$Site==s,"fit.CO2"])/npp
			# -----------------------

		}
	}
}
summary(dat.ecosys)
summary(ci.terms)
summary(wt.terms)

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


indices.wt <- wt.terms$Site=="PHA" & wt.terms$Resolution=="t.001"
indices.dat <- dat.ecosys$Site=="PHA" & dat.ecosys$Resolution=="t.001"

# pdf(file.path(fig.dir, "NPP_Drivers_Time_PHA_1800-2010.pdf"), width=11, height=8.5)
# print(
# ggplot(data=wt.terms[indices.wt,]) + facet_wrap(~Model.Order) +
 	# geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="clm.cn",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="clm.cn","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="ed2","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="ed2.lu",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="ed2.lu","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.stat",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="jules.stat","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="jules.triffid",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="jules.triffid","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.guess",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="lpj.guess","weight.precipf"])), size=3) +
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="lpj.wsl",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="lpj.wsl","weight.precipf"])), size=3)+
	# geom_line(data=wt.terms[indices.wt & wt.terms$Model=="sibcasa",], aes(x=Year, y=fit.full),
	          # color=rgb(abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.tair"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.CO2"]),
                        # abs(wt.terms[indices.wt & wt.terms$Model=="sibcasa","weight.precipf"])), size=3) +
	# scale_x_continuous(limits=c(1700,2010), expand=c(0,0), breaks=c(1750, 1850, 1950)) +
	# scale_y_continuous(name=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")"))), expand=c(0,0)) +
	# theme(axis.line=element_line(color="black", size=0.5), 
	      # panel.grid.major=element_blank(), 
	      # panel.grid.minor=element_blank(), 
	      # panel.border=element_blank(), 
	      # panel.background=element_blank(), 
	      # axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      # axis.text.y=element_text(color="black", size=rel(2.5)), 
	      # axis.title.x=element_text(face="bold", size=rel(2.25), vjust=-0.5),  
	      # axis.title.y=element_text(face="bold", size=rel(2.25), vjust=1)) +
	# theme(strip.text=element_text(size=rel(1.5)))
# )
# dev.off()

# Merging wt.terms & dat.ecosys to get the relativized NPP lined up
wt.terms2 <- merge(wt.terms, dat.ecosys[,c("Model", "Model.Order", "Site", "Year", "Resolution", "Extent", "NPP", "NPP.rel")], all.x=T, all.y=F)
summary(wt.terms2)

indices.wt2 <- wt.terms2$Site=="PHA" & wt.terms2$Resolution=="t.001" & wt.terms2$Extent=="0850-2010" 
indices.wt2b <- wt.terms2$Site=="PHA" & wt.terms2$Resolution=="t.010" & wt.terms2$Extent=="0850-2010" 


pdf(file.path(fig.dir, "NPP_Rel_Drivers_Time_PHA_1800-2010.pdf"), width=11, height=7.5)
# png(file.path(fig.dir, "NPP_Rel_Drivers_Time_PHA_1800-2010.png"), width=11, height=8.5, units="in", res=600)
print( 
ggplot() + facet_grid(Model.Order~., scales="free_y") +
 	# geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="clm.cn",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="clm.cn","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="ed2",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="ed2.lu",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="ed2.lu","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="jules.stat",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.stat","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="jules.triffid",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="jules.triffid","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.guess",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.guess","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.wsl",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="lpj.wsl","weight.precipf"])), size=3)+
	geom_line(data=wt.terms2[indices.wt2 & wt.terms2$Model=="sibcasa",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2 & wt.terms2$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms2[indices.wt2 & wt.terms2$Model=="sibcasa","weight.precipf"])), size=3) +
	scale_x_continuous(limits=c(1800,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("% Mean NPP"))), expand=c(0,0)) +
	ggtitle("NPP Controlling Factor -- Annual") + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(2.5))) + 
	theme(strip.text=element_text(size=rel(0.75), face="bold")) +
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
print( 
ggplot() + facet_grid(Model.Order~., scales="free_y") +
 	# geom_line(data= dat.ecosys[indices.dat,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="clm.cn",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="clm.cn","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="clm.cn","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="clm.cn","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="ed2",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="ed2.lu",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2.lu","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2.lu","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="ed2.lu","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="jules.stat",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.stat","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.stat","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.stat","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="jules.triffid",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.triffid","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.triffid","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="jules.triffid","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.guess",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.guess","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.guess","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.guess","weight.precipf"])), size=3) +
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.wsl",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.wsl","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.wsl","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="lpj.wsl","weight.precipf"])), size=3)+
	geom_line(data=wt.terms2[indices.wt2b & wt.terms2$Model=="sibcasa",], aes(x=Year, y= NPP.rel*100),
	          color=rgb(abs(wt.terms2[indices.wt2b & wt.terms2$Model=="sibcasa","weight.tair"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="sibcasa","weight.CO2"]),
                        abs(wt.terms2[indices.wt2b & wt.terms2$Model=="sibcasa","weight.precipf"])), size=3) +
	scale_x_continuous(limits=c(1800,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("% Mean NPP"))), expand=c(0,0)) +
	ggtitle("NPP Controlling Factor -- Decadal") + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(2.5))) + 
	theme(strip.text=element_text(size=rel(0.75), face="bold")) +
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
dev.off()



# ---------------------------------
# Looking at the ensemble mean of Drivers of change and sensitivities through time
# ---------------------------------
summary(wt.terms2)
ensemble.base.wts   <- wt.terms2[wt.terms2$Resolution=="t.010",]
ensemble.base.terms <- ci.terms[ci.terms$Resolution=="t.010",]

summary(ensemble.base.wts)
summary(ensemble.base.terms)

ensemble.wts1 <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Year), FUN=mean)
names(ensemble.wts1)[1:2] <- c("Site", "Year") 
summary(ensemble.wts1)


ensemble.wts.lo <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Year), FUN=quantile, 0.025)
names(ensemble.wts.lo) <- c("Site", "Year", paste0(names(ensemble.wts.lo[3:ncol(ensemble.wts.lo)]), ".lo")) 
summary(ensemble.wts.lo)

ensemble.wts.hi <- aggregate(ensemble.base.wts[,c("fit.full", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "NPP.rel")], by=list(ensemble.base.wts$Site, ensemble.base.wts$Year), FUN=quantile, 0.975)
names(ensemble.wts.hi) <- c("Site", "Year", paste0(names(ensemble.wts.hi[3:ncol(ensemble.wts.hi)]), ".hi")) 
summary(ensemble.wts.hi)

dim(ensemble.wts1); dim(ensemble.wts.lo); dim(ensemble.wts.hi)

ensemble.wts <- cbind(ensemble.wts1, ensemble.wts.lo[,2:ncol(ensemble.wts.lo)], ensemble.wts.hi[,2:ncol(ensemble.wts.hi)])
summary(ensemble.wts)

ensemble.wts2 <- stack(ensemble.wts[,c("fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel")])[,c(2,1)]
names(ensemble.wts2) <- c("Driver", "fit.rel")
ensemble.wts2$Driver <- as.factor(substr(ensemble.wts2$Driver, 5, nchar(paste(ensemble.wts2$Driver))-4 ))
levels(ensemble.wts2$Driver) <- c("CO2", "Precip", "Temp")
ensemble.wts2[,c("Site", "Year")] <- ensemble.wts[,c("Site", "Year")]
ensemble.wts2$fit.rel.lo <- stack(ensemble.wts[,c("fit.tair.rel.lo", "fit.precipf.rel.lo", "fit.CO2.rel.lo")])[,1]
ensemble.wts2$fit.rel.hi <- stack(ensemble.wts[,c("fit.tair.rel.hi", "fit.precipf.rel.hi", "fit.CO2.rel.hi")])[,1]
summary(ensemble.wts2)


pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_AllSites_1800-2010_Decadal.pdf"), width=11, height=7.5)
# png(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_AllSites_1800-2010_Decadal.png"), width=1100, height=750)
print( 
ggplot() + facet_grid(Site~., scales="free_y") +
 	geom_ribbon(data= ensemble.wts, aes(x=Year, ymin=NPP.rel.lo*100, ymax=NPP.rel.hi*100), alpha=0.5) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PHO",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PHO","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHO","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHO","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PHA",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PMB",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PMB","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PMB","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PMB","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PUN",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PUN","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PUN","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PUN","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PBL",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PBL","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PBL","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PBL","weight.precipf"])), size=4) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PDL",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PDL","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PDL","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PDL","weight.precipf"])), size=4) +
	scale_x_continuous(limits=c(1800,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("% Mean NPP"))), expand=c(0,0)) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(2.5))) + 
	theme(strip.text=element_text(size=rel(0.75), face="bold")) +
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_blank(), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
dev.off()

plot.ensemble.npp <- ggplot() +
 	geom_ribbon(data= ensemble.wts[ensemble.wts$Site=="PHA",], aes(x=Year, ymin=NPP.rel.lo*100, ymax=NPP.rel.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts[ensemble.wts$Site=="PHA",], aes(x=Year, y=NPP.rel*100),
	          color=rgb(abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.tair"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.CO2"]),
                        abs(ensemble.wts[ensemble.wts$Site=="PHA","weight.precipf"])), size=3) +
	geom_hline(y=100, linetype="dashed") +
	scale_x_continuous(limits=c(1800,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_blank(),
		  axis.text.y=element_text(size=rel(2), color="black"), 
		  axis.title.x=element_blank(),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(1,1,-0.9,0), "lines"))


# NEED TO ADD LEGEND!!!
plot.ensemble.weights <- ggplot() +
 	geom_ribbon(data= ensemble.wts2[ensemble.wts2$Site=="PHA",], aes(x=Year, ymin=fit.rel.lo*100, ymax= fit.rel.hi*100, fill=Driver), alpha=0.35) +
 	geom_line(data= ensemble.wts2[ensemble.wts2$Site=="PHA",], aes(x=Year, y=fit.rel*100, color=Driver, linetype=Driver), size=1.5) +
	scale_x_continuous(limits=c(1800,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0), breaks=c(100, 125)) +
	scale_fill_manual(values=c("green","blue", "red")) +
	scale_color_manual(values=c("green4", "blue3", "red3")) +
	scale_linetype_manual(values=c("dashed", "solid", "longdash")) +
	# ggtitle("NPP Controlling Factor") + 
	theme(legend.text=element_text(size=rel(1.5)), 
	      legend.title=element_text(size=rel(1.5)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines"),
	      legend.position=c(0.25, 0.8),
	      legend.direction="horizontal") + 
	      # legend.key.width=unit(2, "lines")) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(),
	      panel.margin.y=unit(0.5, "lines"))  +
	theme(axis.text.x=element_text(color="black", size=rel(2.5)),
		  axis.text.y=element_text(color="black", size=rel(2)), 
		  axis.title.x=element_text(size=rel(2.5), face="bold"),  
		  axis.title.y=element_text(size=rel(2), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines")) +
	      theme(plot.margin=unit(c(0,1,0.5,0), "lines"))
	      
	      
# pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_PHA_1800-2010.pdf"), width=11, height=7.5)
# print(plot.ensemble.npp)
# dev.off()

# pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Weights_Time_PHA_1800-2010.pdf"), width=11, height=7.5)
# print( plot.ensemble.weights)
# dev.off()

# 2-panel Figure!!
pdf(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_PHA_1800-2010_Decadal.pdf"), width=11, height=7.5)
# png(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_PHA_1800-2010_Decadal.png"), width=800, height=600)
# svg(file.path(fig.dir, "NPP_Ensemble_Drivers_Time_PHA_1800-2010_Decadal.svg"), width=8, height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
print(plot.ensemble.npp    , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.ensemble.weights, vp = viewport(layout.pos.row = 2, layout.pos.col=1))
dev.off()
# ----------------------------------------






