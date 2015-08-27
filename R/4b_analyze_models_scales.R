# ----------------------------------------
# Temporal Scaling Analyses
# Looking at whether different model responses are correlated & how those correlations change with scale
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Are NPP, dAGB, and NEE tightly correlated in the models like we think they are in nature? 
#           Are those correlation scale-dependent
# Rationale: When trying to validate the model output with empirical data, the data available varies by 
#            temporal scale (extent & resolution) and things we'd like to equate may not be the same 
#            because of interactions among ecosystem components as we scale up in levels of organization
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
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
out.dir <- "Data/analysis_models_scales"
fig.dir <- "Figures/analysis_models_scales"

if(!dir.exists(out.dir)) dir.create(out.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load gamm results files
# ----------------------------------------
load(file.path(path.data, "EcosysData_Raw.Rdata"))
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]

# Figure out what models we have to work with
response.all <- c("NPP", "AGB.diff", "NEE")
models <- unique(ecosys$Model)
f.res <- dir(file.path(in.base, in.res))
f.ext <- dir(file.path(in.base, in.ext))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")

# Put all responses & models & scales into single data frames
for(r in 1:length(response.all)){
	response <- response.all[r]
	response.res <- grep(response, f.res)
	response.ext <- grep(response, f.ext)
for(i in 1:length(models)){
	# First narrow to the models
	fmod <- response.res[grep(models2[i], f.res[response.res])]
	if(!length(fmod)>0) next

	load(file.path(in.base, in.res, f.res[fmod]))

	if(r==1 & i==1) {
		ci.terms   <- cbind(mod.out$ci.terms, Response=response)
		dat.ecosys <- cbind(mod.out$data[,!names(mod.out$data)==response], mod.out$ci.response[,c("mean", "lwr", "upr")], Response=response)
		# sim.terms <- mod.out$sim.terms
	} else {
		ci.terms  <- rbind(ci.terms, cbind(mod.out$ci.terms, Response=response))
		dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data[,!names(mod.out$data)==response], mod.out$ci.response[,c("mean", "lwr", "upr")], Response=response))
		# sim.terms <- rbind(sim.terms, mod.out$sim.terms)
	}

	# loop through by extent
	fmod <- response.ext[grep(models2[i], f.ext[response.ext])]
	load(file.path(in.base, in.ext, f.ext[fmod]))
	
	# Note: because we're lumping everything together, let's not mess with reiterating the base level
	ci.terms  <- rbind(ci.terms, cbind(mod.out$ci.terms[!(mod.out$ci.terms$Resolution=="t.001" & substr(mod.out$ci.terms$Extent,1,3)=="850"),], Response=response))
	dat.ecosys <- rbind(dat.ecosys,
				  cbind(mod.out$data[!(mod.out$data$Resolution=="t.001" & substr(mod.out$data$Extent,1,3)=="850"),!names(mod.out$data)==response], 
						mod.out$ci.response[!(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850"),c("mean", "lwr", "upr")], Response=response))
		# # sim.terms <- rbind(sim.terms, mod.out$sim.terms)
		# !(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850")
	# }

	# Clear the mod.out to save space
	rm(mod.out)
}
	# ci.terms$Response <- response.all[r]
	
}# Fix Extent Labels for consistency
ci.terms$Extent <- as.factor(ifelse(ci.terms$Extent=="850-2010", "0850-2010", paste(ci.terms$Extent)))
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
summary(ci.terms)
summary(dat.ecosys)

# Get rid of Linkages, because it's just weird
ci.terms   <- ci.terms[!ci.terms$Model=="linkages",]
dat.ecosys <- dat.ecosys[!dat.ecosys$Model=="linkages",]
ecosys     <- ecosys[!ecosys$Model=="linkages",]
# ----------------------------------------



# ----------------------------------------
# Calculate changes in Drivers across scales
# ----------------------------------------

col.Y <- which(names(dat.ecosys)=="Y")
cols.main <- which(!names(dat.ecosys)=="Response") 

dat.NPP  <- dat.ecosys[dat.ecosys$Response=="NPP",]
names(dat.NPP)[col.Y:ncol(dat.NPP)] <- c("NPP", paste0(names(dat.NPP[(col.Y+1):ncol(dat.NPP)]), ".NPP"))

dat.dAGB <- dat.ecosys[dat.ecosys$Response=="AGB.diff",]
names(dat.dAGB)[col.Y:ncol(dat.dAGB)] <- c("dAGB", paste0(names(dat.dAGB[(col.Y+1):ncol(dat.dAGB)]), ".dAGB"))

dat.NEE  <- dat.ecosys[dat.ecosys$Response=="NEE",]
names(dat.NEE)[col.Y:ncol(dat.NEE)] <- c("NEE", paste0(names(dat.NEE[(col.Y+1):ncol(dat.NEE)]), ".NEE"))

dat.ecosys2 <- merge(dat.NPP[,cols.main], dat.NEE[,cols.main], all.x=T, all.y=T)
dat.ecosys2 <- merge(dat.ecosys2, dat.dAGB[,cols.main], all.x=T, all.y=T)
summary(dat.ecosys2)

# Trimming outliers
dat.ecosys2 <- dat.ecosys2[dat.ecosys2$dAGB>quantile(dat.ecosys2$dAGB, 0.05, na.rm=T) & dat.ecosys2$dAGB<quantile(dat.ecosys2$dAGB, 0.95, na.rm=T) ,]
dat.ecosys2 <- dat.ecosys2[!is.na(dat.ecosys2$Resolution) & !is.na(dat.ecosys2$Extent),]
summary(dat.ecosys2)

pdf(file.path(fig.dir, "ModelVars_Correlation_Resolution_dAGB_NPP_t001.pdf"))
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Extent=="0850-2010" & dat.ecosys2$Resolution=="t.001",]) + 
	facet_wrap( ~ Model.Order, scales="free") +
	geom_point(aes(y=dAGB, x=NPP, color=Site), size=0.8) +
	# scale_y_continuous(limits=c(-1,1)) +
	ggtitle("dAGB vs. NPP") +
	guides(color=guide_legend(nrow=2)) +
	theme(legend.position=c(0.5, 0.2)) +
	theme(plot.title=element_text(face="bold", size=rel(1))) + 
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      plot.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(1))) + 
	theme(strip.text=element_text(size=rel(1), face="bold")) +
	theme(legend.text=element_text(size=rel(1)), 
	      legend.title=element_text(size=rel(1)),
	      legend.key=element_blank(),
	      legend.key.size=unit(1, "lines")) #+ 
	      # legend.key.width=unit(2, "lines")) + 
	# theme(axis.text.x=element_text(color="black", size=rel(1)),
		  # axis.text.y=element_text(color="black", size=rel(1)), 
		  # axis.title.x=element_blank(),  
		  # axis.title.y=element_text(size=rel(1), face="bold"),
		  # axis.ticks.length=unit(-0.5, "lines"),
	      # axis.ticks.margin=unit(1.0, "lines"))
)
dev.off()

pdf(file.path(fig.dir, "ModelVars_Correlation_Resolution.pdf"))
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Extent=="0850-2010",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Resolution ~ Model.Order, scales="free") +
	geom_point(aes(y=dAGB, x=NPP, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("dAGB vs. NPP") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Extent=="0850-2010",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Resolution ~ Model.Order, scales="free") +
	geom_point(aes(y=NEE, x=dAGB, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("NEE vs. dAGB") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Extent=="0850-2010",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Resolution ~ Model.Order, scales="free") +
	geom_point(aes(y=NEE, x=NPP, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("NEE vs. NPP") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
dev.off()

pdf(file.path(fig.dir, "ModelVars_Correlation_Extemt.pdf"))
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Resolution=="t.001",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Extent ~ Model.Order, scales="free") +
	geom_point(aes(y=dAGB, x=NPP, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("dAGB vs. NPP") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Resolution=="t.001",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Extent ~ Model.Order, scales="free") +
	geom_point(aes(y=NEE, x=dAGB, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("NEE vs. dAGB") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
print(
ggplot(data=dat.ecosys2[dat.ecosys2$Resolution=="t.001",]) + 
	# facet_grid(Model.Order ~ Resolution, scales="free") +
	facet_grid(Extent ~ Model.Order, scales="free") +
	geom_point(aes(y=NEE, x=NPP, color=Site), size=0.8) +
	scale_y_continuous(limits=c(-1,1)) +
	ggtitle("NEE vs. NPP") +
	# theme(legend.position="top") #+
	guides(color=F) +
	theme_bw()
)
dev.off()

# Printing/Saving model correlations by site/scale
res.corr.npp.agb <- data.frame(Model=unique(dat.ecosys2$Model))
ext.corr.npp.agb <- data.frame(Model=unique(dat.ecosys2$Model))
for(m in unique(dat.ecosys2$Model)){
	for(r in unique(dat.ecosys2$Resolution)){
		lm.temp <- lm(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution==r & dat.ecosys2$Extent=="0850-2010","dAGB"] ~ dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution==r & dat.ecosys2$Extent=="0850-2010","NPP"])
		res.corr.npp.agb[res.corr.npp.agb$Model==m, r] <- summary(lm.temp)$r.squared
	}
	for(e in unique(dat.ecosys2$Extent)){
		lm.temp <- lm(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Resolution=="t.001","dAGB"] ~ dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution=="t.001" & dat.ecosys2$Extent==e,"NPP"])
		ext.corr.npp.agb[ext.corr.npp.agb$Model==m, e] <- summary(lm.temp)$r.squared
	}
}

res.corr.npp.nee <- data.frame(Model=unique(dat.ecosys2$Model))
ext.corr.npp.nee <- data.frame(Model=unique(dat.ecosys2$Model))
for(m in unique(dat.ecosys2$Model)){
	for(r in unique(dat.ecosys2$Resolution)){
		lm.temp <- lm(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution==r & dat.ecosys2$Extent=="0850-2010","NEE"] ~ dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution==r & dat.ecosys2$Extent=="0850-2010","NPP"])
		res.corr.npp.nee[res.corr.npp.nee$Model==m, r] <- summary(lm.temp)$r.squared
	}
	for(e in unique(dat.ecosys2$Extent)){
		lm.temp <- lm(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Extent==e & dat.ecosys2$Resolution=="t.001","NEE"] ~ dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Resolution=="t.001" & dat.ecosys2$Extent==e,"NPP"])
		ext.corr.npp.nee[ext.corr.npp.nee$Model==m, e] <- summary(lm.temp)$r.squared
	}
}
# ----------------------------------------
