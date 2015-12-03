# ----------------------------------------
# Generate the model output table ("ecosys") that is be basis of most analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 2 June 2015 
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(car)
library(zoo)

# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
# setwd("..")

inputs    <- "~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/phase1a_output_variables"
path.data <- "Data"
fig.dir   <- "Figures"

# Making sure the appropriate file paths exist
if(!dir.exists(path.data)) dir.create(path.data)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------


# ----------------------------------------
# Note: Commented out because saved as EcosysData.RData 1 June 2015
#       (with an increasing number of models, running this every time became cumbersome)
# ----------------------------------------
# Load Data Sets
# ----------------------------------------
# Ecosystem Model Outputs
ecosys <- read.csv(file.path(inputs, "MIP_Data_Ann_2015.csv"))
ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
summary(ecosys)

# Note: Now we're going to merge in the raw annual and growing season data from the drivers
met.yr <- read.csv(file.path(path.data, "analysis_met_drivers", "Drivers_Year_GrowingSeason.csv"))
met.yr$CO2 <- met.yr$CO2.yr 

# Merge the met drivers into ecosys; 
# NOTE: this merge screws up the year/site/model ordering, so we need to resort
ecosys <- merge(ecosys, met.yr, all.x=T, all.y=T)
ecosys <- ecosys[order(ecosys$Model, ecosys$Site, ecosys$Year),]
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("raw_inputs/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors

model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors
# ----------------------------------------


# ----------------------------------------
# Convert variables to more user-friendly units
# ----------------------------------------
# Fluxes: change s-1 to yr-1
sec2yr <- 1*60*60*24*365.25 # 1 sec * 60 sec/min * 60 min/hr * 24 hr/day * 365.25 days/yr
vars.flux <- c("precipf", "precipf.yr", "precipf.gs", "GPP", "NPP", "NEE", "AutoResp", "HeteroResp", "Evap")
ecosys[,vars.flux] <- ecosys[,vars.flux]*sec2yr

# Growing Season Precipitation should be in total precip for those months (may - sept = 5 months)
ecosys$precipf.gs = ecosys$precipf.gs/12*5

# Carbon: change kgC/m2 to MgC/ha
kgm2_2_Mgha <- 1*0.001*10000 # 1 kg * 0.001 kg/Mg * 10000 m2/ha
vars.carbon <- c("GPP", "NPP", "NEE", "AutoResp", "HeteroResp", "AGB", "SoilCarb")
ecosys[, vars.carbon] <- ecosys[, vars.carbon]*kgm2_2_Mgha

# Calculate dAGB because it could come in handy later
for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){

		# -----------------------
		# AGB 1st difference		
		# -----------------------
		ecosys[ecosys$Site==s & ecosys$Model==m,"AGB.diff"] <- c(NA, diff(ecosys[ecosys$Site==s & ecosys$Model==m,"AGB"]))
	}
}

save(ecosys, model.colors, file=file.path(path.data, "EcosysData_Raw.Rdata"))
# ----------------------------------------

# ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1990,]) +
	# geom_line(aes(x=Year, y=AGB.diff, color=Model.Order)) +
	# # scale_color_manual(values=model.colors$color) +
	# # scale_y_continuous(limits=c(-1,1)) +
	# # scale_x_continuous(limits=c(1990,2010)) +
	# theme_bw()

# ----------------------------------------
# Standardizations:
# Calculate Deviations from desired reference point
# Reference Point: 0850-0869 (spinup climate)
# ----------------------------------------
vars.response <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Evergreen","Deciduous", "Grass")
vars.climate <- c("tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair", "CO2")
vars.climate2 <- c(vars.climate, paste0(vars.climate, ".yr"), paste0(vars.climate, ".gs"))
vars <- c(vars.response, vars.climate2)
# vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

ref.window <- 850:869

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){

		# -----------------------
		# AGB 1st difference		
		# -----------------------
		ecosys[ecosys$Site==s & ecosys$Model==m,"AGB.diff"] <- c(NA, diff(ecosys[ecosys$Site==s & ecosys$Model==m,"AGB"]))

		# -----------------------
		# Model Variabiles -- Relative Change
		# Deviation = percent above or below the mean for the reference window
		#             (observed-ref.mean)/ref.mean 
		# -----------------------
		for(v in unique(vars)){
			ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".dev")] <- (ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean)/ref.mean

		}
		# -----------------------

		# -----------------------
		# Climate Drivers -- Absolute, not relative change
		# Deviation = absolute deviation from reference window
		#             observed - ref.mean
		# -----------------------
		for(v in unique(vars.climate2)){
			ref.mean <- mean(ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= min(ref.window) & ecosys$Year<=max(ref.window), v], na.rm=T)
			ecosys[ecosys$Site==s & ecosys$Model==m, paste0(v, ".abs.dev")] <- ecosys[ecosys$Site==s & ecosys$Model==m, v] - ref.mean
		}
		# -----------------------
	}
}

summary(ecosys)
# ----------------------------------------


# ----------------------------------------
# Perform Temporal Smoothing on Data
# Note: Smoothing is performed over the PREVIOUS X years becuase ecosystems 
#       cannot respond to what they have not yet experienced
# ----------------------------------------
# vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair", "CO2")
vars.dev <- c("AGB.diff", paste0(vars, ".dev"), paste0(vars.climate2, ".abs.dev"))

vars.all <- c(vars, vars.dev)

# A lazy way of adding a Scale factor (this probably could be done more simply with a merge, but this works too)
ecosys.010 <- ecosys.050 <- ecosys.100 <- ecosys.250 <- ecosys
ecosys$Resolution     <- as.factor("t.001")
ecosys.010$Resolution <- as.factor("t.010")
ecosys.050$Resolution <- as.factor("t.050")
ecosys.100$Resolution <- as.factor("t.100")
ecosys <- rbind(ecosys, ecosys.010, ecosys.050, ecosys.100)
summary(ecosys)

# Doing the scale by model & site
for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		for(v in vars.all){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Resolution=="t.001", v]

			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Resolution=="t.010", v] <- rollapply(temp, FUN=mean, width=10, align="center", fill=NA)
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Resolution=="t.050", v] <- rollapply(temp, FUN=mean, width=50, align="center", fill=NA)
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Resolution=="t.100", v] <- rollapply(temp, FUN=mean, width=100, align="center", fill=NA)
		}
		# -----------------------
	}
}
summary(ecosys)
save(ecosys, model.colors, file=file.path(path.data, "EcosysData.Rdata"))
# ----------------------------------------



# ----------------------------------------
# Making Some general Figures that are handy
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata"))

# Note: CLM-BGC is wrong, so lets exclude it for now
ecosys <- ecosys[!ecosys$Model=="linkages",]
col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])
# col.model <- paste(model.colors[,"color"])

# ---------------------------
# Plotting NPP by model & site
# ---------------------------
pdf(file.path(fig.dir, "NPP_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_Century_AllSites_NoLINKAGES.pdf"))
ggplot(data=ecosys[!ecosys$Model=="linkages",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & !ecosys$Model=="linkages",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & !ecosys$Model=="linkages",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	scale_color_manual(values=col.model) + labs(y="NPP MgC ha-1 yr-1") +
	theme_bw()
dev.off()


pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting AGB by model & site
# ---------------------------
pdf(file.path(fig.dir, "AGB_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=AGB, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=AGB, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_Century_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color=Model), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()


ggplot(data=ecosys[ecosys$Site=="PHA",]) + #facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color= Model.Order), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=rel(3))) + theme(legend.position=c(0.6,0.9), legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) + labs(color="Model", y=expression(bold(paste("Aboveground Biomass (Mg C ha"^"-1",")"))), title="Comparison of Model Aboveground Biomass") +
	guides(color=guide_legend(ncol=3)) +
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=rel(2)), axis.text.y=element_text(color="black", size=rel(2)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.margin=unit(c(0.1,0.5,0.5,0.1), "lines"))

# ---------------------------


# ---------------------------
# Plotting Fraction Evergreen by model & site
# ---------------------------
pdf(file.path(fig.dir, "Evergreen_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Evergreen, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Evergreen_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Evergreen, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Evergreen_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Evergreen, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=Evergreen, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Evergreen_Annual_Century_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=Evergreen, color=Model), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=Evergreen, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

# ggplot(data=ecosys[ecosys$Site=="PHA",]) + #facet_wrap(~Site) +
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Evergreen, color=Model.Order), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.050",], aes(x=Year, y=Evergreen, color= Model.Order), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=rel(3))) + theme(legend.position="bottom", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) + labs(color="Model", y="Fraction Evergreen") +
	guides(color=guide_legend(ncol=3)) +
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=rel(2)), axis.text.y=element_text(color="black", size=rel(2)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.margin=unit(c(0.1,0.5,0.5,0.1), "lines"))

# ---------------------------

# ---------------------------
# Plotting Fraction Deciduous by model & site
# ---------------------------
pdf(file.path(fig.dir, "Deciduous_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Deciduous, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Deciduous_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Deciduous, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Deciduous_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Deciduous, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=Deciduous, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Deciduous_Annual_Century_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=Deciduous, color=Model), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=Deciduous, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Deciduous, color=Model.Order), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.050",], aes(x=Year, y=Deciduous, color= Model.Order), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=rel(3))) + theme(legend.position="bottom", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) + labs(color="Model", y="Fraction Deciduous") +
	guides(color=guide_legend(ncol=3)) +
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=rel(2)), axis.text.y=element_text(color="black", size=rel(2)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.margin=unit(c(0.1,0.5,0.5,0.1), "lines"))
# ---------------------------

# ---------------------------
# Plotting Fraction Grass by model & site
# ---------------------------
pdf(file.path(fig.dir, "Grass_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Grass, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Grass_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=Grass, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Grass_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Grass, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100",], aes(x=Year, y=Grass, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "Grass_Annual_Century_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=Grass, color=Model), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=Grass, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001",], aes(x=Year, y=Grass, color=Model.Order), size=1, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.050",], aes(x=Year, y=Grass, color= Model.Order), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=rel(3))) + theme(legend.position="bottom", legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) + labs(color="Model", y="Fraction Grass") +
	guides(color=guide_legend(ncol=3)) +
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=rel(2)), axis.text.y=element_text(color="black", size=rel(2)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.margin=unit(c(0.1,0.5,0.5,0.1), "lines"))
# ---------------------------

# ---------------------------
# Plotting tair by model & site
# ---------------------------
pdf(file.path(fig.dir, "tair_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=tair, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=tair, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting precipf by model & site
# ---------------------------
pdf(file.path(fig.dir, "precipf_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=precipf, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=precipf, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting CO2 by model & site
# ---------------------------
pdf(file.path(fig.dir, "CO2_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=CO2, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=CO2, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Resolution=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Resolution=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------
