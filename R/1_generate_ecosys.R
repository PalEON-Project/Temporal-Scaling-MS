# ----------------------------------------
# Generate the model output table ("ecosys") that is be basis of most analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 2 June 2015 
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ncdf4)
library(lme4)
library(R2jags)
library(ggplot2); library(grid)
library(car)
library(zoo)
# library(mvtnorm)
# library(MCMCpack)

sec2yr <- 1*60*60*24*365.25 # 1 sec * 60 sec/min * 60 min/hr * 24 hr/day * 365.25 days/yr
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Dropbox/PalEON CR/paleon_mip_site")
inputs    <- "phase1a_output_variables"
path.data <- "Analyses/Temporal-Scaling/Data"
fig.dir   <- "Analyses/Temporal-Scaling/Figures"
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
ecosys$precipf <- ecosys$precipf*sec2yr # convert precip to mm/yr just to help people figure it out
summary(ecosys)

# CO2 Record
nc.co2 <- nc_open("env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)

# Merging CO2 into Model Outputs
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors

model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors
# ----------------------------------------

# ----------------------------------------
# Standardizations:
# Calculate Deviations from desired reference point
# Reference Point: 0850-0869 (spinup climate)
# ----------------------------------------
vars.response <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap")
vars.climate <- c("tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair", "CO2")

vars <- c(vars.response, vars.climate)
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
		for(v in unique(vars.climate)){
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
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair", "CO2")
vars.dev <- c("AGB.diff", paste0(vars[1:(length(vars)-3)], ".dev"), paste0(vars.climate, ".abs.dev"))

vars.all <- c(vars, vars.dev)

# A lazy way of adding a Scale factor (this probably could be done more simply with a merge, but this works too)
ecosys.010 <- ecosys.050 <- ecosys.100 <- ecosys.250 <- ecosys
ecosys$Scale     <- as.factor("t.001")
ecosys.010$Scale <- as.factor("t.010")
ecosys.050$Scale <- as.factor("t.050")
ecosys.100$Scale <- as.factor("t.100")
ecosys.250$Scale <- as.factor("t.250")
ecosys <- rbind(ecosys, ecosys.010, ecosys.050, ecosys.100, ecosys.250)
summary(ecosys)

# Doing the scale by model & site
for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		for(v in vars.all){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Scale=="t.001", v]

			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Scale=="t.010", v] <- rollapply(temp, FUN=mean, width=10, align="right", fill=NA)
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Scale=="t.050", v] <- rollapply(temp, FUN=mean, width=50, align="right", fill=NA)
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Scale=="t.100", v] <- rollapply(temp, FUN=mean, width=100, align="right", fill=NA)
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Scale=="t.250", v] <- rollapply(temp, FUN=mean, width=250, align="right", fill=NA)
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
ecosys <- ecosys[!ecosys$Model=="clm.bgc",]
col.model <- paste(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])
# col.model <- paste(model.colors[,"color"])

# ---------------------------
# Plotting NPP by model & site
# ---------------------------
pdf(file.path(fig.dir, "NPP_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=NPP, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting AGB by model & site
# ---------------------------
pdf(file.path(fig.dir, "AGB_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001",], aes(x=Year, y=AGB, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100",], aes(x=Year, y=AGB, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Site=="PHA",], aes(x=Year, y=AGB, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting tair by model & site
# ---------------------------
pdf(file.path(fig.dir, "tair_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=tair, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=tair, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "tair_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=tair, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting precipf by model & site
# ---------------------------
pdf(file.path(fig.dir, "precipf_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=precipf, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=precipf, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "precipf_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=precipf, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------

# ---------------------------
# Plotting CO2 by model & site
# ---------------------------
pdf(file.path(fig.dir, "CO2_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=CO2, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=CO2, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "CO2_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(data=ecosys[ecosys$Scale=="t.001" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=0.25, alpha=0.3) +
	geom_line(data=ecosys[ecosys$Scale=="t.100" & ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], aes(x=Year, y=CO2, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()
# ---------------------------
