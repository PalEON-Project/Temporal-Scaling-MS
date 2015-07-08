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
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/Dropbox/PalEON CR/paleon_mip_site")
inputs <- "phase1a_output_variables"
path.data <- "~/Desktop/Dropbox/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/Data"
fig.dir <- "~/Desktop/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
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
# Something's weird with Transpiration, so lets just get rid of it
ecosys <- ecosys[,!(names(ecosys)=="Transp")]
summary(ecosys)

# CO2 Record
nc.co2 <- nc_open("~/Desktop/Dropbox/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)

# Merging CO2 into Model Outputs
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("~/Desktop/Dropbox/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
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
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap")
vars.climate <- c("Temp", "Precip", "CO2")

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
# Note: Smoothing is performed over the PREVIOUS 100 years becuase ecosystems 
#       cannot respond to what they have not yet experienced
# ----------------------------------------
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Temp", "Precip", "CO2")
vars.dev <- c("AGB.diff", paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		# -----------------------
		# 10-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollapply(temp, FUN=mean, width=10, align="right", fill=NA)
		}

		## Standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollapply(temp, FUN=mean, width=10, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 50-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollapply(temp, FUN=mean, width=50, align="right", fill=NA)
		}

		## Standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".50")] <- rollapply(temp, FUN=mean, width=50, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 100-yr Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollapply(temp, FUN=mean, width=100, align="right", fill=NA)
		}

		## Standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollapply(temp, FUN=mean, width=100, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# 250-Year Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollapply(temp, FUN=mean, width=250, align="right", fill=NA)
		}

		## Standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".250")] <- rollapply(temp, FUN=mean, width=250, align="right", fill=NA)
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

# Plotting NPP by model & site
pdf(file.path(fig.dir, "NPP_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(aes(x=Year, y=NPP.100, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "NPP_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=NPP, color=Model), size=0.25, alpha=0.3) +
	geom_line(aes(x=Year, y=NPP.100, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

# ---------------------------
# Plotting AGB by model & site
pdf(file.path(fig.dir, "AGB_Annual_AllSites_AllModels.pdf"))
ggplot(data=ecosys) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model)) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_Century_AllSites_AllModels.pdf"))
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line (aes(x=Year, y=AGB, color=Model), size=0.5, alpha=0.3) +
	geom_line(aes(x=Year, y=AGB.100, color=Model), size=1.5) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()

pdf(file.path(fig.dir, "AGB_Annual_PHA_AllModels.pdf"))
ggplot(data=ecosys[ecosys$Site=="PHA",]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model), size=0.5, alpha=0.3) +
	geom_line(aes(x=Year, y=AGB.100, color=Model), size=2) +
	scale_color_manual(values=col.model) +
	theme_bw()
dev.off()