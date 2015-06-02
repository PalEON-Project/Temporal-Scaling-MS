# ----------------------------------------
# Temporal Scaling Analyses
# Non-Stationary Drivers
# Christy Rollinson, crollinson@gmail.com
# Date Created: 7 May 2015
# Last Modified: 7 May 2015 
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
# Define constants
# ----------------------------------------
sec2yr <- 1*60*60*24*365
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
setwd("~/Desktop/PalEON CR/paleon_mip_site")
inputs <- "phase1a_output_variables"
fig.dir <- "~/Desktop/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures"
# ----------------------------------------

# ----------------------------------------
# Load Data Sets
# ----------------------------------------
# Ecosystem Model Outputs
ecosys <- read.csv(file.path(inputs, "MIP_Data_Ann_2015.csv"))
ecosys$Model.Order <- recode(ecosys$Model, "'clm.bgc'='01'; 'clm.cn'='02'; 'ed2'='03'; 'ed2.lu'='04';  'jules.stat'='05'; 'jules.triffid'='06'; 'linkages'='07'; 'lpj.guess'='08'; 'lpj.wsl'='09'; 'sibcasa'='10'")
levels(ecosys$Model.Order) <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
summary(ecosys)

# CO2 Record
nc.co2 <- nc_open("~/Desktop/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)

# Merging CO2 into Model Outputs
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# Colors used for graphing
model.colors <- read.csv("~/Desktop/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors $Model.Order <- recode(model.colors$Model, "'CLM4.5-BGC'='01'; 'CLM4.5-CN'='02'; 'ED2'='03'; 'ED2-LU'='04';  'JULES-STATIC'='05'; 'JULES-TRIFFID'='06'; 'LINKAGES'='07'; 'LPJ-GUESS'='08'; 'LPJ-WSL'='09'; 'SiBCASA'='10'")
levels(model.colors$Model.Order)[1:10] <- c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")
model.colors

model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors
# ----------------------------------------

# ----------------------------------------
# Calculate Deviations from desired reference point
# Reference Point: 0850-0869 (spinup climate)
# ----------------------------------------
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp")
vars.climate <- c("Temp", "Precip", "CO2")

# vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

ref.window <- 850:869

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
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
vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp", "Temp", "Precip", "CO2")
vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		# # -----------------------
		# # Differencing
		# # -----------------------
		# ## Non-standardized
		# for(v in vars[2]){
			# temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			# ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".d1")] <- c(NA, diff(temp, lag=1,))
		# }

		# # -----------------------


		# -----------------------
		# Decadal Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="right", fill=NA)
		}
		# -----------------------

		# -----------------------
		# Centennial Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="right", fill=NA)
		}
		# -----------------------
	}
}
summary(ecosys)

# -----------------------
# Some exploratory Graphing
# -----------------------
colors.use <- 
ggplot(data=ecosys[,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=AGB, color=Model.Order), size=1, alpha=0.6) +
	geom_line(aes(x=Year, y=AGB.100, color=Model.Order), size=1.5) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model %in% unique(ecosys$Model.Order),"color"])) +
	theme_bw()
# -----------------------

# -----------------------
# Subsetting individual models & sites for model prototyping
# -----------------------
ed2.pha <- ecosys[ecosys$Model=="ed2" & ecosys$Site=="PHA",]
summary(ed2.pha)

lpj.g.pha <- ecosys[ecosys$Model=="lpj.guess" & ecosys$Site=="PHA",]
summary(lpj.g.pha)

lpj.g <- ecosys[ecosys$Model=="lpj.guess",]
summary(lpj.g)
# -----------------------
# ----------------------------------------



# ----------------------------------------
# Model approach: AGB ~ 3 non-interactive temporal smoothers: AGB, Temp, Precip
# ----------------------------------------
library(mgcv)

# summary(lpj.g.pha)

# gam.lpj.g.pha <- gamm(AGB ~ s(Temp, k=3) + s(Precip, k=3) + s(CO2, k=3), data=lpj.g.pha, correlation=corCAR1(form=~Year))
# summary(gam.lpj.g.pha$gam)
# summary(gam.lpj.g.pha$lme)
# gam.lpj.g.pha0 <- gam(AGB ~ s(Temp, k=3) + s(Precip, k=3) + s(CO2, k=3), data=lpj.g.pha)
# summary(gam.lpj.g.pha0)

# acf(gam.lpj.g.pha0$resid)
# acf(gam.lpj.g.pha$gam$resid)
# ------------------------------------------------
# All Sites: (for 1 site, see model selection script)
# ------------------------------------------------
source('~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/predict.gamm.model.R', chdir = TRUE)
# gam.lpj.g <- gamm(AGB ~ s(Temp, by=Site, k=3) + s(Precip, by=Site, k=3) + s(CO2, by=Site, k=3) + Site - 1, data=lpj.g, correlation=corCAR1(form=~1|Site))
outdir="~/Desktop/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling/R/"

# ------------------------
# LPJ-GUESS
# ------------------------
gam.lpj.guess <- model.gam(data=ecosys, model="lpj.guess", response="AGB", k=4, outdir)
summary(gam.lpj.guess)
summary(gam.lpj.guess$data)
summary(gam.lpj.guess$gam)
summary(gam.lpj.guess$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.lpj.guess$gam)
dev.off()

lm.lpj.g <- lm(fit.gam ~ response, data=gam.lpj.guess$data)
summary(lm.lpj.g)

ggplot(data=gam.lpj.guess$data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$weights) + facet_wrap(~Site) +
	geom_line(data=gam.lpj.guess$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$data[gam.lpj.guess$data$Site=="PHA" & gam.lpj.guess$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-GUESS_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.guess$data[gam.lpj.guess$data$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"temp"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"co2"]),
			  			abs(gam.lpj.guess$weights[gam.lpj.guess$weights$Site=="PHA" & gam.lpj.guess$data$Year>=1350 & gam.lpj.guess$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-GUESS, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# ------------------------

# ------------------------
# LPJ-WSL
# ------------------------
gam.lpj.wsl <- model.gam(data=ecosys, model="lpj.wsl", response="AGB", k=4, outdir)
summary(gam.lpj.wsl)
summary(gam.lpj.wsl$data)
summary(gam.lpj.wsl$gam)
summary(gam.lpj.wsl$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.lpj.wsl$gam)
dev.off()

lm.lpj.wsl <- lm(fit.gam ~ response, data=gam.lpj.wsl$data)
summary(lm.lpj.wsl)

ggplot(data=gam.lpj.wsl$data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$weights) + facet_wrap(~Site) +
	geom_line(data=gam.lpj.wsl$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$data[gam.lpj.wsl$data$Site=="PHA" & gam.lpj.wsl$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LPJ-WSL_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.lpj.wsl$data[gam.lpj.wsl$data$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"temp"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"co2"]),
			  			abs(gam.lpj.wsl$weights[gam.lpj.wsl$weights$Site=="PHA" & gam.lpj.wsl$data$Year>=1350 & gam.lpj.wsl$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LPJ-WSL, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------


# ------------------------
# LINKAGES
# ------------------------
# note: linkages didn't return temp or precip, but was updated, so we're going to assume it's the same as LPJ-Guess
data.linkages <- ecosys[ecosys$Model=="linkages",]
data.linkages$Temp <- ecosys[ecosys$Model=="lpj.guess", "Temp"]
data.linkages$Precip <- ecosys[ecosys$Model=="lpj.guess", "Precip"]
summary(data.linkages)

gam.linkages <- model.gam(data=data.linkages, model="linkages", response="AGB", k=3, outdir)
summary(gam.linkages)
summary(gam.linkages$data)
summary(gam.linkages$gam)
summary(gam.linkages$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.linkages$gam)
dev.off()

lm.linkages <- lm(fit.gam ~ response, data=gam.linkages$data)
summary(lm.linkages)

ggplot(data= gam.linkages $data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, Linkages, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$weights) + facet_wrap(~Site) +
	geom_line(data=gam.linkages$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$data[gam.linkages$data$Site=="PHA" & gam.linkages$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_LINKAGES_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.linkages$data[gam.linkages$data$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"temp"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"co2"]),
			  			abs(gam.linkages$weights[gam.linkages$weights$Site=="PHA" & gam.linkages$data$Year>=1350 & gam.linkages$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, LINKAGES, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------


# ------------------------
# ED2
# ------------------------
gam.ed2 <- model.gam(data=ecosys, model="ed2", response="AGB", k=3, outdir)
summary(gam.ed2)
summary(gam.ed2$data)
summary(gam.ed2$gam)
summary(gam.ed2$weights)

pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_Splines.pdf"), width=11, height=8.5)
par(mfrow=c(3,6), mar=c(5,5,0.5, 0.5))
plot(gam.ed2$gam)
dev.off()

lm.ed2 <- lm(fit.gam ~ response, data=gam.ed2$data)
summary(lm.ed2)

ggplot(data= gam.ed2 $data) + # facet_wrap(~Site) +
	geom_point(aes(x=response, y=fit.gam, color=Site), size=2) +
	geom_abline(intercept=0, slope=1, color="gray50", linetype="dashed") +
	labs(x="Observed", y="Modeled", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))


pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_0850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$weights) + facet_wrap(~Site) +
	geom_line(data=gam.ed2$data, aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHO",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHO","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PUN",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PUN","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PBL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PBL","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PDL",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PDL","precip"])), size=4) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PMB",], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PMB","precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()

# Just post-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_1850-2010.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$data[gam.ed2$data$Site=="PHA" & gam.ed2$data$Year>=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$weights$Year>= 1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()


# 600 years pre-1850 at PHA
pdf(file.path(fig.dir, "Non-StationaryDrivers_ED2_AGB_Annual_1350-1850.pdf"), width=11, height=8.5)
ggplot(data=gam.ed2$data[gam.ed2$data$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,]) + facet_wrap(~Site) +
	geom_line(aes(x=Year, y=response), color="gray50", size=2) +
	geom_line(data=gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,], aes(x=Year, y=fit), 
			  color=rgb(abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"temp"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"co2"]),
			  			abs(gam.ed2$weights[gam.ed2$weights$Site=="PHA" & gam.ed2$data$Year>=1350 & gam.ed2$data$Year<=1850,"precip"])), size=4) +
	labs(x="Year", y="AGB kg/m2", title="Non-Stationary Drivers of AGB, ED2, Annual") +
	theme_bw() + theme(axis.text.x=element_text(angle=0, color="black", size=rel(1.25)), axis.text.y=element_text(color="black", size=rel(1.25)), axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1), plot.title=element_text(face="bold", size=rel(2)))
dev.off()
# ------------------------