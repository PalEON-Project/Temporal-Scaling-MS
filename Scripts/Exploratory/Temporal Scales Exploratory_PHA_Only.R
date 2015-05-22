# Exploratory analyses for the Carbon Paper
library(zoo)
library(ncdf4)
library(ggplot2)
fig.path <- "Analyses/Temporal Scaling/Figures"

setwd("~/Dropbox/PalEON CR/paleon_mip_site")
outputs <- "phase1a_output_variables"
large.axes <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=24, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

# -----------------------------------------------------------
# Getting Yearly Data
# -----------------------------------------------------------



# # # -----------------------------------------------------------
# # Going from Monthly to Yearly Climate variables
# # -----------------------------------------------------------
# met.month <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Monthly.csv")
# met.month <- met.month[met.month$Site=="PHA",]
# summary(met.month)

# vars.met <- names(met.month[,5:ncol(met.month)])

# # Creating an object of yearly means
# met.year <- aggregate(met.month[,vars.met], by=list(met.month$Site, met.month$Year), FUN=mean)
# names(met.year) <- c("Site", "Year", vars.met)
# summary(met.year)

# # This should be modified based on a priori list of temporally varying env factors that drive ecosystem function
# # This step should include the facotrs for the Carbon Paper and Resilience paper
# for(s in unique(met.year$Site)){
	# for(y in min(met.year$Year):max(met.year$Year)){
		# met.year[met.year$Site==s & met.year$Year==y ,"precipf.JJA"] <- mean(met.month[met.month$Site==s & met.month$Year==y & met.month$Month>=6 & met.month$Month<=8, "precipf"], na.rm=T)
		# met.year[met.year$Site==s & met.year$Year==y ,"tair.JJA"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=6 & met.month $Month<=8, "tair"])
		# met.year[met.year$Site==s & met.year$Year==y ,"tair.MAM"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=3 & met.month $Month<=5, "tair"])
		# met.year[met.year$Site==s & met.year$Year==y ,"tair.SON"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=10 & met.month $Month<=11, "tair"])
	# }
# }
# summary(met.year)


# write.csv(met.year, "phase1a_output_variables/PalEON_MIP_Drivers_Year_PHA.csv", row.names=F)
# # -----------------------------------------------------------


# -----------------------------------------------------------
# Formatting Met & Ecosystem Data
# -----------------------------------------------------------
# met.year <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Year_PHA.csv")
# summary(met.year)
# met.col <- ncol(met.year)
library(ncdf4)
nc.co2 <- nc_open("~/Dropbox/PalEON CR/paleon_mip_site/env_drivers/phase1a_env_drivers_v4/paleon_co2/paleon_annual_co2.nc")
co2.ann <- data.frame(CO2=ncvar_get(nc.co2, "co2"), Year=850:2010)
nc_close(nc.co2)
# plot(CO2~Year, data=co2.ann, type="l")

ecosys <- read.csv(file.path(outputs, "MIP_Data_Ann_2015.csv"))
ecosys <- merge(ecosys, co2.ann)
summary(ecosys)

# write.csv(file.path(outputs, "MIP_Data_Ann_2015.csv"), row.names=F)

ecosys <- ecosys[ecosys$Site=="PHA",]
summary(ecosys)

# Calculating the percent deviations in the ecosystem variables
for(s in unique(ecosys$Site)){
	# ---------------------------
	# Driver Deviations
	# ---------------------------
	# for(v in names(met.year[,3:met.col])){
		# tmp <- mean(met.year[met.year$Site==s, v], na.rm=T)
		# met.year[met.year$Site==s, paste(v, "dev", sep=".")] <- met.year[met.year$Site==s, v] - tmp
	# }

	# ---------------------------
	
	# ---------------------------
	# Ecosystem Deviations
	# Update: centering on 1990-2010
	# ---------------------------
	for(m in unique(ecosys$Model)){
		mean.subset <- ecosys[ecosys$Site==s & ecosys$Model==m & ecosys$Year>= 1990 & ecosys$Year<=2010, ]
		x.gpp <- mean(mean.subset$GPP, na.rm=T)
		x.agb <- mean(mean.subset$AGB, na.rm=T)
		x.lai <- mean(mean.subset$LAI, na.rm=T)
		x.npp <- mean(mean.subset$NPP, na.rm=T)
		x.nee <- mean(mean.subset$NEE, na.rm=T)
		x.ra <- mean(mean.subset$AutoResp, na.rm=T)
		x.rh <- mean(mean.subset$HeteroResp, na.rm=T)
		x.soilcarb <- mean(mean.subset$SoilCarb, na.rm=T)
		x.soilmoist <- mean(mean.subset$SoilMoist, na.rm=T)
		x.evap <- mean(mean.subset$Evap, na.rm=T)
		x.transp <- mean(mean.subset$Transp, na.rm=T)
		x.temp <- mean(mean.subset$Temp, na.rm=T)-273.15
		x.precip <- mean(mean.subset$Precip, na.rm=T)
		x.co2 <- mean(mean.subset$CO2, na.rm=T)

		ecosys[ecosys$Site==s & ecosys$Model==m, "GPP.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "GPP"] - x.gpp)/x.gpp
		ecosys[ecosys$Site==s & ecosys$Model==m, "AGB.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "AGB"] - x.agb)/x.agb
		ecosys[ecosys$Site==s & ecosys$Model==m, "LAI.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "LAI"] - x.lai)/x.lai
		ecosys[ecosys$Site==s & ecosys$Model==m, "NPP.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "NPP"] - x.npp)/x.npp
		ecosys[ecosys$Site==s & ecosys$Model==m, "NEE.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "NEE"] - x.nee)/x.nee
		ecosys[ecosys$Site==s & ecosys$Model==m, "AutoResp.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "AutoResp"] - x.ra)/x.ra
		ecosys[ecosys$Site==s & ecosys$Model==m, "HeteroResp.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "HeteroResp"] - x.rh)/x.rh
		ecosys[ecosys$Site==s & ecosys$Model==m, "SoilCarb.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "SoilCarb"] - x.soilcarb)/x.soilcarb
		ecosys[ecosys$Site==s & ecosys$Model==m, "SoilMoist.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "SoilMoist"] - x.soilmoist)/x.soilmoist
		ecosys[ecosys$Site==s & ecosys$Model==m, "Evap.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Evap"] - x.evap)/x.evap
		ecosys[ecosys$Site==s & ecosys$Model==m, "Transp.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Transp"] - x.transp)/x.transp
		ecosys[ecosys$Site==s & ecosys$Model==m, "Temp.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Temp"]-273.15 - x.temp)/x.temp
		ecosys[ecosys$Site==s & ecosys$Model==m, "Precip.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Precip"] - x.precip)/x.precip
		ecosys[ecosys$Site==s & ecosys$Model==m, "CO2.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "CO2"] - x.co2)/x.co2

		ecosys[ecosys$Site==s & ecosys$Model==m, "Temp.abs.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Temp"]-273.15 - x.temp)
		ecosys[ecosys$Site==s & ecosys$Model==m, "Precip.abs.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "Precip"] - x.precip)
		ecosys[ecosys$Site==s & ecosys$Model==m, "CO2.abs.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "CO2"] - x.co2)

	}
	# ---------------------------
}
# summary(met.year)
summary(ecosys)

sites <- read.csv("env_drivers/phase1a_env_drivers_v4/PalEON_Phase1a_sites.csv")
row.names(sites) <- sites$code
sites <- sites[,3:ncol(sites)]
summary(sites)

sites2 <- data.frame(t(sites))
sites2$Site <- as.factor(row.names(sites2))
summary(sites2)
dim(sites2)
write.csv(sites2, "env_drivers/phase1a_env_drivers_v4/PalEON_Phase1a_sites2.csv", row.names=F)

sites <- read.csv("env_drivers/phase1a_env_drivers_v4/PalEON_Phase1a_sites2.csv")
summary(sites)

# -----------------------------------
# Reshaping for ordination
environ <- merge(met.year[,c("Site", "Year", "CO2")], sites2, all.x=T, all.y=F)
summary(environ)

library(reshape2)
library(reshape)
summary(ecosys)

ecosys$Year <- as.factor(ecosys$Year)
ecosys2 <- melt(ecosys)
summary(ecosys2)
npp.models <- recast(ecosys[,c("Year", "Site", "Model", "NPP.dev")], Site+Year~Model)
summary(npp.models)
# -----------------------------------------------------------


# -----------------------------------------------------------
# Model correlation across temporal scales
# -----------------------------------------------------------
summary(ecosys)
summary(ecosys)
ecosys$Year <- as.numeric(paste(ecosys$Year))
summary(ecosys)

## Don't need this now, just use rollmean in zoo
# smooth <- function(x, y, index, years, window, na.rm=T){
	# for(i in (min(years)+window):(max(years)-window)){
		# mean(x[x[,index]>=(i-window) & x[,index]<=(i+window),y], na.rm=T)
	# }
# }

vars <- c("GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp", "Temp", "Precip", "CO2")
vars.dev <- c(paste0(vars[1:(length(vars)-3)], ".dev"), "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")

for(s in unique(ecosys$Site)){
	for(m in unique(ecosys$Model)){
		# -----------------------
		# Decadal Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]

			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="center", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".10")] <- rollmean(temp, k=10, align="center", fill=NA)
		}
		# -----------------------

		# -----------------------
		# Centennial Smoothing
		# -----------------------
		## Non-standardized
		for(v in vars){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="center", fill=NA)
		}

		## Non-standardized
		for(v in vars.dev){
			temp <- ecosys[ecosys$Model==m & ecosys$Site==s, v]
			ecosys[ecosys$Model==m & ecosys$Site==s, paste0(v, ".100")] <- rollmean(temp, k=100, align="center", fill=NA)
		}
		# -----------------------
	}
}
summary(ecosys)
#summary(met.year)
write.csv(ecosys, "phase1a_output_variables/PalEON_MIP_Yearly_2_PHA.csv", row.names=F)
# write.csv(met.year, "phase1a_output_variables/PalEON_MIP_Drivers_Year_2_PHA.csv", row.names=F)

library(ggplot2)
ggplot(data=ecosys[!ecosys$Model=="jules.stat",]) +
	geom_line(aes(x=Year, y=AGB, color=Model, size=Updated), size=1, alpha=0.5) +
	geom_line(aes(x=Year, y=AGB.100, color=Model, size=Updated)) +
	scale_size_manual(values=c(0.5,2)) +
	theme_bw() + guides(Updated=F)
# --------------------------------------------------------











# --------------------------------------------------------
# Some exploratory Graphing
# --------------------------------------------------------
library(ggplot2)
library(dplR)
library(grid)

# Model Data sets
ecosys <- read.csv("phase1a_output_variables/PalEON_MIP_Yearly_2_PHA.csv")
library(car)
ecosys$Model.Order <- recode(ecosys$Model, "'ed2'='1'; 'ed2.lu'='2'; 'clm45'='3'; 'lpj.wsl'='4'; 'lpj.guess'='5'; 'jules.stat'='6'; 'SiB'='7'; 'linkages'='8'")
levels(ecosys$Model.Order) <- c("ED2", "ED2-LU", "CLM4.5", "LPJ-WSL", "LPJ-GUESS", "JULES", "SiBCASA", "LINKAGES")
summary(ecosys)

# Using only what's been updated
ecosys <- ecosys[ecosys$Updated=="Yes",]
summary(ecosys)


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_Temperature_0850-2010.pdf", width=10, height=3)
ggplot(data=ecosys[ecosys$Model=="ed2",]) + large.axes +
	geom_line(aes(x=Year, y=Temp-273.15), size=0.5, color="red", alpha=0.5) +
	geom_line(aes(x=Year, y=Temp.100-273.15), size=2, color="red") +		
	scale_y_continuous(name=expression(bold(paste("Temp ("^"o","C)")))) +
	scale_x_continuous(name="Year")  +
	theme(axis.text.x=element_text(angle=0, color="black", size=18), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=28, vjust=-1),  axis.title.y=element_text(face="bold", size=20, vjust=2))
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_Precipitation_0850-2010.pdf", width=10, height=3)
ggplot(data=ecosys[ecosys$Model=="ed2",]) + large.axes +
	geom_line(aes(x=Year, y=Precip), size=0.5, color="blue2", alpha=0.5) +
	geom_line(aes(x=Year, y=Precip.100), size=2, color="blue2") +		scale_y_continuous(name=expression(bold(paste("Precip (m yr"^"-1",")"))), breaks=c(0.8,1,1.2, 1.4, 1.6)) +
	scale_x_continuous(name="Year")  +
	theme(axis.text.x=element_text(angle=0, color="black", size=18), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=28, vjust=-1),  axis.title.y=element_text(face="bold", size=20, vjust=2))
dev.off()

# Model color sc
model.colors <- read.csv("~/Dropbox/PalEON CR/PalEON_MIP_Site/Model.Colors.csv")
model.colors$Model.Order <- recode(model.colors$Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'; 'LINKAGES'='8'")
levels(model.colors$Model.Order)[1:8] <- c("ED2", "ED2-LU", "CLM4.5", "LPJ-WSL", "LPJ-GUESS", "JULES", "SiBCASA", "LINKAGES")
model.colors <- model.colors[order(model.colors$Model.Order),]
model.colors

# met.year <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Year_2_PHA.csv")
# summary(met.year)

# Validation Data sets
data.flux <- read.csv("Validation Data/Flux_All_Year.csv")
data.nbcd <- read.csv("Validation Data/NBCD_MIP_Sites.csv")
data.setveg <- read.csv("Validation Data/SetVegBiomass_MIP_Sites.csv")
quru.chron <- read.csv("Validation Data/Lyford_QURU_Chronology.csv")
quru.rwl <- read.rwl("~/Dropbox/PalEON CR/PalEON_MIP_Site/Validation Data/TreeRings/Lyford_Data_13m/RW/Combined/Lyford_QURU.rw")
data.agb <- read.csv("~/Dropbox/PalEON CR/paleon_mip_site/Validation Data/HarvardForest_LTER/HF_C_Synthesis_ConfInt.csv")

# Detrending QURU to get CIs
quru.detrend <- detrend(quru.rwl, method="Spline", make.plot=F)
summary(quru.detrend)
row.names(quru.detrend)

quru.detrend.10 <- quru.detrend
for(j in 1:ncol(quru.detrend.10)){
	temp <- quru.detrend.10[!is.na(quru.detrend[,j]),j]
	quru.detrend.10[!is.na(quru.detrend[,j]), j] <- rollmean(temp, k=10, align="center", fill=NA)
}
summary(quru.detrend.10)



quru.stack <- stack(quru.detrend)
names(quru.stack) <- c("RWI", "TreeID")
quru.stack$RWI <- quru.stack$RWI-1
quru.stack$Year <- as.factor(row.names(quru.detrend))
quru.stack$RWI.Smooth <- stack(quru.detrend.10)[,1]-1
summary(quru.stack)



quru.ci <- data.frame(Year=unique(quru.stack$Year))

lm.quru <- lm(RWI ~ Year, data=quru.stack)

quru.predict.ann <- predict(lm.quru, interval="confidence", newdata=quru.ci)
quru.predict.ann <- cbind(quru.ci, quru.predict.ann)
summary(quru.predict.ann)
ggplot(data=quru.predict.ann) +
	geom_ribbon(aes(x=as.numeric(Year), ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=as.numeric(Year), y=fit))


lm.quru.10 <- lm(RWI.Smooth ~ Year, data=quru.stack)
quru.ci2 <- data.frame(Year=quru.ci[which(quru.ci$Year %in% years.10),])

years.10 <- unique(quru.stack[!is.na(quru.stack$RWI.Smooth), "Year"])
quru.predict.dec <- predict(lm.quru.10, interval="confidence", newdata=quru.ci2)
quru.predict.dec <- cbind(quru.ci2, quru.predict.dec)
summary(quru.predict.dec)
ggplot(data=quru.predict.dec) +
	geom_ribbon(aes(x=as.numeric(Year), ymin=lwr, ymax=upr), alpha=0.5) +
	geom_line(aes(x=as.numeric(Year), y=fit))

ggplot(data=data.agb) + 
	geom_ribbon(aes(x=Year, ymin=CI.Lwr, ymax=CI.Upr), alpha=0.5) +
	geom_line(aes(x=Year, y=Mean))

summary(data.flux)
summary(data.flux[data.flux$Site=="HarvardForest1",])
summary(data.nbcd)
summary(data.setveg)
summary(quru.chron)

large.axes <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=16), axis.text.y=element_text(color="black", size=12), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=24, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

# model.colors <- c("black", "gray50", "cyan4", "hotpink3", "green3", "orange3", "steelblack3")
# model.colors2 <- c("black", "cyan4", "green3", "orange3", "steelblack3")

# SetVeg Range: 0 - 217.42
setveg <- data.frame(Year = 1781:1851, Min=0, Max=217.52/2*.1)
summary(setveg)
modernbm <- data.frame(Year=1990:2010, Min=min(data.nbcd$nbcd.mean, na.rm=T)/2*.1, Max=max(data.nbcd$nbcd.mean, na.rm=T)/2*.1)
sec2yr <- 1*60*60*24*365

pseduo.pollen1 <- data.frame(Year1=seq(1800, 1850, length.out=50), Year2=seq(900, 1200, length.out=50), Year3=seq(1400, 1500, length.out=50), Min=0, Max=15)

# ------------------------------
# Validation data sets
# ------------------------------
# Empirical Data Only
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_EmpiricalOnly_1800-2010.pdf", height=6, width=8)
ggplot() +
	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3, fill="darkgreen")+
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="darkgreen") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="darkgreen", size=2) +
	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=5, color="darkgreen") +
	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=2, color="darkgreen") +
	geom_line(data=quru.chron, aes(x=X, y=xxxstd*6), size=1) +
	geom_line(data=data.flux[data.flux$Site=="HarvardForest1" & data.flux$Year>1991,], aes(x=Year, y=-NEE*3+1), color="black", size=2) + 
	scale_y_continuous(breaks=0, name="") + 
	scale_x_continuous(limits=c(1800, 2015), name="") + 
	large.axes 
dev.off()
	
# Empirical Data Only -- Full Time
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_EmpiricalOnly_0850-2010.pdf", height=6, width=8)
ggplot() +
	geom_ribbon(data=pseduo.pollen1, aes(x=Year1, ymin=Min, ymax=Max), alpha=0.3, fill="red3")+
	geom_ribbon(data=pseduo.pollen1, aes(x=Year2, ymin=Min, ymax=Max), alpha=0.3, fill="red3")+
	geom_ribbon(data=pseduo.pollen1, aes(x=Year3, ymin=Min, ymax=Max), alpha=0.3, fill="red3")+
	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.6, fill="darkgreen")+
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="darkgreen") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="darkgreen", size=2) +
	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=5, color="darkgreen") +
	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=2, color="darkgreen") +
	geom_line(data=quru.chron, aes(x=X, y=xxxstd*6), size=1) +
	geom_line(data=data.flux[data.flux$Site=="HarvardForest1" & data.flux$Year>1991,], aes(x=Year, y=-NEE*3+1), color="black", size=2) +
	scale_y_continuous(breaks=-1, name="") + 
	scale_x_continuous(limits=c(850, 2015), name="Year") + 
	large.axes 	
dev.off()

#colors1 <- as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])
# AGB
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_ModelsOnly_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES",]) + large.axes +
#	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
	# geom_ribbon(data=modernbm, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
#	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
#	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="gray50") +
	geom_line(aes(x=Year, y=AGB, color= Model.Order, size=Updated)) +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
	scale_size_manual(values=c(1,3)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=2), alpha=F, size=F) + 
 	labs(color="Models")
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_ModelsOnly_0850-2010_withCentennial.pdf", height=6, width=9)
ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES",]) + large.axes +
	geom_line(aes(x=Year, y=AGB, color= Model.Order), size=1, alpha=0.5) +
	geom_line(aes(x=Year, y=AGB.100, color= Model.Order), size=2, alpha=1) +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
 	theme(legend.position=c(0.30,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=3), alpha=F, size=F) + 
 	labs(color="Models")
dev.off()


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_ModelsOnly_0850-2010_withCentennial_2.pdf", height=6, width=9)
ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES",]) + large.axes +
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="black") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="black", size=3) +
	geom_line(aes(x=Year, y=AGB, color= Model.Order), size=1, alpha=0.5) +
	geom_line(aes(x=Year, y=AGB.100, color= Model.Order), size=2, alpha=1) +
	# geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
	# geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="black") +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
 	theme(legend.position=c(0.30,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=3), alpha=F, size=F) + 
 	labs(color="Models")
dev.off()


# AGB
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_ModelsOnly_1800-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1800 & !ecosys$Model.Order=="JULES",]) + large.axes +
#	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
	# geom_ribbon(data=modernbm, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
#	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
#	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="gray50") +
	geom_line(aes(x=Year, y=AGB, color=Model.Order, size=Updated)) +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
	scale_size_manual(values=c(1,3)) +
 	theme(legend.position=c(0.75,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=2), alpha=F, size=F)+ 
 	labs(color="Models")
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_1800-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1800 & !ecosys$Model.Order=="JULES",]) + large.axes +
	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3, fill="black")+
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="black") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="black", size=3) +

	# geom_ribbon(data=modernbm, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
	geom_line(aes(x=Year, y=AGB, color=Model.Order, size=Updated)) +
	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="black") +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
	scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.85,0.13), legend.text=element_text(size=10), legend.title=element_text(size=12), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)+ 
 	labs(color="Models")
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1900 & !ecosys$Model.Order=="JULES",]) + large.axes +
	# geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3, fill="black")+
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="black") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="black", size=3) +
	# geom_ribbon(data=modernbm, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
	geom_line(aes(x=Year, y=AGB, color=Model.Order, size=Updated)) +
	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="black") +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
	scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.52,0.85), legend.text=element_text(size=10), legend.title=element_text(size=12), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + 
 	guides(col=guide_legend(ncol=2), alpha=F, size=F) +
 	labs(color="Models")
dev.off()


# AGB
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES",]) + large.axes +
	# geom_ribbon(data=modernbm, aes(x=Year, ymin=Min, ymax=Max), alpha=0.3)+
	geom_line(aes(x=Year, y=AGB, color=Model.Order, size=Updated)) +
	geom_ribbon(data=setveg, aes(x=Year, ymin=Min, ymax=Max), alpha=0.5, fill="black")+
	geom_ribbon(data=data.agb, aes(x=Year, ymin=CI.Lwr*.1/2, ymax=CI.Upr*.1/2), alpha=0.5, fill="black") +
	geom_line(data=data.agb, aes(x=Year, y=Mean*.1/2), color="black", size=2) +
	geom_errorbar(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, ymin=(nbcd.mean-nbcd.sd)/2*.1, ymax=(nbcd.mean+nbcd.sd)/2*.1), size=3, color="black", alpha=0.5) +
	geom_point(data=data.nbcd[data.nbcd$Site=="PHA",], aes(x=2000, y=nbcd.mean/2*.1), size=7, color="black") +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd)) +
	ylab(expression(bold(paste("Aboveground Biomass (kgC m"^"-2",")")))) +
	scale_y_continuous(limits=c(0, max(ecosys$AGB, na.rm=T)))+
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"])) +
	scale_size_manual(values=c(1,3)) +
 	theme(legend.position=c(0.25,0.13), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=F, alpha=F, size=F) +
 	labs(color="Models")
dev.off()

# ----------------------------------
# GPP
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_GPP_ModelsOnly_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	# geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=GPP), color="gray50", size=3) +
	geom_line(aes(x=Year, y=GPP*sec2yr, color= Model.Order, size=Updated)) +
	ylab(expression(bold(paste("GPP (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_y_continuous(limits=c(0, max(ecosys$GPP*sec2yr, na.rm=T)))+
	labs(color="Models") +
	scale_size_manual(values=c(0.25,0.7)) +
 	theme(legend.position=c(0.35,0.15), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_GPP_ModelsOnly_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1900,]) + large.axes +
	# geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=GPP), color="gray50", size=3) +
	geom_line(aes(x=Year, y=GPP*sec2yr, color= Model.Order, size=Updated)) +
	ylab(expression(bold(paste("GPP (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_y_continuous(limits=c(0, max(ecosys$GPP*sec2yr, na.rm=T)))+
	labs(color="Models") +
 	theme(legend.position=c(0.85,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_GPP_Validation_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1900,]) + large.axes +
	geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=GPP), color="black", size=2) +
	geom_line(aes(x=Year, y=GPP*sec2yr, color= Model.Order, size=Updated)) +
	# geom_line(data=quru.chron, aes(x=X, y=xxxstd-0.5), color="gray50", size=3) +
	ylab(expression(bold(paste("GPP (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_y_continuous(limits=c(0, max(ecosys$GPP*sec2yr, na.rm=T)))+
	labs(color="Models") +
	scale_size_manual(values=c(0.8,2)) +
 	theme(legend.position=c(0.85,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()




# ----------------------------------
# NPP
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_ModelsOnly_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	# geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=GPP), color="gray50", size=3) +
	geom_line(aes(x=Year, y=NPP*sec2yr, color= Model.Order), size=1) +
	ylab(expression(bold(paste("NPP (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_y_continuous(limits=c(0, max(ecosys$NPP*sec2yr, na.rm=T)))+
	labs(color="Models") +
#	scale_size_manual(values=c(0.75,2)) +
 	theme(legend.position=c(0.25,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + 
 	guides(col=F, alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_ModelsOnly_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	# geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=GPP), color="gray50", size=3) +
	geom_line(aes(x=Year, y=NPP*sec2yr, color= Model.Order), size=2) +
	ylab(expression(bold(paste("NPP (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_y_continuous(limits=c(0, max(ecosys$NPP*sec2yr, na.rm=T)))+
	labs(color="Models") +
#	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.25,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=F, alpha=F, size=F)
dev.off()

# -------------------
# NEE
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NEE_Validation_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>1900,]) + large.axes +
	geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=-NEE), color="black", size=3) +
	geom_line(aes(x=Year, y=NEE*sec2yr, color= Model.Order), size=2) +
	ylab(expression(bold(paste("NEE (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
#	scale_y_continuous(limits=c(0, max(ecosys$NEE*sec2yr, na.rm=T)))+
	labs(color="Models") +
#	scale_size_manual(values=c(0.8,2)) +
 	theme(legend.position=c(0.35,0.1), legend.text=element_text(size=12), legend.title=element_text(size=12), legend.key=element_rect(fill="white"), legend.key.width=unit(1.5, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NEE_Validation_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	geom_line(data=data.flux[substr(data.flux$Site,1,7)=="Harvard" &data.flux$Year>1992,], aes(x=Year, y=-NEE), color="black", size=2) +
	geom_line(aes(x=Year, y=NEE*sec2yr, color= Model.Order), size=1) +
	ylab(expression(bold(paste("NEE (kgC m"^"-2"," yr"^"-1",")")))) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
#	scale_y_continuous(limits=c(0, max(ecosys$NEE*sec2yr, na.rm=T)))+
	labs(color="Models") +
	# scale_size_manual(values=c(0.5,1)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=F, alpha=F, size=F)
dev.off()
# ------------------------------


# --------------------------------------------------------
# Annual NPP Deviation, PHA
# --------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color= Model.Order), size=0.5) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	# scale_size_manual(values=c(0.25,0.75)) +
 	theme(legend.position=c(0.35,0.10), legend.text=element_text(size=14), legend.title=element_text(size=16), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_1800-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1800,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color= Model.Order, size=Updated)) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.35,0.10), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color= Model.Order, size=Updated)) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.35,0.10), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1.5, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

quru.predict.ann$Year <- as.numeric(paste(quru.predict.ann$Year))
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_1900-2010_QURU_Chron.pdf", height=6, width=8)
quru.predict.ann$Year <- as.numeric(paste(quru.predict.ann$Year))
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color=Model.Order, size=Updated), alpha=0.5) +
	# geom_ribbon(aes(x=1900:2010, ymin=min(ecosys$NPP.dev, na.rm=T), ymax=max(ecosys$NPP.dev, na.rm=T)), alpha=0.3) +
	geom_ribbon(data=quru.predict.ann[quru.predict.ann$Year>1900  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), ymin=lwr, ymax=upr), alpha=0.5, fill="black") +
	geom_line(data= quru.predict.ann[quru.predict.ann$Year>1900  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), y=fit), color="black", size=2) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.35,0.10), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_1900-2010_QURU_Chron2.pdf", height=6, width=8)
quru.predict.ann$Year <- as.numeric(paste(quru.predict.ann$Year))
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color=Model.Order, size=Updated), alpha=0.5) +
	geom_ribbon(aes(x=1900:2010, ymin=min(ecosys$NPP.dev, na.rm=T), ymax=max(ecosys$NPP.dev, na.rm=T)), alpha=0.3) +
	geom_ribbon(data=quru.predict.ann[quru.predict.ann$Year>1900  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), ymin=lwr, ymax=upr), alpha=0.5, fill="black") +
	geom_line(data= quru.predict.ann[quru.predict.ann$Year>1900  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), y=fit), color="black", size=1.5) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.3,0.15), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Annual_1990-2010_QURU_Chron.pdf", height=6, width=12)
quru.predict.ann$Year <- as.numeric(paste(quru.predict.ann$Year))
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1985,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev, color=Model.Order), size=3, alpha=0.6) +
	geom_ribbon(data=quru.predict.ann[quru.predict.ann$Year>1985  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), ymin=lwr, ymax=upr), alpha=0.5, fill="black") +
	geom_line(data= quru.predict.ann[quru.predict.ann$Year>1985  & quru.predict.ann$Year<=2010,], aes(x=as.numeric(paste(Year)), y=fit), color="black", size=3) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
#	scale_size_manual(values=c(1,2)) +
 	theme(legend.position=c(0.8,0.10), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F) +	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank()) +
	theme(axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=0),  axis.title.y=element_text(face="bold", size=24, vjust=1)) 
dev.off()


# --------------------------------------------------------
# Decadal NPP Deviation, PHA
# --------------------------------------------------------
# for(y in (6):(nrow(quru.chron)-5)){
	# quru.chron[y, "Smooth.10"] <- mean(quru.chron[(y-5):(y+5), "xxxstd"] )
	# }
# summary(quru.chron)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Decadal_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=AGB.dev.10, color=Model.Order), size=3) +
	ylab("Relative AGB") +
	xlab("Year") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
#	scale_size_manual(values=c(3,4)) +
 	theme(legend.position=c(0.8,0.9), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)
dev.off()


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Decadal_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_line(aes(x=Year, y=AGB.dev.10, color=Model.Order), size=3) +
	ylab("Relative AGB") +
	xlab("Year") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	# scale_size_manual(values=c(3,4)) +
 	theme(legend.position=c(0.8,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Decadal_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev.10, color=Model.Order), size=2) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	# scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=F, alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Decadal_1900-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev.10, color= Model.Order), size=3) +
	ylab("Relative NPP") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	# scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.8,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)
dev.off()


quru.predict.dec$Year <- as.numeric(paste(quru.predict.dec$Year))
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Decadal_1900-2010_QURU_Chron.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & ecosys$Year>=1900,]) + large.axes +
	geom_ribbon(data=quru.predict.dec[quru.predict.dec$Year>1900  & quru.predict.dec$Year<=2010,], aes(x=as.numeric(paste(Year)), ymin=lwr, ymax=upr), alpha=0.3, fill="black") +
	geom_line(aes(x=Year, y=NPP.dev.10, color= Model.Order, size=Updated), alpha=0.5) +
	geom_ribbon(data=quru.predict.dec[quru.predict.dec$Year>1900  & quru.predict.dec$Year<=2010,], aes(x=as.numeric(paste(Year)), ymin=lwr, ymax=upr), alpha=0.3, fill="black") +
	geom_line(data= quru.predict.dec[quru.predict.dec$Year>1900  & quru.predict.dec$Year<=2010,], aes(x=as.numeric(paste(Year)), y=fit), color="black", size=3) +
	ylab("Relative NPP") +
	xlab("Year") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.75,0.1), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(1, "line")) + guides(col=guide_legend(ncol=3), alpha=F, size=F)
dev.off()



# --------------------------------------------------------

# --------------------------------------------------------
# Centennial NPP Deviation, PHA
# --------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Centennial_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" ,]) + large.axes +
	geom_line(aes(x=Year, y=AGB.dev.100, color=Model.Order, size=Updated)) +
	ylab("Relative AGB") +
	xlab("Year") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(3,4)) +
 	theme(legend.position=c(0.8,0.9), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)
dev.off()


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Centennial_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" ,]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev.100, color=Model.Order), size=2) +
	ylab("Relative NPP") +
	xlab("Year") +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
#	scale_size_manual(values=c(3,4)) +
 	theme(legend.position=c(0.82,0.9), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2), alpha=F, size=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Centennial_0850-2010_with_Temp.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=NPP.dev.100, color=Model.Order, size=Updated)) +
	geom_line(aes(x=Year, y=Temp.abs.dev.100), color="red2", size=2) +
	ylab("Relative NPP") +
	xlab("Year") +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Models") +
	scale_size_manual(values=c(2,3)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=F, alpha=F, size=F)
dev.off()

# --------------------------------------------------------



# -----------------------------------------------------------
# -----------------------------------------------------------
# Correlation between 2 variables
# -----------------------------------------------------------
# -----------------------------------------------------------

#agb.df <- data.frame(AGB=rep(seq(0, max(ecosys$AGB, na.rm=T), length.out=250), 5), Model.Order=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("CLM4.5", 250), rep("SiBCASA", 250)))
agb.df <- data.frame(AGB=rep(seq(0, max(ecosys$AGB, na.rm=T), length.out=250), 5), Model.Order=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)))
summary(agb.df)
#model.colors2b <- c("cyan4", "black", "green3", "orange3", "steelblack3")

lm.lai <- lm(LAI ~ AGB*Model.Order -1, data=ecosys[,])
summary(lm.lai)

lai.ci <- predict(lm.lai, interval="confidence", newdata=agb.df)
summary(lai.ci)
lai.ci <- cbind(agb.df, lai.ci)
summary(lai.ci)

ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES", ]) + large.axes +
	geom_point(aes(x=AGB, y=LAI, color=Model.Order, size=Updated) ) +
	geom_ribbon(data=lai.ci, aes(x=AGB, ymin=lwr, ymax=upr, fill= Model.Order), alpha=0.5) +
	geom_line(data=lai.ci, aes(x=AGB, y=fit, color=Model.Order)) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_fill_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"][c(1,2,4,3)])) +
	labs(color="Models") +
	scale_size_manual(values=c(0.5,1)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=2), alpha=F, size=F, fill=F)

lm.npp <- lm(NPP ~ AGB*Model.Order -1, data=ecosys[!ecosys$Model.Order=="JULES",])
summary(lm.npp)

npp.ci <- predict(lm.npp, interval="confidence", newdata=agb.df)
summary(npp.ci)
npp.ci <- cbind(agb.df, npp.ci)
summary(npp.ci)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_vs_AGB_0850-2010.pdf", height=6, width=8)
ggplot(data=ecosys[ecosys$Site=="PHA" & !ecosys$Model.Order=="JULES", ]) + large.axes +
	geom_point(aes(x=AGB, y=NPP, color=Model.Order, size=Updated) ) +
	geom_ribbon(data=npp.ci, aes(x=AGB, ymin=lwr, ymax=upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=npp.ci, aes(x=AGB, y=fit, color=Model.Order)) +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_fill_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"][c(1,2,4,3)])) +
	labs(color="Models", fill="Models") +
	scale_size_manual(values=c(0.5,1)) +
 	theme(legend.position=c(0.35,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=guide_legend(ncol=2), alpha=F, size=F, fill=F)
dev.off()


dev.scales <- stack(ecosys[,which(substr(names(ecosys), 1, 7)=="NPP.dev")])
names(dev.scales) <- c("NPP.dev", "scale")
levels(dev.scales$scale) <- c("annual", "decadal", "centennial")
dev.scales$Year <- ecosys[,"Year"]
dev.scales$Model <- ecosys[,"Model.Order"]
dev.scales$Site <- as.factor("PHA")
dev.scales$LAI.dev <- stack(ecosys[,which(substr(names(ecosys), 1, 7)=="LAI.dev")])[,1]
dev.scales$AGB.dev <- stack(ecosys[,which(substr(names(ecosys), 1, 7)=="AGB.dev")])[,1]
dev.scales$NEE.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 7)=="NEE.dev")])[,1]
dev.scales$SoilCarb.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 12)=="SoilCarb.dev")])[,1]
dev.scales$SoilMoist.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 13)=="SoilMoist.dev")])[,1]
# dev.scales$Temp.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 8)=="Temp.dev")])[,1]
# dev.scales$Precip.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 10)=="Precip.dev")])[,1]
# dev.scales$CO2.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 7)=="CO2.dev")])[,1]
dev.scales$Temp.abs.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 12)=="Temp.abs.dev")])[,1]
dev.scales$Precip.abs.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 14)=="Precip.abs.dev")])[,1]
dev.scales$CO2.abs.dev <-stack(ecosys[,which(substr(names(ecosys), 1, 11)=="CO2.abs.dev")])[,1]
summary(dev.scales)

agb.df.ann <- data.frame(AGB.dev=rep(seq(min(ecosys$AGB.dev, na.rm=T), max(ecosys$AGB.dev, na.rm=T), length.out=250), 5), Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="annual")
agb.df.dec <- data.frame(AGB.dev=rep(seq(min(ecosys$AGB.dev, na.rm=T), max(ecosys$AGB.dev, na.rm=T), length.out=250), 5), Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="decadal")
agb.df.cent <- data.frame(AGB.dev=rep(seq(min(ecosys$AGB.dev, na.rm=T), max(ecosys$AGB.dev, na.rm=T), length.out=250), 5), Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="centennial")

agb.df2 <- rbind(agb.df.ann, agb.df.dec, agb.df.cent)
summary(agb.df2)


lm.npp2 <- lm(NPP.dev ~ AGB.dev*Model*scale -1, data=dev.scales[])
summary(lm.npp2)

npp.ci2 <- predict(lm.npp2, interval="confidence", newdata=agb.df2)
summary(npp.ci2)
npp.ci2 <- cbind(agb.df2, npp.ci2)
summary(npp.ci2)




pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP.dev_vs_AGB.dev_Scales.pdf", height=6, width=8)
ggplot(data=dev.scales[dev.scales $Site=="PHA" & ! dev.scales $Model=="JULES", ]) + large.axes +
	facet_grid(scale ~ .) +
	geom_point(aes(x=AGB.dev, y=NPP.dev, color=Model), size=0.5 ) +
	geom_ribbon(data=npp.ci2, aes(x=AGB.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci2, aes(x=AGB.dev, y=fit, color=Model)) +
	scale_y_continuous(limits=c(-0.5, 0.5), name="AGB Deviation") +
	scale_x_continuous(limits=c(-0.25, 0.25), name="NPP Deviation") +
	scale_color_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_fill_manual(values=as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	labs(color="Model") +
	# scale_size_manual(values=c(0.5,1)) +
 	theme(legend.position=c(0.35,0.25), legend.text=element_text(size=12), legend.title=element_text(size=18), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	guides(col=F, alpha=F, size=F, fill=F)
dev.off()








# -----------------------------------------------------------
# -----------------------------------------------------------
# Graphing responses averaged over sites
# -----------------------------------------------------------
# -----------------------------------------------------------

# # -----------------------------------------------------------
# # Driver sensitivities Fig 4
# #	NOTE: these currently have done nothing about temporal autocorrelation
# # -----------------------------------------------------------
# met.stack <- stack(met.year[, which(substr(names(met.year), 1, 8)=="tair.dev" | substr(names(met.year), 1, 11)=="precipf.dev" | substr(names(met.year), 1, 7)=="CO2.dev" | substr(names(met.year),3,10)=="down.dev")])


# met.stack <- stack(met.year[, which(substr(names(met.year), 1, 8)=="tair.dev")])
# names(met.stack) <- c("tair.dev", "scale")
# levels(met.stack$scale) <- c("annual", "decadal", "centennial")
# met.stack$Year <- met.year$Year
# met.stack$Site <- as.factor("PHA")
# met.stack$precip.dev <- stack(met.year[, which(substr(names(met.year), 1, 11)=="precipf.dev")])[,1]
# met.stack$CO2.dev <- stack(met.year[, which(substr(names(met.year), 1, 7)=="CO2.dev")])[,1]
# met.stack$lwdown.dev <- stack(met.year[, which(substr(names(met.year), 1, 10)=="lwdown.dev")])[,1]
# met.stack$swdown.dev <- stack(met.year[, which(substr(names(met.year), 1, 10)=="swdown.dev")])[,1]
# summary(met.stack)

# sensitivity <- stack(ecosys[ecosys$Updated=="Yes" | ecosys$Model=="ED2", which(substr(names(ecosys), 4, 7)==".dev")])
# names(sensitivity) <- c("percent.mean", "var.scale")
# sensitivity[,3:5] <- ecosys[ecosys$Updated=="Yes" | ecosys$Model=="ED2",1:3]
# sensitivity$var <- as.factor(substr(sensitivity$var.scale, 1, 3))
# sensitivity$scale <- as.factor(gsub(".*dev", "t", sensitivity$var.scale))
# levels(sensitivity$scale) <- c("annual", "decadal", "centennial")
# summary(sensitivity)

summary(quru.stack)
quru.stack2 <- stack(quru.stack[,c("RWI", "RWI.Smooth")])
names(quru.stack2) <- c("NPP.dev", "scale")
levels(quru.stack2$scale) <- c("annual", "decadal")
quru.stack2$Year <- as.numeric(paste(quru.stack$Year))
quru.stack2$Model <- as.factor("Tree Rings")
#quru.stack2$var <- as.factor("NPP")
quru.stack2$Site <- as.factor("PHA")
summary(quru.stack2)

quru.stack3 <- merge(quru.stack2, dev.scales[dev.scales$Model=="ED2",c("Year", "scale", "Temp.abs.dev", "Precip.abs.dev", "CO2.abs.dev")], all.x=T, all.y=F)
summary(quru.stack3)
dim(quru.stack2); dim(quru.stack3)
# quru.chron2 <- stack(quru.chron[,c("xxxstd", "Smooth.10")])
# names(quru.chron2) <- c("percent.mean", "scale")
# quru.chron2$Year <- quru.chron$X
# levels(quru.chron2$scale) <- c("decadal", "annual")
# quru.chron2$var <- as.factor("NPP")
# quru.chron2$Model <- as.factor("Tree Rings")
# quru.chron2$Site<-as.factor("PHA")
# quru.chron2$percent.mean <- quru.chron2$percent.mean - 1
# summary(quru.chron2)

sensitivity2 <- merge(dev.scales, quru.stack3, all.x=T, all.y=T)
#sensitivity2$Updated <- as.factor(ifelse(substr(sensitivity2$Model,1,3)=="ED2" | substr(sensitivity2$Model, 1, 3)=="LPJ" | sensitivity2$Model=="Tree Rings", "Yes", "No"))
#sensitivity2[sensitivity2$var=="NPP" & sensitivity2$Model=="Tree Rings",]

sensitivity2$Model.Order <- recode(sensitivity2$Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'; 'LINKAGES'='8'; 'Tree Rings'='9'")
levels(sensitivity2$Model.Order) <- c("ED2", "ED2-LU", "LPJ-WSL", "LPJ-GUESS", "JULES", "LINKAGES", "Tree Rings")
summary(sensitivity2)

# dim(sensitivity); 
dim(sensitivity2)
# Giving linkages Temp & Precip data
sensitivity2[sensitivity2$Model=="LINKAGES", "Temp.abs.dev"] <- sensitivity2[sensitivity2$Model=="LPJ-GUESS", "Temp.abs.dev"] 
sensitivity2[sensitivity2$Model=="LINKAGES", "Precip.abs.dev"] <- sensitivity2[sensitivity2$Model=="LPJ-GUESS", "Precip.abs.dev"] 
sensitivity2[sensitivity2$Model=="LINKAGES", "CO2.abs.dev"] <- sensitivity2[sensitivity2$Model=="LPJ-GUESS", "CO2.abs.dev"] 
summary(sensitivity2)


# sensitivity3  <- merge(sensitivity2, met.stack, all.x=T, all.y=T)
# summary(sensitivity3)
# summary(sensitivity3[sensitivity3$var=="NPP" & sensitivity3$Model=="Tree Rings",])

# model.colors <- c("black", "cyan4", "hotpink3", "green3", "orange3", "steelblack3")

npp.ci.ann <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 7), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 7), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 7), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250), rep("JULES", 250), rep("Tree Rings", 250)), scale="annual")

npp.ci.dec <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 7), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 7), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 7), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250), rep("JULES", 250), rep("Tree Rings", 250)), scale="decadal")

npp.ci.cent <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 6), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 6), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 6), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250), rep("JULES", 250)), scale="centennial")

npp.ci <- rbind(npp.ci.ann, npp.ci.dec, npp.ci.cent)
npp.ci$Model.Order <- recode(npp.ci $Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'; 'LINKAGES'='8'; 'Tree Rings'='9'")
levels(npp.ci$Model.Order) <- c("ED2", "ED2-LU", "LPJ-WSL", "LPJ-GUESS", "JULES", "LINKAGES", "Tree Rings")
summary(npp.ci)

#############################################################################################################################################################################################
# NPP Sensitivity
#############################################################################################################################################################################################
# Temperature
###############################################################
lm.temp <- lm(NPP.dev ~ Temp.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.temp)

npp.ci2 <- predict(lm.temp, interval="confidence", newdata=npp.ci)
colnames(npp.ci2) <- paste("Temp", colnames(npp.ci2), sep=".")
summary(npp.ci2)
npp.ci <- cbind(npp.ci, npp.ci2)

# Precip
lm.precip <- lm(NPP.dev ~ Precip.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.precip)

npp.ci2 <- predict(lm.precip, interval="confidence", newdata=npp.ci)
colnames(npp.ci2) <- paste("Precip", colnames(npp.ci2), sep=".")
summary(npp.ci2)
npp.ci <- cbind(npp.ci, npp.ci2)
summary(npp.ci)

# Precip
lm.co2 <- lm(NPP.dev ~ CO2.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.co2)

npp.ci2 <- predict(lm.co2, interval="confidence", newdata=npp.ci)
colnames(npp.ci2) <- paste("CO2", colnames(npp.ci2), sep=".")
summary(npp.ci2)
npp.ci <- cbind(npp.ci, npp.ci2)
summary(npp.ci)

npp.ci$Updated <- as.factor(ifelse(substr(npp.ci$Model,1,3)=="ED2" | substr(npp.ci$Model, 1, 3)=="LPJ" | npp.ci$Model=="Tree Rings", "Yes", "No"))
npp.ci$Empirical <- as.factor(ifelse(npp.ci$Model=="Tree Rings", "Yes", "No"))
summary(npp.ci)

sensitivity2$Empirical <- as.factor(ifelse(sensitivity2$Model=="Tree Rings", "Yes", "No"))
summary(sensitivity2)


colors.senstivity <- c("dodgerblue4", "deepskyblue2", "navajowhite4", "darkorange2", "deeppink4", "deeppink2", "black")

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Temperature.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model.Order =="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=NPP.dev, color=Model.Order, size=Empirical), size=0.5, alpha=0.5) +
	geom_ribbon(data=npp.ci[!npp.ci$Model =="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill= Model.Order), alpha=0.3) +
	geom_line(data=npp.ci[!npp.ci$Model =="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color= Model.Order), size=1, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$Model =="Tree Rings",], aes(x=Temp.abs.dev, y=NPP.dev, color= Model.Order, size=Empirical), size=0.5) +
	geom_ribbon(data=npp.ci[npp.ci$Model =="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill= Model.Order), alpha=0.5) +
	geom_line(data=npp.ci[npp.ci$Model =="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color= Model.Order), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-3, 1.5), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-1, 1), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



npp.ci2 <- npp.ci
npp.ci2$Temp.lwr <- ifelse(npp.ci2$Temp.lwr > 0.3, 0.3, npp.ci2$Temp.lwr)
npp.ci2$Temp.lwr <- ifelse(npp.ci2$Temp.lwr < -0.5, -0.5, npp.ci2$Temp.lwr)
npp.ci2$Temp.upr <- ifelse(npp.ci2$Temp.upr > 0.3, 0.3, npp.ci2$Temp.upr)
npp.ci2$Temp.upr <- ifelse(npp.ci2$Temp.upr < -0.5, -0.5, npp.ci2$Temp.upr)
summary(npp.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Temperature_2.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=NPP.dev, color= Model.Order, size=Empirical), size=0.5, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=NPP.dev, color= Model.Order, size=Empirical), size=0.5) +
	geom_ribbon(data=npp.ci2[!npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill= Model.Order), alpha=0.3) +
	geom_line(data=npp.ci2[!npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model.Order), size=1, alpha=0.5) +
	geom_ribbon(data=npp.ci2[npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill= Model.Order), alpha=0.5) +
	geom_line(data=npp.ci2[npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color= Model.Order), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.5, 1.0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.3), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Precipitation.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	geom_ribbon(data=npp.ci[!npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci[!npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	geom_ribbon(data=npp.ci[npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci[npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-1, 1), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



npp.ci2 <- npp.ci
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr > 0.5, 0.5, npp.ci2$Precip.lwr)
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr < -0.5, -0.5, npp.ci2$Precip.lwr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr > 0.5, 0.5, npp.ci2$Precip.upr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr < -0.5, -0.5, npp.ci2$Precip.upr)
summary(npp.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Precipitation_2.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	geom_ribbon(data=npp.ci2[!npp.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci2[!npp.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	geom_ribbon(data=npp.ci2[npp.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci2[npp.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, 100), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.5), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()
# -----------------------------------------------------------


# -----------------------------------------------------------
# Block by Scale
# -----------------------------------------------------------

# -----------------------------------------------------------
# Annual
# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Temperature_Annual.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=NPP.dev, color=Model.Order, size=Empirical), size=0.5, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=NPP.dev, color= Model.Order, size=Empirical), size=0.5, alpha=0.05) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="annual" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci2[npp.ci2$scale=="annual" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color= Model.Order), size=1, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="annual" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill= Model.Order), alpha=0.5) +
	geom_line(data=npp.ci2[npp.ci2$scale=="annual" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color= Model.Order), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.5, 1.0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.3), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Precipitation_Annual.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	geom_ribbon(data=npp.ci[npp.ci2$scale=="annual" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci[npp.ci2$scale=="annual" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci[npp.ci2$scale=="annual" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci[npp.ci2$scale=="annual" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, 200), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.5), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

npp.ci2 <- npp.ci
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr > 0.5, 0.5, npp.ci2$Precip.lwr)
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr < -0.5, -0.5, npp.ci2$Precip.lwr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr > 0.5, 0.5, npp.ci2$Precip.upr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr < -0.5, -0.5, npp.ci2$Precip.upr)
summary(npp.ci2)

# -----------------------------------------------------------

# -----------------------------------------------------------
# Decadal
# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Temperature_Decadal.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="decadal" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$scale=="decadal" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="decadal" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci2[npp.ci2$scale=="decadal" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="decadal" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci2[npp.ci2$scale=="decadal" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2, 0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.25), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()


npp.ci2 <- npp.ci
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr > 0.3, 0.3, npp.ci2$Precip.lwr)
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr < -0.5, -0.5, npp.ci2$Precip.lwr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr > 0.3, 0.3, npp.ci2$Precip.upr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr < -0.3, -0.3, npp.ci2$Precip.upr)
summary(npp.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Precipitation_Decadal.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="decadal" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	geom_point(data=sensitivity2[sensitivity2$scale=="decadal" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="decadal" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.3) +
	geom_line(data=npp.ci2[npp.ci2$scale=="decadal" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="decadal" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	geom_line(data=npp.ci[npp.ci2$scale=="decadal" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values= colors.senstivity) +
	scale_fill_manual(values= colors.senstivity) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, -50), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.2), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()
# -----------------------------------------------------------

# -----------------------------------------------------------
# Centennial
# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Temperature_Centennial.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=NPP.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci[npp.ci$scale=="centennial" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=npp.ci[npp.ci$scale=="centennial" & !npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=npp.ci2[npp.ci2$scale=="centennial" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=npp.ci2[npp.ci2$scale=="centennial" & npp.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_fill_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.0, -0.5), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.25), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

npp.ci2 <- npp.ci
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr > 0.25, 0.25, npp.ci2$Precip.lwr)
npp.ci2$Precip.lwr <- ifelse(npp.ci2$Precip.lwr < -0.25, -0.25, npp.ci2$Precip.lwr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr > 0.25, 0.25, npp.ci2$Precip.upr)
npp.ci2$Precip.upr <- ifelse(npp.ci2$Precip.upr < -0.25, -0.25, npp.ci2$Precip.upr)
summary(npp.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_NPP_Deviation_Sensitivity_Precipitation_Centennial.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=npp.ci2[npp.ci2$scale=="centennial" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.3) +
	geom_line(data=npp.ci2[npp.ci2$scale=="centennial" & !npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=NPP.dev, color=Model, size=Empirical), size=0.5) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=npp.ci2[npp.ci2$scale=="centennial" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=npp.ci2[npp.ci2$scale=="centennial" & npp.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_fill_manual(values= as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order),"color"])) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-200, -100), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.4, 0.25), name="NPP Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()
# -----------------------------------------------------------





#############################################################################################################################################################################################
# AGB Sensitivity
#############################################################################################################################################################################################
agb.ci.ann <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 5), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 5), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 5), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="annual")

agb.ci.dec <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 5), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 5), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 5), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="decadal")

agb.ci.cent <- data.frame(Temp.abs.dev=rep(seq(min(sensitivity2$Temp.abs.dev, na.rm=T), max(sensitivity2$Temp.abs.dev, na.rm=T), length.out=250), 5), 
                         Precip.abs.dev=rep(seq(min(sensitivity2$Precip.abs.dev, na.rm=T), max(sensitivity2$Precip.abs.dev, na.rm=T), length.out=250), 5), 
                         CO2.abs.dev=rep(seq(min(sensitivity2$CO2.abs.dev, na.rm=T), max(sensitivity2$CO2.abs.dev, na.rm=T), length.out=250), 5), 
                         Model=c(rep("LPJ-GUESS", 250), rep("LPJ-WSL", 250), rep("ED2", 250), rep("ED2-LU", 250), rep("LINKAGES", 250)), scale="centennial")

agb.ci <- rbind(agb.ci.ann, agb.ci.dec, agb.ci.cent)
agb.ci$Model.Order <- recode(agb.ci$Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'; 'LINKAGES'='8'; 'Tree Rings'='9'")
levels(agb.ci$Model.Order) <- c("ED2", "ED2-LU", "LPJ-WSL", "LPJ-GUESS", "LINKAGES")
summary(agb.ci)

summary(agb.ci)



# Temperature
lm.temp <- lm(AGB.dev ~ Temp.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.temp)

agb.ci2 <- predict(lm.temp, interval="confidence", newdata=agb.ci)
colnames(agb.ci2) <- paste("Temp", colnames(agb.ci2), sep=".")
summary(agb.ci2)
agb.ci <- cbind(agb.ci, agb.ci2)

# Precip
lm.precip <- lm(AGB.dev ~ Precip.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.precip)

agb.ci2 <- predict(lm.precip, interval="confidence", newdata=agb.ci)
colnames(agb.ci2) <- paste("Precip", colnames(agb.ci2), sep=".")
summary(agb.ci2)
agb.ci <- cbind(agb.ci, agb.ci2)
summary(agb.ci)

# Precip
lm.co2 <- lm(AGB.dev ~ CO2.abs.dev*Model*scale - 1, data=sensitivity2[,])
summary(lm.co2)

agb.ci2 <- predict(lm.co2, interval="confidence", newdata=agb.ci)
colnames(agb.ci2) <- paste("CO2", colnames(agb.ci2), sep=".")
summary(agb.ci2)
agb.ci <- cbind(agb.ci, agb.ci2)
summary(agb.ci)

agb.ci$Updated <- as.factor(ifelse(substr(agb.ci$Model,1,3)=="ED2" | substr(agb.ci$Model, 1, 3)=="LPJ" | agb.ci$Model=="Tree Rings", "Yes", "No"))
agb.ci$Empirical <- as.factor(ifelse(agb.ci$Model=="Tree Rings", "Yes", "No"))
summary(agb.ci)

sensitivity2$Empirical <- as.factor(ifelse(sensitivity2$Model=="Tree Rings", "Yes", "No"))
summary(sensitivity2)


pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings" & !sensitivity2$Model=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci[!agb.ci$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci[!agb.ci$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1) +
	# geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-3, 1.5), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-1, 1), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



agb.ci2 <- agb.ci
agb.ci2$Temp.lwr <- ifelse(agb.ci2$Temp.lwr > 0.3, 0.3, agb.ci2$Temp.lwr)
agb.ci2$Temp.lwr <- ifelse(agb.ci2$Temp.lwr < -0.5, -0.5, agb.ci2$Temp.lwr)
agb.ci2$Temp.upr <- ifelse(agb.ci2$Temp.upr > 0.3, 0.3, agb.ci2$Temp.upr)
agb.ci2$Temp.upr <- ifelse(agb.ci2$upr < -0.3, -0.3, agb.ci2$Temp.upr)
summary(agb.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_2.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[,], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[,], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.5, 1.0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.3), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci[!agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci[!agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=1) +
	# geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-1, 1), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()



agb.ci2 <- agb.ci
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr > 0.5, 0.5, agb.ci2$Precip.lwr)
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr < -0.5, -0.5, agb.ci2$Precip.lwr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr > 0.5, 0.5, agb.ci2$Precip.upr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr < -0.5, -0.5, agb.ci2$Precip.upr)
summary(agb.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_2.pdf", height=6, width=8)
ggplot(data=sensitivity2[!sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[!agb.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[!agb.ci2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=1) +
	# geom_point(data=sensitivity2[sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, 200), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.5, 0.5), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()
# -----------------------------------------------------------


# -----------------------------------------------------------
# Block by Scale
# -----------------------------------------------------------
large.axes2 <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=16), axis.text.y=element_text(color="black", size=16), axis.title.x=element_text(face="bold", size=18, vjust=-1),  axis.title.y=element_text(face="bold", size=18, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

# -----------------------------------------------------------
# Annual
# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Annual.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="annual" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="annual" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1, alpha=1) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="annual" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="annual" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.5, 1.0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Annual.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci[agb.ci2$scale=="annual" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci[agb.ci2$scale=="annual" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci[agb.ci2$scale=="annual" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci2$scale=="annual" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, 200), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.20), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

# --------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Annual_2.pdf", height=5, width=8)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="annual" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="annual" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1, alpha=1) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="annual" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="annual" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.5, 1.0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F, color=F)
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Annual_2.pdf", height=5, width=8)
ggplot(data=sensitivity2[sensitivity2$scale=="annual" & !sensitivity2$Model=="Tree Rings"& !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci[agb.ci2$scale=="annual" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci[agb.ci2$scale=="annual" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="annual" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci[agb.ci2$scale=="annual" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci2$scale=="annual" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, 200), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.20), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F, color=F)
dev.off()
# -----------------------------------------------------------

# -----------------------------------------------------------
# Decadal
# -----------------------------------------------------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Decadal.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="decadal" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="decadal" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="decadal" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="decadal" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1, alpha=1) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="decadal" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="decadal" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2, 0), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()


agb.ci2 <- agb.ci
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr > 0.2, 0.2, agb.ci2$Precip.lwr)
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr < -0.3, -0.3, agb.ci2$Precip.lwr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr > 0.2, 0.2, agb.ci2$Precip.upr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr < -0.3, -0.3, agb.ci2$Precip.upr)
summary(agb.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Decadal.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="decadal" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="decadal" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="decadal" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1, alpha=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="decadal" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="decadal" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci[agb.ci2$scale=="decadal" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-250, -50), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()
# -----------------------------------------------------------

# -----------------------------------------------------------
# Centennial
# -----------------------------------------------------------

agb.ci2 <- agb.ci
agb.ci2$Temp.lwr <- ifelse(agb.ci2$Temp.lwr > 0.3, 0.3, agb.ci2$Temp.lwr)
agb.ci2$Temp.lwr <- ifelse(agb.ci2$Temp.lwr < -0.3, -0.3, agb.ci2$Temp.lwr)
agb.ci2$Temp.upr <- ifelse(agb.ci2$Temp.upr > 0.3, 0.3, agb.ci2$Temp.upr)
agb.ci2$Temp.upr <- ifelse(agb.ci2$Temp.upr < -0.3, -0.3, agb.ci2$Temp.upr)
summary(agb.ci2)
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Centennial.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.0, -0.5), name="Temperature Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()

agb.ci2 <- agb.ci
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr > 0.2, 0.2, agb.ci2$Precip.lwr)
agb.ci2$Precip.lwr <- ifelse(agb.ci2$Precip.lwr < -0.3, -0.3, agb.ci2$Precip.lwr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr > 0.2, 0.2, agb.ci2$Precip.upr)
agb.ci2$Precip.upr <- ifelse(agb.ci2$Precip.upr < -0.3, -0.3, agb.ci2$Precip.upr)
summary(agb.ci2)

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Centennial.pdf", height=3, width=9)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.3) +
	geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-200, -100), name="Precipitation Deviation") +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation") +
	theme_bw() + guides(alpha=F, fill=F)
dev.off()


# ----------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Centennial_2.pdf", height=5, width=8)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.0, -0.5), name=expression(bold(paste("Temperature Deviation ("^"o","C)")))) +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation (Fraction)") +
	theme_bw() + guides(alpha=F, fill=F, color=F) + 
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank()) +
	theme(axis.text.x=element_text(angle=0, color="black", size=14), axis.text.y=element_text(color="black", size=14), axis.title.x=element_text(face="bold", size=14, vjust=0),  axis.title.y=element_text(face="bold", size=14, vjust=1))
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Centennial_2.pdf", height=5, width=8)
ggplot(data=sensitivity2[sensitivity2$scale=="centennial" & !sensitivity2$Model=="Tree Rings" & !sensitivity2$Model.Order=="JULES",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.3) +
	geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-200, -100), name=expression(bold(paste("Precipitation Deviation (mm yr"^"-1",")")))) +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation (Fraction)") +
	theme_bw() + guides(alpha=F, fill=F, color=F) + 
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank()) +
	theme(axis.text.x=element_text(angle=0, color="black", size=14), axis.text.y=element_text(color="black", size=14), axis.title.x=element_text(face="bold", size=14, vjust=0),  axis.title.y=element_text(face="bold", size=14, vjust=1))
dev.off()


# -----------------------------------------------------------axis.text.x=element_text(angle=0, 


# ----------------
pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Temperature_Annual_Centennial.pdf", height=11.5, width=9.5)
ggplot(data=sensitivity2[(sensitivity2$scale=="centennial" | sensitivity2$scale=="annual") & !sensitivity2$Model.Order=="JULES" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5, alpha=0.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Temp.abs.dev, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[(agb.ci2$scale=="centennial" | agb.ci2$scale=="annual") & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=Temp.lwr, ymax=Temp.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[(agb.ci2$scale=="centennial" | agb.ci2$scale=="annual")  & !agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=Temp.fit, color=Model), size=1.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, ymin=lwr, ymax=upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci2$Model=="Tree Rings",], aes(x= Temp.abs.dev, y=fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-2.0, -0.5), name=expression(bold(paste("Temperature Deviation ("^"o","C)")))) +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation (Fraction)") +
	theme_bw() + guides(alpha=F, fill=F, color=F) + 
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank()) +
	theme(axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=0),  axis.title.y=element_text(face="bold", size=24, vjust=1)) + 
	theme(strip.text=element_text(size=24))
dev.off()

pdf("~/Dropbox/PalEON CR/paleon_mip_site/Analyses/Temporal-Scaling/HF_Talk_2015.04/PHA_AGB_Deviation_Sensitivity_Precipitation_Annual_Centennial.pdf", height=11.5, width=9.5)
ggplot(data=sensitivity2[(sensitivity2$scale=="centennial" | sensitivity2$scale=="annual") & !sensitivity2$Model.Order=="JULES" & !sensitivity2$Model=="Tree Rings",]) + 
	facet_grid(scale ~ .) +
	geom_point(aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# geom_smooth(method=lm, alpha=0.5, size=1) + 
	geom_ribbon(data=agb.ci2[(agb.ci2$scale=="centennial" | agb.ci2$scale=="annual") & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model.Order), alpha=0.5) +
	geom_line(data=agb.ci2[(agb.ci2$scale=="centennial" | agb.ci2$scale=="annual") & !agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	# geom_point(data=sensitivity2[sensitivity2$scale=="centennial" & sensitivity2$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=AGB.dev, color=Model, size=Empirical), size=0.5) +
	# # geom_smooth(method=lm, alpha=0.5, size=1) + 
	# geom_ribbon(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, ymin=Precip.lwr, ymax=Precip.upr, fill=Model), alpha=0.5) +
	# geom_line(data=agb.ci2[agb.ci2$scale=="centennial" & agb.ci$Model=="Tree Rings",], aes(x=Precip.abs.dev*sec2yr, y=Precip.fit, color=Model), size=1.5) +
	scale_color_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_fill_manual(values=c(as.vector(model.colors[model.colors$Model.Order %in% unique(ecosys$Model.Order) & !model.colors$Model.Order=="JULES","color"]))) +
	scale_size_manual(values=c(0.1, 10)) + 
	scale_alpha_manual(values=c(0.2, 0.5)) + 
	scale_x_continuous(limits=c(-200, -100), name=expression(bold(paste("Precipitation Deviation (mm yr"^"-1",")")))) +
	scale_y_continuous(limits=c(-0.3, 0.2), name="AGB Deviation (Fraction)") +
	theme_bw() + guides(alpha=F, fill=F, color=F) + 
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank()) +
	theme(axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=0),  axis.title.y=element_text(face="bold", size=24, vjust=1)) + 
	theme(strip.text=element_text(size=24))
dev.off()

