# ----------------------------------------------------
# Graphs & Analysis for Poster at 2015 NACP Meeting
# Christine R. Rollinson
# crollinson@gmail.com 
# 18 Jan 2015
# ----------------------------------------------------

library(ncdf4)
library(ggplot2)
library(grid)
library(nlme)

setwd("~/Dropbox/PalEON CR/paleon_mip_site")
outputs <- "phase1a_output_variables"
years <- 850:2010

large.axes <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=24, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

# -------------------------------------------------------------------
# Reading in .nc files of variables created using script in MIP_formatting_ModelLoop_Yr.R
# -------------------------------------------------------------------
gpp <- nc_open(file.path(outputs, "GPP.annual.nc"))
agb <- nc_open(file.path(outputs, "AGB.annual.nc"))
lai <- nc_open(file.path(outputs, "LAI.annual.nc"))
npp <- nc_open(file.path(outputs, "NPP.annual.nc"))
nee <- nc_open(file.path(outputs, "NEE.annual.nc"))
temp <- nc_open(file.path(outputs, "Temp.annual.nc"))
precip <- nc_open(file.path(outputs, "Precip.annual.nc"))
auto.resp <- nc_open(file.path(outputs, "RespirationAuto.annual.nc"))
hetero.resp <- nc_open(file.path(outputs, "RespirationHetero.annual.nc"))
soilmoist <- nc_open(file.path(outputs, "SoilMoist.annual.nc"))
soilcarb <- nc_open(file.path(outputs, "TotSoilCarb.annual.nc"))
transp <- nc_open(file.path(outputs, "Transp.annual.nc"))
evap <- nc_open(file.path(outputs, "Evap.annual.nc"))

gpp.sib <- read.csv("phase1a_model_output/Version1_OldMet/Sib-CASA.v0/Sib-CASA_GPP_Month.csv") #kg/m2/s
agb.sib <- read.csv("phase1a_model_output/Version1_OldMet/Sib-CASA.v0/Sib-CASA_AGB_Month.csv") #kg/m2
lai.sib <- read.csv("phase1a_model_output/Version1_OldMet/Sib-CASA.v0/Sib-CASA_LAI_Month.csv") #kg/m2
nee.sib <- read.csv("phase1a_model_output/Version1_OldMet/Sib-CASA.v0/Sib-CASA_NEE_Month.csv") #kg/m2
npp.sib <- read.csv("phase1a_model_output/Version1_OldMet/Sib-CASA.v0/Sib-CASA_NPP_Month.csv") #kg/m2

summary(agb.sib) 

yr.rows <- seq(1, nrow(gpp.sib)-12, by=12)

gpp.sib.yr <- agb.sib.yr <- lai.sib.yr <- npp.sib.yr <- nee.sib.yr <- data.frame(array(dim=c(nrow(gpp.sib)/12, ncol(gpp.sib))))
names(gpp.sib.yr) <- names(agb.sib.yr) <- names(lai.sib.yr) <- names(npp.sib.yr) <- names(nee.sib.yr) <- names(gpp.sib)

for(s in 1:ncol(gpp.sib)){
	gpp.temp <-	agb.temp <- lai.temp <- npp.temp <- nee.temp <- vector()
	for(i in 1:length(yr.rows)){
		gpp.temp <- c(gpp.temp, mean(gpp.sib[yr.rows[i]:(yr.rows[i]+11),s],na.rm=T))
		agb.temp <- c(agb.temp, mean(agb.sib[yr.rows[i]:(yr.rows[i]+11),s], na.rm=T))
		lai.temp <- c(lai.temp, mean(lai.sib[yr.rows[i]:(yr.rows[i]+11),s], na.rm=T))
		npp.temp <- c(npp.temp, mean(npp.sib[yr.rows[i]:(yr.rows[i]+11),s], na.rm=T))
		nee.temp <- c(nee.temp, mean(nee.sib[yr.rows[i]:(yr.rows[i]+11),s], na.rm=T))
	}
	gpp.sib.yr[,s] <- c(gpp.temp, NA)
	agb.sib.yr[,s] <- c(agb.temp, NA)
	lai.sib.yr[,s] <- c(lai.temp, NA)
	npp.sib.yr[,s] <- c(npp.temp, NA)
	nee.sib.yr[,s] <- c(nee.temp, NA)
}
summary(gpp.sib.yr)
summary(agb.sib.yr)
summary(lai.sib.yr)
summary(npp.sib.yr)
summary(nee.sib.yr)


sites <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB")
model.names <- c(ncvar_get(gpp, "ModelNames"), "SiB")
#sites <- "PHA"
GPP <- AGB <- LAI <- NPP <- NEE <- TEMP <- PRECIP <- RA <- RH <- SOILMOIST <- SOILCARB <- EVAP <- TRANSP <- list()
for(i in 1:length(sites)){
	GPP[[i]] <- data.frame(cbind(t(ncvar_get(gpp, names(gpp$var)[i])), gpp.sib.yr[,i]))
	AGB[[i]] <- data.frame(cbind(t(ncvar_get(agb, names(agb$var)[i])), agb.sib.yr[,i]))
	LAI[[i]] <- data.frame(cbind(t(ncvar_get(lai, names(lai$var)[i])), lai.sib.yr[,i]))
	NPP[[i]] <- data.frame(cbind(t(ncvar_get(npp, names(npp$var)[i])), npp.sib.yr[,i]))
	NEE[[i]] <- data.frame(cbind(t(ncvar_get(nee, names(nee$var)[i])), nee.sib.yr[,i]))
	TEMP[[i]] <- data.frame(cbind(t(ncvar_get(temp, names(temp$var)[i])), ncvar_get(temp, names(temp$var)[i])[1,]))
	PRECIP[[i]] <- data.frame(cbind(t(ncvar_get(precip, names(precip$var)[i])), ncvar_get(precip, names(precip$var)[i])[1,]))
	RA[[i]] <- data.frame(cbind(t(ncvar_get(auto.resp, names(auto.resp$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))
	RH[[i]] <- data.frame(cbind(t(ncvar_get(hetero.resp, names(hetero.resp$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))
	SOILCARB[[i]] <- data.frame(cbind(t(ncvar_get(soilcarb, names(soilcarb$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))
	SOILMOIST[[i]] <- data.frame(cbind(t(ncvar_get(soilmoist, names(soilmoist$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))
	EVAP[[i]] <- data.frame(cbind(t(ncvar_get(evap, names(evap$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))
	TRANSP[[i]] <- data.frame(cbind(t(ncvar_get(transp, names(transp$var)[i]))), SiB=rep(NA, length(nee.sib.yr[,i])))

	names(GPP[[i]]) <- names(AGB[[i]]) <- names(LAI[[i]]) <- names(NPP[[i]]) <- names(NEE[[i]]) <- names(TEMP[[i]]) <- names(PRECIP[[i]]) <- names(RA[[i]]) <- names(RH[[i]]) <- names(SOILCARB[[i]]) <- names(SOILMOIST[[i]]) <- names(EVAP[[i]]) <- names(TRANSP[[i]]) <- model.names
	
	row.names(GPP[[i]]) <- row.names(AGB[[i]]) <- row.names(TEMP[[i]]) <- row.names(PRECIP[[i]]) <- row.names(RA[[i]]) <- row.names(RH[[i]]) <- row.names(SOILCARB[[i]]) <- row.names(SOILMOIST[[i]]) <- row.names(EVAP[[i]]) <- row.names(TRANSP[[i]]) <- years

}
names(GPP) <- names(AGB) <- names(LAI) <- names(NPP) <- names(NEE) <- names(TEMP) <- names(PRECIP) <- names(RA) <- names(RH) <- names(SOILCARB) <- names(SOILMOIST) <- names(EVAP) <- names(TRANSP) <- sites

nc_close(gpp); nc_close(agb); nc_close(lai); nc_close(npp); nc_close(nee); nc_close(temp); nc_close(precip); nc_close(auto.resp); nc_close(hetero.resp); nc_close(soilcarb); nc_close(soilmoist); nc_close(evap); nc_close(transp)

df1 <- data.frame(array(dim=c(0,10)))
for(i in 1:length(sites)){
	gpp.df <- stack(GPP[[i]])
	names(gpp.df) <- c("GPP", "Model")
	gpp.df$Site <- as.factor(names(GPP)[i])
	gpp.df$Year <- years

	agb.df <- stack(AGB[[i]])
	names(agb.df) <- c("AGB", "Model")
	agb.df$Site <- as.factor(names(AGB)[i])

	lai.df <- stack(LAI[[i]])
	names(lai.df) <- c("LAI", "Model")
	lai.df$Site <- as.factor(names(LAI)[i])

	npp.df <- stack(NPP[[i]])
	names(npp.df) <- c("NPP", "Model")
	npp.df$Site <- as.factor(names(NPP)[i])

	nee.df <- stack(NEE[[i]])
	names(nee.df) <- c("NEE", "Model")
	nee.df$Site <- as.factor(names(NEE)[i])

	ra.df <- stack(RA[[i]])
	names(ra.df) <- c("AutoResp", "Model")
	ra.df$Site <- as.factor(names(RA)[i])

	rh.df <- stack(RH[[i]])
	names(rh.df) <- c("HeteroResp", "Model")
	rh.df$Site <- as.factor(names(RH)[i])

	soilcarb.df <- stack(SOILCARB[[i]])
	names(soilcarb.df) <- c("SoilCarb", "Model")
	soilcarb.df$Site <- as.factor(names(SOILCARB)[i])

	soilmoist.df <- stack(SOILMOIST[[i]])
	names(soilmoist.df) <- c("SoilMoist", "Model")
	soilmoist.df$Site <- as.factor(names(SOILMOIST)[i])

	evap.df <- stack(EVAP[[i]])
	names(evap.df) <- c("Evap", "Model")
	evap.df$Site <- as.factor(names(EVAP)[i])

	transp.df <- stack(TRANSP[[i]])
	names(transp.df) <- c("Transp", "Model")
	transp.df$Site <- as.factor(names(TRANSP)[i])
	

	temp.new <- data.frame(Year = years, Temp = TEMP[[i]][,"lpj.guess"], Model=c(rep("ed2", nrow(TEMP[[i]])), rep("ed2.lu", nrow(TEMP[[i]])), rep("lpj.guess", nrow(TEMP[[i]])), rep("lpj.wsl", nrow(TEMP[[i]]))))
	precip.new <- data.frame(Year = years, Precip = PRECIP[[i]][,"lpj.guess"], Model=c(rep("ed2", nrow(PRECIP[[i]])), rep("ed2.lu", nrow(TEMP[[i]])), rep("lpj.guess", nrow(PRECIP[[i]])), rep("lpj.wsl", nrow(PRECIP[[i]]))))
	

	temp.old <- data.frame(Year = years, Temp = TEMP[[i]][,"jules.stat"], Model=c(rep("clm45", nrow(TEMP[[i]])), rep("jules.stat", nrow(TEMP[[i]])), rep("SiB", nrow(TEMP[[i]]))))
	precip.old <- data.frame(Year = years, Precip = PRECIP[[i]][,"jules.stat"], Model=c(rep("clm45", nrow(PRECIP[[i]])), rep("jules.stat", nrow(PRECIP[[i]])), rep("SiB", nrow(PRECIP[[i]]))))

	temp.df <- rbind(temp.new, temp.old)
	precip.df <- rbind(precip.new, precip.old)
	
	climate.df <- temp.df
	climate.df$Site <- sites[i]
	climate.df$Precip <- precip.df$Precip

	df.temp <- cbind(gpp.df, agb.df$AGB, lai.df$LAI, npp.df$NPP, nee.df$NEE, ra.df$AutoResp, rh.df$HeteroResp, soilcarb.df$SoilCarb, soilmoist.df$SoilMoist, evap.df$Evap, transp.df$Transp)
	df.temp <- merge(df.temp, climate.df, all.x=T, all.y=T)

	df1 <- rbind(df1, df.temp)
	}
names(df1) <- c("Model", "Site", "Year", "GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "SoilCarb", "SoilMoist", "Evap", "Transp", "Temp", "Precip")
summary(df1)


models.updated <- c("lpj.guess", "lpj.wsl", "ed2", "ed2.lu", "jules.stat", "linkages")
df1$Updated <- as.factor(ifelse(df1$Model %in% models.updated, "Yes", "No"))
summary(df1)

write.csv(df1, file.path(outputs, "MIP_Data_Ann_2015.csv"), row.names=F)
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Graphing Variables at PHA for Poster
# -------------------------------------------------------------------
df1 <- read.csv(file.path(outputs, "MIP_Data_Ann_2015.csv"))
#df1 <- df1[df1$Site=="PHA", 2:ncol(df1)]
summary(df1)

levels(df1$Model) <- c("CLM4.5", "ED2", "ED2-LU", "JULES", "LPJ-GUESS", "LPJ-WSL", "SiBCASA")

sec_2_yr <- 60*60*24*365
kgm2_2_MgHa <- (1/1000) * 10000

MgCHayr <- sec_2_yr*kgm2_2_MgHa

model.colors <- c("black", "gray50", "cyan4", "hotpink3", "green3", "orange3", "steelblue3")

library(car)
df1$Model.Order <- recode(df1$Model, "'ED2'='1'; 'ED2-LU'='2'; 'CLM4.5'='3'; 'LPJ-WSL'='4'; 'LPJ-GUESS'='5'; 'JULES'='6'; 'SiBCASA'='7'")
levels(df1$Model.Order) <- c("ED2", "ED2-LU", "CLM4.5", "LPJ-WSL", "LPJ-GUESS", "JULES", "SibCASA")
summary(df1)

#sec_2_mo
# Precip Conversion: right now we have mean monthly rate/mo = mm/yr
#	A seconds to year = 60 sec/min * 60 min/hr * 24 hr/day * 365 day/yr = s/y
#	B kg to mm3 = 1000 g/kg = 1000 cm3 * 10*10*10 mm3/cm3 = mm3/kg
#	C m2 to mm2 = 1e6 mm2/m2 = mm2/m2
# kg m-2 s-1 = B * A * 1/c 

precip.conversion <- 1 * ((60*60*24*365)) * (1e3 * 1e3) * (1e-6)

summary(df1$Precip*60*60*24*30.4) # so this is mean monthly precip 
precip.month <- df1$Precip*60*60*24*30.4 # kg/m2/mo = cm/mo

precip.month <- df1$Precip*precip.conversion
summary(precip.month)
summary(precip.month*12/10)




# pdf("Figures/GPP_PHA.pdf", width=10, height=6)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=GPP*MgCHayr, color=Model.Order), size=1.5) +
	scale_y_continuous(name=expression(bold(paste("GPP (MgC ha"^"-1"," yr"^"-1",")")))) +
	scale_x_continuous(name="Year") +
	scale_color_manual(values=model.colors) + labs(color="Models") #+
 	# theme(legend.position=c(0.25,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2))
# dev.off()

# pdf("Figures/NACP2015_AGB_PHA.pdf", width=10, height=6)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=AGB*kgm2_2_MgHa, color=Model.Order), size=1.5) +
	scale_y_continuous(name=expression(bold(paste("AGB (MgC ha"^"-1",")")))) +
	scale_x_continuous(name="Year") +
	scale_color_manual(values=model.colors) + labs(color="Models")  + 
 	theme(axis.title.y=element_text(vjust=1)) +
	guides(color=F)
# dev.off()	

# pdf("Figures/NACP2015_LAI_PHA.pdf", width=10, height=6)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=LAI, color=Model.Order), size=1.5) +
	scale_y_continuous(name=expression(bold(paste("LAI (m"^"2"," m"^"-2",")")))) +
	scale_x_continuous(name="Year") +
	scale_color_manual(values=model.colors) + labs(color="Models")  + 
 	theme(axis.title.y=element_text(vjust=3)) +
	guides(color=F)
# dev.off()	

# pdf("Figures/NACP2015_NPP_PHA.pdf", width=10, height=6)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=NPP*MgCHayr, color=Model.Order), size=1.5) +
	scale_y_continuous(name=expression(bold(paste("NPP (MgC ha"^"-1"," yr"^"-1",")")))) +
	scale_x_continuous(name="Year") +
	scale_color_manual(values=model.colors) + labs(color="Models") +
 	theme(legend.position=c(0.25,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + 
 	theme(axis.title.y=element_text(vjust=3)) +
 	guides(col=guide_legend(ncol=2))
# dev.off()	

# pdf("Figures/NACP2015_NEE_PHA.pdf", width=10, height=6)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=NEE* MgCHayr, color=Model.Order), size=1.5) +
	scale_y_continuous(name=expression(bold(paste("NEE (MgC ha"^"-1"," yr"^"-1",")")))) +
	scale_x_continuous(name="Year") +
	scale_color_manual(values=model.colors) + labs(color="Models") +
	guides(color=F)
# dev.off()	




# pdf("Figures/NACP2015_Temperature_ED2_PHA.pdf", width=10, height=3)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=Temp-273.15, color=Updated), size=0.8) +
	scale_y_continuous(name=expression(bold(paste("Temp ("^"o","C)")))) +
	scale_x_continuous(name="Year")  +
	theme(axis.text.x=element_text(angle=0, color="black", size=18), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=28, vjust=-1),  axis.title.y=element_text(face="bold", size=20, vjust=2))
# dev.off()	

# pdf("Figures/NACP2015_Precip_Rate_ED2_PHA.pdf", width=10, height=3)
ggplot(data=df1[df1$Site=="PHA",]) + large.axes +
	geom_line(aes(x=Year, y=Precip*precip.conversion*1e-3, color=Updated), size=0.8) +
	scale_y_continuous(name=expression(bold(paste("Precip (m yr"^"-1",")"))), breaks=c(0.8,1,1.2, 1.4, 1.6)) +
	# scale_y_continuous(name=expression(bold(paste("Precip (m/yr)")))) +
	scale_x_continuous(name="Year") +
	theme(axis.text.x=element_text(angle=0, color="black", size=18), axis.text.y=element_text(color="black", size=18), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=20, vjust=1.5))
# dev.off()	
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Finding Responses to Extreme Events
# -------------------------------------------------------------------
df1 <- read.csv(file.path(outputs, "MIP_Data_Ann_2015.csv"))
summary(df1)

# Finding the 5th and 95th percentiles for extreme events
precip.extremes <- temp.extremes <- data.frame(array(NA, dim=c(2,6)))
names(precip.extremes) <- names(temp.extremes) <- unique(df1$Site)

for(i in names(temp.extremes)){
	temp.extremes[,i] <- quantile(df1[df1$Site==i, "Temp"], c(0.05, 0.95))
	precip.extremes[,i] <- quantile(df1[df1$Site==i, "Precip"], c(0.05, 0.95))
}
temp.extremes
precip.extremes


# Mean Centering the GPP & AGB Data
for(i in unique(df1$Model)){
	for(j in unique(df1$Site)){
		mean.gpp <- mean(df1[df1$Model==i & df1$Site==j,"GPP"], na.rm=T)
		mean.agb <- mean(df1[df1$Model==i & df1$Site==j,"AGB"], na.rm=T)
		mean.temp <- mean(df1[df1$Model==i & df1$Site==j,"Temp"], na.rm=T)
		mean.precip <- mean(df1[df1$Model==i & df1$Site==j,"Precip"], na.rm=T)
		df1[df1$Model==i & df1$Site==j,"GPP.centered"] <- (df1[df1$Model==i & df1$Site==j,"GPP"] - mean.gpp)/mean.gpp
		df1[df1$Model==i & df1$Site==j,"AGB.centered"] <- (df1[df1$Model==i & df1$Site==j,"AGB"] - mean.agb)/mean.agb
		df1[df1$Model==i & df1$Site==j,"Temp.centered"] <- (df1[df1$Model==i & df1$Site==j,"Temp"] - mean.temp)/mean.temp
		df1[df1$Model==i & df1$Site==j,"Precip.centered"] <- (df1[df1$Model==i & df1$Site==j,"Precip"] - mean.precip)/mean.precip
	}
}
summary(df1)
write.csv(df1, file.path(outputs, "MIP_Data_Ann_2015.csv"))

# ggplot(data=df1) + large.axes + facet_grid(Model ~ .) +
	# geom_point(aes(x=Temp.centered, y=GPP.centered, color=Site), size=0.8)

# ggplot(data=df1) + large.axes + facet_grid(Model ~ .) +
	# geom_point(aes(x=Precip.centered, y=GPP.centered, color=Site), size=0.8)

# ggplot(data=df1) + large.axes + facet_grid(Model ~ .) +
	# geom_point(aes(x=Temp.centered, y=AGB.centered, color=Site))

# ggplot(data=df1) + large.axes + facet_grid(Model ~ .) +
	# geom_point(aes(x=Precip.centered, y=AGB.centered, color=Site))


temp.max <- data.frame(array(dim=c(0,ncol(df1))))
temp.min <- data.frame(array(dim=c(0,ncol(df1))))
precip.max <- data.frame(array(dim=c(0,ncol(df1))))
precip.min <- data.frame(array(dim=c(0,ncol(df1))))
for(s in unique(df1$Site)){
	temp.max <- rbind(temp.max, df1[df1$Site==s & df1$Temp > temp.extremes[2,s],])
	temp.min <- rbind(temp.min, df1[df1$Site==s & df1$Temp < temp.extremes[1,s],])
	precip.max <- rbind(precip.max, df1[df1$Site==s & df1$Precip > precip.extremes[2,s],])
	precip.min <- rbind(precip.min, df1[df1$Site==s & df1$Precip < precip.extremes[1,s],])
}
temp.max$Anomaly <- as.factor("Tmax")
temp.min$Anomaly <- as.factor("Tmin")
precip.max$Anomaly <- as.factor("Pmax")
precip.min$Anomaly <- as.factor("Pmin")
summary(temp.max)
summary(temp.min)
summary(precip.max)
summary(precip.min)

# Some quick t-tests
t.test(df1$GPP.centered, temp.max$GPP.centered)
t.test(df1$GPP.centered, precip.min$GPP.centered)
t.test(df1$GPP.centered, precip.max$GPP.centered)
t.test(df1$GPP.centered, temp.min$GPP.centered)


lm.tmax <- lm(GPP.centered ~ Model -1, data=temp.max)
lm.tmin <- lm(GPP.centered ~ Model -1, data=temp.min)
summary(lm.tmax)
summary(lm.tmin)

lm.pmax <- lm(GPP.centered ~ Model - 1, data=precip.max)
lm.pmin <- lm(GPP.centered ~ Model - 1, data=precip.min)
summary(lm.pmax)
summary(lm.pmin)

t.test(abs(precip.min$GPP.centered), abs(temp.max$GPP.centered))

# Comparing sensitivities of AGB & Precip
t.test(abs(precip.min$AGB.centered), abs(precip.min$GPP.centered), paired=T)
t.test(abs(temp.max$AGB.centered), abs(temp.max$GPP.centered), paired=T)
t.test(abs(precip.max$AGB.centered), abs(precip.max$GPP.centered), paired=T)
t.test(abs(temp.min$AGB.centered), abs(temp.min$GPP.centered), paired=T)


extremes <- rbind(temp.max, temp.min, precip.max, precip.min)
summary(extremes)

ggplot(data=extremes) + large.axes + facet_grid(Anomaly ~ .) +
	geom_boxplot(aes(x=Model, y=GPP.centered), color="red", fill="red")


extremes2 <- aggregate(extremes[,c("GPP.centered", "AGB.centered")], by=list(extremes$Anomaly, extremes$Model), FUN=mean, na.rm=T)
names(extremes2) <- c("Anomaly", "Model", "GPP", "AGB")
summary(extremes2)


extremes2sd <- aggregate(extremes[,c("GPP.centered", "AGB.centered")], by=list(extremes$Anomaly, extremes$Model), FUN=sd, na.rm=T)
names(extremes2sd) <- c("Anomaly", "Model", "GPP", "AGB")
summary(extremes2sd)


extremes2$GPP.sd <- extremes2sd$GPP
extremes2$AGB.sd <- extremes2sd$AGB

ggplot(data=extremes2) + large.axes + facet_grid(Anomaly ~ .) +
	geom_point(aes(x=Model, y=GPP), size=5) +
	geom_hline(aes(yintercept=0), size=0.1) +
	geom_errorbar(aes(x=Model, ymax=GPP+GPP.sd, ymin=GPP-GPP.sd), width=0.25)

levels(extremes2$Anomaly) <- c("Max Temp", "Min Temp", "Max Precip", "Min Precip")

# pdf("Figures/NACP2015_GPP_Deviation.pdf", width=12, height=7)
ggplot(data=extremes2[extremes2$Anomaly=="Max Temp" | extremes2$Anomaly=="Min Precip",]) + large.axes + facet_grid(Anomaly ~ .) +
	geom_point(aes(x=Model, y=GPP*100, color=Anomaly), size=5) +
	geom_hline(aes(yintercept=0), size=0.1) +
	geom_errorbar(aes(x=Model, ymax=(GPP+GPP.sd)*100, ymin=(GPP-GPP.sd)*100, color=Anomaly), width=0.25, size=1.5) +
	scale_y_continuous(name="Percent Change in GPP") +
	scale_x_discrete(name="Model") + 
	scale_color_manual(values=c("red", "blue")) + guides(color=F) +
	theme(strip.text.y=element_text(size=18, face="bold"), strip.background=element_rect(fill="gray80"))
# dev.off()



by(temp.max$GPP.centered, temp.max$Model, FUN=mean)
by(temp.min$GPP.centered, temp.max$Model, FUN=mean)


dim(temp.max)
dim(temp.min)
dim(precip.max)
dim(precip.min)



ggplot() + large.axes +
	geom_boxplot(data=temp.max, aes(x=Model, y=GPP.centered), color="red", fill="red") +
	geom_boxplot(data=temp.min, aes(x=Model, y=GPP.centered), color="blue", fill="NA")

ggplot() + large.axes +
	geom_boxplot(data=precip.min, aes(x=Model, y=GPP.centered), color="red", fill="red") +
	geom_boxplot(data=precip.max, aes(x=Model, y=GPP.centered), color="blue", fill="NA")

t.test(abs(precip.min$GPP.centered), abs(precip.max$GPP.centered))
t.test(abs(temp.min$GPP.centered), abs(temp.max$GPP.centered))
