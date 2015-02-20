# Exploratory analyses for the Carbon Paper
setwd("../../..")
fig.path <- "Analyses/Temporal Scaling/Figures"

# -----------------------------------------------------------
# Going from Monthly to Yearly Climate variables
# -----------------------------------------------------------
met.month <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Monthly.csv")
summary(met.month)

vars.met <- names(met.month[,5:ncol(met.month)])

# Creating an object of yearly means
met.year <- aggregate(met.month[,vars.met], by=list(met.month$Site, met.month$Year), FUN=mean)
names(met.year) <- c("Site", "Year", vars.met)
summary(met.year)

# This should be modified based on a priori list of temporally varying env factors that drive ecosystem function
# This step should include the facotrs for the Carbon Paper and Resilience paper
for(s in unique(met.year$Site)){
	for(y in min(met.year$Year):max(met.year$Year)){
		met.year[met.year$Site==s & met.year$Year==y ,"precipf.JJA"] <- mean(met.month[met.month$Site==s & met.month$Year==y & met.month$Month>=6 & met.month$Month<=8, "precipf"], na.rm=T)
		met.year[met.year$Site==s & met.year$Year==y ,"tair.JJA"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=6 & met.month $Month<=8, "tair"])
		met.year[met.year$Site==s & met.year$Year==y ,"tair.MAM"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=3 & met.month $Month<=5, "tair"])
		met.year[met.year$Site==s & met.year$Year==y ,"tair.SON"] <- mean(met.month[met.month$Site==s & met.month $Year==y & met.month $Month>=10 & met.month $Month<=11, "tair"])
	}
}
summary(met.year)


write.csv(met.year, "phase1a_output_variables/PalEON_MIP_Drivers_Year.csv", row.names=F)
# -----------------------------------------------------------


# -----------------------------------------------------------
# Formatting Met & Ecosystem Data
# -----------------------------------------------------------
met.year <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Year.csv")
summary(met.year)
met.col <- ncol(met.year)

ecosys <- read.csv("phase1a_output_variables/PalEON_MIP_AGB_LAI_NPP_Yearly.csv")
summary(ecosys)

# Calculating the percent deviations in the ecosystem variables
for(s in unique(ecosys$Site)){
	# ---------------------------
	# Driver Deviations
	# ---------------------------
	for(v in names(met.year[,3:met.col])){
		tmp <- mean(met.year[met.year$Site==s, v], na.rm=T)
		met.year[met.year$Site==s, paste(v, "dev", sep=".")] <- met.year[met.year$Site==s, v] - tmp
	}

	# ---------------------------
	
	# ---------------------------
	# Ecosystem Deviations
	# ---------------------------
	# for(m in unique(ecosys$Model)){
		# x.npp <- mean(ecosys[ecosys$Site==s & ecosys$Model==m, "NPP"], na.rm=T)
		# x.agb <- mean(ecosys[ecosys$Site==s & ecosys$Model==m, "AGB"], na.rm=T)
		# x.lai <- mean(ecosys[ecosys$Site==s & ecosys$Model==m, "LAI"], na.rm=T)

		# ecosys[ecosys$Site==s & ecosys$Model==m, "NPP.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "NPP"] - x.npp)/x.npp
		# ecosys[ecosys$Site==s & ecosys$Model==m, "NPP.rel"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "NPP"])/x.npp
		# ecosys[ecosys$Site==s & ecosys$Model==m, "AGB.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "AGB"] - x.agb)/x.agb
		# ecosys[ecosys$Site==s & ecosys$Model==m, "AGB.rel"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "AGB"])/x.agb
		# ecosys[ecosys$Site==s & ecosys$Model==m, "LAI.dev"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "LAI"] - x.lai)/x.lai
		# ecosys[ecosys$Site==s & ecosys$Model==m, "LAI.rel"] <- (ecosys[ecosys$Site==s & ecosys$Model==m, "LAI"])/x.lai
	# }
	# ---------------------------
}
write.csv(met.year, "phase1a_output_variables/PalEON_MIP_Drivers_Year.csv", row.names=F)

summary(met.year)
summary(ecosys)


sites <- read.csv("env_drivers/phase1a_env_drivers_v2/PalEON_Phase1a_sites.csv")
row.names(sites) <- sites$code
sites <- sites[,3:ncol(sites)]
summary(sites)

sites2 <- data.frame(t(sites))
sites2$Site <- as.factor(row.names(sites2))
summary(sites2)
dim(sites2)
write.csv(sites2, "env_drivers/phase1a_env_drivers_v2/PalEON_Phase1a_sites2.csv", row.names=F)

sites <- read.csv("env_drivers/phase1a_env_drivers_v2/PalEON_Phase1a_sites2.csv")
summary(sites)

# -----------------------------------
# Reshaping for ordination
environ <- merge(met.year, sites2)
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

for(m in unique(ecosys$Model)){
	for(s in unique(ecosys$Site)){
		# Decadal Smoothing
		for(y in (min(ecosys$Year)+5):(max(ecosys$Year)-5)){
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "NPP.dev.10"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "NPP.dev"] )
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "AGB.dev.10"] <-  mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "AGB.dev"] )
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "LAI.dev.10"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "LAI.dev"] ) 

			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "NPP.rel.10"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "NPP.rel"] )
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "AGB.rel.10"] <-  mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "AGB.rel"] )
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "LAI.rel.10"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-5 & ecosys$Year<=y+5, "LAI.rel"] ) 
		}
		# Centennial Smoothing
		for(y in (min(ecosys$Year)+50):(max(ecosys$Year)-50)){
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "NPP.dev.100"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "NPP.dev"] )
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "AGB.dev.100"] <-  mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "AGB.dev"] )
			# ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "LAI.dev.100"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "LAI.dev"] )  

			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "NPP.rel.100"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "NPP.rel"] )
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "AGB.rel.100"] <-  mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "AGB.rel"] )
			ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year==y, "LAI.rel.100"] <- mean(ecosys[ecosys$Model==m & ecosys$Site==s & ecosys$Year>=y-50 & ecosys$Year<=y+50, "LAI.rel"] )  
		}
	}
}
summary(ecosys)
write.csv(ecosys, "phase1a_output_variables/PalEON_MIP_AGB_LAI_NPP_Yearly_2.csv", row.names=F)

# --------------------------------------------------------
ecosys <- read.csv("phase1a_output_variables/PalEON_MIP_AGB_LAI_NPP_Yearly_2.csv")
summary(ecosys)

library(ggplot2)
library(grid)

large.axes <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=16), axis.text.y=element_text(color="black", size=12), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=24, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

model.colors <- c("black", "cyan4", "hotpink3", "green3", "orange3", "steelblue3")


ggplot(data=ecosys[,]) + large.axes + facet_grid(Site ~ .) +
	geom_line(aes(y=NPP.dev, x=Year, color=Model)) +
	scale_color_manual(values=model.colors) + labs(color="Models") +
 	theme(legend.position=c(0.25,0.9), legend.text=element_text(size=18), legend.title=element_text(size=20), legend.key=element_rect(fill="white"), legend.key.width=unit(2, "line")) + guides(col=guide_legend(ncol=2))

ggplot(data=ecosys) + large.axes + facet_grid(Site ~ .) +
	geom_line(aes(y=NPP.dev.10, x=Year, color=Model)) +
	scale_y_continuous(limits=c(-5,5))

ggplot(data=ecosys) + facet_grid(Site ~ .) +
	geom_line(aes(y=NPP.dev.100, x=Year, color=Model)) +
	scale_y_continuous(limits=range(ecosys$NPP.dev.100, na.rm=T))


summary(ecosys)

# --------------------------------------------------------
# Annual NPP Deviation, PHA
# --------------------------------------------------------
par(mar=c(5,5,4,1))
plot(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",],  range=range(ecosys$NPP.dev, na.rm=T), type="l", col="black", lwd=2, cex.axis=1.5, font.lab=2, cex.lab=2, ylab="% deviation NPP", main="Annual NPP Deviation through time", cex.main=2)
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], lwd=2, col="orange3")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",], lwd=2, col="black")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.wsl",], lwd=1.5, col="blue")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.guess",], lwd=1.5, col="green3")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="clm45",], lwd=1.5, col="red")
legend("topleft", legend=c("ed2", "clm45", "lpj.guess", "lpj.wsl", "jules"), col=c("black", "red", "green3", "blue", "orange3"), lwd=3, bty="l", cex=1.25, bg="white")

par(mar=c(5,5,4,1))
plot(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",],  ylim=range(ecosys[ecosys$Site=="PHA", "NPP.dev"], na.rm=T), type="l", col="black", lwd=3, cex.axis=1.5, font.lab=2, cex.lab=2, ylab="% deviation NPP", main="Annual NPP Deviation 1800-1900", xlim=c(1800, 1900), cex.main=2)
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], lwd=3, col="orange3")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",], lwd=2, col="black")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.wsl",], lwd=3, col="blue")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.guess",], lwd=3, col="green3")
lines(NPP.dev ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="clm45",], lwd=3, col="red")
legend("topleft", legend=c("ed2", "clm45", "lpj.guess", "lpj.wsl", "jules"), col=c("black", "red", "green3", "blue", "orange3"), lwd=3, bty="n", cex=1.5)
# --------------------------------------------------------

# --------------------------------------------------------
# Decadal NPP Deviation, PHA
# --------------------------------------------------------
par(mar=c(5,5,4,1))
plot(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",],  ylim=range(ecosys[ecosys$Site=="PHA", "NPP.dev.10"], na.rm=T), type="l", col="black", lwd=2, cex.axis=1.5, font.lab=2, cex.lab=2, ylab="% deviation NPP", main="NPP Deviation through time", cex.main=2)
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], lwd=2, col="orange3")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",], lwd=2, col="black")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.wsl",], lwd=2, col="blue")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.guess",], lwd=2, col="green3")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="clm45",], lwd=2, col="red")
legend("topleft", legend=c("ed2", "clm45", "lpj.guess", "lpj.wsl", "jules"), col=c("black", "red", "green3", "blue", "orange3"), lwd=3, bty="n")

par(mar=c(5,5,4,1))
plot(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",],  ylim=range(ecosys[ecosys$Site=="PHA", "NPP.dev.10"], na.rm=T), type="l", col="black", lwd=3, cex.axis=1.5, font.lab=2, cex.lab=2, ylab="% deviation NPP", main="Decadal NPP Deviation 1800-1900", xlim=c(1800, 1900), cex.main=2)
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], lwd=3, col="orange3")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",], lwd=2, col="black")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.wsl",], lwd=3, col="blue")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.guess",], lwd=3, col="green3")
lines(NPP.dev.10 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="clm45",], lwd=3, col="red")
legend("topleft", legend=c("ed2", "clm45", "lpj.guess", "lpj.wsl", "jules"), col=c("black", "red", "green3", "blue", "orange3"), lwd=3, bty="n", cex=1.5)
# --------------------------------------------------------

# --------------------------------------------------------
# Centennial NPP Deviation, PHA
# --------------------------------------------------------
par(mar=c(5,5,4,1))
plot(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",],  ylim=range(ecosys[ecosys$Site=="PHA", "NPP.dev.100"], na.rm=T), type="l", col="black", lwd=3, cex.axis=1.5, font.lab=2, cex.lab=2, ylab="% deviation NPP", main="Centennial NPP Deviation through time", cex.main=2)
lines(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="jules.stat",], lwd=3, col="orange3")
lines(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="ed2",], lwd=3, col="black")
lines(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.wsl",], lwd=3, col="blue")
lines(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="lpj.guess",], lwd=3, col="green3")
lines(NPP.dev.100 ~ Year, data=ecosys[ecosys$Site=="PHA" & ecosys$Model=="clm45",], lwd=3, col="red")
legend("topleft", legend=c("ed2", "clm45", "lpj.guess", "lpj.wsl", "jules"), col=c("black", "red", "green3", "blue", "orange3"), lwd=3, bty="n", cex=1.5)
# --------------------------------------------------------

# -----------------------------------------------------------
# Graphing responses averaged over sites
# -----------------------------------------------------------
# -----------------------------------------------------------

# -----------------------------------------------------------
# Driver sensitivities Fig 4
#	NOTE: these currently have done nothing about temporal autocorrelation
# -----------------------------------------------------------
sensitivity <- stack(ecosys[,which(substr(names(ecosys), 4, 7)==".rel")])
names(sensitivity) <- c("percent.mean", "var.scale")
sensitivity[,3:5] <- ecosys[,1:3]
sensitivity$var <- as.factor(substr(sensitivity$var.scale, 1, 3))
sensitivity$scale <- as.factor(gsub(".*rel", "t", sensitivity$var.scale))
levels(sensitivity$scale) <- c("annual", "decadal", "centennial")
summary(sensitivity)

summary(met.year)
sensitivity <- merge(sensitivity, met.year)
summary(sensitivity)

model.colors <- c("black", "cyan4", "hotpink3", "green3", "orange3", "steelblue3")

pdf(file.path(fig.path, "Sensitivity_Tair_AllSites.pdf"))
ggplot(aes(x=tair.dev, y=percent.mean, color=Model), data=sensitivity[,]) + 
	facet_grid(scale ~ var) +
	geom_point(size=1) +
	geom_smooth(method=lm, alpha=0.5, size=1) + 
	scale_color_manual(values=model.colors) + 
	theme_bw()
dev.off()

pdf(file.path(fig.path, "Sensitivity_Tair_NPP_AllSites.pdf"))
ggplot(aes(x=tair.dev, y=percent.mean, color=Model), data=sensitivity[sensitivity$var=="NPP",]) + 
	facet_grid(scale ~ Site) +
	geom_point(size=1) +
	geom_smooth(method=lm, alpha=0.5, size=1) + 
	scale_color_manual(values=model.colors) + 
	theme_bw()
dev.off()


pdf(file.path(fig.path, "Sensitivity_Tair_PHA.pdf"))
ggplot(aes(x=tair.dev, y=percent.mean, color=Model), data=sensitivity[sensitivity$Site=="PHA",]) + 
	facet_grid(scale ~ var) +
	geom_point(size=1) +
	geom_smooth(method=lm, alpha=0.5, size=1) + 
	scale_color_manual(values=model.colors) + 
	theme_bw()
dev.off()

# -----------------------------------------------------------






# -----------------------------------------------------------
# -----------------------------------------------------------

# -----------------------------------------------------------
# Comparison of Model NPP deviations (Annual)
# 1. Model correlation matrix
# 2. Oridination comparion responses & explanatory variables
#	## NOTE: THIS DOESN'T WORK LIKE I EXPECTED
#	## NEW IDEA: ordination of driver sensitivities: % of mean per change in Temp, Precip...
# -----------------------------------------------------------
ecosys <- read.csv("phase1a_output_variables/PalEON_MIP_AGB_LAI_NPP_Yearly_2.csv")
summary(ecosys)

met.year <- read.csv("phase1a_output_variables/PalEON_MIP_Drivers_Year.csv")
summary(met.year)


library(reshape2)
ecosysb <- ecosys
ecosysb$Year <- as.factor(ecosysb$Year)

npp.1 <- recast(ecosysb[,c("Year", "Site", "Model", "NPP.dev")], Site+Year~Model)
summary(npp.1)
npp.10 <- recast(ecosysb[,c("Year", "Site", "Model", "NPP.dev.10")], Site+Year~Model)
summary(npp.10)
npp.100 <- recast(ecosysb[,c("Year", "Site", "Model", "NPP.dev.100")], Site+Year~Model)
summary(npp.100)

agb.1 <- recast(ecosysb[,c("Year", "Site", "Model", "AGB.dev")], Site+Year~Model)
agb.1 <- agb.1[,c(1:4, 6:7)]
summary(agb.1)
agb.10 <- recast(ecosysb[,c("Year", "Site", "Model", "AGB.dev.10")], Site+Year~Model)
agb.10 <- agb.10[,c(1:4, 6:7)]
summary(agb.10)
agb.100 <- recast(ecosysb[,c("Year", "Site", "Model", "AGB.dev.100")], Site+Year~Model)
agb.100 <- agb.100[,c(1:4, 6:7)]
summary(agb.100)


summary(npp.models[,3:ncol(npp.models)])
summary(met.year)

#---------------------------
# NPP
#---------------------------
npp.1b <- merge(npp.1, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(npp.1b)
npp.10b <- merge(npp.10, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(npp.10b)
npp.100b <- merge(npp.100, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(npp.100b)

npp.corr.1 <- cor(npp.1b[,c(3:ncol(npp.1b))], use="complete.obs") 
npp.corr.1[,1:5]
npp.corr.10 <- cor(npp.10b[,c(3:ncol(npp.10b))], use="complete.obs") 
npp.corr.10[,1:5]
npp.corr.100 <- cor(npp.100b[,c(3:ncol(npp.100b))], use="complete.obs") 
npp.corr.100[,1:5]
#---------------------------


#---------------------------
# AGB
#---------------------------
agb.1b <- merge(agb.1, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(agb.1b)
agb.10b <- merge(agb.10, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(agb.10b)
agb.100b <- merge(agb.100, met.year[,c(1:2, 15:ncol(met.year))], all.x=T, all.y=T)
summary(agb.100b)

agb.corr.1 <- cor(agb.1b[,c(3:ncol(agb.1b))], use="complete.obs") 
agb.corr.1[,1:4]
agb.corr.10 <- cor(agb.10b[,c(3:ncol(agb.10b))], use="complete.obs") 
agb.corr.10[,1:4]
agb.corr.100 <- cor(agb.100b[,c(3:ncol(agb.100b))], use="complete.obs") 
agb.corr.100[,1:4]
#---------------------------

#---------------------------
# AGB vs. NPP
#---------------------------
temp <- cor(agb.1[,3:ncol(agb.1)], npp.1[,c(3:4,6:7)])
temp
#---------------------------


#---------------------------
library(vegan)

npp.models2 <- recast(ecosysb[,c("Year", "Site", "Model", "NPP.rel")], Site+Year~Model)
summary(npp.models2)

test <- npp.models2[npp.models2$Site=="PHA",3:ncol(npp.models2)]
row.names(test) <- npp.models2[npp.models$Site=="PHA","Year"]
test <- test[complete.cases(test) & test$jules.stat>=0,]
summary(test)
dim(test)

npp.ann <- metaMDS(test, distance="euclidean", wascores=T)
npp.ann

npp.ord <- ordiplot(npp.ann, type="none")
points(npp.ord, "sites", pch=19, cex=0.2)
text(npp.ord, "species", col="blue", cex=1.5, font=2)

summary(environ)
environ2 <- environ[environ$Site=="PHA" & environ$Year %in% row.names(test),]
npp.env <- envfit(npp.ann, environ2[,c(1, 15:ncol(environ2))])
npp.env$vectors

points(npp.ord, "sites", pch=19, cex=0.2)
text(npp.ord, "species", col="blue", cex=1.5, font=2)
plot(npp.env, p.max=0.01, col="red")
# -----------------------------------------------------------

