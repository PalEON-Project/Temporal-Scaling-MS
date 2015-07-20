# ----------------------------------------
# Temporal Scaling Analyses
# Looking at changes in PalEON driver correllation & distribution with temporal scales
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: How are the relationships among PalEON met drivers affected by temporal scale (resolution + extent)
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Raw Met drivers at annual scale (exist in Data/EcosysData.Rdata)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) compare changes in driver distribution from changes in scale 
#    1) 
# B) analyze changes in correlations among drivers with changes in scale
#    1) Create correlation matrices for each scale
#    2) to test statistical change in correlation: corr[pair] ~ scale
#    3) 
# -------------------------
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(nlme)
library(car)
library(zoo)
# library(mvtnorm)
# library(MCMCpack)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
# setwd("~/Desktop/Research/PalEON CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("..")
path.data <- "Data"
dat.dir <- "Data/analysis_drivers"
fig.dir <- "Figures/analysis_drivers"

if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------
load(file.path(path.data, "EcosysData.Rdata")) # loads ecosys & model.colors
summary(ecosys)

# subset just the met data (lets just use the drivers returned by JULES since it had the fewest issues)
vars.met <- c("tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind", "CO2")
drivers <- ecosys[ecosys$Model=="jules.stat", c("Year", "Site", "Scale", vars.met)]
summary(drivers)
# ----------------------------------------

# ----------------------------------------
# Look at driver distributions
# ----------------------------------------
# Stack the drivers together to make them easy to deal with in ggplot
drivers.stack <- stack(drivers[,vars.met])[,c(2,1)]
names(drivers.stack) <- c("driver", "value")
drivers.stack$Year   <- drivers$Year
drivers.stack$Site   <- drivers$Site
drivers.stack$Scale  <- drivers$Scale
drivers.stack$Extent <- as.factor("0850-2010")
summary(drivers.stack)

# ---------------
# Finding some summary statistics by temporal resolution
# ---------------
mean.res <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=mean, na.rm=T)
names(mean.res) <- c("driver", "Scale", "mean")
mean.res$sd     <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=sd, na.rm=T)[,3]
mean.res$median <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=median, na.rm=T)[,3]
mean.res$min    <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=min, na.rm=T)[,3]
mean.res$max    <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=max, na.rm=T)[,3]
mean.res$p.025  <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=quantile, 0.025, na.rm=T)[,3]
mean.res$p.975  <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Scale), FUN=quantile, 0.975, na.rm=T)[,3]
mean.res <- mean.res[order(mean.res$driver),]
mean.res

write.csv(mean.res, file.path(dat.dir, "Drivers_Summary_byResolution.csv"), row.names=F)
# ---------------

# ---------------
# Finding some summary statistics by temporal extent
# ---------------
extents <- data.frame(Start=c(0850, 1900, 1990), End=c(2010, 2010, 2010)) 

mean.extent <- merge(unique(drivers.stack$driver), extents)
names(mean.extent)[1] <- "driver"
mean.extent <- mean.extent[order(mean.extent$driver),]
mean.extent$Extent <- paste(ifelse(mean.extent$Start==850, paste0("0", mean.extent$Start), mean.extent$Start), mean.extent$End, sep="-")
mean.extent

for(d in unique(mean.extent$driver)){
	for(e in unique(mean.extent$Start)){
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "mean"]   <- mean(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "sd"]     <- sd(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "median"] <- median(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "min"]    <- min(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "max"]    <- max(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "p.025"]  <- quantile(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], 0.025, na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "p.0975"] <- quantile(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Scale=="t.001","value"], 0.975, na.rm=T)
	}
}
mean.extent
write.csv(mean.res, file.path(dat.dir, "Drivers_Summary_byExtent.csv"), row.names=F)

# ---------------

# ---------------
# Adding in the other extents
# ---------------
drivers.stack2 <- drivers.stack[drivers.stack$Year>=1900,]
drivers.stack2$Extent <- as.factor("1900-2010")
drivers.stack3 <- drivers.stack[drivers.stack$Year>=1990,]
drivers.stack3$Extent <- as.factor("1990-2010")

drivers.stack <- rbind(drivers.stack, drivers.stack2, drivers.stack3)
drivers.stack$n <- length(unique(drivers$Site)) * (2010-as.numeric(substr(drivers.stack$Extent,1,4))+1)
summary(drivers.stack)
# ---------------


# n <- length(unique(drivers$Year))*length(unique(drivers$Site))
# Make some simple plots of the different distributions
pdf(file.path(fig.dir, "Driver_Distributions_by_Scale.pdf"), width=8.5, height=11)
# By Temporal resolution
print(
ggplot(data=drivers.stack[drivers.stack$Extent=="0850-2010",]) + facet_grid(Scale ~ driver, scales="free") +
	geom_histogram(aes(x=value, fill=Site, weight=1/n)) + 
	geom_vline(data=mean.res, aes(xintercept=mean), size=0.5) +
	# geom_vline(xintercept=3, size=2) +
	labs(y="p.obs",title="Distribution of Met Drivers by Temporal Resolution, 0850-2010") +
	theme_bw()
)

# By Temporal Extent
print(
ggplot(data=drivers.stack[drivers.stack$Scale=="t.001",]) + facet_grid(Extent ~ driver, scales="free") +
	geom_histogram(aes(x=value, fill=Site, weight=1/n)) + 
	geom_vline(data=mean.extent, aes(xintercept=mean), size=0.5) +
	# geom_vline(xintercept=3, size=2) +
	labs(y="p.obs",title="Distribution of Met Drivers by Temporal Extent, Annual") +
	theme_bw()
)
dev.off()

# ---------------
# Looking for Extremes or anomalies by site
# ---------------
drivers.basic <- drivers.stack[drivers.stack$Extent=="0850-2010" & drivers.stack$Scale=="t.001",]
summary(drivers.basic)

# recoding the drivers to graph in a sensible order
drivers.basic$driver <- recode(drivers.basic$driver, "'CO2'='1'; 'tair'='2'; 'precipf'='3'; 'swdown'='4'; 'lwdown'='5'; 'qair'='6'; 'psurf'='7'; 'wind'='8'")
levels(drivers.basic$driver) <- c("CO2", "tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind")
summary(drivers.basic)

mean.site <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=mean, na.rm=T)
names(mean.site) <- c("driver", "Site", "mean")
mean.site$sd     <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=sd, na.rm=T)[,3]
mean.site$median <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=median, na.rm=T)[,3]
mean.site$min    <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=min, na.rm=T)[,3]
mean.site$max    <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=max, na.rm=T)[,3]
mean.site$p.025  <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=quantile, 0.025, na.rm=T)[,3]
mean.site$p.975  <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=quantile, 0.975, na.rm=T)[,3]
mean.site$p.005  <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=quantile, 0.005, na.rm=T)[,3]
mean.site$p.995  <- aggregate(drivers.basic$value, by=list(drivers.basic$driver, drivers.basic$Site), FUN=quantile, 0.995, na.rm=T)[,3]
mean.site <- mean.site[order(mean.site$driver),]
mean.site
write.csv(mean.site, file.path(dat.dir, "Drivers_Summary_bySite.csv"), row.names=F)

# Making a dummy dataframe to graph 95% CI 
mean.site2 <- merge(mean.site, data.frame(Year=850:2010))

pdf(file.path(fig.dir, "Drivers_Time_by_Site.pdf"), width=8.5, height=11)
ggplot(data=drivers.basic) + facet_grid(driver~Site, scales="free") +
	geom_line(aes(x=Year, y=value, color=driver), size=1.25) +
	geom_ribbon(data=mean.site2, aes(x=Year, ymin=p.005, ymax=p.995), alpha=0.5) +
	geom_line(data=drivers.stack[drivers.stack$Extent=="0850-2010" & drivers.stack$Scale=="t.100",], aes(x=Year, y=value, color=driver), size=1.25) +
	scale_x_continuous(breaks=c(1250,1750)) +
	scale_color_manual(values=c("green3", "red2", "blue2", "goldenrod3", "goldenrod4", "purple3", "navajowhite4", "salmon4")) +
	guides(color=F) +
	ggtitle("Annual Met Drivers with Centennial Smoothing and Shaded 99% Range") +
	theme_bw()
dev.off()
# ---------------

# ----------------------------------------

# ----------------------------------------
# Look at driver correlations
# ----------------------------------------
# Make some simple plots of correlation
# Note: because of the number of data points, these need to be saved in a "flattened" format
for(s in unique(drivers$Scale)){
	png(file.path(fig.dir, paste0("Driver_Correlations_",s,"0850-2010.png")))
	plot(drivers[drivers$Scale==s, vars.met], main=paste0("Met Driver Correlations 0850-2010, ", s))
	dev.off()
}

# run the correlation matrices
scales <- unique(drivers$Scale)
vars.met

# Need to replae vars.met with something that lists the relevant pairs without doing the redundnat pairs
for(s in 1:length(scales)){
	data.temp <- drivers[drivers$Scale==scales[s],]

	# to get proper pairing we need each variable matched with the ones we haven't done yet so
	for(i in 1:(length(vars.met)-1)){
		for(j in (i+1):length(vars.met)){
			# putting the correlation info into a data frame
			cor.temp <- data.frame(pair   = paste(vars.met[i], vars.met[j], sep=":"), 
						 v1     = vars.met[i], v2 = vars.met[j],
						 scale  = scales[s],
						 cor.est= cor(data.temp[,vars.met[c(i,j)]], use="pairwise.complete.obs")[1,2],
						 cor.lo = cor.test(data.temp[,vars.met[i]], data.temp[,vars.met[j]], use="pairwise.complete.obs")$conf.int[1],
						 cor.hi = cor.test(data.temp[,vars.met[i]], data.temp[,vars.met[j]], use="pairwise.complete.obs")$conf.int[2]
						 )
			# If this is the first run through, core.scales is that dataframe, otherwise, rbind it
			if(s==1 & i==1 & j==2){ 
				cor.res <- cor.temp 
			} else {
				cor.res <- rbind(cor.res, cor.temp)
			}
		}
	}
}
# Making a sign color to make the graphs; note there's some tricks to getting everything to graph okay
cor.res$sign <- as.factor(ifelse(cor.res$cor.est>=0 | (cor.res$cor.hi>=0 & !cor.res$scale=="t.250"), "pos", "neg"))
summary(cor.res)

# graph the changes in correlation
# Note: graphing the absolute value so that going down is always a decrease in correlation
pdf(file.path(fig.dir, "Driver_Correlations_Across_Scale.pdf"))
ggplot(data=cor.res[,]) + facet_grid(v2 ~ v1) +
	# geom_ribbon(aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi)), fill="red2", alpha=0.7) +
	geom_hline(data=cor.res[cor.res$scale=="t.001",], aes(yintercept=abs(cor.est)), linetype="dashed", size=0.25) +
	geom_ribbon(data=cor.res[cor.res$sign=="neg",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	geom_ribbon(data=cor.res[cor.res$sign=="pos",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	# geom_line(aes(x=as.numeric(scale), y=cor.est), size=0.25) +
	scale_fill_manual(values=c("blue2", "red2")) +
	guides(fill=guide_legend(title="Corr. Sign")) +
	labs(x="Scale", y="Absolute Value of Correlation (R2)", title="Change in Driver Correlation Across Temporal Resolution") +
	theme_bw()
dev.off()
# ----------------------------------------
