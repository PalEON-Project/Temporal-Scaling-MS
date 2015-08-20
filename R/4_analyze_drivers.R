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
in.base <- "Data/gamms"
in.res  <- "Big4_GS_byResolution"
in.ext  <- "Big4_GS_byExtent"
dat.dir <- "Data/analysis_drivers_big4"
fig.dir <- "Figures/analysis_drivers_big4"

if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# Load data files & function scripts
# ----------------------------------------

load(file.path(path.data, "EcosysData_Raw.Rdata"))
ecosys <- ecosys[!ecosys$Model=="clm.bgc" & !ecosys$Model=="linkages",]

# Figure out what models we have to work with
#response.all <- c("NPP", "AGB.diff", "NEE")
response.all <- c("NPP")
models <- unique(ecosys$Model)

f.res <- dir(file.path(in.base, in.res))
f.ext <- dir(file.path(in.base, in.ext))

# Need to recode the normal ed so that it will only return one model
models2 <- recode(models, "'ed2'='ed2_'")


# -----------------------
# Put all responses & models & scales into single data frames
# This is the same loop as in other scripts, just parsed down for clarity
# -----------------------
# for(r in 1:length(response.all)){
response="NPP"
response.res <- grep(response, f.res)
response.ext <- grep(response, f.ext)
for(i in 1:length(models)){
    # First narrow to the models
	fmod <- response.res[grep(models2[i], f.res[response.res])]

    if(!length(fmod)>0) next
    load(file.path(in.base, in.res, f.res[fmod]))
  
    if(i==1) {
      dat.ecosys <- cbind(mod.out$data[,], mod.out$ci.response[,c("mean", "lwr", "upr")])
    } else {
      dat.ecosys <- rbind(dat.ecosys, cbind(mod.out$data, mod.out$ci.response[,c("mean", "lwr", "upr")]))
    }
    
    # loop through by extent
    fmod <- response.ext[grep(models2[i], f.ext[response.ext])]
    load(file.path(in.base, in.ext, f.ext[fmod]))
    
    # Note: because we're lumping everything together, let's not mess with reiterating the base level
    dat.ecosys <- rbind(dat.ecosys,
                        cbind(mod.out$data[!(mod.out$data$Resolution=="t.001" & substr(mod.out$data$Extent,1,3)=="850"),], 
                              mod.out$ci.response[!(mod.out$ci.response$Resolution=="t.001" & substr(mod.out$ci.response$Extent,1,3)=="850"),c("mean", "lwr", "upr")]))
    
    # Clear the mod.out to save space
    rm(mod.out)
}
  
# Fix Extent Labels for consistency
dat.ecosys$Extent <- as.factor(ifelse(dat.ecosys$Extent=="850-2010", "0850-2010", paste(dat.ecosys$Extent)))
summary(dat.ecosys)

# subset just the met data (lets just use the drivers returned by JULES since it had the fewest issues)
vars.met <- c("tair", "precipf", "swdown", "CO2")
drivers <- dat.ecosys[dat.ecosys$Model=="ed2", c("Year", "Site", "Resolution", "Extent", vars.met)]
drivers$Met.Source <- as.factor(ifelse(drivers$Year<1850, "CCSM4_p1000", ifelse(drivers$Year>=1901, "CRUNCEP", "CCSM4_historical")))
summary(drivers)


# ----------------------------------------
# Plotting Drivers through time at all sites & where the data came from
# ----------------------------------------
drivers.stack <- stack(drivers[,vars.met])[c(2,1)]
names(drivers.stack) <- c("Driver", "Value")
drivers.stack$Year   <- drivers$Year
drivers.stack$Site   <- drivers$Site
drivers.stack$Extent <- drivers$Extent
drivers.stack$Resolution <- drivers$Resolution
drivers.stack$Source <- drivers$Met.Source
summary(drivers.stack)

# Adding in some dummy max/min by driver to make some prettier graphs
for(v in vars.met){
	drivers.stack[drivers.stack$Driver==v, "Val.min"] <- min(drivers.stack[drivers.stack$Driver==v, "Value"])
	drivers.stack[drivers.stack$Driver==v, "Val.max"] <- max(drivers.stack[drivers.stack$Driver==v, "Value"])
}
drivers.stack[drivers.stack$Driver=="CO2","Source"] <- NA
summary(drivers.stack)

# Plotting by driver
pdf(file.path(fig.dir, "MetDrivers_byDriver_Source.pdf"))
for(d in unique(drivers.stack$Driver)){
print(
ggplot(data=drivers.stack[drivers.stack$Driver==d & drivers.stack$Resolution=="t.001",]) + facet_wrap(~Site) +
	geom_ribbon(aes(x=Year, ymin=Val.min, ymax=Val.max, fill=Source), alpha=0.5) +	
	geom_line(aes(x=Year, y=Value)) +
	geom_line(data=drivers.stack[drivers.stack$Driver==d & drivers.stack$Resolution=="t.050",], aes(x=Year, y=Value), color="gray50", size=1.5) +
	geom_point(data=drivers.stack[drivers.stack$Driver==d & drivers.stack$Resolution=="t.050",], aes(x=Year, y=Value), color="gray50", size=2) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	labs(title=d) +
	theme(plot.title=element_text(face="bold", size=rel(2))) + 
	theme(legend.position="top") +
	# theme(legend.position=c(0.85,0.2)) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()
	      )
)
}
dev.off()

pdf(file.path(fig.dir, "MetDrivers_bySite_Source.pdf"))
for(s in unique(drivers.stack$Site)){
print(
ggplot(data=drivers.stack[drivers.stack$Site==s,]) + facet_wrap(~Driver, scales="free_y") +
	geom_ribbon(aes(x=Year, ymin=Val.min, ymax=Val.max, fill=Source), alpha=0.5) +	
	geom_line(aes(x=Year, y=Value)) +
	geom_line(data=drivers.stack[drivers.stack$Site==s & drivers.stack$Resolution=="t.050",], aes(x=Year, y=Value), color="gray50", size=1.5) +
	geom_point(data=drivers.stack[drivers.stack$Site==s & drivers.stack$Resolution=="t.050",], aes(x=Year, y=Value), color="gray50", size=2) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	labs(title=s, y="Value") +
	theme(plot.title=element_text(face="bold", size=rel(2))) + 
	theme(legend.position="top") +
	# theme(legend.position=c(0.85,0.2)) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()
	      )
)
}
dev.off()

for(d in unique(drivers.stack$Driver)){
pdf(file.path(fig.dir, paste0("MetDrivers_", d, "_Source.pdf")))
print(
ggplot(data=drivers.stack[drivers.stack$Driver==d,]) + facet_wrap(~Site) +
	geom_ribbon(aes(x=Year, ymin=Val.min, ymax=Val.max, fill=Source), alpha=0.5) +	
	geom_line(aes(x=Year, y=Value)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	labs(title=d) +
	theme(plot.title=element_text(face="bold", size=rel(2))) + 
	theme(legend.position="top") +
	# theme(legend.position=c(0.85,0.2)) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank()
	      )
)
dev.off()
}
dev.off()

# ------------------
# Plotting just PHA for talks, etc
# ------------------
summary(drivers.stack)

# Subset PHA & rename/reorder drivers
drivers.pha <- drivers.stack[drivers.stack$Site=="PHA",]
drivers.pha$Driver.Order <- recode(drivers.pha$Driver, "'tair'='1'; 'precipf'='2'; 'swdown'='3'; 'CO2'='4'")
levels(drivers.pha$Driver.Order) <- c("Temp", "Precip", "SW Rad", "CO2")

# Conver tair into celcius
drivers.pha[drivers.pha$Driver=="tair","Value"] <- drivers.pha[drivers.pha$Driver=="tair","Value"]-273.15


pdf(file.path(fig.dir, "Drivers_Big4_0850-2010_ggplot_raw.pdf"), width=10, height=7.5)
ggplot(data=drivers.pha) + facet_grid(Driver.Order~., scales="free_y") +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.001" & drivers.pha$Extent=="0850-2010",], 
			  aes(x=Year, y=Value, color=Driver.Order), alpha=0.5) +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.050" & drivers.pha$Extent=="0850-2010",], 
			  aes(x=Year, y=Value, color=Driver.Order), size=1.5) +
	geom_point(data=drivers.pha[drivers.pha$Resolution=="t.050" & drivers.pha$Extent=="0850-2010",], 
			   aes(x=Year, y=Value, color=Driver.Order), size=3) +
	scale_y_continuous(name="", expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0), breaks=c(1000, 1250, 1500, 1750, 2000)) +
	scale_color_manual(values=c("red", "blue", "darkgoldenrod2", "green2")) +
	guides(color=F) +
	theme(plot.title=element_text(face="bold", size=rel(3))) + 
	theme(strip.text=element_text(size=rel(1.5))) +
	theme(legend.text=element_text(size=rel(1.75)), 
	      legend.title=element_text(size=rel(2)),
	      legend.key=element_blank(),
	      legend.key.size=unit(2, "lines")) + 
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
dev.off()

# customizing a plot for each driver so thave independent y-axes

plot.temp <- ggplot() + 
	facet_grid(Driver.Order~., scales="free_y") +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.001" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="tair",], 
			  aes(x=Year, y=Value, color=Driver.Order), alpha=0.5) +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="tair",], 
			  aes(x=Year, y=Value, color=Driver.Order), size=1.5) +
	geom_point(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                            drivers.pha$Extent=="0850-2010" & 
	                            drivers.pha$Driver=="tair",], 
			   aes(x=Year, y=Value, color=Driver.Order), size=3) +
	geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_y_continuous(name=expression(paste("Temp ("^"o","C)")), expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0), breaks=c(1000, 1250, 1500, 1750, 2000)) +
	scale_color_manual(values=c("red")) +
	guides(color=F) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_blank(), 
	      axis.text.y=element_text(color="black", size=rel(1.25)), 
	      axis.title.x=element_blank(),  
	      axis.title.y=element_text(face="bold", size=rel(1.25), vjust=2),
	      axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1, "lines")) +
	theme(plot.margin=unit(c(1,1,-1,0.9), "lines"))

plot.precip <- ggplot() + 
	facet_grid(Driver.Order~., scales="free_y") +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.001" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="precipf",], 
			  aes(x=Year, y=Value, color=Driver.Order), alpha=0.5) +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="precipf",], 
			  aes(x=Year, y=Value, color=Driver.Order), size=1.5) +
	geom_point(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                            drivers.pha$Extent=="0850-2010" & 
	                            drivers.pha$Driver=="precipf",], 
			   aes(x=Year, y=Value, color=Driver.Order), size=3) +
	geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_y_continuous(name=expression(paste("Precip (mm yr"^"-1",")")), expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0), breaks=c(1000, 1250, 1500, 1750, 2000)) +
	scale_color_manual(values=c("blue")) +
	guides(color=F) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_blank(), 
	      axis.text.y=element_text(color="black", size=rel(1)), 
	      axis.title.x=element_blank(),  
	      axis.title.y=element_text(face="bold", size=rel(1.25), vjust=0.6),
	      axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1, "lines")) +
	theme(plot.margin=unit(c(0,1,-1,0.15), "lines"))

plot.swdown <- ggplot() + 
	facet_grid(Driver.Order~., scales="free_y") +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.001" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="swdown",], 
			  aes(x=Year, y=Value, color=Driver.Order), alpha=0.5) +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="swdown",], 
			  aes(x=Year, y=Value, color=Driver.Order), size=1.5) +
	geom_point(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                            drivers.pha$Extent=="0850-2010" & 
	                            drivers.pha$Driver=="swdown",], 
			   aes(x=Year, y=Value, color=Driver.Order), size=3) +
	geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_y_continuous(name=expression(paste("Rad (W m"^"-2",")")), expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0), breaks=c(1000, 1250, 1500, 1750, 2000)) +
	scale_color_manual(values=c("darkgoldenrod3")) +
	guides(color=F) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_blank(), 
	      axis.text.y=element_text(color="black", size=rel(1.25)), 
	      axis.title.x=element_blank(),  
	      axis.title.y=element_text(face="bold", size=rel(1.25), vjust=1),
	      axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1, "lines")) +
	theme(plot.margin=unit(c(0,1,-1,0.3), "lines"))

plot.co2 <- ggplot() + 
	facet_grid(Driver.Order~., scales="free_y") +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.001" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="CO2",], 
			  aes(x=Year, y=Value, color=Driver.Order), alpha=0.5) +
	geom_line(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                           drivers.pha$Extent=="0850-2010" & 
	                           drivers.pha$Driver=="CO2",], 
			  aes(x=Year, y=Value, color=Driver.Order), size=1.5) +
	geom_point(data=drivers.pha[drivers.pha$Resolution=="t.050" & 
	                            drivers.pha$Extent=="0850-2010" & 
	                            drivers.pha$Driver=="CO2",], 
			   aes(x=Year, y=Value, color=Driver.Order), size=3) +
	geom_vline(xintercept=1901, linetype="dashed", size=1.5) +
	scale_y_continuous(name=expression(paste("CO2 (ppm)")), expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0), breaks=c(1000, 1250, 1500, 1750, 2000)) +
	scale_color_manual(values=c("green4")) +
	guides(color=F) +
	theme(strip.text=element_text(size=rel(1.5), face="bold")) +
	theme(axis.line=element_line(color="black", size=0.5), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(2.5)), 
	      axis.text.y=element_text(color="black", size=rel(1.25)), 
	      axis.title.x=element_text(face="bold", angle=0, color="black", size=rel(2.5)),  
	      axis.title.y=element_text(face="bold", size=rel(1.25), vjust=1),
	      axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1, "lines")) +
	theme(plot.margin=unit(c(0,1,1,0.56), "lines"))



pdf(file.path(fig.dir, "Drivers_Big4_0850-2010.pdf"), width=10, height=7.5)
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,1)))
print(plot.temp,   vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.precip, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.swdown, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot.co2,    vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
dev.off()


# ------------------


# ----------------------------------------

# ----------------------------------------
# Look at driver distributions
# ----------------------------------------
# Stack the drivers together to make them easy to deal with in ggplot
drivers.stack <- stack(drivers[,vars.met])[,c(2,1)]
names(drivers.stack) <- c("driver", "value")
drivers.stack$Year   <- drivers$Year
drivers.stack$Site   <- drivers$Site
drivers.stack$Resolution  <- drivers$Resolution
drivers.stack$Extent <- as.factor("0850-2010")
summary(drivers.stack)

# ---------------
# Finding some summary statistics by temporal resolution
# ---------------
mean.res <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=mean, na.rm=T)
names(mean.res) <- c("driver", "Resolution", "mean")
mean.res$sd     <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=sd, na.rm=T)[,3]
mean.res$median <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=median, na.rm=T)[,3]
mean.res$min    <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=min, na.rm=T)[,3]
mean.res$max    <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=max, na.rm=T)[,3]
mean.res$p.025  <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=quantile, 0.025, na.rm=T)[,3]
mean.res$p.975  <- aggregate(drivers.stack$value, by=list(drivers.stack$driver, drivers.stack$Resolution), FUN=quantile, 0.975, na.rm=T)[,3]
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
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "mean"]   <- mean(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "sd"]     <- sd(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "median"] <- median(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "min"]    <- min(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "max"]    <- max(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "p.025"]  <- quantile(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], 0.025, na.rm=T)
		mean.extent[mean.extent$driver==d & mean.extent$Start==e, "p.0975"] <- quantile(drivers.stack[drivers.stack$driver==d & drivers.stack$Year>=e & drivers.stack$Resolution=="t.001","value"], 0.975, na.rm=T)
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
pdf(file.path(fig.dir, "Driver_Distributions_by_Resolution.pdf"), width=8.5, height=11)
# By Temporal resolution
print(
ggplot(data=drivers.stack[drivers.stack$Extent=="0850-2010",]) + facet_grid(Resolution ~ driver, scales="free") +
	geom_histogram(aes(x=value, fill=Site, weight=1/n)) + 
	geom_vline(data=mean.res, aes(xintercept=mean), size=0.5) +
	# geom_vline(xintercept=3, size=2) +
	labs(y="p.obs",title="Distribution of Met Drivers by Temporal Resolution, 0850-2010") +
	theme_bw()
)

# By Temporal Extent
print(
ggplot(data=drivers.stack[drivers.stack$Resolution=="t.001",]) + facet_grid(Extent ~ driver, scales="free") +
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
drivers.basic <- drivers.stack[drivers.stack$Extent=="0850-2010" & drivers.stack$Resolution=="t.001",]
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
	geom_line(data=drivers.stack[drivers.stack$Extent=="0850-2010" & drivers.stack$Resolution=="t.100",], aes(x=Year, y=value, color=driver), size=1.25) +
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
for(s in unique(drivers$Resolution)){
	png(file.path(fig.dir, paste0("Driver_Correlations_",s,"0850-2010.png")))
	plot(drivers[drivers$Resolution==s, vars.met], main=paste0("Met Driver Correlations 0850-2010, ", s))
	dev.off()
}

# run the correlation matrices
scales <- unique(drivers$Resolution)
vars.met

# Need to replae vars.met with something that lists the relevant pairs without doing the redundnat pairs
for(s in 1:length(scales)){
	data.temp <- drivers[drivers$Resolution==scales[s],]

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
pdf(file.path(fig.dir, "Driver_Correlations_Across_Resolution.pdf"))
ggplot(data=cor.res[,]) + facet_grid(v2 ~ v1) +
	# geom_ribbon(aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi)), fill="red2", alpha=0.7) +
	geom_hline(data=cor.res[cor.res$scale=="t.001",], aes(yintercept=abs(cor.est)), linetype="dashed", size=0.25) +
	geom_ribbon(data=cor.res[cor.res$sign=="neg",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	geom_ribbon(data=cor.res[cor.res$sign=="pos",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	# geom_line(aes(x=as.numeric(scale), y=cor.est), size=0.25) +
	scale_fill_manual(values=c("blue2", "red2")) +
	guides(fill=guide_legend(title="Corr. Sign")) +
	labs(x="Resolution", y="Absolute Value of Correlation (R2)", title="Change in Driver Correlation Across Temporal Resolution") +
	theme_bw()
dev.off()
# ----------------------------------------
