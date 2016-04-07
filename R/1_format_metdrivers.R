# ----------------------------------------
# Objective: Pull & format the raw met drivers for use in downstream analyses
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Workflow
# -------------------------
# 1) Extract raw met drivers at MONTHLY scale -- extract from ED & SiBCASA
# 2) Reorganize it a bit
# 3) Aggregate it into annual & growing season means
# 4) Some simple stats comparing annual vs. growing season met
# -------------------------
#
# ----------------------------------------

# ----------------------------------------
# Load Libaries
# ----------------------------------------
library(ggplot2); library(grid)
library(nlme)
library(car)
library(zoo)
library(ncdf4)
# library(mvtnorm)
# library(MCMCpack)
# ----------------------------------------

# ----------------------------------------
# Set Directories
# ----------------------------------------
# setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
path.data <- "Data"
model.dir <- "~/Dropbox/PalEON_CR/paleon_mip_site/phase1a_model_output"
dat.dir <- "Data/analysis_met_drivers"
fig.dir <- "Figures/analysis_met_drivers"

if(!dir.exists(dat.dir)) dir.create(dat.dir)
if(!dir.exists(fig.dir)) dir.create(fig.dir)
# ----------------------------------------

# ----------------------------------------
# 1. Extract raw met
# ----------------------------------------
met.vars <- c("tair", "precipf", "swdown", "lwdown", "qair", "psurf", "wind", "CO2")
site.list <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB")
met <- list() 
# -----------------------------------
# 1.1 ED2 -- all met vars except swdown
# -----------------------------------
for(s in 1:length(site.list)){
  dir.mod <- file.path(model.dir, "ED2.v7", site.list[s])
  files.mod <- dir(dir.mod)
  metvars.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.mod)){
    ncMT <- nc_open(file.path(dir.mod, files.mod[i]))
    for(v in which(!met.vars=="swdown")){
      if(i == 1) temp <- vector() else temp <- as.vector(metvars.list[[v]])
      temp <- c(temp, ncvar_get(ncMT, met.vars[v]))  
      metvars.list[[v]] <- temp
    }    
    nc_close(ncMT)      
  }
  names(metvars.list) <- met.vars[which(!met.vars=="swdown")]
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in which(!met.vars=="swdown")){
    if(s == 1){
      met[[v]] <- data.frame(metvars.list[[v]]) 
    } else {
      met[[v]][,s] <- metvars.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(met) <- met.vars
for(i in which(!met.vars=="swdown")){
  names(met[[i]]) <- site.list
}

# -----------------------------------
# 1.2 SiBCASA -- swdown
# -----------------------------------
for(s in 1:length(site.list)){
  dir.mod <- file.path(model.dir, "SiBCASA.v1", paste0(site.list[s], "_SiBCASA"))
  files.mod <- dir(dir.mod)
  metvars.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.mod)){
    if(i == 1) temp <- vector() 
    ncMT <- nc_open(file.path(dir.mod, files.mod[i]))
    temp <- c(temp, ncvar_get(ncMT, "swdown"))  
    nc_close(ncMT)      
  }
  #-----------------------------------
  # Adding variable groups to master model list
  if(s == 1){
      met[["swdown"]] <- data.frame(temp) 
    } else {
      met[["swdown"]][,s] <- temp
    }
} # Close the model loop
# Adding site label to each variable
names(met) <- met.vars
for(i in which(met.vars=="swdown")){
  names(met[[i]]) <- site.list
}
# -----------------------------------

# -----------------------------------
# 2. Stack Met to make something easier to work with
# -----------------------------------
year <- 850:2010
mo <- 1:12
year.mo <- merge(mo, year); names(year.mo) <- c("Month", "Year")
head(year.mo)

met.mo <- data.frame(Site    = stack(met$tair)[,2], 
				   Year    = year.mo$Year,
				   Month   = year.mo$Month,
                   tair    = stack(met$tair)[,1], 
                   precipf = stack(met$precipf)[,1], 
                   swdown  = stack(met$swdown)[,1], 
                   lwdown  = stack(met$lwdown)[,1],
                   qair    = stack(met$qair)[,1],
                   psurf   = stack(met$psurf)[,1],
                   wind    = stack(met$wind)[,1],
                   CO2     = stack(met$CO2)[,1]
                   )

summary(met.mo)
head(met.mo)
# -----------------------------------

# -----------------------------------
# 3. Aggregate met to growing season & year to compare
# -----------------------------------
met.yr <- aggregate(met.mo[,met.vars], by=list(met.mo$Site, met.mo$Year), FUN=mean)
names(met.yr) <- c("Site", "Year", paste(met.vars, "yr", sep="."))
met.yr[,paste(met.vars, "gs", sep=".")] <- aggregate(met.mo[met.mo$Month>=5 & met.mo$Month<=9,met.vars], by=list(met.mo[met.mo$Month>=5 & met.mo$Month<=9, "Site"], met.mo[met.mo$Month>=5 & met.mo$Month<=9, "Year"]), FUN=mean)[,3:(length(met.vars)+2)]
summary(met.yr)

write.csv(met.yr, file.path(dat.dir, "Drivers_Year_GrowingSeason.csv"), row.names=F)
# -----------------------------------
# ----------------------------------------

# ----------------------------------------
# 4. Compare growing season & yearly met
# ----------------------------------------
# Graphically Showing Correlation
pdf(file.path(fig.dir, "Driver_Correlations_Year_GrowingSeason_t001.pdf"))
par(mfrow=c(3,3), mar=c(2, 2, 2, 2))
for(v in unique(met.vars)){
	lm1 <- lm(met.yr[,paste0(v, ".gs")]~met.yr[,paste0(v, ".yr")])
print(
	plot(met.yr[,paste0(v,".gs")]~met.yr[,paste0(v,".yr")], main=paste0(v, ", R2=", round(summary(lm1)$r.squared, 2)), xlab="Year", ylab="Growing Season", 
	     xlim=range(c(met.yr[,paste0(v,".yr")], met.yr[,paste0(v,".gs")])),
	     ylim=range(c(met.yr[,paste0(v,".yr")], met.yr[,paste0(v,".gs")])))
)
print(
	abline(lm1, col="red")
	# abline(a=0, b=1, col="blue")
)
}
dev.off()

# Making a table of the correlations between 
# cor.met <- data.frame(var=as.factor(met.vars), cor.pearson=NA, r.squared=NA)
cor.met <- data.frame(var=as.factor(met.vars))
for(v in 1:length(cor.met$var)){
	var <- cor.met$var[v] 
	cor.met[v, "cor.pearson"] <- cor(met.yr[,paste0(var, ".yr")], met.yr[,paste0(var, ".gs")])
	cor.met[v, "r.squared"  ] <- summary(lm(met.yr[,paste0(var, ".gs")]~met.yr[,paste0(var, ".yr")]))$r.squared
}
cor.met
write.csv(cor.met, file.path(dat.dir, "Drivers_Year_GrowingSeason_Stats.csv"), row.names=F)
# ----------------------------------------

# ----------------------------------------
# Making a supplemental figure that graphs tair, precip, & CO2 in a similar manner as the NPP figure
# ----------------------------------------
sec2yr <- 1*60*60*24*365.25

# Read in met data
met.mod <- read.csv(file.path(dat.dir, "Drivers_Year_GrowingSeason.csv"))
met.mod$Data <- as.factor("PalEON")
met.mod[,c("tair.yr", "tair.gs")] <- met.mod[,c("tair.yr", "tair.gs")] - 273.15
met.mod$precipf.yr <- met.mod$precipf.yr*sec2yr
met.mod$precipf.gs <- met.mod$precipf.gs*sec2yr*5/12
summary(met.mod)

# Read in PRISM met data from tree rings
met.tr <- read.csv("Data/TreeRing_PRISM_Climate.csv")
met.tr$Site <- ifelse(substr(met.tr$PlotID,1,2) %in% c("LF", "TP"), "PHA",
                      ifelse(substr(met.tr$PlotID,1,2) %in% c("HO", "ME"), "PHO",
                             ifelse(substr(met.tr$PlotID,1,2)=="MN", "PDL", "PUN")))
met.tr$Site <- as.factor(met.tr$Site)
met.tr$Data <- as.factor("PRISM")
met.tr[,c("tair.yr", "tair.gs")] <- met.tr[,c("tair.yr", "tair.gs")] - 273.15
summary(met.tr)
unique(met.tr$PlotID)

# Merge the two-sites together
met.all <- merge(met.mod, met.tr, all.x=T, all.y=T)
summary(met.all)

summary(met.all[,paste0(vars.agg, ".gs")])
summary(met.all[,c("Site", "Data", "Year")])

# Aggregate to get region-level stats
vars.agg <- c("tair", "precipf", "CO2")
met.all.agg <- aggregate(met.all[,paste0(vars.agg, ".gs")], by=met.all[,c("Site", "Data", "Year")], FUN=mean, na.rm=T )
names(met.all.agg)[4:ncol(met.all.agg)] <- paste0(vars.agg, ".mean")
met.all.agg[,paste0(vars.agg, ".lo")] <- aggregate(met.all[,paste0(vars.agg, ".gs")], by=met.all[,c("Site", "Data", "Year")], FUN=quantile, 0.025, na.rm=T)[,paste0(vars.agg, ".gs")]
met.all.agg[,paste0(vars.agg, ".hi")] <- aggregate(met.all[,paste0(vars.agg, ".gs")], by=met.all[,c("Site", "Data", "Year")], FUN=quantile, 0.975, na.rm=T)[,paste0(vars.agg, ".gs")]
summary(met.all.agg)

png(file.path(fig.dir, "SuppFig1_Temperature_GrowingSeason.png"), width=8, height=6, units="in", res=120)
{
ggplot(data=met.all.agg) + facet_wrap(~Site) +
  geom_ribbon(aes(x=Year, ymin=tair.lo, ymax=tair.hi, fill=Data), alpha=0.5) +
  geom_line(aes(x=Year, y=tair.mean, color=Data), alpha=0.8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(name=expression(bold(paste("Growing Season Temperature ("^"o", "C)"))), expand=c(0,0)) + 
  scale_color_manual(values=c("black", "red3")) +
  scale_fill_manual(values=c("black", "red3")) +
  theme(legend.position=c(0.1, 0.37),
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=10),
        legend.key=element_blank(),
        legend.key.size=unit(1.5, "lines"),
        legend.background=element_blank()) +
  theme(strip.text.x=element_text(size=12, face="bold"),
        strip.text.y=element_text(size=12, face="bold")) + 
  theme(axis.line=element_line(color="black", size=0.5), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_rect(fill=NA, color="black", size=0.5), 
        panel.background=element_blank(),
        panel.margin.x=unit(0.5, "lines"),
        panel.margin.y=unit(0.5, "lines"))  +
  theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
        axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
        axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
        axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
        axis.ticks.length=unit(-0.5, "lines"))
}
dev.off()


png(file.path(fig.dir, "SuppFig2_Precipitation_GrowingSeason.png"), width=8, height=6, units="in", res=120)
{
ggplot(data=met.all.agg) + facet_wrap(~Site) +
  geom_ribbon(aes(x=Year, ymin=precipf.lo, ymax=precipf.hi, fill=Data), alpha=0.5) +
  geom_line(aes(x=Year, y=precipf.mean, color=Data), alpha=0.8) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(name=expression(bold(paste("Precipitation (mm yr"^"-1", ")"))), expand=c(0,0)) +
  scale_color_manual(values=c("black", "red3")) +
  scale_fill_manual(values=c("black", "red3")) +
  theme(legend.position=c(0.5, 0.37),
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=10),
        legend.key=element_blank(),
        legend.key.size=unit(1.5, "lines"),
        legend.background=element_blank()) +
  theme(strip.text.x=element_text(size=12, face="bold"),
        strip.text.y=element_text(size=12, face="bold")) + 
  theme(axis.line=element_line(color="black", size=0.5), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_rect(fill=NA, color="black", size=0.5), 
        panel.background=element_blank(),
        panel.margin.x=unit(0.5, "lines"),
        panel.margin.y=unit(0.5, "lines"))  +
  theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
        axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
        axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
        axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
        axis.ticks.length=unit(-0.5, "lines"))
}
dev.off()

# ----------------------------------------
