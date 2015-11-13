# ----------------------------------------
# Temporal Scaling Analyses
# Looking at changes in PalEON driver correllation & distribution with temporal scales
# Christy Rollinson, crollinson@gmail.com
# Date Created: 15 July 2015
# ----------------------------------------
# -------------------------
# Objectives & Overview
# -------------------------
# Question: Does annual climate correlate closely with the growing season?
#            i.e. do we need to re-run the model analyses trying for only growing season temps
# -------------------------
#
# -------------------------
# Input Data/Results:
# -------------------------
# 1) Raw Met drivers at MONTHLY scale -- extract from CLM-BGC & ED (lwdown)
# -------------------------
#
# -------------------------
# Interpretation Analyses:
# -------------------------
# A) correlations
# -------------------------
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
setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
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
# Compare growing season & yearly met
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
