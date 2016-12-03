# # ----------------------------------------
# Objective: Put all of the numbers in the manuscript in one place
# Christy Rollinson, crollinson@gmail.com
# Date Created: 21 April 2016
#
# This script pulls on numbers/analyses originally performed in scripts:
#   - 7_analysis_response_baseline_ensembles.R
#   - 8_analysis_sensitivity_TemporalExtent.R
# ----------------------------------------
rm(list=ls())

setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
load("Data/EcosysData.Rdata")


# ----------------------------------------
# Load libraries we may need
library(mgcv); library(ggplot2)
library(nlme)
# ----------------------------------------

# ----------------------------------------
# Evaluating correlation between biomasss and LAI for validity of using
# LAI as a biomass proxy in JULES
# ----------------------------------------
{
agb.lai <- data.frame(Model=unique(ecosys$Model))
for(m in agb.lai$Model){
  if(!max(ecosys[ecosys$Model==m & ecosys$Resolution=="t.001", "LAI"], na.rm=T)>0 | 
     !max(ecosys[ecosys$Model==m & ecosys$Resolution=="t.001", "AGB"], na.rm=T)>0 ) next
  agb.lai[agb.lai$Model==m, "corr"] <- cor(ecosys[!is.na(ecosys[ecosys$Model==m, "AGB"]) & ecosys$Model==m & ecosys$Resolution=="t.001", "LAI"], ecosys[!is.na(ecosys[ecosys$Model==m, "AGB"]) & ecosys$Model==m & ecosys$Resolution=="t.001", "AGB"])
}
mean(agb.lai$corr, na.rm=T); sd(agb.lai$corr, na.rm=T)
agb.lai
}
# ----------------------------------------



# ----------------------------------------
# Comparing raw NPP & model stats (Results Paragraph 1)
# from 7_analysis_response_baseline_ensembles.R
# ----------------------------------------
{
load("Data/analyses/analysis_baseline/post-process_baseline.RData")

summary(dat.ecosys)

# ------------------------
# Comparison of ensemble mean NPP across sites
# ------------------------
{
fac.ecosys <- c("Y", "Y.rel", "tair", "precipf", "CO2")
dat.ecosys2 <- aggregate(dat.ecosys[,fac.ecosys],  
                         by=dat.ecosys[,c("Y.type", "data.type", "Site", "Year")],
                         FUN=mean, na.rm=T)
dat.ecosys2[,paste0(fac.ecosys, ".sd")] <- aggregate(dat.ecosys[,fac.ecosys],  
                                                     by=dat.ecosys[,c("Y.type", "data.type", "Site", "Year")],
                                                     FUN=sd, na.rm=T)[,fac.ecosys]
summary(dat.ecosys2[dat.ecosys2$data.type=="Model",])


sites.stat <- aggregate(dat.ecosys2[,c("Y", "Y.rel", "Y.sd", "Y.rel.sd")],
                        by=dat.ecosys2[,c("Y.type","data.type", "Site")],
                        FUN=mean, na.rm=T)
sites.stat <- sites.stat[sites.stat$data.type=="Model",]

# Comparison between models & Tree rings at Havard and Howland
sites.stat2 <- aggregate(dat.ecosys2[dat.ecosys2$Year>=1980, c("Y", "Y.rel", "Y.sd", "Y.rel.sd")],
                         by=dat.ecosys2[dat.ecosys2$Year>=1980,c("Y.type", "data.type", "Site")],
                         FUN=mean, na.rm=T)

sites.stat2[sites.stat2$data.typ=="Model",]
mean(dat.ecosys[dat.ecosys$Year>=1980 & dat.ecosys$Model=="jules.stat" & dat.ecosys$Site=="PHA", "Y"]); sd(dat.ecosys[dat.ecosys$Year>=1980 & dat.ecosys$Model=="jules.stat" & dat.ecosys$Site=="PHA", "Y"])
mean(dat.ecosys[dat.ecosys$Year>=1980 & dat.ecosys$Model=="jules.stat" & dat.ecosys$Site=="PHO", "Y"]); sd(dat.ecosys[dat.ecosys$Year>=1980 & dat.ecosys$Model=="jules.stat" & dat.ecosys$Site=="PHO", "Y"])

mean(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Tree Rings" & dat.ecosys2$Y.type=="NPP" & dat.ecosys2$Site=="PHA", "Y"]); sd(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Tree Rings" & dat.ecosys2$Y.type=="NPP" & dat.ecosys2$Site=="PHA", "Y"])
mean(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Tree Rings" & dat.ecosys2$Y.type=="NPP" & dat.ecosys2$Site=="PHO", "Y"]); sd(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Tree Rings" & dat.ecosys2$Y.type=="NPP" & dat.ecosys2$Site=="PHO", "Y"])



# Mean and SD for site with lowest NPP
sites.stat[which(sites.stat$Y==min(sites.stat$Y)), c("Site", "Y", "Y.sd")]

# Mean & SD for site with highest NPP
sites.stat[which(sites.stat$Y==max(sites.stat$Y)), c("Site", "Y", "Y.sd")]


# Statistically comparing
models.stat <- aggregate(dat.ecosys[dat.ecosys$data.type=="Model",fac.ecosys],  
                         by=dat.ecosys[dat.ecosys$data.type=="Model",c("Model","Site")],
                         FUN=mean, na.rm=T)
summary(models.stat)


sites.diff.lme <- lme(Y ~ Site, random=list(Model=~1), data=models.stat[,], na.action=na.omit)
sites.diff <- aov(Y ~ Site, data=models.stat[,], na.action=na.omit)
TukeyHSD(sites.diff)

summary(sites.diff.lme)

sites.diff.lme2 <- lme(Y ~ Site, random=list(Model=~1, Year=~1), data=dat.ecosys[dat.ecosys$data.type=="Model",], na.action=na.omit)
summary(sites.diff.lme2)
anova(sites.diff)

# Comparing models
models.stat2 <- aggregate(dat.ecosys[dat.ecosys$data.type=="Model",c("Y", "Y.rel")],  
                          by=dat.ecosys[dat.ecosys$data.type=="Model",c("Model", "Model.Order")],
                          FUN=mean, na.rm=T)
models.stat2[,paste0(c("Y", "Y.rel"),".sd")] <- aggregate(dat.ecosys[dat.ecosys$data.type=="Model",c("Y", "Y.rel")],  
                                                          by=dat.ecosys[dat.ecosys$data.type=="Model",c("Model", "Model.Order")],
                                                          FUN=sd, na.rm=T)[c("Y", "Y.rel")]
models.stat2

# Mean and SD for model with lowest NPP
models.stat2[which(models.stat2$Y==min(models.stat2$Y)), c("Model", "Y", "Y.sd")]

# Mean & SD for model with highest NPP
models.stat2[which(models.stat2$Y==max(models.stat2$Y)), c("Model", "Y", "Y.sd")]

# Comparing model variability through time
models.stat3a <- aggregate(dat.ecosys[,c("Y", "Y.rel")],  
                          by=dat.ecosys[,c("data.type", "Model", "Model.Order", "Year")],
                          FUN=mean, na.rm=T)
models.stat3 <- aggregate(models.stat3a[,c("Y", "Y.rel")],  
                           by=models.stat3a[,c("data.type", "Model", "Model.Order")],
                           FUN=mean, na.rm=T)
models.stat3[,paste0(c("Y", "Y.rel"),".sd")] <- aggregate(models.stat3a[,c("Y", "Y.rel")],  
                                                          by=models.stat3a[,c("data.type", "Model", "Model.Order")],
                                                          FUN=sd, na.rm=T)[,c("Y", "Y.rel")]
models.stat3

# Mean and SD for model with highest absolute NPP variability
models.stat3[which(models.stat3$Y.sd==max(models.stat3$Y.sd)), c("Model", "Y.sd", "Y.rel.sd")]

# Mean & SD for model with highest relative NPP variability
models.stat3[which(models.stat3$Y.rel.sd==max(models.stat3$Y.rel.sd)), c("Model", "Y.sd", "Y.rel.sd")]

# Mean and SD for model with lowest absolute NPP variability
models.stat3[which(models.stat3$Y.sd==min(models.stat3$Y.sd)), c("Model", "Y.sd", "Y.rel.sd")]

# Mean & SD for model with lowest relative NPP variability
models.stat3[which(models.stat3$Y.rel.sd==min(models.stat3$Y.rel.sd)), c("Model", "Y.sd", "Y.rel.sd")]

# Getting the mean & sd for models alone
mean(models.stat3[models.stat3$data.type=="Model", "Y.rel.sd"]); sd(models.stat3[models.stat3$data.type=="Model", "Y.rel.sd"])

}
# ------------------------

# ------------------------
# Generating table 2
# ------------------------
{
  
  # Condensing model variability across space and time to get general model characteristics

  vars.agg <- c("Y", "Y.rel", "Biomass")
  mod.site                           <- aggregate(dat.ecosys[,vars.agg], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                                  FUN=mean)
  mod.site[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys[,vars.agg], 
                                                  by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site")], 
                                                  FUN=sd)[,vars.agg]
  mod.site$Biomass.sd.per <- mod.site$Biomass.sd/mod.site$Biomass
  summary(mod.site)
  
  mod.agg                           <- aggregate(mod.site[,c(vars.agg, paste0(vars.agg, ".sd"), "Biomass.sd.per")], 
                                                 by=mod.site[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=mean)
  mod.agg
  
  # Finding the change in key variables in the modern era
  for(v in vars.agg){
    mod.agg[,paste0("dModern.", v)] <- NA  
  }
  
  for(m in unique(mod.agg$Model)){
    mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1990 & dat.ecosys$Year<=2010, vars.agg], na.rm=T) -
      colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1830 & dat.ecosys$Year<=1850, vars.agg], na.rm=T)
  }
  mod.agg
  
  # Add in the vegetation scheme
  mod.agg$veg.scheme <- as.factor(ifelse(mod.agg$Model %in% c("clm.bgc", "clm.cn", "sibcasa", "jules.stat"), "Static", "Dynamic"))
  
  # Finding mean fire return
  output.all <- read.csv("../../phase1a_output_variables/PalEON_MIP_Yearly.csv")
  summary(output.all)
  
  agg.fire <- aggregate(output.all[,c("Fire")], by=output.all[,c("Model", "Updated")], FUN=mean, na.rm=T)
  names(agg.fire)[3] <- "Fire"
  agg.fire[is.na(agg.fire$Fire),"Fire"] <- 0
  agg.fire$Fire <- agg.fire$Fire*(60*60*24*365)*(100*100)*1e-3
  agg.fire
  
  mod.agg <- merge(mod.agg, agg.fire[,c("Model", "Fire")], all.x=T, all.y=F)
  mod.agg$Fire <- mod.agg$Fire # changing fire from KgC/m2/s to MgC/HA/yr (same units as NPP)
  mod.agg$fire.scheme <- as.factor(ifelse(mod.agg$Fire>0, "Yes", "No"))
  
  summary(mod.agg)
  
  # --------
  # Need to aggregate to relativize Biomass & Evergreen so that we can get the sd across 
  #   space & time that isn't confounded by baseline differences in model biomass, etc
  # --------
  # First need to get the composition variables in there
  dat.ecosys <- merge(dat.ecosys, ecosys[ecosys$Resolution=="t.001",c("Model", "Site", "Year", "Evergreen", "Deciduous", "Grass")], all.x=T, all.y=F)
  summary(dat.ecosys)
  
  for(m in unique(dat.ecosys$Model)){
    #   biomass.mean <- mean(dat.ecosys[dat.ecosys$Model==m, "Biomass"], na.rm=T)
    #   dat.ecosys[dat.ecosys$Model==m, "Biomass.rel"] <- dat.ecosys[dat.ecosys$Model==m, "Biomass"]/biomass.mean
    for(p in unique(dat.ecosys$PlotID)){
      biomass.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Biomass"], na.rm=T)
      dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Biomass.rel"] <- dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Biomass"]/biomass.mean  
    }
    
    #   evergreen.mean <- mean(dat.ecosys[dat.ecosys$Model==m, "Evergreen"], na.rm=T)
    if(is.na(mean(dat.ecosys[dat.ecosys$Model==m, "Evergreen"], na.rm=T))) next # Skip things that we don't have composition shift handy for
    for(p in unique(dat.ecosys$PlotID)){
      evergreen.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Evergreen"], na.rm=T)
      dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Evergreen.rel"] <- dat.ecosys[dat.ecosys$Model==m & dat.ecosys$PlotID==p, "Evergreen"]/evergreen.mean  
    }
  }
  summary(dat.ecosys)
  
  for(m in unique(mod.agg[,"Model"])){  
    mod.agg[mod.agg$Model==m,"Y.sd"]             <- sd(dat.ecosys[dat.ecosys$Model==m, "Y"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Y.rel.sd"]         <- sd(dat.ecosys[dat.ecosys$Model==m, "Y.rel"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Biomass.mean"]     <- mean(dat.ecosys[dat.ecosys$Model==m, "Biomass"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Biomass.rel.sd"]   <- sd(dat.ecosys[dat.ecosys$Model==m, "Biomass.rel"], na.rm=T)
    
    if(is.na(mean(dat.ecosys[dat.ecosys$Model==m, "Evergreen"], na.rm=T))) next
    mod.agg[mod.agg$Model==m,"Evergreen.mean"]   <- mean(dat.ecosys[dat.ecosys$Model==m, "Evergreen"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Deciduous.mean"]   <- mean(dat.ecosys[dat.ecosys$Model==m, "Deciduous"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Grass.mean"]       <- mean(dat.ecosys[dat.ecosys$Model==m, "Grass"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Evergreen.sd"]     <- sd(dat.ecosys[dat.ecosys$Model==m, "Evergreen"], na.rm=T)
    mod.agg[mod.agg$Model==m,"Evergreen.rel.sd"] <- sd(dat.ecosys[dat.ecosys$Model==m, "Evergreen.rel"], na.rm=T)
  }
  mod.agg
  
  table2 <- mod.agg[mod.agg$data.type=="Model", c("Model.Order", "veg.scheme", "Evergreen.mean", "Evergreen.sd", "fire.scheme", "Fire", "Biomass.rel.sd")]
  table2[,c("Evergreen.mean", "Evergreen.sd", "Fire", "Biomass.rel.sd")] <- round(table2[,c("Evergreen.mean", "Evergreen.sd", "Fire", "Biomass.rel.sd")], 2)
  names(table2) <- c("Model", "Vegetation Scheme", "Mean Fraction Evergreen", "Composition Variability", "Fire Occurrence", "Mean Fire Magnitude", "Biomass Variability")
  
  table2
  write.csv(table2, file.path("Data/analyses/analysis_baseline", "Table2_ModelCharacterstics.csv"), row.names=F)

# Finding who had the least & most variable composition
table2[which(table2[,"Composition Variability"]==min(table2[,"Composition Variability"], na.rm=T)), c("Model", "Composition Variability")]
table2[which(table2[,"Composition Variability"]==max(table2[,"Composition Variability"], na.rm=T)), c("Model", "Composition Variability")]

# Finding who had the least & most variable biomass
table2[which(table2[,"Biomass Variability"]==min(table2[,"Biomass Variability"], na.rm=T)), c("Model", "Biomass Variability")]
table2[which(table2[,"Biomass Variability"]==max(table2[,"Biomass Variability"], na.rm=T)), c("Model", "Biomass Variability")]

}  
# ------------------------
}
# ----------------------------------------


# ----------------------------------------
# Comparing NPP Sensitivity to climate & CO2 at long & short scales (most of results)
# from 8_analysis_sensitivity_TemporalExtent.R
# ----------------------------------------
{
load("Data/analyses/analysis_TempExtent/post-process_TempExtent.RData")


# --------
# 5.b.0. Adding some model-level stats to compare the relative sensitivities
# --------
summary(ci.terms)
{
# 1. Quantify with-in model sensitivity shifts due to change in temporal extent
# calculate mean.rel.cent
{
  # i. classify met by which extent it's observed in 
  for(v in c("tair", "precipf", "CO2")){ # go through met (same for all models/data)
    range.1980 <- range(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$data.type=="Model", v], na.rm=T)
    range.1901 <- range(dat.ecosys2[dat.ecosys2$Year>=1901 & dat.ecosys2$data.type=="Model", v], na.rm=T)
    
    met.fill <-  as.factor(ifelse(ci.terms[ci.terms$Effect==v,"x"]>=range.1980[1] & 
                                    ci.terms[ci.terms$Effect==v,"x"]<=range.1980[2],
                                  "1980-2010", 
                                  ifelse(ci.terms[ci.terms$Effect==v,"x"]>=range.1901[1] & 
                                           ci.terms[ci.terms$Effect==v,"x"]<=range.1901[2],
                                         "1901-2010",
                                         "850-2010"
                                  ))
    )  
    ci.terms[ci.terms$Effect==v,"Met"] <- as.factor(met.fill)  
  }
  summary(ci.terms)
  
  # ii.  Adjusting the sensitivity so that everywhere has the same x-intercept as 1980-2010
  for(v in c("tair", "precipf", "CO2")){
    if(v=="CO2") r=0 else if(v=="tair") r=1 else r=-1
    
    for(m in unique(ci.terms$Model)){
      var.mid <- mean(ci.terms[ci.terms$Effect==v & ci.terms$Met=="1980-2010", "x"], na.rm=T)
      offset.1980 <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="1980-2010","mean.rel"])
      offset.1901 <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="1901-2010","mean.rel"])
      offset.850  <- mean(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & round(ci.terms$x, r)==round(var.mid,r) & ci.terms$Extent=="850-2010","mean.rel"])
      
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","mean.rel"] - offset.1980
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","mean.rel"] - offset.1901
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"mean.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"mean.rel"] - offset.850
      
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","upr.rel"] - offset.1980
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","upr.rel"] - offset.1901
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"upr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"upr.rel"] - offset.850
      
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010","lwr.rel"] - offset.1980
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010","lwr.rel"] - offset.1901
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"lwr.rel.cent"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010" ,"lwr.rel"] - offset.850
      
    }  
  }

  summary(ci.terms)
}

# 2. Condensing model variability across space and time to get general model characteristics
# calculate mod.agg
{
  vars.agg <- c("Y", "Y.rel", "Biomass")
  mod.agg                           <- aggregate(dat.ecosys2[,vars.agg], 
                                                 by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=mean)
  mod.agg[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys2[,vars.agg], 
                                                 by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=sd)[,vars.agg]
  mod.agg
  
  
  # Finding the change in key variables in the modern era
  for(v in vars.agg){
    mod.agg[mod.agg$Model==m,paste0("dModern.", v)] <- NA  
  }
  
  for(m in unique(mod.agg$Model)){
    mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=1990 & dat.ecosys2$Year<=2010, vars.agg], na.rm=T) -
      colMeans(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=1830 & dat.ecosys2$Year<=1850, vars.agg], na.rm=T)
  }
  mod.agg
  
  # Add in the vegetation scheme
  mod.agg$veg.scheme <- as.factor(ifelse(mod.agg$Model %in% c("clm.bgc", "clm.cn", "sibcasa", "jules.stat"), "Static", "Dynamic"))
  
  # Finding mean fire return
{
    output.all <- read.csv("../../phase1a_output_variables/PalEON_MIP_Yearly.csv")
    summary(output.all)
    
    agg.fire <- aggregate(output.all[,c("Fire")], by=output.all[,c("Model", "Updated")], FUN=mean, na.rm=T)
    names(agg.fire)[3] <- "Fire"
    agg.fire[is.na(agg.fire$Fire),"Fire"] <- 0
    agg.fire$Fire <- agg.fire$Fire*(60*60*24*365)*(100*100)*1e-3
    agg.fire
    
    mod.agg <- merge(mod.agg, agg.fire[,c("Model", "Fire")], all.x=T, all.y=F)
    mod.agg$Fire <- mod.agg$Fire # changing fire from KgC/m2/s to MgC/HA/yr (same units as NPP)
    mod.agg$fire.scheme <- as.factor(ifelse(mod.agg$Fire>0, "Yes", "No"))
    
    summary(mod.agg)
  }
}

ci.terms <- merge(ci.terms, mod.agg, all.x=T, all.y=F)
summary(ci.terms)

# 3. Getting the variability of slow processes at each time scale to use as explanatory factors
#    Slow Processes = Composition, Biomass
# NOTE: This overwrites a lot of what we already calcualted in mod.ag
# Add to mod.agg with the appropriate scales
{
  # --------
  # Need to aggregate to relativize Biomass & Evergreen so that we can get the sd across 
  #   space & time that isn't confounded by baseline differences in model biomass, etc
  # --------
  # First need to get the composition variables in there
  dat.ecosys <- merge(dat.ecosys, ecosys[ecosys$Resolution=="t.001",c("Model", "Site", "Year", "Evergreen", "Deciduous", "Grass")], all.x=T, all.y=F)
  dat.ecosys2 <- merge(dat.ecosys2, ecosys[ecosys$Resolution=="t.001",c("Model", "Site", "Year", "Evergreen", "Deciduous", "Grass")], all.x=T, all.y=F)
  summary(dat.ecosys)
  summary(dat.ecosys2)
  
  for(m in unique(dat.ecosys2[dat.ecosys2$data.type=="Model", "Model"])){
    if(m == "TreeRingNPP") ext.full <- "1980-2010" else if(m=="TreeRingRW") ext.full <- "1901-2010" else ext.full <- "850-2010"
    #   biomass.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Biomass"], na.rm=T)
    #   dat.ecosys2[dat.ecosys2$Model==m, "Biomass.rel"] <- dat.ecosys2[dat.ecosys2$Model==m, "Biomass"]/biomass.mean
    for(p in unique(dat.ecosys$PlotID)){
      biomass.mean <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass"], na.rm=T)
      dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass.rel"] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Biomass"]/biomass.mean  
    }
    
    if(is.na(mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Evergreen"], na.rm=T))) next # Skip things that we don't have composition shift handy for
    #   evergreen.mean <- mean(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Extent==ext.full, "Evergreen"], na.rm=T)
    #   dat.ecosys2[dat.ecosys2$Model==m, "Evergreen.rel"] <- dat.ecosys2[dat.ecosys2$Model==m, "Evergreen"]/evergreen.mean  
    for(p in unique(dat.ecosys$PlotID)){
      evergreen.mean <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen"], na.rm=T)
      dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen.rel"] <- dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$PlotID==p, "Evergreen"]/evergreen.mean  
    } 
  }
  summary(dat.ecosys2)
  
  for(m in unique(ci.terms[ci.terms$data.type=="Model","Model"])){
    for(e in unique(ci.terms$Extent)){
      yr.min <- as.numeric(strsplit(e, "-")[[1]][1])
      
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Y.sd"]             <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Y"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Y.rel.sd"]         <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Y.rel"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Biomass.mean"]     <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Biomass"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Biomass.rel.sd"]   <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Biomass.rel"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.mean"]   <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Deciduous.mean"]   <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Deciduous"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Grass.mean"]       <- mean(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Grass"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.sd"]     <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen"], na.rm=T)
      ci.terms[ci.terms$Model==m & ci.terms$Extent==e,"Evergreen.rel.sd"] <- sd(dat.ecosys2[dat.ecosys2$Model==m & dat.ecosys2$Year>=yr.min, "Evergreen.rel"], na.rm=T)
      
    }
  }
}
# --------
summary(ci.terms)


# 4. Condensing the full sensitivity curves to the values at the 25, 50, and 75% for 
#    analysis of continuous characteristics of models (i.e. NPP, modern change, etc)
{
  summary(ci.terms)
  summary(dat.ecosys2)
  
  factors.analy <- c("Y.sd", "Y.rel.sd","Biomass.sd", "Evergreen.sd")
  factors.agg <- c("Effect", "Extent", "veg.scheme", "fire.scheme", "Model", "x", "Quantile")
  
  co2.driver     <- dat.ecosys2[dat.ecosys2$Model=="ed2","CO2"    ]
  tair.driver    <- dat.ecosys2[dat.ecosys2$Model=="ed2","tair"   ]
  precipf.driver <- dat.ecosys2[dat.ecosys2$Model=="ed2","precipf"]
  
  ci.terms[,"Quantile"] <- "Other"
  for(e in c("tair", "precipf", "CO2")){
    # Get the distirbution of met drivers for each driver & specify to what precision we want to round
    effect.driver <- dat.ecosys2[dat.ecosys2$Model=="ed2", e]
    if(e=="CO2") r=0 else if(e=="tair") r=1 else r=-1
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.05, na.rm=T),r), "Quantile"] <- "q05"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.25, na.rm=T),r), "Quantile"] <- "q25"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.50, na.rm=T),r), "Quantile"] <- "q50"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.75, na.rm=T),r), "Quantile"] <- "q75"
    ci.terms[ci.terms$Effect==e & round(ci.terms$x, r)==round(quantile(effect.driver, 0.95, na.rm=T),r), "Quantile"] <- "q95"
  }
  ci.terms$Quantile <- as.factor(ci.terms$Quantile)
  summary(ci.terms)
  
  factors.agg <- c("x", factors.analy, "Fire", "mean.rel", "lwr.rel", "upr.rel")
  ci.terms.agg <- aggregate(ci.terms[,factors.agg], 
                            by=ci.terms[,c("Effect", "Extent", "data.type", "Y.type", "Model", "Model.Order", "veg.scheme", "fire.scheme", "Quantile")],
                            FUN=mean, na.rm=T)
  summary(ci.terms.agg)
}

# 5. Find the model-relative shift in sensitivity by comparing the slopes (1st derivatives) relative to a base extent
#    -- Base Extent == "1980-2010
# mean.cent.deriv = deriv(mean.rel.cent)/deriv(x) ==> units are %change NPP per unit X (C, mm yr-1, or ppm)
{
# Trying to relativize the sensitivity to look at the change from the base extent
# Using absolute value so that positive indicates more extreme even though this 
#   leads to some very weird shapes
ext.base <- "1980-2010"
for(m in unique(ci.terms[!is.na(ci.terms$mean.rel.cent),"Model"])){
  for(v in c("tair", "precipf", "CO2")){
    mod.ext <- unique(ci.terms[ci.terms$Model==m & !is.na(ci.terms$mean.rel.cent),"Extent"]) 
    for(e in mod.ext){
      dif.x   <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "x"            ], lag=1))
      dif.y   <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"], lag=1))
      # dif.lwr <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "lwr.rel.cent" ], lag=1))
      # dif.upr <- c(NA, diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "upr.rel.cent" ], lag=1))
      
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"]   <- dif.y/dif.x
      # ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "lwr.cent.deriv" ] <- dif.lwr/dif.x
      # ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "upr.cent.deriv" ] <- dif.upr/dif.x
    }
    
    # Looking at the change in slope relative to a baseline -- looking at % change r
    if(!length(mod.ext)>1) next # Skip models/datasets where we only have one extent
    deriv.base <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==ext.base, "mean.cent.deriv"]
    cent.base  <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==ext.base, "mean.rel.cent"]
    for(e in mod.ext[which(!mod.ext==ext.base)]){
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev"]     <-     ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"] -      deriv.base
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev.abs"] <- abs(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"]) - abs(deriv.base)
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.deriv.dev.per"] <- (ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.deriv"] - deriv.base)/deriv.base
      
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev"]     <-     ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"] -      cent.base
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev.per"] <- (ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"] - cent.base)/cent.base
      ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.cent.dev.abs"] <-  abs(ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e, "mean.rel.cent"]) - abs(cent.base)
      
    }
  }
}

  
# Making sure we're not comparing the parts that are extrapolated beyond the range of temperatures
for(v in c("tair", "precipf", "CO2")){
  if(v=="CO2") r=0 else if(v=="tair") r=1 else r=-1
  
  for(m in unique(ci.terms$Model)){
    # Trimming the range of conditions in each dataset
    range.1980 <- range(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$Model==m, v], na.rm=T)
    range.1901 <- range(dat.ecosys2[dat.ecosys2$Year>=1901 & dat.ecosys2$Model==m, v], na.rm=T)
    range.0850 <- range(dat.ecosys2[dat.ecosys2$Model==m, v], na.rm=T)
    
    rows.1980 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010" & ci.terms$x>=min(range.1980) & ci.terms$x<=max(range.1980))
    rows.1901 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010" & ci.terms$x>=min(range.1901) & ci.terms$x<=max(range.1901))
    rows.0850 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010"  & ci.terms$x>=min(range.0850) & ci.terms$x<=max(range.0850))
    
    ci.terms[rows.1980,"cent.deriv2"] <- ci.terms[rows.1980,"mean.cent.deriv"]
    ci.terms[rows.1901,"cent.deriv2"] <- ci.terms[rows.1901,"mean.cent.deriv"]
    ci.terms[rows.0850,"cent.deriv2"] <- ci.terms[rows.0850 ,"mean.cent.deriv"]
  }  
}
summary(ci.terms)


summary(ci.terms)
summary(ci.terms[ci.terms$Effect=="precipf",])

}

summary(ci.terms)
}
# --------


# --------
# Getting stats using mean.cent.deriv (d%NPP/dx)
# --------
summary(ci.terms)
{  
# Putting "NONE" in for veg scheme & fire scheme for tree rings
ci.terms$veg.scheme <- as.factor(ifelse(ci.terms$data.type=="Tree Rings", "NONE", ci.terms$veg.scheme))
ci.terms$fire.scheme <- as.factor(ifelse(ci.terms$data.type=="Tree Rings", "NONE", ci.terms$fire.scheme))

# calculating the model and ensemble average derivative for each factor
mod.deriv <- aggregate(ci.terms[, c("mean.rel.cent", "mean.cent.deriv")],
                       by=ci.terms[, c("Effect", "Extent", "data.type", "Model", "Model.Order", "veg.scheme", "fire.scheme")],
                       FUN=mean, na.rm=T)
mod.deriv[,paste0(c("mean.rel.cent", "mean.cent.deriv"), ".sd")] <- aggregate(ci.terms[, c("mean.rel.cent", "mean.cent.deriv")],
                                                                              by=ci.terms[, c("Effect", "Extent", "data.type", "Model", "Model.Order", "veg.scheme", "fire.scheme")],
                                                                              FUN=sd, na.rm=T)[,c("mean.rel.cent", "mean.cent.deriv")]
mod.deriv <- mod.deriv[!mod.deriv$Effect=="Biomass",]
summary(mod.deriv)

ens.deriv <- aggregate(mod.deriv[mod.deriv$data.type=="Model",c("mean.rel.cent", "mean.cent.deriv")], 
                       by=mod.deriv[mod.deriv$data.type=="Model",c("Effect", "Extent")],
                       FUN=mean, na.rm=T)
ens.deriv[,paste0(c("mean.rel.cent", "mean.cent.deriv"), ".sd")] <- aggregate(mod.deriv[,c("mean.rel.cent", "mean.cent.deriv")], 
                                                                              by=mod.deriv[,c("Effect", "Extent")],
                                                                              FUN=sd, na.rm=T)[,c("mean.rel.cent", "mean.cent.deriv")]

ens.deriv
round(ens.deriv$mean.cent.deriv.sd*100, 2)
round(mod.deriv[mod.deriv$data.type=="Tree Rings", "mean.cent.deriv.sd"]*100, 2)

mod.deriv[mod.deriv$data.type=="Tree Rings", ]
}
# -----------

# -----------
# Ensemble & Model Stats (850-2010)
# -----------
{  
# Tair
ens.deriv[ens.deriv$Extent=="850-2010" & ens.deriv$Effect=="tair", c("mean.cent.deriv", "mean.cent.deriv.sd")]*100
mod.deriv[which(mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="tair" & mod.deriv[, "mean.cent.deriv"]==min(mod.deriv[mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="tair", "mean.cent.deriv"])), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]
mod.deriv[which(mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="tair" & mod.deriv[, "mean.cent.deriv"]==max(mod.deriv[mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="tair", "mean.cent.deriv"])), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]

# Precip
ens.deriv[ens.deriv$Extent=="850-2010" & ens.deriv$Effect=="precipf", c("mean.cent.deriv", "mean.cent.deriv.sd")]*100
mod.deriv[which(mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="precipf" & mod.deriv[,"mean.cent.deriv"]==min(mod.deriv[mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="precipf", "mean.cent.deriv"])), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]
mod.deriv[which(mod.deriv$Effect=="precipf" & mod.deriv[,"mean.cent.deriv"]==max(mod.deriv[mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="precipf", "mean.cent.deriv"])), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]

# CO2
ens.deriv[ens.deriv$Extent=="850-2010" & ens.deriv$Effect=="CO2", c("mean.cent.deriv", "mean.cent.deriv.sd")]*100
mod.deriv[mod.deriv$Extent=="850-2010" & !mod.deriv$Model=="linkages" & mod.deriv$Effect=="CO2" & mod.deriv[, "mean.cent.deriv"]==min(mod.deriv[!mod.deriv$Model=="linkages" & mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="CO2", "mean.cent.deriv"]), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]
mod.deriv[mod.deriv$Extent=="850-2010" & !mod.deriv$Model=="linkages" & mod.deriv$Effect=="CO2" & mod.deriv[, "mean.cent.deriv"]==max(mod.deriv[mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="CO2", "mean.cent.deriv"]), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]
# LINKAGES CO2
mod.deriv[which(mod.deriv$Extent=="850-2010" & mod.deriv$Effect=="CO2" &mod.deriv$Model=="linkages"), c("Effect", "Model", "veg.scheme", "fire.scheme", "mean.cent.deriv", "mean.cent.deriv.sd")]
}
# -----------

# -----------
# Calculating Changes in Ensemble agreement across scales
# -----------
{
  source("R/0_Calculate_GAMM_Posteriors.R")

  # 6. Making the data frame so we can graph & compare the curves quantitatively
  {
    factors.analy <- c("Y.sd", "Y.rel.sd","Biomass.sd", "Evergreen.sd")
    factors.agg <- c("Effect", "Extent", "veg.scheme", "fire.scheme", "Model", "x", "Quantile")
    df.co2 <- aggregate(ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",factors.analy], 
                        by=ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",factors.agg],
                        FUN=mean, na.rm=T)
    summary(df.co2)
    
    
    df.tair <- aggregate(ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",factors.analy], 
                         by=ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",factors.agg],
                         FUN=mean, na.rm=T)
    summary(df.tair)
    
    df.precipf <- aggregate(ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",factors.analy], 
                            by=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",factors.agg],
                            FUN=mean, na.rm=T)
    summary(df.precipf)
  }

  co2.ext     <- gam(mean.cent.deriv ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.ext    <- gam(mean.cent.deriv ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.ext <- gam(mean.cent.deriv ~ s(x, by=Extent), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  # Plot the different sensitivity curves by characteristic
  factors.agg <- c("Effect", "Extent", "veg.scheme", "fire.scheme", "Model", "x", "Quantile")
  
  co2.ext.post <- post.distns(model.gam=co2.ext, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.ext.post[,factors.agg] <- df.co2[,factors.agg]
  summary(co2.ext.post)
  
  tair.ext.post <- post.distns(model.gam=tair.ext, model.name="tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.ext.post[,factors.agg] <- df.tair[,factors.agg]
  summary(tair.ext.post)
  
  precipf.ext.post <- post.distns(model.gam=precipf.ext, model.name="precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.ext.post[,factors.agg] <- df.precipf[,factors.agg]  
  summary(precipf.ext.post)
  
  ext.post <- rbind(tair.ext.post, precipf.ext.post, co2.ext.post)
  summary(ext.post)
  
  for(m in unique(ci.terms[ci.terms$data.type=="Model", "Model"])){
    for(v in c("tair", "precipf", "CO2")){
      for(e in unique(ci.terms[ci.terms$Model==m, "Extent"])){
        mean.ens <- aggregate(ci.terms[ci.terms$Effect==v & ci.terms$Extent==e,"mean.rel.cent"], 
                              by=list(ci.terms[ci.terms$Effect==v & ci.terms$Extent==e,"x"]), 
                              FUN=mean)
        
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"dev.ensemble"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"mean.rel.cent"] - ext.post[ext.post$Extent==e & ext.post$Effect==v & ext.post$Model==m , "mean"]
        ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"dev.ensemble2"] <- ci.terms[ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent==e,"mean.rel.cent"] - mean.ens$x
      }
    }
  }
  
  # Making sure we're not comparing the parts that are extrapolated beyond the range of temperatures
  for(v in c("tair", "precipf", "CO2")){
    if(v=="CO2") r=0 else if(v=="tair") r=1 else r=-1
    
    for(m in unique(ci.terms$Model)){
      # Trimming the range of conditions in each dataset
      range.1980 <- range(dat.ecosys2[dat.ecosys2$Year>=1980 & dat.ecosys2$Model==m, v], na.rm=T)
      range.1901 <- range(dat.ecosys2[dat.ecosys2$Year>=1901 & dat.ecosys2$Model==m, v], na.rm=T)
      range.0850 <- range(dat.ecosys2[dat.ecosys2$Model==m, v], na.rm=T)
      
      rows.1980 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1980-2010" & ci.terms$x>=min(range.1980) & ci.terms$x<=max(range.1980))
      rows.1901 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="1901-2010" & ci.terms$x>=min(range.1901) & ci.terms$x<=max(range.1901))
      rows.0850 <- which(ci.terms$Model==m & ci.terms$Effect==v & ci.terms$Extent=="850-2010"  & ci.terms$x>=min(range.0850) & ci.terms$x<=max(range.0850))
      
      ci.terms[rows.1980,"dev.ensemble3a"] <- abs(ci.terms[rows.1980,"dev.ensemble"])
      ci.terms[rows.1901,"dev.ensemble3a"] <- abs(ci.terms[rows.1901,"dev.ensemble"])
      ci.terms[rows.0850,"dev.ensemble3a"] <- abs(ci.terms[rows.0850,"dev.ensemble"])
      
      ci.terms[rows.1980,"dev.ensemble3b"] <- abs(ci.terms[rows.1980,"dev.ensemble2"])
      ci.terms[rows.1901,"dev.ensemble3b"] <- abs(ci.terms[rows.1901,"dev.ensemble2"])
      ci.terms[rows.0850,"dev.ensemble3b"] <- abs(ci.terms[rows.0850,"dev.ensemble2"])
      
    }  
  }
  
  
  summary(ci.terms)  
  
  
  mod.agree <- aggregate(ci.terms[ci.terms$data.type=="Model",c("mean.cent.deriv", "dev.ensemble", "cent.deriv2", "dev.ensemble2", "dev.ensemble3a", "dev.ensemble3b")],
                         by=ci.terms[ci.terms$data.type=="Model", c("Effect", "Extent", "Model")],
                         FUN=mean, na.rm=T)
  mod.agree[,paste0(c("mean.cent.deriv", "dev.ensemble", "cent.deriv2", "dev.ensemble2", "dev.ensemble3a", "dev.ensemble3b"), ".sd")] <- aggregate(ci.terms[ci.terms$data.type=="Model",c("mean.cent.deriv", "dev.ensemble", "cent.deriv2", "dev.ensemble2", "dev.ensemble3a", "dev.ensemble3b")],
                                                                               by=ci.terms[ci.terms$data.type=="Model", c("Effect", "Extent", "Model")],
                                                                               FUN=sd, na.rm=T)[,c("mean.cent.deriv", "dev.ensemble", "cent.deriv2", "dev.ensemble2", "dev.ensemble3a", "dev.ensemble3b")]
  summary(mod.agree)
  
  # Looking at the mean deviation from the ensemble mean
  # 3b finds the model mean analytically rather than through the gam
  ext.co2     <- lm(dev.ensemble3b ~ Extent, data=mod.agree[mod.agree$Effect=="CO2"    ,])
  ext.precipf <- lm(dev.ensemble3b ~ Extent, data=mod.agree[mod.agree$Effect=="precipf",])
  ext.tair    <- lm(dev.ensemble3b ~ Extent, data=mod.agree[mod.agree$Effect=="tair"   ,])
  summary(ext.co2)
  summary(ext.precipf)
  summary(ext.tair)

  ext.co2.sd     <- lm(dev.ensemble.3b.sd ~ Extent, data=mod.agree[mod.agree$Effect=="CO2"    ,])
  ext.precipf.sd <- lm(dev.ensemble.3b.sd ~ Extent, data=mod.agree[mod.agree$Effect=="precipf",])
  ext.tair.sd    <- lm(dev.ensemble.3b.sd ~ Extent, data=mod.agree[mod.agree$Effect=="tair"   ,])
  summary(ext.co2.sd)
  summary(ext.precipf.sd)
  summary(ext.tair.sd)
}
# -----------

# -----------
# Comparing Models & tree rings (1900-2010, 1980-2010)
# -----------
{
mod.deriv2 <- aggregate(ci.terms[,c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")],
                       by=ci.terms[,c("Effect", "Extent", "data.type", "Y.type", "Model", "veg.scheme", "Evergreen.sd","fire.scheme", "Biomass.rel.sd"),],
                       FUN=mean, na.rm=T)
mod.deriv2[,paste0(c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2"), ".sd")] <- aggregate(ci.terms[, c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")],
                                                                                              by=ci.terms[,c("Effect", "Extent", "data.type", "Y.type", "Model", "veg.scheme", "Evergreen.sd","fire.scheme", "Biomass.rel.sd"),],
                                                                                              FUN=sd, na.rm=T)[,c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")]
mod.deriv2 <- mod.deriv2[!mod.deriv2$Effect=="Biomass",]
summary(mod.deriv2)

ens.deriv2 <- aggregate(mod.deriv2[,c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")], 
                       by=mod.deriv2[,c("Effect", "Extent", "Y.type", "data.type")],
                       FUN=mean, na.rm=T)
ens.deriv2[,paste0(c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2"), ".sd")] <- aggregate(mod.deriv2[,c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")], 
                                                                              by=mod.deriv2[,c("Effect", "Extent", "Y.type", "data.type")],
                                                                              FUN=sd, na.rm=T)[,c("mean.rel.cent", "mean.cent.deriv", "cent.deriv2")]
ens.deriv2

# 1980-2010
ens.deriv2[ens.deriv2$Extent=="1980-2010" & ens.deriv2$Effect=="tair",]
ens.deriv2[ens.deriv2$Extent=="1980-2010" & ens.deriv2$Effect=="precipf",]
ens.deriv2[ens.deriv2$Extent=="1980-2010" & ens.deriv2$Effect=="CO2",]
}
# -----------


ci.terms.agg <- aggregate(ci.terms[,c("mean.cent.deriv", "mean.rel.cent", "cent.deriv2", "Y", "Evergreen.sd", "Biomass.rel.sd")], 
                          by=ci.terms[,c("Effect", "Extent", "Model", "veg.scheme", "fire.scheme", "data.type", "Y.type")],
                          FUN=mean, na.rm=T)
summary(ci.terms.agg)
# -----------
# Analyzing Precipitation Sensitivity
# Note: looking at the absolute value of the slope so that smaller numebrs are closer to 0
# -----------
{
# Correlation with Vegetation Scheme
precipf.veg.lme <- lme(abs(cent.deriv2) ~ veg.scheme*(Extent-1) - veg.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",], na.action=na.omit)
summary(precipf.veg.lme)

# precipf.veg.lm2 <- lm(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
precipf.veg.lm2 <- lm(abs(cent.deriv2) ~ veg.scheme*(Extent-1) - veg.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="precipf",])
summary(precipf.veg.lm2)

# Correlation with Composition Variability
precipf.evg.lme <- lme(abs(cent.deriv2) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",], na.action=na.omit)
summary(precipf.evg.lme)

# precipf.evg.lm2 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
# precipf.evg.lm2 <- lm(abs(cent.deriv2) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
precipf.evg.lm2 <- lm(abs(cent.deriv2) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="precipf",])
summary(precipf.evg.lm2)

# Correlation with Fire 
precipf.fire.lme <- lme(abs(mean.cent.deriv) ~ fire.scheme*(Extent-1) - fire.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",], na.action=na.omit)
summary(precipf.fire.lme)

# precipf.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
# precipf.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
precipf.fire.lm2 <- lm(abs(cent.deriv2) ~ fire.scheme*(Extent-1) - fire.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="precipf",])
summary(precipf.fire.lm2) 

# Correlation with Biomass
precipf.bm.lme <- lme(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",], na.action=na.omit)
summary(precipf.bm.lme)

# precipf.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="precipf",])
# precipf.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf",])
precipf.bm.lm2 <- lm(abs(cent.deriv2) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="precipf",])
summary(precipf.bm.lm2)
}
# -----------

# -----------
# Analyzing CO2 Sensitivity
# Note: looking at the absolute value of the slope so that smaller numebrs are closer to 0
# -----------
{
  # Correlation with Vegetation Scheme
  co2.veg.lme <- lme(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",], na.action=na.omit)
  summary(co2.veg.lme)
  
  #   co2.veg.lm2 <- lm(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  # co2.veg.lm2 <- lm(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2",])
  co2.veg.lm2 <- lm(abs(cent.deriv2) ~ veg.scheme*(Extent-1) - veg.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="CO2",])
  summary(co2.veg.lm2)
  
  # Correlation with Composition Variability
  co2.evg.lme <- lme(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",], na.action=na.omit)
  summary(co2.evg.lme)

#   co2.evg.lm2 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  # co2.evg.lm2 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2",])
  co2.evg.lm2 <- lm(abs(cent.deriv2) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="CO2",])
  summary(co2.evg.lm2)
  
  # Correlation with Fire 
  co2.fire.lme <- lme(abs(mean.cent.deriv) ~ fire.scheme*(Extent-1) - fire.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",], na.action=na.omit)
  summary(co2.fire.lme)
  
  #   co2.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  # co2.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2",])
  co2.fire.lm2 <- lm(abs(cent.deriv2) ~ fire.scheme*(Extent-1) - fire.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="CO2",])
  summary(co2.fire.lm2) 
  
  # Correlation with Biomass
  co2.bm.lme <- lme(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",], na.action=na.omit)
  summary(co2.bm.lme)
  
  #   co2.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="CO2",])
  # co2.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2",])
  co2.bm.lm2 <- lm(abs(cent.deriv2) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="CO2",])
  summary(co2.bm.lm2)
}
# -----------


# -----------
# Analyzing Temperature Sensitivity
# Note: looking at the absolute value of the slope so that smaller numebrs are closer to 0
# -----------
{ 
  # Correlation with Vegetation Scheme
  tair.veg.lme <- lme(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="tair",], na.action=na.omit)
  summary(tair.veg.lme)
  
  #   tair.veg.lm2 <- lm(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  # tair.veg.lm2 <- lm(abs(mean.cent.deriv) ~ veg.scheme*(Extent-1) - veg.scheme, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair",])
  tair.veg.lm2 <- lm(abs(cent.deriv2) ~ veg.scheme*(Extent-1) - veg.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="tair",])
  summary(tair.veg.lm2)
  
  # Correlation with Composition Variability
  tair.evg.lme <- lme(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="tair",], na.action=na.omit)
  summary(tair.evg.lme)
  
  #   tair.evg.lm2 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  # tair.evg.lm2 <- lm(abs(mean.cent.deriv) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair",])
  tair.evg.lm2 <- lm(abs(cent.deriv2) ~ Evergreen.sd*(Extent-1) - Evergreen.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="tair",])
  summary(tair.evg.lm2)
  
  # Correlation with Fire 
  tair.fire.lme <- lme(abs(mean.cent.deriv) ~ fire.scheme*(Extent-1) - fire.scheme, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="tair",], na.action=na.omit)
  summary(tair.fire.lme)
  
  #   tair.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms[!ci.terms$Quantile=="Other"  & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  # tair.fire.lm2 <- lm(abs(mean.cent.deriv) ~ (Extent-1)*fire.scheme - fire.scheme, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair",])
  tair.fire.lm2 <- lm(abs(cent.deriv2) ~ fire.scheme*(Extent-1) - fire.scheme, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="tair",])
  summary(tair.fire.lm2) 
  
  # Correlation with Biomass
  tair.bm.lme <- lme(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, random=list(Model=~1, x=~1), data=ci.terms[ci.terms$data.type=="Model" & ci.terms$Effect=="tair",], na.action=na.omit)
  summary(tair.bm.lme)
  
  #   tair.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms[!ci.terms$Quantile=="Other"  & !is.na(ci.terms$mean.rel.cent) & ci.terms$data.type=="Model" & ci.terms$Effect=="tair",])
  # tair.bm.lm2 <- lm(abs(mean.cent.deriv) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair",])  
  tair.bm.lm2 <- lm(abs(cent.deriv2) ~ Biomass.rel.sd*(Extent-1) - Biomass.rel.sd, data=mod.deriv2[mod.deriv2$data.type=="Model" & mod.deriv2$Effect=="tair",])
  summary(tair.bm.lm2)
}
# -----------


# -----------
# Remaking Table 3 with the derivatives
# -----------
{
# vegetation scheme
veg.precip <- data.frame(char="veg.scheme (static)", 
                         effect="precip", 
                         extent=c("1980-2010","1901-2010", "850-2010"),
                         estimate=summary(precipf.veg.lm2)$coefficients[4:6,1],
                         std.err = summary(precipf.veg.lm2)$coefficients[4:6,2],
                         t.stat  = summary(precipf.veg.lm2)$coefficients[4:6,3],
                         p.value = summary(precipf.veg.lm2)$coefficients[4:6,4])
veg.tair <- data.frame(char="veg.scheme (static)", 
                       effect="tair", 
                       extent=c("1980-2010","1901-2010", "850-2010"),
                       estimate=summary(tair.veg.lm2)$coefficients[4:6,1],
                       std.err = summary(tair.veg.lm2)$coefficients[4:6,2],
                       t.stat  = summary(tair.veg.lm2)$coefficients[4:6,3],
                       p.value = summary(tair.veg.lm2)$coefficients[4:6,4])
veg.co2 <- data.frame(char="veg.scheme (static)", 
                      effect="co2", 
                      extent=c("1980-2010","1901-2010", "850-2010"),
                      estimate=summary(co2.veg.lm2)$coefficients[4:6,1],
                      std.err = summary(co2.veg.lm2)$coefficients[4:6,2],
                      t.stat  = summary(co2.veg.lm2)$coefficients[4:6,3],
                      p.value = summary(co2.veg.lm2)$coefficients[4:6,4])
veg.stats <- rbind(veg.tair, veg.precip, veg.co2)

# Composition variability
evg.precip <- data.frame(char="evg.sd", 
                         effect="precip", 
                         extent=c("1980-2010","1901-2010", "850-2010"),
                         estimate=summary(precipf.evg.lm2)$coefficients[4:6,1],
                         std.err = summary(precipf.evg.lm2)$coefficients[4:6,2],
                         t.stat  = summary(precipf.evg.lm2)$coefficients[4:6,3],
                         p.value = summary(precipf.evg.lm2)$coefficients[4:6,4])
evg.tair <- data.frame(char="evg.sd", 
                       effect="tair", 
                       extent=c("1980-2010","1901-2010", "850-2010"),
                       estimate=summary(tair.evg.lm2)$coefficients[4:6,1],
                       std.err = summary(tair.evg.lm2)$coefficients[4:6,2],
                       t.stat  = summary(tair.evg.lm2)$coefficients[4:6,3],
                       p.value = summary(tair.evg.lm2)$coefficients[4:6,4])
evg.co2 <- data.frame(char="evg.sd", 
                      effect="co2", 
                      extent=c("1980-2010","1901-2010", "850-2010"),
                      estimate=summary(co2.evg.lm2)$coefficients[4:6,1],
                      std.err = summary(co2.evg.lm2)$coefficients[4:6,2],
                      t.stat  = summary(co2.evg.lm2)$coefficients[4:6,3],
                      p.value = summary(co2.evg.lm2)$coefficients[4:6,4])

evg.stats <- rbind(evg.tair, evg.precip, evg.co2)

# Fire
fire.precip <- data.frame(char="fire (yes)", 
                          effect="precip", 
                          extent=c("1980-2010","1901-2010", "850-2010"),
                          estimate=summary(precipf.fire.lm2)$coefficients[4:6,1],
                          std.err = summary(precipf.fire.lm2)$coefficients[4:6,2],
                          t.stat  = summary(precipf.fire.lm2)$coefficients[4:6,3],
                          p.value = summary(precipf.fire.lm2)$coefficients[4:6,4])
fire.tair <- data.frame(char="fire (yes)", 
                        effect="tair", 
                        extent=c("1980-2010","1901-2010", "850-2010"),
                        estimate=summary(tair.fire.lm2)$coefficients[4:6,1],
                        std.err = summary(tair.fire.lm2)$coefficients[4:6,2],
                        t.stat  = summary(tair.fire.lm2)$coefficients[4:6,3],
                        p.value = summary(tair.fire.lm2)$coefficients[4:6,4])
fire.co2 <- data.frame(char="fire (yes)", 
                       effect="co2", 
                       extent=c("1980-2010","1901-2010", "850-2010"),
                       estimate=summary(co2.fire.lm2)$coefficients[4:6,1],
                       std.err = summary(co2.fire.lm2)$coefficients[4:6,2],
                       t.stat  = summary(co2.fire.lm2)$coefficients[4:6,3],
                       p.value = summary(co2.fire.lm2)$coefficients[4:6,4])

fire.stats <- rbind(fire.tair, fire.precip, fire.co2)

# Biomass
bm.precip <- data.frame(char="bm.sd", 
                        effect="precip", 
                        extent=c("1980-2010","1901-2010", "850-2010"),
                        estimate=summary(precipf.bm.lm2)$coefficients[4:6,1],
                        std.err = summary(precipf.bm.lm2)$coefficients[4:6,2],
                        t.stat  = summary(precipf.bm.lm2)$coefficients[4:6,3],
                        p.value = summary(precipf.bm.lm2)$coefficients[4:6,4])
bm.tair <- data.frame(char="bm.sd", 
                      effect="tair", 
                      extent=c("1980-2010","1901-2010", "850-2010"),
                      estimate=summary(tair.bm.lm2)$coefficients[4:6,1],
                      std.err = summary(tair.bm.lm2)$coefficients[4:6,2],
                      t.stat  = summary(tair.bm.lm2)$coefficients[4:6,3],
                      p.value = summary(tair.bm.lm2)$coefficients[4:6,4])
bm.co2 <- data.frame(char="bm.sd", 
                     effect="co2", 
                     extent=c("1980-2010","1901-2010", "850-2010"),
                     estimate=summary(co2.bm.lm2)$coefficients[4:6,1],
                     std.err = summary(co2.bm.lm2)$coefficients[4:6,2],
                     t.stat  = summary(co2.bm.lm2)$coefficients[4:6,3],
                     p.value = summary(co2.bm.lm2)$coefficients[4:6,4])

bm.stats <- rbind(bm.tair, bm.precip, bm.co2)




model.stats <- rbind(veg.stats, evg.stats, fire.stats, bm.stats)
model.stats$effect.ext <- as.factor(paste(model.stats$effect, model.stats$extent, sep="."))
model.stats[,c("estimate", "std.err", "t.stat")] <- round(model.stats[,c("estimate", "std.err", "t.stat")]*100, 2)
model.stats[,"p.value"] <- round(model.stats[,"p.value"], 2)
model.stats$sig     <- ifelse(model.stats$p.value<0.05, "*", "")
model.stats$est.err <- paste(model.stats$estimate, "+/-", model.stats$std.err, " ", model.stats$sig)
model.stats
write.csv(model.stats, file.path("Data/analyses/analysis_TempExtent", "ModelCharStats_full.csv"), row.names=F)


library(reshape2)
table3 <- merge(unique(model.stats$char), unique(model.stats$effect))
names(table3) <- c("Character", "Effect")
table3 <- table3[,c(2,1)]
table3

for(i in unique(table3$Character)){
  for(j in unique(table3$Effect)){
    for(e in unique(model.stats$extent)){
      table3[table3$Character==i & table3$Effect==j, e] <- model.stats[model.stats$char==i & model.stats$effect==j & model.stats$extent==e, "est.err"]
    }
  }
}
table3

write.csv(table3, file.path("Manuscript/Submission3/", "Table3_ModelCharStats.csv"), row.names=F)
}
# -----------

}
# ----------------------------------------


# ----------------------------------------
# Stats on shifting influence of factors through time (last paragrph of results)
# from 7_analysis_response_baseline_ensembles.R
# ----------------------------------------
{
load("Data/analyses/analysis_baseline/post-process_baseline.RData")
load("Data/EcosysData.Rdata")

# ---------------------
# 4.a. Aggregate to region level
#  -- Note: going to site first for tree ring products
# ---------------------
{
  summary(dat.ecosys)
  summary(wt.terms)
  
  fac.ecosys <- c("Y", "Y.10", "Y.rel", "Y.rel.10", "tair", "precipf", "CO2", "Time")
  dat.ecosys2 <- aggregate(dat.ecosys[,fac.ecosys],  
                           by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type", "Site", "Year")],
                           FUN=mean, na.rm=T)
  summary(dat.ecosys2)
  
  fac.wt <- c("fit.full", "fit.tair", "fit.precipf", "fit.CO2", "fit.tair.10", "fit.precipf.10", "fit.CO2.10", "fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel", "fit.tair.rel.10", "fit.precipf.rel.10", "fit.CO2.rel.10", "weight.tair", "weight.precipf", "weight.CO2", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
  wt.terms2 <- aggregate(wt.terms[,fac.wt],  
                         by=wt.terms[,c("Model","Site", "Year")],
                         FUN=mean, na.rm=T)
  summary(wt.terms2)
  
  
  # ----------
  # Adjusting CO2 Effect
  # ----------
  # Note: because the gam makes the smoother cross 0 at the MEAN CO2 (which is in the 1800s), 
  # it's saying the region is pretty CO2-limited at periods where that doesn't really make 
  # sense, so we're going to off relativize it to whatever the starting point for the run is
  # ----------
{
  for(m in unique(wt.terms2$Model)){
    yr1       <- min(wt.terms2[wt.terms2$Model==m, "Year"]) # find the minimum year
    yr2       <- min(wt.terms2[wt.terms2$Model==m & !is.na(wt.terms2$weight.CO2.10), "Year"]) # find the minimum year
    co2.base <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr1+5), "fit.CO2"], na.rm=T) # mean of the 5 years around the starting 
    co2.base.10 <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr2+5),"fit.CO2.10"],na.rm=T) # mean of the 5 years around the starting point
    
    co2.rel.base <- 1-mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr1+5), "fit.CO2.rel"], na.rm=T) # mean of the 5 years around the starting 
    co2.rel.base.10 <- 1-mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Year<=(yr2+5),"fit.CO2.rel.10"],na.rm=T) # mean of the 5 years around the starting point
    
    wt.terms2[wt.terms2$Model==m, "fit.CO2.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2"] - co2.base
    wt.terms2[wt.terms2$Model==m, "fit.CO2.10.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.10"] - co2.base.10
    wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.rel"] + co2.rel.base # The adjustment needs to be relative to 1
    wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.10.adj"] <- wt.terms2[wt.terms2$Model==m, "fit.CO2.rel.10"] + co2.rel.base.10 # The adjustment needs ot be relative to 1
  }
  
  wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")] <- abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")])/rowSums(abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")]))
  wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")] <- abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")])/rowSums(abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")]))
  
  summary(wt.terms2)
  
  wt.check <- rowSums(wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")])
  wt.check.10 <- rowSums(wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")])
  summary(wt.check)
  summary(wt.check.10)
}
# ----------

dat.ecosys2 <- merge(dat.ecosys2, wt.terms2, all.x=T, all.y=T)
summary(dat.ecosys2)

fac.wt <- c(fac.wt, "fit.CO2.rel.adj", "fit.CO2.rel.10.adj")
fac.wt.2 <- c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")
fac.wt.ci <- c("fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel.adj", "fit.tair.rel.10", "fit.precipf.rel.10", "fit.CO2.rel.10.adj")
factors.agg <- c(fac.ecosys, fac.wt.2, fac.wt.ci)

dat.region                             <- aggregate(dat.ecosys2[,factors.agg], 
                                                    by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                    FUN=mean, na.rm=T)
dat.region[,paste0(factors.agg, ".lo")] <- aggregate(dat.ecosys2[,factors.agg], 
                                                     by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                     FUN=quantile, 0.025, na.rm=T)[,factors.agg]
dat.region[,paste0(factors.agg, ".hi")] <- aggregate(dat.ecosys2[,factors.agg], 
                                                     by=dat.ecosys2[,c("Model", "Model.Order", "Y.type", "data.type", "Year")], 
                                                     FUN=quantile, 0.975, na.rm=T)[,factors.agg]
summary(dat.region)


# Ensemble-level stats
# Aggregate models, them merge in tree ring region level
dev.region <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg], 
                        by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                        FUN=mean, na.rm=T)
dev.region[,paste0(factors.agg, ".lo")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg], 
                                                     by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                                     FUN=quantile, 0.025, na.rm=T)[,factors.agg]
dev.region[,paste0(factors.agg, ".hi")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.agg],
                                                     by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                                     FUN=quantile, 0.975, na.rm=T)[,factors.agg]

# Adding in the standard deviation for a few other variables we want to use to describe the models through time
factors.sd <- c("Y.rel", "fit.tair.rel", "fit.precipf.rel", "fit.CO2.rel.adj")
dev.region[,paste0(factors.sd, ".sd")] <- aggregate(dat.region[dat.region$data.type=="Model",factors.sd],
                                                    by= dat.region[dat.region$data.type=="Model",c("Y.type", "data.type", "Year")], 
                                                    FUN=sd, na.rm=T)[,factors.sd]
summary(dev.region)

dev.region <- merge(dev.region, dat.region[dat.region$data.type=="Tree Rings",], all.x=T, all.y=T)

dev.region[is.na(dev.region$weight.tair.10.adj   ),"weight.tair.10.adj"   ] <- 0
dev.region[is.na(dev.region$weight.precipf.10.adj),"weight.precipf.10.adj"] <- 0
dev.region[is.na(dev.region$weight.CO2.10.adj    ),"weight.CO2.10.adj"    ] <- 0
}
# ---------------------

# ---------------------
# 4.c. Some summary statistics on model agreement through time
# ---------------------
{
  summary(dat.region)
  
  # --------
  # Model agreement in terms of relative NPP
  # --------
{
  # pre-1900
  summary(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"])
  mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"])
  
  # 1980-2010
  summary(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"])
  mean(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"]); sd(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "Y.rel.sd"])
  
  # 1900-2010
  summary(dev.region[dev.region$Year>=1900 & dev.region$data.type=="Model", "Y.rel.sd"])
  mean(dev.region[dev.region$Year>=1901 & dev.region$data.type=="Model", "Y.rel.sd"]); sd(dev.region[dev.region$Year>=1901 & dev.region$data.type=="Model", "Y.rel.sd"])
  
  
  t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "Y.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "Y.rel.sd"])
}
# --------

# --------
# Comparing Factor weights and the shift in their agreement 
# --------
{
  # tair
  t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
  mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
  mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.tair.rel.sd"])
  apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.tair.rel.lo", "fit.tair.rel.hi", "fit.tair.rel")],2, mean)
  
  # precipf
  t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
  mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
  mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.precipf.rel.sd"])
  apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.precipf.rel.lo", "fit.precipf.rel.hi", "fit.precipf.rel")],2, mean)
  
  # CO2
  t.test(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"], dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
  summary(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
  mean(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"]); sd(dev.region[dev.region$Year<1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
  mean(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"]); sd(dev.region[dev.region$Year>1900 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
  mean(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"]); sd(dev.region[dev.region$Year>1980 & dev.region$data.type=="Model", "fit.CO2.rel.adj.sd"])
  summary(dat.region[dat.region$Year>=1980 & dat.region$data.type=="Model", c("Model", "Year", "fit.CO2.rel.adj")])
  summary(dat.region[!dat.region$Model=="linkages" & dat.region$Year>=1980 & dat.region$data.type=="Model", c("Model", "Year", "fit.CO2.rel.adj")])
  summary(dat.region[dat.region$Model=="linkages" & dat.region$Year>=1980 & dat.region$data.type=="Model", c("Model", "Year", "fit.CO2.rel.adj")])
  summary(wt.terms)
  summary(wt.terms2)
  
  apply(dev.region[dev.region$Year>=1980 & dev.region$data.type=="Model", c("fit.CO2.rel.adj.lo", "fit.CO2.rel.adj.hi", "fit.CO2.rel.adj")],2, mean)
}
# --------
}
# ---------------------

}
# ----------------------------------------
