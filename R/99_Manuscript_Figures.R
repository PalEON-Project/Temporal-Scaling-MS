# # ----------------------------------------
# Objective: Put all of the Figures in the manuscript in one place
# Christy Rollinson, crollinson@gmail.com
# Date Created: 5 May 2016
#
# This script pulls on numbers/analyses originally performed in scripts:
#   - 7_analysis_response_baseline_ensembles.R
#   - 8_analysis_sensitivity_TemporalExtent.R
# ----------------------------------------
rm(list=ls())

setwd("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Temporal-Scaling")
load("Data/EcosysData.Rdata")


# ----------------------------------------
# Load libraries we may need
library(ggplot2); library(grid); library(scales); 
library(gridExtra)
library(car); library(zoo)
# ----------------------------------------

# ----------------------------------------
# Figure 1: Model NPP variability across space & time
# ----------------------------------------
{
  load("Data/analyses/analysis_baseline/post-process_baseline.RData")
  
  dat.ecosys$Site.Order <- recode(dat.ecosys$Site, "'PDL'='1'; 'PBL'='2'; 'PUN'='3'; 'PMB'='4'; 'PHA'='5'; 'PHO'='6'")
  levels(dat.ecosys$Site.Order) <- c("Demming", "Billy's", "UNDERC", "Minden", "Harvard", "Howland")
  
  dat.ecosys$Model.Order <- factor(dat.ecosys$Model.Order, levels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "SiBCASA", "Tree Ring NPP", "Tree Ring RW"))
  models.use <- unique(dat.ecosys[,"Model.Order"])
  
  model.colors$Model.Order <- factor(model.colors$Model.Order, levels=c("CLM-BGC", "CLM-CN", "ED2", "ED2-LU", "JULES-STATIC", "JULES-TRIFFID", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "SiBCASA", "CLM4.5", "CLM4.5-DGVM"))
  model.colors <- model.colors[order(model.colors$Model.Order),]
  colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray30"))
  
  png(file.path("Figures/analyses/analysis_baseline", "Fig1_NPP_Raw_AllSites_0850-2010_Simple_Models.png"), height=5, width=8, units="in", res=180)
  {
    ggplot(data=dat.ecosys[!dat.ecosys$Model %in% c("TreeRingRW", "TreeRingBAI", "TreeRingNPP"),])  + 
      facet_wrap(~Site.Order) +
      geom_line(aes(x=Year, y=Y, color=Model.Order), size=0.1, alpha=0.3) + 
      geom_line(aes(x=Year, y=Y.10, color=Model.Order), size=0.75, alpha=1) + 
      scale_x_continuous(limits=c(0850, 2010), expand=c(0,0), breaks=seq(min(dat.ecosys$Year), max(dat.ecosys$Year), by=250)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_manual(values=c("black", "gray50")) +
      scale_color_manual(values=colors.use) +
      labs(color="Model", x="Year", y=expression(bold(paste("NPP (Mg C ha"^"-1"," yr"^"-1",")")))) +
      guides(col=guide_legend(nrow=2, title="Model"), fill=F) +
      theme(legend.position="top") +
      # 	theme(plot.title=element_text(face="bold", size=rel(3))) + 
      theme(legend.text=element_text(size=8), 
            legend.title=element_text(size=10),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines"),
            legend.position=c(0.5, 0.9)) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(), 
            axis.text.x=element_text(angle=0, color="black", size=10), 
            axis.text.y=element_text(color="black", size=10), 
            axis.title.x=element_text(face="bold", vjust=-0.5, size=12),  
            axis.title.y=element_text(face="bold", vjust=1, size=12))
  }
  dev.off()
  
}
# ----------------------------------------


# ----------------------------------------
# Figure 2: NPP sensitivity to climate and CO2 at different temporal scales
# ----------------------------------------
{
  load(file.path("Data/analyses/analysis_TempExtent", "post-process_TempExtent.RData"))
  
  models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
  colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))
  
  # -----------------
  # Creating a cheat data frame that lets values go off the graph
  # -----------------
  {
  ci.terms.graph <- ci.terms
  ci.terms.graph[!is.na(ci.terms.graph$mean.rel) & ci.terms.graph$mean.rel<(-0.65),"mean.rel"] <- NA 
  ci.terms.graph[!is.na(ci.terms.graph$lwr.rel) & ci.terms.graph$lwr.rel<(-0.65),"lwr.rel"] <- -0.65 
  ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$upr.rel<(-0.65),"upr.rel"] <- -0.65 
  ci.terms.graph[which(ci.terms.graph$mean.rel>1.00),"mean.rel"] <- NA 
  ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$lwr.rel>(1.00),"lwr.rel"] <- 1.00
  ci.terms.graph[!is.na(ci.terms.graph$upr.rel) & ci.terms.graph$upr.rel>(1.00),"upr.rel"] <- 1.00 
  ci.terms.graph[ci.terms.graph$Effect=="tair", "x"] <- ci.terms.graph[ci.terms.graph$Effect=="tair", "x"]-273.15
  
  ci.terms.graph <- merge(ci.terms.graph, models.df, all.x=T, all.y=F)
  summary(ci.terms.graph)
  
  # Grouping the kind and source of the data
  ci.terms.graph$Y.type <- as.factor(ifelse(ci.terms.graph$Model=="TreeRingRW", "RW", "NPP"))
  ci.terms.graph$data.type <- as.factor(ifelse(substr(ci.terms.graph$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
  summary(ci.terms.graph)
  # summary(ci.terms.graph[ci.terms.graph$Model=="linkages",])
  
  # Creating a mask for values outside of the model drivers for that time period
  for(e in unique(ci.terms.graph$Extent)){
    yr <- as.numeric(strsplit(paste(e), "-")[[1]][1])
    
    tair    <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"tair"   ], na.rm=T) - 273.15
    precipf <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"precipf"], na.rm=T)
    co2     <- range(dat.ecosys2[dat.ecosys2$Model=="ed2" & dat.ecosys2$Year>=yr,"CO2"    ], na.rm=T)
    
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="tair"   , "line.min"] <- tair   [1]
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="precipf", "line.min"] <- precipf[1]
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="CO2"    , "line.min"] <- co2    [1]
    
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="tair"   , "line.max"] <- tair   [2]
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="precipf", "line.max"] <- precipf[2]
    ci.terms.graph[ci.terms.graph$Extent==e & ci.terms.graph$Effect=="CO2"    , "line.max"] <- co2    [2]
  }
  ci.terms.graph$x.min <- ifelse(ci.terms.graph$x<ci.terms.graph$line.min, ci.terms.graph$x, ci.terms.graph$line.min)
  ci.terms.graph$x.max <- ifelse(ci.terms.graph$x>ci.terms.graph$line.max, ci.terms.graph$x, ci.terms.graph$line.max)
  ci.terms.graph$mask.min <- min(ci.terms.graph$lwr.rel, na.rm=T)
  ci.terms.graph$mask.max <- max(ci.terms.graph$upr.rel, na.rm=T)
  
  # Playing with the extent labels a bit so that "850-2010" is "all data" and "1901-2010" is left alone
  ci.terms.graph$Extent2 <- as.factor(ifelse(ci.terms.graph$Extent=="850-2010" | 
                                               (ci.terms.graph$Extent=="1980-2010" & ci.terms.graph$Model=="TreeRingNPP") |
                                               (ci.terms.graph$Extent=="1901-2010" & ci.terms.graph$Model=="TreeRingRW")
                                             , "All Data", paste(ci.terms.graph$Extent)))
  
  ci.terms.graph[ci.terms.graph$Extent2=="All Data", c("x.min", "x.max", "mask.min", "mask.max", "line.min", "line.max")] <- NA
  summary(ci.terms.graph)
  
  
  tree.rings.1901 <- ci.terms.graph[ci.terms.graph$Model=="TreeRingRW" & ci.terms.graph$Extent=="1901-2010",]
  tree.rings.1901$Extent2 <- as.factor("1901-2010")
  summary(tree.rings.1901)
  
  tree.rings.npp <- ci.terms.graph[ci.terms.graph$Model=="TreeRingNPP" & ci.terms.graph$Extent=="1980-2010",]
  tree.rings.npp$Extent2 <- as.factor("1980-2010")
  summary(tree.rings.npp)
  
  
  ci.terms.graph <- rbind(ci.terms.graph, tree.rings.1901, tree.rings.npp)
  ci.terms.graph$Extent3 <- recode(ci.terms.graph$Extent2, "'All Data'='0'; '1901-2010'='1'; '1980-2010'='2'")
  levels(ci.terms.graph$Extent3) <- c("0850-2010", "1901-2010", "1980-2010")
  summary(ci.terms.graph)
  
  ci.terms.graph <- ci.terms.graph[!(ci.terms.graph$Extent3=="0850-2010" & ci.terms.graph$data.type=="Tree Rings"),]
  levels(ci.terms.graph$Effect) <- c("Temperature", "Precipitation", "CO2", "Biomass")
  summary(ci.terms.graph)
  }
  # -----------------
  
  fig3.tair <- {
    ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "Temperature",]) + 
      facet_grid(Extent3~Effect, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
      geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
      # Lower Shaded Region
      geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.min), linetype="dashed") +
      # Upper Shaded Region
      geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.max), linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name=expression(bold(paste("Temperature ("^"o", "C)"))), breaks=c(10, 12.5, 15.0, 17.5)) +
      scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
      guides(fill=F, color=F, linetype=F) +
      scale_fill_manual(values=colors.use) +
      scale_color_manual(values=colors.use) +
      scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_blank()) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0, "lines"))  +
      theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.65,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"))
    
  }
  fig3.precip <- {
    ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "Precipitation",]) + 
      facet_grid(Extent3~Effect, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
      geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
      # Lower Shaded Region
      geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.min), linetype="dashed") +
      # Upper Shaded Region
      geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.max), linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name=expression(bold(paste("Precipitation (mm yr"^"-1", ")")))) +
      scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
      guides(fill=F, color=F, linetype=F) +
      scale_fill_manual(values=colors.use) +
      scale_color_manual(values=colors.use) +
      scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
      theme(strip.text.x=element_text(size=12, face="bold"),
            strip.text.y=element_blank()) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0, "lines"))  +
      theme(axis.text.y=element_blank(),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_blank(),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0.0,0.5,0.5), "lines"))
    
  }
  fig3.co2 <- {
    ggplot(data=ci.terms.graph[ci.terms.graph$Effect == "CO2",]) + 
      facet_grid(Extent3~Effect, scales="free_x") +
      geom_ribbon(aes(x=x, ymin=lwr.rel*100, ymax=upr.rel*100, fill=Model.Order), alpha=0.3) +
      geom_line(aes(x=x, y=mean.rel*100, color=Model.Order, linetype=Model.Order), size=1) +
      # Lower Shaded Region
      geom_ribbon(aes(x=x.min, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.min), linetype="dashed") +
      # Upper Shaded Region
      geom_ribbon(aes(x=x.max, ymin=mask.min*100, ymax=mask.max*100), alpha=0.3) +
      geom_vline(aes(xintercept=line.max), linetype="dashed") +
      scale_x_continuous(expand=c(0,0), name=expression(bold(paste("CO" ["2"], " (ppm)")))) +
      scale_y_continuous(name="NPP Contribution (% mean)", expand=c(0,0)) +
      guides(fill=guide_legend(title="Model"), 
             color=guide_legend(title="Model"), 
             linetype=guide_legend(title="Model")) +
      scale_fill_manual(values=colors.use) +
      scale_color_manual(values=colors.use) +
      scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
      theme(legend.title=element_text(size=12, face="bold"),
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
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0, "lines"))  +
      theme(axis.text.y=element_blank(),
            axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
            axis.title.y=element_blank(),  
            axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.79,0,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"))
    
  }
  
  png(file.path("Figures/analyses/analysis_TempExtent", "Fig3_Sensitivity_Rel_extent.png"), width=8, height=6, units="in", res=180)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3, widths=c(1.3,1,2))))
  print(fig3.tair  , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
  print(fig3.precip, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
  print(fig3.co2   , vp = viewport(layout.pos.row = 1, layout.pos.col=3))
  dev.off()
  
}
# ----------------------------------------


# ----------------------------------------
# Figure 3: Drivers of change in NPP over the past 300 years
# ----------------------------------------
{
  load("Data/analyses/analysis_baseline/post-process_baseline.RData")
 
  models.use <- unique(dat.ecosys[,"Model.Order"])
  colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.use, "color"]), "black", "gray30"))
  
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
  
  
  # --------
  # 4.b.1. Raw NPP by model
  # --------
  annotate.npp <- data.frame(Y.type="NPP", x=1715, y=15, label="a)")
  plot.npp <- {
    ggplot(data=dat.region) + 
      scale_x_continuous(limits=c(1700,2004), expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name=expression(bold(paste("NPP MgC HA"^"-1"," yr"^"-1")))) +
      facet_grid(Y.type~., scales="free_y", space="free") +
      geom_ribbon(aes(x=Year, ymin=Y.10.lo, ymax=Y.10.hi, fill=Model.Order), alpha=0.5) +
      geom_line(aes(x=Year, y=Y.10, color=Model.Order, linetype=Model.Order), size=1.5) +
      # annotate(geom="text", label="a)", x=1725, y=15) +
      geom_text(data=annotate.npp, aes(x=x, y=y, label=label), fontface="bold", size=8) +
      scale_fill_manual(values=colors.use) +
      scale_color_manual(values=colors.use) +
      scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
      guides(color=guide_legend(title="Model", nrow=3),
             fill =guide_legend(title="Model", nrow=3),
             linetype =guide_legend(title="Model", nrow=3)) +
      theme(legend.title=element_text(size=8, face="bold"),
            legend.text=element_text(size=8),
            legend.position=c(0.35, 0.83),
            legend.key=element_blank(),
            legend.key.size=unit(0.75, "lines"),
            legend.background=element_blank()) +
      theme(strip.text=element_text(size=11, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0, "lines"))  +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")), 
            axis.title.x=element_blank(),  
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(1.5,1,0.7,0.8), "lines"))
  }
  plot.npp <- ggplotGrob(plot.npp )
  plot.npp $heights[[5]] <- unit(5, "null")
  plot(plot.npp)
  # --------
  
  # --------
  # 4.b.2. NPP deviation by data stream
  # --------
  dat.plot.dev <- dev.region[dev.region$Year >=1700 & dev.region$Year<2004,] 
  dat.plot.dev$Mode.Plot <- factor(ifelse(dat.plot.dev$data.type=="Model", "Model", ifelse(dat.plot.dev$Y.type=="NPP", "NPPtr", "RW")), levels=c("Model", "NPPtr", "RW"))
  summary(dat.plot.dev)
  
  annotate.dev <- data.frame(Mode.Plot=c("Model"), x=1715, y=175, label="b)")
  plot.dev <- {
    ggplot() + 
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="% NPP") +
      facet_grid(Mode.Plot~., scales="free_y", space="free") +
      geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Model",], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.5) +
      geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Model",], aes(x=Year, y=Y.rel.10*100), size=2,
                color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.tair.10.adj"    ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.CO2.10.adj"     ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Model","weight.precipf.10.adj" ])), 
                size=3) +
      
      geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP",],
                  aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.5) +
      geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings"  & dat.plot.dev$Y.type=="NPP",], aes(x=Year, y=Y.rel.10*100), size=2,
                color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.tair.10.adj"    ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.CO2.10.adj"     ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="NPP","weight.precipf.10.adj" ])), 
                size=3) +
      
      geom_ribbon(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW",],
                  aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.5) +
      geom_line(data=dat.plot.dev[dat.plot.dev$data.type=="Tree Rings"  & dat.plot.dev$Y.type=="RW",], aes(x=Year, y=Y.rel.10*100), size=2,
                color=rgb(abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.tair.10.adj"    ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.CO2.10.adj"     ]),
                          abs(dat.plot.dev[dat.plot.dev$data.type=="Tree Rings" & dat.plot.dev$Y.type=="RW","weight.precipf.10.adj" ])), 
                size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      #   annotate(geom="text", label="b)", x=1725, y=180, size=14) +
      geom_text(data=annotate.dev, aes(x=x, y=y, label=label), fontface="bold", size=8) +
      scale_linetype_manual(values=c(rep("solid", length(colors.use)-1), "dashed")) +
      theme(legend.title=element_text(size=rel(1), face="bold"),
            legend.text=element_text(size=rel(1)),
            # legend.position=c(0.2, 0.18),
            legend.key=element_blank(),
            legend.key.size=unit(1.5, "lines")) +
      theme(strip.text=element_text(size=11, face="bold")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_rect(fill=NA, color="black", size=0.5), 
            panel.background=element_blank(),
            panel.margin.x=unit(0, "lines"),
            panel.margin.y=unit(0, "lines"))  +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")), 
            axis.title.x=element_blank(),  
            axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),
            axis.ticks.length=unit(-0.5, "lines")) +
      theme(plot.margin=unit(c(0,1,0.7,1.05), "lines"))
  }
  plot.dev <- ggplotGrob(plot.dev )
  plot.dev $heights[[3]] <- unit(150, "null")
  plot.dev $heights[[5]] <- unit(50, "null")
  plot.dev $heights[[7]] <- unit(50, "null")
  plot(plot.dev)
  # --------
  
  # --------
  # 4.b.3. Factor weights
  # --------
  {
    fit.stack <- stack(dev.region[,c("fit.CO2.rel.10.adj","fit.tair.rel.10","fit.precipf.rel.10")])
    names(fit.stack) <- c("fit.mean", "Effect")
    fit.stack$Effect <- as.factor(ifelse(substr(fit.stack$Effect,5,7)=="CO2", "CO2", ifelse(substr(fit.stack$Effect,5,8)=="tair", "tair", "precipf")))
    fit.stack[,c("data.type", "Y.type", "Year")] <- dev.region[,c("data.type", "Y.type", "Year")]
    fit.stack[,"ci.lo"] <- stack(dev.region[,c("fit.CO2.rel.10.adj.lo","fit.tair.rel.10.lo","fit.precipf.rel.10.lo")])[,1]
    fit.stack[,"ci.hi"] <- stack(dev.region[,c("fit.CO2.rel.10.adj.hi","fit.tair.rel.10.hi","fit.precipf.rel.10.hi")])[,1]
    
    fit.stack[!is.na(fit.stack$ci.lo) & fit.stack$ci.lo<0.44, "ci.lo"] <- 0.44
    fit.stack[!is.na(fit.stack$ci.lo) & fit.stack$ci.hi<0.44, "ci.hi"] <- 0.44
    fit.stack[!is.na(fit.stack$ci.hi) & fit.stack$ci.hi>1.95, "ci.hi"] <- 1.95
    fit.stack[!is.na(fit.stack$fit.mean) & (fit.stack$fit.mean<0.44 | fit.stack$fit.mean>1.95),"fit.mean"] <- NA
    summary(fit.stack)
    
    fit.stack$Mode.Plot <- factor(ifelse(fit.stack$data.type=="Model", "Model", ifelse(fit.stack$Y.type=="NPP", "NPPtr", "RW")), levels=c("Model", "NPPtr", "RW"))
    
    levels(fit.stack$Effect) <- c("CO2", "Precip", "Tair")

    annotate.wts <- data.frame(Mode.Plot="Model", x=1715, y=175, label="c)")
    
    plot.wts <- {
      ggplot(fit.stack[,]) + 
        scale_x_continuous(limits=c(1700,2010), expand=c(0,0), name="Year") +
        scale_y_continuous(expand=c(0,0), name="% NPP") +
        facet_grid(Mode.Plot~., scales="free_y", space="free") +
        geom_ribbon(aes(x=Year, ymin=ci.lo*100, ymax=ci.hi*100, fill=Effect), alpha=0.5) +
        geom_line(aes(x=Year, y=fit.mean*100, color=Effect), size=2) +
        geom_hline(yintercept=100, linetype="dashed") +
        #   annotate("text", label="c)", x=1725, y=180, size=14) +
        geom_text(data=annotate.wts, aes(x=x, y=y, label=label), fontface="bold", size=8) +
        scale_color_manual(values=c("green3", "blue", "red2")) +
        scale_fill_manual(values=c("green3", "blue", "red2")) +
        guides(color=guide_legend(nrow=1), 
               fill =guide_legend(nrow=1))+
        theme(legend.title=element_text(size=10, face="bold"),
              legend.text=element_text(size=10),
              legend.position=c(0.23, 0.85),
              legend.key=element_blank(),
              legend.key.size=unit(1, "lines")) +
        theme(strip.text=element_text(size=11, face="bold")) + 
        theme(axis.line=element_line(color="black", size=0.5), 
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_rect(fill=NA, color="black", size=0.5), 
              panel.background=element_blank(),
              panel.margin.x=unit(0, "lines"),
              panel.margin.y=unit(0, "lines"))  +
        theme(axis.text.x=element_text(color="black", size=12, margin=unit(c(1.5,0,0,0), "lines")),
              axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")), 
              axis.title.x=element_text(size=12, face="bold"),  
              axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),
              axis.ticks.length=unit(-0.5, "lines")) +
        theme(plot.margin=unit(c(0,1,0.5,1.05), "lines"))
    }
    plot.wts <- ggplotGrob(plot.wts )
    plot.wts $heights[[3]] <- unit(150, "null")
    plot.wts $heights[[5]] <- unit(50, "null")
    plot.wts $heights[[7]] <- unit(50, "null")
    plot(plot.wts)
  }
  # --------
  
  # --------
  # Putting NPP & change through time in context
  # --------
  png(file.path("Figures/analyses/analysis_baseline", "Fig3_NPP_Dev_1700-2010_NPP_Rel_Weight.png"), height=8, width=8, units="in", res=180)
  grid.arrange(plot.npp, plot.dev, plot.wts, ncol=1)
  dev.off()
  
  
  # --------

}
# ----------------------------------------


# ----------------------------------------
# Supplemental Figures
# ----------------------------------------
{
  load(file.path("Data/analyses/analysis_TempExtent", "post-process_TempExtent.RData"))
  
  models.df <- data.frame(Model=unique(dat.ecosys[,"Model"]), Model.Order=unique(dat.ecosys[,"Model.Order"]))
  colors.use <- as.vector(c(paste(model.colors[model.colors$Model.Order %in% models.df$Model.Order, "color"]), "black", "gray30"))
  
  ci.terms <- merge(ci.terms, models.df, all.x=T, all.y=F)

  # -------------------
  # Some formatting
  # -------------------
  {
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
  # -------------------
  
  
  png(file.path("Figures/analyses/analysis_TempExtent", "SuppFig1_Sensitivity_byModel_tair.png"), height=8, width=8, units="in", res=180)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="tair",]) +
    facet_wrap(~Model.Order, scales="free_x") +
    geom_ribbon(aes(x=x-273.15, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
    geom_line(aes(x=x-273.15, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
    scale_x_continuous(name=expression(bold(paste("Temperature, May-Sep ("^"o","C)"))), expand=c(0,0), breaks=c(10, 12.5, 15, 17.5)) +
    scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    theme(legend.title=element_text(size=10, face="bold"),
          legend.text=element_text(size=9),
          legend.key=element_blank(),
          legend.key.size=unit(1.5, "lines"),
          legend.background=element_blank(), 
          legend.position=c(0.83, 0.2)) +
    theme(strip.text.x=element_text(size=12, face="bold"),
          strip.text.y=element_text(size=12, face="bold")) + 
    theme(axis.line=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(fill=NA, color="black", size=0.5), 
          panel.background=element_blank(),
          panel.margin.x=unit(0, "lines"),
          panel.margin.y=unit(0.25, "lines"))  +
    theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
          axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
          axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
          axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
          axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
  dev.off()
  
  png(file.path("Figures/analyses/analysis_TempExtent", "SuppFig2_Sensitivity_byModel_precipf.png"), height=8, width=8, units="in", res=180)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="precipf",]) +
    facet_wrap(~Model.Order, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
    geom_line(aes(x=x, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
    scale_x_continuous(name=expression(bold(paste("Precipitation, May-Sep (mm yr"^"-1",")"))), expand=c(0,0)) +
    scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    theme(legend.title=element_text(size=10, face="bold"),
          legend.text=element_text(size=9),
          legend.key=element_blank(),
          legend.key.size=unit(1.5, "lines"),
          legend.background=element_blank(), 
          legend.position=c(0.92, 0.09)) +
    theme(strip.text.x=element_text(size=12, face="bold"),
          strip.text.y=element_text(size=12, face="bold")) + 
    theme(axis.line=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(fill=NA, color="black", size=0.5), 
          panel.background=element_blank(),
          panel.margin.x=unit(0, "lines"),
          panel.margin.y=unit(0.25, "lines"))  +
    theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
          axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
          axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
          axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
          axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
  dev.off()
  
  png(file.path("Figures/analyses/analysis_TempExtent", "SuppFig3_Sensitivity_byModel_CO2.png"), height=8, width=8, units="in", res=180)
  {
    ggplot(data=ci.terms[ci.terms$Effect=="CO2",]) +
    facet_wrap(~Model.Order, scales="free_x") +
    geom_ribbon(aes(x=x, ymin=lwr.rel.cent*100, ymax=upr.rel.cent*100, fill=Extent),alpha=0.5) +
    geom_line(aes(x=x, y=mean.rel.cent*100, color=Extent, linetype=Extent), size=2) +
    #       geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_x_continuous(name=expression(bold(paste("CO"["2"], " (ppm)"))), expand=c(0,0)) +
    scale_y_continuous(name="NPP Effect (%)", expand=c(0,0)) +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    scale_color_manual(values=c("blue3", "red3", "green3")) +
    theme(legend.title=element_text(size=12, face="bold"),
          legend.text=element_text(size=10),
          legend.key=element_blank(),
          legend.key.size=unit(1.5, "lines"),
          legend.background=element_blank(), 
          legend.position=c(0.85, 0.1)) +
    theme(strip.text.x=element_text(size=12, face="bold"),
          strip.text.y=element_text(size=12, face="bold")) + 
    theme(axis.line=element_line(color="black", size=0.5), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(fill=NA, color="black", size=0.5), 
          panel.background=element_blank(),
          panel.margin.x=unit(0, "lines"),
          panel.margin.y=unit(0.25, "lines"))  +
    theme(axis.text.y=element_text(color="black", size=10, margin=unit(c(0,1.5,0,0), "lines")),
          axis.text.x=element_text(color="black", size=10, margin=unit(c(1.5,0,0,0), "lines")), 
          axis.title.y=element_text(size=12, face="bold", margin=unit(c(0,0.5,0,0), "lines")),  
          axis.title.x=element_text(size=12, face="bold", margin=unit(c(0.5,0,0,0), "lines")),
          axis.ticks.length=unit(-0.5, "lines")) +
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  }
  dev.off()
}
# ----------------------------------------
