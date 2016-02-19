

# -----------------------
# 5.b. Analysis: 
# -----------------------
# --------
# 5.b.0. Getting the first derivative (first diference) of each line so we can take the mean slope
# --------
{
  # First make sure the effects are sorted by x to make this easier
  ci.terms <- ci.terms[order(ci.terms$Model,ci.terms$Extent, ci.terms$Effect, ci.terms$x),]
  sim.terms <- sim.terms[order(sim.terms$Model, sim.terms$PFT, sim.terms$Effect, sim.terms$x),]
  ci.terms[1:20,1:11]
  dim(ci.terms)
  summary(ci.terms)
  ci.terms$Y.type     <- as.factor(ifelse(ci.terms $Model=="TreeRingRW", "RW", "NPP"))
  sim.terms$Y.type    <- as.factor(ifelse(sim.terms$Model=="TreeRingRW", "RW", "NPP"))
  ci.terms$data.type  <- as.factor(ifelse(substr(ci.terms $Model,1,8)=="TreeRing", "Tree Rings", "Model"))
  sim.terms$data.type <- as.factor(ifelse(substr(sim.terms$Model,1,8)=="TreeRing", "Tree Rings", "Model"))
  
  # Looking at the range of the 95% CI
  ci.terms$rel.ci.range  <- ci.terms$upr.rel - ci.terms$lwr.rel
  
  # Making a new dataframe dedicated to the derivatives
  cols.sims <- which(substr(names(sim.terms),1,1)=="X")
  sim.deriv <- sim.terms[,]
  sim.deriv[,cols.sims] <- NA
  
  for(e in unique(ci.terms$Effect)){
    for(m in unique(ci.terms$Model)){
      for(p in unique(ci.terms[ci.terms$Model==m & ci.terms$Effect==e, "Extent"])){
        x.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "x"], lag=1), NA)
        y.dif <- c(diff(ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "mean.rel"], lag=1), NA)
        ci.terms[ci.terms$Model==m & ci.terms$Effect==e & ci.terms$Extent==p, "deriv"] <- y.dif/x.dif
        
        # For the full simiulation for robust analysis
        y.dif2 <- rbind(apply(sim.terms[sim.terms$Model==m & sim.terms$Effect==e & sim.terms$Extent==p, cols.sims], 2, FUN=diff), NA)
        x.dif <- c(diff(sim.terms[sim.terms$Model==m & sim.terms$Effect==e  & sim.terms$Extent==p, "x"], lag=1), NA)
        sim.deriv[sim.deriv$Model==m & sim.deriv$Effect==e & sim.terms$Extent==p, cols.sims] <- apply(y.dif2, 2, FUN=function(y){y/x.dif})
      }
    } 
  }
  summary(ci.terms)
  
  # Stacking and aggregating the simulations
  deriv.stack <- stack(sim.deriv[,cols.sims])
  names(deriv.stack)  <- c("deriv", "sim")
  deriv.stack[,names(sim.deriv)[which(!(1:ncol(sim.deriv)) %in% cols.sims)]] <- sim.deriv[,which(!(1:ncol(sim.deriv)) %in% cols.sims)]
  deriv.stack$Model  <- as.factor(deriv.stack$Model)
  deriv.stack$Effect <- as.factor(deriv.stack$Effect)
  summary(deriv.stack)
  
  sim.stack <- stack(sim.terms[,cols.sims])
  names(sim.stack)  <- c("Y", "sim")
  sim.stack[,names(sim.terms)[which(!(1:ncol(sim.terms)) %in% cols.sims)]] <- sim.terms[,which(!(1:ncol(sim.terms)) %in% cols.sims)]
  sim.stack$Model  <- as.factor(sim.stack$Model)
  sim.stack$Effect <- as.factor(sim.stack$Effect)
  summary(sim.stack)
  
  # Getting the mean derivative for each model & effet & Extent
  # -- Need to do this in two bins: overlap climate, extremes
  for(e in unique(ci.terms$Effect)){
    range.modern <- range(dat.ecosys[dat.ecosys$Extent=="1901-2010", e], na.rm=T)
    ci.terms[ci.terms$Effect==e, "Climate"] <- as.factor(ifelse(ci.terms[ci.terms$Effect==e, "x"]>=range.modern[1] & 
                                                                  ci.terms[ci.terms$Effect==e, "x"]<=range.modern[2],
                                                                "Overlap", "Extended"))
    sim.stack[sim.stack$Effect==e, "Climate"] <- as.factor(ifelse(sim.stack[sim.stack$Effect==e, "x"]>=range.modern[1] & 
                                                                    sim.stack[sim.stack$Effect==e, "x"]<=range.modern[2],
                                                                  "Overlap", "Extended"))
    deriv.stack[deriv.stack$Effect==e, "Climate"] <- as.factor(ifelse(deriv.stack[deriv.stack$Effect==e, "x"]>=range.modern[1] & 
                                                                        deriv.stack[deriv.stack$Effect==e, "x"]<=range.modern[2],
                                                                      "Overlap", "Extended"))
  }
  summary(ci.terms)
  summary(sim.stack)
  summary(deriv.stack)
  
  vars.deriv <- c("deriv", "x")
  deriv.agg                             <- aggregate(deriv.stack[,vars.deriv], 
                                                     by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                     FUN=mean, na.rm=T)
  deriv.agg[,paste0(vars.deriv, ".sd")] <- aggregate(deriv.stack[,vars.deriv], 
                                                     by=deriv.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                     FUN=sd, na.rm=T)[,vars.deriv]
  summary(deriv.agg)
  
  vars.sim <- c("Y", "x")
  sim.agg                            <- aggregate(sim.stack[,vars.sim], 
                                                  by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                  FUN=mean, na.rm=T)
  sim.agg[,paste0(vars.sim, ".sd")]  <- aggregate(sim.stack[,vars.sim], 
                                                  by=sim.stack[,c("Model", "Extent", "Y.type", "data.type", "Extent", "Effect", "Climate", "sim")], 
                                                  FUN=sd, na.rm=T)[,vars.sim]
  summary(sim.agg)
  
  
  vars.agg <- c("x", "mean.rel", "lwr.rel", "upr.rel", "rel.ci.range", "deriv")
  ci.terms.agg <- aggregate(ci.terms[,vars.agg],
                            by=ci.terms[,c("Model", "Y.type", "data.type", "Extent", "Effect", "Climate")],
                            FUN=mean, na.rm=T)
  
  # Ignoring 1980-2010 for analyses & presentation
  ci.terms.agg <- ci.terms.agg[!ci.terms.agg$Extent=="1980-2010",]
  sim.agg      <- sim.agg[!sim.agg$Extent=="1980-2010",]
  deriv.agg    <- deriv.agg[!deriv.agg$Extent=="1980-2010",]
  summary(ci.terms.agg)
  summary(sim.agg)
  summary(deriv.agg)
  
}
# --------

# --------
# 5.b.1. Looking at changes in sensitivity based on extent
# --------
library(nlme)
{
  summary(ci.terms.agg)
  ci.terms <- ci.terms[!ci.terms$Extent=="1980-2010",]
  sim.stack <- sim.stack[!sim.stack$Extent=="1980-2010",]
  
  # Tair
  tair.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
  summary(tair.overlap); anova(tair.overlap)
  tair.overlap2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
  summary(tair.overlap2); anova(tair.overlap2)
  
  tair.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
  summary(tair.overlap3); anova(tair.overlap3)
  tair.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Overlap",])
  summary(tair.overlap4); anova(tair.overlap4)
  
  
  tair.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Extended",])
  summary(tair.extended); anova(tair.extended)
  tair.extended2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="tair" & ci.terms.agg$Climate=="Extended",])
  summary(tair.extended2); anova(tair.extended2)
  
  # Precipf
  precipf.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
  summary(precipf.overlap); anova(precipf.overlap)
  precipf.overlap2 <- lm(rel.ci.range ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
  summary(precipf.overlap2); anova(precipf.overlap2)
  
  precipf.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
  summary(precipf.overlap3); anova(precipf.overlap3)
  precipf.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Overlap",])
  summary(precipf.overlap4); anova(precipf.overlap4)
  
  precipf.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="precipf" & ci.terms.agg$Climate=="Extended",])
  summary(precipf.extended); anova(precipf.extended)
  
  # CO2
  CO2.overlap <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
  summary(CO2.overlap); anova(CO2.overlap)
  
  CO2.overlap3 <- lme(deriv ~ Extent, random=list(Model=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
  summary(CO2.overlap3); anova(CO2.overlap3)
  CO2.overlap4 <- lme(deriv ~ Model, random=list(Extent=~1), data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Overlap",])
  summary(CO2.overlap4); anova(CO2.overlap4)
  
  CO2.extended <- lm(deriv ~ Extent, data=ci.terms.agg[ci.terms.agg$data.type=="Model" & ci.terms.agg$Effect=="CO2" & ci.terms.agg$Climate=="Extended",])
  summary(CO2.extended); anova(CO2.extended)
}
# --------

# --------
# 5.b.2. Looking explicitly at change in sensitivity among models
# --------
{
  summary(deriv.stack)
  deriv.stack2 <- deriv.stack[deriv.stack$Extent=="850-2010",]
  summary(deriv.stack2)
  
  for(m in unique(deriv.stack2$Model)){
    deriv.stack2[deriv.stack2$Model==m, "dDeriv"] <- deriv.stack[deriv.stack$Model==m & deriv.stack$Extent=="850-2010","deriv"] - deriv.stack[deriv.stack$Model==m & deriv.stack$Extent=="1901-2010","deriv"]
  }
  summary(deriv.stack2)
  
  # Aggregating a bit; ignoring the periods of climate overlap
  deriv.agg2 <- aggregate(deriv.stack2[,c("deriv", "dDeriv", "x")], 
                          by=deriv.stack2[,c("Model", "Y.type", "data.type", "Effect", "sim")], 
                          FUN=mean, na.rm=T)
  
  mean(abs(deriv.agg2[deriv.agg2$Effect=="CO2","dDeriv"]), na.rm=T); sd(abs(deriv.agg2[deriv.agg2$Effect=="CO2","dDeriv"]), na.rm=T)
  mean(deriv.agg2[deriv.agg2$Effect=="tair","dDeriv"], na.rm=T); sd(deriv.agg2[deriv.agg2$Effect=="tair","dDeriv"], na.rm=T)
  mean(deriv.agg2[deriv.agg2$Effect=="precipf","dDeriv"], na.rm=T); sd(deriv.agg2[deriv.agg2$Effect=="precipf","dDeriv"], na.rm=T)
  
  
  
  deriv.agg3 <- aggregate(deriv.agg2[,c("deriv", "dDeriv", "x")], 
                          by=deriv.agg2[,c("Model", "Y.type", "data.type", "Effect")], 
                          FUN=mean, na.rm=T)
  deriv.agg3
  mean(deriv.agg3[deriv.agg3$Effect=="CO2","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="CO2","dDeriv"], na.rm=T)
  mean(deriv.agg3[deriv.agg3$Effect=="tair","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="tair","dDeriv"], na.rm=T)
  mean(deriv.agg3[deriv.agg3$Effect=="precipf","dDeriv"], na.rm=T); sd(deriv.agg3[deriv.agg3$Effect=="precipf","dDeriv"], na.rm=T)
  
  d.co2.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="CO2",])
  d.tair.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="tair",])
  d.precipf.mod <- lm(abs(dDeriv) ~ Model-1, data=deriv.agg2[deriv.agg2$Effect=="precipf",])
  summary(d.co2.mod)
  summary(d.tair.mod)
  summary(d.precipf.mod)
  
  d.co2.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="CO2",])
  d.tair.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="tair",])
  d.precipf.mod2 <- lm(abs(dDeriv) ~ 1, data=deriv.agg2[deriv.agg2$Effect=="precipf",])
  summary(d.co2.mod2)
  summary(d.tair.mod2)
  summary(d.precipf.mod2)
}
# --------

# --------
# 5.b.3. Attributing change in sensitivity to change in biomass and composition
# --------
{
  summary(ecosys)
  vars.agg <- c("NPP", "AGB", "LAI", "Evergreen", "Deciduous", "Grass")
  # Calculating precent change in time in the two reference periods
  mod.agg1 <- aggregate(ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=850,vars.agg], 
                        by=ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=850,c("Model", "Model.Order")],
                        FUN=mean, na.rm=T)
  mod.agg2 <- aggregate(ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=1901,vars.agg], 
                        by=ecosys[ecosys$Resolution=="t.001" & ecosys$Year>=1901,c("Model", "Model.Order")],
                        FUN=mean, na.rm=T)
  mod.agg1$Extent <- "850-2010"
  mod.agg2$Extent <- "1901-2010"
  
  mod.agg <- mod.agg1[,1:2]
  mod.agg[,vars.agg] <- (mod.agg1[,vars.agg]-mod.agg2[,vars.agg])/mod.agg1[,vars.agg]
  
  mod.agg[mod.agg$Model=="jules.stat", "AGB"] <- mod.agg[mod.agg$Model=="jules.stat", "LAI"] 
  mod.agg[is.na(mod.agg)] <- 0
  mod.agg
  
  # Merging the percent change in time with 
  deriv.agg3 <- merge(deriv.agg3, mod.agg, all.x=T, all.y=T)
  deriv.agg3
  
  effect.AGB <- lm(abs(dDeriv) ~ AGB*(Effect-1)-Effect-AGB, data=deriv.agg3[,])
  summary(effect.AGB)
  
  effect.evg <- lm(abs(dDeriv) ~ Evergreen*(Effect-1)-Effect-Evergreen, data=deriv.agg3[,])
  summary(effect.evg)
}
# --------


}
# -----------------------
# ----------------------------------------

















# ----------------------------------------
# 6. Graphing & Analyzing Emulated Change in NPP
# ----------------------------------------
# -----------------------
# 6.a. Graphing
# -----------------------
wt.terms2 <- merge(wt.terms, dat.ecosys2[dat.ecosys2$Extent=="850-2010",c("Model", "Model.Order", "Y.type", "Site", "data.type", "Year", "Y", "Y.rel", "Y.10", "Y.rel.10")], all.x=T, all.y=F)
{
  # summary(wt.terms2)
  dim(wt.terms); dim(wt.terms2)
  # --------------------------
  # Adjusting CO2 Effect
  # --------------------------
  # Note: because the gam makes the smoother cross 0 at the MEAN CO2 (which is in the 1800s), 
  # it's saying the region is pretty CO2-limited at periods where that doesn't really make 
  # sense, so we're going to off relativize it to whatever the starting point for the run is
  # --------------------------
{
  for(m in unique(wt.terms2$Model)){
    for(e in unique(wt.terms2[wt.terms2$Model==m, "Extent"])){
      yr1       <- min(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "Year"]) # find the minimum year
      yr2       <- min(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & !is.na(wt.terms2$weight.CO2), "Year"]) # find the minimum year
      co2.base <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & wt.terms2$Year<=(yr1+5), "fit.CO2"], na.rm=T) # mean of the 5 years around the starting 
      co2.base.10 <- mean(wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e & wt.terms2$Year<=(yr2+5),"fit.CO2.10"],na.rm=T) # mean of the 5 years around the starting point
      wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.adj"] <- wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2"] - co2.base
      wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.10.adj"] <- wt.terms2[wt.terms2$Model==m & wt.terms2$Extent==e, "fit.CO2.10"] - co2.base.10
    }
  }
  summary(wt.terms2)
  
  wt.terms2[,c("weight.tair.adj", "weight.precipf.adj","weight.CO2.adj")] <- abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")])/rowSums(abs(wt.terms2[,c("fit.tair", "fit.precipf", "fit.CO2.adj")]))
  wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")] <- abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")])/rowSums(abs(wt.terms2[,c("fit.tair.10", "fit.precipf.10", "fit.CO2.10.adj")]))
  
  wt.terms2[is.na(wt.terms2$weight.tair.10), c("weight.tair.10", "weight.precipf.10", "weight.CO2.10")] <- 0
  wt.terms2[is.na(wt.terms2$weight.tair.10.adj), c("weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")] <- 0
  
  summary(rowSums(wt.terms2[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")]))
  
  summary(wt.terms2)
}
# -----------
# Aggregating to the region ensemble-level
# -----------
{
  # factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
  factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "Y.rel", "Y.rel.10", "weight.tair.adj", "weight.precipf.adj", "weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")
  
  # aggregate across plots
  ensemble.wts0 <- aggregate(wt.terms2[,factors.aggregate], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=mean, na.rm=T)
  ensemble.wts0[,paste0(factors.aggregate, ".lo")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=quantile, 0.025, na.rm=T)[,factors.aggregate[1:5]]
  ensemble.wts0[,paste0(factors.aggregate, ".hi")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], by=wt.terms2[,c("Model", "data.type", "Year", "Extent")], FUN=quantile, 0.975, na.rm=T)[,factors.aggregate[1:5]]
  summary(ensemble.wts0)
  
  ensemble.wts1 <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=mean, na.rm=T)
  summary(ensemble.wts1)
  
  ensemble.wts.lo <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=quantile, 0.025, na.rm=T)
  names(ensemble.wts.lo)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.lo[5:ncol(ensemble.wts.lo)]), ".lo")) 
  summary(ensemble.wts.lo)
  
  ensemble.wts.hi <- aggregate(ensemble.wts0[,factors.aggregate], by= ensemble.wts0[,c("data.type", "Year", "Extent")], FUN=quantile, 0.975, na.rm=T)
  names(ensemble.wts.hi)[5:ncol(ensemble.wts.lo)] <- c(paste0(names(ensemble.wts.hi[5:ncol(ensemble.wts.hi)]), ".hi")) 
  summary(ensemble.wts.hi)
  
  ensemble.wts.final <- cbind(ensemble.wts1, ensemble.wts.lo[5:ncol(ensemble.wts.lo)], ensemble.wts.hi[5:ncol(ensemble.wts.hi)])
  summary(ensemble.wts.final)
}
# -----------

# -----------
# Aggregating to the site ensemble-level
# -----------
{
  # factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "fit.tair", "fit.tair.rel", "weight.tair", "fit.precipf", "fit.precipf.rel", "weight.precipf", "fit.CO2", "fit.CO2.rel", "weight.CO2", "Y.rel", "Y.10", "Y.rel.10", "weight.tair.10", "weight.precipf.10", "weight.CO2.10")
  factors.aggregate <- c("fit.full", "fit.full.rel", "fit.full.rel.10", "Y.rel", "Y.rel.10", "weight.tair.adj", "weight.precipf.adj", "weight.CO2.adj", "weight.tair.10.adj", "weight.precipf.10.adj", "weight.CO2.10.adj")
  
  # aggregate across plots
  ensemble.wts.site                                    <- aggregate(wt.terms2[,factors.aggregate], 
                                                                    by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                    FUN=mean, na.rm=T)
  ensemble.wts.site[,paste0(factors.aggregate, ".lo")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], 
                                                                    by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                    FUN=quantile, 0.025, na.rm=T)[,factors.aggregate[1:5]]
  ensemble.wts.site[,paste0(factors.aggregate, ".hi")] <- aggregate(wt.terms2[,factors.aggregate[1:5]], 
                                                                    by=wt.terms2[,c("Site", "data.type", "Year", "Extent")], 
                                                                    FUN=quantile, 0.975, na.rm=T)[,factors.aggregate[1:5]]
  summary(ensemble.wts.site)
  
}
# -----------

# -----------
# Re-normalizing factor weights
# -----------
{
  summary(rowSums(ensemble.wts.final[,c("weight.tair.10.adj", "weight.precipf.10.adj","weight.CO2.10.adj")]))
  
  
  # {
  # wts.sum.10 <- abs(ensemble.wts.final$weight.tair.10) + abs(ensemble.wts.final$weight.precipf.10) + abs(ensemble.wts.final$weight.CO2.10)
  # ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")] <- ensemble.wts.final[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")]/wts.sum.10
  # ensemble.wts.final[is.na(ensemble.wts.final$weight.tair.10   ),"weight.tair.10"   ] <- 0
  # ensemble.wts.final[is.na(ensemble.wts.final$weight.precipf.10),"weight.precipf.10"] <- 0
  # ensemble.wts.final[is.na(ensemble.wts.final$weight.CO2.10    ),"weight.CO2.10"    ] <- 0
  
  
  # wts.sum <- abs(ensemble.wts.final$weight.tair) + abs(ensemble.wts.final$weight.precipf) + abs(ensemble.wts.final$weight.CO2)
  # ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")] <- ensemble.wts.final[,c("weight.tair","weight.precipf", "weight.CO2")]/wts.sum
  # summary(ensemble.wts.final)
  # }
  
  summary(ensemble.wts.final)
  
  ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
  ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
  ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
  # ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$ci.min)
  # ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
}
# -----------

# -----------
# Region figures
# -----------
{
  pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_1500-2010_Decadal_AllExtent.pdf"), width=11, height=8.5)
{
    print(
      ggplot(ensemble.wts.final) + facet_grid(Extent~., scales="fixed") +
        scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
        geom_ribbon(data= ensemble.wts.final[,], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
        geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, y=fit.graph*100),
                  color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1980-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="1980-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="1980-2010","weight.precipf.10.adj"])), size=3) +
        geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                  color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
        geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                  color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        geom_hline(yintercept=100, linetype="dashed") +
        scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
        # ggtitle("NPP Controlling Factor") + 
        theme(legend.text=element_text(size=rel(1)), 
              legend.title=element_text(size=rel(1)),
              legend.key=element_blank(),
              legend.key.size=unit(1, "lines")) + 
        # legend.key.width=unit(2, "lines")) + 
        theme(axis.line=element_line(color="black", size=0.5), 
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(), 
              panel.background=element_blank(),
              panel.margin.y=unit(0.5, "lines"))  +
        theme(axis.text.x=element_text(size=rel(1), color="black"),
              axis.text.y=element_text(size=rel(1), color="black"), 
              axis.title.x=element_text(size=rel(1), face="bold"),  
              axis.title.y=element_text(size=rel(1), face="bold"),
              # axis.ticks.length=unit(-0.5, "lines"),
              axis.ticks.margin=unit(1.0, "lines"))
    )
  }
dev.off()

pdf(file.path(fig.dir, "Fig5_Ensemble_Drivers_Time_Region_1500-2010_Decadal_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0.5,0.5,ensemble.wts.final$fit.full.rel.10.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      # ggtitle("NPP Controlling Factor") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
dev.off()

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_0850-2010_Decadal_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.10.hi>2,2,ensemble.wts.final$fit.full.rel.10.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.10.lo<0,0,ensemble.wts.final$fit.full.rel.10.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel.10<0,NA,ensemble.wts.final$fit.full.rel.10)
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(850,2010), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      # ggtitle("NPP Controlling Factor") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(900,1300), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      ggtitle("Medieval Warm Period") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1400,1800), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.10.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      ggtitle("Little Ice Age") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
dev.off()

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_Annual_AllExtent.pdf"), width=11, height=8.5)
ensemble.wts.final$ci.max <- ifelse(ensemble.wts.final$fit.full.rel.hi>2,2,ensemble.wts.final$fit.full.rel.hi)
ensemble.wts.final$ci.min <- ifelse(ensemble.wts.final$fit.full.rel.lo<0,0,ensemble.wts.final$fit.full.rel.lo)
ensemble.wts.final$fit.graph <- ifelse(ensemble.wts.final$fit.full.rel<0,NA,ensemble.wts.final$fit.full.rel)
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1850,2010), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      # ggtitle("NPP Controlling Factor") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
{
  print(
    ggplot(ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
      scale_x_continuous(limits=c(1750,1850), expand=c(0,0)) +
      geom_ribbon(data= ensemble.wts.final[!ensemble.wts.final$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="1901-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="1901-2010","weight.precipf.adj"])), size=3) +
      geom_line(data= ensemble.wts.final[ensemble.wts.final$Extent=="850-2010",], aes(x=Year, y=fit.full.rel*100),
                color=rgb(abs(ensemble.wts.final[ensemble.wts.final$Extent=="850-2010","weight.tair.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.CO2.adj"]),
                          abs(ensemble.wts.final[ensemble.wts.final$Extent =="850-2010","weight.precipf.adj"])), size=3) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
      # ggtitle("NPP Controlling Factor") + 
      theme(legend.text=element_text(size=rel(1)), 
            legend.title=element_text(size=rel(1)),
            legend.key=element_blank(),
            legend.key.size=unit(1, "lines")) + 
      # legend.key.width=unit(2, "lines")) + 
      theme(axis.line=element_line(color="black", size=0.5), 
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(), 
            panel.background=element_blank(),
            panel.margin.y=unit(0.5, "lines"))  +
      theme(axis.text.x=element_text(size=rel(1), color="black"),
            axis.text.y=element_text(size=rel(1), color="black"), 
            axis.title.x=element_text(size=rel(1), face="bold"),  
            axis.title.y=element_text(size=rel(1), face="bold"),
            # axis.ticks.length=unit(-0.5, "lines"),
            axis.ticks.margin=unit(1.0, "lines"))
  )
}
dev.off()

}
# -----------

# -----------
# Site Figures
# -----------
{
  pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Sites_1500-2010_Decadal.pdf"), width=11, height=8.5)
  ensemble.wts.site$ci.max    <- ifelse(ensemble.wts.site$fit.full.rel.10.hi > 2  , 2  , ensemble.wts.site $fit.full.rel.10.hi)
  ensemble.wts.site$ci.min    <- ifelse(ensemble.wts.site$fit.full.rel.10.lo < 0.25, 0.25, ensemble.wts.site $fit.full.rel.10.lo)
  ensemble.wts.site$fit.graph <- ifelse(ensemble.wts.site$fit.full.rel.10   < 0.25   , NA , ensemble.wts.site $fit.full.rel.10)
  for(s in unique(ensemble.wts.site$Site)){
    print(
      ggplot(ensemble.wts.site[ensemble.wts.site$Site==s & !ensemble.wts.site$Extent=="1980-2010",]) + facet_grid(Extent~., scales="fixed") +
        scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
        geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Site==s & !ensemble.wts.site$Extent=="1980-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="1901-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="1901-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="1901-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site==s & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        geom_hline(yintercept=100, linetype="dashed") +
        scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
        ggtitle(s) + 
        theme(legend.text=element_text(size=rel(1)), 
              legend.title=element_text(size=rel(1)),
              legend.key=element_blank(),
              legend.key.size=unit(1, "lines")) + 
        # legend.key.width=unit(2, "lines")) + 
        theme(axis.line=element_line(color="black", size=0.5), 
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(), 
              panel.background=element_blank(),
              panel.margin.y=unit(0.5, "lines"))  +
        theme(axis.text.x=element_text(size=rel(1), color="black"),
              axis.text.y=element_text(size=rel(1), color="black"), 
              axis.title.x=element_text(size=rel(1), face="bold"),  
              axis.title.y=element_text(size=rel(1), face="bold"),
              # axis.ticks.length=unit(-0.5, "lines"),
              axis.ticks.margin=unit(1.0, "lines"))
    )
  }
  dev.off()
  
  
  pdf(file.path(fig.dir, "Ensemble_Drivers_Time_AllSites_1500-2010_Decadal.pdf"), width=11, height=8.5)
{
    print(
      ggplot(ensemble.wts.site[ensemble.wts.site$Extent=="850-2010",]) + 
        facet_grid(Site~., scales="fixed") +
        scale_x_continuous(limits=c(1500,2010), expand=c(0,0)) +
        geom_ribbon(data= ensemble.wts.site[ensemble.wts.site$Extent=="850-2010",], aes(x=Year, ymin=ci.min*100, ymax=ci.max*100), alpha=0.35) +
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        
        geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent=="850-2010",], aes(x=Year, y= fit.graph*100),
                  color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent=="850-2010","weight.tair.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent =="850-2010","weight.CO2.10.adj"]),
                            abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
        geom_hline(yintercept=100, linetype="dashed") +
        scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
        theme(legend.text=element_text(size=rel(1)), 
              legend.title=element_text(size=rel(1)),
              legend.key=element_blank(),
              legend.key.size=unit(1, "lines")) + 
        # legend.key.width=unit(2, "lines")) + 
        theme(axis.line=element_line(color="black", size=0.5), 
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(), 
              panel.background=element_blank(),
              panel.margin.y=unit(0.5, "lines"))  +
        theme(axis.text.x=element_text(size=rel(1), color="black"),
              axis.text.y=element_text(size=rel(1), color="black"), 
              axis.title.x=element_text(size=rel(1), face="bold"),  
              axis.title.y=element_text(size=rel(1), face="bold"),
              # axis.ticks.length=unit(-0.5, "lines"),
              axis.ticks.margin=unit(1.0, "lines"))
    )
  }
dev.off()
}
# -----------

}

# -----------------------

# -----------------------
# 6.b. Analysis: 
# -----------------------
{
  source("R/0_Calculate_GAMM_Posteriors.R")
  library(mgcv)
  # --------
  # 6.b.0. Adding some model-level stats to compare the relative sensitivities
  # --------
{
  # 1. Quantify with-in model sensitivity shifts due to change in temporal extent
{
  # i. classify met by extent as "observed" in that extent and "extrapolated" 
  summary(ci.terms)
  # ii.  
}

# 2. Condensing model variability across space and time to get general model characteristics
{
  vars.agg <- c("Y", "Y.rel", "Biomass")
  mod.agg                           <- aggregate(dat.ecosys[,vars.agg], 
                                                 by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=mean)
  mod.agg[,paste0(vars.agg, ".sd")] <- aggregate(dat.ecosys[,vars.agg], 
                                                 by=dat.ecosys[,c("Model", "Model.Order", "Y.type", "data.type")], 
                                                 FUN=sd)[,vars.agg]
  mod.agg
  
  
  # Finding the change in key variables in the modern era
  for(v in vars.agg){
    mod.agg[mod.agg$Model==m,paste0("dModern.", v)] <- NA  
  }
  
  for(m in unique(mod.agg$Model)){
    mod.agg[mod.agg$Model==m, paste0("dModern.", vars.agg)] <- colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1990 & dat.ecosys$Year<=2010, vars.agg], na.rm=T) -
      colMeans(dat.ecosys[dat.ecosys$Model==m & dat.ecosys$Year>=1830 & dat.ecosys$Year<=1850, vars.agg], na.rm=T)
  }
  mod.agg
  
  # Add in the vegetation scheme
  mod.agg$veg.scheme <- as.factor(ifelse(mod.agg$Model %in% c("clm.bgc", "clm.cn", "sibcasa", "jules.stat", "jules.triffid"), "Static", "Dynamic"))
  
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
}

ci.terms <- merge(ci.terms, mod.agg, all.x=T, all.y=F)
summary(ci.terms)

# 3. Making the data frame so we can graph & compare the curves quantitatively
{
  factors.analy <- c("Y", "Y.rel", "Y.sd", "Y.rel.sd", "dModern.Y", "dModern.Y.rel", "dModern.Biomass", "Biomass", "Biomass.sd")
  
  df.co2 <- aggregate(ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",c(factors.analy)], 
                      by=ci.terms[ci.terms$Effect=="CO2" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                      FUN=mean, na.rm=T)
  summary(df.co2)
  
  df.tair <- aggregate(ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",c(factors.analy)], 
                       by=ci.terms[ci.terms$Effect=="tair" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                       FUN=mean, na.rm=T)
  summary(df.tair)
  
  df.precipf <- aggregate(ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",c(factors.analy)], 
                          by=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",c("x", "veg.scheme", "fire.scheme", "Model")],
                          FUN=mean, na.rm=T)
  summary(df.precipf)
}

# 4. setting up some null models for the full sensitivity curves
{
  co2.null     <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="CO2"     & ci.terms$data.type=="Model",])
  tair.null    <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="tair"    & ci.terms$data.type=="Model",])
  precipf.null <- gam(mean.rel ~ s(x), data=ci.terms[ci.terms$Effect=="precipf" & ci.terms$data.type=="Model",])
  
  summary(co2.null)
  summary(tair.null)
  summary(precipf.null)
  
  # Plot the different sensitivity curves by characteristic
  co2.null.post <- post.distns(model.gam=co2.null, model.name="CO2", n=50, newdata=df.co2, vars="x", terms=F)$ci
  co2.null.post$null.scheme <- df.co2$null.scheme
  co2.null.post$Effect <- as.factor("CO2")
  summary(co2.null.post)
  
  tair.null.post <- post.distns(model.gam=tair.null, model.name="Tair", n=50, newdata=df.tair, vars="x", terms=F)$ci
  tair.null.post$null.scheme <- df.tair$null.scheme
  tair.null.post$Effect <- as.factor("tair")
  summary(tair.null.post)
  
  precipf.null.post <- post.distns(model.gam=precipf.null, model.name="Precipf", n=50, newdata=df.precipf, vars="x", terms=F)$ci
  precipf.null.post$null.scheme <- df.precipf$null.scheme
  precipf.null.post$Effect <- as.factor("precipf")
  summary(precipf.null.post)
  
  null.post <- rbind(tair.null.post, precipf.null.post, co2.null.post)
  summary(null.post)
}

# 5. Condensing the full sensitivity curves to the values at the 25, 50, and 75% for 
#    analysis of continuous characteristics of models (i.e. NPP, modern change, etc)
{
  summary(ci.terms)
  summary(dat.ecosys)
  
  co2.driver     <- dat.ecosys[dat.ecosys$Model=="ed2","CO2"    ]
  tair.driver    <- dat.ecosys[dat.ecosys$Model=="ed2","tair"   ]
  precipf.driver <- dat.ecosys[dat.ecosys$Model=="ed2","precipf"]
  
  for(e in c("tair", "precipf", "CO2")){
    # Get the distirbution of met drivers for each driver & specify to what precision we want to round
    effect.driver <- dat.ecosys[dat.ecosys$Model=="ed2", e]
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
}
# --------


# --------
# 6.b.1 Ensemble-level Cross-Scale shifts in sensitivity
#   Hypothesis: Model sensitivities are more similar at short temporal scales
#               because it's the differences in long-term feedbacks that cuase
#               them to diverge
# --------
{
  
  # Calculate ensemeble mean sensitivitiy & deviation of individual models from 
  #   the ensemble at each temporal scale
{  
}

# Calculate ensemble Biome sensitivity & model deviation from the ensemble   
{
  
}

# Statistical test!
{
  
}

}
# --------

# --------
# 6.b.2 Finding patterns in Model Shifts: dynamic or static veg
#  Hypothesis: Climate sensitivities of models with static veg are more consistent
#              across temporal scales because they less capacity for ecosystem-level
#              adaptation to climate variaibility & directional change
# --------
{
}
# --------

# --------
# 6.b.3 Slow processes & sensitivity change: Biomass
#  Hypothesis: Models with greater changes in biomass have greater differences
#              among scales because change in biomass is a slow stabilizing process 
#              that can't happen at short time scales
# --------
{
}
# --------

# --------
# 6.b.3 Slow processes & sensitivity change: Composition (Fraction Evergreen)
#  Hypothesis: Like Biomass, shifts in composition is a slow process that can mediate
#              climate sensitivity at more long time scales but has limited influence
#              at short time scales
# --------
{
  
}
# --------
}
# -----------------------
# ----------------------------------------


















summary(ensemble.wts0)
summary(ensemble.wts.final)

weight.CO2     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010",])
weight.tair    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010",])
weight.precipf <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010",])
summary(weight.CO2)
summary(weight.tair)
summary(weight.precipf)


weight.CO22     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year>=1850,])
weight.tair2    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year>=1850,])
weight.precipf2 <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year>=1850,])
summary(weight.CO22)
summary(weight.tair2)
summary(weight.precipf2)

weight.CO22     <- lme(abs(weight.CO2.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year<=1850,])
weight.tair2    <- lme(abs(weight.tair.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year<=1850,])
weight.precipf2 <- lme(abs(weight.precipf.adj) ~ Extent, random=list(Year=~1), data=ensemble.wts0[!ensemble.wts0$Extent=="1980-2010" & ensemble.wts0$Year<=1850,])
summary(weight.CO22)
summary(weight.tair2)
summary(weight.precipf2)

summary(ensemble.wts0)
mean(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])

mean(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1990,"fit.full.rel"])

mean(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
quantile(ensemble.wts0[ensemble.wts0$Extent=="1901-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel.10"], c(0.025, 0.975))


mean(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
sd(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel"])
quantile(ensemble.wts0[ensemble.wts0$Extent=="850-2010" & ensemble.wts0$Year>=1830 & ensemble.wts0$Year<=1850,"fit.full.rel.10"], c(0.025, 0.975))

vars.compare <- c("Model", "Year", "fit.full.rel", "weight.CO2.adj", "weight.tair.adj", "weight.precipf.adj")
ensemble.ext1 <- ensemble.wts0[ensemble.wts0$Extent=="850-2010",vars.compare]
names(ensemble.ext1)[3:length(vars.compare)] <- paste0(vars.compare[3:length(vars.compare)],".850")
ensemble.ext2 <- ensemble.wts0[ensemble.wts0$Extent=="1901-2010",vars.compare]
names(ensemble.ext2)[3:length(vars.compare)] <- paste0(vars.compare[3:length(vars.compare)],".1901")

ensemble.comparison <- merge(ensemble.ext1, ensemble.ext2, all.x=T, all.y=T)
summary(ensemble.comparison)

ggplot(data=ensemble.comparison) + 
  facet_wrap(~Model, scales="fixed") +
  geom_point(aes(x= fit.full.rel.850, y= fit.full.rel.1901)) +
  coord_cartesian() +
  theme_bw()


plot(abs(weight.CO2.adj) ~ Extent, data=ensemble.wts0[,])
# -----------------------
# ----------------------------------------
