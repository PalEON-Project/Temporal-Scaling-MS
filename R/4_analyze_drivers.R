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
out.dir <- "Data/analysis_drivers"
fig.dir <- "Figures/analysis_drivers"

if(!dir.exists(out.dir)) dir.create(out.dir)
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
summary(drivers.stack)

n <- length(unique(drivers$Year))*length(unique(drivers$Site))
# Make some simple plots of the different distributions
pdf(file.path(fig.dir, "Driver_Distributions_by_Scale.pdf"))
print(
ggplot(data=drivers.stack) + facet_grid(Scale ~ driver, scales="free") +
	geom_histogram(aes(x=value, fill=Site, weight=1/n)) + 
	labs(y="p.obs",title="Distribution of PalEON Met Drivers with Temporal Resolution") +
	theme_bw()
)
dev.off()
# ----------------------------------------

# ----------------------------------------
# Look at driver correlations
# ----------------------------------------
# Make some simple plots of correlation
# pdf(file.path(fig.dir, "Driver_Correlations_by_Scale.pdf"))
# for(s in unique(drivers$Scale)){
	# print(
	# plot(drivers[drivers$Scale==s, vars.met], main=paste0("Met Driver Correlations, ", s))
	# )
# }
# dev.off()

for(s in unique(drivers$Scale)){
	png(file.path(fig.dir, paste0("Driver_Correlations_",s,".png")))
	plot(drivers[drivers$Scale==s, vars.met], main=paste0("Met Driver Correlations, ", s))
	dev.off()
}

# run the correlation matrices
scales <- unique(drivers$Scale)
vars.met

test <- merge(vars.met, vars.met)
data.temp <- drivers[driver$scale==scales[s],]
cor.scales <- data.frame(pair   = paste(vars.met[1], vars.met[2], sep=":"), 
						 scale  = scales[s],
						 cor.est= cor(data.temp[,vars.met[1:2]], use="pairwise.complete.obs")[2,2],
						 cor.lo = cor.test(data.temp[,vars.met[1]], data.temp[,vars.met[2]], use="pairwise.complete.obs")$conf.int[1],
						 cor.hi = cor.test(data.temp[,vars.met[1]], data.temp[,vars.met[2]], use="pairwise.complete.obs")$conf.int[2],
						 )
for(s in 1:length(scales)){
  data.temp <- drivers[drivers$Scale==scales[s],]
  cor.temp <- data.frame(pair   = paste(vars.met[1], vars.met[2], sep=":"), 
						 scale  = scales[s],
						 cor.est= cor(data.temp[,vars.met[1:2]], use="pairwise.complete.obs")[1,2],
						 cor.lo = cor.test(data.temp[,vars.met[1]], data.temp[,vars.met[2]], use="pairwise.complete.obs")$conf.int[1],
						 cor.hi = cor.test(data.temp[,vars.met[1]], data.temp[,vars.met[2]], use="pairwise.complete.obs")$conf.int[2]
						 )
	for(i in 2:length(vars.met)){
		# cor.temp2 <- data.frame(pair = paste(vars.met[i-1], vars.met[i], sep=":"), 
						   # scale  = scales[s],
						   # cor.est= cor(data.temp[,vars.met[(i-1):i]], use="pairwise.complete.obs")[2,2],
						   # cor.lo = cor.test(data.temp[,vars.met[(i-1)]], data.temp[,vars.met[i]], use="pairwise.complete.obs")$conf.int[1],
						   # cor.hi = cor.test(data.temp[,vars.met[(i-1)]], data.temp[,vars.met[i]], use="pairwise.complete.obs")$conf.int[2]
						   # )
		cor.temp <- rbind(cor.temp, data.frame(pair = paste(vars.met[i-1], vars.met[i], sep=":"), 
						   scale  = scales[s],
						   cor.est= cor(data.temp[,vars.met[(i-1):i]], use="pairwise.complete.obs")[1,2],
						   cor.lo = cor.test(data.temp[,vars.met[(i-1)]], data.temp[,vars.met[i]], use="pairwise.complete.obs")$conf.int[1],
						   cor.hi = cor.test(data.temp[,vars.met[(i-1)]], data.temp[,vars.met[i]], use="pairwise.complete.obs")$conf.int[2]
						   ))
	}
    # cor.temp <- cor(drivers[drivers$Scale==scales[s], vars.met], use="pairwise.complete.obs")
	# v1 <- rownames(cor.temp)
	# v2 <- stack(data.frame(cor.temp))[,2]
	# cor.temp2 <- data.frame(scale=scales[s], pair=paste(v1, v2, sep=":"), v1=v1, v2=v2, var.cor=stack(data.frame(cor.temp))[,1])
	if(s==1){
		cor.scales <- cor.temp
	} else {
		cor.scales <- rbind(cor.scales, cor.temp)
	}
}
summary(cor.scales)

# graph the changes in correlation
pdf(file.path(fig.dir, "Driver_Correlations_Across_Scales.pdf"))
ggplot(data=cor.scales[!cor.scales$v1==cor.scales$v2,]) + facet_grid(v1 ~ v2) +
	geom_line(aes(x=as.numeric(scale), y=var.cor), size=2) +
	geom_hline(data=cor.scales[!cor.scales$v1==cor.scales$v2 & cor.scales$scale=="t.001",], aes(yintercept=var.cor), color="red", size=0.5, linetype="dashed") +
	labs(x="Scale", y="Correlation (R2)", title="Change in Driver Correlation Across Scales") +
	theme_bw()
dev.off()

	cor2 <- cor.test(drivers[drivers$Scale=="t.001", vars.met[1]], drivers[drivers$Scale=="t.001", vars.met[2]], use="pairwise.complete.obs")
cor2$conf.int[2]
# ----------------------------------------
