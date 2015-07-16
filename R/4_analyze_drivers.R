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
# Note: because of the number of data points, these need to be saved in a "flattened" format
for(s in unique(drivers$Scale)){
	png(file.path(fig.dir, paste0("Driver_Correlations_",s,".png")))
	plot(drivers[drivers$Scale==s, vars.met], main=paste0("Met Driver Correlations, ", s))
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
				cor.scales <- cor.temp 
			} else {
				cor.scales <- rbind(cor.scales, cor.temp)
			}
		}
	}
}
# Making a sign color to make the graphs; note there's some tricks to getting everything to graph okay
cor.scales$sign <- as.factor(ifelse(cor.scales$cor.est>=0 | (cor.scales$cor.hi>=0 & !cor.scales$scale=="t.250"), "pos", "neg"))
summary(cor.scales)

# graph the changes in correlation
# Note: graphing the absolute value so that going down is always a decrease in correlation
pdf(file.path(fig.dir, "Driver_Correlations_Across_Scales.pdf"))
ggplot(data=cor.scales[,]) + facet_grid(v2 ~ v1) +
	# geom_ribbon(aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi)), fill="red2", alpha=0.7) +
	geom_hline(data=cor.scales[cor.scales$scale=="t.001",], aes(yintercept=abs(cor.est)), linetype="dashed", size=0.25) +
	geom_ribbon(data=cor.scales[cor.scales$sign=="neg",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	geom_ribbon(data=cor.scales[cor.scales$sign=="pos",], aes(x=as.numeric(scale), ymin=abs(cor.lo), ymax=abs(cor.hi), fill=sign), alpha=0.7) +
	# geom_line(aes(x=as.numeric(scale), y=cor.est), size=0.25) +
	scale_fill_manual(values=c("blue2", "red2")) +
	guides(fill=guide_legend(title="Corr. Sign")) +
	labs(x="Scale", y="Absolute Value of Correlation (R2)", title="Change in Driver Correlation Across Scales") +
	theme_bw()
dev.off()
# ----------------------------------------
