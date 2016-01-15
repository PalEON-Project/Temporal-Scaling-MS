# Climate Means Mapping Locally

# Script to make some maps of climate averages over the paleon domain

# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(raster)
library(ncdf4); library(ggplot2); library(grid)

met.dir  <- "~/Dropbox/PalEON_CR/paleon_mip_site/Analyses/Temporal-Scaling/Data/climate_means"
veg.dir  <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/biome/biome_potential_vegtype_biome.nc"
fig.dir  <- "~/Dropbox/PalEON_CR/paleon_mip_site/Analyses/Temporal-Scaling/Figures/phase2_met_regional_climate_means"
if(!dir.exists(fig.dir)) dir.create(fig.dir)

met.slices <- data.frame(start=c(0851, 1851, 1981), end=c(0880, 1880, 2010))
sites <- data.frame(Site=c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB"), name=c("Harvard", "Howland", "UNDERC", "Billy's", "Demming", "Minden"), lon=c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83), lat=c(42.54, 45.25, 46.22, 46.28, 47.17, 43.68))

met.vars <- c("tair", "precipf")
paleon.states <- map_data("state")
# ----------------------------------------------

# ----------------------------------------------
# Mapping Climate
# ----------------------------------------------
tair   <- raster(file.path(met.dir, "PalEON_Met_Normals_tair_1981-2010"))
precip <- raster(file.path(met.dir, "PalEON_Met_Normals_precipf_1981-2010"))

tair.x   <- data.frame(rasterToPoints(tair))
precipf.x <- data.frame(rasterToPoints(precip))
names(tair.x)   <- c("lon", "lat", "tair")
names(precipf.x) <- c("lon", "lat", "precipf")

summary(tair.x)
summary(precipf.x)

# Convert to celcius & mm/yr
tair.x$tair       <- tair.x$tair-273.15
precipf.x$precipf <- precipf.x$precipf*1*60*60*24*365.25 # 1 sec * 60 sec/min * 60 min/hr * 24 hr/day * 365.25 days/yr

tiff(file.path(fig.dir, "Tair_1981-2010_withSites.tiff"), height=900, width=1600)
ggplot(data=tair.x) +
		geom_raster(aes(x=lon, y=lat, fill=tair)) +
		scale_fill_gradientn(name=expression(bold(paste("Mean Ann Temp ("^"o","C)"))), colours=c("blue4", "red3")) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=1) +
		geom_point(data=sites, aes(x=lon, y=lat), size=20, color="darkgoldenrod3") +
		scale_x_continuous(limits=range(tair.x$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(tair.x$lat), expand=c(0,0), name="Latitude") +
		coord_equal(ratio=1) +
		# ggtitle(paste("Tair", year, mon, sep=" - ")) +
		guides(fill=guide_colorbar(title.position="top")) +
        theme(plot.margin=unit(c(1,1,1,1), "lines")) +
        theme(legend.position=c(0.83, 0.2), 
              legend.direction="horizontal",
              legend.key.size=unit(4, "lines"),
              legend.text=element_text(color="black", size=rel(3.5)),
              legend.title=element_text(color="black", size=rel(3.5))) +
        theme(axis.line=element_blank(), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(4)), 
	      axis.text.y=element_text(color="black", size=rel(4)), 
	      axis.title.x=element_text(face="bold", size=rel(4), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(4), vjust=1)) 
dev.off()

# tiff(file.path(fig.dir, "Tair_1981-2010_PHA_PHO_PDL.tiff"), height=1800, width=3200)
ggplot(data=tair.x) +
		geom_raster(aes(x=lon, y=lat, fill=tair)) +
		scale_fill_gradientn(name=expression(bold(paste("Mean Ann Temp ("^"o","C)"))), colours=c("blue4", "red3")) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=1) +
		# geom_point(data=sites, aes(x=lon, y=lat), size=8, color="darkgoldenrod3") +
		geom_point(data=sites[sites$Site %in% c("PHA", "PDL", "PUN"),], aes(label=name, x=lon, y=lat), size=8, color="darkgoldenrod3") +
		geom_text(data=sites[sites$Site %in% c("PHA", "PDL", "PUN"),], aes(label=name, x=lon-c(0.5,0.5,0.25), y=lat-c(-1,1,1)), size=8, color="darkgoldenrod3", face="bold") +
		scale_x_continuous(limits=range(tair.x$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(tair.x$lat), expand=c(0,0), name="Latitude") +
		coord_equal(ratio=1) +
		# ggtitle(paste("Tair", year, mon, sep=" - ")) +
		guides(fill=guide_colorbar(title.position="top")) +
        theme(plot.margin=unit(c(1,1,1,1), "lines")) +
        theme(legend.position=c(0.8, 0.2), 
              legend.direction="horizontal",
              legend.key.size=unit(2, "lines"),
              legend.text=element_text(color="black", size=rel(1.5)),
              legend.title=element_text(color="black", size=rel(1.5))) +
        theme(axis.line=element_blank(), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(1.5)), 
	      axis.text.y=element_text(color="black", size=rel(1.5)), 
	      axis.title.x=element_text(face="bold", size=rel(1.5), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(1.5), vjust=1)) 
# dev.off()


tiff(file.path(fig.dir, "Tair_1981-2010_PHA_PHO_PDL.tiff"), height=900, width=1600)
ggplot(data=tair.x) +
		geom_raster(aes(x=lon, y=lat, fill=tair)) +
		scale_fill_gradientn(name=expression(bold(paste("Mean Ann Temp ("^"o","C)"))), colours=c("blue4", "red3")) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=1) +
		# geom_point(data=sites, aes(x=lon, y=lat), size=8, color="darkgoldenrod3") +
		geom_point(data=sites[sites$Site %in% c("PHA", "PDL", "PUN"),], aes(label=name, x=lon, y=lat), size=15, color="darkgoldenrod3") +
		geom_text(data=sites[sites$Site %in% c("PHA", "PDL", "PUN"),], aes(label=name, x=lon-c(0.75,0.5,0.25), y=lat-c(-1,1,1)), size=15, color="darkgoldenrod3", fontface="bold") +
		scale_x_continuous(limits=range(tair.x$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(tair.x$lat), expand=c(0,0), name="Latitude") +
		coord_equal(ratio=1) +
		# ggtitle(paste("Tair", year, mon, sep=" - ")) +
		guides(fill=guide_colorbar(title.position="top")) +
        theme(plot.margin=unit(c(1,1,1,1), "lines")) +
        theme(legend.position=c(0.83, 0.2), 
              legend.direction="horizontal",
              legend.key.size=unit(4, "lines"),
              legend.text=element_text(color="black", size=rel(3.5)),
              legend.title=element_text(color="black", size=rel(3.5))) +
        theme(axis.line=element_blank(), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(4)), 
	      axis.text.y=element_text(color="black", size=rel(4)), 
	      axis.title.x=element_text(face="bold", size=rel(4), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(4), vjust=1)) 
dev.off()
# ----------------------------------------------



# ----------------------------------------------
# Mapping Biomes
# ----------------------------------------------
veg   <- raster(veg.dir)

biomes   <- data.frame(rasterToPoints(veg))
names(biomes)   <- c("lon", "lat", "biome")
biomes$biome <- as.factor(biomes$biome)
levels(biomes$biome) <- c("Forest, Hardwood", "Forest, Mixed", "Forest, Conifer", "Savanna & Shrublands", "Grasslands")
summary(biomes)


biome.colors <- c("green3", "darkseagreen3", "darkgreen", "darkorange2", "goldenrod2")

tiff(file.path(fig.dir, "PotVeg_Biomes_withSites.tiff"), height=900, width=1600)
ggplot(data=biomes) +
		geom_raster(aes(x=lon, y=lat, fill=biome)) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group), size=1) +
		geom_point(data=sites, aes(x=lon, y=lat), size=20, color="black") +
		scale_x_continuous(limits=range(biomes$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(biomes$lat), expand=c(0,0), name="Latitude") +
		coord_equal(ratio=1) +
		# ggtitle(paste("Tair", year, mon, sep=" - ")) +
		scale_fill_manual(values=biome.colors, name="Biome") +
		guides(fill=guide_legend(ncol=1)) +
        theme(plot.margin=unit(c(1,1,1,1), "lines")) +
        theme(legend.position=c(0.82, 0.2), 
              # legend.direction="horizontal",
              legend.key.size=unit(2.5, "lines"),
              legend.text=element_text(color="black", size=rel(2)),
              legend.title=element_text(color="black", size=rel(2))) +
        theme(axis.line=element_blank(), 
	      panel.grid.major=element_blank(), 
	      panel.grid.minor=element_blank(), 
	      panel.border=element_blank(), 
	      panel.background=element_blank(), 
	      axis.text.x=element_text(angle=0, color="black", size=rel(4)), 
	      axis.text.y=element_text(color="black", size=rel(4)), 
	      axis.title.x=element_text(face="bold", size=rel(4), vjust=-0.5),  
	      axis.title.y=element_text(face="bold", size=rel(4), vjust=1)) 
dev.off()