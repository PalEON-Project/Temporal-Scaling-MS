# Script to make some maps of climate averages over the paleon domain

# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(raster)
library(ncdf4); library(ggplot2); library(grid)

met.dir  <- "/projectnb/dietzelab/paleon/met_regional/phase2_met_regional_v2_monthly"
met.out  <- "/projectnb/dietzelab/paleon/met_regional/phase2_met_regional_climate_means"
if(!dir.exists(met.out)) dir.create(met.out)
setwd(met.out)

met.slices <- data.frame(start=c(0851, 1851, 1981), end=c(0880, 1880, 2010))

met.vars <- c("tair", "precipf")
paleon.states <- map_data("state")

# ----------------------------------------------


# ----------------------------------------------
# Finding and graphing climate means
# ----------------------------------------------

tair <- stack(file.path(met.dir, "tair.nc"))
precip <- stack(file.path(met.dir, "precipf.nc"))

yrs <- unlist(strsplit(names(tair), "[.]"))
yrs <- yrs[seq(1, length(yrs), by=3)]
yrs <- as.numeric(substr(yrs, 2, nchar(paste(yrs))))


for(i in 1:nrow(met.slices)){
tair.mean   <- mean(tair  [[which(yrs>=met.slices[i,"start"] & yrs<=met.slices[i,"end"])]])
# filename=paste0("PalEON_Met_Normals_tair_", met.slices[i,"start"], "-", met.slices[i,"end"]))
tair.mean   <- calc(tair  [[which(yrs>=met.slices[i,"start"] & yrs<=met.slices[i,"end"])]], fun=mean, filename=file.path(met.out, paste0("PalEON_Met_Normals_tair_", met.slices[i,"start"], "-", met.slices[i,"end"])))

precip.mean <- calc(precip[[which(yrs>=met.slices[i,"start"] & yrs<=met.slices[i,"end"])]], fun=mean, filename=file.path(met.out, paste0("PalEON_Met_Normals_precipf_", met.slices[i,"start"], "-", met.slices[i,"end"])))


tair.x   <- data.frame(rasterToPoints(tair.mean))
precipf.x <- data.frame(rasterToPoints(precip.mean))
names(tair.x)   <- c("lon", "lat", "tair")
names(precipf.x) <- c("lon", "lat", "precipf")

# Convert to celcius & mm/yr
tair.x$tair       <- tair.x$tair-273.15
precipf.x$precipf <- precipf.x$precipf*1*60*60*24*365.25 # 1 sec * 60 sec/min * 60 min/hr * 24 hr/day * 365.25 days/yr


	plot.tair <- ggplot(data=tair.x) +
		geom_raster(aes(x=lon, y=lat, fill=tair)) +
		scale_fill_gradientn(colours=c("gray50", "red3")) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
		scale_x_continuous(limits=range(tair.x$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(tair.x$lat), expand=c(0,0), name="Latitude") +
		# ggtitle(paste("Tair", year, mon, sep=" - ")) +
        theme(plot.margin=unit(c(0,0,0,0), "lines")) +
        # theme(legend.position="bottom", 
              # legend.direction="horizontal") +
        theme(panel.background=element_blank(), 
              axis.title.y=element_blank(),
              axis.title.x=element_blank()) +
		coord_equal(ratio=1)

	plot.precipf <- ggplot(data=precipf.x) +
		geom_raster(aes(x=lon, y=lat, fill=precipf)) +
		scale_fill_gradientn(colours=c("gray50", "blue3")) +
		geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
		scale_x_continuous(limits=range(tair.x$lon), expand=c(0,0), name="Longitude") +
		scale_y_continuous(limits=range(tair.x$lat), expand=c(0,0), name="Latitude") +
		# ggtitle(paste("Precipf", year, mon, sep=" - ")) +
        theme(plot.margin=unit(c(0,0,0,0), "lines")) +
        # theme(legend.position=c(0.75, 0.1), legend.direction="horizontal") +
        theme(panel.background=element_blank()) +
		coord_equal(ratio=1)


# Plot tair
	pdf(file.path(met.out, paste0("PalEON_Met_Normals_tair_", met.slices[i,"start"], "-", met.slices[i,"end"], ".pdf")), width=16, height=9)
		print(plot.tair)
	dev.off()

	tiff(file.path(met.out, paste0("PalEON_Met_Normals_tair_", met.slices[i,"start"], "-", met.slices[i,"end"], ".tiff")), width=3200, height=1800)
		print(plot.tair)
	dev.off()

# Plot precipf
	pdf(file.path(met.out, paste0("PalEON_Met_Normals_precipf_", met.slices[i,"start"], "-", met.slices[i,"end"], ".pdf")), width=16, height=9)
		print(plot.precipf)
	dev.off()

	tiff(file.path(met.out, paste0("PalEON_Met_Normals_precipf_", met.slices[i,"start"], "-", met.slices[i,"end"], ".tiff")), width=3200, height=1800)
		print(plot.precipf)
	dev.off()
}
# ----------------------------------------------
