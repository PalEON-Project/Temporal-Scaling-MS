# Handy figures for plotting the gamm outputs

plot.time.ci <- function(df, model.name, var, scale, yrs=c(850,2010)){
	col.model <- model.colors[model.colors$Model.Order %in% unique(df$data$Model.Order),"color"]
	ggplot(data=df$ci.response) + facet_wrap(~Site) + theme_bw() +
		geom_line(data=df$data, aes(x=Year, y=response)) +
		geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), alpha=0.5, fill=col.model) +
		geom_line(aes(x=Year, y=mean), size=0.35, color= col.model) +
		scale_x_continuous(limits=yrs) +
		# scale_fill_manual(values=col.model) +
		# scale_color_manual(values=col.model) +		
		labs(title=paste0(var, ", ", ifelse(scale=="", "  1", scale), "yr: ", model.name), x="Year", y=var)
}

rgb.line <- function(df, size, y, red, green, blue){
	df.temp       <- df[, c("Year", "Site", "Resolution", y, red, green, blue)]
	names(df.temp) <- c("Year", "Site", "Resolution", "y", "red", "green", "blue")
	 
	geom_line(data=df.temp, aes(x=Year, y=y),
	          color=rgb(abs(df.temp$red),
                        abs(df.temp$green),
                        abs(df.temp$blue)), size=size)

}


plot.weights.time <- function(df, y="fit.full", red="weight.tair", green="weight.CO2", blue="weight.precipf", xmin=min(df$Year), xmax=max(df$Year), breaks=seq(min(df$Year), max(df$Year), by=(max(df$Year)-min(df$Year))/5)[2:5], plot.labs=labs(x="Year"), formats=theme_bw()){
ggplot(data=df) + facet_grid(Site~Resolution, scales="free") +
 	geom_line(data= mod.out$data[,], aes(x=Year, y=NPP), alpha=0.5, size=1.5) +
	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PHA",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PHA",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PHA",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PHA",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PHA",], size=2, y=y, red=red, green=green, blue=blue) +

	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PHO",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PHO",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PHO",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PHO",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PHO",], size=2, y=y, red=red, green=green, blue=blue) +

	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PUN",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PUN",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PUN",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PUN",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PUN",], size=2, y=y, red=red, green=green, blue=blue) +

	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PBL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PBL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PBL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PBL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PBL",], size=2, y=y, red=red, green=green, blue=blue) +

	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PDL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PDL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PDL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PDL",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PDL",], size=2, y=y, red=red, green=green, blue=blue) +

	rgb.line(df=df[df$Resolution=="t.001" & df$Site=="PMB",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.010" & df$Site=="PMB",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.050" & df$Site=="PMB",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.100" & df$Site=="PMB",], size=2, y=y, red=red, green=green, blue=blue) +
	rgb.line(df=df[df$Resolution=="t.250" & df$Site=="PMB",], size=2, y=y, red=red, green=green, blue=blue) +
	scale_x_continuous(limits=c(xmin,xmax), breaks=breaks) +
	# scale_y_continuous(limits=quantile(mod.out$data$response, c(0.01, 0.99), na.rm=T)) +
	plot.labs + formats
}

