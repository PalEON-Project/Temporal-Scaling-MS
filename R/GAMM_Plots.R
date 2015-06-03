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

rgb.line <- function(df, site, size){
	geom_line(data=df$weights[df$weights$Site==site, ], aes(x=Year, y=fit),
	          color=rgb(abs(df$weights[df$weights$Site==site, "weight.temp"]),
                        abs(df$weights[df$weights$Site==site, "weight.co2"]),
                        abs(df$weights[df$weights$Site==site, "weight.precip"])), size=size)

}
