print( 
ggplot() + facet_grid(Data.Type ~ Site, scales="free", space="free") +
	scale_x_continuous(limits=c(1900,2010), expand=c(0,0), breaks=c(1800, 1850, 1900, 1950, 2000)) +
 	geom_ribbon(data= ensemble.wts.site[,], aes(x=Year, ymin=Y.rel.10.lo*100, ymax=Y.rel.10.hi*100), alpha=0.5) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4)	+
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4)	+                        
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Data.Type=="TreeRingRW",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="TreeRingRW","weight.precipf.10"])), size=4)	+


	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Data.Type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="TreeRingNPP",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="TreeRingNPP","weight.precipf.10"])), size=4) +

	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHO" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHO"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PHA" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PHA"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PUN" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PUN"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4) +
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PBL" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PBL"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4)	+
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PDL" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PDL"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4)	+                        
	geom_line(data= ensemble.wts.site[ensemble.wts.site$Site=="PMB" & ensemble.wts.site$Data.Type=="Model",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="Model","weight.tair.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="Model","weight.CO2.10"]),
                        abs(ensemble.wts.site[ensemble.wts.site$Site=="PMB"& ensemble.wts.site$Data.Type=="Model","weight.precipf.10"])), size=4)	+
	geom_hline(y=100, linetype="dashed") +
	scale_y_continuous(name=expression(bold(paste("% Mean NPP"))), expand=c(0,0)) +
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
	      panel.background=element_blank()) +
	theme(plot.title=element_text(face="bold", size=rel(2.5))) + 
	theme(strip.text=element_text(size=rel(0.75), face="bold")) +
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
	theme(axis.text.x=element_text(color="black", size=rel(1)),
		  axis.text.y=element_text(color="black", size=rel(1)), 
		  axis.title.x=element_text(size=rel(1), face="bold"),  
		  axis.title.y=element_text(size=rel(1), face="bold"),
		  axis.ticks.length=unit(-0.5, "lines"),
	      axis.ticks.margin=unit(1.0, "lines"))
)
