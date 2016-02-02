wts.sum <- abs(wt.terms2$weight.tair.10) + abs(wt.terms2$weight.precipf.10) + abs(wt.terms2$weight.CO2.10)
wt.terms2[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")] <- wt.terms2[,c("weight.tair.10","weight.precipf.10", "weight.CO2.10")]/wts.sum
wt.terms2[is.na(wt.terms2$weight.tair.10   ),"weight.tair.10"   ] <- 0
wt.terms2[is.na(wt.terms2$weight.precipf.10),"weight.precipf.10"] <- 0
wt.terms2[is.na(wt.terms2$weight.CO2.10    ),"weight.CO2.10"    ] <- 0

{
print(
ggplot(data= wt.terms2[wt.terms2$Model=="lpj.wsl" ,]) + facet_grid(Site~Model, scales="free_y") +
 	# geom_ribbon(data= wt.terms2[,], aes(x=Year, ymin=fit.full.rel.lo*100, ymax=fit.full.rel.hi*100), alpha=0.35) +
	geom_line(data= wt.terms2[wt.terms2$Model=="lpj.wsl" & wt.terms2$Site=="PHA",], aes(x=Year, y=Y.rel.10*100),
	          color=rgb(abs(wt.terms2[wt.terms2$Model=="lpj.wsl" & wt.terms2$Site=="PHA","weight.tair.10"]),
                        abs(wt.terms2[wt.terms2$Model=="lpj.wsl" & wt.terms2$Site=="PHA","weight.CO2.10"]),
                        abs(wt.terms2[wt.terms2$Model=="lpj.wsl" & wt.terms2$Site=="PHA","weight.precipf.10"])), size=3) +
 	# geom_hline(yintercept=100, linetype="dashed") +
	scale_x_continuous(limits=c(1850,2010), expand=c(0,0), breaks=seq(round(min(wt.terms2$Year), -2), round(max(wt.terms2$Year), -2), by=100)) +
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