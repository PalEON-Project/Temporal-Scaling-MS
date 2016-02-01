ensemble.wts.graph <- ensemble.wts0
ensemble.wts.graph <- ensemble.wts.graph[!ensemble.wts.graph$Extent=="1985-2010",]
summary(ensemble.wts.graph)

# wts.sum <- abs(ensemble.wts.graph$weight.tair.10.adj) + abs(ensemble.wts.graph$weight.precipf.10.adj) + abs(ensemble.wts.graph$weight.CO2.10.adj)
# ensemble.wts.graph[,c("weight.tair.10.adj","weight.precipf.10.adj", "weight.CO2.10.adj")] <- ensemble.wts.graph[,c("weight.tair.10.adj","weight.precipf.10.adj", "weight.CO2.10.adj")]/wts.sum
ensemble.wts.graph[is.na(ensemble.wts.graph$weight.tair.10.adj   ),"weight.tair.10.adj"   ] <- 0
ensemble.wts.graph[is.na(ensemble.wts.graph$weight.precipf.10.adj),"weight.precipf.10.adj"] <- 0
ensemble.wts.graph[is.na(ensemble.wts.graph$weight.CO2.10.adj    ),"weight.CO2.10.adj"    ] <- 0
ensemble.wts.graph$ci.max <- ifelse(ensemble.wts.graph$fit.full.rel.10.hi>1.5,1.5,ensemble.wts.graph$fit.full.rel.10.hi)
ensemble.wts.graph$ci.min <- ifelse(ensemble.wts.graph$fit.full.rel.10.lo<0.5,0.5,ensemble.wts.graph$fit.full.rel.10.lo)
ensemble.wts.graph$fit.graph <- ifelse(ensemble.wts.graph$fit.full.rel.10<0.5,NA,ensemble.wts.graph$fit.full.rel.10)

pdf(file.path(fig.dir, "Ensemble_Drivers_Time_Region_1500-2010_Decadal_byModel.pdf"), width=11, height=8.5)
for(m in unique(ensemble.wts.graph$Model)){
print(
ggplot(data= ensemble.wts.graph[ensemble.wts.graph$Model==m ,]) + facet_grid(Extent~Model, scales="free_y") +
  scale_x_continuous(limits=c(1500,2010), expand=c(0,0), breaks=seq(round(min(ensemble.wts.graph$Year), -2), round(max(ensemble.wts.graph$Year), -2), by=100)) +
  geom_ribbon(data= ensemble.wts.graph[ensemble.wts.graph$Model==m,], aes(x=Year, ymin=fit.full.rel.10.lo*100, ymax=fit.full.rel.10.hi*100), alpha=0.35) +
	geom_line(data= ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent=="1901-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent=="1901-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent =="1901-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent =="1901-2010","weight.precipf.10.adj"])), size=3) +
	geom_line(data= ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent=="850-2010",], aes(x=Year, y=fit.full.rel.10*100),
	          color=rgb(abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent=="850-2010","weight.tair.10.adj"]),
                        abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent =="850-2010","weight.CO2.10.adj"]),
                        abs(ensemble.wts.graph[ensemble.wts.graph$Model==m & ensemble.wts.graph$Extent =="850-2010","weight.precipf.10.adj"])), size=3) +
 	geom_hline(yintercept=100, linetype="dashed") +
# 	scale_y_continuous(name=expression(bold(paste("Relative NPP (%)"))), expand=c(0,0)) +
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