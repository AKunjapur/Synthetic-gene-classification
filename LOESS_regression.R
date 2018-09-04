f <- read.table("CrossKingdom_LOWESS_wo_CRISPR.csv", header = TRUE, sep = ";")
df$Category <- factor(df$FeatureOrigin, c("Natural", "Synthetic"))
p<-ggplot(data=df, aes(x=PartLength, y=GeneticDistance, group=Category, shape=Category, color=Category)) +
	theme(	panel.grid.major.x = element_blank(),
			panel.grid.major.y = element_line(size=0.1, color="black"),
			panel.background = element_rect(fill = "white"),
			plot.background = element_rect(fill = "white", colour = "white"),
			axis.text.x = element_text(color="black", size=8, angle=0, family="Arial"),
			axis.text.y = element_text(color="black", size=8, angle=0, family="Arial"),
			axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), family="Arial"),
			axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), family="Arial"),
			#axis.line = element_line(arrow = arrow(angle = 25, length = unit(0.2, "cm"),ends = "last", type = "closed")),
			axis.line = element_line(size=0.1, color="black"),
			axis.ticks.y = element_blank(),
			legend.position="top",
			legend.text = element_text(color="black", size=8, angle=0, family="Arial"),
			legend.title = element_blank(),
			legend.background = element_rect(fill = "white"),
			legend.key = element_rect(fill = "white", colour = "white"),
			plot.margin=unit(c(0,0.4,0.1,0.1),"cm")
		) +

	scale_color_manual(values=c("#15a5d5", "#e84925")) +	

	scale_fill_manual(values = c("#8edbf4","#f39f8d")) +

	scale_x_continuous(	name="Gene Length",
					breaks=seq(0,70000,1000),
					labels=seq(0,70000,1000),
					limits=c(0,70000),
					position="bottom”) +

	scale_y_continuous(	name="Genetic Distance",
					breaks=seq(0,1,0.2),
					labels=seq(0,1,0.2),
					limits=c(0,10))  +

	stat_smooth(aes(fill=Category), method="loess", formula = y ~ x, level = 0.68, span = 0.9) +
	coord_cartesian(xlim=c(0, 7000), ylim=c(0, 1))
