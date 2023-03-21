apr=read.table("~/Desktop/apr24_v3.plot",sep="\t",header=T,stringsAsFactors = F)
head(apr)
apr=apr[,-c(3)]
mm=melt(apr)

ggplot(mm, aes(FISH_STATUS,value,color=FISH_STATUS))+geom_boxplot(outlier.size=0)+geom_point(aes(alpha=0.3,color=FISH_STATUS), size=3)+facet_wrap(~variable,scales="free")+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  theme( legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
         panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
         panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")