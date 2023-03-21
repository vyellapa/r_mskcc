###Venn
library(VennDiagram)
venn.diagram(list(MSK = c(1:942,1001:1041), NYGC = c(1:942,2001:2351)),fill = c("#0000FF", "#FF00CC"),alpha = c(0.50, 0.50),fontface=7 ,cex = 1.5,cat.cex=1.5,cat.fontface = 8,lty =0,cat.pos=180,filename = "~/Desktop/T1.png");
venn.diagram(list(MSK = c(1:1507,10001:10071), NYGC = c(1:1507,2001:2676)),fill = c("#0000FF", "#FF00CC"),alpha = c(0.50, 0.50),fontface=7 ,cex = 1.5,cat.cex=1.5,cat.fontface = 8,lty =0,cat.pos=180,filename = "~/Desktop/T2.png");
venn.diagram(list(MSK = c(1:40,101:295), NYGC = c(1:40,501:644)),fill = c("#0000FF", "#FF00CC"),alpha = c(0.50, 0.50),fontface=7 ,cex = 1.5,cat.cex=1.5,cat.fontface = 8,lty =0,cat.pos=180,filename = "~/Desktop/T1_Indel.png");
venn.diagram(list(MSK = c(1:52,101:455), NYGC = c(1:52,501:734)),fill = c("#0000FF", "#FF00CC"),alpha = c(0.50, 0.50),fontface=7 ,cex = 1.5,cat.cex=1.5,cat.fontface = 8,lty =0,cat.pos=180,filename = "~/Desktop/T2_Indel.png");



library(ggplot2)
all=read.table("/Users/yellapav/Desktop/elli/all.plot",sep="\t",header=T)
all=all[all$chr!="X",]
all$MMQS_DIFF=abs(all$count-all$count.1)

ggplot(all, aes(MAPQ.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+coord_cartesian(c(55,60))+xlab("MAPQ of ALT supporting reads")
ggplot(all, aes(BASEQ.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()
ggplot(all, aes(Plus_strand.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Reads on Pos Strand")
ggplot(all, aes(Minus_strand.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Reads on Neg Strand")
ggplot(all, aes(pos_as_fraction.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Pos of variant as fraction of read")
ggplot(all, aes(mismatches_as_fraction.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+coord_cartesian(c(0,0.07))+xlab("Mismatches as fraction of read")
ggplot(all, aes(q2_containing_reads.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+coord_cartesian(c(0,200))+xlab("# Q2 containing reads")
ggplot(all, aes(avg_clipped_length.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+coord_cartesian(c(130,150))
ggplot(all, aes(avg_distance_to_3p_end.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Distance to 3' end of the read")
ggplot(all, aes(count.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Distance to 3' end of the read")
ggplot(all, aes(MMQS.1,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("Alt MMQS")+coord_cartesian(c(0,150))
ggplot(all, aes(MMQS_DIFF,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("MMQS Diff")+coord_cartesian(c(0,500))
ggplot(all, aes(RLEN_DIFF,fill=GROUP))+geom_density(alpha=0.5, linetype="blank")+theme_bw()+xlab("RLEN Diff")

u=all[all$GROUP=="UNIQUE",]
nrow(u[u$count.1>5,])
nrow(u)

ggplot(all, aes(count,count.1,color=GROUP))+geom_point(alpha=0.7,size=0.75, linetype="blank")+theme_bw()+xlab("Ref supporting read")+ylab("Alt supporting reads")+coord_cartesian(c(0,300),c(0,300))