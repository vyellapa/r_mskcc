library(reshape2)
library(ggplot2)
library(VennDiagram)
setwd("~/Desktop/SOHN_163/NYGC_download/Project_SOH_11516_B01_SOM_WGS.2016-06-30/call_compare/elli2/")


b=read.table("a.cnv",sep="\t",header=F)
b$V22=abs(b$V22)

b$CNV="NEUTRAL"
b[which(b$V20<1.5 | b$V20>2.5 | b$V21<1.5 | b$V21>2.5),c("CNV")]="ABERRANT"


b$V11=as.character(b$V11)
b[b$V12=="translocation",]$V11="TRA"
b[b$V12=="deletion",]$V11="DEL"
b[b$V12=="tandem-duplication",]$V11="DUP"
b[b$V12=="inversion",]$V11="INV"

b[grep("=",b$V11),]$V11=b[grep("=",b$V11),]$V7
#b[b$V11==127,]$V11="INS"


b[b$V7=="CTX",]$V11="TRA"
b[b$V7=="ITX",]$V11="TRA"
b[b$V7=="TRA",]$V11="TRA"
b[b$V7=="INV",]$V11="INV"
b[b$V7=="DEL",]$V11="DEL"
b[b$V7=="DUP",]$V11="DUP"
b[b$V7=="INS",]$V11="INS"


table(b$V11,b$V19,b$CNV)
p=as.data.frame(table(b$V11,b$V19,b$CNV))

temp=b[b$V14=="BRASS",]
x=as.data.frame(table(temp$V11,temp$CNV))
x$Caller="BRASS"

temp=b[b$V15=="crest",]
y=as.data.frame(table(temp$V11,temp$CNV))
y$Caller="crest"
x=rbind(x,y)

temp=b[b$V16=="breakdancer",]
y=as.data.frame(table(temp$V11,temp$CNV))
y$Caller="breakdancer"
x=rbind(x,y)

temp=b[b$V17=="delly",]
y=as.data.frame(table(temp$V11,temp$CNV))
y$Caller="delly"
x=rbind(x,y)


temp=b[b$V18=="NOVOBREAK",]
y=as.data.frame(table(temp$V11,temp$CNV))
y$Caller="NOVOBREAK"
x=rbind(x,y)

m=melt(x)

ggplot(m,aes(Caller,value,fill=Var1))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("Count")+
  facet_wrap(~Var2)




b$LEN=0
b[b$V11=="TRA",]$LEN=0
b$V1=as.character(b$V1)
b$V4=as.character(b$V4)
b[which(b$V1==b$V4 & b$V11!="TRA"),]$LEN=abs(b[which(b$V1==b$V4 & b$V11!="TRA"),]$V5-b[which(b$V1==b$V4 & b$V11!="TRA"),]$V2)


ggplot(b,aes(fill=as.factor(V19),LEN))+geom_histogram(bins=100)+facet_wrap(V11~CNV,scales="free")+
theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("Count")


ggplot(b,aes(V11,LEN,fill=as.factor(V19)))+geom_boxplot()+geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+facet_wrap(~CNV,scales="free")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("Count")+coord_cartesian(ylim=c(0,1000000))


#temp=b[which(b$V14=="BRASS" & !(is.na(b$V14))),]
#x=as.data.frame(table(b$V11,b$LEN,b$CNV,b$V19))
#x$Caller="BRASS"





ggplot(p,aes(Var1,Freq,fill=factor(Var2)))+geom_bar(stat="identity")+
 theme_bw()+ 
 scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
   panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
       xlab("")+ylab("Count")+
       facet_wrap(~Var3)


venn.diagram(list(T1.groups.filtered.bedpe= c(1:25,1001:7606), NYGC.SV = c(1:25, 101:169)), main="E-H-109099-T1",fill = c("#66CC99","#99CCFF"),alpha = c(0.55, 0.55),fontface=3 ,cex = 1.4,cat.cex=1.4,cat.fontface = 1,lty =0,cat.pos=353,filename = "~/Desktop/SOHN_163/NYGC_download/Project_SOH_11516_B01_SOM_WGS.2016-06-30/call_compare/elli2/figs/E-H-109099-T1.sv.png")
venn.diagram(list(T1.groups.filtered.bedpe= c(1:40,1001:11235), NYGC.SV = c(1:40, 101:297)), main="E-H-109099-T2",fill = c("#66CC99","#99CCFF"),alpha = c(0.55, 0.55),fontface=3 ,cex = 1.4,cat.cex=1.4,cat.fontface = 1,lty =0,cat.pos=353,filename = "~/Desktop/SOHN_163/NYGC_download/Project_SOH_11516_B01_SOM_WGS.2016-06-30/call_compare/elli2/figs/E-H-109099-T2.sv.png")
venn.diagram(list(T1.annot.bedpe= c(1:21,101:202), NYGC.SV.all = c(1:21, 1001:2468)), main="CUL-65-Diagnosis",fill = c("#66CC99","#99CCFF"),alpha = c(0.55, 0.55),fontface=3 ,cex = 1.4,cat.cex=1.4,cat.fontface = 1,lty =0,cat.pos=353,filename = "~/Desktop/SOHN_163/NYGC_download/Project_SOH_11516_B01_SOM_WGS.2016-06-30/call_compare/elli2/figs/CUL65-Diagnosis.sv.png")
venn.diagram(list(T2.annot.bedpe= c(1:31,101:354), NYGC.SV.all = c(1:31, 1001:3488)), main="CUL-65-Relapse",fill = c("#66CC99","#99CCFF"),alpha = c(0.55, 0.55),fontface=3 ,cex = 1.4,cat.cex=1.4,cat.fontface = 1,lty =0,cat.pos=353,filename = "~/Desktop/SOHN_163/NYGC_download/Project_SOH_11516_B01_SOM_WGS.2016-06-30/call_compare/elli2/figs/CUL65-Relapse.sv.png")


e5=read.table("E5.plot",sep="\t",header=F)
#e5=e5[grep("099",e5$V1,invert=T),]
ggplot(e5,aes(x=reorder(V1,-V2),y=V2))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("")+geom_hline(yintercept = 6,color="blue", aes(color="blue"))+annotate("text",x=5,y=8,label="y=6")






merged=read.table("files/merged.bedpe",sep="\t",header=F)

colnames(merged)=c("Sname","Chr1","Start1","End1","Chr2","Start2","End2","NA","NA1","Starnd1","Strand2","SV_TYPE","NA2","NA3","BRASS","crest","breakdancer","delly","NOVOBREAK","NUM_CALLERS","BP1_AVG_MAPQ","BP1_AVG_MATCH_LEN","BP1_SUPPORT_READS","BP2_AVG_MAPQ","BP2_AVG_MATCH_LEN","BP2_SUPPORT_READS","BP1_CNV","BP2_CNV")

ggplot(merged,aes(Sname,BP1_AVG_MAPQ,fill=factor(NUM_CALLERS)))+geom_boxplot(alpha=0.8)+
  geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")


ggplot(merged,aes(Sname,BP2_AVG_MAPQ,fill=factor(NUM_CALLERS)))+geom_boxplot(alpha=0.8)+
  geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+geom_hline(yintercept = 10,color="blue")+annotate("text",x=1,y=10,label="y=10")

ggplot(merged,aes(Sname,BP2_SUPPORT_READS,fill=factor(NUM_CALLERS)))+geom_boxplot(alpha=0.8)+
  geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+coord_cartesian(ylim=c(0,1000))+geom_hline(yintercept = 800,color="blue")+annotate("text",x=1,y=809,label="y=800")


ggplot(merged,aes(Sname,BP2_AVG_MATCH_LEN,fill=factor(NUM_CALLERS)))+geom_boxplot(alpha=0.8)+
  geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+geom_hline(yintercept = 50,color="blue")+annotate("text",x=1,y=57,label="y=50")


ggplot(merged,aes(Sname,BP2_CNV,fill=factor(NUM_CALLERS)))+geom_boxplot(alpha=0.8)+
  geom_point(position = position_dodge(width=0.75),alpha=0.25)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")




svnum=read.table("files/sv_number.plot",sep=" ",header=F)
svnum=melt(svnum)
ggplot(svnum,aes(V2,value,fill=V3))+geom_bar(stat="identity",position = position_dodge())+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")


prePost=read.table("files/pre_post.plot",sep=" ",header=T)
m=melt(prePost)
ggplot(m,aes(SAMPLE,value,fill=variable))+geom_bar(stat="identity", position = position_dodge())+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")


ggplot(plot, aes(V2,V1))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("Number")


ggplot(a,aes(x=reorder(V1,-V4),y=V4))+geom_bar(stat="identity", aes(fill=factor(1)))+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+
  xlab("")+ylab("Frequency %")


col = RColorBrewer::brewer.pal(n = 12, name = 'Paired')
ggplot(m, aes(Translocation,value,fill=Sample))+geom_bar(stat="identity", position = position_dodge())+facet_wrap(~variable,scales = "free")+
  theme_bw(base_size = 18)+ 
  scale_fill_brewer(palette = "Set2")+ scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())

dplot=read.table("~/Desktop/dplot", header = TRUE, stringsAsFactors = F, sep="\t")
dm=melt(dplot)
ggplot(dm, aes(Translocation,value,fill=SAMPLE))+geom_bar(stat="identity", position = position_dodge())+facet_wrap(~variable,scales = "free")+
  theme_bw(base_size = 18)+ 
  scale_fill_brewer(palette = "Set2")+ scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("Delly")


b=read.table("E-H-131564-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131565-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131568-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131569-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131570-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131571-T1-1-D1-1_flt.plot",sep="\t",header=F)
b=read.table("E-H-131572-T1-1-D1-1_flt.plot",sep="\t",header=F)


ggplot(b,aes(V2,V9))+geom_point(color="blue", alpha=0.5)+facet_wrap(~V1,scales="free_x")+
  theme_bw(base_size = 10)+ 
  scale_colour_brewer(palette = "Set1")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")

ggplot(m,aes(variable,value,fill=ANNOT))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=6,position = position_dodge(width=0.5),aes(group=ANNOT, color=factor(ANNOT)), alpha=0.3)+facet_wrap(~variable, scales="free")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme( legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+geom_hline(yintercept = 30)+geom_hline(yintercept = 60)

ggplot(m,aes(variable,value,fill=ANNOT))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=6,position = position_dodge(width=0.5),aes(group=ANNOT, color=factor(ANNOT)), alpha=0.3)+facet_grid(variable~Translocation, scales="free")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
         axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
         panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+geom_hline(yintercept = 30)+geom_hline(yintercept = 60)


t15=read.table("~/Desktop/p219/delly/p219_trans_all.plot",header = F)
colnames(t15)=colnames(all)
mt15=melt(t15)

ggplot(mt15,aes(variable,value,fill=ANNOT))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=6,position = position_dodge(width=0.5),aes(group=ANNOT, color=factor(ANNOT)), alpha=0.3)+facet_grid(variable~Translocation, scales="free")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+geom_hline(yintercept = 30)+geom_hline(yintercept = 60)


ggplot(m222,aes(CHR2.CHR1,value,fill=variable))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=4,position = position_dodge(width=0.5),aes(group=variable), alpha=0.3)+facet_grid(variable~CHR2.CHR1, scales="free")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+geom_hline(yintercept = 30)+geom_hline(yintercept = 60)




t15=read.table("/Users/yellapav/Desktop/SV_test//Users/yellapav/Desktop/SV_test/p251_50_plot_non14.txt",header = T) 
ggplot(m,aes("",value,fill=ANNOT))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=4,position = position_dodge(width=0.8),aes(group=ANNOT), alpha=0.3)+facet_wrap(  ~  variable, scales="free", nrow=1)+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+geom_hline(yintercept = 30)+geom_hline(yintercept = 60)

ggplot(plot, aes(V1,V2))+geom_bar(stat="identity")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("Median Coverage (X)")

