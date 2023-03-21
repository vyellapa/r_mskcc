args=commandArgs(TRUE)
s1=args[1]
s2=args[2]
library(ggplot2)
library(reshape)
setwd(s2)

n=read.table(s1,sep=' ',header=T)
m=melt(n)

p=ggplot(m,aes(Sample,value,fill=factor(variable)))+geom_bar(stat="identity",alpha=0.85)+theme_bw(base_size=40)+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.text = element_text(size=37))

#p=ggplot(m,aes(reorder(Sample,value,sum),value,fill=factor(variable)))+geom_bar(stat="identity",alpha=0.85)+theme_bw(base_size=40)+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.text = element_text(size=37))+geom_hline(yintercept=500)
