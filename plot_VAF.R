library(ggplot2)
for(i in list.files(pattern=".caveman.tsv.gz", recursive=TRUE)) {
print(i);


rhoFile=(paste(strsplit(i,"/")[[1]][1],"/*/*_rho_and_psi.txt")); 
rhoFile=gsub(" ","",rhoFile); 
r=list.files(pattern=rhoFile, recursive=T); 
print(r);
caveman=read.table(gzfile(i), sep="\t", header=T, stringsAsFactors = FALSE)
caveman=caveman[caveman$FILTER=="PASS",]
caveman=caveman[caveman$ASRD>0.93,]
med=median(caveman$TARGET_VAF)
#med =0.5
rho=read.table(r[1], sep="\t", header=T)
pur1=rho$rho[2];
pur2=rho$rho[2]/2;


name=gsub(".caveman.tsv.gz","",strsplit(i,"/")[[1]][3])
outFile=paste(name,"pdf",sep=".")
p=ggplot(caveman, aes(caveman$TARGET_VAF))+geom_density(fill="blue", alpha=0.6)+ggtitle(name)+annotate("text", label=paste("Tumor.Cellularity =", pur1),x=pur2,y=2)+geom_vline(xintercept = pur2, color="red")+annotate("text", label=paste("N =",nrow(caveman)),x=0.2,y=3)+annotate("text", label=paste("Ploidy =", rho$psi[2]),x=pur2,y=1)+
theme_bw()+
scale_fill_brewer(palette = "Set1")+
theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(),legend.title=element_blank(),
 panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
 panel.grid.major.x = element_blank())+
xlab("TUMOR_VAF")+ylab("")

#p=ggplot(caveman, aes(caveman$TARGET_VAF))+geom_density(fill="blue", alpha=0.6)+ggtitle(name)+annotate("text", label=paste("frac_genome_rho =", pur1),x=pur2,y=2)+geom_vline(xintercept = pur2, color="red")+annotate("text", label=paste("N =",nrow(caveman)),x=0.1,y=3)+annotate("text", label=paste("frac_genome_ploidy =", rho$psi[2]),x=pur2,y=1)
ggsave(p, file=outFile, width=6, height=6)





}

