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


rho=read.table(r, sep="\t", header=T)
pur1=rho[rho$is.best=="TRUE",]$rho[2]
pur2=pur1/2;

ploidy=rho[rho$is.best=="TRUE",]$ploidy[2]
label=rownames(rho[rho$is.best=="TRUE",])[2]



name=gsub(".caveman.tsv.gz","",strsplit(i,"/")[[1]][3])
outFile=paste(name,"v2","pdf",sep=".")
p=ggplot(caveman, aes(caveman$TARGET_VAF))+geom_density(fill="blue", alpha=0.6)+ggtitle(name)+annotate("text", label=paste("frac_genome_rho =", pur1),x=pur2,y=2)+geom_vline(xintercept = pur2, color="red")+annotate("text", label=paste("N =",nrow(caveman)),x=0.1,y=3)+annotate("text", label=paste("frac_genome_ploidy =", ploidy, x=pur2,y=1))+annotate("text", label=label,x=0.1,y=4)
ggsave(p, file=outFile, width=6, height=6)





}

