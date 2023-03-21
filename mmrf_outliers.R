library(tidyr)
library(dplyr)
library(GenomicRanges)
a=read.table("/Users/yellapav/Desktop/record/MMRF/MMRF_CoMMpass_IA13a_E74GTF_Salmon_V7.2_Filtered_Gene_TPM.txt", header=T,sep="\t")
aa=a
rownames(a)=a[,c(1)]
a=a[,-c(1)]
anno=read.table("~/local/resources/Homo_sapiens.GRCh37.75.ensg.gname.bed", header = F,sep="\t")
cgc=read.table("~/local/resources/cancer_genes.txt", header = F,sep="\t")

a$max=apply(a,1,max)

a=a[a$max>2,]
dim(a)

a=a[,!(colnames(a) %in% c("max"))]
use=log(a+1)
a$sd=apply(a,1,sd)
a=a[a$sd>2,]
a=a[,!(colnames(a) %in% c("sd"))]
dim(a)
l=log(a+1)
l$sd=apply(l,1,sd)
ll=l
l=l[l$sd>0.5,]
l=l[order(-l$sd),]
l$ensg=rownames(l)
m=merge(anno,l,by.x="V1", by.y="ensg", all.y = TRUE)
mm=merge(m,cgc,by.x="V2",by.y="V1")
mmm=mm[mm$sd>0.8,]
mmm[,c("sd","V2","V1")]
mmmm=mmm[,!(colnames(mmm) %in% c("sd","V1"))]


ggplot(mmmm,aes(V2,value))+geom_boxplot(outlier.colour ="blue")+geom_jitter(width=0.25,alpha=0.3,size=0.5)+ylab("log(TPM+1)")+xlab("")+
  theme_bw(base_size = 15)+ 
  scale_fill_brewer(palette = "Set1")+
  theme( legend.position="none", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_blank(), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))

#l=l[,!(colnames(l) %in% c("sd"))]
#l=a[a$sd>2,]
#a=a[,!(colnames(a) %in% c("sd"))]

mmrf=read.table("/Users/yellapav/Downloads/108942/all_lesions.conf_90.txt",sep="\t",header=T)
colnames(mmrf)=gsub(" ","_",colnames(mmrf))
mmrf=mmrf[1:154,]

intersected=read.table("/Users/yellapav/Downloads/108942/insersected.bed",sep="\t",header=F)
intersected1=unique(intersected[,c(1,2,3,4)])
intersected1$V1=gsub("chr","",intersected1$V1)
intersected1$V2=intersected1$V2-2000000
intersected1$V3=intersected1$V3+2000000

intersected1 = intersected1 %>% mutate(V2 = ifelse(V2<0, 1, V2))



genes=read.table("~/local/resources/GRCh37.e75.gene_boundaries.bed",sep="\t",header=F)
genes=separate(genes,V4,sep=";",c("V4","V5","V6"))
head(intersected)

gr.bins = GRanges(seqnames=Rle(genes$V1), IRanges(genes$V2, genes$V3), ensg=genes$V6)
gr.sv = GRanges(seqnames=Rle(as.character(intersected1$V1)), IRanges(intersected1$V2, intersected1$V3))

overlapGenes <- findOverlaps(gr.bins, gr.sv)
results = data.frame(intersected1[subjectHits(overlapGenes),], genes[queryHits(overlapGenes),])
head(results)
m=c(colnames(use), c("Unique.Name","Descriptor","Wide.Peak.Limits","Peak.Limits","Region_Limits","q.values","Residual q values after removing segments shared with higher peaks","Broad.or.Focal","Amplitude.Threshold"))

mmmrf=mmrf[,colnames(mmrf) %in% m]


new <- data.frame(feature=character(), 
                  genus=character(), 
                  ENSG0=character(),
                  meanx=numeric(),
                  xq1=numeric(),
                  xq3=numeric(),
                  meany=numeric(),
                  yq1=numeric(),
                  yq3=numeric(),
                  lenx=numeric(),
                  leny=numeric(),
                  p.value=numeric(), 
                  p10.value=numeric(),
                  p10.wilcox=numeric(),
                  stringsAsFactors=FALSE)

for(i in 1:nrow(mmmrf)){
#  for(i in 1:2){
  sub=mmmrf[i,c(8:641)]
  feature=(gsub(" ","_",mmmrf[i,1]))
  pos=colnames(sub)[sub>0]
  neg=colnames(sub)[sub==0]
  gaga=results[results$V4==feature,]
  for(j in unique(gaga$V5)) {
    V6=unique(gaga[gaga$V5==j,c("V6")])
    sub.sub=use[rownames(use) %in% j,]
    if(nrow(sub.sub)>0){
    sub.sub.pos=as.vector(sub.sub[,colnames(sub.sub) %in% pos])
    sub.sub.neg=as.vector(sub.sub[,colnames(sub.sub) %in% neg])
    t=t.test(sub.sub.pos,sub.sub.neg)
    lenx=length(sub.sub.pos)
    leny=length(sub.sub.neg)
    
    p10=0;
    w10=0;
    for(k in 1:20){
      lenx3=round(0.33*lenx)
      sub.sub.pos3=sub.sub.pos[sample(lenx,lenx3)]
      
      leny3=round(0.33*leny)
      sub.sub.neg3=sub.sub.neg[sample(leny,leny3)]
      
      t10=t.test(sub.sub.pos3,sub.sub.neg3)
      w=wilcox.test(as.numeric(as.vector(sub.sub.neg3)),as.numeric(as.vector(sub.sub.pos3)))
      p10=p10+t10$p.value
      w10=w10+w$p.value
      
    }
    p10=p10/20;
    w10=w10/20;
    
    xq1=as.numeric(strsplit(as.character(as.vector(summary(t(sub.sub.pos))))[2],":")[[1]][2])
    xq3=as.numeric(strsplit(as.character(as.vector(summary(t(sub.sub.pos))))[5],":")[[1]][2])
    yq1=as.numeric(strsplit(as.character(as.vector(summary(t(sub.sub.neg))))[2],":")[[1]][2])
    yq3=as.numeric(strsplit(as.character(as.vector(summary(t(sub.sub.neg))))[5],":")[[1]][2])
    #X is positive
    print(c(feature,V6,rownames(sub.sub),j,t$p.value))
    new <- rbind(new, data.frame(feature=feature, genus=V6, ENSG0=j,meanx=(t$estimate)[1],xq1=xq1,xq3=xq3,meany=(t$estimate)[2],yq1=yq1,yq3=yq3,lenx=lenx,leny=leny,p.value=t$p.value,p10.value=p10, p10.wilcox=w10))
    }
  }
}

unique(gaga$V6)
n=inner_join(new,results,by=c("feature" = "V4"))
n=n[n$genus==n$V6,]
n$qval=p.adjust(n$p.value, method = "fdr")
n.fil=(n[(n$qval<0.05 & (n$meanx>0.75 | n$meany>0.75) & (abs(n$meanx-n$meany)>0.25)),])
n.fil=n.fil[n.fil$genus==n.fil$V6,]
#n.fil$V1.1==14 & n.fil$V2.1<103885201
#108850199
#n.fil=n.fil[n.fil$lenx<100,]
n.fil$diff=abs(n.fil$meanx-n.fil$meany)
nrow(n.fil)
View(n.fil)

write.table(n.fil, "~/Desktop/mmrf_gistic/gene_expression_corr.txt",append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE, quote=FALSE)
