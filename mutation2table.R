setwd("/Users/yellapav/Desktop/SOHN_163/cave_pindel")
cgc=read.table("~/local/resources/cgc_syn.tsv",sep="\t")
ct1DP=read.table("E-H-109099-T1-1-D1-1_vs_E-H-109099-N1-1-D1-1.caveman.tsv.gz",header=T,sep="\t",quote="~")
ct2DP=read.table("E-H-109099-T2-1-D1-1_vs_E-H-109099-N1-1-D1-1.caveman.tsv.gz",header=T,sep="\t",quote="~")
ct1DP=ct1DP[,c(1,11,12,22)]
ct2DP=ct2DP[,c(1,11,12,22)]

ct1=read.table("E-H-109099-T1-1-D1-1.filter.caveman.rc.gm.cosmic.maf",header=T,sep="\t",quote="~")
ct2=read.table("E-H-109099-T2-1-D1-1.filter.caveman.rc.gm.cosmic.maf",header=T,sep="\t",quote="~")


ct1=merge(ct1,ct1DP,by.x="variant_id",by.y="ID_VARIANT",all.x=T)
ct2=merge(ct2,ct2DP,by.x="variant_id",by.y="ID_VARIANT",all.x=T)
#ct1=ct1[which(ct1$FILTER=="PASS" & ct1$ASRD>0.93),]
#ct2=ct2[which(ct2$FILTER=="PASS" & ct2$ASRD>0.93),]



pt1=read.table("E-H-109099-T1-1-D1-1.v2.pindel.gm.maf",header=T,sep="\t",quote="~")
pt2=read.table("E-H-109099-T2-1-D1-1.v2.pindel.gm.maf",header=T,sep="\t",quote="~")
pp1=read.table("E-H-109099-T1-1-D1-1.PPNP",header=T,sep="\t")
pp2=read.table("E-H-109099-T2-1-D1-1.PPNP",header=T,sep="\t")

pt1=merge(pt1,pp1,by.x="variant_id",by.y="ID",all.x=T)
pt2=merge(pt2,pp2,by.x="variant_id",by.y="ID",all.x=T)

pt1$DP=(pt1$GEN.1..PR+pt1$GEN.1..NR)
pt1$VAF=(pt1$GEN.1..PP+pt1$GEN.1..NP)/(pt1$DP)
pt1$VAF=round(pt1$VAF, digits = 2)


pt2$DP=(pt2$GEN.1..PR+pt2$GEN.1..NR)
pt2$VAF=(pt2$GEN.1..PP+pt2$GEN.1..NP)/(pt2$DP)
pt2$VAF=round(pt2$VAF, digits = 2)

ct1=ct1[,c(17,6,7,12,14,122,123,2,38,124,36)]
pt1=pt1[,c(17,6,7,12,14,136,135,2,38,10,36)]
colnames(pt1)=colnames(ct1)
pt1=rbind(pt1,ct1)


ct2=ct2[,c(17,6,7,12,14,122,123,2,38,124,36)]
pt2=pt2[,c(17,6,7,12,14,136,135,2,38,10,36)]
colnames(pt2)=colnames(ct2)
pt2=rbind(pt2,ct2)
pt2=rbind(pt1,pt2)

#pt2$EFFECT=gsub("intron_variant","Intron",pt2$EFFECT)
#pt2$EFFECT=gsub("non_synonymous_codon","Missense",pt2$EFFECT)
pt2=pt2[pt2$Hugo_Symbol %in% cgc$V1,]
pt2[pt2$EFFECT=="intron_variant",]$EFFECT="Intron"
pt2[pt2$EFFECT=="non_synonymous_codon",]$EFFECT="Missense_Mutation"


write.table(pt2,file="T2.tsv", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
