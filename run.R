library(reshape2)
library(ggplot2)
library(VennDiagram)
library("stringr")
library(plyr)
setwd("/Users/yellapav/Desktop/SOHN_163/one_offs/SNP")


z="E-H-109100-T2-1-D1-1_vs_E-H-109100-N1-1-D1-1.snps.ids.vcf.gz"
name=strsplit(z,"_",fixed=TRUE)[[1]][1]
snp=read.table(z,sep="\t",header=F,stringsAsFactors = F)
snv=read.table("E-H-109100-T2-1-D1-1_vs_E-H-109100-N1-1-D1-1.caveman.muts.annot.vcf.gz",sep="\t",header=F,stringsAsFactors = F)


pat<-'MP=[0-9,e,\\-,\\+,\\.]+'
snp$MP=gsub("MP=","",str_extract(snp$V8, pat))

foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V10),':',fixed=TRUE)))
colnames(foo)=paste("NORMAL_",colnames(foo),sep="")
snp=cbind(snp,foo)
foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V11),':',fixed=TRUE)))
colnames(foo)=paste("TUMOR_",colnames(foo),sep="")
snp=cbind(snp,foo)
snp$CALL="SNP"


pat<-'MP=[0-9,e,\\-,\\+,\\.]+'
snv$MP=gsub("MP=","",str_extract(snv$V8, pat))

foo <- data.frame(do.call('rbind', strsplit(as.character(snv$V10),':',fixed=TRUE)))
colnames(foo)=paste("NORMAL_",colnames(foo),sep="")
snv=cbind(snv,foo)
foo <- data.frame(do.call('rbind', strsplit(as.character(snv$V11),':',fixed=TRUE)))
colnames(foo)=paste("TUMOR_",colnames(foo),sep="")
snv=cbind(snv,foo)
snv$CALL="SNV"


snp=rbind(snp,snv)
snp$NORMAL_DP=as.integer(as.character(snp$NORMAL_X2))+as.integer(as.character(snp$NORMAL_X3))+as.integer(as.character(snp$NORMAL_X4))+as.integer(as.character(snp$NORMAL_X5))+as.integer(as.character(snp$NORMAL_X6))+as.integer(as.character(snp$NORMAL_X7))+as.integer(as.character(snp$NORMAL_X8))+as.integer(as.character(snp$NORMAL_X9))
snp$TUMOR_DP=as.integer(as.character(snp$TUMOR_X2))+as.integer(as.character(snp$TUMOR_X3))+as.integer(as.character(snp$TUMOR_X4))+as.integer(as.character(snp$TUMOR_X5))+as.integer(as.character(snp$TUMOR_X6))+as.integer(as.character(snp$TUMOR_X7))+as.integer(as.character(snp$TUMOR_X8))+as.integer(as.character(snp$TUMOR_X9))
snp$TUMOR_X10=as.numeric(as.character(snp$TUMOR_X10))
snp$NORMAL_X10=as.numeric(as.character(snp$NORMAL_X10))
snp$MP=as.numeric(as.character(snp$MP))

snp$TUMOR_X2=as.numeric(as.character(snp$TUMOR_X2))
snp$TUMOR_X3=as.numeric(as.character(snp$TUMOR_X3))
snp$TUMOR_X4=as.numeric(as.character(snp$TUMOR_X4))
snp$TUMOR_X5=as.numeric(as.character(snp$TUMOR_X5))
snp$TUMOR_X6=as.numeric(as.character(snp$TUMOR_X6))
snp$TUMOR_X7=as.numeric(as.character(snp$TUMOR_X7))
snp$TUMOR_X8=as.numeric(as.character(snp$TUMOR_X8))
snp$TUMOR_X9=as.numeric(as.character(snp$TUMOR_X9))
snp$TUMOR_X10=as.numeric(as.character(snp$TUMOR_X10))


snp$A_DP=snp$TUMOR_X2+snp$TUMOR_X6
snp$C_DP=snp$TUMOR_X3+snp$TUMOR_X7
snp$G_DP=snp$TUMOR_X4+snp$TUMOR_X8
snp$T_DP=snp$TUMOR_X5+snp$TUMOR_X9

snp$NORMAL_DP=as.numeric(as.character(snp$NORMAL_DP))
snp$TUMOR_DP=as.numeric(as.character(snp$TUMOR_DP))


snp_ids=snp[snp$CALL=="SNP",]
somatic=snp[snp$CALL=="SNV",]






########### 1) SUBSET 1 #########################
subset1=snp_ids[snp_ids$NORMAL_DP>10 & snp_ids$TUMOR_DP>10 & snp_ids$NORMAL_X10==0 & snp_ids$TUMOR_X10>0.05,]



############ 2) SUBSET 1.2 ###################################
subset1.2_A=subset1[subset1$V5=="A" & subset1$TUMOR_X2>4 & subset1$TUMOR_X6>4,]
subset1.2_C=subset1[subset1$V5=="C" & subset1$TUMOR_X3>4 & subset1$TUMOR_X7>4,]
subset1.2_G=subset1[subset1$V5=="G" & subset1$TUMOR_X4>4 & subset1$TUMOR_X8>4,]
subset1.2_T=subset1[subset1$V5=="T" & subset1$TUMOR_X5>4 & subset1$TUMOR_X9>4,]
subset1.2=rbind(subset1.2_A,subset1.2_C,subset1.2_G,subset1.2_T)
subset1.2=subset1.2[subset1.2$TUMOR_X10>0.08,]

subset1.2$STATUS="SUBSET1.2"

################# SUBSET 1.3 Gunes suggestion to remove multi-allelic calls #####################################
subset1.3_GA=subset1.2[subset1.2$V4=="G" & subset1.2$V5=="A" & (subset1.2$TUMOR_X3+subset1.2$TUMOR_X7+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_GC=subset1.2[subset1.2$V4=="G" & subset1.2$V5=="C" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_GT=subset1.2[subset1.2$V4=="G" & subset1.2$V5=="T" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X3+subset1.2$TUMOR_X7)<1,]
subset1.3_AC=subset1.2[subset1.2$V4=="A" & subset1.2$V5=="C" & (subset1.2$TUMOR_X4+subset1.2$TUMOR_X8+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_AG=subset1.2[subset1.2$V4=="A" & subset1.2$V5=="G" & (subset1.2$TUMOR_X3+subset1.2$TUMOR_X7+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_AT=subset1.2[subset1.2$V4=="A" & subset1.2$V5=="T" & (subset1.2$TUMOR_X3+subset1.2$TUMOR_X7+subset1.2$TUMOR_X4+subset1.2$TUMOR_X8)<1,]
subset1.3_CA=subset1.2[subset1.2$V4=="C" & subset1.2$V5=="A" & (subset1.2$TUMOR_X4+subset1.2$TUMOR_X8+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_CG=subset1.2[subset1.2$V4=="C" & subset1.2$V5=="G" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X5+subset1.2$TUMOR_X9)<1,]
subset1.3_CT=subset1.2[subset1.2$V4=="C" & subset1.2$V5=="T" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X4+subset1.2$TUMOR_X8)<1,]
subset1.3_TA=subset1.2[subset1.2$V4=="T" & subset1.2$V5=="A" & (subset1.2$TUMOR_X3+subset1.2$TUMOR_X7+subset1.2$TUMOR_X4+subset1.2$TUMOR_X8)<1,]
subset1.3_TG=subset1.2[subset1.2$V4=="T" & subset1.2$V5=="G" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X3+subset1.2$TUMOR_X7)<1,]
subset1.3_TC=subset1.2[subset1.2$V4=="T" & subset1.2$V5=="C" & (subset1.2$TUMOR_X2+subset1.2$TUMOR_X6+subset1.2$TUMOR_X4+subset1.2$TUMOR_X8)<1,]

subset1.3=rbind(subset1.3_GA,subset1.3_GC,subset1.3_GT,subset1.3_AC,subset1.3_AG,subset1.3_AT,subset1.3_CA,subset1.3_CG,subset1.3_CT,subset1.3_TA,subset1.3_TG,subset1.3_TC)

subset1.3$STATUS="SUBSET1.3"

################### 3) SUBSET 2 ##############################
subset2=snp_ids[as.numeric(as.character(snp_ids$NORMAL_DP))==0,]


################## 4) REF-REF Calls Gunes ###############################
subset3=snp_ids[as.numeric(as.character(snp_ids$NORMAL_DP))>9,]
subset3=subset3[as.numeric(as.character(subset3$TUMOR_DP))>9,]
subset3=subset3[which(subset3$V4==subset3$V5 & subset3$NORMAL_X10>0.995  & subset3$TUMOR_X10<0.995),]
for(i in 1:nrow(subset3)) {
  #print(i)
  if(as.character(subset3[i,]$V4)=="T" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$G_DP)==subset3[i,]$A_DP) { subset3[i,]$V5="A";subset3[i,]$TUMOR_X10=subset3[i,]$A_DP/subset3[i,]$T_DP  }
  else if(as.character(subset3[i,]$V4)=="T" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$G_DP)==subset3[i,]$C_DP) { subset3[i,]$V5="C";subset3[i,]$TUMOR_X10=subset3[i,]$C_DP/subset3[i,]$T_DP  }
  else if(as.character(subset3[i,]$V4)=="T" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$G_DP)==subset3[i,]$G_DP) { subset3[i,]$V5="G";subset3[i,]$TUMOR_X10=subset3[i,]$G_DP/subset3[i,]$T_DP  }
  
  else if(as.character(subset3[i,]$V4)=="G" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$T_DP)==subset3[i,]$A_DP) { subset3[i,]$V5="A";subset3[i,]$TUMOR_X10=subset3[i,]$A_DP/subset3[i,]$G_DP  }
  else if(as.character(subset3[i,]$V4)=="G" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$T_DP)==subset3[i,]$C_DP) { subset3[i,]$V5="G";subset3[i,]$TUMOR_X10=subset3[i,]$G_DP/subset3[i,]$G_DP  }
  else if(as.character(subset3[i,]$V4)=="G" & max(subset3[i,]$A_DP,subset3[i,]$C_DP,subset3[i,]$T_DP)==subset3[i,]$T_DP) { subset3[i,]$V5="T";subset3[i,]$TUMOR_X10=subset3[i,]$T_DP/subset3[i,]$G_DP  }
  
  else if(as.character(subset3[i,]$V4)=="C" & max(subset3[i,]$A_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$A_DP) { subset3[i,]$V5="A";subset3[i,]$TUMOR_X10=subset3[i,]$A_DP/subset3[i,]$C_DP  }
  else if(as.character(subset3[i,]$V4)=="C" & max(subset3[i,]$A_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$G_DP) { subset3[i,]$V5="G";subset3[i,]$TUMOR_X10=subset3[i,]$G_DP/subset3[i,]$C_DP  }
  else if(as.character(subset3[i,]$V4)=="C" & max(subset3[i,]$A_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$T_DP) { subset3[i,]$V5="T";subset3[i,]$TUMOR_X10=subset3[i,]$T_DP/subset3[i,]$C_DP  }
  
  else if(as.character(subset3[i,]$V4)=="A" & max(subset3[i,]$C_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$C_DP) { subset3[i,]$V5="C";subset3[i,]$TUMOR_X10=subset3[i,]$C_DP/subset3[i,]$A_DP  }
  else if(as.character(subset3[i,]$V4)=="A" & max(subset3[i,]$C_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$G_DP) { subset3[i,]$V5="G";subset3[i,]$TUMOR_X10=subset3[i,]$G_DP/subset3[i,]$A_DP  }
  if(as.character(subset3[i,]$V4)=="A" & max(subset3[i,]$C_DP,subset3[i,]$G_DP,subset3[i,]$T_DP)==subset3[i,]$T_DP) { subset3[i,]$V5="T";subset3[i,]$TUMOR_X10=subset3[i,]$T_DP/subset3[i,]$A_DP  }
}




subset3_GA=subset3[subset3$V4=="G" & subset3$V5=="A" & (subset3$TUMOR_X3+subset3$TUMOR_X7+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_GC=subset3[subset3$V4=="G" & subset3$V5=="C" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_GT=subset3[subset3$V4=="G" & subset3$V5=="T" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X3+subset3$TUMOR_X7)<1,]
subset3_AC=subset3[subset3$V4=="A" & subset3$V5=="C" & (subset3$TUMOR_X4+subset3$TUMOR_X8+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_AG=subset3[subset3$V4=="A" & subset3$V5=="G" & (subset3$TUMOR_X3+subset3$TUMOR_X7+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_AT=subset3[subset3$V4=="A" & subset3$V5=="T" & (subset3$TUMOR_X3+subset3$TUMOR_X7+subset3$TUMOR_X4+subset3$TUMOR_X8)<1,]
subset3_CA=subset3[subset3$V4=="C" & subset3$V5=="A" & (subset3$TUMOR_X4+subset3$TUMOR_X8+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_CG=subset3[subset3$V4=="C" & subset3$V5=="G" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X5+subset3$TUMOR_X9)<1,]
subset3_CT=subset3[subset3$V4=="C" & subset3$V5=="T" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X4+subset3$TUMOR_X8)<1,]
subset3_TA=subset3[subset3$V4=="T" & subset3$V5=="A" & (subset3$TUMOR_X3+subset3$TUMOR_X7+subset3$TUMOR_X4+subset3$TUMOR_X8)<1,]
subset3_TG=subset3[subset3$V4=="T" & subset3$V5=="G" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X3+subset3$TUMOR_X7)<1,]
subset3_TC=subset3[subset3$V4=="T" & subset3$V5=="C" & (subset3$TUMOR_X2+subset3$TUMOR_X6+subset3$TUMOR_X4+subset3$TUMOR_X8)<1,]

subset3_MA=rbind(subset3_GA,subset3_GC,subset3_GT,subset3_AC,subset3_AG,subset3_AT,subset3_CA,subset3_CG,subset3_CT,subset3_TA,subset3_TG,subset3_TC)
subset3_MA=subset3_MA[subset3_MA$TUMOR_X10>0.08,]

subset3_MA$STATUS="SUBSET3"

#################### 5) chi-squre #########################
group=rbind(subset1.2,subset1.3, subset3_MA)
snp_file_subset=snp_ids[!(snp_ids$V3 %in% unique(group$V3)), ]
snp_file_subset=snp_file_subset[snp_file_subset$V4!=snp_file_subset$V5,]

snp_file_subset=snp_file_subset[which(as.numeric(as.character(snp_file_subset$NORMAL_DP))>9 & as.numeric(as.character(snp_file_subset$TUMOR_DP))>9),]

## Remove Multi Alleles
snp_file_subset_GA=snp_file_subset[snp_file_subset$V4=="G" & snp_file_subset$V5=="A" & snp_file_subset$TUMOR_X2>4 & snp_file_subset$TUMOR_X6>4 & (snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_GC=snp_file_subset[snp_file_subset$V4=="G" & snp_file_subset$V5=="C" & snp_file_subset$TUMOR_X3>4 & snp_file_subset$TUMOR_X7>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_GT=snp_file_subset[snp_file_subset$V4=="G" & snp_file_subset$V5=="T" & snp_file_subset$TUMOR_X5>4 & snp_file_subset$TUMOR_X9>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7)<1,]
snp_file_subset_AC=snp_file_subset[snp_file_subset$V4=="A" & snp_file_subset$V5=="C" & snp_file_subset$TUMOR_X3>4 & snp_file_subset$TUMOR_X7>4 & (snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_AG=snp_file_subset[snp_file_subset$V4=="A" & snp_file_subset$V5=="G" & snp_file_subset$TUMOR_X4>4 & snp_file_subset$TUMOR_X8>4 & (snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_AT=snp_file_subset[snp_file_subset$V4=="A" & snp_file_subset$V5=="T" & snp_file_subset$TUMOR_X5>4 & snp_file_subset$TUMOR_X9>4 & (snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7+snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8)<1,]
snp_file_subset_CA=snp_file_subset[snp_file_subset$V4=="C" & snp_file_subset$V5=="A" & snp_file_subset$TUMOR_X2>4 & snp_file_subset$TUMOR_X6>4 & (snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_CG=snp_file_subset[snp_file_subset$V4=="C" & snp_file_subset$V5=="G" & snp_file_subset$TUMOR_X4>4 & snp_file_subset$TUMOR_X8>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X5+snp_file_subset$TUMOR_X9)<1,]
snp_file_subset_CT=snp_file_subset[snp_file_subset$V4=="C" & snp_file_subset$V5=="T" & snp_file_subset$TUMOR_X5>4 & snp_file_subset$TUMOR_X9>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8)<1,]
snp_file_subset_TA=snp_file_subset[snp_file_subset$V4=="T" & snp_file_subset$V5=="A" & snp_file_subset$TUMOR_X2>4 & snp_file_subset$TUMOR_X6>4 & (snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7+snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8)<1,]
snp_file_subset_TG=snp_file_subset[snp_file_subset$V4=="T" & snp_file_subset$V5=="G" & snp_file_subset$TUMOR_X4>4 & snp_file_subset$TUMOR_X8>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X3+snp_file_subset$TUMOR_X7)<1,]
snp_file_subset_TC=snp_file_subset[snp_file_subset$V4=="T" & snp_file_subset$V5=="C" & snp_file_subset$TUMOR_X3>4 & snp_file_subset$TUMOR_X7>4 & (snp_file_subset$TUMOR_X2+snp_file_subset$TUMOR_X6+snp_file_subset$TUMOR_X4+snp_file_subset$TUMOR_X8)<1,]

snp_file_subset_MA=rbind(snp_file_subset_GA,snp_file_subset_GC,snp_file_subset_GT,snp_file_subset_AC,snp_file_subset_AG,snp_file_subset_AT,snp_file_subset_CA,snp_file_subset_CG,snp_file_subset_CT,snp_file_subset_TA,snp_file_subset_TG,snp_file_subset_TC)
snp_file_subset_MA$N_REF=snp_file_subset_MA$NORMAL_DP-(snp_file_subset_MA$NORMAL_X10*snp_file_subset_MA$NORMAL_DP)
snp_file_subset_MA$N_ALT=(snp_file_subset_MA$NORMAL_X10*snp_file_subset_MA$NORMAL_DP)

snp_file_subset_MA$T_REF=snp_file_subset_MA$TUMOR_DP-(snp_file_subset_MA$TUMOR_X10*snp_file_subset_MA$TUMOR_DP)
snp_file_subset_MA$T_ALT=(snp_file_subset_MA$TUMOR_X10*snp_file_subset_MA$TUMOR_DP)

#file_p <- adply(snp_file_subset_MA, 1, function(x) {chisq.test(matrix(c(x$T_ALT, x$N_ALT, x$T_REF, x$N_REF),nrow=2))$p.value})

################## 6) ########################################
group=rbind(subset1.2,subset2, subset3)
snp_file_subset=snp_ids[!(snp_ids$V3 %in% unique(group$V3)), ]


#################### 7) ####################################
pat<-'ASRD=[0-9,\\.]+'

somatic$ASRD=0
somatic$ASRD=gsub("ASRD=","",str_extract(somatic$V8, pat))
somatic$ASRD=as.numeric(as.character(somatic$ASRD))
somatic_pass=somatic[somatic$ASRD>0.9 & somatic$V7=="PASS",]
somatic_fail=somatic[!(somatic$V3 %in% unique(somatic_pass$V3)), ]


somatic_pass=somatic_pass[,1:39]
somatic_fail=somatic_fail[,1:39]
somatic_pass$STATUS="SOMATIC_PASS"
somatic_fail$STATUS="SOMATIC_FAIL"




################### PLOT #####################
plot=rbind(subset1.2,subset1.3,subset3_MA,somatic_pass,somatic_fail)
ggplot(plot,aes(MP,fill=STATUS))+geom_density(linetype="blank",alpha=0.4)+coord_cartesian(ylim=c(0,2.5e+02), xlim=c(0,0.07))+theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  theme( legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())

ggplot(plot,aes(MP,fill=STATUS))+geom_density(alpha=0.7)+theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  theme( legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
         panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
         panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())


