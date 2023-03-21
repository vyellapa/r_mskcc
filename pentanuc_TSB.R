#http://htmlpreview.github.io/?https://github.com/kgori/sigfit/blob/master/inst/doc/sigfit_vignette.html
library(dplyr)
library(BSgenome)
library(Biostrings)
library("seqRFLP")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg19)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

path_input = "/Users/yellapav/Desktop/tsb/"
setwd(path_input)


#Load TxDB
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Read mutations and create pentanucleotide sequence
muts = read.table("cave_fm6_chap_alex_pd26414_5col.txt",sep="\t")
muts$context=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr",muts$V2),muts$V3-2,muts$V3+2))
muts$revcomp=as.character(unlist(lapply(muts$context,revComp)))
muts$RREF=as.character(unlist(lapply(muts$V4,revComp)))
muts$RALT=as.character(unlist(lapply(muts$V5,revComp)))

muts$base5S=paste0(substr(muts$context,1,2),"[",muts$V4,">",muts$V5,"]",substr(muts$context,4,5))
muts$base5R=paste0(substr(muts$revcomp,1,2),"[",muts$RREF,">",muts$RALT,"]",substr(muts$revcomp,4,5))
muts$base5="A"
muts = muts %>% mutate(base5 = ifelse(V4=="G"| V5=="A",base5R,base5S))


#In the temp directory create VCFs
dir.create("./temp")
temp_path=sprintf("%stemp/",path_input)
setwd(path_input)


#MutationalPatterns requires vcf; Create VCFs
vcf_header="##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
write.table(vcf_header,file="./temp/vcf_header", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE,quote=F)
for(i in unique(muts$V1)) {
  sub = muts %>% dplyr::filter(V1 == i) %>% dplyr::mutate(ID=base5,QUAL=".",FILTER=".",INFO=".") %>% dplyr::select(V2,V3,ID,V4,V5,QUAL,FILTER,INFO)
  write.table(vcf_header,file=sprintf("./temp/%s.vcf",i), append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE,quote=F)
  write.table(sub,file=sprintf("./temp/%s.vcf",i), append=TRUE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE,quote=F)
}
#Read vcfs and get trancribed vs untranscribed status
vcf_files = paste0(temp_path,list.files(path = temp_path, pattern = "*.vcf"))
sample_names = gsub(".vcf","",list.files(path = temp_path, pattern = "*.vcf"))
vcfs = read_vcfs_as_granges(vcf_files, sample_names, ref_genome)


strand = mut_strand(vcfs[[1]], genes_hg19)
TSB = data.frame(cbind(data.frame(vcfs[[1]]),data.frame(strand)))
TSB.bkup=TSB

for(i in 2:length(vcfs)) {
  strand = mut_strand(vcfs[[i]], genes_hg19)
  TSB = rbind(TSB, data.frame(cbind(data.frame(vcfs[[i]]),data.frame(strand))))
}

TSB$ALT1=as.character(unlist(lapply(TSB$ALT,function(x) {return(unlist(as.list(as.character(x))))})))
TSB = data.frame(TSB) %>% dplyr::select(seqnames, start, REF, ALT1, strand.1) %>% dplyr::mutate(key=paste(seqnames, start,sep=":"))
#write.table(TSB,file="TSB_results.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote=F)


muts$key = paste(muts$V6,muts$V7,muts$V9,muts$V10,sep=":")
TSB$key= paste(gsub("chr","",TSB$seqnames),TSB$start,TSB$REF,TSB$ALT1,sep=":")
muts.tsb=left_join(muts,TSB,by="key") %>% unique()
muts.tsb = muts.tsb %>% dplyr::select(V2,V6,V7,V9,V10,base5,strand.1)
write.table(muts.tsb,file="muts.tsb_results.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote=F)

mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
mut_mat_df = as.data.frame(mut_mat_s) %>% dplyr::mutate(type = rownames(as.data.frame(mut_mat_s)))

mut_df_melt = melt(mut_mat_df)
s=as.data.frame(matrix(unlist(strsplit(mut_df_melt$type,"-")),ncol=2,byrow=T))
mut_df_melt = cbind(s,mut_df_melt)
mut_df_melt$mut=substr(mut_df_melt$V1,3,5)
colnames(mut_df_melt) = c("trinucl","effect","type","sample","value","mut")


file=(dcast(mut_df_melt,type~sample, value=value))
rownames(file) = gsub("-"," ",file$type) 
file = file %>% dplyr::select(-type)

col_file<-unique(paste(mut_df_melt$trinucl, z$effect))
mm_col <- c(rep("lightskyblue",32), rep("black",32), rep("firebrick2",32),rep("gray88",32),rep("darkolivegreen3",32),
            rep("lightpink1",32))
file2 = file[col_file,]

pdf("strandbias_all_DP.pdf", width = 15, height = 7)
for(i in (1:ncol(file2))){
  
  par(mar=c(5,5,2,2), xpd=T)
  x<- barplot(as.numeric(file2[,i]), col=c("blue4","brown3"), main ="", cex.names=0.7, names.arg = "", las=2, border = NA,
              ylim=c(0, max(as.numeric(file2[,i]))+10), yaxt="n")
  
  x<- format(x, scientific = FALSE)
  
  options(digits=2)
  seque<- seq(1, 191, by=2)
  loc<- NULL
  for(w in (1:length(seque)))
  {
    loc<- c(loc, (as.numeric(x[seque[w],1])+as.numeric(x[seque[w]+1,1]))/2)
  }
  mtext(side=1, at=loc, text = unique(all_sel$trinucl), las=2, line = 1, cex=0.8)
  
  mycol <- t_col("lightskyblue", perc = 50, name = "lightskyblue")
  rect(0, max(as.numeric(file2[,i]))+10, x[32.75], 0, col = mycol, border = NA)
  
  mycol <- t_col("black", perc = 50, name = "black")
  rect(x[32.75], max(as.numeric(file2[,i]))+10, x[64.75], 0, col = mycol, border = NA)
  
  mycol <- t_col("firebrick2", perc = 50, name = "firebrick2")
  rect(x[64.75], max(as.numeric(file2[,i]))+10, x[96.75], 0, col = mycol, border = NA)
  
  mycol <- t_col("gray88", perc = 50, name = "gray88")
  rect(x[96.75], max(as.numeric(file2[,i]))+10, x[128.75], 0, col = mycol, border = NA)
  
  mycol <- t_col("darkolivegreen3", perc = 50, name = "darkolivegreen3")
  rect(x[128.75], max(as.numeric(file2[,i]))+10, x[160.75], 0, col = mycol, border = NA)
  
  mycol <- t_col("lightpink1", perc = 50, name = "lightpink1")
  rect(x[160.75], max(as.numeric(file2[,i]))+10, x[192], 0, col = mycol, border = NA)
  par(new=T)
  x<- barplot(as.numeric(file2[,i]), col=c("blue4","brown3"), main =paste(colnames(file2)[i], sum(file2[,i]), sep=" - Number of SNVs="), cex.names=0.7,
              names.arg = "", las=2, border = NA,
              ylim=c(0, max(as.numeric(file2[,i]))+10), cex.axis=1.5)
  
  legend("topright",legend=c("Transcribed Strand","Untranscribed Strand"), pch=15,border="n",
         col=c("blue4","brown3"),box.lwd = 0, box.col = "white",bg = "white",
         cex=1, pt.cex=1, inset=c(0.035,0.0),x.intersp = 1,y.intersp = 1)
  
  
}

dev.off()




# Fishers test between mm1 & non-mm1
mut_df_melt$code_mm1 = NA
mut_df_melt$code_mm1[mut_df_melt$trinucl  %in% c("A[C>T]C","C[C>T]A","G[C>T]A","G[C>T]T","G[C>T]C")]<-"mm1"
mut_df_melt$code_mm1[! mut_df_melt$trinucl  %in% c("A[C>T]C","C[C>T]A","G[C>T]A","G[C>T]T","G[C>T]C")]<-"no_mm1"
sample_list<- unique(mut_df_melt$sample)
alfa<- list() 

for(i in (1:length(sample_list)))
{
  samz<- mut_df_melt[mut_df_melt$sample == sample_list[i],]
  if(nrow(as.matrix(table(samz$code_mm1, samz$effect)))>1 &
     ncol(as.matrix(table(samz$code_mm1, samz$effect)))>1){
    
    fishers=c(as.data.frame((samz %>% group_by(code_mm1,effect) %>% summarize(sum =sum(value)))) %>% dplyr::select(sum))
    alfa[[i]]<-c(as.character(sample_list[i]), as.vector(fishers$sum), fisher.test(matrix(as.vector(fishers$sum),byrow=T,ncol=2))$p.value)
  }
}
alfa2<- do.call("rbind", alfa)
alfa2<- as.data.frame(alfa2)
colnames(alfa2)<-c("sample","T_MM1","U_MM1","T","U","pvalue")
alfa2$pvalue<- as.numeric(as.character(alfa2$pvalue))
alfa2$q<- p.adjust(alfa2$pvalue, method = "fdr")
alfa2[alfa2$pvalue<0.05,]
alfa2[alfa2$q<0.1,]



