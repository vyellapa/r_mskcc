#splits 2 vcf files by clusters based on dirichlet cluster file
#Rscript ~/local/msk_scripts/cluster2VCF.R E-H-109099-T1-1-D1-1.caveman.filter.vcf E-H-109099-T2-1-D1-1.caveman.filter.vcf ../dirchlet/resultsXF/E-H-109099-T1-1-D1-1_DPoutput_10000iters_1000burnin/E-H-109099-T1-1-D1-1_bestConsensusAssignments.bed `pwd`/signature/ > signature/E99_cluster_signature.snv.input


library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)

args <- commandArgs(trailingOnly = TRUE)
vfile1=args[1]
vfile2=args[2]
clusterFile=args[3]
outDir=args[4]

cat("sample	vcf_file\n")
vcf1=readVcf(vfile1, "hg19")
vcf2=readVcf(vfile2, "hg19")
#direc=dirname(vfile1)
name=gsub(".caveman.filter.vcf","",basename(vfile1))

vcf=rbind(vcf1,vcf2)
#seqlevels(vcf) <- paste("chr", seqlevels(vcf), sep="")
cluster=read.table(clusterFile,header=T,sep="\t")
cluster=cluster[!is.na(cluster$cluster),]

#for(i in unique(cluster$cluster)){
for(i in 1:7){
clusterSub=cluster[cluster$cluster==i,]
g1=GRanges(seqnames=Rle(as.character(clusterSub$chr)), ranges = IRanges(start = clusterSub$start, end=clusterSub$end))
#seqlevels(g1) <- paste("chr", seqlevels(g1), sep="")
#hits=findOverlaps(vcf,g1)
#new=ranges(vcf)[subjectHits(hits)]
rd = rowRanges(vcf)
starts = start(ranges(rd))
new=vcf[paste(seqnames(rd),starts) %in% paste(seqnames(g1),end(ranges(g1))),]
filename=paste(outDir,"/",name,"_cluster",i,".vcf",sep="")
writeVcf(new, file=filename)

sName=paste(name,"_cluster",i,sep="")
cat(paste(sName,filename,sep="\t"))
cat("\n")
}
