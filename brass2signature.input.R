
suppressMessages(library(optparse, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))


option_list = list(make_option(c("-i", "--input"), type = "character",
                               help = "Input file.This is the brass annot.bedpe",
                               metavar = "path", action = "store"),
                   make_option(c("-o", "--output"), type = "character",
                               help = "Output file", metavar = "path"));



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



library("BSgenome.Hsapiens.UCSC.hg19")

a=read.table(opt$i, sep="\t",header=F)
colnames(a)=c("chr1","start1","end1","chr2","start2","end2","id/name","brass_score","strand1","strand2","sample","svclass","bkdist","assembly_score","readpair names","readpair count","bal_trans","inv","occL","occH","copynumber_flag","range_blat","Brass Notation","non-template","micro-homology","assembled readnames","assembled read count","gene1","gene_id1","transcript_id1","strand1","end_phase1","region1","region_number1","total_region_count1","first/last1","gene2","gene_id2","transcript_id2","strand2","phase2","region2","region_number2","total_region_count2","first/last2","fusion_flag")

a$bp1_seq="NN"
a$bp1_index="NN"
a$bp2_seq="NN"
a$bp2_index="NN"

for(i in 1:nrow(a)) {
  bp1_chrom=paste("chr",a[i,"chr1"], sep="")
  bp1_start=a[i,"start1"]
  #bp1_seq_post=getSeq(Hsapiens,bp1_chrom,bp1_start,bp1_start+15)
  #bp1_seq_pre=getSeq(Hsapiens,bp1_chrom,bp1_start-15,bp1_start)
  bp1_seq=as.character(getSeq(Hsapiens,bp1_chrom,bp1_start-15,bp1_start+15))
  bp1_index=paste(a[i,"chr1"], a[i,"start1"], sep="_")
  
  a[i,c("bp1_seq")]=as.character(bp1_seq)
  a[i,]$bp1_index=bp1_index
  
  bp2_chrom=paste("chr",a[i,"chr2"], sep="")
  bp2_start=a[i,"start2"]
  #bp2_seq_post=getSeq(Hsapiens,bp2_chrom,bp2_start,bp2_start+15)
  #bp2_seq_pre=getSeq(Hsapiens,bp2_chrom,bp2_start-15,bp2_start)
  bp2_seq=as.character(getSeq(Hsapiens,bp2_chrom,bp2_start-15,bp2_start+15))
  bp2_index=paste(a[i,"chr2"], a[i,"start2"], sep="_")
  
  a[i,]$bp2_seq=bp2_seq
  a[i,]$bp2_index=bp2_index
}
  cat(c("Total number of SVs process is ",nrow(a),"\n"))

bp1_lengths <- data.frame(ind=character(), 
                 next_r=numeric(), 
                 stringsAsFactors=FALSE) 

bp2_lengths <- data.frame(ind=character(), 
                          next_r=numeric(), 
                          stringsAsFactors=FALSE) 

for(j in unique(a$chr1)){
b=sort(a[which(a$chr1==j),c(2)])
if(length(b) > 1) {
  for(i in 1:(length(b)-1)) {
    len_diff=b[i+1]-b[i]
    zzz=paste(j,b[i],sep="_")
    bp1_lengths <- rbind(bp1_lengths,data.frame(ind=zzz,next_r=len_diff))
  #print(c(zz,dist_bp1))
}}
}

for(j in unique(a$chr2)){
  b=sort(a[which(a$chr2==j),c(5)])
  if(length(b) > 1) {
    for(i in 1:(length(b)-1)) {
      len_diff=b[i+1]-b[i]
      zzz=paste(j,b[i],sep="_")
      bp2_lengths <- rbind(bp2_lengths,data.frame(ind=zzz,next_r=len_diff))
      #print(c(zz,dist_bp1))
    }}
}

merged=merge(a,bp1_lengths,by.x = "bp1_index", by.y = "ind", all.x = T)
merged=merge(merged,bp2_lengths,by.x = "bp2_index", by.y = "ind", all.x = T)


write.table(merged,file=opt$o, append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
