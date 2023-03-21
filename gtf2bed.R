## Usage 
#Rscript --vanilla gtf2bed.R <input.gtf> <out.bed>
args=commandArgs(TRUE)
s1=args[1] 
s2=args[2]
unlink(s2)
#gtf=read.table("~/Desktop/dtop/record/splice/context/context/Homo_sapiens.GRCh37.75.gtf",header=F,sep="\t")
gtf=read.table(s1,header=F,sep="\t")
gtfe=gtf[which(gtf$V3=="exon"),]
library("stringr")
pat<-"gene_id ENSG[0-9]+"
pat2<-"gene_name [A-Z0-9a-z-.]+"
uniqG=unique(gsub("gene_id ","",str_extract(gtfe$V9, pat)))
gtfe$gene_id=gsub("gene_id ","",str_extract(gtfe$V9, pat))
gtfe$gene_name=gsub("gene_name ","",str_extract(gtfe$V9, pat2))

for(i in uniqG) {
  sub=gtfe[which(gtfe$gene_id==i),]
  mn=min(sub$V4)
  mx=max(sub$V5)
  chr=unique(sub$V1)
  gid=unique(sub$gene_id)
  gname=unique(sub$gene_name)
  ##v22=unique(sub$V2) ## These are not uniq and causing redundancy
  v77=unique(sub$V7)
  p=paste(v77, gid,gname,sep=":")
  r=paste(chr,mn,mx,p,sep="\t")
  #print(r)
  write.table(r,file=s2, append=TRUE, quote=F,sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
}
