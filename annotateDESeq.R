library("biomaRt")
setwd("/Users/yellapav/Desktop/MM_CD99")
listMarts(host="www.ensembl.org")
chr=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

all75=getBM(attributes=c('hgnc_symbol','ensembl_gene_id','strand'), mart = ensembl_75)
all75c=getBM(attributes=c('chromosome_name','external_gene_id','uniprot_genename','wikigene_name','hgnc_symbol','ensembl_gene_id','strand'), mart = ensembl_75)
all75c=all75c[all75c$chromosome_name %in% chr, ]


for(i in list.files(pattern = "*.tsv")) {
  out=gsub(".tsv", "_annotated.tsv", i)
  deseqOut=read.table(i,header=T,sep="\t")
  
  merged=merge(deseqOut,all75, by.x="id", by.y="ensembl_gene_id", all.x=TRUE)
  merged=merged[,-c(2)]
  write.table(merged,file=out, append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
  }