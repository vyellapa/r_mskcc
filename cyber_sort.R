library(tximport)
library(tximportData)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(dplyr)


txdb = makeTxDbFromGFF("~/local/resources/gencode.v19.annotation.gtf",format = c("gtf"))
k <- keys(txdb, keytype = "TXNAME")
keytypes(txdb)
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files = list.files("/Users/yellapav/Desktop/record/BMS/salmon/", pattern = "1[.]quant.sf",full.names = T)
names(files)=substr(basename(files),1,15)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$counts)
TPM = (txi$abundance)
write.table(TPM, "~/Desktop/record/BMS/TPM_matrix.txt", quote = F, sep = "\t",eol = "\n")


h = read.table("/Users/yellapav/Desktop/record/BMS/salmon//IID_H196503_N01_01_WT01_1d20ba0f2cf4d748_EA27-301978_1_quant.sf",sep = "\t", header = T)

ensg = sapply((strsplit(as.character(h$Name), split = "[|]")), "[[", 2) %>% as.data.frame
gname = sapply((strsplit(as.character(h$Name), split = "[|]")), "[[", 5) %>% as.data.frame

anno.df = cbind(ensg, gname) 
colnames(anno.df) = c("ENSG","Gene")
anno.df = distinct(anno.df)
anno.df$ENSG = as.character(anno.df$ENSG)


TPM.df = as.data.frame(rbind((TPM)))

TPM.df$ENSG = rownames(TPM.df)
TPM.df = left_join(TPM.df,anno.df, by = "ENSG") 

TPM.df = TPM.df[!duplicated(TPM.df$Gene),]
rownames(TPM.df) = TPM.df$Gene
TPM.df = TPM.df %>% dplyr::select(-c("ENSG","Gene"))


write.table(TPM.df, "~/Desktop/record/BMS/TPM_anno_matrix.txt", quote = F, sep = "\t",eol = "\n")

TPM.df.see = TPM.df[apply(TPM.df,1,max)>2,]
write.table(TPM.df.see, "~/Desktop/record/BMS/TPM_anno_matrix_fpkm2.txt", quote = F, sep = "\t",eol = "\n")
############
cyber = read.table("~/Desktop/record/BMS/cybersort.txt",sep="\t",header = T)






col.ann = as.data.frame(cyber$Mixture)
colnames(col.ann) = c("sample")
col.ann = (col.ann) %>% dplyr::mutate(type= ifelse(sample %in% grep("_T0",col.ann$sample,value=T), "TUMOR","NORMAL"))
types = col.ann$type;color.annot=c("#7D093B","#60B96C","#D31996","#E9E9E9"); names(color.annot) = c("NORMAL","TUMOR","DEL","NULL")
ha = HeatmapAnnotation(df = data.frame(type = types), 
                       FAM46C = data.frame(FAM46C = as.vector(colan$status)), 
                       Translocation = data.frame(Translocation = as.vector(colan$status)),
                       col = list(type = c(color.annot[1:2]), 
                                  FAM46C = c("DEL" = "#D31996","NULL" = "#E9E9E9","MUT" = "#52A20D", "DEL_MUT" = "#3F3451") ) )
ha = HeatmapAnnotation(type = types, 
                       FAM46C = colan$status, 
                       Translocation = colan.t$status,
                       col = list(type = c(color.annot[1:2]), 
                                  FAM46C = c("DEL" = "#D31996","NULL" = "#E9E9E9","MUT" = "#52A20D", "DEL_MUT" = "#3F3451"),
                                  Translocation = c("CCND1" = "#2FAD3B","NULL" = "#E9E9E9","CCND3" = "#CF0535", 
                                                    "HRD" = "#3D1C00", "MAF" = "#FEA014", "MAFB" = "#511A8E", "MMSET" = "#0580E7", "MYC" = "#FA2A00")) )

ha = HeatmapAnnotation(type = types, 
                       col = list(type = c(color.annot[1:2])) )

head(cyber)
rownames(cyber) = as.character(cyber$Mixture)
cyber = cyber %>% dplyr::select(-c(Mixture,P.value,Correlation,RMSE,Absolute.score..sig.score.))

pdf("~/Desktop/hm.pdf",width = 16, height = 10)
Heatmap(t(as.matrix(cyber)), top_annotation = ha, clustering_method_columns = "ward.D",clustering_distance_columns =  "euclidean", column_km = 4)
dev.off()










txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files <- file.path("Users/yellapav/Desktop/record/BMS/salmon/", "salmon", samples$run, "quant.sf")

txi <- tximport("/Users/yellapav/Desktop/record/BMS/salmon/IID_H196349_N02_01_WT01_1d20ba0f2cf4d748_EA27-301937_1_quant.sf", type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
names(txi)

library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)




check = function(j) {
  foreach(p=unique(data$marker), .combine=rbind) %dopar% {
  #  for(j in unique(data$trt_group)) {
  i=as.character(as.vector(unlist(p)))
  sub=data %>% dplyr::filter(marker==i & trt_group!=j & (sample_time=="DAY8" |sample_time=="DAY1" )) %>% summarise(p.vl=t.test(digital_count~sample_time)$p.value) %>% dplyr::mutate(marker=i, treatment=j) %>% as.data.frame() 
  #results_ttest[[counter]]=sub
  #counter=counter+1
  return(sub)
  
}
}

