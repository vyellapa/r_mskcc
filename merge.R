setwd("/Users/yellapav/Desktop/record/images/metadata")
svs = read.table("luke_mrn.txt",sep = "\t", header = T)
pur = read.table("final_purity.txt",sep = "\t", header = T) %>% dplyr::select(-c(X))
p292 = read.table("p292_MRN_pathid_leukgenID.txt",sep = "\t", header = T)
p292 = p292 %>% dplyr::select(-c(Has.Bam,Has.Fastq,Individual.ID,Extraction.ID,Workflow.ID,Primary.Key,Individual.Leukid,Specimen.Leukid)) %>% distinct()

p292 = p292 %>% left_join(pur, by = c("Workflow.Leukid" = "sample"))

svs_p292 = svs %>% dplyr::left_join(p292, by = c("mrn" = "MRN"))
write.table(svs_p292,"image_metadata.txt",append = F,quote = F,sep = "\t",row.names = F, col.names = T)



library("readxl")
library(tidyr)

sv = read_excel("/Users/yellapav/Desktop//record/p292/uploads/p292_supplementary.xlsx", sheet = "FISH.Comparison")
cnv = read_excel("/Users/yellapav/Desktop//record/p292/uploads/p292_supplementary.xlsx", sheet = "Array.Comparison")
snv = read_excel("/Users/yellapav/Desktop//record/p292/uploads/p292_supplementary.xlsx", sheet = "myTYPE.SNVs")
indel = read_excel("/Users/yellapav/Desktop//record/p292/uploads/p292_supplementary.xlsx", sheet = "myTYPE.Indels")


sv = sv %>% dplyr::filter(called.by %in% unique(grep("myTYPE",sv$called.by,value = T))) %>% 
  dplyr::mutate(TARGET_NAME = Sample,VAG_GENE = translocation) %>% 
  dplyr::select(TARGET_NAME, VAG_GENE) %>% distinct()

cnv = cnv %>% dplyr::filter(called.by %in% unique(grep("myTYPE",cnv$called.by,value = T))) %>%
  dplyr::mutate(VAG_GENE=paste(arm,status,sep="_")) %>% 
  dplyr::filter(VAG_GENE %in% c("11q_gain","13q_loss","17p_loss","1q_gain","1p_loss")) %>% 
  dplyr::select(ID,VAG_GENE) %>% dplyr::rename(TARGET_NAME = ID) %>% distinct()


snv = snv %>% dplyr::filter(Annotation %in% c("ONCOGENIC","LIKELY")) %>% 
  dplyr::filter(VAG_GENE %in% c("BRAF","FAM46C","KRAS","NRAS","TP53","DIS3")) %>% 
  dplyr::select(TARGET_NAME,VAG_GENE)


indel = indel %>% dplyr::filter(Annotation %in% c("ONCOGENIC","LIKELY")) %>% 
  dplyr::filter(VAG_GENE %in% c("BRAF","FAM46C","KRAS","NRAS","TP53","DIS3")) %>% 
  dplyr::select(TARGET_NAME,VAG_GENE)

muts = rbind(snv,indel)
muts = rbind(muts, cnv)
muts = rbind(muts, sv)
muts = distinct(muts)


muts$value=1
muts1 = muts %>% tidyr::spread(VAG_GENE,value) %>% as.data.frame()
muts1[is.na(muts1)] = 0

muts1 = muts1[,c("TARGET_NAME","11q_gain","13q_loss","17p_loss","1p_loss","1q_gain","BRAF","DIS3","FAM46C","KRAS","NRAS","TP53","t(11;14)","t(14;16)","t(14;20)","t(4;14)","t(6;14)","t(8;14)")]
write.table(muts1,"mutation_matrix.txt",append = F,quote = F,sep = "\t",row.names = F, col.names = T)

metadata = read.table("image_metadata.txt",sep = "\t",header = T)
meta = dplyr::left_join(svs_p292, muts1, by = c( "Workflow.Leukid" = "TARGET_NAME"))
write.table(meta,"merged_mutation_matrix.txt",append = F,quote = F,sep = "\t",row.names = F, col.names = T)
