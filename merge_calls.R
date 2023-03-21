#OLa request
library(dplyr)
library(reshape2)
library(tidyr)

MRN = read.table("~/Desktop/ola_request/MRN.tsv", sep="\t", header=F)
colnames(MRN) = c("MRN","Name","LEUKID")

nih = read.table("~/Desktop/ola_request/ola_request_nih_ids.txt", sep="\t", header=T)
nih=nih[,c(1,6)]


MRN = MRN %>% left_join(nih, by=c('LEUKID' = 'Individual.ID')) %>% 
  dplyr::mutate(LEUKID = ifelse(is.na(Experiment.System.ID), LEUKID,as.character(Experiment.System.ID))) %>% 
  dplyr::select(-Experiment.System.ID) %>% unique()



mutations = read.table("~/Desktop/ola_request/p222_mutations.csv", sep=",", header=T)
mutations = mutations[,c(1,2,5,6,7,20,21,22,24,25)]

p12.snv = read.table("~/Desktop/ola_request/p212_snv.txt", sep="\t", header=T) %>% 
  filter(MANUAL_ANNOTATION=="LIKELY" | MANUAL_ANNOTATION=="ONCOGENIC") %>% 
  dplyr::select(TARGET_NAME,MANUAL_ANNOTATION,GENE, PROTEIN_CHANGE,TARGET_VAF,EFFECT,CHR,START,REF,ALT)
p12.indel = read.table("~/Desktop/ola_request/p212_indel.txt", sep="\t", header=T) %>% 
  filter(MANUAL_ANNOTATION=="LIKELY" | MANUAL_ANNOTATION=="ONCOGENIC") %>%
  dplyr::select(TARGET_NAME,MANUAL_ANNOTATION,GENE, PROTEIN_CHANGE,TARGET_VAF,EFFECT, CHR, START, REF,ALT) 
p12.mut = rbind(p12.snv,p12.indel) %>% dplyr::rename(VAG_GENE=GENE,VAG_PROTEIN_CHANGE=PROTEIN_CHANGE, LEUKID=TARGET_NAME,VAG_EFFECT=EFFECT)

mutations = rbind(mutations, p12.mut) %>% unique()

#############
cnv = read.table("~/Desktop/ola_request/p222_cnv.csv", sep=",", header=T)
cnv = (cnv) %>% dplyr::select(LEUKID,Cytogenetic_abnormality,SEQ_result) %>% unique() %>% mutate(SEQ_result=ifelse(SEQ_result==1,"YES","NO")) #%>% spread(Cytogenetic_abnormality,SEQ_result)

p12.cna = read.table("~/Desktop/ola_request/p212_cna.txt", sep="\t", header=T) %>% dplyr::select(Sample,del_1p,del_17p,del_16q,del_13q,del_12p,amp_1q)
colnames(p12.cna)=c("Sample","Del(1p)","Del(17p)","Del(16q)","Del(13q)","Del(12p)","Gain(1q)")
p12.c = p12.cna  %>% melt() %>% unique() %>% mutate(value=ifelse(value==1,"YES","NO")) %>% dplyr::rename(Cytogenetic_abnormality=variable,SEQ_result=value)





#################
p12.brass = read.table("~/Desktop/ola_request/p212_brass.txt", sep="\t", header=T) %>% dplyr::select(CHR1,CHR2,SAMPLE) %>% dplyr::filter(CHR1==14 | CHR2==14)
p12.delly = read.table("~/Desktop/ola_request/p212_delly.txt", sep="\t", header=T) %>% dplyr::select(chr1,chr2,sample)
colnames(p12.delly) = colnames(p12.brass)
p12.delly$SAMPLE = gsub("_50","",p12.delly$SAMPLE)
p12.trans = rbind(p12.delly, p12.brass) %>% dplyr::mutate(t= ifelse(as.numeric(as.character(CHR1))<as.numeric(as.character(CHR2)), paste0("t(",CHR1,";",CHR2,")"),paste0("t(",CHR2,";",CHR1,")") ))
p12.trans$SEQ_result="YES"
p12.t = p12.trans %>% dplyr::select(SAMPLE, t,SEQ_result ) %>% dplyr::rename(Sample = SAMPLE, Cytogenetic_abnormality=t ) %>% 
  dplyr::filter( Cytogenetic_abnormality %in% c("t(4;14)","t(6;14)","t(11;14)","t(14;16)","t(14;20)")) %>% unique()

p12.cnv = rbind(p12.c, p12.t) %>% unique() %>% dplyr::rename(LEUKID=Sample)
all.cnv = rbind(cnv, p12.cnv) %>% spread(Cytogenetic_abnormality,SEQ_result)


RESULTS = (MRN) %>% left_join(mutations, by = c("LEUKID" = "LEUKID")) %>% left_join(all.cnv,by = c("LEUKID" = "LEUKID")) 
write.table(RESULTS,"~/Desktop/ola_request/results_p212_p222.txt",sep="\t",eol="\n",col.names = T)









#colnames(p12.delly) = colnames(p12.brass)
#p12.delly$SAMPLE = gsub("_50","",p12.delly$SAMPLE)
#p12.trans = rbind(p12.delly, p12.brass) %>% dplyr::mutate(t= ifelse(as.numeric(as.character(CHR1))<as.numeric(as.character(CHR2)), paste0("t(",CHR1,";",CHR2,")"),paste0("t(",CHR2,";",CHR1,")") ))
#p12.trans$SEQ_result="YES"
#p12.t = p12.trans %>% dplyr::select(SAMPLE, t,SEQ_result ) %>% dplyr::rename(Sample = SAMPLE, Cytogenetic_abnormality=t ) %>% 
#  dplyr::filter( Cytogenetic_abnormality %in% c("t(4;14)","t(6;14)","t(11;14)","t(14;16)","t(14;20)")) %>% unique()

#p12.cnv = rbind(p12.c, p12.t) %>% unique() %>% dplyr::rename(LEUKID=Sample)
#all.cnv = rbind(cnv, p12.cnv) %>% spread(Cytogenetic_abnormality,SEQ_result)






