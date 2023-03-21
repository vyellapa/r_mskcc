library(gsheet)
library(dplyr)
library(reshape2)
library(ggplot2)
names=read.table("/Users/yellapav/Desktop/record/p292/misc/p292_sample_sheet.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
names=names[,c("Workflow.Leukid","Aliquot.ID")]


exclude = c("I-H-135277-T1-1-D1-1","I-H-135351-T1-1-D1-1","I-H-135304-T1-1-D1-1")
clinical <- read.csv(construct_download_url("https://docs.google.com/spreadsheets/d/1vOO7ihNDfv4cElTLJT18rmpjZeJrHzjMixXkAuIgza4", format = "csv", sheetid = "1075620643"), stringsAsFactors=FALSE)
clinical = clinical %>% dplyr::filter(!(CATEG %in% c("WALDENSTROM","SMOLDERING","AMYLOIDOSIS","MGUS","PLASMA.CELL.LEUKEMIA"))) %>% 
  left_join(names,by = c('Pathology.ID' = 'Aliquot.ID')) %>% 
  dplyr::filter(!(Workflow.Leukid %in% exclude))



plot = read.table("/Users/yellapav/Desktop/record/p292/misc/fish.arm.results.plot.txt",sep="\t",header=T,stringsAsFactors = F) %>%
  dplyr::filter(translocation!="11q.gain" & translocation!="12p.loss" & translocation!="BIRC3") %>%
  dplyr::filter(translocation!="Sequenced") %>% mutate(translocation = ifelse(translocation=="13q.loss", "del(13q)",
                                                                        ifelse(translocation=="14q.loss", "del(14q)",
                                                                        ifelse(translocation=="16q.loss", "del(16q)",
                                                                        ifelse(translocation=="17p.loss", "del(17p)",  
                                                                        ifelse(translocation=="1p.loss", "del(1p)",
                                                                        ifelse(translocation=="1q.gain", "gain(1q)",  
                                                                        ifelse(translocation=="6q.loss", "del(6q)", 
                                                                        ifelse(translocation=="8p.loss", "del(8p)", translocation)))))))))

plot.c = read.table("/Users/yellapav/Desktop/record/p292/misc/fish.arm.results.txt",sep="\t",header=T,stringsAsFactors = F) %>%
  dplyr::filter(translocation!="11q.gain" & translocation!="12p.loss" & translocation!="BIRC3") %>%
  dplyr::filter(translocation!="Sequenced") %>% mutate(translocation = ifelse(translocation=="13q.loss", "del(13q)",
                                                                       ifelse(translocation=="14q.loss", "del(14q)",
                                                                        ifelse(translocation=="16q.loss", "del(16q)",
                                                                        ifelse(translocation=="17p.loss", "del(17p)",  
                                                                        ifelse(translocation=="1p.loss", "del(1p)",
                                                                        ifelse(translocation=="1q.gain", "gain(1q)",  
                                                                        ifelse(translocation=="6q.loss", "del(6q)", 
                                                                        ifelse(translocation=="8p.loss", "del(8p)", translocation)))))))))

colors=c("#06A77D","#F1A208","#035d46","#3E37A1","#005377","#D5C67A","#8D37A1")
library(ComplexHeatmap)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  MYELOMA = function(x, y, w, h) {
    grid.rect(x, y, w*0.85, h*0.88, gp = gpar(fill = "#0570b0", col = NA)) 
  }, 
  AMYLOIDOSIS = function(x, y, w, h) {
    grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = "#DC6941", col = NA))
  }, 
  SMOLDERING = function(x, y, w, h) {
    grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = "#74a9cf", col = NA)) 
  },
  PLASMA.CELL.LEUKEMIA = function(x, y, w, h) {
    grid.rect(x, y, w*0.85, h-unit(0.5, "mm"), gp = gpar(fill = "#005377", col = NA))
  },
  MGUS = function(x, y, w, h) {
    grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = "#a6bddb", col = NA))
  },
  inframe_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.85, gp = gpar(fill = "#444da3", col = NA))
  },
  complex_change_in_transcript = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.85, gp = gpar(fill = "#D5C67A", col = NA))
  }, 
  stop_lost = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.85, gp = gpar(fill = "#8D37A1", col = NA))
  }, 
  inframe_codon_gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.85, gp = gpar(fill = "#ff0900", col = NA))
  }, 
  Sequenced = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.85, gp = gpar(fill = "#A9A9A9", col = NA))
  }, 
  initiator_codon_change = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*0.8, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)


#rorder = as.character(unlist(data.frame(table(all.muts$VAG_GENE)) %>% arrange(desc(Freq)) %>% filter(Var1!="sequenced") %>% dplyr::select(Var1)))
col = c("SMOLDERING" = "#F1A208","PLASMA.CELL.LEUKEMIA"="#005377","AMYLOIDOSIS"="#3E37A1","inframe_codon_gain"="#ff0900","MGUS"="#035d46","stop_gained"="#3E37A1","stop_lost"="#8D37A1","MYELOMA"="#06A77D","initiator_codon_change"="#EA3599","Sequenced"="#A9A9A9" )
col = c("SMOLDERING" = "#74a9cf","PLASMA.CELL.LEUKEMIA"="#005377","MGUS"="#a6bddb","MYELOMA"="#0570b0","AMYLOIDOSIS"="#DC6941","Sequenced"="#A9A9A9" )





n=plot.c
colnames(n)=gsub(".T[1,2,3].1.D1.1","",colnames(n))
colnames(n)=gsub(".T1.2.D1.1",".2",colnames(n))
colnames(n)=gsub("I.H.","H",colnames(n))

cln=clinical
cln$sample=cln$Workflow.Leukid
cln$sample=gsub(".T[1,2,3].1.D1.1","",cln$sample)
cln$sample=gsub(".T1.2.D1.1",".2",cln$sample)
cln$sample=gsub("I.H.","H",cln$sample)

n=n[,(colnames(n) %in% c("translocation",unique(cln$sample)))]
mat=n
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)



cccc=clinical
cccc$SAMPLE=gsub("I-H-","H",cccc$Workflow.Leukid)
cccc$SAMPLE=gsub("-T[1,2,3]-[1,2,3]-D[1,2,3]-[1,2,3]","",cccc$SAMPLE)
cccc[cccc$Pathology.ID=="YM40",]$SAMPLE="H135270.2"


tt=data.frame(t(mat))
tt$SAMPLE=rownames(tt)
tt=inner_join(tt,cccc, by=c('SAMPLE'='SAMPLE')) 

head(tt)



corder1 = (tt) %>% mutate(level = ifelse(CATEG=="MGUS", "AMGUS",
                                  ifelse(CATEG=="SMOLDERING", "ABSMOLDERING",
                                  ifelse(CATEG=="MYELOMA", "AAMYELOMA",
                                  ifelse(CATEG=="AMYLOIDOSIS", "ZAMYLOIDOSIS",CATEG))))) %>% 
  mutate(levels = ifelse(del.13q.!="", "AZ",
                         #ifelse(X11q.gain!="", "B",
                          ifelse(gain.1q.!="", "C",
                          ifelse(t.11.14.!="", "AA",
                          ifelse(del.14q.!="", "E",
                          ifelse(del.8p.!="", "F",
                          ifelse(del.6q.!="", "G",
                          ifelse(del.17p.!="", "H",
                          ifelse(HRD!="", "J",
                          ifelse(NRAS!="", "K",
                          ifelse(KRAS!="", "I","Z"))))))))))) %>% 
  arrange((levels)) %>% arrange((level)) #%>% dplyr::select(SAMPLE)
corder=as.character(unique(corder1$SAMPLE))


rorder=c('Sequenced','13q.loss','1q.gain','14q.loss','16q.loss','8p.loss','6q.loss','17p.loss','1p.loss','HRD','t(11;14)','t(4;14)','t(8;14)','t(6;14)','t(14;16)','t(14;20)','KRAS','NRAS','TP53','FAM46C','BRAF','LTB','SP140','DIS3','IRF4','PTPN11','RB1','TRAF3')     
rorder=c('del(13q)','gain(1q)','del(14q)','del(16q)','del(8p)','del(6q)','del(17p)','del(1p)','HRD','t(11;14)','t(4;14)','t(8;14)','t(6;14)','t(14;16)','t(14;20)','KRAS','NRAS','TP53','FAM46C','BRAF','LTB','SP140','DIS3','IRF4','PTPN11','RB1','TRAF3')                                             
df = data.frame(rorder)                           
zz = plot %>% group_by(translocation,CATEG) %>% tally()                           
zz=dcast(zz,translocation ~ CATEG,sum) %>% arrange(desc(MYELOMA)) #%>% arrange(desc(Sequenced))
zz=zz[match(df$rorder,zz$translocation),]
m=zz %>% dplyr::select("AMYLOIDOSIS","MGUS","MYELOMA","PLASMA.CELL.LEUKEMIA","SMOLDERING") %>% dplyr::select("MYELOMA")
m = t(apply(m, 1, function(x) x/sum(x)))    
head(m)  


mat= mat[rorder,]  

nmat = (as.data.frame(t(mat)))

fi=paste("~/Downloads/","malin_review_oncoplot.png",sep="")
png(file = fi, res=300, width=3200, height=2200)
ht = oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, 
               column_title = "Alterations in Myeloma Samples", show_column_names = FALSE, remove_empty_columns=FALSE, column_names_gp = gpar(fontsize = 4),row_names_gp = gpar(fontsize = 12),
               #row_barplot_width = unit(2, "cm"),#bottom_annotation=HeatmapAnnotation(TP53.Loss = c(tt$CATEG)),
               column_order=corder,row_order=rownames(mat),#row_split=rep(c("A", "B"), 6), #row_order = rorder,
               heatmap_legend_param = list(title = "Disease", at = c("MYELOMA","SMOLDERING","MGUS","PLASMA.CELL.LEUKEMIA","AMYLOIDOSIS","Sequenced"), 
                                           labels = c("Multiple myeloma","Smoldering multiple myeloma","MGUS","Plasma cell leukemia","AL amyloidosis","Sequenced"))) #+
  #rowAnnotation(Proportion = row_anno_barplot(m, gp = gpar(fill = c("#DC6941","#a6bddb","#0570b0","#005377","#74a9cf","#A9A9A9")), 
  #                                            axis_param = list(direction = "normal"), row_order=sample(1:30),
  #                                            bar_width = 1, border=F),width = unit(3, "cm"))
ht
dev.off()


l=c("A",rep("CNA",9),rep("ASV",6),rep("SNV",12))
l=c(rep("CNA",9),rep("ASV",6),rep("SNV",12))
#fi=paste("/ifs/res/leukgen/projects/292/RESULTS/misc/fish/figs/","see.png",sep="")
#png(file = fi, res=300, width=3200, height=2600)
draw(ht,split=l)
#ht




