packages <- c("readr", "tidyr", "dplyr", "ggplot2", "car", "broom",
              "reshape2","foreach","doParallel","parallel","DESeq2")
invisible(suppressWarnings(suppressMessages(lapply( packages, library, 
                                                    character.only = TRUE))))
options(scipen=999)
set.seed(42)

onc=read.table("~/Downloads/oncpoplot_adaptive.txt", header=T, sep="\t")



library(ComplexHeatmap)
onc=onc[,c(1,15,17:26)]

l = list()
onc1 = onc[,c(1,2,3)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[1]] = onc1


onc1 = onc[,c(1,2,4)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[2]] = onc1

onc1 = onc[,c(1,2,5)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[3]] = onc1

onc1 = onc[,c(1,2,6)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[4]] = onc1

onc1 = onc[,c(1,2,7)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[5]] = onc1

onc1 = onc[,c(1,2,8)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[6]] = onc1

onc1 = onc[,c(1,2,9)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[7]] = onc1

onc1 = onc[,c(1,2,10)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[8]] = onc1

onc1 = onc[,c(1,2,11)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[9]] = onc1

onc1 = onc[,c(1,2,12)]
colnames(onc1) = c("SAMP","CAPTURE","MUT")
l[[10]] = onc1

plot = do.call(rbind, l)
plot = distinct(plot)

p1 = plot[,c(1,2)]
colnames(p1)[2]="MUT"
plot1 = rbind(p1,plot[,c(1,3)]) #%>% dplyr::filter(MUT!="Failure")
plot1$MUT = gsub("[*]","",plot1$MUT)
plot1$MUT = gsub("_p.+","",plot1$MUT)
plot1 = plot1 %>% dplyr::filter(MUT!="N") %>% dplyr::filter(SAMP!="0")
plot1 = plot1 %>% dplyr::mutate(pass = ifelse(MUT=="Failure",NA,"pass;")) %>% 
  dplyr::mutate(MUT=ifelse((MUT=="Failure" | MUT=="Success"),"Clonotype detection",MUT))
plot1 = distinct(plot1)


plot1 = dcast(plot1, MUT ~ SAMP, value.var = "pass")
mat=plot1
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  pass = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4575b4", col = NA))
  },
  fail = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#d73027", col = NA))
  },
  Translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#f46d43", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#B8DC3C", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#1a9850", col = NA))
  }, 
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#2b8cbe", col = NA))
  }, 
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)

col = c("Translocation" = "#f46d43", "fail" = "#d73027", "pass" = "#4575b4", "Splice_Site"="#B8DC3C", "Missense_Mutation" = "#1a9850", "Frame_Shift_Del"="#a6cee3", "In_Frame_Del" = "#EA3599", "Nonsense_Mutation" = "#a6cee3", "Frame_Shift_Ins" = "#2b8cbe" )

ht = oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in all Samples", show_column_names = FALSE,
          heatmap_legend_param = list(title = "Alternations", at = c("fail", "pass", "Translocation", "Splice_Site", "Missense_Mutation", "Frame_Shift_Del", "In_Frame_Del", "Nonsense_Mutation", "Frame_Shift_Ins"), 
                                      labels = c("Amplification", "Success/Present", "Translocation", "Splice Site", "Missense", "Frame Shift Del", "In-Frame Del", "Nonsense_Mutation", "Frame_Shift_Ins")))



#l=c(rep("CNA",9),rep("ASV",6),rep("SNV",12))
l=rownames(mat) %>% as.data.frame() %>% dplyr::rename(GENE=".") %>% 
  dplyr::mutate(val="MUTATION") %>% dplyr::mutate(i=grepl("Del|Gain",GENE)) %>% 
  dplyr::mutate(val=ifelse(i,"CNA",val)) %>% dplyr::mutate(i=grepl('t',GENE)) %>% 
  dplyr::mutate(val=ifelse(i,"ASV",val)) %>% dplyr::mutate(val=ifelse(GENE=="Clonotype detection","A",val)) %>%
  dplyr::select(val) 
l=as.vector(l$val)


draw(ht,split=l)
