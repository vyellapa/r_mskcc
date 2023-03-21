setwd("/Users/yellapav/Desktop/p292/scripts/data")
library(ComplexHeatmap)
library(dplyr)
library(reshape2)

p1 = read.table("p1.txt",header=T, sep="\t") 
clinical = read.table("clinical2.txt",header=T, sep="\t") 

p1.total = p1 %>% group_by(Sample) %>% tally()
p1.perc = p1 %>% group_by(Sample, ASSAY) %>% tally() %>% left_join(p1.total, by=c('Sample' = 'Sample')) %>% mutate(BOTH.percent = ifelse(ASSAY=="BOTH",n.x/n.y,0)) %>% 
mutate(PATH.percent = ifelse(ASSAY=="Pathology",n.x/n.y,0)) %>% mutate(myTYPE.percent = ifelse(ASSAY=="myTYPE",n.x/n.y,0))

p1.both = p1.perc %>% filter(BOTH.percent==1) %>% arrange(desc(n.x))
p1.mixed = p1.perc %>% filter(!(Sample %in% unique(p1.both$Sample))) %>% arrange(desc(n.x))


### make 1 data frame with aligned PATH and myTYPE percentage
p1.mixed.both = unique(p1.mixed[p1.mixed$ASSAY=="BOTH",c("Sample","BOTH.percent")])
p1.mixed.PATH = unique(p1.mixed[p1.mixed$ASSAY=="Pathology",c("Sample","PATH.percent")])
p1.mixed.myTYPE = unique(p1.mixed[p1.mixed$ASSAY=="myTYPE",c("Sample","myTYPE.percent")])

p1.mixed.combined = left_join(p1.mixed.both,p1.mixed.PATH ) %>% left_join(p1.mixed.myTYPE ) %>% 
mutate(PATH.percent = ifelse(is.na(PATH.percent), 0, PATH.percent), myTYPE.percent = ifelse(is.na(myTYPE.percent), 0, myTYPE.percent))



p1.MP = p1.mixed.combined %>% filter(PATH.percent>0 & myTYPE.percent > 0) %>% arrange(desc(BOTH.percent))
p1.M = p1.mixed.combined %>% filter(PATH.percent==0 & myTYPE.percent > 0) %>% arrange(desc(BOTH.percent))
p1.P = p1.mixed.combined %>% filter(PATH.percent>0 & myTYPE.percent == 0) %>% arrange(desc(BOTH.percent))


mms = (p1) %>% left_join(clinical, by = c('Sample' = 'Workflow.Leukid'))  %>% filter(CATEG == "MYELOMA") %>% dplyr::select(Sample,CATEG) %>% unique()
head(mms)
p1.mm = p1 %>% filter(Sample %in% unique(mms$Sample))
p1.nmm = p1 %>% filter(!(Sample %in% unique(mms$Sample)))
nrow(p1.both)
nrow(p1.mixed)


#p2.both = p1[p1$Sample %in% unique(p1.both$Sample),]
new=dcast(p1.mm,translocation~Sample,paste, collapse = ";", value.var=c("ASSAY"))
#write.table(new, "~/new.txt",append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE, quote=FALSE)

colors=c("#4CAE55","#4C4537","#E7AB1B")
colors=c("#06A77D","#035d46","#F1A208","#3E37A1","#005377","#D5C67A")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  BOTH = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#035d46", col = NA))
  },
  Pathology = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#752F3C", col = NA))
  },
  myTYPE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F1A208", col = NA))
  },
   Sana = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#abb9c9", col = NA))
  },
  Sequenced = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#ffffff", col = NA))
  }
  
)

col = c("myTYPE" = "#F1A208", "Pathology" = "#752F3C", "BOTH" = "#035d46", "Sana" = "#abb9c9", "Sequenced" = "#ffffff" )

order=c('1p','1q','6q','8p','13q','14q','16q','17p','HRD','t(11;14)','t(14;16)','t(14;20)','t(4;14)','t(6;14)','t(8;14)')


n=new
colnames(n)=gsub("-T[1,2,3]-1-D1-1","",colnames(n))
colnames(n)=gsub("-T1-2-D1-1",".2",colnames(n))
colnames(n)=gsub("I-H-","H",colnames(n))

mat=n
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

mat=gsub("myTYPE;BOTH","BOTH",mat)
mat=gsub("Pathology;BOTH","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)
mat=gsub("Pathology;myTYPE","BOTH",mat)
mat=gsub("myTYPE;Pathology","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)

fi=paste("../",as.character("fig1_"),"MM.png",sep="")
png(file = fi, res=300, width=2500, height=4200)
oncoPrint(t(mat), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col, 
          column_title = "", show_column_names = TRUE, show_row_names = FALSE, remove_empty_columns=FALSE, column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 7), column_order=order, 
          show_pct=F, right_annotation = rowAnnotation(width = unit(5, "cm"),rbar = anno_oncoprint_barplot(bar_width=0.8)), row_split=c(rep(c("MM"),115)),column_split = c(rep(c("CNA"), 9), rep(c("Translocations"), 6)),
          heatmap_legend_param = list(title = "Assay", at = c("BOTH", "Pathology", "myTYPE", "Sana"), 
                                      labels = c("myTYPE & Array/FISH", "Array/FISH", "myTYPE", "NA")))

dev.off()


#p2.mixed = p1[p1$Sample %in% unique(p1.mixed$Sample),]

order=c('1p','1q','6q','8p','13q','14q','16q','17p','HRD','t(11;14)','t(14;16)','t(14;20)','t(4;14)','t(8;14)')
new=dcast(p1.nmm,translocation~Sample,paste, collapse = ";", value.var=c("ASSAY"))


nmms = (p1.nmm) %>% inner_join(clinical, by = c('Sample' = 'Workflow.Leukid'))  %>% filter(CATEG != "MYELOMA") %>% dplyr::select(Sample,CATEG) %>% unique() %>% left_join(p1.total) %>% 
mutate(CATEG = ifelse(CATEG=="PLASMA.CELL.LEUKEMIA","ZZ",as.character(CATEG))) %>% arrange(CATEG,desc(n))

#nmms = (p1) %>% left_join(clinical, by = c('Sample' = 'Workflow.Leukid'))  %>% filter(CATEG != "MYELOMA") %>% dplyr::select(Sample,CATEG) %>% unique()
#nmms = (nmms) %>% left_join(p1.total) %>% arrange(CATEG,desc(n))
row_order=gsub("-T[1,2,3]-1-D1-1","",nmms$Sample)
row_order=gsub("-T1-2-D1-1",".2",row_order)
row_order=gsub("I-H-","H",row_order)


n=new
colnames(n)=gsub("-T[1,2,3]-1-D1-1","",colnames(n))
colnames(n)=gsub("-T1-2-D1-1",".2",colnames(n))
colnames(n)=gsub("I-H-","H",colnames(n))

mat=n
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

mat=gsub("myTYPE;BOTH","BOTH",mat)
mat=gsub("Pathology;BOTH","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)
mat=gsub("Pathology;myTYPE","BOTH",mat)
mat=gsub("myTYPE;Pathology","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)


mat = mat[,(row_order)]

fi=paste("../",as.character("fig1_"),"nonMM.png",sep="")
png(file = fi, res=300, width=2400, height=3000)
oncoPrint(t(mat), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col, row_order = row_order,
          column_title = "", show_column_names = TRUE, show_row_names = FALSE, remove_empty_columns=FALSE, column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 7), column_order=order, 
          show_pct=F, right_annotation = rowAnnotation(width = unit(4, "cm"),rbar = anno_oncoprint_barplot(bar_width=0.8)), 
	column_split = c(rep(c("CNA"), 9), rep(c("Translocations"), 5)),row_split = c(rep(c("Amyloidosis"), 9), rep(c("MGUS"), 4), rep(c("SMM"), 17), rep(c("ZPCL"), 1)),
         heatmap_legend_param = list(title = "Assay", at = c("BOTH", "Pathology", "myTYPE", "Sana"),
                                      labels = c("myTYPE & Array/FISH", "Array/FISH", "myTYPE", "NA")))

dev.off()


