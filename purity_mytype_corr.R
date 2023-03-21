library(mclust)
#library(ASCAT)
library(DNAcopy)
library(data.table) 
library(VariantAnnotation) 
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library('gsheet')
library(dplyr)
setwd("/Users/yellapav/Desktop/p292/purity")
#clin<- read.delim("CLASSIFICATION_MSKCC.txt", stringsAsFactors = F)

names<- read.delim("/Users/yellapav/Desktop/p292/purity/data/p292_sample_sheet.txt", stringsAsFactors = F)
names<- names[,c("Workflow.Leukid","Specimen.ID","Aliquot.ID")]
names2 = names
colnames(names2)[1]<-"sample"


exclude = c("I-H-135277-T1-1-D1-1","I-H-135351-T1-1-D1-1","I-H-135304-T1-1-D1-1")
clin <- read.csv(construct_download_url("https://docs.google.com/spreadsheets/d/1vOO7ihNDfv4cElTLJT18rmpjZeJrHzjMixXkAuIgza4", format = "csv", sheetid = "1075620643"), stringsAsFactors=FALSE)
clin = clin %>% dplyr::filter(CATEG!="WALDENSTROM") %>% 
  left_join(names,by = c('Pathology.ID' = 'Aliquot.ID')) %>% 
  dplyr::filter(!(Workflow.Leukid %in% exclude))


purity<- read.delim("data/purity.mytype.ascat.txt", stringsAsFactors = F)
head(purity)


mutations <- read.csv(construct_download_url("https://docs.google.com/spreadsheets/d/1Uif7ICxqHtvqjt0IyAxHx5m16VqY6PDORH0VbPnC1yw", format = "csv", sheetid = "199311900"), stringsAsFactors=FALSE)
mutations=mutations[!(mutations$Annotation %in% c('SNP','SNP/ARTIFACT','ARTIFACT','ARTIFACT, SNP','NOT IN PANEL')),]
mut = mutations
mut3<- mut[!mut$Annotation %in% c("SNP","SNP/ARTIFACT" ,"ARTIFACT","NOT IN PANEL"),]
mut2<- mut3[,c("TARGET_NAME","CHR","START","REF","ALT","TARGET_VAF_MEAN")]

cnv<- read.delim("data/consensus_cnvkit_facets.txt", stringsAsFactors = F)

unique(cnv$sample)[! unique(cnv$sample) %in% unique(mut2$TARGET_NAME) ]



gr0 = with(mut2, GRanges(CHR, IRanges(start=(START), end=(START))))
values(gr0) <- DataFrame(ref = mut2$REF,
                         alt = mut2$ALT, 
                         vaf = mut2$TARGET_VAF_MEAN,
                         sample = mut2$TARGET_NAME)


gr1 = with(cnv, GRanges(chrom, IRanges(start=start, end=end)))
values(gr1) <- DataFrame(tot = cnv$tcn.em,
                         minor = cnv$lcn.em,
                         sample = cnv$sample
                          )

ranges2 <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
ranges2 <- ranges2[with(ranges2, startB <= startA & endB >= endA),]
ranges3<- ranges2[which(ranges2$sampleA == ranges2$sampleB),]
head(ranges3)

file<- ranges3[,c("sampleA","tot","minor","seqnames","startA","ref","alt","vaf")]
file_diplo<- file[file$tot ==2 & file$minor==1,]

sam<- unique(file_diplo$sampleA)
length(sam)
all<- list()
for(i in (1:length(sam)))
{
  if(max(na.omit(as.numeric(as.character(file_diplo[file_diplo$sampleA == sam[i],]$vaf))))>0.15){
  all[[i]]<- c(sam[i],mean(na.omit(as.numeric(as.character(file_diplo[file_diplo$sampleA == sam[i],]$vaf)))))
  }
}
all2<- do.call("rbind", all)


head(all2)
all2<- as.data.frame(all2)
colnames(all2)<-c("sample","vaf")
colnames(purity)[2]<-"sample"

def<- merge(purity, all2, by="sample")
def$vaf <- as.numeric(as.character(def$vaf))
def$pur_my<- def$vaf * 2
def$pur_my[def$pur_my>1]<-1
def$purity.ascat[def$sample == "I-H-135270-T1-2-D1-1"]<-0.9 ### wrong ASCAT
def$pur_my[def$sample=="I-H-135322-T1-1-D1-1"]<- 0.7 ### one clonal with many subclonal
def$pur_my[def$sample=="I-H-135327-T1-1-D1-1"]<- 0.7 ### one clonal with many subclonal
def$pur_my[def$sample=="I-H-135274-T1-1-D1-1"]<- 0.5 ### one clonal with many subclonal

plot(def$purity.ascat, def$pur_my, pch=16,ylim=c(0,1), xlim=c(0,1) )
summary(lm(def$pur_my~def$purity.ascat, method = 'spearman'))

getwd()
write.table(def, "data/purity_fm.txt", sep="\t")


mut<- read.delim("data/VDJ.mutations.minimalCols.txt", stringsAsFactors = F)

names<- read.delim("data/cel_samples.txt", stringsAsFactors = F)
names2<- names[,c("sample","Workflow.Leukid","Aliquot.ID")]
pur <- read.csv(construct_download_url("https://docs.google.com/spreadsheets/d/14BkSwcuGKIISk4eBzsXxp6Hb9tyOi7xrsHxUddItkPU", format = "csv", sheetid = "119655860"), stringsAsFactors=FALSE)
pur = pur %>% filter(include2!="X") %>% mutate(purity = ifelse(include2=="alt", include,ascat.output.aberrantcellfraction))
pur[pur$sample=="20160627_161857_005_CG162828-C,G-05",]$sample="20160627_161857_005_CG16-2828-C,G-05"
#pur<- read.delim("p292.SNParray.myTYPE.unique.Feb26 - ascat.purity.tsv", stringsAsFactors = F)
pur2<- pur[pur$include!="X",]
ascat<- merge(names2,pur2, by="sample")
colnames(ascat)[1:2]<-c("id","sample")

head(mut)
mut_cov<- mut[mut$strelka_TARGET_DEPTH>400,]
mut_cov_igh<- mut_cov[mut_cov$START>106300000,]

mut_cov_igh_high<- mut_cov_igh[mut_cov_igh$TARGET_VAF_MEAN>0.6,]
mut_cov_igh_low<- mut_cov_igh[mut_cov_igh$TARGET_VAF_MEAN<0.6,]

sam<- unique(mut_cov_igh_low$TARGET_NAME)
length(sam)
all_low<- list()
for(i in (1:length(sam)))
{
  if(max(na.omit(as.numeric(as.character(mut_cov_igh_low[mut_cov_igh_low$TARGET_NAME == sam[i],]$TARGET_VAF_MEAN))))>0.15){
    all_low[[i]]<- c(sam[i],mean(na.omit(as.numeric(as.character(mut_cov_igh_low[mut_cov_igh_low$TARGET_NAME == sam[i],]$TARGET_VAF_MEAN)))))
  }
}
all_low2<-do.call("rbind", all_low)
all_low2<- as.data.frame.matrix(all_low2)
colnames(all_low2)<-c("sample","low")

all_high<- list()
for(i in (1:length(sam)))
{
  if(max(na.omit(as.numeric(as.character(mut_cov_igh_high[mut_cov_igh_high$TARGET_NAME == sam[i],]$TARGET_VAF_MEAN))))>0.15){
    all_high[[i]]<- c(sam[i],mean(na.omit(as.numeric(as.character(mut_cov_igh_high[mut_cov_igh_high$TARGET_NAME == sam[i],]$TARGET_VAF_MEAN)))))
  }
}
all_high2<-do.call("rbind", all_high)
all_high2<- as.data.frame.matrix(all_high2)
colnames(all_high2)<-c("sample","high")

vdj<- merge(all_low2,all_high2,by="sample")
vdj$high<- as.numeric(as.character(vdj$high))
vdj$low<- as.numeric(as.character(vdj$low))
vdj$purity<- ((vdj$high/2)+vdj$low)

vdj_final<- merge(ascat, vdj, by="sample")
plot(vdj_final$purity.y, vdj_final$ascat.output.aberrantcellfraction, pch=16,
     ylim=c(0,1), xlim=c(0,1))


head(def)


vdj_final2<- vdj_final[! vdj_final$sample %in% def$sample, ]
vdj_final3<- vdj_final2[vdj_final2$notes!="X",]

plot(vdj_final3$purity.y, vdj_final3$ascat.output.aberrantcellfraction, pch=16,
     ylim=c(0,1), xlim=c(0,1))
summary(lm(vdj_final3$ascat.output.aberrantcellfraction~vdj_final3$purity.y, method = 'spearman'))



vdj_final3_short<- vdj_final3[,c("sample","purity.y","ascat.output.aberrantcellfraction")]
def_short<- def[,c("sample","pur_my","purity.ascat")]
colnames(vdj_final3_short)<-colnames(def_short)

final_purity<- rbind.data.frame(def_short, vdj_final3_short)

plot(final_purity$pur_my, final_purity$purity.ascat, pch=16,
     ylim=c(0,1), xlim=c(0,1))
summary(lm(final_purity$purity.ascat~final_purity$pur_my, method = 'spearman'))



##### outlayers

final_purity[final_purity$pur_my<0.4 & final_purity$purity.ascat>0.5,]

final_purity[final_purity$pur_my>0.6 & final_purity$purity.ascat<0.5,]

vdj_final3[vdj_final3$sample == "I-H-135252-T1-1-D1-1",]
#def2[def2$sample=="I-H-135327-T1-1-D1-1",]


final_purity$purity.ascat[final_purity$sample=="I-H-135270-T1-1-D1-1"]<-0.75
final_purity$purity.ascat[final_purity$sample=="I-H-135307-T1-1-D1-1"]<-0.7
final_purity$purity.ascat[final_purity$sample=="I-H-135252-T1-1-D1-1"]<-0.75


summary(lm(final_purity$purity.ascat~final_purity$pur_my, method = 'spearman'))
plot(final_purity$pur_my, final_purity$purity.ascat, pch=16,
     ylim=c(0,1), xlim=c(0,1))
abline(lm(final_purity$purity.ascat~final_purity$pur_my))



library(ggplot2)
ggplot(final_purity, aes(x=purity.ascat, y=pur_my, color="dodgerblue")) + 
  geom_point( show.legend = FALSE) +theme_bw(base_size=35) +
  geom_smooth(method=lm, data = final_purity, col="red") +
  # geom_smooth(method=lm, data = sig_nhl_all_nhl2, col="dodgerblue4") +
  xlim(0.2,1) + ylim(0.2,1) +labs(x = "SNP Array Purity", y="myTYPE Purity")+ 
  theme(
    # plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(size=32), 
    axis.title.y = element_text(size=32), 
    axis.text.x = element_text(angle = 0, hjust = 1, size=30), 
    axis.text.y = element_text(angle = 0, hjust = 1, size=30), 
  )+ annotate("text", x = 0.3, y = 0.9, label = "italic(r) ^ 2 == 0.61", parse = TRUE,size=10) +
  annotate("text", x = 0.3, y = 0.85, label = "p < 0.001", parse = TRUE,size=10) +
  theme( panel.border = element_blank(),
         panel.grid = element_blank(),
         panel.background = element_blank(), axis.line.x = element_line(size = 0.5, linetype = "solid"),
         axis.line.y = element_line(size = 0.5, linetype = "solid")) 



ggsave("purity_mytype.pdf", width = 10, height = 10)
system("open purity_mytype.pdf")

cor.test(final_purity$pur_my,final_purity$purity.ascat, method = "pearson")

png(file = "./purity_mytype_hist.png", res=300, width=1200, height=1000)
ggplot(final_purity,aes(as.numeric(as.character(pur_my))))+geom_histogram(binwidth=0.02,fill="#d7301f")+xlab("myTYPE Purity")+ theme_bw(base_size = 13)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position="none", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0.2), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+coord_cartesian(xlim=c(0,1))+ylab("Count")
dev.off()


write.table(final_purity,file="final_purity.txt", append=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA, quote=FALSE)

ggplot(final_purity, aes(pur_my)) + 
  geom_histogram( show.legend = FALSE)+
  theme_bw(base_size=18)



# 
# def[def$pur_my>0.3 & def$purity.ascat <0.3,]
# 
# setwd("/Users/mauraf/LUNA/fm6/mytype/ascat/")
# f<- list.files()
# f[grep("135270", f)]
# 
# # ASCAT ==> 90 "I-H-135270-T1-2-D1-1",]
# 
# 
# DEF<- merge(def, names2, by="sample")
# 
# 
# # remove for mytype low purity
# # "I-H-130664-T2-1-D1-1"
# 
# DEF[DEF$pur_my<0.4 & DEF$purity.ascat >0.5,]
# 
# mut2[mut2$TARGET_NAME %in% "I-H-135259-T1-1-D1-1",]
# file[file$sampleA=="I-H-135270-T1-2-D1-1",]
# 
# 
# def_inc<- def[def$pur_my>0.3,]
# 
# plot(def_inc$purity.ascat, def_inc$pur_my, pch=16,ylim=c(0,1), xlim=c(0,1) )
# summary(lm(def_inc$pur_my~def_inc$purity.ascat, method = 'spearman'))
# 
# cnv[cnv$sample=="I-H-106909-T3-1-D1-1",]
# mut2[mut2$TARGET_NAME== "I-H-130664-T2-1-D1-1",]
