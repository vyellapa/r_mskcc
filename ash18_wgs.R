as.data.frame(table(zz$inv.VAF))
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
 
m=read.table("/Users/yellapav/Desktop/subs_indels/caveman_wgs/I-H-106917-T2-1-D1-2.caveman.maf",header=T,sep="\t",stringsAsFactors = FALSE, quote="~")
rem=as.data.frame(m[grep("^rs",m[,c(57)],perl = TRUE),]$variant_id)
colnames(rem)=c("V1")
g=read.table("/Users/yellapav/Desktop/subs_indels/caveman_wgs/I-H-106917-T2-1-D1-2.dnds.vcf",header=F,sep="\t")


#g=g[!(g$V3 %in% rem$V1),]
g$s="I-H-106917-T2-1-D1-2"
g=g[,c("s","V1","V2","V4","V5")]
g=read.table("/Users/yellapav/Desktop/RESULTS/misc/signatures/dn_ds_sub.in",header=F,sep="\t")
g=read.table("/Users/yellapav/Desktop/RESULTS/misc/signatures/dn_ds_clo.in",header=F,sep="\t")
colnames(g)=c("sampleID","chr","pos","ref","mut")
dndsout = dndscv(g)
print(dndsout$globaldnds)

cgc=read.table("~/local/resources/cgc_1col.tsv")
cgc_genes=sel_cv[which(sel_cv$gene_name %in% cgc$V1),]
signif_genes = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv")]
rownames(signif_genes) = NULL
print(signif_genes)


signif_genes = sel_cv[sel_cv$ptrunc_cv<0.1, c("gene_name","ptrunc_cv")]
rownames(signif_genes) = NULL
print(signif_genes)



##################################################################################
################ Clonal and Subclonal 96nt #######################################
##################################################################################

vcf_files=c("/Users/yellapav/Desktop/RESULTS/misc/signatures/I-H-106917_clone.vcf","/Users/yellapav/Desktop/RESULTS/misc/signatures/I-H-106917_subclone.vcf","/Users/yellapav/Desktop/RESULTS/misc/signatures/I-H-130719_clone.vcf","/Users/yellapav/Desktop/RESULTS/misc/signatures/I-H-130719_subclone.vcf")
sample_names=c("I-H-106917_clone","I-H-106917_subclone","I-H-130719_clone","I-H-130719_subclone")

ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)

#vcf=read_vcfs_as_granges("c.vcf", "sample_names", genome = "hg19")

auto = extractSeqlevelsByGroup(species="Homo_sapiens",
                               style="UCSC",
                               group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
mut_type(vcfs)
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(test_matrix, condensed = TRUE, ymax=0.1)
plot_compare_profiles(test_matrix, condensed = TRUE)


##################################################################################
################ Enrichment in Enhancers #########################################
##################################################################################
surveyed_file <- system.file("extdata/callableloci-sample.bed", package = "MutationalPatterns")
enh=read.table("/Users/yellapav/Downloads//Merged.MM.primary.K27ac.remove2.5k.bed", sep="\t",header=FALSE)
#enh$V1=gsub("chr","",enh$V1)
  enh=reduce(GRanges(enh$V1, IRanges(enh$V2, enh$V3)))
  CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds", package="MutationalPatterns"))
  seqlevels(CTCF_g)=paste("chr",seqlevels(CTCF_g), sep="")
  promoter_g <- readRDS(system.file("states/promoter_g_data.rds",package="MutationalPatterns"))
  seqlevels(promoter_g)=paste("chr",seqlevels(promoter_g), sep="")
  flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds", package="MutationalPatterns"))
  seqlevels(flanking_g)=paste("chr",seqlevels(flanking_g), sep="")


  regions = GRangesList(promoter_g, flanking_g, CTCF_g, enh)
  
  #regions = GRangesList(enh)
  names(regions) <- c("Promoter", "Promoter flanking", "CTCF", "Enhancers")
  #names(regions) = c("Enhancers")

library(rtracklayer)
surveyed <- import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"
#seqlevelsStyle(surveyed) <- "NCBI"
## For this example we use the same surveyed file for each sample.
surveyed_list <- rep(list(surveyed), 4)
distr <- genomic_distribution(vcfs, surveyed_list, regions)
distr_test <- enrichment_depletion_test(distr)
plot_enrichment_depletion(distr_test)



#############################################################################

g=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g2=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g3=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g4=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")

ca=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917-T2-1-D1-2_DPoutput_10000iters_1000burnin/I-H-106917-T2-1-D1-2_10000iters_1000burnin_bestConsensusAssignments.bed",header=T,sep="\t")
ca=unique(ca)
g$key=paste(g$chr,g$start,g$end, sep="_")
g2$key=paste(g2$chr,g2$start,g2$end, sep="_")
g3$key=paste(g3$chr,g3$start,g3$end, sep="_")
g4$key=paste(g4$chr,g4$start,g4$end, sep="_")

ca$key=paste(ca$chr,ca$start,ca$end, sep="_")
m=merge(g,ca, by.x="key", by.y="key")
mm=m[m$cluster==11,]
m14=m[m$cluster==14,]
nrow(m14[m14$subclonal.CN<2,])/nrow(m14[m14$subclonal.CN>=2,])
#nrow(m14[m14$subclonal.CN==2,])/nrow(m14[m14$subclonal.CN!=2,])
nrow(mm[mm$subclonal.CN<2,])/nrow(mm[mm$subclonal.CN>=2,])


ggplot(m14,aes(subclonal.CN))+geom_density()+xlab("I-106917-1 subclonal.CN Clone 14")+annotate("text", x = 4.5, y = 1.5,label="nrow(CN==2)/nrow(CN!=2)=0.33")+theme_bw()
ggplot(mm,aes(subclonal.CN))+geom_density()+xlab("I-106917-1 subclonal.CN Clone 11")+annotate("text", x = 4.5, y = 1.5,label="nrow(CN==2)/nrow(CN!=2)=0.34")+theme_bw()


##################################################################################
g=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g2=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g3=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")
g4=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917/I-H-106917-T2-1-D1-2_allDirichletProcessInfo.txt",header=T,sep="\t")

ca=read.table("/Users/yellapav/Desktop/RESULTS/misc/dirichlet_step1/I-H-106917-T2-1-D1-2_DPoutput_10000iters_1000burnin/I-H-106917-T2-1-D1-2_10000iters_1000burnin_bestConsensusAssignments.bed",header=T,sep="\t")
ca=unique(ca)
g$key=paste(g$chr,g$start,g$end, sep="_")
g2$key=paste(g2$chr,g2$start,g2$end, sep="_")
g3$key=paste(g3$chr,g3$start,g3$end, sep="_")
g4$key=paste(g4$chr,g4$start,g4$end, sep="_")


#g=read.table("/Users/yellapav/Desktop/subs_indels/caveman_wgs/I-H-106917-T2-1-D1-2.dnds.vcf",header=F,sep="\t")

am=read.table("/Users/yellapav/Desktop/subs_indels/caveman_wgs/dndn_all_917.txt",header=F,sep="\t")

#g=g[!(g$V3 %in% rem$V1),]
am$s="I-H-106917-T2-1-D1-2"
am=am[,c("s","V1","V2","V3","V4")]
am$key=paste(am$V1, am$V2, sep="_")

m1=m[m$cluster==11,]
m2=m[m$cluster==7,]
m3=m[m$cluster==4,]
m4=m[m$cluster==3,]
mm=rbind(m1,m2,m3,m4)

mm$key=paste(mm$chr.x, mmm$end.x, sep="_")
mm=merge(mm,am, by.x="key", by.y="key")

gg2=mm[mm$subclonal.CN>=2,]
gl2=mm[mm$subclonal.CN<2,]
gg2=gg2[,c("s","V1","V2","V3","V4")]
colnames(gg2)=c("sampleID","chr","pos","ref","mut")

gl2=gl2[,c("s","V1","V2","V3","V4")]
colnames(gl2)=c("sampleID","chr","pos","ref","mut")

dndsout = dndscv(gl2)
print(dndsout$globaldnds)

dndsout = dndscv(gg2)
print(dndsout$globaldnds)

################################################################
############## Summary of Aberrations in WGS ###################
################################################################
ss=read.table("/Users/yellapav/Desktop/RESULTS/misc/signatures/summary_aberrations_wgs.txt",header=F,sep="\t",stringsAsFactors = FALSE, quote="~")
colnames(ss)=c("Sample","Type","Count")
ss$Patient = substr(ss$Sample,5,10)

summary(ss[which(ss$Type=="SNV" & ss$Patient=="106917"),c(3)])
sd(ss[which(ss$Type=="SNV" & ss$Patient=="106917"),c(3)])
summary(ss[which(ss$Type=="Indels" & ss$Patient=="106917"),c(3)])
sd(ss[which(ss$Type=="Indels" & ss$Patient=="106917"),c(3)])
summary(ss[which(ss$Type=="SV" & ss$Patient=="106917"),c(3)])
sd(ss[which(ss$Type=="SV" & ss$Patient=="106917"),c(3)])


summary(ss[which(ss$Type=="SNV" & ss$Patient=="130719"),c(3)])
sd(ss[which(ss$Type=="SNV" & ss$Patient=="130719"),c(3)])
summary(ss[which(ss$Type=="Indels" & ss$Patient=="130719"),c(3)])
sd(ss[which(ss$Type=="Indels" & ss$Patient=="130719"),c(3)])
summary(ss[which(ss$Type=="SV" & ss$Patient=="130719"),c(3)])
sd(ss[which(ss$Type=="SV" & ss$Patient=="130719"),c(3)])


###################################
maf=read.table("/Users/yellapav/Desktop/RESULTS/misc/mafs/I-H-106917-T2-1-D1-2.caveman.maf",header=T,sep="\t",stringsAsFactors = FALSE, quote="~")
maf$key=paste(maf$Chromosome,maf$Start_Position-1, sep="_")
ca$key=paste(ca$chr,ca$start, sep="_")
merged=merge(maf,ca, by.x="key", by.y="key")




######################################################################
############# CNA Clustering #########################################
######################################################################

sc=read.table("/Users/yellapav/Desktop/p220_wgs/misc/I-H-106917-T2-2-D1-2_subclones.txt",sep="\t", header=TRUE)
ggplot(sc[!is.na(sc$nMaj2_A),], aes(frac2_A))+geom_histogram(bins=25)+theme_bw()
ggplot(sc[!is.na(sc$nMaj2_A),], aes(frac2_A))+geom_histogram(bins=25)+theme_bw()+facets_wrap(.~chr)

subc=sc[!is.na(sc$nMaj2_A),]


d <- dist(subc$frac2_A, method = "euclidean")
fit <- kmeans(subc$frac2_A, 3) 
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=3) # cut tree into 3 clusters
# draw dendogram with red borders around the 3 clusters
rect.hclust(fit, k=3, border="red") 
######################################################################

library(ggplot2)
library(reshape2)
a=read.table("numbers",header=F,sep=" ")
ggplot(a,aes(V2,V1))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")+ylab("")

a=read.table("DP_VAF.plot",header=F,sep="\t")

aa=read.table("t",header=T,sep="\t")
m=melt(aa)
ggplot(m,aes(variable,Gene_AA))+geom_tile(aes(fill=value))+xlab("")+ylab("")+
  theme_bw()+ 
  scale_fill_gradientn(colours = hm.palette(100),na.value = '#ffffbf')+geom_text(aes(label=value))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),panel.grid.minor.y = element_blank())+theme(legend.title=element_blank())+xlim("I.H.106917.T2.1","I.H.106917.T2.2","I.H.106917.T2.3","I.H.106917.T2.4","I.H.130718.T1.1","I.H.130718.T1.2","I.H.130718.T1.3","I.H.130718.T1.4","I.H.130718.T1.5","I.H.130718.T1.6","I.H.130718.T1.7","I.H.130718.T1.8","I.H.130718.T1.9","I.H.130718.T1.10","I.H.130719.T1.1","I.H.130719.T1.2","I.H.130719.T1.3","I.H.130719.T1.4","I.H.130719.T1.5","I.H.130719.T1.6","I.H.130720.T1.1","I.H.130720.T1.2","I.H.130720.T1.3","I.H.130720.T1.4","I.H.130720.T1.5","I.H.130720.T1.6","I.H.130720.T1.7","I.H.130720.T1.8")+
  ylim("XPO1_p.P392R","FBXW7_p.E107K","KMT2D_p.P3408Q","MAF_p.P110H","MAP3K1_p.R313L","PIK3CA_p.D578E","PRDM1_p.H705Y","RB1_p.W99L","SETD2_p.R236K","APC_p.M701I","ARID2_p.G1605C","ATM_p.Q1636H","ATM_p.R2244T","FAM46C_p.G217W","BIRC3_p.V461G","KRAS_p.G12D","NRAS_p.Q61R","TP53_p.C176Y","BRAF_p.V600E")


ggplot(m,aes(variable,Gene_AA))+geom_tile(aes(fill=value))+xlab("")+ylab("")+
  theme_bw()+ 
  scale_fill_gradient(low = "#ffffbf", high = "#006837",limits=c(0, 50), breaks=seq(0,5,by=0.25))+geom_text(aes(label=value))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())



