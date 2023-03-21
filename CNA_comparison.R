# CN5 CN2 comparison

a=read.table("/Users/yellapav/Desktop/SOHN_163/VCF_MAF/misc/CN5_2.plot",sep="\t",header = F,stringsAsFactors = F)


colnames(a)=c("CHR","POS","POS1","REF","ALT","POS","REF1","COUNT","MAPQ","BASEQ","SE_MAPQ","POS_STRAND_READS","NEG_STRAND_READS","POS_AS_FRAC","MM_AS_FRAC","MMQS","NUM_Q2","AVG_DIST_Q2","CLIPPED_LEN","DIST_3P","ALT1","ALT_COUNT","ALT_MAPQ","ALT_BASEQ","ALT_SE_MAPQ","ALT_POS_STRAND_READS","ALT_NEG_STRAND_READS","ALT_POS_AS_FRAC","ALT_MM_AS_FRAC","ALT_MMQS","ALT_NUM_Q2","ALT_AVG_DIST_Q2","ALT_CLIPPED_LEN","ALT_DIST_3P","CATEGORY","SAMPLE","MP")
a$STRANDEDNESS=apply(a,1,function(x) {p=as.numeric(as.character(x["ALT_POS_STRAND_READS"]));n=as.numeric(as.character(x["ALT_NEG_STRAND_READS"]));m=min(p,n)/(p+n); return(m)})
a$VAF=apply(a,1,function(x) {p=as.numeric(as.character(x["COUNT"]));n=as.numeric(as.character(x["ALT_COUNT"]));m=n/(p+n); return(m)})

a$DP=apply(a,1,function(x) {p=as.numeric(as.character(x["COUNT"]));n=as.numeric(as.character(x["ALT_COUNT"]));m=(p+n); return(m)})
a$MMQS_DIFF=abs(a$ALT_MMQS-a$MMQS)

a=a[which(a$ALT_POS_AS_FRAC>0.3),]


write.table(a,file="~/Desktop/sig.tsv", append=FALSE, sep="\t", eol="\n",quote = F, row.names=FALSE, col.names=TRUE)

ggplot(a,aes(ALT_COUNT,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+coord_cartesian(xlim=c(0,150),ylim=c(0,0.15))+theme_bw()
ggplot(a,aes(ALT_MAPQ,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(45,60),ylim=c(0,2))

ggplot(a,aes(ALT_BASEQ,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(ylim=c(0,0.4))
ggplot(a,aes(STRANDEDNESS,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggplot(a,aes(ALT_POS_AS_FRAC,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+geom_vline(xintercept = 0.01)
ggplot(a,aes(ALT_MM_AS_FRAC,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(ylim=c(0,300),xlim=c(0.01,0.065))
ggplot(a,aes(ALT_MMQS,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,120),ylim=c(0,0.2))
ggplot(a,aes(ALT_AVG_DIST_Q2,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggplot(a,aes(VAF,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggplot(a,aes(MMQS_DIFF,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,100),ylim=c(0,0.15))

ggplot(a,aes(ALT_CLIPPED_LEN,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(125,150),ylim=c(0,0.75))
ggplot(a,aes(ALT_DIST_3P,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()

ggplot(a,aes(MP,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0.80,1),ylim=c(0,25))

ggplot(a,aes(DP,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,750),ylim=c(0,0.015))




setwd("/Users/yellapav/Desktop/SOHN_163/VCF_MAF/CN5")
vcf_files=c("e99_t1_common.vcf","e99_t1_s1.vcf","e99_t2_common.vcf","e99_t2_s1.vcf","e99_t1_filter.vcf","e99_t2_filter.vcf")
sample_names=c("SHARED_T1","UNIQ_CN5_T1","SHARED_T2","UNIQ_CN5_T2","UNIQ_CN5_T1_FILTER","UNIQ_CN5_T2_FILTER")


ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = "hg19")

#vcf=read_vcfs_as_granges("c.vcf", "sample_names", genome = "hg19")

auto = extractSeqlevelsByGroup(species="Homo_sapiens",
                               style="UCSC",
                               group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
mutation_types(vcfs)
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(test_matrix)


