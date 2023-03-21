#################################################
### Usage: Rscript3.3 make_signature.R /ifs/res/leukgen/local/opt/leukdc/projects/163/analyses/WGS-2-P150-X10-NYGC/yellapav/VCF_MAF/signature/snv.input /ifs/res/leukgen/local/opt/leukdc/projects/163/analyses/WGS-2-P150-X10-NYGC/yellapav/VCF_MAF/signature/snv.pdf

#######################################################

library(MutationalPatterns)
library(BSgenome)


ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


args = commandArgs(TRUE)
outFile = toString(args[2])
inFile = toString(args[1])
#libdir = toString(args[1])

#inFile="/ifs/res/leukgen/local/opt/leukdc/projects/163/analyses/WGS-2-P150-X10-NYGC/yellapav/VCF_MAF/signature/snv.input"
#outFile="/ifs/res/leukgen/local/opt/leukdc/projects/163/analyses/WGS-2-P150-X10-NYGC/yellapav/VCF_MAF/signature/see.pdf"

r=read.table(inFile,sep="\t",header=T,stringsAsFactors=F)


#vcf_files=c("E-H-109099-T1-1-D1-1.caveman.signature.common.vcf","E-H-109099-T1-1-D1-1.caveman.signature.uniq.vcf","E-H-109099-T2-1-D1-1.caveman.signature.uniq.vcf","E-H-109099-T2-1-D1-1.caveman.signature.vcf","E-H-109099-T1-1-D1-1.caveman.signature.vcf")
#sample_names=c("COMMON","UNIQUE_T1","UNIQUE_T2","T2","T1")


########## Specifically for SOHN ##########
r$sample=gsub("H-109","",r$sample)
r$sample=gsub("-1-D1-1","",r$sample)
##########################################
r$sample

vcf_files=r$vcf_file
sample_names=r$sample

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
pdf(outFile, width=10,height=8)
plot_96_profile(test_matrix)
dev.off()
