setwd("~/Desktop")
vcf_files=c("E-H-109099-T1-1-D1-1.caveman.signature.common.vcf","E-H-109099-T1-1-D1-1.caveman.signature.uniq.vcf","E-H-109099-T2-1-D1-1.caveman.signature.uniq.vcf","E-H-109099-T2-1-D1-1.caveman.signature.vcf","E-H-109099-T1-1-D1-1.caveman.signature.vcf")
sample_names=c("COMMON","UNIQUE_T1","UNIQUE_T2","T2","T1")


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
