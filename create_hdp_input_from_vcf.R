#
# Create input for HDP. Mutation counts for 96 nucleotides with 5' and 3' context
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/mutational.signatures/create_hdp_input_from_vcf.R /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/mutational_signatures/368_create_hdp_input_from_vcf_subclones.input mutational_signatures/368_create_hdp_input_from_vcf_subclones.output F
#
#

library(BSgenome.Hsapiens.UCSC.hg19) # human genome
library(VariantAnnotation) # 

UTIL.DIR = '/ifs/work/leukgen/home/gg10/soft/integrated.circos/utils'
SIGN.CONTEX.FILE = '/ifs/work/leukgen/home/gg10/soft/integrated.circos/data/signatures_probabilities.txt'
source(paste(UTIL.DIR,'/toPyr.R',sep=''))
source(paste(UTIL.DIR,'/generateHist.R', sep=''))
source(paste(UTIL.DIR,'/processSubs.R',sep=''))
source(paste(UTIL.DIR,'/merge.with.order.R',sep=''))
source(paste(UTIL.DIR,'/generateHist.R',sep=''))
source(paste(UTIL.DIR,'/getMutTables.R',sep=''))

args <- commandArgs(trailingOnly = TRUE)
sample.info <- args[1]
out <- args[2]
filter <- args[3]

sh <- read.table(SIGN.CONTEX.FILE, header=TRUE, sep='\t')
mut.order <- (sh[,'Somatic.Mutation.Type'])

sample.info = read.table(file=sample.info, stringsAsFactors=F, header=T)
result = NULL
for(i in 1:nrow(sample.info)) {
	sample = as.character(sample.info$sample[i])
	print(sample)
	if(filter=="T") {
		subs.data <- getMutTables(as.character(sample.info$vcf_file[i]), onlyPASSED=TRUE, onlyPassedAsrd=TRUE)
	} else {
		subs.data <- getMutTables(as.character(sample.info$vcf_file[i]), onlyPASSED=FALSE, onlyPassedAsrd=FALSE)
	}
	result = rbind(result, subs.data[['all.hist']])
}
rownames(result) = sample.info$sample
write.table(result, file=out)
result = read.table(file=out, sep=' ')

format.cols <- function(x) {
	return(paste(c(
		unlist(strsplit(colnames(result)[x], split='\\.'))[c(2,3)],
		paste(unlist(strsplit(colnames(result)[x], split='\\.'))[c(1,2,4)], sep='', collapse='')
	), sep='\\.', collapse='.'))
}
colnames(result) <- sapply(1:length(colnames(result)), format.cols)
write.table(result, file=out)
