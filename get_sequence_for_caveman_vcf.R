library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(reshape2)

source('/ifs/work/leukgen/home/gg10/soft/integrated.circos/utils/toPyr.R')

args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_stub = args[2]

onlyPASSED=TRUE; onlyPassedAsrd=TRUE; genome.v="hg19"; genomeSeq=Hsapiens; addContext=TRUE; normalise=TRUE;
nwindow=50;
B35.1.2.vcf <- readVcf(input_file, genome.v)
rd <- rowRanges(B35.1.2.vcf)
fs <- rd$FILTER
fs.passed <- (fs=='PASS')
if (onlyPASSED) {
	B35.1.2.vcf <- B35.1.2.vcf[fs.passed,]
	rd <- rowRanges(B35.1.2.vcf)
	fs <- rd$FILTER
	fs.passed <- (fs=='PASS')
}
# Use soft filter for ASRD>0.90
info.data <- info(B35.1.2.vcf)
if(onlyPassedAsrd) {
	asrd.passed = info.data$ASRD>0.90
	B35.1.2.vcf <- B35.1.2.vcf[asrd.passed,]
	fs.passed <- (fs=='PASS' & asrd.passed)
}

info.data <- info(B35.1.2.vcf)
rd <- rowRanges(B35.1.2.vcf)
fs <- rd$FILTER
fs.passed <- (fs=='PASS')
rgs <- ranges(rd)
starts <- start(rgs)
ends <-  end(rgs)
chroms <- paste('chr', seqnames(rd), sep='')

fxd <- (fixed(B35.1.2.vcf))
wt <- as.character(rd$REF)
mt <- as.character(unlist(rd$ALT))

barcode <- paste(chroms, '-',starts,'-', mt, sep='')

bb <- as.character(getSeq(genomeSeq, chroms, start=starts-nwindow, end=ends-1))
ba <- as.character(getSeq(genomeSeq, chroms, start=starts+1, end=ends+nwindow))
wt.ref <- as.character(getSeq(genomeSeq, chroms, start=starts, end=ends))

mut.table <- data.frame(bbef=as.character(bb), wt=as.character(wt), mt=as.character(mt), baft=as.character(ba), stringsAsFactors=FALSE)

mut.table$pyrwt <- as.character(mut.table$wt)
mut.table$pyrmt <- as.character(mut.table$mt)
mut.table$pyrbbef <- as.character(mut.table$bbef)
mut.table$pyrbaft <- as.character(mut.table$baft)

not.pyr <- ((wt=='G') | (wt=='A'))
mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
mut.table$pyrmt[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))
mut.table$context = paste(mut.table$pyrbbef, mut.table$pyrwt, mut.table$pyrbaft, sep='')

mut.table = cbind(paste(chroms, starts, wt, mt, sep='_'), mut.table)
colnames(mut.table)[1] = 'id'

write.table(mut.table, file=paste(output_stub, '_context.txt', sep=''), row.names=F, sep="\t", quote=F)
command = paste(
	"less", paste(output_stub, '_context.txt', sep=''), "| grep -v id | awk '{print \">\"$1; print $10;}' >",
	paste(output_stub, '_context.fasta', sep=''), sep=' '
)
system(command)

#pyrbef = colsplit(string=mut.table$pyrbbef, pattern="", names=paste('pb',nwindow:1,sep=""))
#pyraft = colsplit(string=mut.table$pyrbaft, pattern="", names=paste('pa',1:nwindow,sep=""))
#mut.table.context <- cbind(pyrbef, as.character(mut.table$pyrwt), as.character(mut.table$pyrmt), pyraft, stringsAsFactors=FALSE)
#colnames(mut.table.context)[nwindow+1] = 'pyrwt'
#colnames(mut.table.context)[nwindow+2] = 'pyrmt'
#
#par(mfrow=c(5,6))
#for(k in 1:nwindow) {
#	cur.context <- paste(mut.table.context[,paste('pb',k,sep='')], '[', mut.table.context$pyrwt, '>',mut.table.context$pyrmt,']', sep='')
#	cur.context.table <- as.data.frame(table(cur.context))
#	total.muts <- sum(cur.context.table[,'Freq'])
#	signHist <- cur.context.table[,'Freq']/total.muts
#	barplot(signHist, names=as.character(cur.context.table$cur.context), las=2, main=k)
#}
