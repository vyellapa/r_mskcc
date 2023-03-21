#!/usr/bin/env Rscript

#
# gundemg@mskcc.org
# Brass for related samples from the same patient
#

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-r", "--rearr", type="character", nargs=1, help="Absolute path to the BEDPE format rearrangement file.")
parser$add_argument("-c", "--cnngs", type="character", nargs=1, help="Absolute path to the NGS CN file: Raw bin counts from Brass.")
parser$add_argument("-o", "--out", type="character", nargs=1, help="Output prefix. Please include a directory and prefix. Multiple images will be created.")
parser$add_argument("-s", "--subc", default=NULL, type="character", nargs=1, help="Optional. Absolute path to the Battenberg subclones file.")
parser$add_argument("-b", "--baf", default=NULL, type="character", nargs=1, help="Optional. Absolute path to the Battenberg logR-BAF file.")
parser$add_argument("-l", "--legend", default=NULL, type="character", nargs=1, help="Optional. Legend position. Not drawn if 'none'.")
parser$add_argument("-t", "--typec", default="type", type="character", nargs=1, help="Optional. Coloring scheme to use for junctions: 'custom' for user-specified color | 'type' for junction type.")
parser$add_argument("-i", "--itype", default="pdf", type="character", nargs=1, help="Optional. Type of the output file: png or pdf ")
parser$add_argument("-f", "--flip", action="store_true", default=FALSE, help="Optional. Flip the strand for the higher bp window. Junction strands should be in alignment orientation.")
parser$add_argument("-q", "--logscale", action="store_true", default=FALSE, help="Optional. Log2-transform the copy number data. Useful with high-level amplifications.")
parser$add_argument("-w", "--win", default=1e5, type="double", nargs=1, help="Optional. Window size for summarizing the raw bin counts data for the specified genonimic region only.")
parser$add_argument("-y", "--yrange", default=NULL, type="character", nargs=1, help="Optional. Range for CN values on y-axis.")
parser$add_argument("-n", "--neutralcns", default=NULL, type="character", nargs=1, help="Optional. Neural copy number state: 1:1, 2:2, 3:3, etc. Default=1:1")

#####################################
source('/ifs/work/leukgen/home/gg10/soft/plot_rearrangements_copynumber_and_battenberg.R')
library(quantsmooth)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19) # human genome
library(VariantAnnotation) # 
#####################################
chr_lens = read.table("/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta.fai", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
temp = rownames(chr_lens)
chr_lens = chr_lens[,1]
names(chr_lens) = temp

###################################################
# MAIN

# Parse arguments and complain
args <- parser$parse_args()
bedpe_file = args$rearr
bg_file = args$cnngs
bb_file = args$subc
baf_file = args$baf
legend_pos = args$legend
rearr_col = args$typec
image_type = args$itype
flip_higher = args$flip
cn_log_scale = args$logscale
cn_win_size = args$win
yrange = args$yrange
outfile = args$out
if(!is.null(yrange)) {
	yrange = unlist(strsplit(yrange, split="_"))
	yrange = c(as.double(yrange[1]), as.double(yrange[2]))
}
if(is.null(args$neutralcns)) {
	neutral_cns = c(1,1)
} else {
	neutral_cns = as.integer(unlist(strsplit(args$neutralcns, split=":")))
}
########################################################################
# Read the breakpoints into a bedpe format data frame
cat(paste('Reading rearrangements...', '\n',sep=''))
count = as.integer(unlist(strsplit(system(paste('wc -l ', bedpe_file, sep=''), intern=TRUE), split=" "))[1])
headered = as.integer(unlist(strsplit(system(paste('grep color ', bedpe_file, ' | wc -l', sep=''), intern=TRUE), split=" "))[1])
if(count>0) {
	if(headered==0) {
		brass = read.table(file=bedpe_file, sep='\t', stringsAsFactors=F, header=F, comment.char="#")
		brass$id = rownames(brass)
		bedpe = cbind(brass[,c(1:7)],'black',brass[,9:10])
		bedpe = cbind(bedpe,'a')
	}
	if(headered==1) {
		brass = read.csv(file=bedpe_file, sep='\t', stringsAsFactors=F, header=T, comment.char="")
		bedpe = brass[,unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id color strand1 strand2', split=' '))]
		bedpe = cbind(bedpe,brass[,'presence.summary'])
	}
	colnames(bedpe) = c('m_chrl','m_l5','m_l3','m_chrh','m_h5','m_h3','id','color','strl','strh','presence.summary')
	bedpe$color = as.character(bedpe$color)
	bedpe$m_l5 = as.integer(bedpe$m_l5)
	bedpe$m_l3 = as.integer(bedpe$m_l3)
	bedpe$m_h5 = as.integer(bedpe$m_h5)
	bedpe$m_h3 = as.integer(bedpe$m_h3)
} else {
	bedpe = NULL
}
cat(paste('Read rearrangements...', '\n',sep=''))
if(flip_higher) {
	strh = bedpe$strh
	strh[bedpe$strh == '-'] = '+'
	strh[bedpe$strh == '+'] = '-'
	bedpe$strh = strh
}
########################################################################
# Read the raw bin counts data, Battenberg segments and raw BAF data
cat(paste('Reading NGS copy number data...', '\n',sep=''))
bg = read.table(file=bg_file, header=F)
########################################################################
# Read the BAF data
cat(paste('Reading BAF data...', '\n',sep=''))
tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
cmd = paste("zless ", baf_file, " | grep -v \"#\" | awk 'OFS=\"\t\"{split($8,a,\";\"); split(a[1],b,\"=\"); split(a[2],d,\"=\"); split(a[3],c,\"=\"); print $1,$2,b[2],c[2],d[2]}' > ", tmp.file.name, sep="")
system(cmd)
baf_data = read.table(tmp.file.name, header=F, stringsAsFactors=F)
system(paste('rm ', tmp.file.name, sep=""))
colnames(baf_data) = c('Chromosome','Position','BAFE','BAFS','BAFP')
# Order BAF data by chromosome: 1, 2, 3, ...
baf_data_sorted = NULL
chr_order = as.character(unique(baf_data$Chromosome))
#chr_order = c(sort(as.integer(chr_order[chr_order!="X"])),"X")
chr_order = sort(as.integer(chr_order[chr_order!="X"]))
prev_chr_len = 0
baf_data$plot_position = 0
for(i in 1:length(chr_order)) {
	baf_data[baf_data$Chromosome==chr_order[i],'plot_position'] = baf_data[baf_data$Chromosome==chr_order[i],'Position'] + prev_chr_len
	chr_len = chr_lens[chr_order[i]]
	prev_chr_len = prev_chr_len + chr_len
	baf_data_sorted = rbind(baf_data_sorted, baf_data[baf_data$Chromosome==chr_order[i],])
}
# Read Battenberg segments
bb_segments = NULL
bb_segments = read.table(file=bb_file, header=T)
bb_segments$length = bb_segments$endpos - bb_segments$startpos
cat(paste('Read BBG segements...', '\n',sep=''))
neutral_cns = c(1,1)
sel_neutral = is.na(bb_segments$nMaj2_A) & bb_segments$nMaj1_A==neutral_cns[1] & bb_segments$nMin1_A==neutral_cns[2]
sel_clonal = !sel_neutral & is.na(bb_segments$nMaj2_A)
sel_subclonal = !sel_neutral & !is.na(bb_segments$nMaj2_A)
clonal_segs = bb_segments[sel_clonal,]
subclonal_segs = bb_segments[sel_subclonal,]
clonal_segs = rbind(clonal_segs, subclonal_segs[which(subclonal_segs$frac1_A<0.11 | subclonal_segs$frac2_A<0.11 | subclonal_segs$length<1e6),])
subclonal_segs = subclonal_segs[which(subclonal_segs$frac1_A>=0.11 & subclonal_segs$frac2_A>=0.11 & subclonal_segs$length>=1e6),]
baf_data_sorted$color = "darkgrey"
for(k in 1:nrow(clonal_segs)) {
	sel = baf_data_sorted$Chromosome==clonal_segs$chr[k] & baf_data_sorted$Position>=clonal_segs$startpos[k] & baf_data_sorted$Position<=clonal_segs$endpos[k]
	baf_data_sorted$color[sel] = "brown"
}
for(k in 1:nrow(subclonal_segs)) {
	sel = baf_data_sorted$Chromosome==subclonal_segs$chr[k] & baf_data_sorted$Position>=subclonal_segs$startpos[k] & baf_data_sorted$Position<=subclonal_segs$endpos[k]
	baf_data_sorted$color[sel] = "salmon"
}

# If no coordinate is specified, plot rearrangements, BAF and BBG segments for all chromosomes seperately
cat(paste('Plotting...', '\n',sep=''))
bg_data = bg[bg$V1 %in% chr_order,]
bg_data$V1 = as.integer(bg_data$V1)
pdf(outfile, width=30, height=12)
par(mfrow=c(2,1), mai=c(0.1,0.82,0.82,0.42))
# Plot BAF
plot(c(), xlim=c(1,max(baf_data_sorted$plot_position)), ylim=c(0,1), type = "n", xaxt = "n", xlab="", ylab="BAF", yaxs="i", xaxs="i")
points(y=baf_data_sorted$BAFP, x=baf_data_sorted$plot_position, col=baf_data_sorted$color, pch=".")
abline(h=0.5,lty=1,col="lightgrey")
prev_chr_len = 0
for(i in 1:length(chr_order)) {
        chr_border_pos = max(baf_data_sorted$plot_position[baf_data_sorted$Chromosome==chr_order[i]])
        abline(v=chr_border_pos, lty=1, col="lightgrey")
}
# Plot logR and rearrangements
if(legend_pos=="none") {legend_pos=NULL}
plot_rearrangements(
	bedpe=bedpe[bedpe$m_chrl %in% chr_order | bedpe$m_chrh %in% chr_order,], chrs=chr_order, cn_bedgraph=bg,
	segments=NULL, lwd=2, yrange=NULL, rearr_col=rearr_col, legend_pos=legend_pos, cn_win_size=10e5, ideogram=FALSE
)
dev.off()

