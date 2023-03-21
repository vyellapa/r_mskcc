#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-g", "--groups", type="character", nargs=1, help="Absolute path to file with new rearrangement breakpoints. In the same format as groups.gz file.")
parser$add_argument("-b", "--bedpe", type="character", nargs=1, help="Bedpe output file from the BRASS pipeline.")
parser$add_argument("-t", "--bam", type="character", nargs=1, help="BAM file for the tumour.")

# Parse arguments and complain
args <- parser$parse_args()
group_file = args$groups
bedpe_file = args$bedpe
bam_file = args$bam

library(reshape2)

# Convert the groups format ro r3 format
groups = read.table(group_file, header=F, sep="\t", stringsAsFactors=F)
groups$id = 1:nrow(groups)
r3_file = gsub('.txt','.r3',group_file)
write.table(groups[,c(1,3,4,5,7,8,11,10,2,6)], row.names=F, col.names=F, quote=F, sep="\t", r3_file)

# Identify the absolute breakpoints using the script from BRASS pipeline
extra_rearr_file = gsub('.txt','.r4',group_file)
cmd = paste(
	"/ifs/work/leukgen/opt/cgp/5.18.4/brass/4.0.5/bin/get_abs_bkpts_from_clipped_reads.pl -fasta /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta -out",
	extra_rearr_file, bam_file, r3_file, sep=" "
)
print("Identifying absolute breakpoints from clipped reads...")
system(cmd)

# Read the bedpe file from the pipeline
bedpe = read.table(bedpe_file, header=F, sep="\t", stringsAsFactors=F)
colnames(bedpe) = unlist(strsplit(gsub('# ','',system(paste0('grep brass_score ', bedpe_file), intern=T)), split="\t"))
colnames(bedpe)[grep('strand1',colnames(bedpe))[2]] = "sstrand1"
colnames(bedpe)[grep('strand2',colnames(bedpe))[2]] = "sstrand2"

# Append the extra rearrangements to those in the bedpe file
extras = colsplit(
	string=system(paste0('less ', extra_rearr_file, ' | cut -f1-12'), intern=T), pattern="\t",
	names=c('chr1','lstart1','lend1','chr2','hstart2','hend2','id/name','brass_score','astr1','astr2','bkpt1','bkpt2')
)
extras$strand1 = extras$astr1
extras$strand2 = "-"
extras$strand2[extras$astr2=="-"] = "+"
extras$svclass = ""
extras$svclass[extras$strand1=="+" & extras$strand2=="+"] = "deletion"
extras$svclass[extras$strand1=="-" & extras$strand2=="-"] = "tandem-duplication"
extras$svclass[extras$strand1=="-" & extras$strand2=="+"] = "inversion"
extras$svclass[extras$strand1=="+" & extras$strand2=="-"] = "inversion"
extras$svclass[extras$chr1!=extras$chr2] = "translocation"
extras$start1 = extras$bkpt1 - 1
extras$end1 = extras$bkpt1
extras$start2 = extras$bkpt2 - 1
extras$end2 = extras$bkpt2
extras$sample = grep('N1',unique(bedpe$sample),invert=T, value=T)
max_id = max(bedpe$"id/name")
extras$"id/name" = c((max_id + 1000):(max_id + 1000 + nrow(extras)-1))
shared_columns = intersect(colnames(bedpe), colnames(extras))
bedpe_specific_columns = setdiff(colnames(bedpe), colnames(extras))
bedpe_specific_data = as.data.frame(mat.or.vec(nrow(extras),length(bedpe_specific_columns)))
colnames(bedpe_specific_data) = bedpe_specific_columns
extras = cbind(extras[,shared_columns], bedpe_specific_data)
bedpe = rbind(bedpe, extras[,colnames(bedpe)])

# Create a new bedpe file
new_bedpe_file = gsub('annot.bedpe','annot.new.bedpe',bedpe_file)
system(paste0('grep "#" ', bedpe_file, ' > ', new_bedpe_file))
write.table(bedpe, new_bedpe_file, row.names=F, col.names=F, quote=F, sep="\t", append=T)
file.remove(r3_file)
file.remove(extra_rearr_file)
