#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-s", "--spfile", type="character", nargs=1, help="Path to the sample2purity file used for the full run that generated the Gibbs sampling data.")
parser$add_argument("-f", "--fdir", type="character", nargs=1, help="Path to the output directory from the full run that contains the Gibbs sampling data.")
parser$add_argument("-d", "--ddir", type="character", nargs=1, help="Path to the directory that contains the DP input files used for the Gibbs sampling.")
parser$add_argument("-o", "--oparallel", type="character", nargs=1, help="Output directory where DP_output directory will created.")
parser$add_argument("-p", "--patient", type="character", nargs=1, help="Patient name to be used as job id while submitting jobs to LSF.")

args <- parser$parse_args()
sample2purity_file = args$spfile
dir_full_run = args$fdir
dir_dp_input = args$ddir
parallel_output = args$oparallel
patient = args$patient
rundp_script = "/ifs/work/leukgen/home/gg10/soft/dpclust_pipeline_v0.2.2/parallel_assignments_gg10/RunDP_gg10_2.R"
rundp_script = "/home/yellapav/local/msk_scripts/RunDP_gg10_2.R"

# Function to create permutations
get.subsample.indices<-function(choose.index,choose.number,choose.from){
	subsample.indices = 1:choose.number
	temp.choose.index = choose.index
	for(i in 1:choose.number){
		last.subtotal = 0
		subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
		while(temp.choose.index>subtotal){
			subsample.indices[i] = subsample.indices[i] + 1
			last.subtotal = subtotal
			subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)
		}
		if(i<choose.number){
			subsample.indices[i+1] = subsample.indices[i]+1
		}
		temp.choose.index = temp.choose.index - last.subtotal
	}
	return(subsample.indices)
}
############################################################################################################
sample2purity = read.table(sample2purity_file, header=T, stringsAsFactors=F)
## Create the subsample info file
choose.number = 3
choose.from = nrow(sample2purity)
perm.info = NULL
no.perms = choose(choose.from,choose.number)
for(choose.index in 1:no.perms){
	perm.info = rbind(perm.info, sample2purity$samplename[get.subsample.indices(choose.index,choose.number,choose.from)])
}
write.table(perm.info, 'permutation_info.txt', row.names=F, col.names=F, quote=F, sep="\t")

for(choose.index in 1:no.perms) {
	cmd = paste(
	'bsub -oo dpp_',choose.index,'.log -W 480 -R "rusage[mem=16]" ', "'docker-circos Rscript ", rundp_script, " ", parallel_output," ",
	choose.index, " ", dir_full_run, " ", sample2purity_file, " ", dir_dp_input,"'", sep=""
	)
	print(cmd)
	system(cmd)
}
