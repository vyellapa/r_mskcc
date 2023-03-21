#'
#' Simple DPClust preprocessing pipeline that takes a VCF file with SNV calls. Then fetches
#' allele counts from the specified bam file. With the counts and the copy number fit it
#' creates a DPClust input file with all the required columns.
#' 
#' Dependencies:
#'  * The alleleCounter C utility must be in $PATH
#' 
#' v1.0 - 2016-02-29 - gundemg@mskcc.org
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/RunDP_pipeline_input.R dpclustering_conf_file dpclustering_wdir [vcf_type:cgp/standard] [do_pileups:T/F]
# Rscript /ifs/work/leukgen/home/gg10/soft/RunDP_pipeline_input.R /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/nd_dp_clustering_input_current/dpcluster.conf JHBone368/nd_dpclustering_vs_JHLiver368 cgp T

library(dpclust3p)
library(hexbin)

args = commandArgs(T)
dp.conf = toString(args[1])
output_dir = toString(args[2])
vcf_type = args[3]
run_pileup = 'F'
if(length(args)>=4) {run_pileup='T'}

dp.conf = read.table(file=dp.conf, header=T, sep="\t", stringsAsFactors=F)

loci_file1 = paste(output_dir, '/loci1.txt', sep='')
loci_file2 = paste(output_dir, '/loci2.txt', sep='')
if(run_pileup=='T') {
	# Filter the vcf file
	for(i in 1:nrow(dp.conf)) {
		cmd = ""
		if(vcf_type=='cgp') {
			cmd = paste(
			"zless ",dp.conf[i,'vcf_file'],
			" | awk -F\"\\t\" '$7==\"PASS\"{split($8,a,\"ASRD=\"); split(a[2],a2,\";\"); {if(a2[1]>0.90) print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}}' >> ",
			loci_file1, sep=""
			)
		} else if (vcf_type=='standard') {
			cmd = paste(
			"zless ",dp.conf[i,'vcf_file'],
			" | awk -F\"\\t\" '$7==\"PASS\"{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}' >> ",
			loci_file1, sep=""
			)
		}
		system(cmd)
		system(paste("sort -u ",loci_file1," | grep -v X | sort -k1 -nk2.2 > ",loci_file2, sep=" "))
		system(paste("sort -u ",loci_file1," | grep X | sort -k1 -nk2.2 >> ",loci_file2, sep=" "))
	}
} 
# Define various temp files and the final output file
for (i in 1:nrow(dp.conf)) {
	dpoutput_file = paste(output_dir, '/', dp.conf[i,'samplename'], "_allDirichletProcessInfo.txt", sep="")
	allelecounts_file = paste(output_dir, "/", dp.conf[i,'samplename'], "_alleleFrequencies.txt", sep="")
	# Fetch allele counts
	alleleCount(locifile=loci_file2, bam=dp.conf[i,'bam_file'], outfile=allelecounts_file, min_baq=20, min_maq=35)
	# Create dpIn file
	runGetDirichletProcessInfo(
		loci_file=loci_file2, allele_frequencies_file=allelecounts_file,
		cellularity_file=dp.conf[i,'rho_and_psi_file'], subclone_file=dp.conf[i,'subclones_file'],
		gender=dp.conf[i,'gender'], SNP.phase.file="NA", mut.phase.file="NA", output_file=dpoutput_file
	)
}

for(i in 1:nrow(dp.conf)) {
	for(k in 1:nrow(dp.conf)) {
		if(i!=k) {
		dpo1 = read.table(paste(output_dir, "/", dp.conf[i,'samplename'], "_allDirichletProcessInfo.txt", sep=""), header=T, stringsAsFactors=F)
		dpo2 = read.table(paste(output_dir, "/", dp.conf[k,'samplename'], "_allDirichletProcessInfo.txt", sep=""), header=T, stringsAsFactors=F)
		dpo1$VAF.uncorrected = dpo1$mut.count / (dpo1$mut.count + dpo1$WT.count)
		dpo2$VAF.uncorrected = dpo2$mut.count / (dpo2$mut.count + dpo2$WT.count)
		df = data.frame(dpo1$subclonal.fraction, dpo2$subclonal.fraction)
		colnames(df) = c(dp.conf[i,'samplename'],dp.conf[k,'samplename'])
		df = df[df[,1]<=2 & df[,2]<=2,]
		h = hexbin(df, xbins=100)
		pdf(paste(output_dir, '/dpcluster__', dp.conf[i,'samplename'], '_vs_', dp.conf[k,'samplename'], '.pdf', sep=''))
		plot(h)
		dev.off()
		}
	}
}

