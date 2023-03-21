#'
#' Simple DPClust preprocessing pipeline that takes a VCF file with SNV calls. Then fetches
#' allele counts from the specified bam file. With the counts and the copy number fit it
#' creates a DPClust input file with all the required columns.
#' 
#' Dependencies:
#'  * The alleleCounter C utility must be in $PATH
#' 
#' v1.0 - 2016-02-29 - gundemg@mskcc.org
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/RunDPclustering.R /ifs/work/leukgen/local/opt/leukdc/workflows/W0000304F/data/W0000304F.bam /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR5/caveman/W0000304F_vs_W0000299F.caveman.muts.annot.vcf.gz /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR5/bbg/W0000304F_rho_and_psi.txt /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR5/bbg/W0000304F_subclones.txt male /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR5/dpclustering/

library(dpclust3p)

args = commandArgs(T)
bam_file = toString(args[1])
vcf_file_raw = toString(args[2])
rho_and_psi_file = toString(args[3])
subclones_file = toString(args[4])
sex = toString(args[5])
output_dir = toString(args[6])
fai_file = '/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta.fai'
if(length(args)>=7) {
	fai_file = toString(args[7])
}
ign_file = '/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/battenberg/ignored-contigs.txt'
if(length(args)>=8) {
	ign_file = toString(args[8])
}

# Filter the vcf file
samplename = unlist(strsplit(basename(vcf_file_raw), split="_vs_"))[1]
loci_file = paste(output_dir, "/", samplename, "_loci.txt", sep="")
cmd = paste("zless ",vcf_file_raw, " | awk -F\"\\t\" '$7==\"PASS\"{split($8,a,\"ASRD=\"); split(a[2],a2,\";\"); {if(a2[1]>0.90) print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5\"\\t\"$3}}' > ", loci_file, sep="")
system(cmd)

# Define various temp files and the final output file
dpplot_file = paste(output_dir, '/', samplename,'_dpplot.pdf', sep='')
dpoutput_file = paste(output_dir, '/', samplename, "_allDirichletProcessInfo.txt", sep="")
allelecounts_file = paste(output_dir, "/", samplename, "_alleleFrequencies.txt", sep="")
# Fetch allele counts
alleleCount(locifile=loci_file, bam=bam_file, outfile=allelecounts_file, min_baq=20, min_maq=35)
# Create dpIn file
runGetDirichletProcessInfo(loci_file=loci_file, allele_frequencies_file=allelecounts_file, cellularity_file=rho_and_psi_file, subclone_file=subclones_file, gender=sex, SNP.phase.file="NA", mut.phase.file="NA", output_file=dpoutput_file)

subclones = read.table(subclones_file, header=T, sep="\t", stringsAsFactors=F)
del.chrs = unique(subclones[subclones$nMin1_A==0,'chr'])
dpo = read.table(dpoutput_file, header=T)
dpo$VAF.uncorrected = dpo$mut.count / (dpo$mut.count + dpo$WT.count)
write.table(dpo, file=dpoutput_file, row.names=F, quote=F, sep="\t")
if(length(del.chrs)>0) {
	pdf(dpplot_file, width=12, height=12)
	par(mfrow=c(2,2));
	plot(density(dpo[!(dpo$chr %in% del.chrs),'subclonal.fraction'], na.rm=T), ylab="", xlab="", main='subclonal fraction\nother chrs')
	plot(density(dpo[(dpo$chr %in% del.chrs),'subclonal.fraction'], na.rm=T), ylab="", xlab="", main='subclonal fraction\tdeleted chrs')
	plot(density(dpo[!(dpo$chr %in% del.chrs),'VAF.uncorrected'], na.rm=T), main='uncorrected VAF\tother chrs',ylab="",xlab="", xlim=c(0,1))
	plot(density(dpo[(dpo$chr %in% del.chrs),'VAF.uncorrected'], na.rm=T), main='uncorrected VAF\tdeleted chrs',ylab="",xlab="", xlim=c(0,1))
} else {
	pdf(dpplot_file, width=12)
	par(mfrow=c(1,2));
	plot(density(dpo[,'subclonal.fraction'], na.rm=T), ylab="", xlab="", main='subclonal fraction')
	plot(density(dpo[,'VAF.uncorrected'], na.rm=T), main='uncorrected VAF',ylab="",xlab="", xlim=c(0,1))
}
rho_and_psi = read.table(file=rho_and_psi_file)
legend("topright",c(paste('ACF=',rho_and_psi['FRAC_GENOME','rho'], sep=''),paste('PLOIDY=',rho_and_psi['FRAC_GENOME','psi'], sep='')))
dev.off()

