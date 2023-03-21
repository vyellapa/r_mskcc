#!/ifs/work/leukgen/bin/R/3.2.3/bin/Rscript
library(reshape2)
#suppressPackageStartupMessages(library(argparse))
#parser <- ArgumentParser()
#parser$add_argument("-d", "--dir", type="character", nargs=1, help="Directory containing symlinks to the output files from the annotation pipeline.")
#parser$add_argument("-s", "--sample", type="character", nargs=1, help="Sample name to be used as prefix for output files.")
#args <- parser$parse_args()
#tc_dir = args$dir
#sample = args$sample
args = commandArgs(trailingOnly=TRUE)
tc_dir = args[1]
sample = args[2]

write_vcf <- function(locidata, output_vcf_file) {
	cmd = paste('cp /ifs/work/leukgen/home/gg10/soft/mutational.signatures/vcf_header.vcf ', output_vcf_file, sep="")
	system(cmd)
	loci = cbind(locidata$CHR, locidata$START, ".", locidata$REF, locidata$ALT, ".","PASS",".",".",".",".")
	loci[loci[,1]==23,1] = "X"
	loci[loci[,1]==24,1] = "Y"
	colnames(loci) = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOUR")
	write.table(loci, file=output_vcf_file, col.names=F, row.names=F, quote=F, sep="\t", append=T)
}

diagnostic_plots <- function(subdata, title, plot_snp_vaf=T) {
	flags = grep('^FLAG_',colnames(merged),value=T)
	pdata = sort(colSums(subdata[,flags]), decreasing=T)/nrow(subdata)
	barplot(pdata[pdata>=0.01], names=sub('FLAG_','',names(pdata[pdata>=0.01])), las=2)
	with_rs = which(!is.na(subdata$VEP_Existing_variation) | !is.na(subdata$VEP_gnomAD_genome))
	no_rs = which(is.na(subdata$VEP_Existing_variation) & is.na(subdata$VEP_gnomAD_genome))
	# Barplot for RS id
	pdata = c(length(with_rs), length(no_rs))
	pdata = pdata / nrow(subdata)
	barplot(pdata, names=c('with RS','no RS'))
	# Density plots for VAF
	plot(density(subdata$TARGET_VAF_MEAN, na.rm=T), lwd=2, main='VAF for all calls', xlim=c(0,0.6))
	lines(density(subdata$REFERENCE_VAF_MEAN, na.rm=T), col="orange", lwd=2)
	legend('topright', legend=c('tumor','normal'), lwd=c(2,2), col=c('black','orange'))
	if(plot_snp_vaf) {
		plot(density(subdata[with_rs,'TARGET_VAF_MEAN'], na.rm=T), main='VAF for calls with RS id', lwd=2, xlim=c(0,0.6))
		lines(density(subdata[with_rs,'REFERENCE_VAF_MEAN'], na.rm=T), col="orange", lwd=2)
	}
}

print("Reading all mutations")
print(list.files(tc_dir, pattern='.snvs.output.annot.tsv.gz$', full.names=T))
#merged = read.delim(list.files(tc_dir, pattern='.snvs.output.annot.tsv.gz$', full.names=T), header=T, sep="\t", stringsAsFactors=F)
cmd = paste('zless ', list.files(tc_dir, pattern='.snvs.output.annot.tsv.gz$', full.names=T), ' | grep -v "##"', sep="")
capture = system(cmd, intern=T)
merged = colsplit(capture[2:length(capture)], pattern="\t", names=unlist(strsplit(capture[1], split="\t")))
merged = merged[which(!grepl('GL',merged$CHR) & !grepl('MT',merged$CHR)),]
sel = which(grepl('C',merged$PASSED_BY) & merged$CAVEMAN_ASRD<=0.9)
merged[sel,'PASSED_BY'] = sub('C','',sub('C,','',merged[sel,'PASSED_BY']))
merged[which(merged$PASSED_BY==""),'PASSED_BY'] = NA
# Get variant subsets
passed_caveman = merged[grepl('C',merged$PASSED_BY),]
passed_strelka = merged[grepl('S',merged$PASSED_BY),]
passed_mutect = merged[grepl('M',merged$PASSED_BY),]
passed_triple = merged[which(merged$PASSED_BY=="C,M,S"),]
passed_twoplus = merged[which(merged$PASSED_BY %in% c('C,M','C,S','M,S','C,M,S')),]
passed_strelka_mutect = merged[which(merged$PASSED_BY=="M,S"),]
passed_cavemanonly = merged[which(merged$PASSED_BY=="C"),]
passed_cavemanonly$TARGET_VAF_MEAN = passed_cavemanonly$caveman_TARGET_VAF
passed_cavemanonly$REFERENCE_VAF_MEAN = passed_cavemanonly$caveman_REFERENCE_VAF
passed_strelkaonly = merged[which(merged$PASSED_BY=="S"),]
passed_strelkaonly$TARGET_VAF_MEAN = passed_strelkaonly$strelka_TARGET_VAF
passed_strelkaonly$REFERENCE_VAF_MEAN = passed_strelkaonly$strelka_REFERENCE_VAF
passed_strelkaonly_pass_rs_and_vum = passed_strelkaonly[which(passed_strelkaonly$FLAG_VUM==0 & is.na(passed_strelkaonly$VEP_Existing_variation) & is.na(passed_strelkaonly$VEP_gnomAD_genome)),]
passed_strelkaonly_fail_rs_and_vum = passed_strelkaonly[which(passed_strelkaonly$FLAG_VUM==1 | !is.na(passed_strelkaonly$VEP_Existing_variation) | !is.na(passed_strelkaonly$VEP_gnomAD_genome)),]
passed_mutectonly = merged[which(merged$PASSED_BY=="M"),]
passed_mutectonly$TARGET_VAF_MEAN = passed_mutectonly$mutect_TARGET_VAF
passed_mutectonly$REFERENCE_VAF_MEAN = passed_mutectonly$mutect_REFERENCE_VAF
passed_mutectonly_pass_rs_and_vum = passed_mutectonly[which(passed_mutectonly$FLAG_VUM==0 & is.na(passed_mutectonly$VEP_Existing_variation) & is.na(passed_mutectonly$VEP_gnomAD_genome)),]
passed_mutectonly_fail_rs_and_vum = passed_mutectonly[which(passed_mutectonly$FLAG_VUM==1 | !is.na(passed_mutectonly$VEP_Existing_variation) | !is.na(passed_mutectonly$VEP_gnomAD_genome)),]
# Go to the new output directory
system(paste('mkdir ', tc_dir, '/diagnostics', sep=""))
setwd(paste(tc_dir, '/diagnostics', sep=""))
# Write VCF files for the different subsets
write_vcf(passed_caveman, paste(sample, '_caveman.vcf',sep=""))
write_vcf(passed_strelka, paste(sample, '_strelka.vcf',sep=""))
write_vcf(passed_mutect, paste(sample, '_mutect.vcf',sep=""))
write_vcf(passed_triple, paste(sample, '_triple.vcf',sep=""))
write_vcf(passed_twoplus, paste(sample, '_twoplus.vcf',sep=""))
write_vcf(passed_strelka_mutect, paste(sample, '_strelkamutect.vcf',sep=""))
write_vcf(passed_strelka_mutect[passed_strelka_mutect$FLAGS_ALL=='PASS',], paste(sample, '_strelkamutect_passfl.vcf',sep=""))
write_vcf(passed_strelka_mutect[passed_strelka_mutect$FLAGS_ALL!='PASS',], paste(sample, '_strelkamutect_failfl.vcf',sep=""))
write_vcf(passed_twoplus[passed_twoplus$FLAGS_ALL=='PASS' & is.na(passed_twoplus$VEP_Existing_variation) & is.na(passed_twoplus$VEP_gnomAD_genome),], paste(sample, '_twoplus_filtered.vcf',sep=""))
write.table(
	passed_twoplus[passed_twoplus$FLAGS_ALL=='PASS' & is.na(passed_twoplus$VEP_Existing_variation) & is.na(passed_twoplus$VEP_gnomAD_genome),],
	row.names=F, sep="\t", quote=F, file=paste(sample, '_twoplus_filtered.tsv',sep="")
)
write_vcf(passed_cavemanonly, paste(sample, '_cavemanonly.vcf',sep=""))
write_vcf(passed_strelkaonly, paste(sample, '_strelkaonly.vcf',sep=""))
write_vcf(passed_mutectonly, paste(sample, '_mutectonly.vcf',sep=""))
write_vcf(passed_strelkaonly_pass_rs_and_vum, paste(sample, '_strelkaonly_pass_rs_and_vum.vcf',sep=""))
write_vcf(passed_strelkaonly_fail_rs_and_vum, paste(sample, '_strelkaonly_fail_rs_and_vum.vcf',sep=""))
write_vcf(passed_mutectonly_pass_rs_and_vum, paste(sample, '_mutectonly_pass_rs_and_vum.vcf',sep=""))
write_vcf(passed_mutectonly_fail_rs_and_vum, paste(sample, '_mutectonly_fail_rs_and_vum.vcf',sep=""))
# Count mutations at 96 trinucleotide types
conf_file = "conf.txt"
cmd = paste('echo "sample__vcf_file" | ', "sed 's/__/\t/g' > ", conf_file, sep="")
system(cmd)
write.table(t(c('caveman',paste(sample, '_caveman.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka',paste(sample, '_strelka.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('mutect',paste(sample, '_mutect.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('triple',paste(sample, '_triple.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('twoplus',paste(sample, '_twoplus.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_mutect',paste(sample, '_strelkamutect.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_mutect_passfl',paste(sample, '_strelkamutect_passfl.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_mutect_failfl',paste(sample, '_strelkamutect_failfl.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('caveman_only',paste(sample, '_cavemanonly.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_only',paste(sample, '_strelkaonly.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_only_pass_rs_and_vum',paste(sample, '_strelkaonly_pass_rs_and_vum.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('strelka_only_fail_rs_and_vum',paste(sample, '_strelkaonly_fail_rs_and_vum.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('mutect_only',paste(sample, '_mutectonly.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('mutect_only_pass_rs_and_vum',paste(sample, '_mutectonly_pass_rs_and_vum.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
write.table(t(c('mutect_only_fail_rs_and_vum',paste(sample, '_mutectonly_fail_rs_and_vum.vcf',sep=""))), file=conf_file, row.names=F, col.names=F, sep="\t", quote=F, append=T)
cmd = paste("/ifs/work/leukgen/home/gg10/soft/mutational.signatures/get_normalised_mutational_profile.R -c ", conf_file, " -o .", sep="")
system(cmd)
tri_context_counts = read.table('normalised_pyrimidine_trinucleotide_counts.txt', header=T)
# Print counts
print('Caveman\tStrelka\tMutect\tTriple\tTwo_plus\tStrelka_Mutect\tStrelka_Mutect_PASS\tStrelka_Mutect_FAIL\tCaveman_only\tStrelka_only\tMutect_only\tStrelka_only_pass_rs_and_vum\tStrelka_only_fail_rs_and_vum\tMutect_only_pass_rs_and_vum\tMutect_only_fail_rs_and_vum')
print(c(
	nrow(passed_caveman), nrow(passed_strelka), nrow(passed_mutect), nrow(passed_triple), nrow(passed_twoplus),
	nrow(passed_strelka_mutect), length(which(passed_strelka_mutect$FLAGS_ALL=='PASS')), length(which(passed_strelka_mutect$FLAGS_ALL!='PASS')),
	nrow(passed_cavemanonly), nrow(passed_strelkaonly), nrow(passed_mutectonly),
	nrow(passed_strelkaonly_pass_rs_and_vum), nrow(passed_strelkaonly_fail_rs_and_vum), nrow(passed_mutectonly_pass_rs_and_vum), nrow(passed_mutectonly_fail_rs_and_vum)
))
# Diagnostic plots
pdf('diagnostics.pdf', width=20, height=40)
par(mfrow=c(9,4))
max_y = max(tri_context_counts[,2:ncol(tri_context_counts)])
barplot(
tri_context_counts[,'caveman'], ylim=c(0,max_y+0.05), main=paste('caveman N=',nrow(passed_caveman), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_caveman, 'caveman', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka'], ylim=c(0,max_y+0.05), main=paste('strelka N=',nrow(passed_strelka), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelka, 'strelka', plot_snp_vaf=F)

barplot(
tri_context_counts[,'mutect'], ylim=c(0,max_y+0.05), main=paste('mutect N=',nrow(passed_mutect), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_mutect, 'mutect', plot_snp_vaf=F)

barplot(
tri_context_counts[,'triple'], ylim=c(0,max_y+0.05), main=paste('triple N=',nrow(passed_triple), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_triple, 'triple', plot_snp_vaf=F)

barplot(
tri_context_counts[,'twoplus'], ylim=c(0,max_y+0.05), main=paste('twoplus N=',nrow(passed_twoplus), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_twoplus, 'twoplus', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka_mutect'], ylim=c(0,max_y+0.05), main=paste('strelka_mutect N=',nrow(passed_strelka_mutect), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelka_mutect, 'strelka_mutect', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka_mutect_passfl'], ylim=c(0,max_y+0.05), main=paste('strelka_mutect_passfl N=',length(which(passed_strelka_mutect$FLAGS_ALL=='PASS')), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelka_mutect[passed_strelka_mutect$FLAGS_ALL=='PASS',], 'strelka_mutect_passfl', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka_mutect_failfl'], ylim=c(0,max_y+0.05), main=paste('strelka_mutect_failfl N=',length(passed_strelka_mutect$FLAGS_ALL!='PASS'), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelka_mutect[passed_strelka_mutect$FLAGS_ALL!='PASS',], 'strelka_mutect_failfl', plot_snp_vaf=F)

barplot(
tri_context_counts[,'caveman_only'], ylim=c(0,max_y+0.05), main=paste('caveman_only N=',nrow(passed_cavemanonly), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_cavemanonly, 'caveman_only', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka_only'], ylim=c(0,max_y+0.05), main=paste('strelka_only N=',nrow(passed_strelkaonly), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelkaonly, 'strelka_only', plot_snp_vaf=F)

barplot(
tri_context_counts[,'mutect_only'], ylim=c(0,max_y+0.05), main=paste('mutect_only N=',nrow(passed_mutectonly), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_mutectonly, 'mutect_only', plot_snp_vaf=F)
dev.off()

pdf('diagnostics2.pdf', width=20, height=20)
par(mfrow=c(4,4))
barplot(
tri_context_counts[,'strelka_only_pass_rs_and_vum'], ylim=c(0,max_y+0.05), main=paste('strelka_only_pass_rs_and_vum N=',nrow(passed_strelkaonly_pass_rs_and_vum), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelkaonly_pass_rs_and_vum, 'mutect_only_pass_rs_and_vum', plot_snp_vaf=F)

barplot(
tri_context_counts[,'strelka_only_fail_rs_and_vum'], ylim=c(0,max_y+0.05), main=paste('strelka_only_fail_rs_and_vum N=',nrow(passed_strelkaonly_fail_rs_and_vum), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_strelkaonly_fail_rs_and_vum, 'mutect_only_pass_rs_and_vum', plot_snp_vaf=F)

barplot(
tri_context_counts[,'mutect_only_pass_rs_and_vum'], ylim=c(0,max_y+0.05), main=paste('mutect_only_pass_rs_and_vum N=',nrow(passed_mutectonly_pass_rs_and_vum), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_mutectonly_pass_rs_and_vum, 'mutect_only_pass_rs_and_vum', plot_snp_vaf=F)

barplot(
tri_context_counts[,'mutect_only_fail_rs_and_vum'], ylim=c(0,max_y+0.05), main=paste('mutect_only_fail_rs_and_vum N=',nrow(passed_mutectonly_fail_rs_and_vum), sep=""),
col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16))
)
axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
diagnostic_plots(passed_mutectonly_fail_rs_and_vum, 'mutect_only_fail_rs_and_vum', plot_snp_vaf=F)

dev.off()
##############
chr_lens = read.table("/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta.fai", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
temp = rownames(chr_lens)
chr_lens = chr_lens[,1]
names(chr_lens) = temp

change = paste(passed_twoplus$REF, ">", passed_twoplus$ALT, sep="")
passed_twoplus$CHANGE[change=="A>C"] = "T>G"
passed_twoplus$CHANGE[change=="A>G"] = "T>C"
passed_twoplus$CHANGE[change=="A>T"] = "T>A"
passed_twoplus$CHANGE[change=="G>A"] = "C>T"
passed_twoplus$CHANGE[change=="G>C"] = "C>G"
passed_twoplus$CHANGE[change=="G>T"] = "C>A"
passed_twoplus$with_rs = "yes"
no_rs = which(is.na(passed_twoplus$VEP_Existing_variation) & is.na(passed_twoplus$VEP_gnomAD_genome))
passed_twoplus$with_rs[no_rs] = "no"

pdf('diagnostics3.pdf', width=6, height=8)
par(mfrow=c(3,1))
counts = table(passed_twoplus$CHR)
barplot(counts / (chr_lens[names(counts)]/1000000), las=2, ylab="TMB (#muts per Mb)")
barplot(prop.table(table(passed_twoplus$with_rs, passed_twoplus$CHR),margin=2), las=2, ylab="% with RS id")
barplot(
	prop.table(table(passed_twoplus$CHANGE, passed_twoplus$CHR),margin=2), las=2,
	col=c('royalblue','black','red','grey','green2','hotpink'), ylab="% pyrimidine nucleotide change"
)
dev.off()
