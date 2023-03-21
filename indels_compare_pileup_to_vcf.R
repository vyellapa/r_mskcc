args <- commandArgs(trailingOnly = TRUE)
pileup_file = args[1]
vcf_file = args[2]
output_file = args[3]

tmp_file_name = tempfile(pattern='tmp_file_', tmpdir='.')
system(paste('zless ', vcf_file, ' | grep -v "#" | awk \'$1!=""\' >', tmp_file_name))
vcf = read.delim(tmp_file_name, header=T, sep="\t")
system(paste('rm', tmp_file_name))
vcf$key = paste(vcf$CHR, vcf$START, vcf$REF, vcf$ALT, sep="_")

pileup = read.table(pileup_file, header=T, sep="\t", stringsAsFactors=T)

vcf_subset = merge(x=vcf, y=pileup, by.x='key', by.y='key')
tumor_vaf_column = paste(gsub('-','.',unlist(strsplit(basename(vcf_file), split="_vs_"))[1]), '_VAF', sep="")
normal_vaf_column = grep('\\.N',grep('_VAF',colnames(vcf_subset),value=T),value=T)
vcf_subset$LENGTH = vcf_subset$END - vcf_subset$START
vcf_subset$color = "black"
vcf_subset$color[vcf_subset[,normal_vaf_column]>0.02] = "red"
prop = length(which(vcf_subset$color=='red')) / nrow(vcf_subset)
pdf(output_file, width=12, height=12)
par(mfrow=c(3,3))
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'mutect_TARGET_VAF'], col=vcf_subset$color, pch=16,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Mutect VAF', main=paste('prop. >0.01 in normal ', prop, sep="")
)
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'strelka_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Strelka VAF', pch=16
)
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'pindel_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Pindel VAF', pch=16
)
vcf_subset$color = "black"
vcf_subset$color[vcf_subset[,'VAG_VT']=="Del"] = "brown"
vcf_subset$color[vcf_subset[,'VAG_VT']=="Ins"] = "blue2"
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'mutect_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Mutect VAF', pch=16
)
legend('topleft', legend=c('Complex','Del','Ins'), fill=c('black','brown','blue2'))
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'strelka_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Strelka VAF', pch=16
)
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'pindel_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Pindel VAF', pch=16
)
vcf_subset$color = "black"
vcf_subset$color[vcf_subset$LENGTH>1 & vcf_subset$LENGTH<4] = "brown"
vcf_subset$color[vcf_subset$LENGTH>=4] = "blue2"
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'mutect_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Mutect VAF', pch=16
)
legend('topleft', legend=c('1bp','1-3bp','>=4bp'), fill=c('black','brown','blue2'))
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'strelka_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Strelka VAF', pch=16
)
plot(
	x=vcf_subset[,tumor_vaf_column], y=vcf_subset[,'pindel_TARGET_VAF'], col=vcf_subset$color,
	xlim=c(0,1), ylim=c(0,1), xlab="Corrected VAF", ylab='Pindel VAF', pch=16
)


dev.off()
