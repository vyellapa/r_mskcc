
system("zless chrs | awk '{print \"grep -v Count_A tumour_allele_freqs/PD4021a_alleleFrequencies_chr\"$1\".txt | awk '$NF>10' >> tumour_allele_freqs.txt\"}' | bash")
system("zless chrs | awk '{print \"grep -v Count_A tumour_allele_freqs/PD4021a_alleleFrequencies_chr\"$1\".txt | awk '$NF>10' >> tumour_allele_freqs.txt\"}' | bash")
normal.data = read.table(file='normal_allele_freqs.txt', header=F, sep="\t")
tumour.data = read.table(file='tumour_allele_freqs.txt', header=F, sep="\t")
tumour.data$id = paste(tumour.data$V1, tumour.data$V2, sep="_")
normal.data$id = paste(normal.data$V1, normal.data$V2, sep="_")
m = merge(x=tumour.data, y=normal.data, by.x='id', by.y='id')
m = m[,c(2,3)]
colnames(m) = c('chr','pos')
rownames(m) = paste('snp',rownames(m), sep='')
write.table(m, file='snp_positions.txt', quote=F, sep="\t")
for(chr in as.character(unique(m$chr))) {
	write.table(m[m$chr==chr,], file=paste('chr',chr,'_snp_positions.txt',sep=''), quote=F, sep="\t")
}
