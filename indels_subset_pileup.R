args <- commandArgs(trailingOnly = TRUE)
pileup_file = args[1]
loci_file = args[2]
output_file = args[3]


tmp_file_name = tempfile(pattern='tmp_file_', tmpdir='.')
system(paste('grep -v "#"', pileup_file, '| awk \'$1!=""\' >', tmp_file_name))
pileup = read.table(tmp_file_name, header=T, sep="\t")
pileup$key = paste(pileup$Chrom, pileup$Pos, pileup$Ref, pileup$Alt, sep="_")
system(paste('rm', tmp_file_name))

loci = read.table(loci_file, header=F, stringsAsFactors=F)
loci$key = paste(loci$V1, loci$V2, loci$V3, loci$V4, sep="_")
pileup = merge(x=pileup, y=loci, by.x='key', by.y='key')
write.table(pileup, output_file, row.names=F, quote=F, sep="\t")
