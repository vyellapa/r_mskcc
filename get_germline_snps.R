#
# Interrogate the caveman SNP output for possible germline mutations that fall within the regions of interest
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/get_germline_snps.R Adrenal-162414T/against_blood_normal/caveman/W0006075F_vs_W0006082F.snps.ids.vcf.gz /ifs/work/leukgen/home/gg10/ref_data/icr_gene_panel_exons.txt outfile.txt
#


args <- commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
exons_file = args[2]
out_file = args[3]

TMP.VCF.FILE = tempfile(pattern = "tmp.", tmpdir = '.')
TMP.VCF.ANNOT.FILE = tempfile(pattern = "tmp.", tmpdir = '.')
TMP.EXAC.FILE = tempfile(pattern = "tmp.", tmpdir = '.')

# Retrieve relevant SNPs using tabix
exons = read.table(exons_file, header=F, sep="\t", stringsAsFactors=F)
cmd = paste('zgrep "#" ', vcf_file, ' > ', TMP.VCF.FILE, sep='')
system(cmd)
for(i in 1:nrow(exons)) {
	cat('Region ',i,'\n')
	cmd = paste('tabix ', vcf_file, ' ', exons[i,'V4'], ':', as.integer(exons[i,'V5'])-20, '-', as.integer(exons[i,'V6'])+20, ' >> ', TMP.VCF.FILE, sep='')
	system(cmd)
}

# Annotate the SNPs using Vagrent
if(as.integer(system(paste('less ', TMP.VCF.FILE, ' | grep -v "#" | wc -l', sep=''), intern=TRUE))>0) {
	cmd = paste('AnnotateVcf.pl -i ', TMP.VCF.FILE, ' -o ', TMP.VCF.ANNOT.FILE, ' -c /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz -sp Human -as GRCh37 -t', sep='')
	system(cmd)

	snps = read.table(file=paste(TMP.VCF.ANNOT.FILE, '.gz', sep=''), comment.char="#", header=F, colClasses=c("character"))
	for(i in 1:nrow(snps)) {
		cmd = paste('tabix /ifs/work/leukgen/ref/homo_sapiens/37/exac_broad/ExAC.r0.3.sites.vep.vcf.gz ', paste(snps[i,'V1'], ':', snps[i,'V2'], '-', snps[i,'V2'], sep=''), ' >> ', TMP.EXAC.FILE, sep='')
	system(cmd)
	}
	if(as.integer(system(paste('less ', TMP.EXAC.FILE, ' | grep -v "#" | wc -l', sep=''), intern=TRUE))>0) {
		exac = read.table(file=TMP.EXAC.FILE, header=F, sep="\t")
		exac$id = paste(exac$V1, exac$V2, sep="_")
		snps$id = paste(snps$V1, snps$V2, sep="_")
		m = merge(x=snps, y=exac, by.x='id', by.y='id', all.x=T)
		m = m[,2:20]
	} else {
		m = cbind(snps, t(rep("NONE",8)))	
	}
	colnames(m) = c('chr','pos','var.id','wt','mt','qual','flag','info','format','normal','tumour','exac.chr','exac.pos','exac.id','exac.wt','exac.mt','exac.qual','exac.flag','exac.info')
	write.table(m, out_file, row.names=F, sep="\t", quote=F)
	system(paste('rm ', TMP.EXAC.FILE, sep=""))
	system(paste('rm ', TMP.VCF.FILE, ' ', TMP.VCF.ANNOT.FILE, '*', sep=''))
} else {
	print('No SNPs found.')
	system(paste('rm ', TMP.VCF.FILE, sep=''))
}
