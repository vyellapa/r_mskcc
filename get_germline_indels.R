#
# Interrogate the Pindel output for possible germline mutations that fall within the regions of interest
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/get_germline_indels.R Adrenal-162414T/against_blood_normal/pindel/W0006075F_vs_W0006082F.flagged.annot.vcf.gz /ifs/work/leukgen/home/gg10/ref_data/icr_gene_panel_exons.txt germline_mutation_check/Adrenal-162414T_icrgene_germline_indels.txt
#


args <- commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
exons_file = args[2]
out_file = args[3]

TMP.VCF.FILE = tempfile(pattern = "tmp.", tmpdir = '.')
TMP.VCF.ANNOT.FILE = tempfile(pattern = "tmp.", tmpdir = '.')
TMP.EXAC.FILE = tempfile(pattern = "tmp.", tmpdir = '.')

# Retrieve relevant indels using tabix
exons = read.table(exons_file, header=F, sep="\t", stringsAsFactors=F)
cmd = paste('zgrep "#" ', vcf_file, ' > ', TMP.VCF.FILE, sep='')
system(cmd)
for(i in 1:nrow(exons)) {
	cat('Region ',i,'\n')
	cmd = paste('tabix ', vcf_file, ' ', exons[i,'V4'], ':', as.integer(exons[i,'V5'])-20, '-', as.integer(exons[i,'V6'])+20, ' >> ', TMP.VCF.FILE, sep='')
	system(cmd)
}
if(as.integer(system(paste('less ', TMP.VCF.FILE, ' | grep -v "#" | wc -l', sep=''), intern=TRUE))>0) {
	cmd = paste('AnnotateVcf.pl -i ', TMP.VCF.FILE, ' -o ', TMP.VCF.ANNOT.FILE, ' -c /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz -sp Human -as GRCh37 -t', sep='')
	print(cmd)
	system(cmd)

	snps = read.table(file=paste(TMP.VCF.ANNOT.FILE, '.gz', sep=''), comment.char="#", header=F, colClasses=c("character"))
	for(i in 1:nrow(snps)) {
		cmd = paste('tabix /ifs/work/leukgen/ref/homo_sapiens/37/exac_broad/ExAC.r0.3.sites.vep.vcf.gz ', paste(snps[i,'V1'], ':', snps[i,'V2'], '-', snps[i,'V2'], sep=''), ' >> ', TMP.EXAC.FILE, sep='')
		print(cmd)
		system(cmd)
	}
	exac = read.table(file=TMP.EXAC.FILE, header=F, sep="\t")
	exac$id = paste(exac$V1, exac$V2, sep="_")
	snps$id = paste(snps$V1, snps$V2, sep="_")
	m = merge(x=snps, y=exac, by.x='id', by.y='id', all.x=T)
	colnames(m) = c('id','chr','pos','var.id','wt','mt','qual','flag','info','format','normal','tumour','exac.chr','exac.pos','exac.id','exac.wt','exac.mt','exac.qual','exac.flag','exac.info')
	write.table(m[,2:20], out_file, row.names=F, sep="\t", quote=F)

	system(paste('rm ', TMP.VCF.FILE, sep=''))
	system(paste('rm ', TMP.EXAC.FILE, sep=""))
	system(paste('rm ', TMP.VCF.FILE, ' ', TMP.VCF.ANNOT.FILE, '*', sep=''))
} else {
	print('No indels found.')
	system(paste('rm ', TMP.VCF.FILE, sep=''))
}
