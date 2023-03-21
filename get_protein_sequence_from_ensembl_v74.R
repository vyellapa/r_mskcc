library(biomaRt)

ENSEMBL_v74_HOST = "dec2013.archive.ensembl.org"
ENSEMBL_v75_HOST = "feb2014.archive.ensembl.org"

args=(commandArgs(trailingOnly=TRUE))

ids = read.table(args[1], stringsAsFactors=F)
ofile = args[2]
host = args[3]

ensembl_host = ENSEMBL_v75_HOST
if(length(args)>=3) {
	version = args[3]
	if (version=="v74") {ensembl_host=ENSEMBL_v74_HOST}
}
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host=ensembl_host, path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
all = getBM(attributes=c(
		'hgnc_symbol','ccds','ensembl_transcript_id','ensembl_peptide_id','gene_biotype',
		'transcript_biotype','transcript_status','chromosome_name','transcript_start',
		'transcript_end'
), mart=ensembl)

ccds_ids = grep('CCDS', ids[,2], value=T)
enst_ids = grep('ENST', ids[,2], value=T)

ccds_data = NULL
if(length(ccds_ids)>0 & length(all[all$ccds %in% ccds_ids,'ensembl_transcript_id'])>0) {
	ccds_data = getSequence(id=all[all$ccds %in% ccds_ids,'ensembl_transcript_id'], type='ensembl_transcript_id', seqType='peptide', mart=ensembl)
	ccds_data = merge(x=all, y=ccds_data, by.x='ensembl_transcript_id', by.y='ensembl_transcript_id')
	ccds_data = merge(x=ids[ids$V2 %in% ccds_ids,], y=ccds_data, by.x='V2', by.y='ccds', all.x=T)
	colnames(ccds_data)[1:2] = c('query_id','var_id')
	ccds_data$ccds = ccds_data$query_id
} else if (length(ccds_ids)>0 & length(all[all$ccds %in% ccds_ids,'ensembl_transcript_id'])==0) {
	write.table('no CCDS', file=ofile, row.names=F, quote=F, sep="\t", col.names=F)
	quit(save="no", status=0)
}
enst_data = NULL
if(length(enst_ids)>0 & length(all[all$ensembl_transcript_id %in% enst_ids,'ensembl_transcript_id'])) {
	enst_data = getSequence(id=all[all$ensembl_transcript_id %in% enst_ids,'ensembl_transcript_id'], type='ensembl_transcript_id', seqType='peptide', mart=ensembl)
	enst_data = merge(x=all, y=enst_data, by.x='ensembl_transcript_id', by.y='ensembl_transcript_id')
	enst_data = merge(x=ids[ids$V2 %in% enst_ids,], y=enst_data, by.x='V2', by.y='ensembl_transcript_id', all.x=T)
	colnames(enst_data)[1:2] = c('query_id','var_id')
	enst_data$ensembl_transcript_id = enst_data$query_id
}
write.table(rbind(ccds_data, enst_data), file=ofile, row.names=F, quote=F, sep="\t")
