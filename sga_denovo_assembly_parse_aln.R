args=(commandArgs(trailingOnly=TRUE))
bam_file = args[1]

matcher <- function(pattern, x) {
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}
doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}
## takes a cigar string and parses it, not very fast but...
cigarsums <- function(cigar, chars=c("M","N","D","I","S","H", "P", "X", "=")) {
   sapply (chars, doone, cigar)
}

cmd = paste0("samtools view ", bam_file, " | awk 'OFS=\"\\t\"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}'")
aln = as.data.frame(do.call(rbind, sapply(system(cmd, intern=T), FUN=strsplit, split="\t")))
rownames(aln) = NULL
colnames(aln)[1:11] = c('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL')
aln$CIGAR = as.character(aln$CIGAR)
cigar_split = t(sapply (aln$CIGAR, cigarsums))
rownames(cigar_split) = NULL
cigar_split = as.data.frame(cigar_split)
colnames(cigar_split) = paste0('CIGAR_',colnames(cigar_split))
cigar_split$N_M = 0
cigar_split$N_N = 0
cigar_split$N_D = 0
cigar_split$N_I = 0
cigar_split$N_S = 0
cigar_split$N_H = 0
cigar_split$N_P = 0
cigar_split$N_X = 0
for(char in c('M','N','D','I','S','H','P','X')) {
	for(i in 1:nrow(cigar_split)) {
		cigar_split[i,paste0('N_',char)] = length(grep(char,unlist(strsplit(aln$CIGAR[i], split=""))))
	}
}
selected = rowSums((cigar_split[,grep('^N_',colnames(cigar_split), value=T)]>1)+0)==0
if(length(which(selected))>0) {
	selected_aln = aln[selected,]
	selected_cigar_split = cigar_split[selected,]
	selected = rowSums((selected_cigar_split[,c('N_N','N_D','N_I','N_S','N_H','N_P','N_X')]>0)+0)>0
	if(length(which(selected))>0) {
		selected_aln = selected_aln[selected,]
		selected_cigar_split = selected_cigar_split[selected,c('CIGAR_M','CIGAR_S','CIGAR_H')]
		selected_cigar_split$IDX_M = 0
		selected_cigar_split$IDX_S = 0
		selected_cigar_split$IDX_H = 0
		for(char in c('M','S','H')) {
			for(i in 1:nrow(selected_cigar_split)) {
				idx = grep(char,unlist(strsplit(selected_aln$CIGAR[i], split="")))
				if(length(idx)==0) {idx=0}
				selected_cigar_split[i,paste0('IDX_',char)] = idx
			}
		}
		selected_aln = cbind(selected_aln, selected_cigar_split)
	
		selected_aln$CIGAR_FIRST = 0
		selected_aln$CIGAR_SECOND = 0
		for(i in 1:nrow(selected_aln)) {
			present = selected_aln[i,c('IDX_M','IDX_S','IDX_H')][,selected_aln[i,c('IDX_M','IDX_S','IDX_H')]>0]
			selected_aln[i,c('CIGAR_FIRST','CIGAR_SECOND')] = gsub('IDX_','',colnames(present[,order(present)]))
		}
		selected_aln$bkpt = 0
		selected_aln$POS = as.integer(as.character(selected_aln$POS))
		selected_aln$SEQ = as.character(selected_aln$SEQ)
		for(i in 1:nrow(selected_aln)) {
			bkpt = 0
			if(selected_aln$CIGAR_FIRST[i]=='M') {bkpt = selected_aln$POS[i] + selected_aln$CIGAR_M[i]}
			else if(selected_aln$CIGAR_FIRST[i] %in% c('H','S')) {bkpt = selected_aln$POS[i]}
			selected_aln$bkpt[i] = bkpt
		}
		write.table(selected_aln[,c('QNAME','FLAG','RNAME','POS','CIGAR','bkpt')], file=gsub('.bam','_bkpt.txt',bam_file), row.names=F, quote=F, sep="\t")
	} else {write.table('No chimeric contig.', file=gsub('.bam','_bkpt.txt',bam_file), row.names=F, col.names=F, quote=F, sep="\t")}
} else {write.table('No chimeric contig.', file=gsub('.bam','_bkpt.txt',bam_file), row.names=F, col.names=F, quote=F, sep="\t")}

