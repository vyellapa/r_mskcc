library(GenomicRanges)

match_sample <- function(sample_name, sample_file, merged) {
	merged_lgr = with(merged, GRanges(chr1, IRanges(start=pmin(start1,end1)-10, end=pmax(start1,end1)+10), strand=strand1,id=id))
	sample1 = unique(read.table(sample_file, header=T, sep="\t", stringsAsFactors=F, comment=""))
	if(nrow(sample1) > 0 ) {
		colnames(sample1)[1] = "chr1"
		sample1_lgr = with(sample1, GRanges(chr1, IRanges(start=pmin(start1,end1)-10, end=pmax(start1,end1)+10), strand=strand1,id=id))
		lower_match = as.data.frame(cbind(as.data.frame(sample1_lgr[subjectHits(findOverlaps(merged_lgr,sample1_lgr)),])$id, as.data.frame(merged_lgr[queryHits(findOverlaps(merged_lgr,sample1_lgr)),])$id))
		colnames(lower_match) = c('sample_id','merged_id')
		lower_match$sample_id = as.character(lower_match$sample_id)
		lower_match$merged_id = as.character(lower_match$merged_id)
		matched = c()
		for(i in 1:nrow(lower_match)){
			sample = sample1[sample1$id==lower_match$sample_id[i],]
			merge =unique( merged[merged$id==lower_match$merged_id[i],])
			mstrand = "+"
			if(merge$strand2=="+") {mstrand = "-"}
			condition = sample$strand1==merge$strand1 & sample$strand2==mstrand & sample$chr2==merge$chr2 & (sample$start2-10)<merge$start2 & merge$start2<(sample$start2+10)
			if(length(condition)>1) {print(lower_match[i,])}
			if(condition) {
				matched = c(matched, TRUE)
			} else {
				matched = c(matched, FALSE)
			}
		}
		final_match = lower_match[matched,]
		sample1_subset = sample1[sample1$id %in% final_match$sample_id,]
		sample1_subset = sample1_subset[match(final_match$sample_id, sample1_subset$id),]

		merged_subset = merged[merged$id %in% final_match$merged_id,]
		merged_subset = merged_subset[match(final_match$merged_id, merged_subset$id),]
		sample1_subset = sample1_subset[,c('numcallers','callers','numreads_brass','numreads_svaba','numreads_gridss')]
		colnames(sample1_subset) = paste(sample_name,colnames(sample1_subset), sep="__")
		merged_subset = cbind(merged_subset, sample1_subset)

		rest_ids = setdiff(merged$id, merged_subset$id)
		rest = as.data.frame(matrix(data=NA, nrow=length(rest_ids), ncol=ncol(sample1_subset)))
		colnames(rest) = colnames(sample1_subset)
		final_merged = rbind(merged_subset, cbind(unique(merged[merged$id %in% rest_ids,]), rest))
	} else {
		final_merged = merged
	}
	return(final_merged)
}

args = commandArgs(trailingOnly=TRUE)

samples_file = args[1]
mbrass_file = args[2]
output_file = args[3]
merged = unique(read.delim(mbrass_file, header=T, sep="\t", stringsAsFactors=F))
conf = read.table(samples_file, header=T, sep="\t", stringsAsFactors=F)

for(i in 1:nrow(conf)) {
	merged = unique(match_sample(conf$sample[i], conf$bedpe[i], merged))
}
write.table(merged[,!grepl('reads$',colnames(merged))], file=output_file, row.names=F, sep="\t", quote=F)
