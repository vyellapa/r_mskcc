#
# grgsundemg@mskcc.org
# Brass for related samples
# args: Configuration file that contains info about all the tumour samples as well as the matched normal and output directory
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/multi_sample_brass2.R /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/multi_sample_brass/multi_sample_brass.conf /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/multi_sample_brass TRUE TRUE
#

library(GenomicRanges)

args=(commandArgs(trailingOnly=TRUE))
#args = unlist(strsplit('multi_sample_brass_with_PDX_samples/multi_sample_brass.conf multi_sample_brass_with_PDX_samples FALSE FALSE', split=' '))

conf_file = args[1]
output_dir = args[2]
preprocess_bams = TRUE
if(length(args)>=3) {
	if(args[3]=='FALSE' | args[3]=='F') { preprocess_bams = FALSE; print('Skipping the BAM file preprocessing step...'); }
}
regroup = TRUE
if(length(args)>=4) {
	if(args[4]=='FALSE' | args[4]=='F') { regroup = FALSE; print('Skipping the BAM file regrouping step...') }
}

MIN_NUMBER_READ_GR = 3
TMP_HEADER_FILE = paste(output_dir, '/h.sam', sep="")
TMP_MERGED_BAM_FILE = paste(output_dir, '/merged.bam', sep="")

reheader_bamfile <- function(bam_file, sample) {
	# Reheader the BRASS-remapped BAM file for tumour 1
	tmp_bam_file = paste(output_dir, '/', sample, '.bam', sep="")
	cmd1 = paste(
		"samtools view -H ", bam_file,
		" | /ifs/work/leukgen/home/gg10/soft/fix_pi_in_bam_header.py > ",
		TMP_HEADER_FILE, sep=""
	)
	cmd2 = paste('samtools reheader', TMP_HEADER_FILE, bam_file, '>', tmp_bam_file, sep=' ')
	cat('Re-headering the Brass-remapped BAM file for ', sample, '...\n')
	print(cmd1)
	system(cmd1)
	print(cmd2)
	system(cmd2)
}

conf = read.table(file=conf_file, header=T, stringsAsFactors=F)
conf$tmpbamfile = rep('',nrow(conf))
if(preprocess_bams) {
	print('Preprocessing BAM files...')
	for(i in 1:nrow(conf)) {
		reheader_bamfile(conf[i,'brmfile'], conf[i,'sample'])
		conf[i,'tmpbamfile'] = paste(output_dir, '/', conf[i,'sample'], '.bam', sep="")
	}
	# Merge BAM files
	cmd1 = paste(
		'java -jar /opt/common/CentOS_6/picard/picard-tools-1.96/MergeSamFiles.jar ',
		paste(paste('INPUT=',as.character(conf[conf$sampletype=='tumour','tmpbamfile']), sep=''), collapse=' '),
		' OUTPUT=', TMP_MERGED_BAM_FILE, sep=''
	)
	# Index the merged BAM file
	cmd2 = paste('samtools index ', TMP_MERGED_BAM_FILE, se='')
	cat('Merging the Brass-remapped BAM files...\n')
	print(cmd1)
	system(cmd1)
	cat('Indexing the merged BAM file...')
	print(cmd2)
	system(cmd2)
}
if(regroup) {
	print('Re-grouping the merge brm BAM file...')
	# Regroup reads in the merged BAM file
	conf$merged_groups_file = rep('',nrow(conf))
	for(i in 1:nrow(conf)) {
		merged_groups_file = paste(output_dir, '/', conf[i,'sample'], '_merged.groups', sep='')
		conf[i,'merged_groups_file'] = merged_groups_file
		cmd1 = paste(
			"bsub -We 180 -J bgroup -o ", conf[i,'sample'], 'bgroupd.log ',
			"'brass-group -I /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz ",
			'-n ', MIN_NUMBER_READ_GR,
			' -F /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/brass/brassRepeats.bed.gz ', TMP_MERGED_BAM_FILE,' ', conf[conf$sampletype=='normal','brmfile'],
			' | brassI_pre_filter.pl -i - -t ', conf[i,'bamid'], ' -o ', merged_groups_file,
			" -n /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/brass/brass_np_v4.groups.gz'", sep=""
		)
		cat('Running brass-regroup for the merged BAM file for sample #', i, '...\n')
		print(cmd1)
		system(cmd1)
	}
Sys.sleep(400)
}



for(i in 1:nrow(conf)) {
merged_groups_file = paste(output_dir, '/', conf[i,'sample'], '_merged.groups', sep='')
                conf[i,'merged_groups_file'] = merged_groups_file

}

#################################################
for(i in 1:nrow(conf)) {
	cat('Intersecting merged groups with sample-level results for sample # ',i, ': ', conf[i,'sample'], '\n')
	if(conf[i,'sampletype']=='tumour') {
	# COMPARE MERGED RESULTS WITH SAMPLE-LEVEL RESULTS
	# Read merged groups file
	# Read index information for the samples from the header
	cmd = paste("zless ", conf[i,'merged_groups_file'], " | awk '$1~/^#SAMPLE/' | sed 's/#//g'", sep="")
	samples.info = t(sapply(system(cmd, intern=TRUE), function(s) unlist(strsplit(s, split="\t"))[2:3]))
	# Retrieve the column index of each sample from the bedpe header
	conf$idx_sample = rep(-1, nrow(conf))
	for(k in 1:nrow(conf)) {
		conf[k,'idx_sample'] = as.integer(which(samples.info[,2]==conf[k,'bamid']))
	}
	# Read rearrangement calls
	merged.gr = read.table(as.character(conf[i,'merged_groups_file']), sep="\t", stringsAsFactors=F)
	# Filter calls that present in >=1 read in the matched normal
	# Keep only the columns for tumour read counts and read read names
	# The first 8 columns are lower and higher breakpoints.
	# V113 to ncol(Merged.gr) are for read names
	idx_normal = conf[conf$sampletype=='normal','idx_sample'] + 8
	idx_tumours = c(1:8, conf[conf$sampletype=='tumour','idx_sample']+8, conf[conf$sampletype=='tumour','idx_sample']+112)
	print(paste('Number of rearrangements before matched normal filtering: ', nrow(merged.gr), sep=''))
	merged.gr = merged.gr[merged.gr[,idx_normal]<1,idx_tumours]
	print(paste('Number of rearrangements after matched normal filtering: ', nrow(merged.gr), sep=''))
	colnames(merged.gr) = c(
		'chr1','strand1','start1','end1','chr2','strand2','start2','end2', c(
		paste(as.character(conf[conf$sampletype=='tumour','sample']),'.count', sep=''),
		paste(as.character(conf[conf$sampletype=='tumour','sample']),'.reads', sep='')
		)
	)
	merged.gr$id = rownames(merged.gr)
	# Create a genomic range object for merge groups
	merged.ranges = with(merged.gr, GRanges(chr1, IRanges(start=start1, end=end1), strand=strand1, id=id))
	# Read rearrangement calls for tumour1 into a genomic ranges object
	tumour.gr = read.table(as.character(conf[i,'bedfile']), sep="\t", stringsAsFactors=F)
	colnames(tumour.gr) = unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id tumour_count strand1 strand2 sample rearr_type bkpt_distance brass_score spanning_read_names spanning_read_count col1 col2 col3 col4 col5 col6 brass_summary col7 col8 contributing_read_names contributing_read_number gene1 gene1_id gene1_transcript_id gene1_strand gene1_end_phase gene1_region gene1_region_number gene1_total_region_count gene1_first_last gene2 gene2_id gene2_transcript_id gene2_strand gene2_phase gene2_region gene2_region_number gene2_total_region_count gene2_first_last fusion_flag', split=" "))
	tumour.ranges = with(tumour.gr, GRanges(chr1, IRanges(start=pmin(start1,end1)-250, end=pmax(start1,end1)+250), strand=strand1,id=id))
	print(paste('Number of rearrangements after called in tumour_vs_normal BRASS run: ', nrow(tumour.gr), sep=''))
	# Intersect tumour ranges with merged group ranges by lower end
	idx.merged.calls = as.data.frame(merged.ranges[subjectHits(findOverlaps(tumour.ranges,merged.ranges)),])$id
	idx.query.hits = as.data.frame(tumour.ranges[queryHits(findOverlaps(tumour.ranges,merged.ranges)),])$id
	merged.calls = cbind(idx.query.hits, merged.gr[idx.merged.calls,])
	# Intersect by higher end
	match = rep(FALSE,nrow(merged.calls))
	flipped_strand = rep('-', nrow(merged.calls))
	flipped_strand[merged.calls$strand2=='-'] = '+'
	for(k in 1:length(match)) {
		match[k]= (
		merged.calls[k,'chr2']==tumour.gr[tumour.gr$id==merged.calls[k,'idx.query.hits'],'chr2'] &
		
		(tumour.gr[tumour.gr$id==merged.calls[k,'idx.query.hits'],'start2']-500) < merged.calls[k,'start2'] &
		(tumour.gr[tumour.gr$id==merged.calls[k,'idx.query.hits'],'start2']+500) > merged.calls[k,'end2'] &
		(tumour.gr[tumour.gr$id==merged.calls[k,'idx.query.hits'],'strand2'] == flipped_strand[k])
		)
	}
	merged.calls = merged.calls[match,]
	# Grouped BRASS file might have >1 junctions very close to each other i.e. within 500 bp of each other
	# Merge such calls in so as to maximize the total number of read support across the tumour samples
	idx.to.keep = c()
	for(m in merged.calls$idx.query.hits[duplicated(merged.calls$idx.query.hits)]) {
		merged.subset = merged.calls[merged.calls$idx.query.hits==m,]
		max.value = max(rowSums(merged.subset[,grep('count$',colnames(merged.calls), value=T)]))
		max.idx = merged.subset[rowSums(merged.subset[,grep('count$',colnames(merged.calls), value=T)])==max.value,'id']
		idx.to.keep = c(idx.to.keep, max.idx)
	}
	merged.calls = rbind(
		merged.calls[!(merged.calls$idx.query.hits %in% merged.calls$idx.query.hits[duplicated(merged.calls$idx.query.hits)]),],
		merged.calls[merged.calls$id %in% idx.to.keep,]
	)
	colnames(merged.calls)[ncol(merged.calls)] = 'mbrass.id'
	print(paste('Number of rearrangements remaining after filtering by single-sample BRASS run: ', nrow(merged.calls), sep=''))
	m1 = merge(x=tumour.gr, y=merged.calls, by.x='id', by.y='idx.query.hits')
	new.str = m1$strand2.y
	new.str[m1$strand2.y=="-"] = "+"
	new.str[m1$strand2.y=="+"] = "-"
	m1 = m1[m1$strand1.x==m1$strand1.y & m1$strand2.x==new.str,]
	# Retreive the rearrangements called in single-sample BRASS analysis but not called in multi-sample BRASS analysis
	notcalled = tumour.gr[!(tumour.gr$id %in% m1$id),]
	notcalled = cbind(notcalled[,'id'], notcalled[,c(1:6,8:ncol(notcalled))])
	colnames(notcalled)[1] = 'id'
	notcalled2 = mat.or.vec(nrow(notcalled),ncol(merged.calls)-1)
	colnames(notcalled2) = colnames(merged.calls)[colnames(merged.calls)!="idx.query.hits"]
	tumour_out_file = paste(output_dir, '/', conf[i,'sample'], '_final_groups.bedpe', sep='')
	write.table(m1, file=tumour_out_file, row.names=F, quote=F, sep="\t")
	write.table(cbind(notcalled, notcalled2), file=tumour_out_file, row.names=F, quote=F, sep="\t", append=T, col.names=F)

	}
}
# Cross-ref results from different multi-sample re-group runs to get a unique set of rearrangements across the samples
all.files = list.files(path=output_dir, pattern='*final_groups.bedpe', full.names=T)
all.output = list()
for(file in all.files) {all.output[[unlist(strsplit(basename(file), split='_final'))[1]]] = read.table(file, sep="\t", header=T, stringsAsFactors=F)}
all.output = do.call('rbind', all.output)
# Unique rearrangements using the lower breakpoint
uniq.data = NULL
max_distance = 400
for(chr in unique(all.output$chr1.x)) {
	all.subset = all.output[all.output$chr1.x==chr,]
#	all.subset = all.output[all.output$chr1.x==chr & all.output$start1.x>10958160 & all.output$start1.x<15808734,]
	if(nrow(all.subset)>1) {
		all.subset = all.subset[order(all.subset$start1.x, all.subset$chr2.x, all.subset$start2.x),]
		rownames(all.subset) = 1:nrow(all.subset)
		all.subset$distance = c(all.subset$start1.x[2:(nrow(all.subset))] - all.subset$start1.x[1:(nrow(all.subset)-1)], -1)
		n_to_unique = length(which(all.subset$distance<max_distance & all.subset$distance!=(-1)))
		while(n_to_unique>0) {
			all.subset$distance = c(all.subset$start1.x[2:(nrow(all.subset))] - all.subset$start1.x[1:(nrow(all.subset)-1)], -1)
			idx.to.keep = c()
			idx.uniqued = c()
			for (i in which(all.subset$distance<max_distance & all.subset$distance!=(-1))) {
				if(
					all.subset[i,'strand1.x'] == all.subset[i+1,'strand1.x'] & all.subset[i,'strand2.x'] == all.subset[i+1,'strand2.x'] &
					abs(all.subset[i+1,'start2.x'] - all.subset[i,'start2.x'])<max_distance & all.subset[i,'chr2.x'] == all.subset[i+1,'chr2.x']
				) {
					idx.to.keep = c(idx.to.keep, i)
					idx.uniqued = c(idx.uniqued, i, i+1)
				}
			}
			if(!is.null(idx.to.keep)) {
				all.subset = rbind(all.subset[idx.to.keep,], all.subset[setdiff(rownames(all.subset), idx.uniqued),])
			} else {
				break
			}
			all.subset = all.subset[order(all.subset$start1.x, all.subset$chr2.x, all.subset$start2.x),]
			rownames(all.subset) = 1:nrow(all.subset)
			all.subset$distance = c(all.subset$start1.x[2:(nrow(all.subset))] - all.subset$start1.x[1:(nrow(all.subset)-1)], -1)
			n_to_unique = length(which(all.subset$distance<100 & all.subset$distance!=(-1)))
		}
		uniq.data = rbind(uniq.data, all.subset)
	} else {
		all.subset$distance = -1
		uniq.data = rbind(uniq.data, all.subset)
	}
}
# Unique rearrangements using the higher breakpoint
uniq.data2 = NULL
for(chr in unique(uniq.data$chr2.x)) {
	all.subset = uniq.data[uniq.data$chr2.x==chr,]
	if(nrow(all.subset)>1) {
		all.subset = all.subset[order(all.subset$start2.x, all.subset$chr1.x, all.subset$start1.x),]
		rownames(all.subset) = 1:nrow(all.subset)
		all.subset$distance = c(all.subset$start2.x[2:(nrow(all.subset))] - all.subset$start2.x[1:(nrow(all.subset)-1)], -1)
		n_to_unique = length(which(all.subset$distance<max_distance & all.subset$distance!=(-1)))
		while(n_to_unique>0) {
			all.subset$distance = c(all.subset$start2.x[2:(nrow(all.subset))] - all.subset$start2.x[1:(nrow(all.subset)-1)], -1)
			idx.to.keep = c()
			idx.uniqued = c()
			for (i in which(all.subset$distance<max_distance & all.subset$distance!=(-1))) {
				if(
					all.subset[i,'strand2.x'] == all.subset[i+1,'strand2.x'] & all.subset[i,'strand1.x'] == all.subset[i+1,'strand1.x'] &
					abs(all.subset[i+1,'start1.x'] - all.subset[i,'start1.x'])<max_distance & all.subset[i,'chr1.x'] == all.subset[i+1,'chr1.x']
				) {
					idx.to.keep = c(idx.to.keep, i)
					idx.uniqued = c(idx.uniqued, i, i+1)
				}
			}
			if(!is.null(idx.to.keep)) {
				all.subset = rbind(all.subset[idx.to.keep,], all.subset[setdiff(rownames(all.subset), idx.uniqued),])
			} else {
				break
			}
			all.subset = all.subset[order(all.subset$start2.x, all.subset$chr1.x, all.subset$start1.x),]
			rownames(all.subset) = 1:nrow(all.subset)
			all.subset$distance = c(all.subset$start2.x[2:(nrow(all.subset))] - all.subset$start2.x[1:(nrow(all.subset)-1)], -1)
			n_to_unique = length(which(all.subset$distance<100 & all.subset$distance!=(-1)))
		}
		uniq.data2 = rbind(uniq.data2, all.subset)
	} else {
	all.subset$distance = -1
	uniq.data2 = rbind(uniq.data2, all.subset)
	}
}
uniq.data2 = uniq.data2[order(uniq.data2$chr1.x, uniq.data2$start1.x, uniq.data2$chr2.x, uniq.data2$start2.x),]
write.table(uniq.data2, paste(output_dir, '/final_consolidated_rearrangements.txt', sep=""), sep="\t")
# Clean up
#system(paste('rm ',  TMP_HEADER_FILE, ' ', output_dir, '/*.bam*', sep=''))
cat('DONE.')
