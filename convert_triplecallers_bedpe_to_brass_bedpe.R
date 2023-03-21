args = commandArgs(trailingOnly=TRUE)
input_bedpe = args[1]
output_bedpe = args[2]


bedpe = read.table(input_bedpe, header=F, sep="\t", stringsAsFactors=F)
#bedpe = bedpe[,(2:ncol(bedpe)-2)]

columns1 = unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id tumour_count strand1 strand2 sample rearr_type bkpt_distance brass_score spanning_read_names spanning_read_count col1 col2 col3 col4 col5 col6 brass_summary col7 col8 contributing_read_names contributing_read_number gene1 gene1_id gene1_transcript_id gene1_strand gene1_end_phase gene1_region gene1_region_number gene1_total_region_count gene1_first_last gene2 gene2_id gene2_transcript_id gene2_strand gene2_phase gene2_region gene2_region_number gene2_total_region_count gene2_first_last fusion_flag', split=" "))

columns2 = c('chr1','start1','end1','chr2','start2','end2','id','brass_score','strand1','strand2','sample','brass_summary','col7','col8','contributing_read_names','gene1','gene1_id','gene1_transcript_id','gene1_strand','gene1_end_phase','gene1_region','gene1_region_number','gene1_total_region_count','gene1_first_last','gene2','gene2_id','gene2_transcript_id','gene2_strand','gene2_phase','gene2_region','gene2_region_number','gene2_total_region_count','gene2_first_last','fusion_flag')

extras = as.data.frame(mat.or.vec(nrow(bedpe), length(setdiff(columns1, columns2))))
colnames(extras) = setdiff(columns1, columns2)
colnames(bedpe) = columns2
bedpe = cbind(bedpe, extras)
bedpe = bedpe[,columns1]
colnames(bedpe)[1] = "# chr1"
write.table(bedpe, file=output_bedpe, row.names=F, quote=F, sep="\t")
