# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/intersect_cn_rearr_merge_outputs.R /ifs/res/papaemme/users/gg10/mycn_amplified_nbl/data/I-H-112908-T1-1-D1-1/brass/extras /ifs/res/papaemme/users/gg10/mycn_amplified_nbl/data/I-H-112908-T1-1-D1-1/brass/pipeline_results/I-H-112908-T1-1-D1-1_vs_I-H-112908-N1-1-D1-1.annot.bedpe

args = commandArgs(trailingOnly=TRUE)

library(reshape2)

td_col         = rgb(238/255, 118/255, 0) # darkorange2
del_col        = rgb(0, 178/255, 238/255) # deepskyblue2
inter_chrs_col = rgb(191/255, 62/255, 255/255) #darkorchid1
tail_tail_col  = rgb(205/255, 205/255, 0/255) # yellow3
head_head_col  = rgb(85/255, 107/255, 47/255) # darkolivegreen

output_dir = args[1]
first_assembly_file = args[2]
REARRS_WITH_SEGS_FILE = paste(output_dir, "/rearrs_with_cnsegs.txt", sep="")

print(output_dir)
second_assembly = NULL
for (file in list.files(path=output_dir, pattern='assembly_output*', full.names=T)) {
	if(file.info(file)$size>0) {
		second_assembly = rbind(second_assembly, read.table(file, sep="\t", header=F))
	}
}
colnames(second_assembly) = unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id_name assembly_score strand1 strand2 sample Brass_Notation non-template micro-homology assembled_read_names', split=' '))


all_rearrs = read.table(REARRS_WITH_SEGS_FILE)
colnames(all_rearrs) = c('cnaseg_chr','cnaseg_start','cnaseg_end','cnaseg_ntot','cnaseg_nmarks','chrl','startl','endl','chrh','starth','endh','id_name','brass_score','strl','strh')
all_rearrs$is_assembled = rep('notassembled',nrow(all_rearrs))
all_rearrs[all_rearrs$id_name %in% second_assembly$id_name,'is_assembled'] = 'assembled'
#all_rearrs = all_rearrs[order(all_rearrs$id_name),]

merged = merge(x=second_assembly, y=all_rearrs, by.x='id_name', by.y='id_name', all.y=T)
merged$strand1 = merged$strl
merged[merged$strh=='-','strand2'] = '+'
merged[merged$strh=='+','strand2'] = '-'
assembeled = merged[merged$is_assembled=='assembled',unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id_name brass_score strand1 strand2 sample Brass_Notation non-template micro-homology assembled_read_names', split=' '))]
notassembeled = merged[merged$is_assembled=='notassembled',unlist(strsplit('chrl startl endl chrh starth endh id_name brass_score strand1 strand2 sample Brass_Notation non-template micro-homology assembled_read_names', split=' '))]
notassembeled$chr1 = notassembeled$chrl
notassembeled$chr2 = notassembeled$chrh
notassembeled$start1 = notassembeled$startl + round((notassembeled$endl - notassembeled$startl)/2)
notassembeled$start2 = notassembeled$starth + round((notassembeled$endh - notassembeled$starth)/2)
notassembeled$end1 = notassembeled$start1
notassembeled$end2 = notassembeled$start2
all_cna = rbind(assembeled, notassembeled[,colnames(assembeled)])
all_cna$assembly_score = colsplit((colsplit(string=all_cna$Brass_Notation, pattern="score ", names=c("Part1", "Part2"))$Part2), pattern="\\)", names=c("Part1","Part2"))$Part1
a = sapply(as.character(all_cna$assembled_read_names), function(x) length(unlist(strsplit(x, ","))))
names(a) = NULL
all_cna$assembled_read_count = a
all_cna$bkdist = all_cna$end2 - all_cna$start1
all_cna[all_cna$chr1 != all_cna$chr2,'bkdist'] = -1
all_cna$svclass = rep('translocation',nrow(all_cna))
all_cna[all_cna$strand1=='+' & all_cna$strand2=='+','svclass'] = 'deletion'
all_cna[all_cna$strand1=='-' & all_cna$strand2=='-','svclass'] = 'tandom-duplication'
all_cna[all_cna$strand1=='-' & all_cna$strand2=='+','svclass'] = 'inversion'
all_cna[all_cna$strand1=='+' & all_cna$strand2=='-','svclass'] = 'inversion'
all_cna$copynumber_flag = rep(1,nrow(all_cna))
all_cna = unique(all_cna)
write.table(all_cna, file=REARRS_WITH_SEGS_FILE, row.names=F, col.names=T, sep='\t', quote=F)

first_assembly = read.table(first_assembly_file, comment.char='#', header=F, sep="\t", stringsAsFactors=F)
colnames(first_assembly) = unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id_name brass_score strand1 strand2 sample svclass bkdist assembly_score readpair_names readpair_count bal_trans inv occL occH copynumber_flag range_blat Brass_Notation non-template micro-homology assembled_read_names assembled_read_count gene1 gene_id1 transcript_id1 grass.strand1 end_phase1 region1 region_number1 total_region_count1 first_last1 gene2 gene_id2 transcript_id2 grass.strand2 phase2 region2 region_number2 total_region_count2 first/last2 fusion_flag', split=' '))



first_assembly$copynumber_flag[first_assembly$id_name %in% all_cna$id_name] = 1
extra_rearrs = unique(all_cna[all_cna$id_name %in% setdiff(all_cna$id_name, first_assembly$id_name),])
missing.cols = setdiff(colnames(first_assembly), colnames(extra_rearrs))
missing.data = as.data.frame(mat.or.vec(nrow(extra_rearrs), length(missing.cols)))
colnames(missing.data) = missing.cols
extra_rearrs = cbind(extra_rearrs, missing.data)
all = rbind(first_assembly, extra_rearrs[,colnames(first_assembly)])
all$confidence = rep('NA', nrow(all))
all[all$id_name %in% first_assembly$id_name,'confidence'] = 'pipeline'
all[!is.na(all$assembly_score) & !(all$id_name %in% first_assembly$id_name),'confidence'] = 'assembled'
all[is.na(all$assembly_score) & !(all$id_name %in% first_assembly$id_name),'confidence'] = 'notassembled'

all = all[order(all$chr1, all$start1, all$chr2, all$start2),]
all$colour = rep('black',nrow(all))
all[all$id_name %in% first_assembly$id_name & all$strand1=='+' & all$strand2=='+','colour'] = del_col
all[all$id_name %in% first_assembly$id_name & all$strand1=='-' & all$strand2=='-','colour'] = td_col
all[all$id_name %in% first_assembly$id_name & all$chr1!=all$chr2,'colour'] = inter_chrs_col
all[all$id_name %in% first_assembly$id_name & all$strand1=='-' & all$strand2=='+','colour'] = tail_tail_col
all[all$id_name %in% first_assembly$id_name & all$strand1=='+' & all$strand2=='-','colour'] = head_head_col
write.table(all, file=paste(output_dir, '/final_output.bedpe', sep=''), row.names=F, col.names=T, quote=F, sep="\t")
