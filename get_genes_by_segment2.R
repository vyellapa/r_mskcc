CGC_FILE = '/ifs/work/leukgen/home/gg10/ref_data/cancer_gene_census_v79_hg19.csv'

args=commandArgs(TRUE)
segment_file = args[1]
ploidy = as.numeric(args[2])
output_stub = args[3]

source('/ifs/res/papaemme/users/gg10/neuroblastoma_braf_mt/mskilab/do_neuroblastoma_manually_curated_junctions_functions.R')

segments = read.table(file=segment_file)
segments = rbind(segments[segments$V4>ploidy,], segments[segments$V4<2,])
gr.segments = GRanges(
	seqnames=Rle(sub('chr','', segments$V1)), IRanges(segments$V2, segments$V3),
	abs_cn=segments$V4, num_marks=segments$V5, logr=segments$V6
)
gt.ge = track.gencode(); ge = unlist(dat(gt.ge)[[1]]); ge[ge$type == 'gene']
hits = findOverlaps(ge, gr.segments)

write.table(as.data.frame(hits), file='hits.txt', row.names=F, col.names=F, quote=F)
write.table(as.data.frame(ge, row.names=1:length(ge)), file='ge.txt', row.names=T, col.names=F, quote=F)
write.table(as.data.frame(gr.segments, row.names=1:length(gr.segments)), file='junctions.txt', row.names=T, col.names=F, quote=F)
system('sort hits.txt > k')
system('mv k hits.txt')
system('sort ge.txt > k')
system('mv k ge.txt')
system('sort junctions.txt > k')
system('mv k junctions.txt')
system('join hits.txt ge.txt > joint1')
system("awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' joint1 | sort > k")
system('mv k joint1')
system('join joint1 junctions.txt > joint2')
out1 = paste(output_stub, '_annotations_all.txt.gz', sep='')
cmd = paste("awk '{print $1,$15,$16,$17,$18,$20,$21,$22,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' joint2 | gzip > ", out1, sep='')
system(cmd)
out2 = paste(output_stub, '_annotations_genes.txt.gz', sep='')
cmd = paste("awk '{if($11==\"gene\") print $1,$15,$16,$17,$18,$20,$21,$22,$3,$4,$5,$7,$8,$9,$10}' joint2 | sort -n | awk 'OFS=\"\\t\"{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | gzip > ", out2, sep="")
system(cmd)
system('rm ge.txt junctions.txt hits.txt joint1 joint2')

gene_annotations = read.table(out2)
cgc = read.csv(file=CGC_FILE)
m = merge(x=gene_annotations, y=cgc, by.x='V12', by.y='Gene.Symbol')
colnames(m)[1:14] = c('segment_chr','segment_start','segment_end','segment_length','segment_abscn','segment_nummark','segment_logr','gene_chr','gene_start','gene_end','gene_strand','gene_symbol','gene_name','ensg_id')
write.table(m, file=paste(output_stub, '_annotations_genes_cgc.txt.gz', sep=''), row.names=F, quote=F, sep="\t")

