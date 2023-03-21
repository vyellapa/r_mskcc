#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-b", "--brassdir", type="character", nargs=1, help="Absolute path to the BRASS output directory. If not specified, will use LK to retrieve it.")
parser$add_argument("-t", "--tumor", type="character", nargs=1, help="Leukid for the tumor. Will use LK to retrieve files from pipeline.")
parser$add_argument("-n", "--normal", type="character", nargs=1, help="Leukid for the normal. Will use LK to retrieve files from pipeline.")
parser$add_argument("-m", "--minreads", default=10, type="double", nargs=1, help="Min number of aberrant reads. Will keep rearrangements supported with at least this many reads in both tumor and normal.")
parser$add_argument("-w", "--gwindow", default=10000, type="double", nargs=1, help="Will extend the genomic search space by this number of nucleotides.")
parser$add_argument("-r", "--regions", type="character", nargs=1, help="Target regions. A headerless tab-delimited file containing these columns: ENSE ENSG GENE CHR START END STRAND")
parser$add_argument("-o", "--outputstub", type="character", nargs=1, help="Output stub.")


# Parse arguments and complain
args <- parser$parse_args()
tumor.name = args$tumor
normal.name = args$normal
window = args$gwindow
regions.file = args$regions
brassdir = args$brassdir
outputstub = args$outputstub
if(is.null(brassdir)) {
	print("BRASS directory not provided. Will get the one from the pipeline.")
	brassdir = system(paste('lk get_outdirs -fi name BRASS -fi as_target__leukid ', tumor.name, sep=""), intern=T)
}
groups.file = list.files(path=brassdir, recursive=T, pattern='groups.gz', full.names=T)
if(is.null(groups.file) | length(groups.file)==0) {
	print('Could not file the groups file for that leukid.')
	quit(save='no')
}
if(is.null(tumor.name)) {
	print("Specify the leukid for the tumor.")
	quit(save='no')
}
if(is.null(normal.name)) {
	print("Specify the leukid for the normal.")
	quit(save='no')
}
if(is.null(regions.file)) {
	print("Specify the regions file.")
	quit(save='no')
}
if(is.null(outputstub)) {
	print("Specify the output stub.")
	quit(save='no')
}

VAGRENT.ENS.FILE = "/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/vagrent/Homo_sapiens.GRCh37.74.vagrent.cache.gz"
GENOME.REF = '/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta'

library(reshape2)
print(paste("Groups file:", groups.file, sep=""))
tumor.bam = system(paste("lk get_data -fi leukid ", tumor.name, sep=""), intern=TRUE)
normal.bam = system(paste("lk get_data -fi leukid ", tumor.name, sep=""), intern=TRUE)
regions = read.table(file=regions.file, header=F, sep="\t", stringsAsFactors=F)
assembled = NULL
all = NULL
for(gene in unique(regions$V3)) {
	# Retrieve rearrangements
	print(paste("Checking rearrangements for ", gene, sep=""))
	assembly.input = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	assembly.output = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	grass.input = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	grass.output = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	cmd = paste(
		"zless ", groups.file, " | awk 'OFS=\"\\t\"{if(",
			"($1==", min(regions[regions$V3==gene,'V4']), " && ",
			"$3>(", min(regions[regions$V3==gene,'V5'])-window, ") && ",
			"$3<(", max(regions[regions$V3==gene,'V6'])+window, ") && ",
			"$9>10 && $10>10) || ",
			"($5==", min(regions[regions$V3==gene,'V4']), " && ",
			"$7>(", min(regions[regions$V3==gene,'V5'])-window, ") && ",
			"$7<(", max(regions[regions$V3==gene,'V6'])+window, ") && ",
			"$9>10 && $10>10)", 
		") print $1,$3,$4,$5,$7,$8,$9,$10,$2,$6}' | awk 'OFS=\"\\t\"{print $1,$2,$3,$4,$5,$6,$7\"_\"$8,NR,$9,$10}'> ", assembly.input, sep=""
	)
	system(cmd)
	n1.rearr = as.integer(unlist(strsplit(system(paste("wc -l ", assembly.input, sep=""), intern=T), split=' '))[1])
	# Assemble rearrangements
	if(n1.rearr > 0) {
		all = rbind(all, cbind(gene,read.table(file=assembly.input, header=F, sep="\t")))
		print(paste("Attempting assembly for ", n1.rearr, " rearrangements for ", gene, sep=""))
		cmd = paste(
		"brass-assemble -X -m mem -O bedpe -r ", GENOME.REF, " -T . -o ", assembly.output, " ", assembly.input,
		" ", tumor.bam, ":", tumor.bam, ".bai ", normal.bam, ":", normal.bam, ".bai", sep=""
		)
		system(cmd)
		n2.rearr = as.integer(unlist(strsplit(system(paste("wc -l ", assembly.output, sep=""), intern=T), split=' '))[1])
		if(n2.rearr > 0) {
			print(paste("Assembled for ", n2.rearr, "/", n1.rearr, " rearrangements for ", gene, sep=""))
			# Annotate with Grass
			cmd = paste("sort -k1,1 -k 2,2n ", assembly.output, " > ", grass.input, sep="")
			system(cmd)
			cmd = paste(
				"grass.pl -genome_cache ", VAGRENT.ENS.FILE, " -ref ", GENOME.REF,
				" -species HUMAN -assembly GRCh37 -platform ILLUMINA-X10 -protocol WGS -tumour ",
				tumor.name, " -normal ", normal.name, " -file ", grass.input
			)
			system(cmd)
			assembled = rbind(assembled, cbind(gene,read.table(file=paste(grass.input, "_ann", sep=""), header=F, sep="\t")))
		} else {
			print(paste('No rearrangement assembled for ', gene, sep=""))
		}
		system(paste("rm ", assembly.output, " ", grass.input, "*", sep=""))
	} else {
		print(paste('No rearrangement found for ', gene, sep=""))
	}
	system(paste("rm ", assembly.input, sep=""))
}
if(!is.null(assembled) && nrow(assembled) > 0) {
	assembled = assembled[,1:35]
	assembled = cbind(colsplit(string=as.character(assembled$V7), pattern="_", names=c("normal.nreads", "tumor.nreads")), assembled)
	colnames(assembled) = unlist(strsplit("normal.nreads tumor.nreads query.gene chr1 start1 end1 chr2 start2 end2 id/name assembly_score strand1 strand2 samples assembly_notation non-template micro-homology assembled_readnames gene1 gene_id1 transcript_id1 strand1 end_phase1 region1 region_number1 total_region_count1 first/last1 gene2 gene_id2 transcript_id2 strand2 phase2 region2 region_number2 total_region_count2 first/last2 fusion_flag", split=" "))
	write.table(assembled, file=paste(outputstub, 'assembled.txt', sep=""), row.names=F, sep="\t", quote=F)

	# Plot rearrangements and copy number
	print(paste('BRASS=', brassdir))
	cov.files = list.files(path=brassdir, recursive=T, pattern='.ngscn.bed.gz', full.names=T)
	tumor.cov.file = grep(tumor.name, cov.files, value=T)
	normal.cov.file = grep(normal.name, cov.files, value=T)

	library(gTrack)
	library(gUtils)
	gt.ge = track.gencode(); ge = unlist(dat(gt.ge)[[1]]); ge[ge$type == 'gene']
	tumor.cov.sub.file = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	normal.cov.sub.file = tempfile(pattern = "germline_rearr_", tmpdir = '.')
	for(gene in assembled$query.gene) {
		win = ge[grep(gene, ge$gene_name)] + window
		cmd = paste(
		"zless ", tumor.cov.file, " | awk '",
		"$1==", as.character(seqnames(win))[1], " && ",
		"$3>", min(start(ranges(win))), " && ",
		"$3<", max(end(ranges(win))), "' > ", tumor.cov.sub.file, sep=""
		)
		print(cmd)
		system(cmd)
		cmd = paste(
		"zless ", normal.cov.file, " | awk '",
		"$1==", as.character(seqnames(win))[1], " && ",
		"$3>", min(start(ranges(win)))," && ",
		"$3<", max(end(ranges(win))), "' > ", normal.cov.sub.file, sep=""
		)
		print(cmd)
		system(cmd)
		t.data = read.table(file=tumor.cov.sub.file, header=F, sep="\t")
		n.data = read.table(file=normal.cov.sub.file, header=F, sep="\t")
		colnames(t.data) = c('seqnames','start','end','width','GC','depth')
		tdata.gr = with(t.data, GRanges(seqnames, IRanges(as.numeric(start), as.numeric(start), names=start), tumor_depth=depth))
		tdata.gt = gTrack(tdata.gr, 'tumor_depth', labels.suppress=T, circles=T, col='blue')
		colnames(n.data) = c('seqnames','start','end','width','GC','depth')
		ndata.gr = with(n.data, GRanges(seqnames, IRanges(as.numeric(start), as.numeric(start), names=start), normal_depth=depth))
		ndata.gt = gTrack(ndata.gr, 'normal_depth', labels.suppress=T, circles=T, col='orange')

		lower = assembled[assembled$query.gene==gene,c(4,5,12)]
		higher = assembled[assembled$query.gene==gene,c(7,8,13)]
		if(assembled[assembled$query.gene==gene,13]=='+') {higher[,3]='-'}
		if(assembled[assembled$query.gene==gene,13]=='-') {higher[,3]='+'}
		colnames(lower) = c('seqnames','start','strand')
		colnames(higher) = c('seqnames','start','strand')
		gr.data = rbind(lower, higher)
		gr = with(gr.data, GRanges(seqnames, IRanges(as.numeric(start), as.numeric(start), names=start), strand=strand))
		grl <- GRangesList('gr1' = gr)
		values(grl)$lty = 1
		values(grl)$lwd = 3
		values(grl)$col = 'purple'
		pdf(paste(outputstub, '_', gene, '.pdf', sep=""), height=10)
		plot(c(gt.ge, tdata.gt, ndata.gt), win, links=grl)
		dev.off()
	}
	system(paste("rm ", tumor.cov.sub.file, " ", tumor.cov.sub.file, sep=""))
} else {
        print("No genomic rearrangements have been found.")
}
print("DONE.")

