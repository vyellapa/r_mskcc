#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-s", "--samples", default="none", type="character", nargs=1, help="Tab-delimited file the following header: id, sample, bam_file")
parser$add_argument("-a", "--lbkpt", default="none", type="character", nargs=1, help="Lower breakpoint window: chr:start:end:strand")
parser$add_argument("-b", "--hbkpt", default="none", type="character", nargs=1, help="Higher breakpoint window: chr:start:end:strand")
parser$add_argument("-d", "--outdir", type="character", nargs=1, help="Output directory. Should exist.")
parser$add_argument("-w", "--window", default=0, type="double", nargs=1, help="Will retrieve reads aligning to start-Window to end+Window.")
parser$add_argument("-o", "--outputstub", type="character", nargs=1, help="Output stub. Will create files for contigs (.fa), alignments (.bam) and breakpoints (_bkpt.txt)")

collate_reads <- function(sample_name, sample_bam, chr, start, end, strand) {
        if(strand1=='-'){
                cmd = paste('samtools view -f 0x10 -b ', sample_bam, ' "', chr, ':', start, '-', end, '" > tmp_', sample, '.bam', sep='')
        } else if (strand1=='+') {
                cmd = paste('samtools view -F 0x10 -b ', sample_bam, ' "', chr, ':', start, '-', end, '" > tmp_', sample, '.bam', sep='')
        } else {
                cmd = paste('samtools view -b ', sample_bam, ' "', chr, ':', start, '-', end, '" > tmp_', sample, '.bam', sep='')
        }
        system(cmd)
        cmd = paste(
                '/ifs/work/leukgen/opt/cgp/5.18.4/pcap/1.14.0/bin/bamtofastq exclude=QCFAIL,SECONDARY,SUPPLEMENTARY ',
                'tryoq=1 gz=1 level=1 outputperreadgroup=1 outputperreadgroupsuffixF=_i.fq outputperreadgroupsuffixF2=_i.fq ',
                'T=./bamtofastq.1 outputdir=. filename=tmp_', sample, '.bam split=20000000', sep=''
        )
        system(cmd)
	system('zcat *000000.gz >> all.reads.fastq')
}

# Parse arguments and complain
args <- parser$parse_args()
N_WINDOW = args$window
if(args$samples=="none") {
	print("Please provide a file sample info.")
	q(save="none")
} else {sample_file = args$samples}
if(args$lbkpt=="none") {
	print("At least one breakpoint must be provided.")
	q(save="none")
} else {
	lbkpt = unlist(strsplit(args$lbkpt, split=":"))
	chr1 = lbkpt[1]
	start1 = as.integer(lbkpt[2]) - N_WINDOW
	end1 = as.integer(lbkpt[3]) + N_WINDOW
	strand1 = lbkpt[4]
}
if(args$hbkpt=="none") {
	print("No higher breakpoint provided. Will collate reads from the lower breakpoint only.")
} else {
	hbkpt = unlist(strsplit(args$hbkpt, split=":"))
	chr2 = hbkpt[1]
	start2 = as.integer(hbkpt[2]) - N_WINDOW
	end2 = as.integer(hbkpt[3]) + N_WINDOW
	strand2 = hbkpt[4]
}
output_dir = args$outdir
output_prefix = args$outputstub

# Create temporary output directory
print(paste('Output dir = ', output_dir, sep=''))
gsa_root = tempfile(tmpdir=output_dir, pattern='tmp_sga_root_')
system(paste('mkdir ', gsa_root, sep=''))
samples = read.table(file=sample_file, header=T, stringsAsFactors=F)
setwd(gsa_root)

## Collect the reads from target regions
print('Collecting reads from all the samples...')
for(i in 1:nrow(samples)) {
	sample = samples[i,1]
	sample_bam = samples[i,3]
	collate_reads(sample, sample_bam, chr1, start1, end1, strand1)
	# If higher breakpoint is defined
	if(args$hbkpt!="none") {
		collate_reads(sample, sample_bam, chr2, start2, end2, strand2)	
	}
}
system('rm tmp*.bam')
system('mv all.reads.fastq reads.1.fastq')
system(paste('zless reads.1.fastq | awk ', "'{if($1~/@/ || $1~/+/ || length($1)==151) {print $0}}' > t", sep=''))
system('mv t reads.1.fastq ')

## Pre-process the reads
print('Preprocessing reads...')
system(paste('mkdir ', gsa_root, '/preprocessed', sep=''))
setwd(paste(gsa_root, '/preprocessed', sep=''))
cmd = 'sga preprocess -o reads.1.fastq.gz --pe-mode 2 ../reads.1.fastq'
print(cmd)
system(cmd)

## Index the raw data
system(paste('mkdir ', gsa_root, '/index_raw', sep=''))
setwd(paste(gsa_root, '/index_raw', sep=''))
system('ln -s ../preprocessed/reads.1.fastq.gz final.fastq.gz')
system('sga index --no-reverse -d 5000000 final.fastq.gz')

## Perform error correction
#print('Performing error correction...')
#cmd = paste('mkdir ', gsa_root, '/corrected', sep='')
#system(cmd)
#setwd(paste(gsa_root, '/corrected', sep=''))
#cmd = 'ln -s ../preprocessed/reads.1.fastq.gz'
#system(cmd)
#cmd = 'ln -s ../index_raw/final.bwt'
#system(cmd)
#cmd = 'ln -s ../index_raw/final.fastq.gz'
#system(cmd)
#cmd = 'sga correct -k 55 --learn -d 256 -t 8 -p final reads.1.fastq.gz -o reads.1.k55.ec.fa --metrics metrics_k55.txt > correct_k55.log'
#print(cmd)
#system(cmd)

## Filter reads
print('Filtering reads...')
system(paste('mkdir ', gsa_root, '/index_corrected', sep=''))
setwd(paste(gsa_root, '/index_corrected', sep=''))
#cmd = 'ln -s ../corrected/reads.1.k55.ec.fa'							# Soft-link to corrected reads file
system('gunzip ../preprocessed/reads.1.fastq.gz')
system('ln -s ../preprocessed/reads.1.fastq reads.1.k55.ec.fa')
cmd = 'sga index -d 5000000 -t 8 reads.1.k55.ec.fa'						# Index corrected reads
print(cmd)
system(cmd)
cmd = 'sga filter -d 256 -x 2 -t 8 reads.1.k55.ec.fa'						# Filter the reads to remove duplicates and low-quality reads remaining after correction
print(cmd)
system(cmd)
cmd = 'sga fm-merge -m 65 -t 8 reads.1.k55.ec.filter.pass.fa'					# Merge together reads that can be unambiguously assembled
print(cmd)
system(cmd)
cmd = 'sga index -d 2000000 -t 8 reads.1.k55.ec.filter.pass.merged.fa'				# Index the merged sequences
print(cmd)
system(cmd)

## Assemble contigs
print('Performing contig assembly...')
cmd = 'sga overlap -e 0.025 -m 65 -t 8 reads.1.k55.ec.filter.pass.merged.fa'			# Use overlap to construct the string graph of the merged reads
print(cmd)
system(cmd)
cmd = 'sga assemble -m 77 -d 0.4 -g 0.1 -r 10 -l 200 reads.1.k55.ec.filter.pass.merged.asqg.gz'	# Perform assembly
print(cmd)
system(cmd)

## Align contigs to hg19
print('Aligning contigs to the reference genome...')
cmd = 'bwa mem -M -T 0 -t 16 /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta default-contigs.fa | samtools view -Sb -o default-contigs.bam -'
system(cmd)
cmd = 'samtools sort -T default-contigs.tmp.sorted -o default-contigs.sorted.bam default-contigs.bam'
system(cmd)
cmd = 'mv default-contigs.sorted.bam default-contigs.bam'
system(cmd)
cmd = 'samtools index default-contigs.bam'
system(cmd)
cmd = 'Rscript /ifs/work/leukgen/home/gg10/soft/sga_denovo_assembly_parse_aln.R default-contigs.bam'
system(cmd)

cmd = paste('mv default-contigs.fa ', output_dir, '/', output_prefix, '.fa', sep='')
system(cmd)
cmd = paste('mv default-contigs.bam ', output_dir, '/', output_prefix, '.bam', sep='')
system(cmd)
cmd = paste('mv default-contigs.bam.bai ', output_dir, '/', output_prefix, '.bam.bai', sep='')
system(cmd)
cmd = paste('mv default-contigs_bkpt.txt ', output_dir, '/', output_prefix, '_bkpt.txt', sep='')
system(cmd)
system(paste('cd ', normalizePath(output_dir), sep=''))
system(paste('rm -r ', gsa_root, sep=''))
