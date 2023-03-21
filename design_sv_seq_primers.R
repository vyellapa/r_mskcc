########################## inter-pretation of the genomic strands and the corresponding VCF annotation
#	GENOMIC_STR_LOW	GENOMIC_STR_HIGH	REARR_TYPE	VCF_REF	VCF_ALT	MEANING
#	+		+			deletion	s	t[p[	piece extending to the right of p is joined after t
#	+		-			inversion	s	t]p]	reverse compliment of piece extending to the left of p is joined after t
#	-		-			tandem_dup	s	]p]t	piece extending to the left of p is joined before t
#	-		+			inversion	s	[p[t	reverse compliment of piece extending to the right of p is joined before t
##### Note that Primer3 accepts 5'>3' orientation (positive genomic strand).


library(Rsamtools)
REF_FASTA = "/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta"
PRIMER3_ROOT = "/ifs/work/leukgen/home/gg10/soft/primer3-2.3.7/src"
DISTANCE_FROM_BKPT = 20

pick_genomic_region <- function(chr1, strand1, start1, chr2, strand2, start2, n1=200, n2=200) {
	if(strand1=="-" & strand2=="+") {
	# lower == - and higher == + ==> reverse compliment of the genomic segment at position > start2 is joined to the genomic segment at position < start1
	# Resulting mutant contig: (start2+n1)-----(start2) | (start1)-----(start1+n2)
		res = rbind(c(chr1, start1, start1+n1, 'lower', 'no',2), c(chr2, start2, start2+n2, 'higher', 'yes', 1))
	} else if (strand1=="-" & strand2=="-") {
	# lower == - and higher == - ==> the genomic segment at position < start2 is joined to the genomic segment at position < start1
	# Resulting mutant contig: (start2-n1)-----(start2) | (start1)-----(start1+n2)
		res = rbind(c(chr1, start1, start1+n1, 'lower', 'no',2), c(chr2, start2-n2, start2, 'higher', 'no', 1))
	} else if (strand1=="+" & strand2=="-") {
	# lower == + and higher == - ==> reverse compliment of the genomic segment at position < start2 is joined to the genomic segment at position > start1
	# Resulting mutant contig: (start1-n1)-----(start1) | (start2)-----(start2-n2)
		res = rbind(c(chr1, start1-n1, start1, 'lower', 'no',1), c(chr2, start2-n2, start2, 'higher', 'yes', 2))
	} else if (strand1=="+" & strand2=="+") {
	# lower == + and higher == + ==> the genomic segment at position > start2 is joined to the genomic segment at position > start1
	# Resulting mutant contig: (start1-n1)-----(start1) | (start2)-----(start2+n2)
		res = rbind(c(chr1, start1-n1, start1, 'lower', 'no',1), c(chr2, start2, start2+n2, 'higher', 'no', 2))
	}
	colnames(res) = c('chr','start','end','bkpt','use_rev_comp','order')
	res = as.data.frame(res)
	res$start = as.integer(as.character(res$start))
	res$end = as.integer(as.character(res$end))
	refseqs <- getSeq(ref_fasta, GRanges(res$chr, IRanges(start = res$start, end = res$end)))
	res$seq = as.data.frame(refseqs)$x
	res$rev_comp_seq = reverse(as.data.frame(complement(refseqs))$x)
	return(res)
}

create_contig <- function(bkpt) {
	seq1 = bkpt[bkpt$order==1,'seq']
	if(bkpt[bkpt$order==1,'use_rev_comp']=='yes') {
		seq1 = bkpt[bkpt$order==1,'rev_comp_seq']
	}
	seq2 = bkpt[bkpt$order==2,'seq']
	if(bkpt[bkpt$order==2,'use_rev_comp']=='yes') {
		seq2 = bkpt[bkpt$order==2,'rev_comp_seq']
	}
	return(list(
		'mt_contig'=paste(seq1, '|', seq2, sep=" "),
		'primer3_contig'=paste(seq1, seq2, sep=""),
		'bkpt_pos_in_contig'=nchar(seq1)
	))
}

write_primer3_input <- function(p3_input_file, rearr_id, rearr_contig, bkpt_pos_in_contig, distance_from_bkpt) {
	write.table(file=p3_input_file, paste('SEQUENCE_ID=', rearr_id, sep=""), row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, paste('SEQUENCE_TEMPLATE=', rearr_contig, sep=""), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, paste('SEQUENCE_TARGET=', bkpt_pos_in_contig-distance_from_bkpt, ',', 2*distance_from_bkpt, sep=""), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_TASK=pick_sequencing_primers'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_PICK_LEFT_PRIMER=1'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_PICK_INTERNAL_OLIGO=0'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_PICK_RIGHT_PRIMER=1'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_OPT_SIZE=20'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_MIN_SIZE=18'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_MAX_SIZE=23'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_MAX_NS_ACCEPTED=1'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_PRODUCT_SIZE_RANGE=100-300'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('P3_FILE_FLAG=1'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('SEQUENCE_INTERNAL_EXCLUDED_REGION=180,40'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('PRIMER_EXPLAIN_FLAG=1'), append=T, row.names=F, col.names=F, quote=F)
	write.table(file=p3_input_file, c('='), append=T, row.names=F, col.names=F, quote=F)
}
args = commandArgs(TRUE)
rearr_file = args[1]
output_dir = normalizePath(args[2])

rearrangements = read.table(rearr_file, header=F, stringsAsFactors=F, sep="\t")
cmd = paste('grep chr1 ', rearr_file, ' | sed \'s/# //g\' ', sep="")
colnames(rearrangements) = unlist(strsplit(as.character(system(cmd, intern=T)), split="\t"))
ref_fasta <- FaFile(file = REF_FASTA)

### Pick genomic regions present in the mutant contig
all_bkpts = NULL
for(i in 1:nrow(rearrangements)) {
	bkpth = NULL
	bkptl = NULL
	bkpt = pick_genomic_region(
		rearrangements[i,'chr1'], rearrangements[i,'strand1'], rearrangements[i,'start1'],
		rearrangements[i,'chr2'], rearrangements[i,'strand2'], rearrangements[i,'start2']
	)
	bkpt = cbind(rearrangements[i,'id/name'], bkpt)
	colnames(bkpt)[1] = 'id'
	all_bkpts = rbind(all_bkpts, bkpt)
}
### Construct the mutant contig
rearrangements$mt_contig = ""
rearrangements$primer3_contig = ""
rearrangements$bkpt_pos_in_contig = 0
for(id in rearrangements$id) {
	bkpt = all_bkpts[all_bkpts$id==id,]
	contig_info = create_contig(bkpt)
	rearrangements[rearrangements$id==id, 'mt_contig'] = contig_info[['mt_contig']] 
	rearrangements[rearrangements$id==id, 'primer3_contig'] = contig_info[['primer3_contig']]
	rearrangements[rearrangements$id==id, 'bkpt_pos_in_contig'] = contig_info[['bkpt_pos_in_contig']]
}
write.table(rearrangements, file=paste(output_dir, '/mutant_contigs.txt', sep=""), row.names=F, quote=F, sep="\t")
### Design primers
FWD = getwd()
rearrangements$left_primer_seq = ""
rearrangements$right_primer_seq = ""
for(i in 1:nrow(rearrangements)) {
	print(paste("Designing primers for rearrangement with id: ", rearrangements[i,'id/name'], rep=""))
	p3_input_file = paste(output_dir, '/', rearrangements[i,'id/name'], '_primer3_input.txt', sep="")
	p3_outout_file = paste(output_dir, '/', rearrangements[i,'id/name'], '_primer3_output.txt', sep="")
	write_primer3_input(p3_input_file, rearrangements[i,'id/name'], rearrangements[i,'primer3_contig'], rearrangements[rearrangements$id==id, 'bkpt_pos_in_contig'], DISTANCE_FROM_BKPT)
	setwd(PRIMER3_ROOT)
	cmd = paste("./primer3_core -p3_settings_file=../primer3web_v3_0_0_default_settings.txt -output=", p3_outout_file, " ", p3_input_file, sep="")
	system(cmd)
	### Parse primer3 output
	p3_output = read.table(p3_outout_file, sep="=", stringsAsFactors=F)
	lprimer = p3_output[p3_output$V1=='PRIMER_LEFT_0_SEQUENCE','V2']
	rprimer = p3_output[p3_output$V1=='PRIMER_RIGHT_0_SEQUENCE','V2']
	### If one of the primers failed, move the window 50bp in one direction and re-attempt primer creation
	if(length(lprimer)==0 | length(rprimer)==0) {
		print(paste("Failed primer design for rearrangement with id: ", rearrangements[i,'id/name'], rep=""))
		print("Second attempt!")
		bkpt = pick_genomic_region(
			rearrangements[i,'chr1'], rearrangements[i,'strand1'], rearrangements[i,'start1'],
			rearrangements[i,'chr2'], rearrangements[i,'strand2'], rearrangements[i,'start2'], n1=150, n2=250
		)
		refseqs <- getSeq(ref_fasta, GRanges(bkpt$chr, IRanges(start = bkpt$start, end = bkpt$end)))
		bkpt$seq = as.data.frame(refseqs)$x
		bkpt$rev_comp_seq = reverse(as.data.frame(complement(refseqs))$x)
		
		contig_info = create_contig(bkpt)
		write_primer3_input(p3_input_file, rearrangements[i,'id/name'], contig_info[['primer3_contig']], contig_info[['bkpt_pos_in_contig']], DISTANCE_FROM_BKPT)
		cmd = paste("./primer3_core -p3_settings_file=../primer3web_v3_0_0_default_settings.txt -output=", p3_outout_file, " ", p3_input_file, sep="")
		system(cmd)
		### Parse second primer3 output
		p3_output = read.table(p3_outout_file, sep="=", stringsAsFactors=F)
		lprimer = p3_output[p3_output$V1=='PRIMER_LEFT_0_SEQUENCE','V2']
		rprimer = p3_output[p3_output$V1=='PRIMER_RIGHT_0_SEQUENCE','V2']
		### If one of the primars failed again, move the window 50bp in the opposite direction and re-attempt primer creation
		if(length(lprimer)==0 | length(rprimer)==0) {
			print(paste("Failed primer design for rearrangement with id: ", rearrangements[i,'id/name'], rep=""))
			print("Third attempt!")
			bkpt = pick_genomic_region(
				rearrangements[i,'chr1'], rearrangements[i,'strand1'], rearrangements[i,'start1'],
				rearrangements[i,'chr2'], rearrangements[i,'strand2'], rearrangements[i,'start2'], n1=250, n2=150
			)
			contig_info = create_contig(bkpt)
			write_primer3_input(p3_input_file, rearrangements[i,'id/name'], contig_info[['primer3_contig']], contig_info[['bkpt_pos_in_contig']], DISTANCE_FROM_BKPT)
			cmd = paste("./primer3_core -p3_settings_file=../primer3web_v3_0_0_default_settings.txt -output=", p3_outout_file, " ", p3_input_file, sep="")
			system(cmd)
			p3_output = read.table(p3_outout_file, sep="=", stringsAsFactors=F)
			lprimer = p3_output[p3_output$V1=='PRIMER_LEFT_0_SEQUENCE','V2']
			rprimer = p3_output[p3_output$V1=='PRIMER_RIGHT_0_SEQUENCE','V2']
		}
	} 
	cmd = paste("./primer3_core -format_output -p3_settings_file=../primer3web_v3_0_0_default_settings.txt -output=", p3_outout_file, " ", p3_input_file, sep="")
	system(cmd)
	if(length(lprimer)==0) {lprimer='failed'}
	if(length(rprimer)==0) {rprimer='failed'}
	rearrangements[i,'left_primer_seq'] = lprimer
	rearrangements[i,'right_primer_seq'] = rprimer
	setwd(FWD)
}
write.table(rearrangements, file=paste(output_dir, '/mutant_contigs.txt', sep=""), row.names=F, quote=F, sep="\t")
