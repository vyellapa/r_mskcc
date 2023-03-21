#!/usr/bin/env Rscript

#
# gundemg@mskcc.org
# Brass for related samples from the same patient
#

#####################################
source('/ifs/work/leukgen/home/gg10/soft/plot_rearrangements_copynumber_and_battenberg.R')
#library(quantsmooth)
#library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19) # human genome
library(VariantAnnotation) # 

chrs = read.table('/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta.fai')
cytoband = read.table('/ifs/work/leukgen/home/gg10/ref_data/ucsc_goldenpath_hg19_cytoband.txt')
UTILS_DIR = '/ifs/work/leukgen/home/gg10/soft/integrated.circos/utils/'
source(paste(UTILS_DIR, 'RCircos.Scatter.Plot.color.R', sep=''))
source(paste(UTILS_DIR, 'my.RCircos.Tile.Plot.R', sep=''))
source(paste(UTILS_DIR, 'calcIntermutDist.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Get.Plot.Data.nosort.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Chromosome.Ideogram.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Track.Outline.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Line.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Link.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Track.Positions.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Gene.Connector.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Heatmap.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Validate.Genomic.Data.my.R', sep=''))
source(paste(UTILS_DIR, 'toPyr.R', sep=''))
source(paste(UTILS_DIR, 'processSubs.R', sep=''))
source(paste(UTILS_DIR, 'RCircos.Gene.Name.Plot.my.R', sep=''))
source(paste(UTILS_DIR, 'merge.with.order.R', sep=''))
source(paste(UTILS_DIR, 'read.ascat.R', sep=''))
source(paste(UTILS_DIR, 'read.brassI.R', sep=''))
source(paste(UTILS_DIR, 'read.brassII.R', sep=''))
source(paste(UTILS_DIR, 'plotDP.R', sep=''))
source(paste(UTILS_DIR, 'generateHist.R', sep=''))
source(paste(UTILS_DIR, 'getMutTables.R', sep=''))
source(paste(UTILS_DIR, 'getPindelVcf.R', sep=''))

##################################### functions
plot_BAF_BBG <- function(baf_data, bb_segments, chr, start=NULL, end=NULL, is_bbg=TRUE) {
	xlim = c(1,chrs[chrs$V1==chr,'V2'])
	if(!is.null(start) & !is.null(end)) {
		xlim = c(start, end)
	}
	plot(c(), ylim=c(0,1), xlim=xlim, bty="n", yaxt="n", xaxt="n", xlab="", ylab="", yaxs="i", xaxs="i")
	axis(2, at=axisTicks(usr=c(0,1), log=F), las=2)
	title(ylab="BAF")
	points(x=baf_data[baf_data$Chromosome==chr,2], y=baf_data[baf_data$Chromosome==chr,3], col='darkgrey', pch=16, cex=0.3)
	sel = bb_segments$chr==chr
	
	if(!is.null(start) & !is.null(end)) {
		sel = (
		(sel & bb_segments$startpos>start & bb_segments$startpos<end) | (sel & bb_segments$endpos>start & bb_segments$endpos<end) |
		(sel & bb_segments$startpos>end & bb_segments$endpos<end) | (sel & bb_segments$startpos<start & bb_segments$endpos>end)
		) & (bb_segments$endpos - bb_segments$startpos>20000)
	}
	segments(
		x0 = bb_segments[sel, 'startpos'] + 1, x1 = bb_segments[sel, 'endpos'],
		y0 = bb_segments[sel, 'BAF'], lwd = 2, col = bb_segments[sel, 'colour']
	)
	if(is_bbg) {
	bb_sub = bb_segments[sel,]
	for (i in 1:nrow(bb_sub)) {
		if(is_bbg) {
			text(
				(as.numeric(bb_sub[i,"startpos"])+as.numeric(bb_sub[i,"endpos"]))/2,
				as.numeric(bb_sub[i,"BAF"])-0.04,
				paste(bb_sub[i,"nMaj1_A"],"+",bb_sub[i,"nMin1_A"],": ",100*round(as.numeric(bb_sub[i,"frac1_A"]),3),"%",sep=""),
				cex = 0.8
			)
			if(!is.na(bb_sub[i,"nMaj2_A"])) {
				text((as.numeric(bb_sub[i,"startpos"])+as.numeric(bb_sub[i,"endpos"]))/2, as.numeric(bb_sub[i,"BAF"])-0.08,
				paste(bb_sub[i,"nMaj2_A"],"+",bb_sub[i,"nMin2_A"],": ",100*round(as.numeric(bb_sub[i,"frac2_A"]),3),"%",sep=""), cex = 0.8)
			}
		} else {
			text((as.numeric(bb_sub[i,"startpos"])+as.numeric(bb_sub[i,"endpos"]))/2, as.numeric(bb_sub[i,"BAF"])-0.04, paste('TCN=',bb_sub[i,"cn"],sep=""), cex=0.8)

		}
	}
	}
}

###################################################
# MAIN

# Parse arguments and complain
args=commandArgs(TRUE)
bedpe_file = args[1]
bg_file = args[2]
bb_file = args[3]
baf_file = args[4]
output_file_stub = args[5]
flip_higher = args[6]
legend_pos = args[7]
rearr_col = args[8]

coord = "none"
plot_chrs = "none"
cn_win_size = 1e5
if(length(args)>=9) {
	if(grepl('chrs',args[9])) {
		plot_chrs = sub('^chrs_','',args[9])
	} else {
		coord = args[9]
		cn_win_size = 5000
	}
}
subs_file = "none"
if(length(args)>=10) {
	subs_file = args[10]
}
yrange = NULL
if(length(args)>=11) {
	yrange = args[11]
	yrange = unlist(strsplit(yrange, split="_"))
	yrange = c(as.double(yrange[1]), as.double(yrange[2]))

}

image_type = "pdf"
cn_log_scale = FALSE
facet = FALSE
image_width = 9

# Read the breakpoints into a bedpe format data frame
cat(paste('Reading rearrangements...', '\n',sep=''))
count = as.integer(unlist(strsplit(system(paste('wc -l ', bedpe_file, sep=''), intern=TRUE), split=" "))[1])
headered = as.integer(unlist(strsplit(system(paste('grep chr1 ', bedpe_file, ' | wc -l', sep=''), intern=TRUE), split=" "))[1])
if(count>0) {
	if(headered==0) {
		brass = read.table(file=bedpe_file, sep='\t', stringsAsFactors=F, header=F, comment.char="#")
		brass$id = rownames(brass)
		bedpe = cbind(brass[,c(1:7)],'black',brass[,9:10])
		bedpe = cbind(bedpe,'a')
	}
	if(headered==1) {
		brass = read.csv(file=bedpe_file, sep='\t', stringsAsFactors=F, header=T, comment.char="")
		if('presence.summary' %in% colnames(brass)) {
			bedpe = brass[,unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id color strand1 strand2', split=' '))]
			bedpe = cbind(bedpe,brass[,'presence.summary'])
		} else {
			brass = cbind("","black",brass[,unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id strand1 strand2', split=' '))])
			colnames(brass)[1:2] = c('presence.summary','color')
			bedpe = brass[,unlist(strsplit('chr1 start1 end1 chr2 start2 end2 id color strand1 strand2 presence.summary', split=' '))]
		}
	}
	colnames(bedpe) = c('m_chrl','m_l5','m_l3','m_chrh','m_h5','m_h3','id','color','strl','strh','presence.summary')
	bedpe$color = as.character(bedpe$color)
	bedpe$m_l5 = as.integer(bedpe$m_l5)
	bedpe$m_l3 = as.integer(bedpe$m_l3)
	bedpe$m_h5 = as.integer(bedpe$m_h5)
	bedpe$m_h3 = as.integer(bedpe$m_h3)
} else {
	bedpe = NULL
}
cat(paste('Read rearrangements...', '\n',sep=''))
if(flip_higher) {
	strh = bedpe$strh
	strh[bedpe$strh == '-'] = '+'
	strh[bedpe$strh == '+'] = '-'
	bedpe$strh = strh
}
# Read the raw bin counts data, Battenberg segments and raw BAF data
cat(paste('Reading NGS copy number data...', '\n',sep=''))
chr = NA
start = NA
end = NA
xlim = NA
if(coord!='none') {
	coordinates = unlist(strsplit(coord,split="_"))
	chr = coordinates[1]
	start = as.integer(coordinates[2])
	end = as.integer(coordinates[3])
	xlim = c(start, end)
}
if(coord!="none") {
	tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
	cmd = paste("less ", bg_file, " | awk '$1==\"", chr, "\" && $2>", start, " && $3<", end, "' > ", tmp.file.name, sep="")
	system(cmd)
	bg = read.table(file=tmp.file.name, header=F)
	system(paste('rm ', tmp.file.name, sep=""))
} else {
	bg = read.table(file=bg_file, header=F)
}
cat(paste('Read NGS copy number data...', '\n',sep=''))
baf_data = NULL
cat(paste('Reading BAF data...', '\n',sep=''))
if(!is.null(baf_file)) {
	if(grepl('mutantBAF',baf_file)) {
		cat(paste('Reading raw BAF data...', '\n',sep=''))
		if(coord!='none') {
			tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
			cmd = paste("zless ", baf_file, " | awk '$2==", chr, " && $3>=", start, " && $3<=", end, "' > ", tmp.file.name, sep="")
			system(cmd)
			baf_data = read.table(tmp.file.name, header=F, stringsAsFactors=F, row.names=1)
			system(paste('rm ', tmp.file.name, sep=""))
		} else {
			baf_data = read.table(baf_file, header=T, stringsAsFactors=F, row.names=1)
		}
		colnames(baf_data) = c('Chromosome','Position','BAF')
	}
	if(grepl('logR_Baf_segmented.vcf.gz', baf_file) | grepl('.BAFsegmented.txt.gz', baf_file)) {
		cat(paste('Reading BAF data from Battenberg...', '\n',sep=''))
		if(grepl('logR_Baf_segmented.vcf.gz', baf_file)) {
			tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
			cmd = paste("zless ", baf_file, " | grep -v \"#\" | awk 'OFS=\"\t\"{split($8,a,\";\"); split(a[2],b,\"=\"); print $1,$2,b[2]}' > ", tmp.file.name, sep="")
			if(coord!="none") {
				cmd = paste(
					"zless ", baf_file, " | grep -v \"#\" | awk 'OFS=\"\t\"{if($1==", chr, " && $2>", start, " && $2<", end,
					") {split($8,a,\";\"); split(a[2],b,\"=\"); print $1,$2,b[2]}}' > ", tmp.file.name, sep=""
				)
			}
			system(cmd)
			baf_data = read.table(tmp.file.name, header=F)
			system(paste('rm ', tmp.file.name, sep=""))
			colnames(baf_data) = c('Chromosome','Position','BAF')
		} else {
			if(coord!="none") {
				tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
				cmd = paste(
					"zless ", baf_file, " | grep -v \"#\" | awk 'OFS=\"\t\"{if($2==", chr, " && $3>=", start, " && $3<=", end,
					") {split($8,a,\";\"); split(a[1],b,\"=\"); print $1,$2,b[2]}}' > ", tmp.file.name, sep=""
				)
				system(cmd)
				baf_data = read.table(tmp.file.name, header=F, stringsAsFactors=F)
				colnames(baf_data) = c('Chromosome','Position','BAF','BAFphased','BAFseg')
				system(paste('rm ', tmp.file.name, sep=""))
			}
			baf_data = read.table(baf_file, header=T, stringsAsFactors=F)
		}
	}
}
cat(paste('Read BAF data...', '\n',sep=''))
bb_segments = NULL
is_bbg = TRUE
cat(paste('Reading BBG segements...', '\n',sep=''))
if(!is.null(bb_file)) {
	bb_segments = read.delim(file=bb_file, header=T)
	if(grepl('_subclones.txt$', bb_file)) {
		bb_segments$colour = rep('darkgreen', nrow(bb_segments))
		bb_segments[bb_segments$frac1_A<1,'colour'] = 'red'
	} else if(grepl('jabba.seg', bb_file)) {
		is_bbg = FALSE
		cat(paste('Calculating BAF for non-BBG segments...', '\n',sep=''))
		bb_segments$BAF = rep(NA,nrow(bb_segments))
		for(i in 1:nrow(bb_segments)) {
		bb_segments[i,'BAF'] = mean(baf_data[baf_data$Chromosome==bb_segments[i,'chrom'] & baf_data$Position>bb_segments[i,'start'] & baf_data$Position<bb_segments[i,'end'],5])
		}
		bb_segments = bb_segments[,c('chrom','start','end','cn','BAF')]
		bb_segments$colour = rep('darkgreen',nrow(bb_segments))
		colnames(bb_segments) = c('chr','startpos','endpos','cn','BAF','colour')
		
	}
}
cat(paste('Read BBG segements...', '\n',sep=''))
# If no coordinate is specified, plot rearrangements, BAF and BBG segments for all chromosomes seperately
cat(paste('Plotting...', '\n',sep=''))
if(coord == 'none') {
	library(quantsmooth)
	chrs_to_plot = c()
	if(plot_chrs!='none') {
		chrs_to_plot = unlist(strsplit(plot_chrs, split="_"))
	} else if(bg_file!="none") {
		chrs_to_plot = as.character(unique(bg$V1))
	} else {
		chrs_to_plot = union(unique(bedpe[,1]), unique(bedpe[,4]))
	}
	if(plot_chrs=='none') {
		for(chr in chrs_to_plot) {
			print(paste('Plotting chr ', chr, sep=''))
			height = 4
			if(!is.null(baf_data)) {height=12}
			if(image_type == 'pdf') {
				pdf(paste(output_file_stub, '_chr_', chr, '.pdf', sep=""), width=image_width, height=height)
			}
			if(image_type == 'svg') {
				svg(paste(output_file_stub, '_chr_', chr, '.svg', sep=""), width=image_width, height=height)
			}
			if(!is.null(baf_data) & chr!='Y') {
				par(mfrow=c(2,1))
				plot_BAF_BBG(baf_data=baf_data, bb_segments=bb_segments, chr=chr, is_bbg=is_bbg)
			}
			plot_rearrangements(
				bedpe=bedpe[bedpe$m_chrl==chr | bedpe$m_chrh==chr,],
				chrs=c(chr), cn_bedgraph=bg, segments=NULL, lwd=2, yrange=yrange,
				rearr_col=rearr_col, legend_pos=legend_pos, cn_win_size=cn_win_size
			)
			dev.off()
		}
	} else {
		library(quantsmooth)
		print(paste('Plotting multiple chromosomes: ', plot_chrs, sep=''))
		height = 4
		if(image_type == 'pdf') {
			pdf(paste(output_file_stub, '_chrs_', gsub(',','_',plot_chrs), '.pdf', sep=""), width=image_width, height=height)
		}
		if(image_type == 'svg') {
			svg(paste(output_file_stub, '_chrs_', gsub(',','_',plot_chrs), '.svg', sep=""), width=image_width, height=height)
		}
		plot_rearrangements(
			bedpe=bedpe[bedpe$m_chrl %in% chrs_to_plot | bedpe$m_chrh %in% chrs_to_plot,],
			chrs=chrs_to_plot, cn_bedgraph=bg, segments=NULL, lwd=2, yrange=yrange,
			rearr_col=rearr_col, legend_pos=legend_pos, cn_win_size=cn_win_size
		)
		dev.off()
	}
}
# If substitutions are available, plot both inter-mutation distance and rearrangements for the specific genomic region
if(subs_file!="none") {
	print('Will plot SUBS, REARRS and CN data.')
	substitutions = NULL
	if(grepl('vcf',subs_file)) {
		subs.data <- getMutTables(subs_file, onlyPASSED=TRUE, onlyPassedAsrd=FALSE)
		substitutions = subs.data$muts
	} else {
		subs.data = read.table(subs_file, header=T, sep="\t", stringsAsFactors=F)
		subs.data$context = paste(subs.data$REF, subs.data$ALT, sep=">")
		subs.data[subs.data$context=='A>C','context'] = 'T>G'
		subs.data[subs.data$context=='A>G','context'] = 'T>C'
		subs.data[subs.data$context=='A>T','context'] = 'T>A'
		subs.data[subs.data$context=='G>A','context'] = 'C>T'
		subs.data[subs.data$context=='G>C','context'] = 'C>G'
		subs.data[subs.data$context=='G>T','context'] = 'C>A'
		subs.data = subs.data[,c('CHR','START','context')]
		colnames(subs.data) = c('chroms','starts','context')
		for(chrom in unique(subs.data$chroms)) {
			ssubs.data = subs.data[subs.data$chroms==chrom,]
			substitutions = rbind(substitutions, ssubs.data[order(ssubs.data$starts),])
		}
		substitutions$chroms = paste0('chr',substitutions$chroms)
	}
	substitutions = droplevels(substitutions[substitutions$chroms!="chrY",])
	substitutions$intermut_dist = -1
	for(chrom in as.character(unique(substitutions$chroms))) {
		substitutions_sub = substitutions[substitutions$chroms==chrom,]
		if(nrow(substitutions_sub)>1) {
			substitutions_sub$intermut_dist = c(substitutions_sub[2:nrow(substitutions_sub),'starts'], 0) - substitutions_sub$starts
		}
		substitutions_sub[substitutions_sub$intermut_dist<0,'intermut_dist'] = mean(substitutions_sub[substitutions_sub$intermut_dist>0,'intermut_dist'])
		substitutions[substitutions$chroms==chrom,'intermut_dist'] = substitutions_sub$intermut_dist
	}
	substitutions$color = "royalblue"					# C>A
	substitutions[grep('C>G',substitutions$context),'color'] = "black"	# C>G
	substitutions[grep('C>T',substitutions$context),'color'] = "red"	# C>T
	substitutions[grep('T>A',substitutions$context),'color'] = "grey"	# T>A
	substitutions[grep('T>C',substitutions$context),'color'] = "green2"	# T>C
	substitutions[grep('T>G',substitutions$context),'color'] = "hotpink"	# T>G
	pdf(paste(output_file_stub, '__zoom_', coord, '.pdf', sep=""), width=9, height=9)
	# Read the BG file for the selected chromosome only
	tmp.file.name = tempfile(pattern = "plot_rearr_", tmpdir = '.')
	cmd = paste("less ", bg_file, " | awk '$1==\"", chr, "\"' > ", tmp.file.name, sep="")
	system(cmd)
	bg_chr = read.table(file=tmp.file.name, header=F)
	system(paste('rm ', tmp.file.name, sep=""))
	# Plot all the mutations on the whole chromosome
	selected_subs2 = substitutions$chroms==paste0('chr',chr)
	par(mfrow=c(2,2))
	xlim1 = c(1,as.integer(as.character(chrs[chrs$V1==chr,'V2'])))
	plot(c(), ylim=c(0,30), xlim=xlim1, bty="n", yaxt="n", xaxt="n", xlab="", ylab="", yaxs="i", xaxs="i")
	axis(2, at=axisTicks(usr=c(0,30), log=F), las=2)
	title(ylab="log(inter-mutation distance)")
	points(x=substitutions[selected_subs2,'starts'], y=log2(substitutions[selected_subs2,'intermut_dist']), pch=16, col=substitutions[selected_subs2,'color'])
	axis(1, at=axisTicks(usr=xlim1, log=F), las=2)
	print('Done with whole chr.')
	# plot only those mutations that fall within the selected genomic region
	selected_subs1 = substitutions$chroms==paste0('chr',chr) & substitutions$starts>start & substitutions$starts<end
	plot(c(), ylim=c(0,30), xlim=xlim, bty="n", yaxt="n", xaxt="n", xlab="", ylab="", yaxs="i", xaxs="i")
	axis(2, at=axisTicks(usr=c(0,30), log=F), las=2)
	title(ylab="log(inter-mutation distance)")
	points(x=substitutions[selected_subs1,'starts'], y=log2(substitutions[selected_subs1,'intermut_dist']), pch=16, col=substitutions[selected_subs1,'color'])
	axis(1, labels=rep("",length(axisTicks(usr=xlim, log=F))), at=axisTicks(usr=xlim, log=F), las=2)
	legend(x="topright", legend=c('C>A','C>G','C>T','T>A','T>C','T>G'), col=c("royalblue","black","red","grey","green2","hotpink"), pch=19, pt.cex=0.6, horiz=F)
	print('Done with genomic region.')
	# Plot rearrangements
	plot_rearrangements(
		bedpe=bedpe[bedpe$m_chrl==chr | bedpe$m_chrh==chr,], chrs=c(chr), cn_bedgraph=bg_chr, segments=NULL, lwd=2,
		yrange=c(0,5), rearr_col=rearr_col, legend_pos=legend_pos, cn_win_size=1e5
	)
	plot_rearrangements(
		bedpe=bedpe, chrs=c(chr), cn_bedgraph=bg, segments=NULL, cn_win_size=cn_win_size, yrange=yrange,
		xlim=xlim, lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col, cn_log_scale=cn_log_scale
	)
	dev.off()
}
# If substitutions are not available, then plot BAF data, BBG and rearrangments for the specified genomic region only
if(coord != "none" & subs_file=="none" & !facet){
	library(quantsmooth)
	height = 7
	if(!is.null(baf_data)) {height=12}
	if(!grepl('chrs:',coord)) {
		if(image_type == 'pdf') {
			pdf(paste(output_file_stub, '_chr_', chr, '__zoom_', coord, '.pdf', sep=""), width=image_width, height=height)
		}
		if(image_type == 'svg') {
			svg(paste(output_file_stub, '_chr_', chr, '__zoom_', coord, '.svg', sep=""), width=image_width, height=height)
		}
		if(!is.null(baf_data)) {
			par(mfrow=c(2,1))
			plot_BAF_BBG(baf_data=baf_data, bb_segments=bb_segments, chr=chr, start=start, end=end, is_bbg=is_bbg)
		}
		chrs_to_plot = c(chr)
		bedpe=bedpe[((bedpe$m_chrl==chr & bedpe$m_l5>start & bedpe$m_l3<end) | (bedpe$m_chrh==chr & bedpe$m_h5>start & bedpe$m_h3<end)),]
	} else if (grepl('chrs:',coord)) {
		chrs_to_plot = unlist(strsplit(unlist(strsplit(coord,split=":"))[2], split=','))
		xlim = NULL
		if(image_type == 'pdf') {
			coord = gsub(',','_',gsub(':','',coord))
			pdf(paste(output_file_stub, '__zoom_', coord, '.pdf', sep=""), width=image_width, height=height)
		} else if (image_type == 'svg') {
			svg(paste(output_file_stub, '__zoom_', coord, '.pdf', sep=""), width=image_width, height=height)
		}
	}
	plot_rearrangements(
		bedpe=bedpe, chrs=chrs_to_plot, cn_bedgraph=bg, segments=NULL, cn_win_size=cn_win_size, yrange=yrange,
		xlim=xlim, lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col, cn_log_scale=cn_log_scale
	)

	dev.off()
} else if (coord !="none" & subs_file=="none" & facet) {
	library(quantsmooth)
	chrs_to_plot = c(chr)
	coord = gsub(',','_',gsub(':','',coord))
	bedpe=bedpe[((bedpe$m_chrl==chr & bedpe$m_l5>start & bedpe$m_l3<end) | (bedpe$m_chrh==chr & bedpe$m_h5>start & bedpe$m_h3<end)),]
	categories = as.character(unique(bedpe$presence.summary))
	print(paste0('Will facet for ', length(categories), ' categories'))
	image_width = 4 * length(categories)
	height = 4
	if(length(categories)>3) {
		image_width = 4 * 2
		height = ceiling(length(categories)/2)
	}
	if(image_type == 'pdf') {
		pdf(paste(output_file_stub, '__zoom_', coord, '.pdf', sep=""), width=image_width, height=height)
	} else if (image_type == 'svg') {
		svg(paste(output_file_stub, '__zoom_', coord, '.pdf', sep=""), width=image_width, height=height)
	}
	if(length(categories)>3) {
		par(mfrow=c(ceiling(length(categories)/2), 2))
	} else {
		par(mfrow=c(1, length(categories)))
	}
	for(category in categories) {
		plot_rearrangements(
			bedpe=bedpe[bedpe$presence.summary==category,], chrs=chrs_to_plot, cn_bedgraph=bg, segments=NULL, cn_win_size=cn_win_size, yrange=yrange,
			xlim=xlim, lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col, cn_log_scale=cn_log_scale
		)
	}
	dev.off()
}

cat(paste('DONE.', '\n',sep=''))
