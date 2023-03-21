# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/run_plot_rearrangements_copynumber_and_battenberg_bedpe2.R merged_brass.txt Cerebellum-167276T/against_blood_normal/bga/intermediates/W0006077F.ngscn.abs_cn.bg.gz none C none topright pdf TRUE 0,10 10000 custom split_plot

source('/ifs/work/leukgen/home/gg10/soft/plot_rearrangements_copynumber_and_battenberg.R')
library(quantsmooth)
library(GenomicRanges)

colour.palette = c(colors()[555], colors()[81], colors()[28], colors()[95], colors()[143], colors()[11], colors()[624], colors()[619], colors()[114], colors()[101], colors()[490])
check_presence <- function(data, n_cutoff) {
	samples = sapply(strsplit(grep('*\\.count', colnames(data), value=T), split='.', fixed=T),"[[",1)
	m = mat.or.vec(nrow(data), length(samples))
	m.str = mat.or.vec(nrow(data), length(samples))
	for(i in 1:length(samples)) {
		sample = samples[i]
		scol = paste(sample,'.count', sep='')
		m[as.numeric(data[,scol]) >= n_cutoff,i] = 1
		m.str[as.numeric(data[,scol]) >= n_cutoff,i] = sample
	}
	m.str[which(m.str=="0")] = "-"
	colnames(m) = paste('presence.', samples, sep="")
	m = as.data.frame(m)
	m$presence.summary = do.call(paste, c(as.data.frame(m.str), sep = "/"))
	return(m)
}

args=(commandArgs(trailingOnly=TRUE))
bedpe_file = args[1]
bg_file = args[2]
bb_file = args[3]
sname = args[4]
coord = args[5]
legend_pos = args[6]
image_type = args[7]
colour_by_presence = args[8]
n_cutoff = 3
yrange = NULL
if(length(args)>=9) {
	yrange = as.numeric(unlist(strsplit(args[9], split=',')))
}
cn_win_size = 100000
if(length(args)>=10) {
	cn_win_size = as.numeric(args[10])
}
rearr_col = 'type'
if(length(args)>=11) {
	rearr_col = 'custom'
}
split_plot = FALSE
if(length(args)>=12) {
	split_plot = TRUE
}
# Read the breakpoints into a bedpe format data frame
count = as.integer(unlist(strsplit(system(paste('wc -l ', bedpe_file, sep=''), intern=TRUE), split=" "))[1])
if(count>0) {
	if(colour_by_presence=="T" || colour_by_presence=='TRUE') {
		brass = read.delim(file=bedpe_file, sep='\t', stringsAsFactors=F, header=T)
		brass.sub = brass[which(brass[,paste(sname,'.count',sep='')]>=n_cutoff),]
		brass.sub = cbind(brass.sub, check_presence(brass.sub, n_cutoff=n_cutoff))
		brass.sub$id = rownames(brass.sub)
		presence.categories = table(brass.sub$presence.summary)
		presence.categories = as.data.frame(presence.categories[order(presence.categories, decreasing=T)])
		presence.categories$category = rownames(presence.categories)
		colnames(presence.categories)[1] = 'count'
		rownames(presence.categories) = NULL
		if(nrow(presence.categories)>length(colour.palette)) {
			presence.categories = cbind(presence.categories, c(colour.palette, rainbow(nrow(presence.categories)-length(colour.palette))))
		} else {
			presence.categories = cbind(presence.categories, colour.palette[1:nrow(presence.categories)])
		}
		colnames(presence.categories)[3] = 'color'
		brass.sub = merge(x=brass.sub, y=presence.categories[,c('category','color')], by.x='presence.summary', by.y='category')
		brass.sub = brass.sub[brass.sub[,paste('presence.',sname,sep='')]>0,]
		brass.sub = brass.sub[,c('chr1','start1','end1','chr2','start2','end2','id','color','str1','str2','presence.summary')]
		bedpe = brass.sub
	} else {
		brass = read.delim(file=bedpe_file, sep='\t', stringsAsFactors=F, header=T, comment.char="#")
		brass$id = rownames(brass)
		bedpe = cbind(brass[,c(1:6,16)],'black',brass[,9:10])
		bedpe[brass[,10]=='-',10] = "+"
		bedpe[brass[,10]=='+',10] = "-"
		bedpe[,8] = as.character(bedpe[,8])
		bedpe[brass$V14=="_",8] = 'red'
		bedpe = cbind(bedpe,'a')
	}
	colnames(bedpe) = c('m_chrl','m_l5','m_l3','m_chrh','m_h5','m_h3','id','color','strl','strh','presence.summary')
} else {
	bedpe = NULL
}
# Read the copy number data
bg = NULL
if(bg_file!='none') {
	bg = read.table(file=bg_file)
}
if(bb_file != "none") {
	# Read the segmentation file from Battenberg algorithm
	bb_segments = read.table(file=bb_file, header=T)
	bb_segments$colour = rep('black', nrow(bb_segments))
}
if(bb_file == 'none') {
	bb_segments=NULL
}
if(coord == 'none') {
	chrs_to_plot = c()
	if(bg_file!="none") {
		chrs_to_plot = as.character(unique(bg$V1))
	} else {
		chrs_to_plot = union(unique(bedpe[,1]), unique(bedpe[,4]))
	}
	for(chr in chrs_to_plot) {
		if(image_type == 'pdf') {
			pdf(paste('chr_', chr, '.pdf', sep=""), width=24, height=15)
		}
		if(image_type == 'svg') {
			svg(paste('chr_', chr, '.svg', sep=""), width=24, height=15)
		}
		bedpe.sub = bedpe[bedpe$m_chrl==chr | bedpe$m_chrh==chr,]
		par(mfrow=c(2,1))
		if(nrow(bedpe.sub)==0) {legend_pos='none'}
		plot_rearrangements(bedpe=bedpe.sub, chrs=c(chr), cn_bedgraph=bg, segments=bb_segments, lwd=2, legend_pos=legend_pos, rearr_col=rearr_col)
		plot_rearrangements(bedpe=bedpe.sub, chrs=c(chr), cn_bedgraph=bg, segments=bb_segments, lwd=2, legend_pos=legend_pos, rearr_col=rearr_col)
		dev.off()
	}
}
if(coord != "none"){
	coordinates = unlist(strsplit(coord,split="_"))
	if(image_type == 'pdf') {
		pdf(paste(sname,'_chr_', coordinates[1], '__zoom_', coord, '.pdf', sep=""), width=12)
	}
	if(image_type == 'svg') {
		svg(paste(sname,'_chr_', coordinates[1], '__zoom_', coord, '.svg', sep=""), width=12)
	}
	bedpe.sub = bedpe[(
		(bedpe$m_chrl==coordinates[1] & bedpe$m_l5 > as.integer(coordinates[2]) & bedpe$m_l3 < as.integer(coordinates[3])) |
		(bedpe$m_chrh==coordinates[1] & bedpe$m_h5 > as.integer(coordinates[2]) & bedpe$m_h3 < as.integer(coordinates[3]))
	),]
	if(split_plot) {
		t = table(bedpe.sub$presence.summary)[order(table(bedpe.sub$presence.summary), decreasing=T)]
		for(cat in names(t[t>=5])) {
			plot_rearrangements(
				bedpe=bedpe.sub[bedpe.sub$presence.summary==cat,], chrs=c(coordinates[1]), cn_bedgraph=bg,
				segments=bb_segments, yrange=yrange, cn_win_size=cn_win_size,
				xlim=c(as.numeric(coordinates[2]), as.numeric(coordinates[3])),
				lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col
			)
		}
		plot_rearrangements(
			bedpe=bedpe.sub[bedpe.sub$presence.summary %in% names(t[t<5]),], chrs=c(coordinates[1]), cn_bedgraph=bg,
			segments=bb_segments, yrange=yrange, cn_win_size=cn_win_size,
			xlim=c(as.numeric(coordinates[2]), as.numeric(coordinates[3])),
			lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col
		)
	} else {
		plot_rearrangements(
			bedpe=bedpe.sub, chrs=c(coordinates[1]), cn_bedgraph=bg,
			segments=bb_segments, yrange=yrange, cn_win_size=cn_win_size,
			xlim=c(as.numeric(coordinates[2]), as.numeric(coordinates[3])),
			lwd=2, cn_cex=0.4, legend_pos=legend_pos, rearr_col=rearr_col
		)
	}
	dev.off()
}
