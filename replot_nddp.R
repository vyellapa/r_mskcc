#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-d", "--datapath", type="character", nargs=1, help="Absolute path to the directory containing all DP input files.")
parser$add_argument("-x", "--xsample", type="character", nargs=1, help="Name of the sample to be plotted on x-axis.")
parser$add_argument("-y", "--ysample", type="character", nargs=1, help="Name of the sample to be plotted on y-axis.")
parser$add_argument("-xl", "--xlabel", type="character", nargs=1, help="Label to show on x-axis. Will default to sample name if not specified.")
parser$add_argument("-yl", "--ylabel", type="character", nargs=1, help="Label to show on y-axis. Will default to sample name if not specified.")
parser$add_argument("-pa", "--passignments", type="character", nargs=1, help="Cluster assignment file from patient-level DP.")
parser$add_argument("-pn", "--pnodes", type="character", nargs=1, help="Node positions file from patient-level DP.")
parser$add_argument("-sn", "--snodes", default=NULL, type="character", nargs=1, help="Node positions file for sample-level DP.")
parser$add_argument("-m", "--max", default=2, type="double", nargs=1, help="Max value to plot on x- and y-axes.")
parser$add_argument("-o", "--output", default=2, type="character", nargs=1, help="Name of the output PNG file.")

args <- parser$parse_args()

# Read cluster info and cluster assignment info
print(paste('Patient-level clusters =', args$pnodes))
print(paste('Patient-level assingment file =', args$passignments))
print(paste('Sample on x-axis =', args$xsample))
print(paste('Sample on y-axis =', args$ysample))
print(paste('Sample-level clusters =', args$snodes))

pnodepos = read.table(args$pnodes, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
pnodepos = pnodepos[order(pnodepos$no.muts.in.cluster, decreasing=T),]
passignment = read.table(args$passignments, header=T, sep="\t")
passignment$mutkey = paste(passignment$chr, passignment$end)
passignment$col = 'black'
vaf_cutoff = 0.05
for(cluster in unique(pnodepos$cluster.no)) {passignment$col[passignment$cluster==cluster] = pnodepos[pnodepos$cluster.no==cluster,'colour']}

# Read subclonal fraction data for the two samples
sfrac_x_file = paste(args$datapath, "/", args$xsample, "_allDirichletProcessInfo.txt", sep="")
sfrac_y_file = paste(args$datapath, "/", args$ysample, "_allDirichletProcessInfo.txt", sep="")
sfrac_x = read.table(sfrac_x_file, header=T, sep="\t", stringsAsFactors=F)
sfrac_x$vaf = sfrac_x$mut.count / (sfrac_x$mut.count + sfrac_x$WT.count)
sfrac_x$mutkey = paste(sfrac_x$chr, sfrac_x$end)
sfrac_y = read.table(sfrac_y_file, header=T, sep="\t", stringsAsFactors=F)
sfrac_y$vaf = sfrac_y$mut.count / (sfrac_y$mut.count + sfrac_y$WT.count)
sfrac_y$mutkey = paste(sfrac_y$chr, sfrac_y$end)
sfrac_x = merge(x=sfrac_x, y=passignment[,c('mutkey','col')], by.x='mutkey', by.y='mutkey')
sfrac_y = merge(x=sfrac_y, y=passignment[,c('mutkey','col')], by.x='mutkey', by.y='mutkey')
sel = !is.na(sfrac_x$subclonal.fraction) & !is.na(sfrac_y$subclonal.fraction) & sfrac_x$col!='black'
sfrac_x = sfrac_x[sel,]
sfrac_y = sfrac_y[sel,]

png(args$output, width = 500, height = 500)
if(is.null(args$snodes)) {
	plot(
		x=sfrac_x$subclonal.fraction, y=sfrac_y$subclonal.fraction, xlim=c(0,as.integer(args$max)), ylim=c(0,as.integer(args$max)),
		col=sfrac_x$col, pch=20, cex=.7, xlab=args$xlabel, ylab=args$ylabel, cex.lab=1.4
	)
	abline(v=c(0,1,2), col='grey')
	abline(h=c(0,1,2), col='grey')
} else {
	snodepos = read.table(args$snodes, sep="\t", comment.char = "", header=T, stringsAsFactors=F)
	layout(mat = matrix(c(2, 1, 0, 3), nrow=2, ncol=2), heights=c(1, 2), widths=c(2, 1))
	# Scatter-plot
	par(mar = c(5, 4, 0, 0))
	plot(
		x=sfrac_x$subclonal.fraction, y=sfrac_y$subclonal.fraction, xlim=c(0,as.integer(args$max)), ylim=c(0,as.integer(args$max)),
		col=sfrac_x$col, pch=20, cex=.7, xlab=args$xlabel, ylab=args$ylabel, cex.lab=1.4
	)
	abline(v=c(0,1,2), col='grey')
	abline(h=c(0,1,2), col='grey')
	# Density plot for the sample on x-axis
	par(mar = c(0, 4, 0, 0))
	d.x = density(sfrac_x[sfrac_x$vaf>=vaf_cutoff,'subclonal.fraction'], na.rm=T)
	plot(d.x, frame=F, main="", axes=F, ann=F, xlim=c(0,as.integer(args$max)))
	polygon(d.x, col="lightgrey", border="lightgrey")
	xnodes = snodepos[snodepos$sample==args$xsample,]
	for(i in 1:nrow(xnodes)) {
		xpos = d.x$x[d.x$x>(xnodes$location[i]-0.01) & d.x$x<(xnodes$location[i]+0.01)][1]
		segments(xpos, 0, xpos, d.x$y[d.x$x==xpos], col="darkgrey")
	}
	# Density plot for the sample on y-axis
	par(mar = c(5, 0, 0, 0))
	d.y = density(sfrac_y[sfrac_y$vaf>=vaf_cutoff,'subclonal.fraction'], na.rm=T)
	plot(d.y$y, d.y$x, xlim=range(d.y$y), ylim=c(0,as.integer(args$max)), frame=F, type="l", main="", axes=F, ann=F)
	polygon(d.y$y, d.y$x, col="lightgrey", border="lightgrey")
	ynodes = snodepos[snodepos$sample==args$ysample,]
	for(i in 1:nrow(ynodes)) {
		xpos = d.y$x[d.y$x>(ynodes$location[i]-0.01) & d.y$x<(ynodes$location[i]+0.01)][1]
		segments(0, xpos, d.y$y[d.y$x==xpos], xpos, col="darkgrey")
	}

}
dev.off()
