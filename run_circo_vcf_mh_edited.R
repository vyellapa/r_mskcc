#
# Plot the integrated Circos plot (inter-mutation distance for substitutions, genomic position for indels, Battenberg segments and rearrangements
#
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/integrated.circos/run_circo_vcf_mh.R SRCCT1 /ifs/res/papaemme/users/gg10/renal/sarcomatoid/SRCCT1 total

library(RCircos);
library(scales)

library(BSgenome.Hsapiens.UCSC.hg19) # human genome
library(VariantAnnotation) # 

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

source(paste(UTILS_DIR, 'generateHist.R', sep=''))
source(paste(UTILS_DIR, 'getMutTables.R', sep=''))
source(paste(UTILS_DIR, 'getPindelVcf.R', sep=''))

getVAF <- function(myFile, onlyPASSED=FALSE) {
	B35.1.2.vcf <- readVcf(myFile, "hg19")
	info(B35.1.2.vcf)
	#filters failed for each variant
	rd <- rowData(B35.1.2.vcf)
	fs <- fixed(B35.1.2.vcf )$FILTER
	fs.passed <- (fs=='PASS')
	if (onlyPASSED) {
	B35.1.2.vcf <- B35.1.2.vcf[fs.passed,]
	rd <- rowData(B35.1.2.vcf)
	fs <- fixed(B35.1.2.vcf )$FILTER
	fs.passed <- (fs=='PASS')
	}
	result <- as.data.frame(geno(B35.1.2.vcf)[['PM']])$TUMOUR
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args)>0) {
	SAMPLE.ID <- as.character(args[1]) # inputs
	ANALYSIS.BASE.DIR <- as.character(args[2])
	copy_number_type <- args[3]  # can be 'major' clone, 'minor' clone or 'total' copy number
	loadData <- TRUE
	forceHistogram <- FALSE
#	cPanelWidth <- as.numeric(as.character(args[length(args)]))
	cPanelWidth <- 0
	clinicalInfo <- FALSE
}

cat(paste( SAMPLE.ID))

# circos parameters 
params.my <- list()
params.my$track.background <- 'white'
params.my$highlight.width <- 0.2
params.my$point.size <- 0.15
params.my$point.type <- 19
params.my$radius.len <- 3
params.my$chr.ideog.pos <- 3.2
params.my$highlight.pos <- 3.35
params.my$chr.name.pos <- 3.45
params.my$plot.radius <- 3.55
params.my$track.in.start <- 3.05
params.my$track.out.start <- 3.2
params.my$text.size <-  0.6
params.my$track.padding <- c(0.07,  0.0, 0.07, 0.0,0.07, 100)
params.my$grid.line.color <- 'lightgrey'
params.my$track.heights <- c(0.85, 0.07, 0.07, 0.1, 0.1,  100)
params.my$track.height <- 0.1
params.my$sub.tracks <- 1
params.my$heatmap.cols <- c(
	alpha('lightcoral', 1), alpha('lightcoral', 0.5), alpha('lightgrey',0.10), alpha('olivedrab2', 0.3),
	alpha('olivedrab2', 0.5), alpha('olivedrab2',.7), alpha('olivedrab2', 0.75), alpha('olivedrab3', 0.9), alpha('olivedrab4', 0.9)
)
params.my$heatmap.ranges <- c(0,1,3,4,8,16, 32,64,1000)
# rearrangement links colors
INV.COLOR <- alpha('lightsteelblue1', 0.8)
DEL.COLOR <- alpha('lightcoral', 0.5)
DUPL.COLOR <-  alpha('olivedrab2', 0.75)
TRANSLOC.COLOR <- alpha('grey', 0.5)

cat(paste(SAMPLE.ID, '\n'))

#FILE.SUBS <-  list.files(path=paste(ANALYSIS.BASE.DIR, '/caveman', sep=''), pattern='*.caveman.muts.annot.vcf.gz$', full.names=T);
FILE.SUBS = list.files(path=ANALYSIS.BASE.DIR, pattern="*.caveman.muts.annot.vcf.gz$", recursive=T, full.names=T)

print(FILE.SUBS)
#FILE.REARR <- list.files(path=paste(ANALYSIS.BASE.DIR, '/bga', sep=''), pattern='.annot.bedpe$', full.names=T);


#### Changed here for it to accept the latest filtered SV file
FILE.REARR <- list.files(path=ANALYSIS.BASE.DIR, pattern='.sv.bedpe$', recursive=T, full.names=T);
print(FILE.REARR)
#FILE.CN <- list.files(path=paste(ANALYSIS.BASE.DIR, '/bbg', sep=''), pattern='.subclones.txt$', full.names=T); 
FILE.CN <- list.files(path=ANALYSIS.BASE.DIR, pattern='^E.*.subclones.txt$', recursive=T, full.names=T); 
print(FILE.CN)
FILE.INDELS <- list.files(path=ANALYSIS.BASE.DIR, pattern='.indeltype.output.txt$', recursive=T, full.names=T); 
print(FILE.INDELS)
SIGN.FILE <- '/lustre/scratch109/sanger/gg10/822_prostate/prostatectomies/mutational_signatures/integrated_circos_plots/822_mutations_singatures_ludmil.csv'    
DRIVER.FILE = ''
SIGN.CONTEX.FILE = '/ifs/work/leukgen/home/gg10/soft/integrated.circos/data/signatures_probabilities.txt'
if (loadData) {
	sh <- read.table(SIGN.CONTEX.FILE, header=TRUE, sep='\t')
	mut.order <- (sh[,'Somatic.Mutation.Type'])
	# indels
	indels <- read.table(FILE.INDELS, header=TRUE, sep='\t', quote = "")
	indels <- subset(indels, ACCESSION!='Y')
	ins <- indels[indels$VARIANT_TYPE=='I',]
	dels <- indels[indels$VARIANT_TYPE=='D' | indels$VARIANT_TYPE=='DI',]
	ins.formatted <- ins[,c('ACCESSION', 'MIN_POSITION', 'MAX_POSITION')]; names( ins.formatted) <- c('Chromosome','chromStart','chromEnd')
	dels.formatted <- dels[,c('ACCESSION', 'MIN_POSITION', 'MAX_POSITION')]; names( dels.formatted) <- c('Chromosome','chromStart','chromEnd')
	if (nrow(ins.formatted)>0) {    ins.formatted$Chromosome <- paste('chr',ins.formatted$Chromosome ,sep='')}
	if (nrow(dels.formatted)>0) {    dels.formatted$Chromosome <- paste('chr',dels.formatted$Chromosome ,sep='')}
	tile.cols <- vector()
	tile.cols[dels$Classification=='Microhomology-mediated'] <- 'firebrick4'
	tile.cols[dels$Classification=='Repeat-mediated'] <- 'firebrick1'
	tile.cols[dels$Classification=='None'] <- 'firebrick3'
	# rearrangements
	rearrs.formatted <- data.frame()
	if (file.exists(FILE.REARR)) {
		rearrs.formatted <- read.brassI(FILE.REARR)
	}
	# substitutions
	if (file.exists(FILE.SUBS)){
		subs.data <- getMutTables(FILE.SUBS, onlyPASSED=TRUE, onlyPassedAsrd=TRUE)
	}
	subs <- data.frame(chr=substring(subs.data$muts$chroms,4,100),
		position = subs.data$muts$starts,
		wt=subs.data$muts$wt,
		mt =  subs.data$muts$mt,
		ref_base_pyrimidine_context =  subs.data$muts$pyrwt,
		mutant_base_pyrimidine_context = subs.data$muts$pyrmut
	)
	subs <- subset(subs, chr!='Y')
	subs <- processSubs(subs)
	scatter.data <- subs$scatter.data 
	scatter.colors <- subs$scatter.colors 
	scatter.data.formatted <- data.frame(
		chromosome=scatter.data$chr, start=scatter.data$position , stop=scatter.data$position,
		logDistPrev=log10(scatter.data$distPrev), color=scatter.colors, mutType=scatter.data$mutType
	)
	cat(paste( dim(scatter.data)[1], ' subs \n'))
	# copy number
	cv.data <- data.frame()
	if (file.exists(FILE.CN)){
		cv.data <- read.table(FILE.CN, header=T, sep='\t')
		if(copy_number_type=='total' ) {
			cv.data = cbind(paste("chr",cv.data$chr, sep=""), cv.data[,c('startpos','endpos','BAF','LogR')])
		} else if (copy_number_type=='major') {
			cv.data = cbind(paste("chr",cv.data$chr, sep=""), cv.data[,c('startpos','endpos','BAF')], rowSums(cv.data[,c('nMaj1_A','nMin1_A')]))
		} else if (copy_number_type=='minor') {
			cv.data = cbind(paste("chr",cv.data$chr, sep=""), cv.data[,c('startpos','endpos','BAF')], rowSums(cv.data[,c('nMaj2_A','nMin2_A')]))
		}
		colnames(cv.data) = c('Chromosome','chromStart','chromEnd','BAF','LogR')
		cv.data$total.copy.number.inTumour = cv.data$LogR
		cv.data$total.copy.number.inNormal = rep(2, nrow(cv.data))
		cv.data$minor.copy.number.inNormal = rep(1, nrow(cv.data))
		cv.data = cv.data[,c('Chromosome','chromStart','chromEnd','total.copy.number.inNormal','minor.copy.number.inNormal','total.copy.number.inTumour','BAF','LogR')]

	} else if (file.exists(FILE.CN)) {
        cat("line 162\n")
	#	cv.data <- read.ascat(FILE.CN.NGS)
	}
	# signatures
	sign.file <- SIGN.FILE
	#sign.data <- try(suppressWarnings(read.table(SIGN.FILE, header=TRUE, check.names=FALSE)))
	#if (!inherits(sign.data, "try-error")) {
	#	sign.data <- read.table(sign.file, sep=',', header=TRUE, check.names=FALSE)
	#	sign.data[,'Sample Names'] <- as.character(sign.data[,'Sample Names'])
	#} else {
		sign.data = NULL
	#}
	driver.genes = NULL
}


#####################################################################################################################
png(file=paste(ANALYSIS.BASE.DIR, '/', SAMPLE.ID, "-integrated-circos.png",sep=''), height=4100, width=(5400*(1/(1-cPanelWidth))), res=550)
data(UCSC.HG19.Human.CytoBandIdeogram);
hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;
RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, tracks.inside=10, tracks.outside=1);
    
params <- RCircos.Get.Plot.Parameters();
params[names(params.my)] <- params.my
params$sub.tracks <- 1
RCircos.Reset.Plot.Parameters(params)
par(mar=c(0.001, 0.001, 0.001, 0.001))
par(fig=c(cPanelWidth,0.75*(1-cPanelWidth)+cPanelWidth,0,1))
 
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot.my();
title(paste(SAMPLE.ID, sep=''), line=-1);

#if (!is.null(driver.genes) & nrow(driver.genes)>0) {
#	RCircos.Gene.Connector.Plot.my(genomic.data=driver.genes, track.num=1, side="out"); # plot the driver genes
#	RCircos.Gene.Name.Plot(gene.data=driver.genes, name.col=4, track.num=2, side="out"); # plot the driver genes
#}
# substitutions
RCircos.Scatter.Plot.color(
	scatter.data=scatter.data.formatted , data.col=4, track.num=1,
	side="in", by.fold=0, scatter.colors = scatter.colors
);
cat('subs plotted \n')

params <- RCircos.Get.Plot.Parameters();
params$line.color <- 'white'
params$highlight.width <- 0.2
params$max.layers <- 5
params$tile.color <- 'darkgreen'
RCircos.Reset.Plot.Parameters(params)

op <- par(lwd = 0.1)	# square plotting region,
			# independent of device size
if (nrow(ins.formatted)>0) {
	my.RCircos.Tile.Plot(tile.data=ins.formatted, track.num=2, side="in");
}

params <- RCircos.Get.Plot.Parameters();
params$tile.color <- 'firebrick4'
RCircos.Reset.Plot.Parameters(params)
if (nrow(dels.formatted)>0) {
	my.RCircos.Tile.Plot(tile.data=dels.formatted, track.num=3, side="in", tile.colors=tile.cols);
}
cat('indels plotted \n')
    
## At end of plotting, reset to previous settings:
par(op)
# copy number
if (exists('cv.data') && (nrow(cv.data)>0)) {
#	heatmap.ranges.logr = seq(from = min(cv.data$LogR), to = max(cv.data$LogR), by = ((max(cv.data$LogR) - min(cv.data$LogR))/(256 - 1)))
	if(max(cv.data$LogR)>0) {
		logrfinite = cv.data$LogR[is.finite(cv.data$LogR)]
		heatmap.ranges.logr = c(seq(from=min(logrfinite), to=0 , by=((abs(min(logrfinite)) - 0)/127)), seq(from=0, to=max(logrfinite), by=((abs(max(logrfinite)) - 0)/128))[-1])
		heatmap.color.logr = colorRampPalette(c('blue','white','red'))(256)
	} else {
		heatmap.ranges.logr = seq(from=min(cv.data$LogR), to=0 , by=((abs(min(cv.data$LogR)) - 0)/128))
		heatmap.color.logr = colorRampPalette(c('blue','white'))(128)
	}
	heatmap.ranges.baf = seq(from = min(cv.data$BAF), to = max(cv.data$BAF), by = ((max(cv.data$BAF) - min(cv.data$BAF))/(256 - 1)))
	heatmap.color.baf = colorRampPalette(c('white','purple'))(256)
	RCircos.Heatmap.Plot.my(
		heatmap.data=cv.data, data.col=6, track.num=4, side="in",
		heatmap.ranges=heatmap.ranges.logr, heatmap.color=heatmap.color.logr
	)
	RCircos.Heatmap.Plot.my(
		heatmap.data=cv.data, data.col=7, track.num=5, side="in",
		heatmap.ranges=heatmap.ranges.baf, heatmap.color=heatmap.color.baf
	)
}
cat('copy number plotted \n')


# rearrangements
link.colors <- vector()
link.data <- rearrs.formatted
link.colors[link.data$pf==1 | link.data$pf==8] <- INV.COLOR
link.colors[link.data$pf==2] <- DEL.COLOR
link.colors[link.data$pf==4] <- DUPL.COLOR
link.colors[link.data$pf==32] <- TRANSLOC.COLOR

if (nrow( rearrs.formatted)>0) {
	RCircos.Link.Plot.my(link.data= rearrs.formatted, track.num=6, by.chromosome=TRUE, link.colors);
	cat('rearrangements plotted \n')
}

# # # # # # # # # # # # # # # # # # # # # 
# SIDE PLOTS
# Barplot/histogram for signatures
par(cex=0.6)
sign.row <- which(sign.data[,'Sample Names']==SAMPLE.ID)
if ((!forceHistogram) && (length(sign.row)==1)) { # if signatures for the sample were found
	par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.70, 0.95), new=TRUE) # subs and signatures
	subs.sign.data <- sign.data[sign.row, 2:ncol(sign.data)]
	subs.sign.lbs <- names(subs.sign.data)
	barplot.values <- rev(as.numeric(subs.sign.data))
	barplot.labels <- rev(subs.sign.lbs)
	mp <- barplot(
		barplot.values, main=paste(nrow(subs.data$muts),'substitutions of signatures:'), axes = FALSE,
		axisnames = FALSE, width=1, horiz=TRUE, border=NA, col=params.my$grid.line.color
	)
	axis(2, at=mp[barplot.values>=0], las=2, labels=barplot.labels[barplot.values>=0], col='grey', tick=FALSE, cex=0.5)
	axis(1, las=2, col='grey')
} else {
	# plot a histogram
	par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.70, 0.95), new=TRUE)
	names(subs.data$passed.hist) <- NA
	barplot(
		subs.data$passed.hist , col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16),
		rep('green2', 16), rep('hotpink',16)), main=paste(nrow(subs.data$muts),'substitutions'), border=NA, xaxt="n"
	)
	names(subs.data$passed.hist) <- mut.order
	axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
}
# Barplot for indels    
par(fig=c(cPanelWidth+0.8*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.49, 0.67), new=TRUE) # indels
indel.data <- c(
	sum(indels$VARIANT_TYPE=='D' & indels$Classification == 'Microhomology-mediated'),
	sum(indels$VARIANT_TYPE=='D' & indels$Classification == 'Repeat-mediated'),
	sum(indels$VARIANT_TYPE=='D' & indels$Classification == 'None'),
	sum(indels$VARIANT_TYPE=='I'),
	sum(indels$VARIANT_TYPE=='DI')
)
indel.col <- c('firebrick4', 'firebrick1', 'firebrick3', 'darkgreen', 'grey')
indel.lbs <- c('deletion \n m-homology',  'deletion repeat', 'deletion other', 'insertion', 'complex')
mp <- barplot(indel.data, main=paste(nrow(indels), 'deletions and insertions'), axes = FALSE, col=indel.col, axisnames = FALSE, width=1 , horiz=TRUE, border=NA)
axis(2, at = mp, las=2,  labels = indel.lbs, col='grey', tick=FALSE, cex=0.5)
axis(1, las=2, col='grey')
# Barplot for rearrangements
par(fig=c(cPanelWidth+0.76*(1-cPanelWidth),cPanelWidth+ 0.98*(1-cPanelWidth),0.26,0.46), new=TRUE) # rearrangements
rearrs.bar <- c(sum(rearrs.formatted$pf==32), sum(rearrs.formatted$pf==1)+sum(rearrs.formatted$pf==8), sum(rearrs.formatted$pf==2), sum(rearrs.formatted$pf==4))
rearrs.col <- c(TRANSLOC.COLOR, INV.COLOR, DEL.COLOR, DUPL.COLOR)


#### Changed here by so that they correctly represent the SV type
#rearrs.lbs <- c('translocation',  'inversion', 'deletion', 't. duplication')
rearrs.lbs <- c('translocation', 'deletion', 'inversion', 't. duplication')
mp <- barplot(rearrs.bar, main=paste(nrow(rearrs.formatted), 'rearrangements'), axes = FALSE, col=rearrs.col, axisnames = FALSE, width=1, horiz=TRUE, border=NA)
axis(2, at = mp[rearrs.bar>0], las=2,  labels =rearrs.lbs[rearrs.bar>0], col='grey', tick=FALSE, cex=0.5)
axis(1, las=2, col='grey')
# Plot subclonal structure
#par(fig=c(cPanelWidth+0.8*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.05, 0.24), new=TRUE) # density plot for VAF
#vafs = getVAF(FILE.SUBS, onlyPASSED=TRUE)
#sp = plot(density(vafs), main='Subclonal structure', axes=FALSE, col='darkgrey')
#axis(2, las=2, col='grey', cex=0.5)
#axis(1, las=2, col='grey', cex=0.5)
# Legend for substitutions
par(fig=c(cPanelWidth+0.06*(1-cPanelWidth),cPanelWidth + 0.15*(1-cPanelWidth),0.1,0.27), new=TRUE) # subs legend
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(
	x='center', legend=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'),
	col=c("royalblue", "black", "red", "grey", "green2", "hotpink"),
	pch=19, pt.cex=0.6, horiz=FALSE, bty='n'
)
# Legend for copy number tracks
par(fig=c(cPanelWidth+0.60*(1-cPanelWidth),cPanelWidth+0.70*(1-cPanelWidth),0.12,0.22), new=TRUE) # LogR/BAF legend
x1 <- seq(1,length(heatmap.ranges.logr),1)
y1 <- rep(1,length(heatmap.ranges.logr))
plot(c(), ylim=c(0,2), xlim=c(0,length(heatmap.ranges.logr)), axes=FALSE, xlab="", ylab="" , main='')
for(i in 0:length(x1)-1) {polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 0, 0), col=heatmap.color.logr[i+1], border=heatmap.color.logr[i+1])}  # First LogR
y1 <- rep(2,length(heatmap.ranges.baf))
for(i in 0:length(x1)-1) {polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 1, 1), col=heatmap.color.baf[i+1], border=heatmap.color.baf[i+1])}
axis(1, at=c(0,length(x1)), labels=c(round(min(cv.data$LogR), digits=2), round(max(cv.data$LogR), digits=2)), las=2, col='grey', cex=0.5)           # Then BAF
axis(3, at=c(0,length(x1)), labels=c(0.5, 1), las=2, col='grey', cex=0.5) 
axis(2, at=c(0.5,1.5), labels=c('LogR','BAF'), las=2, tick=FALSE, col='grey', cex=0.5)
dev.off()
