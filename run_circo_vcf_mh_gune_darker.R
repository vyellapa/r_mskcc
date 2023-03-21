#
# Plot the integrated Circos plot (inter-mutation distance for substitutions, genomic position for indels, Battenberg segments, rearrangements, mutational signaturesi, DP results)
# Usage:
# Rscript /ifs/work/leukgen/home/gg10/soft/integrated.circos/run_circo_vcf_mh.R
# tumour_name		Should be the same as the tumour name in the outout files from the pipeline
# tumour_analysis_root	Should contain at least these directories/files:
#			Directories 'caveman', 'pindel', 'brass' and 'bbg' is required.
#				Each directory should contain soft links to all the output files from the corresponding algorithm.
#			Directory 'subclonality' for Dirichlet process clustering is optional.
#			File 'driver_alterations.txt' is optional. This file is tab-delimited and in the following format:
#				Chromosome	chromStart	chromEnd	Gene
#				12	11802788	12048325	ETV6
#				15	88419988	88799962	NTRK3
#				9	21967751	21975132	CDKN2A
#			Directory 'mutational_signatures' is optional. If it does not exist, Will plot counts for 96 mutation types instead.
#				Should contain the output from deconstrucSigs.
# copy_number_type	One of these options:
#			major: Copy number for the major clone from Battenberg
#			minor: Copy number for the major clone from Battenberg
#			total: nTot and BAF columns from Battenberg
#			integer: Caveman copy number input from Battenberg (preferred)
# Example: Rscript /ifs/work/leukgen/home/gg10/soft/integrated.circos/run_circo_vcf_mh.R SRCCT1 /ifs/res/papaemme/users/gg10/renal/sarcomatoid/SRCCT1 integer
#
.libPaths(c("/home/gundemg/R_libs/","/home/gundemg/R/x86_64-pc-linux-gnu-library/3.2","/ifs/work/leukgen/bin/R/3.2.3/lib64/R/library", "/home/yellapav/R/x86_64-pc-linux-gnu-library/3.2"))
#print(.libPaths())
library(RCircos);
library(scales)

library(BSgenome.Hsapiens.UCSC.hg19) # human genome
library(VariantAnnotation) # 

library(DPClust)

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
source("/ifs/work/leukgen/home/gg10/soft/RunDP.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)>0) {
	SAMPLE.ID <- as.character(args[1]) # inputs
	ANALYSIS.BASE.DIR <- as.character(args[2])
	copy_number_type <- args[3]
	loadData <- TRUE
	forceHistogram <- 'none'
	imageType = 'png'
	if(length(args)>=4) {
		if(args[4]=="pdf") {imageType = 'png'}
	}
	if(length(args)>=5) {
		if(args[5]=="96nt") {forceHistogram = '96nt'}
		if(args[5]=="6nt") {forceHistogram = '6nt'}
		if(args[5]=="none") {forceHistogram = 'none'}
	}
	titleOn <- TRUE
	if(length(args)>=6) {
		if(args[6]=="FALSE") {titleOn = FALSE}
	}
	cn.max = 4
	if(length(args)>=7) {
		cn.max = as.integer(args[7])
	}
	driver.genes.track.num = -2
	if(length(args)>=8) {
		driver.genes.track.num = as.integer(args[8])
	}
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
params.my$text.size <-  0.7
params.my$track.padding <- c(0.07,  0.0, 0.07, 0.0,0.07, 100)
params.my$grid.line.color <- '#737373'
params.my$track.heights <- c(0.85, 0.07, 0.07, 0.1, 0.1,  100)
params.my$track.height <- 0.1
params.my$sub.tracks <- 1
params.my$heatmap.cols <- c(
	alpha('lightcoral', 1), alpha('lightcoral', 0.5), alpha('lightgrey',0.10), alpha('olivedrab2', 0.3),
	alpha('olivedrab2', 0.5), alpha('olivedrab2',.7), alpha('olivedrab2', 0.75), alpha('olivedrab3', 0.9), alpha('olivedrab4', 0.9)
)
params.my$heatmap.ranges <- c(0,1,3,4,8,16, 32,64,1000)
# rearrangement links colors
INV.COLOR <- alpha('#33a02c', 0.9)
DEL.COLOR <- alpha('#1f78b4', 0.9)
DUPL.COLOR <-  alpha('#e31a1c', 0.9)
TRANSLOC.COLOR <- alpha('#737373', 0.75)
# DP-clustering parameters
is.male = T
is.vcf = F

cat(paste(SAMPLE.ID, '\n'))
### Input files 
FILE.SUBS <-  list.files(path=paste(ANALYSIS.BASE.DIR, '/caveman', sep=''), pattern='*.caveman.muts.annot.vcf.gz$', full.names=T);
if(file.exists(FILE.SUBS)) {
	print(paste('Found subs file: ', FILE.SUBS, sep=''))
} else {
	print(paste("Couldn't find subs file : ", FILE.SUBS, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
}
FILE.REARR <- list.files(path=paste(ANALYSIS.BASE.DIR, '/brass', sep=''), pattern='.annot.bedpe$', full.names=T);
if(file.exists(FILE.REARR)) {
	print(paste('Found rearrangements file: ', FILE.REARR, sep=''))
} else {
	print(paste("Couldn't find rearrangements file : ", FILE.REARR, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
}
FILE.SUBCLONES <- list.files(path=paste(ANALYSIS.BASE.DIR, '/bbg', sep=''), pattern='.subclones.txt$', full.names=T);
if(file.exists(FILE.SUBCLONES)) {
	print(paste('Found BBG subclones file: ', FILE.SUBCLONES, sep=''))
} else {
	print(paste("Couldn't find BBG subclones file : ", FILE.SUBCLONES, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
} 
FILE.INTEGER.CN <- list.files(path=paste(ANALYSIS.BASE.DIR, '/bbg/', sep=''), pattern='_battenberg_cn.vcf.gz$', full.names=T); 
if(file.exists(FILE.INTEGER.CN)) {
	print(paste('Found integer CN file: ', FILE.INTEGER.CN, sep=''))
} else {
	print(paste("Couldn't find integer CN file : ", FILE.INTEGER.CN, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
}
FILE.INDELS <- list.files(path=paste(ANALYSIS.BASE.DIR, '/pindel', sep=''), pattern='.indeltype.output.txt$', full.names=T);
if(file.exists(FILE.INDELS)) {
	print(paste('Found indels file: ', FILE.INDELS, sep=''))
} else {
	print(paste("Couldn't find indels file : ", FILE.INDELS, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
}
FILE.SIGNATURES <- paste(ANALYSIS.BASE.DIR, '/mutational_signatures/deconstrucsigs_summary.csv', sep='')
if(file.exists(FILE.SIGNATURES)) {
	print(paste('Found signature file: ', FILE.SIGNATURES, sep=''))
} else {
	print(paste("Couldn't signature file : ", FILE.SIGNATURES, sep=''))
}
FILE.DRIVER.GENES = paste(ANALYSIS.BASE.DIR, '/driver_alterations.txt', sep='')
if(file.exists(FILE.DRIVER.GENES)) {
	print(paste('Found drivers file: ', FILE.DRIVER.GENES, sep=''))
} else {
	print(paste("Couldn't drivers file : ", FILE.DRIVER.GENES, sep=''))
}
FILE.RHO.PSI = list.files(path=paste(ANALYSIS.BASE.DIR, '/bbg', sep=''), pattern='*_rho_and_psi.txt$', full.names=T)
if(file.exists(FILE.RHO.PSI)) {
	print(paste('Found BBG rho-and-psi file: ', FILE.RHO.PSI, sep=''))
} else {
	print(paste("Couldn't find BBG rho-and-psi file : ", FILE.RHO.PSI, sep=''))
	quit(save = "no", status=0, runLast=FALSE)
}
FILE.DP.INPUT = list.files(path=paste(ANALYSIS.BASE.DIR, '/subclonality', sep=''), pattern='*_allDirichletProcessInfo.txt$', full.names=T)
if(length(FILE.DP.INPUT)>0 && file.exists(FILE.DP.INPUT)) {
	print(paste('Found DP input file: ', FILE.DP.INPUT, sep=''))
} else {
	print(paste("Couldn't DP input file : ", FILE.DP.INPUT, sep=''))
}
FILE.DP.DENSITY = list.files(path=paste(ANALYSIS.BASE.DIR, '/subclonality', sep=''), pattern='*_DirichletProcessplotdensity.txt', recursive=T, full.names=T)
if(length(FILE.DP.DENSITY)>0 && file.exists(FILE.DP.DENSITY)) {
	print(paste('Found DP density file: ', FILE.DP.DENSITY, sep=''))
} else {
	print(paste("Couldn't DP density file : ", FILE.DP.DENSITY, sep=''))
}
FILE.DP.POLYGON.DATA = list.files(path=paste(ANALYSIS.BASE.DIR, '/subclonality', sep=''), pattern='*_DirichletProcessplotpolygonData.txt', recursive=T, full.names=T)
if(length(FILE.DP.POLYGON.DATA)>0 && file.exists(FILE.DP.POLYGON.DATA)) {
	print(paste('Found DP polygon data file: ', FILE.DP.POLYGON.DATA, sep=''))
	dp.polygon.data = read.table(FILE.DP.POLYGON.DATA, header=T)
} else {
	print(paste("Couldn't DP polygon data file : ", FILE.DP.POLYGON.DATA, sep=''))
}
SIGN.CONTEX.FILE = '/ifs/work/leukgen/home/gg10/soft/integrated.circos/data/signatures_probabilities.txt'

## Load all input files
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
	if (nrow(ins.formatted)>0) {ins.formatted$Chromosome <- paste('chr',ins.formatted$Chromosome ,sep='')}
	if (nrow(dels.formatted)>0) {dels.formatted$Chromosome <- paste('chr',dels.formatted$Chromosome ,sep='')}
	tile.cols <- vector()
	tile.cols[dels$Classification=='Microhomology-mediated'] <- 'firebrick4'
	tile.cols[dels$Classification=='Repeat-mediated'] <- 'firebrick1'
	tile.cols[dels$Classification=='None'] <- 'firebrick3'
	# rearrangements
	rearrs.formatted <- data.frame()
	if (file.exists(FILE.REARR)) {
		rearrs.formatted <- read.brassII(FILE.REARR)
	}
	# substitutions
        ## Changed the ASRD filter here
	if (file.exists(FILE.SUBS)){
		subs.data <- getMutTables(FILE.SUBS, onlyPASSED=TRUE)
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
	if (copy_number_type!="integer" & file.exists(FILE.SUBCLONES)){
		cv.data <- read.table(FILE.SUBCLONES, header=T, sep='\t')
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

	} else if (copy_number_type=="integer" && file.exists(FILE.INTEGER.CN)) {
		cmd = paste(
			'zless ', FILE.INTEGER.CN,
			' | grep -v "#" | awk -F"\\t" \'BEGIN{OFS=","; ln=1}{split($8,a,"END="); split($10,n,":"); split($11,t,":"); ',
			'print ln,$1,$2,a[2],n[2], n[3], t[2]+t[3],t[3]; ln=ln+1}\' > bbg.ngs.tmp.txt', sep=''
		)
		system(cmd)
		cv.data <- read.ascat('bbg.ngs.tmp.txt')
		cv.data[cv.data$major.copy.number.inTumour>=cn.max,'major.copy.number.inTumour'] = cn.max
		cv.data[cv.data$minor.copy.number.inTumour>=cn.max,'minor.copy.number.inTumour'] = cn.max
		system('rm bbg.ngs.tmp.txt')
	}
	# signatures
	sign.data = NULL
	if (file.exists(FILE.SIGNATURES)) {
		sign.data <- read.table(FILE.SIGNATURES, sep=',', header=TRUE, check.names=FALSE)
		sign.data[,'Sample Names'] <- as.character(sign.data[,'Sample Names'])
	} else {
		cat("No signature results was found.\n")
	}
	# driver alterations
	driver.genes = NULL
	dp.driver.genes = NULL
	if(file.exists(FILE.DRIVER.GENES)) {
		driver.genes = read.table(file=FILE.DRIVER.GENES, header=T, sep="\t", stringsAsFactors=F)
		if(nrow(driver.genes[driver.genes$AnnotateDP=='yes',])>0) {
			dp.driver.genes = driver.genes[driver.genes$AnnotateDP=='yes',]
		}
		cat(paste( dim(driver.genes)[1], ' driver alterations \n'))
	} else {
		cat("No driver file was found.\n")
	}
	# DP clustering results
	cellularity = NULL
	ploidy = NULL
	if(file.exists(FILE.RHO.PSI)) {
		cellularity = read.table(FILE.RHO.PSI, header=T)['FRAC_GENOME','rho']
		ploidy = read.table(FILE.RHO.PSI, header=T)['FRAC_GENOME','psi']
	}
	dp.datafile = NULL
	if(length(FILE.DP.INPUT)>0 && file.exists(FILE.DP.INPUT)) {
		dp.datafile = FILE.DP.INPUT
	}
	dp.density = NULL
	if(length(FILE.DP.DENSITY)>0 && file.exists(FILE.DP.DENSITY)) {
		dp.density = read.table(FILE.DP.DENSITY, header=T)
	}
	dp.polygon.data = NULL
	if(length(FILE.DP.POLYGON.DATA)>0 && file.exists(FILE.DP.POLYGON.DATA)) {
		dp.polygon.data = read.table(FILE.DP.POLYGON.DATA, header=T)
	}
	dp.dataset = NULL
	if(!is.null(cellularity) & !is.null(dp.datafile) & !is.null(dp.density) & !is.null(dp.polygon.data)) {
		dp.dataset = load.data(
			dp.datafile, cellularity=cellularity, Chromosome="chr", position="end", WT.count="WT.count",
			mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut",
			mutation.copy.number="mutation.copy.number", subclonal.fraction="subclonal.fraction",
			is.male=is.male, is.vcf=is.vcf, ref.genome.version="hg19"
		)
	}
}

#####################################################################################################################
if(imageType=='png') {
	png(file=paste(ANALYSIS.BASE.DIR, '/', SAMPLE.ID, "-integrated-circo.png",sep=''), height=4100, width=(5400*(1/(1-cPanelWidth))), res=550)
} else if (imageType=='pdf') {
	## changed as image was too big pdf(file=paste(ANALYSIS.BASE.DIR, '/', SAMPLE.ID, "-integrated-circos.pdf",sep=''), height=410, width=(540*(1/(1-cPanelWidth))))
	pdf(file=paste(ANALYSIS.BASE.DIR, '/', SAMPLE.ID, "-integrated-circo.pdf",sep=''), height=16, width=20)
}
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
title.circos = SAMPLE.ID
if(titleOn & !is.null(cellularity) && !is.null(ploidy)) {
	title.circos = paste(SAMPLE.ID, '\npurity=', round(as.numeric(cellularity), digits=2), ' ploidy=', round(as.numeric(ploidy), digits=2), sep='')
} else if (!is.null(cellularity) && !is.null(ploidy)) {
	title.circos = paste('purity=', round(as.numeric(cellularity), digits=2), ' ploidy=', round(as.numeric(ploidy), digits=2), sep='')
}
title(title.circos, line=-2);

if (!is.null(driver.genes)) {
#	RCircos.Gene.Connector.Plot.my(genomic.data=driver.genes, track.num=1, side="out"); # plot the driver genes
	RCircos.Gene.Name.Plot.my2(gene.data=driver.genes, name.col=4, track.num=driver.genes.track.num, side="out"); # plot the driver genes
}
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
if (exists('cv.data') & copy_number_type!="integer" & nrow(cv.data)>0) {
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
} else if (exists('cv.data') & copy_number_type=="integer" & nrow(cv.data)>0) {
	maxcn = max(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
	mincn = min(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
	if(mincn==1) {mincn=0}
	heatmap.ranges = seq(from=mincn, to=maxcn, by=(1/128))
	heatmap.color.minorcn = colorRampPalette(c('white','deeppink4'))(128*(maxcn-mincn)+1)
	heatmap.color.majorcn = colorRampPalette(c('white','darkorange4'))(128*(maxcn-mincn)+1)
	RCircos.Heatmap.Plot.my(heatmap.data=cv.data, data.col=7, track.num=4, side="in", heatmap.ranges=heatmap.ranges, heatmap.color=heatmap.color.minorcn)
	RCircos.Heatmap.Plot.my(heatmap.data=cv.data, data.col=8, track.num=5, side="in", heatmap.ranges=heatmap.ranges, heatmap.color=heatmap.color.majorcn)
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
##########################
####### SIDE PLOTS #######
# Barplot/histogram for signatures
par(cex=0.6)
sign.row <- which(sign.data[,'Sample Names']==SAMPLE.ID)
if (forceHistogram=='none' && (length(sign.row)==1)) { # if signatures for the sample were found
	par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.70, 0.95), new=TRUE) # subs and signatures
	sign.data = sign.data[,order(sign.data, decreasing=T)[1:10]]
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
} else if (forceHistogram=='96nt') {
	# plot a histogram
	par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.70, 0.95), new=TRUE)
	names(subs.data$passed.hist) <- NA
	barplot(
		subs.data$passed.hist , col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16),
		rep('green2', 16), rep('hotpink',16)), main=paste(nrow(subs.data$muts),'substitutions'), border=NA, xaxt="n"
	)
	names(subs.data$passed.hist) <- mut.order
	axis(1, at=c(10, 29, 48, 68, 86, 105), labels=c('C>A','C>G','C>T','T>A','T>C','T>G'), las=2)
} else if (forceHistogram=='6nt') {
	b = c(
		sum(subs.data$passed.hist[grep('C>A', names(subs.data$passed.hist), value=T)]),
		sum(subs.data$passed.hist[grep('C>G', names(subs.data$passed.hist), value=T)]),
		sum(subs.data$passed.hist[grep('C>T', names(subs.data$passed.hist), value=T)]),
		sum(subs.data$passed.hist[grep('T>A', names(subs.data$passed.hist), value=T)]),
		sum(subs.data$passed.hist[grep('T>C', names(subs.data$passed.hist), value=T)]),
		sum(subs.data$passed.hist[grep('T>G', names(subs.data$passed.hist), value=T)])
	)
	names(b) = c('C>A','C>G','C>T','T>A','T>C','T>G')
	par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.70, 0.95), new=TRUE)
	barplot(b, col=c('royalblue','black','red','grey','green2','hotpink'), las=2, border=NA, main=paste(nrow(subs.data$muts),'substitutions'))
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
par(fig=c(cPanelWidth+0.76*(1-cPanelWidth),cPanelWidth+ 0.98*(1-cPanelWidth),0.30,0.46), new=TRUE) # rearrangements
rearrs.bar <- c(sum(rearrs.formatted$pf==32), sum(rearrs.formatted$pf==1)+sum(rearrs.formatted$pf==8), sum(rearrs.formatted$pf==2), sum(rearrs.formatted$pf==4))
rearrs.col <- c(TRANSLOC.COLOR, INV.COLOR, DEL.COLOR, DUPL.COLOR)
rearrs.lbs <- c('translocation',  'inversion', 'deletion', 't. duplication')
mp <- barplot(rearrs.bar, main=paste(nrow(rearrs.formatted), 'rearrangements'), axes = FALSE, col=rearrs.col, axisnames = FALSE, width=1, horiz=TRUE, border=NA)
axis(2, at = mp[rearrs.bar>0], las=2,  labels =rearrs.lbs[rearrs.bar>0], col='grey', tick=FALSE, cex=0.5)
axis(1, las=2, col='grey')
# Plot subclonal structure
if(!is.null(dp.dataset)) {
	par(fig=c(cPanelWidth+0.76*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.05, 0.28), new=TRUE) # density plot for corrected VAF
	plot1D_DP(
		dp.density, dp.polygon.data[,1], pngFile=NA, dp.driver.genes, density.from=0, x.max=1.5, #y.max=6,
		mutationCopyNumber=dp.dataset$mutation.copy.number, no.chrs.bearing.mut=dp.dataset$copyNumberAdjustment, samplename=SAMPLE.ID
	)
}
#######################
####### LEGENDS #######
# Legend for substitutions
par(fig=c(cPanelWidth+0.02*(1-cPanelWidth),cPanelWidth + 0.15*(1-cPanelWidth),0.05,0.27), new=TRUE) # subs legend
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(
	x='center', legend=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'),
	col=c("royalblue", "black", "red", "grey", "green2", "hotpink"),
	pch=19, pt.cex=0.6, horiz=FALSE, bty='n'
)
# Legend for copy number tracks
if(copy_number_type!="integer") {
	par(fig=c(cPanelWidth+0.60*(1-cPanelWidth),cPanelWidth+0.70*(1-cPanelWidth),0.12,0.22), new=TRUE) # LogR/BAF legend
	x1 <- seq(1,length(heatmap.ranges.logr),1)
	y1 <- rep(1,length(heatmap.ranges.logr))
	plot(c(), ylim=c(0,2), xlim=c(0,length(heatmap.ranges.logr)), axes=FALSE, xlab="", ylab="" , main='')
	# First LogR
	for(i in 0:length(x1)-1) {
		polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 0, 0), col=heatmap.color.logr[i+1], border=heatmap.color.logr[i+1])
	}
	y1 <- rep(2,length(heatmap.ranges.baf))
	# Then BAF
	for(i in 0:length(x1)-1) {
		polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 1, 1), col=heatmap.color.baf[i+1], border=heatmap.color.baf[i+1])
	}
	axis(1, at=c(0,length(x1)), labels=c(round(min(cv.data$LogR), digits=2), round(max(cv.data$LogR), digits=2)), las=2, col='grey', cex=0.5)
	axis(3, at=c(0,length(x1)), labels=c(0.5, 1), las=2, col='grey', cex=0.5) 
	axis(2, at=c(0.5,1.5), labels=c('LogR','BAF'), las=2, tick=FALSE, col='grey', cex=0.5)
} else if(copy_number_type=="integer") {
	par(fig=c(cPanelWidth+0.60*(1-cPanelWidth),cPanelWidth+0.70*(1-cPanelWidth),0.12,0.22), new=TRUE) # Major/minor allele legend
	x1 <- seq(1,length(heatmap.ranges),1)
	y1 <- rep(1,length(heatmap.ranges))
	plot(c(), ylim=c(0,2), xlim=c(0,length(heatmap.ranges)), axes=FALSE, xlab="", ylab="" , main='')
	# First minor CN
	for(i in 0:length(x1)-1) {
		polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 0, 0), col=heatmap.color.minorcn[i+1], border=heatmap.color.minorcn[i+1])
	}
	y1 <- rep(2,length(heatmap.ranges))
	# Then major CN
	for(i in 0:length(x1)-1) {
		polygon(c(i,x1[i+1],x1[i+1],i), c(y1[i+1], y1[i+1], 1, 1), col=heatmap.color.majorcn[i+1], border=heatmap.color.majorcn[i+1])
	}
	axis.labels = seq(from=round(mincn, digits=2), to=round(maxcn, digits=2), by=1)
	if(length(axis.labels)>4) {
		axis(1, at=c(0,length(x1)), labels=c(round(mincn, digits=2), round(maxcn, digits=2)), las=2, col='grey', cex=0.5)
		axis(3, at=c(0,length(x1)), labels=c(round(mincn, digits=2), round(maxcn, digits=2)), las=2, col='grey', cex=0.5)
		axis(2, at=c(0.5,1.5), labels=c('Minor allele','Major allele'), las=2, tick=FALSE, col='grey', cex=0.5)
	} else {
		axis(1, at=x1[heatmap.ranges %in% axis.labels]-1, labels=axis.labels, las=2, col='grey', cex=0.5)
		axis(3, at=x1[heatmap.ranges %in% axis.labels]-1, labels=axis.labels, las=2, col='grey', cex=0.5)
		axis(2, at=c(0.5,1.5), labels=c('Minor allele','Major allele'), las=2, tick=FALSE, col='grey', cex=0.5)
	}
}
dev.off()
