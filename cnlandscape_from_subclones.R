plot_genomewide_cn_landscape <- function(all.subclones, pdf_name) {
	sGain.segs = NULL
	cGain.segs = NULL
	sLoss.segs = NULL
	cLoss.segs = NULL
	highGain = NULL
	hGain.segs = NULL
	dLoss.segs = NULL
	for(chr in unique(all.subclones$chr)) {
		datasub = all.subclones[all.subclones$chr==chr,]
		all.positions = sort(unique(c(datasub$startpos, datasub$endpos)))
		all.positions = cbind(all.positions[1:length(all.positions)-1], all.positions[2:length(all.positions)])
		colnames(all.positions) = c('start','end')
		for(i in 1:nrow(all.positions)) {
			sel = datasub$startpos<=all.positions[i,'start'] & datasub$endpos>=all.positions[i,'end']
			clonal.cn = datasub[datasub$frac1_A==1 & sel,c('sample','nMaj1_A','nMin1_A','ntot')]
			n.gain = 0
			n.loss = 0
			if (nrow(clonal.cn)>0) {
				n.gain = length(which(rowSums(clonal.cn[,c('nMaj1_A','nMin1_A')]) > sample.info[clonal.cn$sample,'ploidy']))
	#			n.loss = length(which(rowSums(clonal.cn[,c('nMaj1_A','nMin1_A')]) < 2))
				n.loss = length(which(clonal.cn$nMin1_A==0))
			}
			cGain.segs = rbind(cGain.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.gain))
			cLoss.segs = rbind(cLoss.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.loss))
			subclonal.cn = datasub[datasub$frac1_A<1 & sel,c('sample','nMaj1_A','nMin1_A','ntot')]
			n.gain = 0
			n.loss = 0
			if (nrow(subclonal.cn)>0) {
				n.gain = length(which(rowSums(subclonal.cn[,c('nMaj1_A','nMin1_A')]) > sample.info[subclonal.cn$sample,'ploidy']))
	#			n.loss = length(which(rowSums(subclonal.cn[,c('nMaj1_A','nMin1_A')]) < 2))
				n.loss = length(which(subclonal.cn$nMin1_A==0))
			}
			sGain.segs = rbind(sGain.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.gain))
			sLoss.segs = rbind(sLoss.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.loss))
			n.high.gain = length(which(datasub[sel,'ntot']>8))
			hGain.segs = rbind(hGain.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.high.gain))
			n.deep.loss = length(which(datasub[sel,'ntot']<1))
			dLoss.segs = rbind(dLoss.segs, c(chr, all.positions[i,'start'], all.positions[i,'end'],n.deep.loss))
			
		}
	}
	colnames(cGain.segs) = c('chr','sp','ep','val')
	cGain.segs = as.data.frame(cGain.segs)
	cGain.segs$ep = as.integer(as.character(cGain.segs$ep))
	cGain.segs$sp = as.integer(as.character(cGain.segs$sp))
	cGain.segs$val = as.integer(as.character(cGain.segs$val))
	colnames(sGain.segs) = c('chr','sp','ep','val')
	sGain.segs = as.data.frame(sGain.segs)
	sGain.segs$ep = as.integer(as.character(sGain.segs$ep))
	sGain.segs$sp = as.integer(as.character(sGain.segs$sp))
	sGain.segs$val = as.integer(as.character(sGain.segs$val))
	colnames(cLoss.segs) = c('chr','sp','ep','val')
	cLoss.segs = as.data.frame(cLoss.segs)
	cLoss.segs$ep = as.integer(as.character(cLoss.segs$ep))
	cLoss.segs$sp = as.integer(as.character(cLoss.segs$sp))
	cLoss.segs$val = as.integer(as.character(cLoss.segs$val))
	colnames(sLoss.segs) = c('chr','sp','ep','val')
	sLoss.segs = as.data.frame(sLoss.segs)
	sLoss.segs$ep = as.integer(as.character(sLoss.segs$ep))
	sLoss.segs$sp = as.integer(as.character(sLoss.segs$sp))
	sLoss.segs$val = as.integer(as.character(sLoss.segs$val))
	colnames(dLoss.segs) = c('chr','sp','ep','val')
	dLoss.segs = as.data.frame(dLoss.segs)
	dLoss.segs$ep = as.integer(as.character(dLoss.segs$ep))
	dLoss.segs$sp = as.integer(as.character(dLoss.segs$sp))
	dLoss.segs$val = as.integer(as.character(dLoss.segs$val))
	colnames(hGain.segs) = c('chr','sp','ep','val')
	hGain.segs = as.data.frame(hGain.segs)
	hGain.segs$ep = as.integer(as.character(hGain.segs$ep))
	hGain.segs$sp = as.integer(as.character(hGain.segs$sp))
	hGain.segs$val = as.integer(as.character(hGain.segs$val))
	# Plot
	chrs = read.table(file='~/Documents/projects/gr37.fasta.fai', header=F, stringsAsFactors=F)
	maxchr = chrs$V2
	names(maxchr) = chrs$V1
	loci = rep(0,23)
	names(loci) = names(maxchr[1:23])
	for (i in 1:(length(loci)-1)) {
		chr = names(loci)[i]
		loci[i+1] = loci[i] + maxchr[chr]
	}
	chrg = NULL
	for (i in 2:length(loci)) {
		chrg = c(chrg, (loci[i]+loci[i-1])/2)
	}
	pdf(pdf_name, width=18, height=12)
	par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
	clonal.data = cGain.segs[cGain.segs$chr!='X',]
	clonal.data$chrlen = rep(0,nrow(clonal.data))
	for(chr in unique(clonal.data$chr)) {clonal.data[clonal.data$chr==chr,'chrlen'] = loci[chr]}
	clonal.data[,'sp'] = clonal.data[,'sp'] + clonal.data[,'chrlen']
	clonal.data[,'ep'] = clonal.data[,'ep'] + clonal.data[,'chrlen']
	subclonal.data = sGain.segs[sGain.segs$chr!='X',]
	subclonal.data$chrlen = rep(0,nrow(subclonal.data))
	for(chr in unique(subclonal.data$chr)) {subclonal.data[subclonal.data$chr==chr,'chrlen'] = loci[chr]}
	subclonal.data[,'sp'] = subclonal.data[,'sp'] + subclonal.data[,'chrlen']
	subclonal.data[,'ep'] = subclonal.data[,'ep'] + subclonal.data[,'chrlen']
	plot(c(), ylim=c(-0.5,1), xlim=c(1,max(clonal.data$ep)), bty="n", xlab="", ylab="", xaxt="n")
	segments(x0=clonal.data[,'sp'], x1=clonal.data[,'ep'], y0=clonal.data[,'val']/132, lwd=2, col='orangered3')
	segments(x0=subclonal.data[,'sp'], x1=subclonal.data[,'ep'], y0=subclonal.data[,'val']/132, lwd=2, col='orange')

	clonal.data = cLoss.segs[cLoss.segs$chr!='X',]
	clonal.data$chrlen = rep(0,nrow(clonal.data))
	for(chr in unique(clonal.data$chr)) {clonal.data[clonal.data$chr==chr,'chrlen'] = loci[chr]}
	clonal.data[,'sp'] = clonal.data[,'sp'] + clonal.data[,'chrlen']
	clonal.data[,'ep'] = clonal.data[,'ep'] + clonal.data[,'chrlen']
	subclonal.data = sLoss.segs[sLoss.segs$chr!='X',]
	subclonal.data$chrlen = rep(0,nrow(subclonal.data))
	for(chr in unique(subclonal.data$chr)) {subclonal.data[subclonal.data$chr==chr,'chrlen'] = loci[chr]}
	subclonal.data[,'sp'] = subclonal.data[,'sp'] + subclonal.data[,'chrlen']
	subclonal.data[,'ep'] = subclonal.data[,'ep'] + subclonal.data[,'chrlen']
	segments(x0=clonal.data[,'sp'], x1=clonal.data[,'ep'], y0=(clonal.data[,'val']/(-132)), lwd=2, col='navyblue')
	segments(x0=subclonal.data[,'sp'], x1=subclonal.data[,'ep'], y0=(subclonal.data[,'val']/(-132)), lwd=2, col='lightsteelblue3')
	abline(v=loci, lty = 2)
	text(x = chrg, y = rep(1, 22), labels=1:22)
	legend(0,0.9, horiz=T, lty=c(1,1), lwd=c(2.5,2.5), c('clonal gain','subclonal gain','clonal loss','subclonal loss'), col = c('orangered3','orange','navyblue','lightsteelblue3'), bg="white")
	### High level amplification and deep deletions
	highampdata = hGain.segs[cLoss.segs$chr!='X',]
	highampdata$chrlen = rep(0,nrow(highampdata))
	for(chr in unique(highampdata$chr)) {highampdata[highampdata$chr==chr,'chrlen'] = loci[chr]}
	highampdata[,'sp'] = highampdata[,'sp'] + highampdata[,'chrlen']
	highampdata[,'ep'] = highampdata[,'ep'] + highampdata[,'chrlen']
	deeplossdata = dLoss.segs[cLoss.segs$chr!='X',]
	deeplossdata$chrlen = rep(0,nrow(deeplossdata))
	for(chr in unique(deeplossdata$chr)) {deeplossdata[deeplossdata$chr==chr,'chrlen'] = loci[chr]}
	deeplossdata[,'sp'] = deeplossdata[,'sp'] + deeplossdata[,'chrlen']
	deeplossdata[,'ep'] = deeplossdata[,'ep'] + deeplossdata[,'chrlen']

	plot(c(), ylim=c(-0.1,0.3), xlim=c(1,max(highampdata$ep)), bty="n", xlab="", ylab="", xaxt="n")
	for(i in 1:nrow(highampdata)) {segments(x0=highampdata[i,'sp'], x1=highampdata[i,'ep'], y0=highampdata[i,'val']/132, lwd=2, col='orangered3')}
	for(i in 1:nrow(deeplossdata)) {segments(x0=deeplossdata[i,'sp'], x1=deeplossdata[i,'ep'], y0=(-1)*deeplossdata[i,'val']/132, lwd=2, col='navyblue')}
	abline(v=loci, lty = 2)
	text(x = chrg, y = rep(1, 22), labels=1:22)
	legend(0,0.3, horiz=T, lty=c(1,1), lwd=c(2.5,2.5), c('high amplification','deep loss'), col = c('orangered2','navyblue'), bg="white")
	dev.off()
}
#########################################
bands = read.table('~/Documents/projects/ucsc_goldenpath_hg19_cytoband.txt', header=F, sep="\t")
chrs = read.table('~/Documents/projects/gr37.fasta.fai', header=F, sep="\t")
sample.info = read.table(file='subclonality_analysis/purity_ploidy_estimates.txt', header=T)
rownames(sample.info) = sample.info$sample
all.files = list.files('subclonality_analysis/subclones', pattern='*subclones.txt', full.names=T)
all.subclones = NULL
for(file in all.files) {
	s = read.table(file, sep="\t", header=T, row.names=1)
	s$sample = unlist(strsplit(basename(file), split="_"))[1]
	all.subclones = rbind(all.subclones, s)
}
all.subclones$chr = as.character(all.subclones$chr)
#########################################
####### ANNOTATED GENOMIC EVENTS ########
#########################################
annotations = read.table(file='annotated_genomics_events.txt', header=T, sep="\t", stringsAsFactors=F)
#### chr1p loss events
sel.samples = annotations[grepl('^LOH_p', annotations$chr1p_events),'sample']
subclones.sub = cbind(all.subclones[all.subclones$sample %in% sel.samples & all.subclones$chr==1 & all.subclones$endpos<=min(bands[bands$V1=='chr1' & bands$V5=='acen','V2']) & all.subclones$nMin1_A==0,c('sample','chr','startpos','endpos')],'darkgreen')
manual.segs = cbind(all.subclones[all.subclones$sample %in% c('TARGET-30-PASSRS-01A-01D','TARGET-30-PAPCTS-01A-01D') & all.subclones$chr==1 & all.subclones$endpos<150000000,c('sample','chr','startpos','endpos')], 'darkgreen')
colnames(subclones.sub)[5] =  'color'
colnames(manual.segs) = colnames(subclones.sub)
subclones.sub = rbind(subclones.sub, manual.segs)
subclones.sub$color = as.character(subclones.sub$color)
subclones.sub$startpos = as.integer(subclones.sub$startpos)
subclones.sub$endpos = as.integer(subclones.sub$endpos)
subclones.sub = subclones.sub[order(subclones.sub$endpos, decreasing=T),]
vals = (1:length(unique(subclones.sub$sample))) / 2
names(vals) = unique(subclones.sub$sample)
subclones.sub$val = vals[subclones.sub$sample]
plot(c(), ylim=range(subclones.sub$val), xlim=c(1,chrs[chrs$V1==1,2]), bty="n", xlab="", ylab="", xaxt="n")
segments(x0=subclones.sub$startpos, x1=subclones.sub$endpos, y0=as.numeric(subclones.sub$val), lwd=3, col=subclones.sub$color)
text.data = unique(subclones.sub[,c('sample','val')])
text(x=rep(178000000, nrow(text.data)), y=text.data$val, labels=text.data$sample, cex=0.6)
##############
chr1p.loss = rep(0, nrow(annotations))
chr1p.loss[which(annotations$sample %in% subclones.sub$sample)] = 1
#### chr3p loss events
sel.samples = annotations[grepl('^LOH_p', annotations$chr3p_events),'sample']
chr = 3
subclones.sub = cbind(all.subclones[all.subclones$sample %in% sel.samples & all.subclones$chr==chr & all.subclones$endpos<=min(bands[bands$V1==paste('chr',chr, sep='') & bands$V5=='acen','V2']) & all.subclones$nMin1_A==0,c('sample','chr','startpos','endpos')],'darkgreen')
colnames(subclones.sub)[5] = 'color'
subclones.sub = subclones.sub[!(subclones.sub$startpos %in% c(61495)),]
manual.segs = rbind(
	c('TARGET-30-PASKJX-01A-01D', 3, 61495, 36794480, 'darkgreen'), c('TARGET-30-PASKJX-01A-01D', 3, 36794480, 37189247, 'hotpink'), c('TARGET-30-PASKJX-01A-01D', 3, 36794480, 59638534, 'darkgreen')
)
colnames(manual.segs) = colnames(subclones.sub)
subclones.sub = rbind(subclones.sub, manual.segs)
subclones.sub$color = as.character(subclones.sub$color)
subclones.sub$startpos = as.integer(subclones.sub$startpos)
subclones.sub$endpos = as.integer(subclones.sub$endpos)
subclones.sub = subclones.sub[order(subclones.sub$endpos, decreasing=T),]
vals = (1:length(unique(subclones.sub$sample))) / 2
names(vals) = unique(subclones.sub$sample)
subclones.sub$val = vals[subclones.sub$sample]
plot(c(), ylim=range(subclones.sub$val), xlim=c(1,chrs[chrs$V1==chr,2]), bty="n", xlab="", ylab="", xaxt="n")
segments(x0=subclones.sub$startpos, x1=subclones.sub$endpos, y0=as.numeric(subclones.sub$val), lwd=3, col=subclones.sub$color)
text.data = unique(subclones.sub[,c('sample','val')])
text(x=rep(143000000, nrow(text.data)), y=text.data$val, labels=text.data$sample, cex=0.6)
##############
chr3p.loss = rep(0, nrow(annotations))
chr3p.loss[which(annotations$sample %in% subclones.sub$sample)] = 1
#### chr4p loss events
sel.samples = annotations[grepl('^LOH_p', annotations$chr4p_event),'sample']
chr = 4
parm_end = min(bands[bands$V1==paste('chr',chr, sep='') & bands$V5=='acen','V2']) + 1e7
subclones.sub = cbind(all.subclones[all.subclones$sample %in% sel.samples & all.subclones$chr==chr & all.subclones$endpos<=parm_end & all.subclones$nMin1_A==0,c('sample','chr','startpos','endpos')],'darkgreen')
colnames(subclones.sub)[5] = 'color'
subclones.sub$color = as.character(subclones.sub$color)
subclones.sub$startpos = as.integer(subclones.sub$startpos)
subclones.sub$endpos = as.integer(subclones.sub$endpos)
subclones.sub = subclones.sub[order(subclones.sub$endpos, decreasing=T),]
vals = (1:length(unique(subclones.sub$sample))) / 2
names(vals) = unique(subclones.sub$sample)
subclones.sub$val = vals[subclones.sub$sample]
plot(c(), ylim=range(subclones.sub$val), xlim=c(1,chrs[chrs$V1==chr,2]), bty="n", xlab="", ylab="", xaxt="n")
segments(x0=subclones.sub$startpos, x1=subclones.sub$endpos, y0=as.numeric(subclones.sub$val), lwd=3, col=subclones.sub$color)
text.data = unique(subclones.sub[,c('sample','val')])
text(x=rep(parm_end+3e7, nrow(text.data)), y=text.data$val, labels=text.data$sample, cex=0.6)
##############
chr4p.loss = rep(0, nrow(annotations))
chr4p.loss[which(annotations$sample %in% subclones.sub$sample)] = 1
#### chr9p loss events
sel.samples = annotations[grepl('LOH_p',annotations$chr9p_events) | grepl('LOH_SMALL',annotations$chr9p_events),'sample']
chr = 9
subclones.sub = cbind(all.subclones[all.subclones$sample %in% sel.samples & all.subclones$chr==chr & all.subclones$endpos<=min(bands[bands$V1==paste('chr',chr, sep='') & bands$V5=='acen','V2']) & all.subclones$nMin1_A==0,c('sample','chr','startpos','endpos')],'darkgreen')
colnames(subclones.sub)[5] =  'color'
subclones.sub = subclones.sub[!(subclones.sub$startpos %in% c(203811, 273233, 203761, 206512)),]
manual.segs = rbind(
	c('TARGET-30-PASUYG-01A-01D', 9, 203811, 8354222, 'darkgreen'), c('TARGET-30-PASUYG-01A-01D', 9, 8354222, 8652977, 'hotpink'),
	c('TARGET-30-PASUYG-01A-01D', 9, 8652977, 11339418, 'darkgreen'), c('TARGET-30-PASUYG-01A-01D', 9, 11339418, 14249755, 'hotpink'),
	c('TARGET-30-PASUYG-01A-01D', 9, 14249755, 27183340, 'darkgreen'), c('TARGET-30-PARAMT-01A-01D', 9, 8870273, 9231109, 'hotpink'),
	c('TARGET-30-PASXHE-01A-01D', 9, 19713156, 26606272, 'hotpink'), c('TARGET-30-PASWVY-01A-01D', 9, 273233, min(bands[bands$V1==paste('chr',chr, sep='') & bands$V5=='acen','V2']), 'darkgreen'),
	c('TARGET-30-PARKGJ-01A-01D', 9, 20151950, 21959236, 'darkgreen'), c('TARGET-30-PARKGJ-01A-01D', 9, 21959236, 22086093, 'hotpink'),
	c('TARGET-30-PARKGJ-01A-01D', 9, 22086093, 38771831, 'darkgreen'), c('TARGET-30-PAPSEI-01A-01D', 9, 8691072, 9143521, 'hotpink'),
	c('TARGET-30-PAMNLH-01A-01D', 9, 20703584, 22514268, 'darkgreen'), c('TARGET-30-PAMNLH-01A-01D', 9, 8562775, 9059566, 'darkgreen'),
	c('TARGET-30-PATNKP-02A-01D', 9, 8696292, 9326249, 'hotpink'), c('TARGET-30-PATCFL-01A-01D', 9, 203761, min(bands[bands$V1==paste('chr',chr, sep='') & bands$V5=='acen','V2']), 'darkgreen'),
	c('TARGET-30-PASBEN-01A-01D', 9, 9072655, 9172128, 'darkgreen'), c('TARGET-30-PATBMM-01A-01D', 9, 203937, 19149599, 'darkgreen'),
	c('TARGET-30-PARZIP-01A-01D', 9, 9514929, 9566864, 'darkgreen'), c('TARGET-30-PARDCK-01A-01D', 9, 203761, 33838690, 'darkgreen'),
	c('TARGET-30-PARBLH-01A-01D', 9, 8470168, 8533182, 'darkgreen'), c('TARGET-30-PARBGP-01A-01D', 9, 8378504, 8509808, 'darkgreen'),
	c('TARGET-30-PAPTIP-01A-01D', 9, 203761, 29383088, 'darkgreen')
)
colnames(manual.segs) = colnames(subclones.sub)
subclones.sub = rbind(subclones.sub, manual.segs)
subclones.sub$color = as.character(subclones.sub$color)
subclones.sub$startpos = as.integer(subclones.sub$startpos)
subclones.sub$endpos = as.integer(subclones.sub$endpos)
subclones.sub[subclones.sub$startpos==21834543,'color'] = 'hotpink'
subclones.sub = subclones.sub[order(subclones.sub$endpos, decreasing=T),]
vals = (1:length(unique(subclones.sub$sample))) / 2
names(vals) = unique(subclones.sub$sample)
subclones.sub$val = vals[subclones.sub$sample]
plot(c(), ylim=range(subclones.sub$val), xlim=c(1,chrs[chrs$V1==chr,2]), bty="n", xlab="", ylab="", xaxt="n")
segments(x0=subclones.sub$startpos, x1=subclones.sub$endpos, y0=as.numeric(subclones.sub$val), lwd=4, col=subclones.sub$color)
abline(v=9e6, col='blue', lwd=1.5, lty=3)
abline(v=2.2e7, col='blue', lwd=1.5, lty=3)
text.data = unique(subclones.sub[,c('sample','val')])
text(x=rep(65000000, nrow(text.data)), y=text.data$val, labels=text.data$sample, cex=0.8)
##############
chr9p.loss = rep(0, nrow(annotations))
chr9p.loss[which(annotations$sample %in% subclones.sub$sample)] = 1
#### chr11p loss events
chr11q.loss.events = c(unique(grep('^LOH_q', annotations$chr11q_events, value=T)), '11qtel', 'SMALL_CNLOH', 'CCND1_AMP+SMALL_LOH', 'CNLOH_small')
chr11q.loss = rep(0, nrow(annotations))
chr11q.loss[annotations$chr11q_events %in% chr11q.loss.events] = 1

segs.loss = cbind(all.subclones[all.subclones$sample %in% annotations[chr11q.loss==1,'sample'] & all.subclones$chr==11 & all.subclones$startpos>=min(bands[bands$V1=='chr11' & bands$V5=='acen','V2']) & all.subclones$nMin1_A==0,c('sample','chr','startpos','endpos')], 'darkgreen')
segs.loss = segs.loss[!(segs.loss$sample %in% c('TARGET-30-PARSEA-01A-01D')),]
segs.loss = segs.loss[!(segs.loss$startpos %in% c(69010651)),]
manual.segs = rbind(
	c('TARGET-30-PARSEA-01A-01D', 11, 65998757, 70764889, 'orangered3'), c('TARGET-30-PARSEA-01A-01D', 11, 70764889, 134946184, 'darkgreen'),
	c('TARGET-30-PATCFL-01A-01D',11, 71633530, 134946184, 'darkgreen'), c('TARGET-30-PAPVRN-01A-01D',11,78820297, 134946184, 'darkgreen'),
	c('TARGET-30-PASKJX-01A-01D', 11, 69838610, 134944160, 'darkgreen')
)
colnames(segs.loss)[5] = 'color'
colnames(manual.segs) = colnames(segs.loss)
segs.loss = rbind(segs.loss, manual.segs)
sel1 = all.subclones$sample %in% annotations[chr11q.loss==1,'sample'] & all.subclones$chr==11 & all.subclones$startpos>=3e7 & all.subclones$nMaj1_A>1
sel2 = all.subclones$sample %in% annotations[chr11q.loss==1,'sample'] & all.subclones$chr==11 & all.subclones$startpos>=3e7 & all.subclones$frac1_A<1 & all.subclones$nMaj2_A>1
segs.gain = cbind(all.subclones[sel1 | sel2,c('sample','chr','startpos','endpos')], 'orangered3')
colnames(segs.gain)[5] = 'color'
segs.gain = segs.gain[!(segs.gain$startpos %in% c(37889834, 48347267, 49483346, 48173404, 48182109, 48346729, 47778125, 55585438, 48332139, 48441008, 48346999, 68335797)),]
manual.segs = rbind(
	c('TARGET-30-PAPUWY-01A-01D', 11, 71236934, 71275195, 'orangered3'), c('TARGET-30-PAPWUC-01A-01D', 11, 66244008, 69466771, 'orangered3'),
	c('TARGET-30-PARBGP-01A-01D', 11, 68878492, 69030545, 'orangered3'), c('TARGET-30-PASFGG-01A-01D', 11, 65663132, 72726168, 'orangered3'),
	c('TARGET-30-PATBMM-01A-01D', 11, 68859812, 69194695, 'orangered3'), c('TARGET-30-PASKJX-01A-01D', 11, 68335797, 69010886, 'gold'),
	c('TARGET-30-PASKJX-01A-01D', 11, 69010886, 69838610, 'orangered3')
)
colnames(manual.segs) = colnames(segs.gain)
segs.gain = rbind(segs.gain, manual.segs)
segs.gain = segs.gain[!(segs.gain$startpos %in% segs.loss$startpos),]
subclones.sub = rbind(segs.loss, segs.gain)

subclones.sub$patient = grep('^P', unlist(strsplit(subclones.sub$sample, split='-')), value=T)
subclones.sub$startpos = as.integer(subclones.sub$startpos)
subclones.sub$endpos = as.integer(subclones.sub$endpos)
subclones.sub$width = subclones.sub$endpos - subclones.sub$startpos
subclones.sub$color = as.character(subclones.sub$color)
subclones.sub[subclones.sub$sample=='TARGET-30-PATBMM-01A-01D' & subclones.sub$startpos==68859812,'color'] = 'gold'
subclones.sub[subclones.sub$sample=='TARGET-30-PARFWB-01A-01D' & subclones.sub$startpos==67419144,'color'] = 'gold'
subclones.sub[subclones.sub$sample=='TARGET-30-PARDCK-01A-01D' & subclones.sub$startpos==65263001,'color'] = 'gold'
#subclones.sub = subclones.sub[order(subclones.sub$sample, decreasing=T),]
subclones.sub = subclones.sub[order(subclones.sub$startpos, decreasing=F),]
vals = (1:length(unique(subclones.sub$sample))) / 2
names(vals) = unique(subclones.sub$sample)
subclones.sub$val = vals[subclones.sub$sample]
plot(c(), ylim=range(subclones.sub$val), xlim=c(1,chrs[chrs$V1==1,2]), bty="n", xlab="", ylab="", xaxt="n")
segments(x0=subclones.sub$startpos, x1=subclones.sub$endpos, y0=as.numeric(subclones.sub$val), lwd=3, col=subclones.sub$color)
abline(v=68604629, col='blue', lwd=1.5, lty=3)
abline(v=71809473, col='blue', lwd=1.5, lty=3)
text.data = unique(subclones.sub[,c('sample','val')])
text(x=rep(145000000, nrow(text.data)), y=text.data$val, labels=text.data$sample, cex=0.3)
##############
mycn.amp = rep(0, nrow(annotations))
mycn.amp[annotations$MYCN_AMP=='HIGH'] = 1
cn.changes = cbind(chr1p.loss, chr3p.loss, chr9p.loss, chr11q.loss, mycn.amp)
rownames(cn.changes) = annotations[,'sample']
