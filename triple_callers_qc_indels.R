library(reshape2)

setwd("/ifs/res/papaemme/users/gg10/neuroblastoma_analysis_across_projects/triple_callers")

samples = read.table('all_samples.txt', header=T, stringsAsFactors=F, sep="\t")
samples$passed_pindel = 0
samples$passed_strelka = 0
samples$passed_mutect = 0
samples$passed_triple = 0
samples$passed_twoplus_lcc = 0
for(i in 1:nrow(samples)) {
	sample = samples[i,'sample']
	print(sample)
	triple_dir = paste(samples[i,'dir'], sample, "triple_callers", sep="/")
	indel_file = list.files(path=triple_dir, pattern="*.indels.output.annot.tsv.gz", full.names=T)

	cmd = paste("less ", indel_file, " | awk '$1!~/#/'", sep="")
	indels = system(cmd, intern = TRUE)
	names = unlist(strsplit(indels[1], split="\t"))
	indels = colsplit(string=indels[2:length(indels)], pattern="\t", names=names)

	passed_pindel = indels[which(grepl('P',indels$PASSED_BY)),]
	passed_strelka = indels[which(grepl('S',indels$PASSED_BY)),]
	passed_mutect = indels[which(grepl('M',indels$PASSED_BY)),]
	passed_triple = indels[which(indels$PASSED_BY=="M,P,S"),]
	passed_twoplus_lcc = indels[which(indels$ANY2LCC==1),]
	samples$passed_pindel[i] = nrow(passed_pindel)
	samples$passed_strelka[i] = nrow(passed_strelka)
	samples$passed_mutect[i] = nrow(passed_mutect)
	samples$passed_triple[i] = nrow(passed_triple)
	samples$passed_twoplus_lcc[i] = nrow(passed_twoplus_lcc)
}
pdf('mut_counts.pdf', width=25, height=10)
par(mar=c(8.1,4.1,4.1,2.1), mfrow=c(2,1))
barplot(samples$passed_twoplus_lcc, names=samples$sample, las=2)
barplot(samples$number_subs, names=samples$sample, las=2)
dev.off()

densm = density(passed_mutect$TARGET_VAF_MEAN)
denss = density(passed_strelka$TARGET_VAF_MEAN)
densp = density(passed_pindel$TARGET_VAF_MEAN)
denslcc = density(passed_twoplus_lcc$TARGET_VAF_MEAN)
maxy = mean(c(max(denslcc$y), max(densm$y), max(denss$y), max(densp$y)))

par(mfrow=c(1,2))
plot(density(passed_twoplus_lcc$TARGET_VAF_MEAN), ylim=c(0,maxy), lwd=2, main=sample, xlab="", ylab="")
abline(v=purity/2)
lines(density(passed_pindel$TARGET_VAF_MEAN), col="red", lwd=2)
lines(density(passed_strelka$TARGET_VAF_MEAN), col="darkgreen", lwd=2)
lines(density(passed_mutect$TARGET_VAF_MEAN), col="orange", lwd=2)
barplot(c(nrow(passed_twoplus_lcc),nrow(passed_pindel),nrow(passed_strelka),nrow(passed_mutect)), col=c("black","red","darkgreen","orange"), names=c("any2lcc","Pindel","Strelka","Mutect"))

