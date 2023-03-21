all.pileup.files = list.files(path='results/', full.names=F)
genotypes = as.data.frame(mat.or.vec(97, length(all.pileup.files)))
colnames(genotypes) = sub('_calls.txt','',sub('1202_','',all.pileup.files))

s = read.table(file='sequenom_loci.txt', header=T, stringsAsFactors=F)
all.pileup.files = list.files(path='results/', full.names=T)
for(sample in colnames(genotypes)) {
	d = read.table(file=grep(paste(sample,'_calls.txt',sep=''), all.pileup.files, value=T), header=T, sep="\t", stringsAsFactors=F)
	m = merge(x=s, y=d, by.x='START', by.y='POS')
	m$ALT.n = rep(-1,nrow(m)); m$REF.n = rep(-1,nrow(m));
	for(i in 1:nrow(m)) {m[i,'REF.n']=m[i,m[i,'REF']]; m[i,'ALT.n']=m[i,m[i,'ALT']]}
	m$ALT.prop = m$ALT.n / m$DEPTH; m$REF.prop = m$REF.n / m$DEPTH;
	m$genotype = rep('none',nrow(m))
	m[m$DEPTH>20 & m$ALT.prop>0.90,'genotype'] = 'alt'
	m[m$DEPTH>20 & m$ALT.prop>=0.1 & m$ALT.prop<=0.90,'genotype'] = 'het'
	m[m$DEPTH>20 & !(m$genotype %in% c('alt','het')),'genotype'] = 'ref'
	genotypes[,sample] = m$genotype
}
genotypes.comparison = mat.or.vec(ncol(genotypes), ncol(genotypes))
colnames(genotypes.comparison) = colnames(genotypes)
rownames(genotypes.comparison) = colnames(genotypes)
for(s1 in colnames(genotypes.comparison)) {
	for(s2 in colnames(genotypes.comparison)) {
		informative = genotypes[,s1]!="none" & genotypes[,s2]!="none"
		genotypes.comparison[s1,s2] = length(which(genotypes[informative,s1]==genotypes[informative,s2])) / length(which(informative))
	}
}
write.table(genotypes.comparison, file='genotypes.comparison.txt', quote=F, sep="\t")
