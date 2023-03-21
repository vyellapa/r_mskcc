library(hdp)
library(deconstructSigs)
library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
direc <- args[1] #Current directory
inFile <- args[2] #Input file


setwd(direc)
type.cols = c(rgb(86/255,180/255,233/255), rgb(0/255,0/255,0/255), rgb(255/255,0/255,0/255), rgb(180/255,180/255,180/255), rgb(107/255,228/255,66/255), rgb(255/255,175/255,178/255))

# Normalized counts for 96 subtypes
nb.counts = read.table(file=inFile)
format.cols <- function(x) {
        return(paste(
                unlist(strsplit(unlist(strsplit(colnames(nb.counts)[x], split="\\."))[3], split=''))[1],
                "[",
                unlist(strsplit(colnames(nb.counts)[x], split="\\."))[1],
                ">",
                unlist(strsplit(colnames(nb.counts)[x], split="\\."))[2],
                "]",
                unlist(strsplit(unlist(strsplit(colnames(nb.counts)[x], split="\\."))[3], split=''))[3],
                sep='')
        )
}
colnames(nb.counts) = sapply(1:length(colnames(nb.counts)), format.cols)
# Plot the normalized counts for 96 mutation subtypes
mt.96.nuc.context = NULL
for(i in 1:nrow(nb.counts)) {
        s = whichSignatures(tumor.ref=nb.counts, signatures.ref = signatures.cosmic, sample.id=rownames(nb.counts)[i], contexts.needed=TRUE, tri.counts.method="genome")
        mt.96.nuc.context = rbind(s$tumor, mt.96.nuc.context)
}
mt.96.nuc.context.melt = melt(mt.96.nuc.context)
colnames(mt.96.nuc.context.melt) = c('Sample','Subtype','Normalized_prop')
mt.96.nuc.context.melt$Type = rep('',nrow(mt.96.nuc.context.melt))
mt.96.nuc.context.melt[grep('C>A',mt.96.nuc.context.melt$Subtype),'Type'] = 'C>A'
mt.96.nuc.context.melt[grep('C>G',mt.96.nuc.context.melt$Subtype),'Type'] = 'C>G'
mt.96.nuc.context.melt[grep('C>T',mt.96.nuc.context.melt$Subtype),'Type'] = 'C>T'
mt.96.nuc.context.melt[grep('T>A',mt.96.nuc.context.melt$Subtype),'Type'] = 'T>A'
mt.96.nuc.context.melt[grep('T>C',mt.96.nuc.context.melt$Subtype),'Type'] = 'T>C'
mt.96.nuc.context.melt[grep('T>G',mt.96.nuc.context.melt$Subtype),'Type'] = 'T>G'
p = ggplot(mt.96.nuc.context.melt) + geom_bar(aes(x=factor(Subtype), y=Normalized_prop, fill=Type), stat="identity")
p = p +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()+facet_grid(Sample ~ .) + theme(legend.position="none")
p = p + scale_fill_manual(values=type.cols)
p= p + theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background  = element_blank(), panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80")) + theme(strip.text.y = element_text(size = 9)) +xlab("")
pdf('hdp_normalised_96_mt_type_raw.pdf', height=20, width=20)
print(p)
dev.off()
pdf('hdp_normalised_96_mt_type_zoom.pdf', height=20, width=20)
p = p + scale_y_continuous(limits = c(0, summary(mt.96.nuc.context.melt$Normalized_prop)['3rd Qu.'] + 0.1))
print(p)
dev.off()
nb.count = nb.counts[,colnames(mt.96.nuc.context)]

#############################################
# Run HDP

nb.counts = read.table(file=inFile, header=T)
nb.counts = nb.counts[,sort(colnames(nb.counts))]
colnames(nb.counts) = colnames(lusc_count)

# Read the COSMIC signatures
cosmic.sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', header=TRUE, sep='\t')
cosmic.sigs = cosmic.sigs[,1:33]
colnames(cosmic.sigs)[1:3] = c('Substitution.Type','Trinucleotide','Somatic.Mutation.Type')
cosmic.sigs <- cosmic.sigs[order(cosmic.sigs$Substitution.Type, cosmic.sigs$Trinucleotide),]
sigs <- as.matrix(cosmic.sigs[,grep('Signature', colnames(cosmic.sigs))])

##################
n_lung = nrow(lusc_count)
n_nbl = nrow(nb.counts)
sign_number=30
# The following is the custom HDP structure for this dataset
# Top DP node = 0
# 32 children DP nodes for 
	# A single node for each signature 1...31
	# Parent node for all the lung samples 32
	# Parent node for all the neuroblastoma smaples 33
		# Nodes for lung samples 34...133
		# Nodes for neuroblastoma sample 134 142

# Initialise the HDP structure conditioning on the COSMIC signatures
prior_pseudoc <- rep(100, ncol(sigs)) # The value 100 is arbitrary.
full_prior <- hdp_prior_init(sigs, prior_pseudoc, hh=rep(1, nrow(sigs)), alphaa=c(1, 1), alphab=c(1, 2))
# Add concentration parameters: one for the top node, one for the signatures plus the parents for the lungs and the nbls, one for lung and one for nbls
full_prior <- hdp_addconparam(full_prior, c(1,1,1,1), c(1,1,1,1))
full_prior <- hdp_adddp(full_prior, (1 + n_lung + 1 + n_nbl), # nodes for the lung parent, the lungs, the nbl parent and the nbl samples
                  ppindex=c(1, 1, rep(1+sign_number+1, n_lung), rep(1+sign_number+1+1, n_nbl)), # index=n_sigs+2 for the lungs and n_sigs+3 for the the nbls
                  cpindex=c(3, 3, rep(4, n_lung + n_nbl)))
# Attach data to the nodes
full_prior <- hdp_setdata(full_prior, (1+sign_number+1+1)+1:n_lung, lusc_count)
full_prior <- hdp_setdata(full_prior, (1+sign_number+1+1+n_lung)+1:n_nbl, nb.counts)
# Run a set of HDP chains
chlist <- vector("list", 100)
for (i in 1:length(chlist)){
        nb_pr <- dp_activate(full_prior, (1+sign_number+1)+0:(n_lung+n_nbl+1), initcc=sign_number)
        chlist[[i]] <- hdp_posterior(nb_pr, burnin=1500, n=50, space=50, cpiter=3, seed=i*1e6)
}
saveRDS(nb_pr, 'hdp_chains.rds')
# Extract components from the HDP chain
full_pr_multi <- hdp_multi_chain(chlist)
full_pr_ec <- hdp_extract_components(full_pr_multi)
saveRDS(full_pr_ec, 'hdp_components.rds')

pdf('hdp_all_samples.pdf')
plot_dp_comp_exposure(full_pr_ec, 1+sign_number+1+1+(1:(n_lung+n_nbl)), incl_nonsig = TRUE, col=rainbow(31))
dev.off() 
# Get the exposure matrix for the neuroblastoma samples only
nb.indeces = 1+sign_number+1+1+((n_lung+1):(n_lung+n_nbl))
d.exposure = comp_dp_distn(full_pr_ec)$mean[nb.indeces,]
rownames(d.exposure) = rownames(nb.counts)
d.exposure.melt = melt(d.exposure)
d.exposure.melt = d.exposure.melt[d.exposure.melt$value>0.05,]
# Plot all signatures per sample
p = ggplot(d.exposure.melt) + geom_bar(aes(x=factor(Var1), y=value, fill=Var2), stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC", "#66CC99", "red", "blue", "green"))
pdf('hdp_ALL_facet_by_samples.pdf')
print(p)
dev.off()
# Plot all samples per signature
p = ggplot(d.exposure.melt) + geom_bar(aes(x=factor(Var1), y=value, fill=Var2), stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC", "#66CC99", "red", "blue", "green"))
p = p + facet_grid(Var2 ~.)
pdf('hdp_ALL_facet_by_signatures.pdf', height=10)
print(p)
dev.off()
# Get the exposure matrix for the neuroblastoma samples only
lung.indeces = 1+sign_number+1+1+(1:n_lung)
d.exposure = comp_dp_distn(full_pr_ec)$mean[lung.indeces,]
rownames(d.exposure) = rownames(lusc_count)
d.exposure.melt = melt(d.exposure)
d.exposure.melt = d.exposure.melt[d.exposure.melt$value>0.05,]
p = ggplot(d.exposure.melt) + geom_bar(aes(x=factor(Var1), y=value, fill=Var2), stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC", "#66CC99", "red", "blue", "green"))
pdf('hdp_lusc_facet_by_samples.pdf', width=12)
print(p)
dev.off()

