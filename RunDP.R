#library("dpclust3p")
library("DPClust")

#patch to cope with zero allele frequency in some samples
library(R.utils)
source("/well/wedge/shared/software/subclone_Dirichlet_Gibbs_sampler_nD_binomial_DCW.R")
source("/well/wedge/shared/software/plotnD.DCW.R")
reassignInPackage("plotnD",pkgName="DPClust",plotnD.DCW)
reassignInPackage("subclone.dirichlet.gibbs",pkgName="DPClust",subclone.dirichlet.gibbs.DCW)
reassignInPackage("Gibbs.subclone.density.est",pkgName="DPClust",Gibbs.subclone.density.est.DCW)
reassignInPackage("Gibbs.subclone.density.est.nd",pkgName="DPClust",Gibbs.subclone.density.est.nd.DCW)

args = commandArgs(T)
run=as.numeric(args[1])

setwd("/well/wedge/david/data/organoids")

input.files = list.files(pattern=".snv.readcounts.filtered.txt.gz")
samplenames = gsub(".snv.readcounts.filtered.txt.gz","",input.files)

info = read.table("tumour_normal_pairs.txt",sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
info[,1] = gsub("/","_",info[,1])
info[,2] = gsub("-",".",info[,2])
info[,3] = gsub("-",".",info[,3])


#THIS DOESN'T WORK FOR MULTIPLE SAMPLES
#interface:
#(loci_file, allele_frequencies_file, cellularity_file,
#    subclone_file, gender, SNP.phase.file, mut.phase.file, output_file)
#runGetDirichletProcessInfo()

donor = samplenames[run]
print(donor)

subsamplenames = info[info[,1]==donor,2]
no.samples = length(subsamplenames)
subsamplenames = paste(donor,subsamplenames,sep="_")

#get cellularities
cellularity = NULL
for(i in 1:no.samples){
	rho.psi = read.table(paste("Battenberg_",subsamplenames[i],"/tmpBattenberg/",subsamplenames[i],"_T_rho_and_psi.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)			
	cellularity = c(cellularity,rho.psi["FRAC_GENOME","rho"])
}

#get DP info
DP.files = paste("DP_info/",subsamplenames,"_DPinfo_withoutChrLosses.txt",sep="")

mutCount = NULL
WTCount = NULL
totalCopyNumber = NULL
copyNumberAdjustment = NULL
mutation.copy.number = NULL
for(i in 1:length(subsamplenames)){
	DP.info = read.table(DP.files[i],sep="\t",header=T, stringsAsFactors=F)
	if(i==1){
		mutCount = DP.info$mut.count
		WTCount = DP.info$WT.count
		totalCopyNumber = DP.info$subclonal.CN
		copyNumberAdjustment = DP.info$no.chrs.bearing.mut
		mutation.copy.number = DP.info$mutation.copy.number
	}else{
		mutCount = cbind(mutCount,DP.info$mut.count)
		WTCount = cbind(WTCount,DP.info$WT.count)
		totalCopyNumber = cbind(totalCopyNumber,DP.info$subclonal.CN)
		copyNumberAdjustment = cbind(copyNumberAdjustment,DP.info$no.chrs.bearing.mut)
		mutation.copy.number = cbind(mutation.copy.number,DP.info$mutation.copy.number)	
	}
}
#we need to filter the mutations
MCN.sum = rowSums(mutation.copy.number)
adj.sum = rowSums(copyNumberAdjustment)
exclude.indices = which(is.na(MCN.sum)|MCN.sum==0|is.na(adj.sum)|adj.sum==0)
mutCount = mutCount[-exclude.indices,]
WTCount = WTCount[-exclude.indices,]
totalCopyNumber = totalCopyNumber[-exclude.indices,]
copyNumberAdjustment = copyNumberAdjustment[-exclude.indices,]
mutation.copy.number = mutation.copy.number[-exclude.indices,]

no.iters=1000
if(nrow(mutation.copy.number)>50000){
	no.iters=250
}

output_folder = "DP_output"
if(!file.exists(output_folder)){
	dir.create(output_folder)
}

#This function does everything. It may be better to run separate function to perform Gibbs sampling and mutation assignment
DirichletProcessClustering(mutCount = mutCount, WTCount = WTCount, totalCopyNumber = totalCopyNumber, copyNumberAdjustment = copyNumberAdjustment,
mutation.copy.number = mutation.copy.number, cellularity = cellularity, output_folder = output_folder, no.iters = no.iters, no.iters.burn.in = ceiling(no.iters/5),
subsamplesrun = subsamplenames,samplename=donor, conc_param = 1, cluster_conc = 5, mut.assignment.type = 1, most.similar.mut = NA, mutationTypes = NA, min.frac.snvs.cluster = NA)
