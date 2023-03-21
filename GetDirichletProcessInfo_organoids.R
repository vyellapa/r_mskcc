args = commandArgs(T)

source("/well/wedge/shared/software/GetDirichletProcessInfo_multipleSamples.R")

setwd("/gpfs2/well/wedge/david/data/organoids")

input.files = list.files(pattern=".snv.readcounts.filtered.txt.gz")
samplenames = gsub(".snv.readcounts.filtered.txt.gz","",input.files)

nucleotides = c("A","C","G","T")

#runs=1:length(all.donors)
#remaining donors have already run
runs=1:length(samplenames)
if(length(args)>0){
	runs=as.numeric(args[1])
}

info = read.table("tumour_normal_pairs.txt",sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
info[,1] = gsub("/","_",info[,1])
info[,2] = gsub("-",".",info[,2])
info[,3] = gsub("-",".",info[,3])

for(run in runs){
	donor = samplenames[run]
	print(donor)
	mut.data = read.table(gzfile(input.files[run]),header=T,sep="\t",stringsAsFactors=F)
	mut.data[,1] = gsub("chr","",mut.data[,1])

	subsamplenames = info[info[,1]==donor,2]
	no.samples = length(subsamplenames)
	
	ref.indices = match(paste0("Ref.count.",subsamplenames),names(mut.data))
	alt.indices = ref.indices+1

	subsamplenames = paste(donor,subsamplenames,sep="_")
	
	#get mutation allele counts and cellularity
	mut.count = mut.data[,alt.indices]
	wt.count = mut.data[,ref.indices]	
	cellularity = NULL
	for(i in 1:no.samples){
		rho.psi = read.table(paste("Battenberg_",subsamplenames[i],"/tmpBattenberg/",subsamplenames[i],"_T_rho_and_psi.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)			
		cellularity = c(cellularity,rho.psi["FRAC_GENOME","rho"])
	}

	no.muts = nrow(wt.count)

	#get subclones files
	subcl.files = paste("Battenberg_",subsamplenames,"/tmpBattenberg/",subsamplenames,"_T_subclones.txt",sep="")

	#get DP info
	GetDirichletProcessInfo_multipleSamples(subsamplenames, cellularity, mut.count,wt.count,subclone.files=subcl.files,out.files = cbind(paste("DP_info/",subsamplenames,"_DPinfo.txt",sep=""),paste("DP_info/",subsamplenames,"_DPinfo_withoutChrLosses.txt",sep="")),info = mut.data[,c(1,2)],keep.muts.not.explained.by.CN=T)	
}
#check warnings
print(warnings())

#q(save="no")
