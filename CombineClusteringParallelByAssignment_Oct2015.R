args=commandArgs(TRUE)
run = as.integer(args[1])

noiters=1000
conc_param=1
cluster_conc = 5
density.smooth = 0.01
burn.in = 200
basedir = "/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/"
setwd(paste(basedir,"prostate",sep=""))

samplenames=c("PD12337","PD13412","PD11328","PD11329","PD11330","PD11331","PD11332","PD11333","PD11334","PD11335")
subsamples = list()
subsamples[[1]] = c("a","c","d","e","f","g","h","i","j")
subsamples[[2]] = c("a","c","d","e","f")
subsamples[[3]] = c("a","c","d","e")
subsamples[[4]] = c("a","c","d","e","f","g","h","i","j","k")
subsamples[[5]] = c("a","c")
subsamples[[6]] = c("a","c","d","e","f")
subsamples[[7]] = c("a","c","d","e","f")
subsamples[[8]] = c("a","c","d")
subsamples[[9]] = c("a","c","d","e")
subsamples[[10]] = c("a","c","d","e")

cellularities = list()
cellularities[[1]]=c(0.8879469,0.8793653,0.9081555,0.868476,0.8707937,0.9065765,0.878091,0.7842113,0.8164292)
cellularities[[2]]=c(0.8976817,0.8788655,0.8725436,0.8206648,0.8784499)
cellularities[[3]] = c(0.930,0.856,0.912,0.826)
#cellularities[[4]] = c(0.800,0.675,0.772,0.890,0.808,0.704,0.763,0.672,0.747,0.817)
#25June2014 - new Bberg for PD11329j
cellularities[[4]] = c(0.800,0.675,0.772,0.890,0.808,0.704,0.763,0.672,0.798,0.817)
cellularities[[5]]=c(0.327,0.434)
#cellularities[[6]]=c(0.896,0.496,0.810,0.860,0.778)
#25June2014 - new Bberg for PD11331c
cellularities[[6]]=c(0.896,0.625,0.810,0.860,0.778)
#210714 - revert to old ploidy for PD11331c
#cellularities[[6]]=c(0.896,0.496,0.810,0.860,0.778)
cellularities[[7]]=c(0.824,0.907,0.895,0.883,0.918)
cellularities[[8]]=c(0.772,0.619,0.810)
cellularities[[9]] = c(0.753,0.802,0.849,0.863)
cellularities[[10]]=c(0.794,0.840,0.851,0.865)

samplename = samplenames[run]
cellularity = cellularities[[run]]

no.subsamples = length(subsamples[[run]])
choose.number = 3
choose.from = no.subsamples
no.perms = choose(no.subsamples,3)

new_output_folder = "/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/prostate/clustering_parallel_October2015"

get.subsample.indices<-function(choose.index,choose.number,choose.from){
	subsample.indices = 1:choose.number
	temp.choose.index = choose.index
	for(i in 1:choose.number){
		last.subtotal = 0
		subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
		while(temp.choose.index>subtotal){
			subsample.indices[i] = subsample.indices[i] + 1
			last.subtotal = subtotal
			subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)			
		}
		if(i<choose.number){
			subsample.indices[i+1] = subsample.indices[i]+1
		}
		temp.choose.index = temp.choose.index - last.subtotal
	}
	return(subsample.indices)
}


data=list()
for(s in 1:no.subsamples){
	data[[s]]=read.table(paste(samplename,subsamples[[run]][s],"_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_Oct2015.txt",sep=""),sep="\t",header=T)
}

library(MASS)
library(MCMCpack)
library(mvtnorm)
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/multidimensionalDensityEstimator.R")
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/interconvertMutationBurdens.R")

WTCount = data[[1]]$WT.count
mutCount = data[[1]]$mut.count
totalCopyNumber = data[[1]]$subclonal.CN
copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
normalCopyNumber = data[[1]]$normalCopyNumber
for(s in 2:length(subsamples[[run]])){
	WTCount = cbind(WTCount,data[[s]]$WT.count)
	mutCount = cbind(mutCount,data[[s]]$mut.count)
	totalCopyNumber = cbind(totalCopyNumber,data[[s]]$subclonal.CN)
	copyNumberAdjustment = cbind(copyNumberAdjustment,data[[s]]$no.chrs.bearing.mut)
	normalCopyNumber = cbind(normalCopyNumber,data[[1]]$normalCopyNumber)
}

non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
normalCopyNumber = normalCopyNumber[non.zero,]
for(s in 1:length(subsamples[[run]])){
	data[[s]]=data[[s]][non.zero,]
}
mutation.copy.number = array(NA,dim(totalCopyNumber))
for(i in 1:length(subsamples[[run]])){
	mutation.copy.number[,i] = mutationBurdenToMutationCopyNumber(mutCount[,i] / (mutCount[,i]+WTCount[,i]) , totalCopyNumber[,i], cellularity[i], normalCopyNumber[,i])
}

no.muts = nrow(data[[1]])

node.assignments=NULL
for(choose.index in 1:no.perms){
	print(choose.index)
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	print(subsample.indices)
	perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and cluster_info_0.01.txt",sep=""),header=T,sep="\t")

	if(is.null(node.assignments)){
		node.assignments = array(data.matrix(perm.assignments[ncol(perm.assignments)]),c(no.muts,1))
	}else{
		node.assignments = cbind(node.assignments,data.matrix(perm.assignments[ncol(perm.assignments)]))
	}
}

#get different combinations of assignments from the different permutations
unique.assignments = unique(node.assignments)
write.table(unique.assignments,paste(new_output_folder,"/",samplename,"_matchedClustersInParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

no.consensus.nodes = nrow(unique.assignments)
print(paste("#consensus nodes = ",no.consensus.nodes,sep=""))
consensus.assignments = vector(mode="numeric",length = no.consensus.nodes)	
for(n in 1:no.consensus.nodes){
	print(n)
	consensus.assignments[sapply(1:no.muts,function(a,u,i){all(u==a[i,])},a=node.assignments,u=unique.assignments[n,])]=n
}

#get probabilities of assignment to each set of assignments
prob.consensus.assignments = array(1,c(no.muts,no.consensus.nodes)) 
#prob.consensus.assignments = array(0,c(no.muts,no.consensus.nodes)) #get average probability, rather than product
for(choose.index in 1:no.perms){
	print(choose.index)
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and cluster_info_0.01.txt",sep=""),header=T,sep="\t")
	for(n in 1:no.consensus.nodes){
		prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]+2]
		#prob.consensus.assignments[,n] = prob.consensus.assignments[,n] + perm.assignments[,unique.assignments[n,choose.index]+2]
	}
}
#prob.consensus.assignments = prob.consensus.assignments / no.perms

out = data[[1]][,1:2]
for(s in 1:no.subsamples){
	out = cbind(out,data[[s]]$subclonal.fraction)
}
out = cbind(out,prob.consensus.assignments)
names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),paste("prob.cluster",1:no.consensus.nodes,sep=""))
write.table(out,paste(new_output_folder,"/",samplename,"_allClusterProbabilitiesFromParallelRuns_16Oct2014.txt",sep=""),sep="\t",row.names=F,quote=F)

##############################################################################################################################################
all.CIs = array(NA,c(no.consensus.nodes,no.subsamples,no.perms,2))
for(choose.index in 1:no.perms){
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	print(subsample.indices)
	CIs = data.matrix(read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_confInts_",density.smooth,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F))
	for(n in 1:no.consensus.nodes){
		all.CIs[n,subsample.indices,choose.index,1] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)-1]
		all.CIs[n,subsample.indices,choose.index,2] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)]
	}
}
median.CIs = array(NA,c(no.consensus.nodes,no.subsamples,2))
for(n in 1:no.consensus.nodes){
	for(s in 1:no.subsamples){
		for(c in 1:2){
			median.CIs[n,s,c] = median(all.CIs[n,s,,c],na.rm=T)
		}
	}
}
median.CIs.2D = t(sapply(1:no.consensus.nodes,function(m,i){as.vector(t(m[i,,]))},m=median.CIs))
write.table(cbind(1:no.consensus.nodes,median.CIs.2D,table(consensus.assignments)),sep="\t",row.names=F,quote=F,col.names=c("cluster.no",paste(rep(paste(samplename,subsamples[[run]],sep=""),each=2),rep(c("lowerCI","upperCI"),times=no.subsamples),sep="_"),"no.muts"),paste(new_output_folder,"/",samplename,"_consensusClustersByParallelNodeAssignment_16Oct2014.txt",sep=""))
out = data[[1]][,1:2]
for(s in 1:no.subsamples){
	out = cbind(out,data[[s]]$subclonal.fraction)
}
out = cbind(out,consensus.assignments)
names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),"cluster.no")
write.table(out,paste(new_output_folder,"/",samplename,"_allClusterassignmentsFromParallelRuns_16Oct2014.txt",sep=""),sep="\t",row.names=F,quote=F)

pdf(paste(new_output_folder,"/","consensus_cluster_assignment_",samplename,"_combined_",density.smooth,"_16Oct2014.pdf",sep=""),height=4,width=4)
#its hard to distinguish more than 8 different colours
max.cols = 8
cols = rainbow(min(max.cols,no.consensus.nodes))
plot.data = mutation.copy.number/copyNumberAdjustment
plot.data[is.na(plot.data)]=0
for(i in 1:(no.subsamples-1)){
	for(j in (i+1):no.subsamples){
		plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,(subsamples[[run]][subsample.indices])[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[[run]][j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
		for(n in 1:no.consensus.nodes){
			pch=20 + floor((n-1)/max.cols)
			#pch is not implmeneted above 25
			if(pch>25){
				pch=pch-20
			}
			points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=pch,cex=0.5)
		}
		pch=20 + floor((0:(no.consensus.nodes-1))/max.cols)
		pch[pch>25] = pch[pch>25]-20
		legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.consensus.nodes,col=cols[(0:(no.consensus.nodes-1)) %% max.cols + 1],pch=pch,cex=1)

	}
}	
dev.off()

q(save="no")
