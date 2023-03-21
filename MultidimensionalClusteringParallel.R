args=commandArgs(TRUE)
run = as.integer(args[1])

if(length(args)>1){
	#subsample.indices = as.numeric(args[-1])
	choose.from = as.numeric(args[2])
	choose.number = as.numeric(args[3])
	choose.index = as.numeric(args[4])
	subsample.indices = 1:choose.number
	for(i in 1:choose.number){
		last.subtotal = 0
		subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
		while(choose.index>subtotal){
			subsample.indices[i] = subsample.indices[i] + 1
			last.subtotal = subtotal
			subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)			
		}
		if(i<choose.number){
			subsample.indices[i+1] = subsample.indices[i]+1
		}
		choose.index = choose.index - last.subtotal
	}
}

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

basedir = "/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/"
setwd(paste(basedir,"prostate",sep=""))

no.subsamples = length(subsamples[[run]])
#no.subsamples = length(subsample.indices)

data=list()
for(s in 1:no.subsamples){
	data[[s]]=read.table(paste(samplename,subsamples[[run]][s],"_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_Oct2015.txt",sep=""),sep="\t",header=T)
	#data[[s]]=read.table(paste(samplename,subsamples[[run]][subsample.indices[s]],"_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_Dec2013.txt",sep=""),sep="\t",header=T)
}

library(MASS)
library(MCMCpack)
library(mvtnorm)
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/multidimensionalDensityEstimator.R")
source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/interconvertMutationBurdens.R")

output_folder = paste("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/prostate/metastases_Oct2015_outputs_2D_DirichletProcess_",samplename,"_adjustedToSubclonalProportion",sep="")

new_output_folder = "/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/prostate/clustering_parallel_October2015"
if(!file.exists(new_output_folder)){
	dir.create(new_output_folder)
}

noiters=1000
conc_param=1
cluster_conc = 5

#140514
#if(no.subsamples>=5){
#	density.smooth = 0.005
#}else{
	density.smooth = 0.01
#}

burn.in = 200

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
	write.table(cbind(data[[i]],mutation.copy.number[,i]),paste(samplename,subsamples[[run]][i],"_muts_withAllSubclonalCNinfo_andMutationCopyNumber.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c(names(data[[i]]),"MCN"))
}

mutation.copy.number = mutation.copy.number[,subsample.indices]
copyNumberAdjustment = copyNumberAdjustment[,subsample.indices]
temp.data = list()
for(i in 1:length(subsample.indices)){
	temp.data[[i]] = data[[subsample.indices[i]]]
}
data = temp.data


S.i = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_states.csv",sep=""),row.names=1)
V.h = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_stickbreaking_weights.csv",sep=""),row.names=1)
pi.h = read.csv(paste(output_folder,"/",samplename,"_2D_iters",noiters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_discreteMutationCopyNumbers.csv",sep=""),row.names=1)
pi.h = array(data.matrix(pi.h),c(noiters,30,no.subsamples))

no.subsamples = length(subsample.indices)
pi.h=pi.h[,,subsample.indices]

GS.data = list(S.i=S.i,V.h=V.h,pi.h=pi.h)

density.out = Gibbs.subclone.density.est(mutation.copy.number/copyNumberAdjustment,GS.data,density.smooth,burn.in+1,noiters,max.burden = 1.5)
save(density.out,file=paste(samplename,"_",paste(subsample.indices,collapse="_"),"_multidimensionalDensityInfo_",density.smooth,".RData",sep=""))

range = density.out$range
gridsize = density.out$gridsize
median.density = density.out$median.density
lower.CI = density.out$lower.CI

getHypercubeIndices<-function(gridsize,lastMin,hypercube.size){
	indices = array(0,(2*hypercube.size+1)^length(lastMin))
	pos.within.hypercube = array(0,c((2*hypercube.size+1)^length(lastMin),length(lastMin)))
	for(i in 1:length(lastMin)){
		pos.within.hypercube[,i]=rep(0:(2*hypercube.size),each=(2*hypercube.size+1)^(i-1), times = (2*hypercube.size+1)^(length(lastMin)-i))
	}
	indices = pos.within.hypercube[,1] + lastMin[1]
	for(i in 2:length(lastMin)){
		indices = indices + (lastMin[i]+pos.within.hypercube[,i]-1) * prod(gridsize[1:(i-1)])
	}
	return(indices)
}

hypercube.size = 2
localMins = array(NA,c(0,no.subsamples))
lastMin = rep(1,no.subsamples)
if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
	localMins = rbind(localMins,lastMin)
}

getNextHyperCube<-function(gridsize,lastMin,hypercube.size){
	current.dimension = length(gridsize)
	lastMin[current.dimension] = lastMin[current.dimension] + 1
	while(T){
		if(lastMin[current.dimension]==gridsize[current.dimension] - 2 * hypercube.size + 1){
			if(current.dimension==1){
			return(NULL)
		}
		lastMin[current.dimension] = 1
		current.dimension = current.dimension-1
		lastMin[current.dimension] = lastMin[current.dimension] + 1
		}else{
			return(lastMin)
		}
	}
}

localMins = array(NA,c(0,no.subsamples))
above95confidence = NULL
lastMin = rep(1,no.subsamples)
if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
	localMins = rbind(localMins,lastMin)
	above95confidence = c(above95confidence,lower.CI[rbind(lastMin+hypercube.size)]>0)
}

while(!is.null(lastMin)){
	if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
		localMins = rbind(localMins,lastMin)
		print("local maximum indices:")
		print(lastMin)
		above95confidence = c(above95confidence,lower.CI[rbind(lastMin+hypercube.size)]>0)
	}
	lastMin = getNextHyperCube(gridsize,lastMin,hypercube.size)
}

localMins = localMins + hypercube.size
localOptima = array(rep(range[,1],each=no.subsamples),dim(localMins))
localOptima = localOptima + array(rep((range[,2] - range[,1])/(gridsize-1),each=no.subsamples),dim(localMins)) * localMins

print("localMins")
print(localMins)
print("localOptima")
print(localOptima)

write.table(cbind(localOptima,above95confidence),paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_localMultidimensionalOptima_",density.smooth,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = c(paste(samplename,subsamples[[run]][subsample.indices],sep=""),"above95percentConfidence"))
write.table(localOptima[above95confidence,],paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_localHighConfidenceMultidimensionalOptima_",density.smooth,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = paste(samplename,subsamples[[run]][subsample.indices],sep=""))

no.optima = nrow(localOptima)

if(no.optima>1){	
	stepsize = (range[,2]-range[,1])/(gridsize-1)
	peak.indices = round((localOptima - rep(range[,1],each = nrow(localOptima))) / rep(stepsize,each = nrow(localOptima)))
	
	#peak.heights are not currently used, but they could be used to construct a mixture model
	#and then estimate posterior probs of belonging to each cluster
	peak.heights = median.density[peak.indices]
	
	#if T, j is 'above' i, relative to the plane through the origin
	vector.direction = array(NA,c(no.optima,no.optima))
		
	boundary = array(NA,c(no.optima,no.optima))
	vector.length = array(NA,c(no.optima,no.optima))
	plane.vector = array(NA,c(no.optima,no.optima,no.subsamples+1))
	for(i in 1:(no.optima-1)){
		for(j in (i+1):no.optima){
			#coefficients describe a plane ax + by + cz + d = 0,
			#perpendicular to the line between optimum i and optimum j, passing through optimum i
			plane.vector[i,j,1:no.subsamples] = localOptima[j,] - localOptima[i,]
			plane.vector[i,j,no.subsamples+1] = -sum(plane.vector[i,j,1:no.subsamples] * localOptima[i,])
			
			#not sure this is needed - it may always be true
			vector.direction[i,j] = sum(plane.vector[i,j,1:no.subsamples] * localOptima[j,]) > sum(plane.vector[i,j,1:no.subsamples] * localOptima[i,])
			#this is needed, in order to normalise the distance
			vector.length[i,j] = sqrt(sum(plane.vector[i,j,1:no.subsamples]^2))
			
			longest.dimension = which.max(abs(plane.vector[i,j,1:no.subsamples])/stepsize)			
			no.steps = max(abs(plane.vector[i,j,1:no.subsamples])/stepsize) - 1
			#how long is each step?
			Euclidean.stepsize = sqrt(sum((plane.vector[i,j,1:no.subsamples])^2)) / (no.steps + 1)
			
			step.coords = array(NA,c(no.steps,no.subsamples))
			step.coords = t(sapply(1:no.steps,function(o,x){o[i,] + x * (o[j,] - o[i,]) / (no.steps + 1)},o=localOptima))
			step.indices = round((step.coords - rep(range[,1],each = nrow(step.coords))) / rep(stepsize,each = nrow(step.coords)))
			densities.on.line = median.density[step.indices]
			min.indices = which(densities.on.line == min(densities.on.line))
			#what distance along the line between a pair of optima do we have to go to reach the minimum density
			boundary[i,j] = (max(min.indices) + min(min.indices))/2 * Euclidean.stepsize
			#boundary[j,i] = Euclidean.stepsize * (no.steps+1) - boundary[i,j]
		}
	}
	#boundary = boundary - plane.vector[,,no.subsamples + 1] # make distance relative to a plane through the origin
	boundary = boundary - plane.vector[,,no.subsamples + 1] / vector.length # 020714 - normalise adjustment
	
	no.muts = nrow(mutation.copy.number)
	
	mutation.preferences = array(0,c(no.muts,no.optima))
	
	sampledIters = (burn.in + 1) : noiters
	#don't use the intitial state
	sampledIters = sampledIters[sampledIters!=1]		
	if(length(sampledIters) > 1000){
		sampledIters=floor(post.burn.in.start + (1:1000) * (noiters - burn.in)/1000)			
	}

	if(choose.index==1){
		save.image(file="parallelClusteringBeforeClusterAssignment.RData")
	}

	S.i = data.matrix(S.i)
	for(s in sampledIters){
		temp.preferences = array(0,c(no.muts,no.optima))
		for(c in unique(S.i[s,])){
			for(i in 1:(no.optima-1)){
				for(j in (i+1):no.optima){
					distance.from.plane = sum(pi.h[s,c,] * plane.vector[i,j,1:no.subsamples]) / vector.length[i,j]
					if(distance.from.plane<=boundary[i,j]){
						if(vector.direction[i,j])
						{
							temp.preferences[S.i[s,]==c,i] = temp.preferences[S.i[s,]==c,i] + 1
						}else{
							temp.preferences[S.i[s,]==c,j] = temp.preferences[S.i[s,]==c,j] + 1
						}
					}else{
						if(vector.direction[i,j])
						{
							temp.preferences[S.i[s,]==c,j] = temp.preferences[S.i[s,]==c,j] + 1
						}else{
							temp.preferences[S.i[s,]==c,i] = temp.preferences[S.i[s,]==c,i] + 1
						}
					}
				}
			}			
		}
		iter.preferences = t(sapply(1:no.muts,function(p,i){as.integer(p[i,]==max(p[i,])) / sum(p[i,]==max(p[i,]))},p=temp.preferences))
		mutation.preferences = mutation.preferences + iter.preferences
	}
	mutation.preferences = mutation.preferences / length(sampledIters)
	most.likely.cluster = sapply(1:no.muts,function(m,i){which.max(m[i,])},m=mutation.preferences)

	#get confidence intervals and median
	#subclonal.fraction = data.matrix(mutation.copy.number/copyNumberAdjustment)
	#no.perms = 10000
	#sampled.vals = array(0,c(no.perms,no.optima,no.subsamples))
	#no.muts.per.cluster = array(0,c(no.perms,no.optima))
	#for(p in 1:no.muts){
	#	print(p)
	#	sampled.cluster = sample(1:no.optima,no.perms,mutation.preferences[p,],replace=T)
	#	for(c in unique(sampled.cluster)){
	#		sampled.vals[sampled.cluster == c,c,] = sampled.vals[sampled.cluster == c,c,] + rep(subclonal.fraction[p,],each = sum(sampled.cluster == c))
	#		no.muts.per.cluster[sampled.cluster == c,c] = no.muts.per.cluster[sampled.cluster == c,c] + 1
	#	}
	#}
	#sampled.vals = sampled.vals/rep(no.muts.per.cluster,times = no.subsamples)
	#quantiles = array(NA, c(no.optima,no.subsamples,3))
	#for(c in 1:no.optima){
	#	for(s in 1:no.subsamples){
	#		quantiles[c,s,] = quantile(sampled.vals[,c,s],probs=c(0.025,0.5,0.975),na.rm=T)
	#	}
	#}
	#new method - 180714
 	#quantiles = array(NA, c(no.optima,no.subsamples,3))
	#sampled.thetas = list()
	#totals = table(factor(most.likely.cluster,levels = 1:no.optima))
	#for(i in 1:no.optima){
	#	sampled.thetas[[i]] = array(NA,c(length(sampledIters)*totals[i],no.subsamples))
	#	for(s in 1:length(sampledIters)){
	#		sampled.thetas[[i]][((s-1)*totals[i]+1):(s*totals[i]),] = pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],]
	#	}
	#	for(s in 1:no.subsamples){
	#		quantiles[i,s,] = quantile(sampled.thetas[[i]][,s],probs=c(0.025,0.5,0.975),na.rm=T)
	#	}
	#}
	#new method - 210714 - should be intermediate between previous methods
 	quantiles = array(NA, c(no.optima,no.subsamples,3))
	sampled.thetas = list()
	totals = table(factor(most.likely.cluster,levels = 1:no.optima))
	for(i in 1:no.optima){
		sampled.thetas[[i]] = array(NA,c(length(sampledIters),totals[i],no.subsamples))
		for(s in 1:length(sampledIters)){
		if(sum(most.likely.cluster==i)>1){
			sampled.thetas[[i]][s,,] = pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],]
		}else{
			sampled.thetas[[i]][s,,] = array(pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],],c(totals[i],no.subsamples))
		}
		}
		for(s in 1:no.subsamples){
			median.sampled.vals = sapply(1:length(sampledIters),function(x){median(sampled.thetas[[i]][x,,s])})
			quantiles[i,s,] = quantile(median.sampled.vals,probs=c(0.025,0.5,0.975),na.rm=T)
		}
	}
	
	CIs = array(quantiles[,,c(1,3)],c(no.optima,no.subsamples*2))
	CIs = CIs[,rep(1:no.subsamples,each=2) + rep(c(0,no.subsamples),no.subsamples)]
	out = cbind(data[[1]][,1:2],mutation.preferences,most.likely.cluster)
	names(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")

	write.table(cbind(quantiles[,,2],colSums(mutation.preferences),table(factor(most.likely.cluster,levels = 1:no.optima))),paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_optimaInfo_",density.smooth,".txt",sep=""),col.names = c(paste(samplename,subsamples[[run]][subsample.indices],sep=""),"estimated.no.of.mutations","no.of.mutations.assigned"),row.names=F,sep="\t",quote=F)		
	write.table(out,paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and cluster_info_",density.smooth,".txt",sep=""),sep="\t",row.names=F,quote=F)
	write.table(CIs,paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_confInts_",density.smooth,".txt",sep=""),col.names = paste(rep(paste(samplename,subsamples[[run]][subsample.indices],sep=""),each=2),rep(c(".lower.CI",".upper.CI"),no.subsamples),sep=""),row.names=F,sep="\t",quote=F)
}else{
	most.likely.cluster = rep(1,no.muts)
}


pdf(paste(new_output_folder,"/","most_likely_cluster_assignment_",samplename,"_",paste(subsample.indices,collapse="_"),"_",density.smooth,".pdf",sep=""),height=4,width=4)
#its hard to distinguish more than 8 different colours
max.cols = 8
cols = rainbow(min(max.cols,no.optima))
plot.data = mutation.copy.number/copyNumberAdjustment
plot.data[is.na(plot.data)]=0
for(i in 1:(no.subsamples-1)){
	for(j in (i+1):no.subsamples){
		plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,(subsamples[[run]][subsample.indices])[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[[run]][subsample.indices][j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
		for(n in 1:no.optima){
			pch=20 + floor((n-1)/max.cols)
			#pch is not implmeneted above 25
			if(pch>25){
				pch=pch-20
			}
			points(plot.data[,i][most.likely.cluster==n],plot.data[,j][most.likely.cluster==n],col=cols[(n-1) %% max.cols + 1],pch=pch)
		}
		#legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
		#legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.optima,col=cols[(0:(no.optima-1)) %% max.cols + 1],pch=20 + floor((0:(no.optima-1))/max.cols),cex=1)
		pch=20 + floor((0:(no.optima-1))/max.cols)
		pch[pch>25] = pch[pch>25]-20
		legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.optima,col=cols[(0:(no.optima-1)) %% max.cols + 1],pch=pch,cex=1)

	}
}	
dev.off()

q(save="no")
