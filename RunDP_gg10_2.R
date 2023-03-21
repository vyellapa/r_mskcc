.libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library","/home/yellapav/R/OS7"))

library("DPClust")

#patch to cope with zero allele frequency in some samples
library(R.utils)
source('/ifs/work/leukgen/home/gg10/soft/dpclust_pipeline_v0.2.2/parallel_assignments_gg10/multidimensionalDensityEstimator.R')

args = commandArgs(T)
wdir = args[1]
perm.index = as.numeric(args[2])
dir_full_run = args[3]
sample2purity = args[4]
dp_input_dir = args[5]

setwd(wdir)
perm.info = read.table("permutation_info.txt", sep="\t", header=F, row.names=NULL, stringsAsFactors=F)

donor = perm.info[perm.index,1]
print(donor)
subsamplenames = as.vector(perm.info[perm.index,])
print('Current triple:')
print(subsamplenames)
no.samples = length(subsamplenames)

# Set number of iterations
no.iters=10000
burn.in = 1000
conc_param = 1
cluster_conc = 5
density.smooth = 0.01
# Create DP output directory
output_folder = "DP_output/"
if(!file.exists(output_folder)){
        dir.create(output_folder)
}
output_folder = paste(output_folder, perm.index, sep="/")
if(!file.exists(output_folder)){
        dir.create(output_folder)
}

########################## Read the same dataset that contains the full sample set and that was used to do Gibbs sampling
all_info = read.table(sample2purity, header=T, stringsAsFactors=F)
cellularity = all_info$cellularity
subsample.indeces = NULL
for(name in subsamplenames) {
        subsample.indeces = c(subsample.indeces, match(gsub('\\.','-',as.character(name)), all_info$samplename))
}
dp_info_files = paste(dp_input_dir,all_info$datafile, sep="/")
dataset = load.data(
        dp_info_files, cellularity=cellularity, Chromosome="chr", position="end", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN",
        no.chrs.bearing.mut="no.chrs.bearing.mut", mutation.copy.number="mutation.copy.number", subclonal.fraction="subclonal.fraction", phase="phase",
        is.male=TRUE, is.vcf=FALSE, ref.genome.version="hg19", min.mutreads=0, min.depth=1
)
print("Dataset: #mutations vs #samples")
print(dim(dataset$chromosome))
mutation.copy.number = dataset$mutation.copy.number[,subsample.indeces]
copyNumberAdjustment = dataset$copyNumberAdjustment[,subsample.indeces]
########################## Load Gibbs sampling from the complete run and run multi-dimensional density estimator
load(list.files(path=dir_full_run, pattern="_gsdata.RData", full.names=T))
S.i = GS.data$S.i
V.h = GS.data$V.h
pi.h = GS.data$pi.h
print("Gibbs sampling data:")
print(paste("Iterations: ", nrow(S.i), sep=""))
print(paste("Mutations: ", ncol(S.i), sep=""))
print(paste("Samples: ", dim(pi.h)[3], sep=""))
pi.h = array(data.matrix(pi.h),c(no.iters,30,length(all_info$samplename)))
pi.h=pi.h[,,subsample.indeces]
new.GS.data = list(S.i=S.i,V.h=V.h,pi.h=pi.h)
post.burn.in.start = burn.in + 1
density.out = Gibbs.subclone.density.est.multi(
	burden=mutation.copy.number/copyNumberAdjustment, GS.data=new.GS.data, density.smooth=density.smooth,
	post.burn.in.start=post.burn.in.start, post.burn.in.stop=no.iters,  max.burden=1.5
)
save(density.out, file=paste(output_folder, "/", paste(gsub('-',',',subsamplenames), collapse="_"), "__multidimensionalDensityInfo__", density.smooth,".RData",sep=""))
########################## Identify number of nodes, node positions and assign mutations
no.subsamples = length(subsamplenames)
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

outfile_prefix = paste(output_folder, "/", paste(gsub('-','.',subsamplenames),collapse="_"), sep="")
write.table(
	cbind(localOptima,above95confidence), file=paste(outfile_prefix,"__localMultidimensionalOptima__",density.smooth,".txt",sep=""),
	quote=F, sep="\t", row.names=F, col.names = c(subsamplenames,"above95percentConfidence")
)
write.table(
	localOptima[above95confidence,], file=paste(outfile_prefix,"__localHighConfidenceMultidimensionalOptima__",density.smooth,".txt",sep=""),
	quote=F, sep="\t", row.names=F, col.names = subsamplenames
)

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

        sampledIters = (burn.in + 1) : no.iters
        #don't use the intitial state
        sampledIters = sampledIters[sampledIters!=1]
        if(length(sampledIters) > 1000){
                sampledIters = floor(post.burn.in.start + (1:1000) * (no.iters - burn.in)/1000)
		sampledIters = sampledIters[sampledIters<nrow(S.i)]
        }

        S.i = data.matrix(S.i)
        for(s in sampledIters){
		print(paste("sampledIter = ", s, sep=""))
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
        #       print(p)
        #       sampled.cluster = sample(1:no.optima,no.perms,mutation.preferences[p,],replace=T)
        #       for(c in unique(sampled.cluster)){
        #               sampled.vals[sampled.cluster == c,c,] = sampled.vals[sampled.cluster == c,c,] + rep(subclonal.fraction[p,],each = sum(sampled.cluster == c))
        #               no.muts.per.cluster[sampled.cluster == c,c] = no.muts.per.cluster[sampled.cluster == c,c] + 1
        #       }
        #}
        #sampled.vals = sampled.vals/rep(no.muts.per.cluster,times = no.subsamples)
        #quantiles = array(NA, c(no.optima,no.subsamples,3))
        #for(c in 1:no.optima){
        #       for(s in 1:no.subsamples){
        #               quantiles[c,s,] = quantile(sampled.vals[,c,s],probs=c(0.025,0.5,0.975),na.rm=T)
        #       }
        #}
        #new method - 180714
        #quantiles = array(NA, c(no.optima,no.subsamples,3))
        #sampled.thetas = list()
        #totals = table(factor(most.likely.cluster,levels = 1:no.optima))
        #for(i in 1:no.optima){
        #       sampled.thetas[[i]] = array(NA,c(length(sampledIters)*totals[i],no.subsamples))
        #       for(s in 1:length(sampledIters)){
        #               sampled.thetas[[i]][((s-1)*totals[i]+1):(s*totals[i]),] = pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],]
        #       }
        #       for(s in 1:no.subsamples){
        #               quantiles[i,s,] = quantile(sampled.thetas[[i]][,s],probs=c(0.025,0.5,0.975),na.rm=T)
        #       }
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
        out = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], mutation.preferences, most.likely.cluster)
        colnames(out) = c("chr","start","end",paste("prob.cluster.",1:no.optima,sep=""),"most.likely.cluster")
        write.table(
		cbind(quantiles[,,2],colSums(mutation.preferences),table(factor(most.likely.cluster,levels = 1:no.optima))),
		file=paste(outfile_prefix, "__optimaInfo__",density.smooth,".txt",sep=""),
		col.names = c(subsamplenames,"estimated.no.of.mutations","no.of.mutations.assigned"), row.names=F, sep="\t", quote=F
	)
        write.table(out, paste(outfile_prefix, "__DP_and_cluster_info__",density.smooth,".txt",sep=""), sep="\t", row.names=F, quote=F)
        write.table(
		CIs, file=paste(outfile_prefix, "__confInts__", density.smooth,".txt", sep=""),
		col.names=paste(rep(subsamplenames,each=2),rep(c(".lower.CI",".upper.CI"),no.subsamples),sep=""), row.names=F, sep="\t", quote=F)
}else{
        most.likely.cluster = rep(1,no.muts)
}

q(save="no")


