setwd("C:/dog")
source("interconvertMutationBurdens.R")

#############################################24T################################################
cellularity.24T = 0.87
CN.data.24T = read.table("ASCAT_CopyNumberSegments_1June2013_24T.txt",header=T,sep="\t",stringsAsFactors=F)
data.24T = read.table("LOH_variants_with_timing_24T.txt",sep="\t",header=T,stringsAsFactors=F)

nucleotides = c("A","C","G","T")

total.snp.rate.24T = 0.5232 #unknown SNPs per KB
het.snp.rate.24T = 0.2993 #unknown heterozygous SNPs per KB
snp.rate.single.chr.24T = total.snp.rate.24T - het.snp.rate.24T/2

LOH.amplification.timings = list()

for(r in 1:nrow(CN.data.24T)){
	print(r)
	inds = which(data.24T$X.CHR==CN.data.24T$chr[r] & data.24T$POS>=CN.data.24T$start.pos[r] & data.24T$POS<=CN.data.24T$end.pos[r])
	if(CN.data.24T$segmented.minor.CN[r]==0){ #LOH
		if(abs(CN.data.24T$segmented.major.CN[r] - round(CN.data.24T$segmented.major.CN[r]))<=0.2 & round(CN.data.24T$segmented.major.CN[r])>1){
			maxCN = round(CN.data.24T$segmented.major.CN[r])
			totals = array(0,maxCN)
			rounded.vals = round(data.24T$MCN[inds])
			for(c in 1:maxCN){
				totals[c] = sum(rounded.vals==c)
			}
			#discount muts that have occurred on unamplified branches
			if(maxCN>2){
				totals[1] = totals[1] - sum(totals[(maxCN-1):2] * 1:(maxCN-2))
			}
			#040513 - discount SNPs (found on all copies)
			totals[maxCN] = max(totals[maxCN] - (CN.data.24T$end.pos[r] - CN.data.24T$end.pos[r] + 1) * snp.rate.single.chr.24T/1000,0)
			
			#divide mutations that have occurred after all amplifications by the total number of copies
			totals[1] = max(totals[1]/maxCN,0)
			timings = cumsum(totals[maxCN:2])/sum(totals)
			LOH.amplification.timings[[length(LOH.amplification.timings)+1]] = c(CN.data.24T$chr[r],CN.data.24T$start.pos[r],CN.data.24T$end.pos[r],maxCN,timings)
			
			#280513 - get variants within each timed segment
			timings=c(0,timings,1)
			for(t in 1:(length(timings)-1)){
				if(sum(rounded.vals==maxCN-t+1)>0){
					LOH.variants.with.timing = rbind(LOH.variants.with.timing,cbind((data.24T[inds,])[rounded.vals==maxCN-t+1,],timings[t],timings[t+1]))
				}
			}
		}
	}else if(CN.data.24T$segmented.minor.CN[r]>0.8 & CN.data.24T$segmented.minor.CN[r]<1.2){ #one unamplified chromosome
		if(abs(CN.data.24T$segmented.major.CN[r] - round(CN.data.24T$segmented.major.CN[r]))<=0.2 & round(CN.data.24T$segmented.major.CN[r])>1){
			maxCN = round(CN.data.24T$segmented.major.CN[r])
			totals = array(0,maxCN)
			rounded.vals = round(data.24T$MCN[inds])
			for(c in 1:maxCN){
				totals[c] = sum(rounded.vals==c)
			}

			#discount het SNPs on unamplified chr		
			totals[1] = totals[1] - (CN.data.24T$end.pos[r] - CN.data.24T$end.pos[r] + 1) * het.snp.rate.24T/1000
			
			#040513 - discount SNPs found on all copies
			totals[maxCN] = max(totals[maxCN] - (CN.data.24T$end.pos[r] - CN.data.24T$end.pos[r] + 1) * het.snp.rate.24T/1000,0)

			#discount muts that have occurred on unamplified branches
			totals[1] = max(totals[1] - sum(totals[maxCN:2] * 1:(maxCN-1)),0)
						
			#divide mutations that have occurred after all amplifications by the total number of copies
			totals[1] = max(totals[1]/(maxCN+1),0)
			
			timings = cumsum(totals[maxCN:2])/sum(totals)
			amplification.timings.single.chromosome[[length(amplification.timings.single.chromosome)+1]] = c(CN.data.24T$chr[r],CN.data.24T$start.pos[r],CN.data.24T$end.pos[r],maxCN,timings)
		}
	}else if(abs(CN.data.24T$segmented.major.CN[r] - round(CN.data.24T$segmented.major.CN[r]))<=0.2 & round(CN.data.24T$segmented.major.CN[r])>1
	& abs(CN.data.24T$segmented.minor.CN[r] - round(CN.data.24T$segmented.minor.CN[r]))<=0.2 & round(CN.data.24T$segmented.minor.CN[r])>1
	& round(CN.data.24T$segmented.minor.CN[r]) == round(CN.data.24T$segmented.major.CN[r])){
		maxCN = round(CN.data.24T$segmented.major.CN[r])
		totals = array(0,maxCN)
		rounded.vals = round(data.24T$MCN[inds])
		for(c in 1:maxCN){
			totals[c] = sum(rounded.vals==c)
		}
		#discount muts that have occurred on unamplified branches
		if(maxCN>2){
			totals[1] = totals[1] - sum(totals[(maxCN-1):2] * 1:(maxCN-2))
		}
		#040513 - discount SNPs (found on all copies)
		totals[maxCN] = max(totals[maxCN] - (CN.data.24T$end.pos[r] - CN.data.24T$end.pos[r] + 1) * het.snp.rate.24T/1000,0)

		#divide mutations that have occurred after all amplifications by the total number of copies
		totals[1] = max(totals[1]/maxCN,0)
		timings = cumsum(totals[maxCN:2])/sum(totals)
		amplification.timings.both.chromosomes[[length(amplification.timings.both.chromosomes)+1]] = c(CN.data.24T$chr[r],CN.data.24T$start.pos[r],CN.data.24T$end.pos[r],maxCN,timings)
	}
}
#280513 - write LOH variants, with timing
names(LOH.variants.with.timing)[(ncol(LOH.variants.with.timing)-1):ncol(LOH.variants.with.timing)]= c("start.time","end.time")
write.table(LOH.variants.with.timing,"LOH_variants_with_timing_24T.txt",sep="\t",quote=F,row.names=F)

save(LOH.amplification.timings,amplification.timings.single.chromosome,amplification.timings.both.chromosomes,file="timings.24T.Rdata")

lens = sapply(1:length(LOH.amplification.timings),function(t,i){length(t[[i]])},t=LOH.amplification.timings)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-4])
	}
	write.table(arr,paste("LOH.amplification.timings.24T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))

	png(paste("LOH.amplification.timings.24T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}
#combine all first duplications, second duplications...
for(l in unique(lens)){
	indices = which(lens %in% l:max(unique(lens)))
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		if(lens[i]>l){
			arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-c(4,(l+1):lens[i])])
		}else{
			arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-4])
		}
	}

	png(paste("LOH.all.amplification.timings.24T.CN",l-3,".png",sep=""),width=600,height = 400)
	i=l-1
	wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
	wts=wts/sum(wts)
	dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
	h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
	normalised.y = max(h$counts) * dens$y/max(dens$y)
	lines(dens$x,normalised.y,col="red",lwd=3)
	dev.off()
}

lens = sapply(1:length(amplification.timings.single.chromosome),function(t,i){length(t[[i]])},t=amplification.timings.single.chromosome)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,amplification.timings.single.chromosome[[indices[i]]][-4])
	}
	write.table(arr,paste("single.chromosome.amplification.timings.24T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))
	png(paste("single.chromosome.amplification.timings.24T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}

lens = sapply(1:length(amplification.timings.both.chromosomes),function(t,i){length(t[[i]])},t=amplification.timings.both.chromosomes)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,amplification.timings.both.chromosomes[[indices[i]]][-4])
	}
	write.table(arr,paste("both.chromosomes.amplification.timings.24T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))
	png(paste("both.chromosomes.amplification.timings.24T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}

######################################79T################################################################
cellularity.79T = 0.86
CN.data.79T = read.table("/nfs/dog_n_devil/dog/CNVs/complete_copy_number_profiles_1June2013/ASCAT_CopyNumberSegments_1June2013_79T.txt",header=T,sep="\t",stringsAsFactors=F)
#data.79T = read.table("/nfs/dog_n_devil/dog/substitutions/caveman/79T_muts_flagged_PASS_caveman_matchedTo24T_withCN.txt",sep="\t",header=T,stringsAsFactors=F)
data.79T = read.table("79T_muts_flagged_PASS_caveman_withCN.txt",sep="\t",header=T,stringsAsFactors=F)

nucleotides = c("A","C","G","T")
data.79T$REF_BASE = toupper(data.79T$REF_BASE)
data.79T$ALT_BASE = sapply(1:nrow(data.79T),function(d,i){substr(gsub(d$REF_BASE[i],"",(strsplit(d$TOP_GENOTYPE[i],"/")[[1]])[2]),1,1)},d=data.79T)
ref.indices = match(data.79T$REF_BASE,nucleotides)
alt.indices = match(data.79T$ALT_BASE,nucleotides)
ref.count = as.numeric(data.79T[cbind(1:nrow(data.79T),14+ref.indices)])
alt.count = as.numeric(data.79T[cbind(1:nrow(data.79T),14+alt.indices)])
data.79T$allele.frequency = alt.count/(ref.count+alt.count)
data.79T$MCN = mutationBurdenToMutationCopyNumber(data.79T$allele.frequency,data.79T$major.CN + data.79T$minor.CN,cellularity.79T)


#total.snp.rate.79T = 0.3543 #unknown SNPs per KB
total.snp.rate.79T = 0.4826 #unknown SNPs per KB
het.snp.rate.79T = 0.2566 #unknown heterozygous SNPs per KB
snp.rate.single.chr.79T = total.snp.rate.79T - het.snp.rate.79T/2

LOH.amplification.timings = list()
amplification.timings.single.chromosome = list()
amplification.timings.both.chromosomes = list() # for chromosomes that have the same amplified CN
#280513
LOH.variants.with.timing = cbind(data.79T[0,],array(NA,c(0,2)))
for(r in 1:nrow(CN.data.79T)){
	print(r)
	inds = which(data.79T$X.CHR==CN.data.79T$chr[r] & data.79T$POS>=CN.data.79T$start.pos[r] & data.79T$POS<=CN.data.79T$end.pos[r])
	if(CN.data.79T$segmented.minor.CN[r]==0){ #LOH
		if(abs(CN.data.79T$segmented.major.CN[r] - round(CN.data.79T$segmented.major.CN[r]))<=0.2 & round(CN.data.79T$segmented.major.CN[r])>1){
			maxCN = round(CN.data.79T$segmented.major.CN[r])
			totals = array(0,maxCN)
			rounded.vals = round(data.79T$MCN[inds])
			for(c in 1:maxCN){
				totals[c] = sum(rounded.vals==c)
			}
			#discount muts that have occurred on unamplified branches
			if(maxCN>2){
				totals[1] = totals[1] - sum(totals[(maxCN-1):2] * 1:(maxCN-2))
			}
			#040513 - discount SNPs (found on all copies)
			totals[maxCN] = max(totals[maxCN] - (CN.data.79T$end.pos[r] - CN.data.79T$end.pos[r] + 1) * snp.rate.single.chr.79T/1000,0)
			
			#divide mutations that have occurred after all amplifications by the total number of copies
			totals[1] = max(totals[1]/maxCN,0)
			timings = cumsum(totals[maxCN:2])/sum(totals)
			LOH.amplification.timings[[length(LOH.amplification.timings)+1]] = c(CN.data.79T$chr[r],CN.data.79T$start.pos[r],CN.data.79T$end.pos[r],maxCN,timings)
			
			#280513 - get variants within each timed segment
			timings=c(0,timings,1)
			for(t in 1:(length(timings)-1)){
				if(sum(rounded.vals==maxCN-t+1)>0){
					LOH.variants.with.timing = rbind(LOH.variants.with.timing,cbind((data.79T[inds,])[rounded.vals==maxCN-t+1,],timings[t],timings[t+1]))
				}
			}			
		}
	}else if(CN.data.79T$segmented.minor.CN[r]>0.8 & CN.data.79T$segmented.minor.CN[r]<1.2){ #one unamplified chromosome
		if(abs(CN.data.79T$segmented.major.CN[r] - round(CN.data.79T$segmented.major.CN[r]))<=0.2 & round(CN.data.79T$segmented.major.CN[r])>1){
			maxCN = round(CN.data.79T$segmented.major.CN[r])
			totals = array(0,maxCN)
			rounded.vals = round(data.79T$MCN[inds])
			for(c in 1:maxCN){
				totals[c] = sum(rounded.vals==c)
			}

			#discount het SNPs on unamplified chr		
			totals[1] = totals[1] - (CN.data.79T$end.pos[r] - CN.data.79T$end.pos[r] + 1) * het.snp.rate.79T/1000
			
			#040513 - discount SNPs found on all copies
			totals[maxCN] = max(totals[maxCN] - (CN.data.79T$end.pos[r] - CN.data.79T$end.pos[r] + 1) * het.snp.rate.79T/1000,0)

			#discount muts that have occurred on unamplified branches
			totals[1] = max(totals[1] - sum(totals[maxCN:2] * 1:(maxCN-1)),0)
						
			#divide mutations that have occurred after all amplifications by the total number of copies
			totals[1] = max(totals[1]/(maxCN+1),0)
			
			timings = cumsum(totals[maxCN:2])/sum(totals)
			amplification.timings.single.chromosome[[length(amplification.timings.single.chromosome)+1]] = c(CN.data.79T$chr[r],CN.data.79T$start.pos[r],CN.data.79T$end.pos[r],maxCN,timings)
		}
	}else if(abs(CN.data.79T$segmented.major.CN[r] - round(CN.data.79T$segmented.major.CN[r]))<=0.2 & round(CN.data.79T$segmented.major.CN[r])>1
	& abs(CN.data.79T$segmented.minor.CN[r] - round(CN.data.79T$segmented.minor.CN[r]))<=0.2 & round(CN.data.79T$segmented.minor.CN[r])>1
	& round(CN.data.79T$segmented.minor.CN[r]) == round(CN.data.79T$segmented.major.CN[r])){
		maxCN = round(CN.data.79T$segmented.major.CN[r])
		totals = array(0,maxCN)
		rounded.vals = round(data.79T$MCN[inds])
		for(c in 1:maxCN){
			totals[c] = sum(rounded.vals==c)
		}
		#discount muts that have occurred on unamplified branches
		if(maxCN>2){
			totals[1] = totals[1] - sum(totals[(maxCN-1):2] * 1:(maxCN-2))
		}
		#040513 - discount SNPs (found on all copies)
		totals[maxCN] = max(totals[maxCN] - (CN.data.79T$end.pos[r] - CN.data.79T$end.pos[r] + 1) * het.snp.rate.79T/1000,0)

		#divide mutations that have occurred after all amplifications by the total number of copies
		totals[1] = max(totals[1]/maxCN,0)
		timings = cumsum(totals[maxCN:2])/sum(totals)
		amplification.timings.both.chromosomes[[length(amplification.timings.both.chromosomes)+1]] = c(CN.data.79T$chr[r],CN.data.79T$start.pos[r],CN.data.79T$end.pos[r],maxCN,timings)
	}
}
#280513 - write LOH variants, with timing
names(LOH.variants.with.timing)[(ncol(LOH.variants.with.timing)-1):ncol(LOH.variants.with.timing)]= c("start.time","end.time")
write.table(LOH.variants.with.timing,"LOH_variants_with_timing_79T.txt",sep="\t",quote=F,row.names=F)

save(LOH.amplification.timings,amplification.timings.single.chromosome,amplification.timings.both.chromosomes,file="timings.79T.Rdata")

lens = sapply(1:length(LOH.amplification.timings),function(t,i){length(t[[i]])},t=LOH.amplification.timings)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-4])
	}
	write.table(arr,paste("LOH.amplification.timings.79T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))

	png(paste("LOH.amplification.timings.79T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}
#combine all first duplications, second duplications...
for(l in unique(lens)){
	indices = which(lens %in% l:max(unique(lens)))
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		if(lens[i]>l){
			arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-c(4,(l+1):lens[i])])
		}else{
			arr = rbind(arr,LOH.amplification.timings[[indices[i]]][-4])
		}
	}

	png(paste("LOH.all.amplification.timings.79T.CN",l-3,".png",sep=""),width=600,height = 400)
	i=l-1
	wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
	wts=wts/sum(wts)
	dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
	h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
	normalised.y = max(h$counts) * dens$y/max(dens$y)
	lines(dens$x,normalised.y,col="red",lwd=3)
	dev.off()
}

lens = sapply(1:length(amplification.timings.single.chromosome),function(t,i){length(t[[i]])},t=amplification.timings.single.chromosome)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,amplification.timings.single.chromosome[[indices[i]]][-4])
	}
	write.table(arr,paste("single.chromosome.amplification.timings.79T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))
	png(paste("single.chromosome.amplification.timings.79T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}

lens = sapply(1:length(amplification.timings.both.chromosomes),function(t,i){length(t[[i]])},t=amplification.timings.both.chromosomes)
for(l in unique(lens)){
	indices = which(lens %in% l)
	arr = array(NA,c(0,l-1))
	for(i in 1:length(indices)){
		arr = rbind(arr,amplification.timings.both.chromosomes[[indices[i]]][-4])
	}
	write.table(arr,paste("both.chromosomes.amplification.timings.79T.CN",l-3,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=c("chr","start.pos","end.pos",paste("timepoint",1:(l-4),sep="")))
	png(paste("both.chromosomes.amplification.timings.79T.CN",l-3,".png",sep=""),width=600,height = 400*(l-4))
	par(mfrow = c(l-4,1))
	for(i in 4:(l-1)){
		wts = (as.numeric(arr[,i-1]) - as.numeric(arr[,i-2]) + 1)[!is.nan(as.numeric(arr[,i]))]
		wts=wts/sum(wts)
		dens=density(as.numeric(arr[,i]),weights = wts, from=0,to=1,na.rm=T,width=0.1)
		h=hist(as.numeric(arr[,i]),breaks=seq(0,1,0.01),col="blue",xlab="time",main=paste("timepoint ",i-3,sep=""))
		normalised.y = max(h$counts) * dens$y/max(dens$y)
		lines(dens$x,normalised.y,col="red",lwd=3)
	}
	dev.off()
}

