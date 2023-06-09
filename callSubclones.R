##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpBattenberg.
#
# cgpBattenberg is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

args=commandArgs(TRUE)
lib_path<-toString(args[1])
impute_input_file<-toString(args[2])
is.male<-as.logical(args[3])
impute.info = read.table(impute_input_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
chr_names=unique(impute.info[,1])

source(paste(lib_path,"orderEdges.R",sep="/"))

start.file = toString(args[4])
rho_psi_info = read.table(paste(start.file,"_rho_and_psi.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

#rho = rho_psi_info$rho[which(rho_psi_info$is.best)] # rho = tumour percentage (called tp in previous versions)
#psit = rho_psi_info$psi[which(rho_psi_info$is.best)] # psi of tumour cells
#always use best solution from grid search - reference segment sometimes gives strange results
rho = rho_psi_info$rho[rownames(rho_psi_info)=="FRAC_GENOME"] # rho = tumour percentage (called tp in previous versions)
psit = rho_psi_info$psi[rownames(rho_psi_info)=="FRAC_GENOME"] # psi of tumour cells

psi = rho*psit + 2 * (1-rho) # psi of all cells

#gamma is the shrinkage parameter for a particular platform
gamma=1

#segmentation.gamma is the gamma used in PCF
segmentation.gamma=NA
if(length(args)>=5){
	segmentation.gamma = as.integer(args[5])
	print(paste("segmentation.gamma=",segmentation.gamma,sep=""))
	if(length(args)>=6){
		gamma = as.numeric(args[6])
		print(paste("platform gamma=",gamma,sep=""))
	}
}

siglevel = 0.05
maxdist = 0.02 # This is distance between the observed BAF and the test BAF level it is compared.
	       # Setting this to 0.02 means that the observed BAF has to be >(closest BAF level) + 0.02
	       # If this is not the case, the segment will be assumed to be clonal
noperms = 1000

BAFvals = read.table(paste(start.file,".BAFsegmented.txt",sep=""),sep="\t",header=T, row.names=1)


BAF = BAFvals[,3]
names(BAF)=rownames(BAFvals)
BAFphased = BAFvals[,4]
names(BAFphased)=rownames(BAFvals)
BAFseg = BAFvals[,5]
names(BAFseg)=rownames(BAFvals)

SNPpos = BAFvals[,c(1,2)]

#assume filename ends with mutantLogR.tab
LogRfile = paste(start.file,"_mutantLogR.tab",sep="")

LogRvals = read.table(LogRfile,sep="\t", header=T, row.names=1)

ctrans = c(1:length(chr_names))
names(ctrans)=chr_names
ctrans.logR = c(1:length(chr_names))
names(ctrans.logR)=chr_names

#note the as.vector to avoid problems
LogRpos = as.vector(ctrans.logR[as.vector(LogRvals[,1])]*1000000000+LogRvals[,2])
names(LogRpos) = rownames(LogRvals)
BAFpos = as.vector(ctrans[as.vector(BAFvals[,1])]*1000000000+BAFvals[,2])
names(BAFpos) = rownames(BAFvals)

#DCW 240314
switchpoints = c(0,which(BAFseg[-1] != BAFseg[-(length(BAFseg))] | BAFvals[-1,1] != BAFvals[-nrow(BAFvals),1]),length(BAFseg))
BAFlevels = BAFseg[switchpoints[-1]]

pval = NULL

BAFpvals = vector(length=length(BAFseg))

subcloneres = NULL

for (i in 1:length(BAFlevels)) {
  l = BAFlevels[i]

  #280814 - make sure that BAF>=0.5, otherwise nMajor and nMinor may be the wrong way around
  l=max(l,1-l)

  #DCW 240314
  BAFke = BAFphased[(switchpoints[i]+1):switchpoints[i+1]]

  startpos = min(BAFpos[names(BAFke)])
  endpos = max(BAFpos[names(BAFke)])
  chrom = names(ctrans[floor(startpos/1000000000)])
  LogR = mean(LogRvals[LogRpos>=startpos&LogRpos<=endpos & !is.infinite(LogRvals[,3]),3],na.rm=T)
  #if we don't have a value for LogR, fill in 0
  if (is.na(LogR)) {
    LogR = 0
  }
  nMajor = (rho-1+l*psi*2^(LogR/gamma))/rho
  nMinor = (rho-1+(1-l)*psi*2^(LogR/gamma))/rho

  # to make sure we're always in a positive square:
  #if(nMajor < 0) {
  #  nMajor = 0.01
  #}
  #if(nMinor < 0) {
  #  nMinor = 0.01
  #}
	#DCW - increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
	if(nMinor<0){
		if(l==1){
			#avoid calling infinite copy number
			nMajor = 1000
		}else{
			nMajor = nMajor + l * (0.01 - nMinor) / (1-l)
		}
		nMinor = 0.01
	}

  # note that these are sorted in the order of ascending BAF:
  nMaj = c(floor(nMajor),ceiling(nMajor),floor(nMajor),ceiling(nMajor))
  nMin = c(ceiling(nMinor),ceiling(nMinor),floor(nMinor),floor(nMinor))
  x = floor(nMinor)
  y = floor(nMajor)

  # total copy number, to determine priority options
  ntot = nMajor + nMinor

  levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
  #problem if rho=1 and nMaj=0 and nMin=0
  levels[nMaj==0 & nMin==0] = 0.5

  #whichclosestlevel = which.min(abs(levels-l))
  ## if 0.5 and there are multiple options, finetune, because a random option got chosen
  #if (levels[whichclosestlevel]==0.5 && levels[2]==0.5 && levels[3]==0.5) {
  #  whichclosestlevel = ifelse(ntot>x+y+1,2,3)
  #}

  #DCW - just test corners on the nearest edge to determine clonality
  #If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
  all.edges = orderEdges(levels, l, ntot,x,y)
  nMaj.test = all.edges[1,c(1,3)]
  nMin.test = all.edges[1,c(2,4)]
  test.levels = (1-rho+rho*nMaj.test)/(2-2*rho+rho*(nMaj.test+nMin.test))
  whichclosestlevel.test = which.min(abs(test.levels-l))

  #270713 - problem caused by segments with constant BAF (usually 1 or 2)
  if(sd(BAFke)==0){
	  pval[i]=0
  }else{
  	#pval[i] = t.test(BAFke,alternative="two.sided",mu=levels[whichclosestlevel])$p.value
  	pval[i] = t.test(BAFke,alternative="two.sided",mu=test.levels[whichclosestlevel.test])$p.value
  }
  #if(min(abs(l-levels))<maxdist) {
  #if(min(abs(l-test.levels))<maxdist) {
  if(abs(l-test.levels[whichclosestlevel.test])<maxdist) {
    pval[i]=1
  }

  #BAFpvals[BAFseg==l]=pval[i]
  #DCW 240314
  BAFpvals[(switchpoints[i]+1):switchpoints[i+1]]=pval[i]

  # if the difference is significant, call subclone level
  if(pval[i] <= siglevel) {

	all.edges = orderEdges(levels, l, ntot,x,y)
	#switch order, so that negative copy numbers are at the end
	na.indices = which(is.na(rowSums(all.edges)))
	if(length(na.indices)>0){
		all.edges = rbind(all.edges[-na.indices,],all.edges[na.indices,])
	}
	nMaj1 = all.edges[,1]
	nMin1 = all.edges[,2]
	nMaj2 = all.edges[,3]
	nMin2 = all.edges[,4]

    tau = (1 - rho + rho * nMaj2 - 2 * l * (1 - rho) - l * rho * (nMin2 + nMaj2)) / (l * rho * (nMin1 + nMaj1) - l * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2)
    sdl = sd(BAFke,na.rm=T)/sqrt(sum(!is.na(BAFke)))
    sdtau = abs((1 - rho + rho * nMaj2 - 2 * (l+sdl) * (1 - rho) - (l+sdl) * rho * (nMin2 + nMaj2)) / ((l+sdl) * rho * (nMin1 + nMaj1) - (l+sdl) * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2) - tau) / 2 +
            abs((1 - rho + rho * nMaj2 - 2 * (l-sdl) * (1 - rho) - (l-sdl) * rho * (nMin2 + nMaj2)) / ((l-sdl) * rho * (nMin1 + nMaj1) - (l-sdl) * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2) - tau) / 2

    sdtaubootstrap = vector(length=length(tau),mode="numeric")
    tau25 = vector(length=length(tau),mode="numeric")
    tau975 = vector(length=length(tau),mode="numeric")

    #bootstrapping, adapted from David Wedge
    for (option in 1:length(tau)) {
      nMaj1o = nMaj1[option]
      nMin1o = nMin1[option]
      nMaj2o = nMaj2[option]
      nMin2o = nMin2[option]
      set.seed(as.integer(Sys.time()))
      permFraction=vector(length=noperms,mode="numeric")
      for(j in 1:noperms){
        permBAFs=sample(BAFke,length(BAFke),replace=T)
        permMeanBAF=mean(permBAFs)
        permFraction[j] = (1 - rho + rho * nMaj2o - 2 * permMeanBAF * (1 - rho) - permMeanBAF * rho * (nMin2o + nMaj2o)) / (permMeanBAF * rho * (nMin1o + nMaj1o) - permMeanBAF * rho * (nMin2o + nMaj2o) - rho * nMaj1o + rho * nMaj2o)
      }
      orderedFractions=sort(permFraction)
      sdtaubootstrap[option] = sd(permFraction)
      tau25[option] = orderedFractions[25]
      tau975[option] = orderedFractions[975]
    }

    subcloneres = rbind(subcloneres,c(chrom,startpos-floor(startpos/1000000000)*1000000000,
      endpos-floor(endpos/1000000000)*1000000000,l,pval[i],LogR,ntot,
      nMaj1[1],nMin1[1],tau[1],nMaj2[1],nMin2[1],1-tau[1],sdtau[1],sdtaubootstrap[1],tau25[1],tau975[1],
      nMaj1[2],nMin1[2],tau[2],nMaj2[2],nMin2[2],1-tau[2],sdtau[2],sdtaubootstrap[2],tau25[2],tau975[2],
      nMaj1[3],nMin1[3],tau[3],nMaj2[3],nMin2[3],1-tau[3],sdtau[3],sdtaubootstrap[3],tau25[3],tau975[3],
      nMaj1[4],nMin1[4],tau[4],nMaj2[4],nMin2[4],1-tau[4],sdtau[4],sdtaubootstrap[4],tau25[4],tau975[4],
      nMaj1[5],nMin1[5],tau[5],nMaj2[5],nMin2[5],1-tau[5],sdtau[5],sdtaubootstrap[5],tau25[5],tau975[5],
      nMaj1[6],nMin1[6],tau[6],nMaj2[6],nMin2[6],1-tau[6],sdtau[6],sdtaubootstrap[6],tau25[6],tau975[6]))
  }else {
    #if called as clonal, use the best corner from the nearest edge
	subcloneres = rbind(subcloneres,c(chrom,startpos-floor(startpos/1000000000)*1000000000,
	endpos-floor(endpos/1000000000)*1000000000,l,pval[i],LogR,ntot,
	nMaj.test[whichclosestlevel.test],nMin.test[whichclosestlevel.test],1,rep(NA,57)))

  }
  colnames(subcloneres) = c("chr","startpos","endpos","BAF","pval","LogR","ntot",
    "nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","SDfrac_A","SDfrac_A_BS","frac1_A_0.025","frac1_A_0.975",
    "nMaj1_B","nMin1_B","frac1_B","nMaj2_B","nMin2_B","frac2_B","SDfrac_B","SDfrac_B_BS","frac1_B_0.025","frac1_B_0.975",
    "nMaj1_C","nMin1_C","frac1_C","nMaj2_C","nMin2_C","frac2_C","SDfrac_C","SDfrac_C_BS","frac1_C_0.025","frac1_C_0.975",
    "nMaj1_D","nMin1_D","frac1_D","nMaj2_D","nMin2_D","frac2_D","SDfrac_D","SDfrac_D_BS","frac1_D_0.025","frac1_D_0.975",
    "nMaj1_E","nMin1_E","frac1_E","nMaj2_E","nMin2_E","frac2_E","SDfrac_E","SDfrac_E_BS","frac1_E_0.025","frac1_E_0.975",
    "nMaj1_F","nMin1_F","frac1_F","nMaj2_F","nMin2_F","frac2_F","SDfrac_F","SDfrac_F_BS","frac1_F_0.025","frac1_F_0.975")
}

write.table(subcloneres,paste(start.file,"_subclones.txt",sep=""),quote=F,col.names=NA,row.names=T,sep="\t")


#DCW 211012
sample.name = strsplit(start.file,"/")
sample.name = sample.name[[1]][length(sample.name[[1]])]
sample.name = gsub("_","",sample.name)

for (chr in chr_names) {
  pos = SNPpos[SNPpos[,1]==chr,2]
  #if no points to plot, skip
  if(length(pos)==0){next}
  BAFchr = BAF[SNPpos[,1]==chr]
  BAFsegchr = BAFseg[SNPpos[,1]==chr]
  BAFpvalschr = BAFpvals[SNPpos[,1]==chr]
  LogRchr = LogRvals[LogRvals[,1]==chr,3]
  LogRposke = LogRvals[LogRvals[,1]==chr,2]

  png(filename = paste(start.file,"_subclones_chr",chr,".png",sep=""), width = 2000, height = 2000, res = 200)
  par(mar = c(2.5,2.5,2.5,0.25), cex = 0.4, cex.main=1.5, cex.axis = 1, cex.lab = 1, mfrow = c(2,1))
  plot(c(min(pos)/1000000,max(pos)/1000000),c(-1,1),pch=".",type = "n",
       main = paste(sample.name,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "LogR")
  points(LogRposke/1000000,LogRchr,pch=".",col="grey")
  plot(c(min(pos)/1000000,max(pos)/1000000),c(0,1),pch=".",type = "n",
       main = paste(sample.name,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "BAF (phased)")
  points(pos/1000000,BAFchr,pch=".",col="grey")
  points(pos/1000000,BAFsegchr,pch=19,cex=0.5,col=ifelse(BAFpvalschr>siglevel,"darkgreen","red"))
  points(pos/1000000,1-BAFsegchr,pch=19,cex=0.5,col=ifelse(BAFpvalschr>siglevel,"darkgreen","red"))
  for (i in 1:dim(subcloneres)[1]) {
    if(subcloneres[i,1]==chr) {
      text((as.numeric(subcloneres[i,"startpos"])+as.numeric(subcloneres[i,"endpos"]))/2/1000000,as.numeric(subcloneres[i,"BAF"])-0.04,
        paste(subcloneres[i,"nMaj1_A"],"+",subcloneres[i,"nMin1_A"],": ",100*round(as.numeric(subcloneres[i,"frac1_A"]),3),"%",sep=""),cex = 0.8)
      if(!is.na(subcloneres[i,"nMaj2_A"])) {
        text((as.numeric(subcloneres[i,"startpos"])+as.numeric(subcloneres[i,"endpos"]))/2/1000000,as.numeric(subcloneres[i,"BAF"])-0.08,
          paste(subcloneres[i,"nMaj2_A"],"+",subcloneres[i,"nMin2_A"],": ",100*round(as.numeric(subcloneres[i,"frac2_A"]),3),"%",sep=""), cex = 0.8)
      }
    }
  }
  dev.off()
}

q(save="no")




