# I would think it would also be better to do GC correction before building the reference normal profile
# tried that, but didn't help much!

source("/nfs/team78pc2/pvl/ASCAT/GCcorrection/ascat.R")
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt")
ascat.bc = ascat.GCcorrect(ascat.bc, "GC_targetedpulldown.txt")
Tumor_LogR = ascat.bc$Tumor_LogR
Tumor_BAF = ascat.bc$Tumor_BAF

tumorfiles = "Tumor"
samples = colnames(Tumor_LogR)
ch = ascat.bc$ch
chrs = ascat.bc$chrs

    for (i in 1:dim(Tumor_LogR)[2]) {
      png(filename = paste(tumorfiles,samples[i],".png",sep=""), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = 20)
      plot(c(1,dim(Tumor_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
      points(Tumor_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
#      points(Tumor_BAF[,i],col="red")
      # plot de-mirrored BAF!!
      points(ifelse(runif(length(Tumor_BAF[,i]))<0.5,Tumor_BAF[,i],1-Tumor_BAF[,i]),col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }


Tumor_LogR_GCcorrected = cbind(ascat.bc$SNPpos,ascat.bc$Tumor_LogR)
write.table(Tumor_LogR_GCcorrected,"Tumor_LogR_GCcorrected.txt",sep="\t",col.names=NA,row.names=T,quote=F)

# do germline genotype prediction:
#ascat.gg = ascat.bc$Tumor_BAF < 0.1 | ascat.bc$Tumor_BAF > 0.9

#ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)

#ascat.plotSegmentedData(ascat.bc)

# changed ascat.R such that ACF> 50%, ploidy always around 2
#ascat.output = ascat.runAscat(ascat.bc, gamma = 1)
