Tumor_LogR_file = "Normal_LogR.txt"
Tumor_BAF_file = "Normal_BAF.txt"

Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")

SNPpos <- Tumor_LogR[,1:2]
chrs=as.vector(unique(SNPpos[,1]))

Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]


  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))  
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]


tumorfiles = "Normal"
samples = colnames(Tumor_LogR)

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

