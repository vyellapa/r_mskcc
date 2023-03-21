
setwd("/Volumes/p/AMoritz/SV_project/")
setwd("~/Desktop/bolli_targeted/commpass_cnv/")
cnv<- read.delim("commpass_cnv_new_fm6.txt", sep="\t", stringsAsFactors = F)

setwd("~/Desktop/bolli_targeted/SV_project//")
setwd("/Volumes/p/AMoritz/SV_project/")
setwd("~/Desktop/bolli_targeted/SV_project/")
sv_all<- read.delim("delly_mapq_60_all.txt", stringsAsFactors = F, header=T, sep="\t")
head(sv_all)
sv_art<- sv_all[!sv_all$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv<- sv_art[sv_art$code_driver!="artefact_intronic",]
# FERRAN<- head(sv)
# setwd("~/Desktop/")
# write.table(FERRAN, "ferran_delly.txt", sep="\t", quote=F, row.names = F)
# max(table(sv$Specimen_ID))
setwd("/Volumes/p/PTCL_NOS/PTCL_NOS_WGS/CN_RR_PLOTS/")
setwd("~/Desktop/bolli_targeted/SV_project//")
gene<- read.delim("ptcl_GENE_LIST.txt",sep="\t", stringsAsFactors = F)

setwd("/Volumes/p/AMoritz/SV_project/plots/")
setwd("~/Desktop/bolli_targeted/SV_project/cn_plot/")
library(plotrix)
library(igraph)
library(quantsmooth)
sample_list<- intersect(unique(sv$sample) , unique(cnv$IDA))
chrom<- unique(cnv$seqnames)

which(sample_list == "MMRF_1700_1_BM")
which(sample_list == "MMRF_1700_2_BM")
sv$pos1<- as.numeric(as.character(sv$pos1))
sv$pos2<- as.numeric(as.character(sv$pos2))

sample_list<- c("MMRF_1700_1_BM","MMRF_1700_2_BM" )
setwd("/Volumes/GoogleDrive/My Drive/MSKCC/GRANT/SCOR/figure/new_jan/figure_2/cn_plots_fig_2/")
for(w in (1:length(sample_list)))
{
  
  pdf(sprintf("%s_CN_COMMPASS_plots.pdf", sample_list[w]), height=4, width=12)
  ascat_sample2<- cnv[cnv$IDA == sample_list[w],]
  ascat_sample<- ascat_sample2[,c(1:6)]
  colnames(ascat_sample)<- c("sample","chrom","start","end","tot","min")
  ascat_sample$diff<- ascat_sample$end - ascat_sample$start
  ascat_sample<- ascat_sample[order(ascat_sample$chrom, ascat_sample$start),]
  kk<-sample_list[w]
  brass_sample2<- sv[sv$sample == kk,]
  brass_sample<- brass_sample2[,c("sample","chrom1","pos1","chrom2","pos2","SVTYPE",
                                  "RGENUPS_Strand","RGENDNS_Strand")]
  
  
  for(jj in (1:23))
  {  
    i<- chrom[jj]
    gene_chr<- gene[gene$Chromosome.Name==i,]
    brass_sample_chr <- brass_sample[brass_sample$chrom1 == i | brass_sample$chrom2 ==i,]
    ascat_sample_filter<- ascat_sample[ascat_sample$chrom ==i,]
    if(nrow(ascat_sample_filter)>0){
      ascat_sample_filter_chr<- ascat_sample_filter
    ascat_sample_filter_chr$tot[is.na(ascat_sample_filter_chr$tot)]<-0
    ascat_sample_filter_chr$min[is.na(ascat_sample_filter_chr$min)]<-0
    iArrows <- igraph:::igraph.Arrows
    
    # setwd("/Volumes/fm6/chapman_data/copy_number_plots/all_copy_number_plot_batt//")
    name<-sample_list[w]
    par(mar=c(5,7,2,5), xpd=F)
    par(mfrow=c(1,1))
    plot(1, type="n", xlab="", xlim=c(0, max(ascat_sample_filter_chr$end)), 
         ylim=c(0, max(ascat_sample_filter_chr$tot)+3), 
         ylab="", las=2,bty ="n",
         xaxt="n", yaxt="n", main=paste(kk, "-","Chromosome", i, sep=" "))
    
    
    axis(2, line=0, at = c(0:(max(ascat_sample_filter_chr$tot) +2)),  col.axis = "black",
         labels=c(0:(max(ascat_sample_filter_chr$tot) +2)) ,lwd.ticks =1, cex.axis = 1, las=2)
    mtext("Copy Number", 2, +2.5, at = median(c(0:(max(ascat_sample_filter_chr$tot) +1))), cex=1)
    par(new=TRUE, xpd=T)
    paintCytobands(i, pos = c(0, -0.4), units = "bases", width = 0.3, cex.leg = 0.1,
                   bands = "major", legend=F, length.out =max(ascat_sample_filter_chr$end)) 
    
    rect(par("usr")[1],par("usr")[3],par("usr")[2],0.5,col = "gray96", border = FALSE, lwd = 0.001)
    
    col<- rep(c("gray90","gray96"), max(ascat_sample_filter_chr$tot)+2)
    for(z in (1:(max(ascat_sample_filter_chr$tot)+2)))
    {
      
      rect(par("usr")[1],z-0.5,par("usr")[2],z+0.5,col = col[z], border = FALSE, lwd = 0.001)
    }
    #rect(par("usr")[1],1.5,par("usr")[2],2.5,col = "gray90", border = FALSE, lwd = 0.001)
    #rect(par("usr")[1],2.5,par("usr")[2],3.5,col = "gray96", border = FALSE, lwd = 0.001)
    #rect(par("usr")[1],3.5,par("usr")[2],4.5,col = "gray90", border = FALSE, lwd = 0.001)
    #rect(par("usr")[1],4.5,par("usr")[2],5.5,col = "gray96", border = FALSE, lwd = 0.001)
    #rect(par("usr")[1],5.5,par("usr")[2],6.5,col = "gray90", border = FALSE, lwd = 0.001)
    #rect(par("usr")[1],6.5,par("usr")[2],max(ascat_sample_filter_chr$tot) +1.5,col = "gray96", border = FALSE, lwd = 0.001)
    
    segments(ascat_sample_filter_chr$start, ascat_sample_filter_chr$tot, ascat_sample_filter_chr$end,ascat_sample_filter_chr$tot,
             lwd=2.5, col="black")
    segments(ascat_sample_filter_chr$start, ascat_sample_filter_chr$min, ascat_sample_filter_chr$end,ascat_sample_filter_chr$min,
             lwd=2.5, col="darkgoldenrod1", lty = 2)
    
    if(nrow(brass_sample_chr)>0)
    {  
      for(j in (1:nrow(brass_sample_chr)))
      {
        if(brass_sample_chr$SVTYPE[j] == "DEL")
        {
          segments(brass_sample_chr$pos1[j], 0, brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1, col="firebrick2", lwd=1)
          segments(brass_sample_chr$pos2[j], 0, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1, col="firebrick2", lwd=1)
          grid.curve(brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1)
          iArrows(brass_sample_chr$pos1[j],max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos2[j],
                  max(ascat_sample_filter_chr$tot)+1,h.lwd=1, sh.lwd=1, sh.col="firebrick2", curve=0.00000005 , width=1, size=0.002)
        }else
        {
          if(brass_sample_chr$SVTYPE[j] == "DUP")
          {
            segments(brass_sample_chr$pos1[j], 0, brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1.5, col="green", lwd=1)
            segments(brass_sample_chr$pos2[j], 0, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1.5, col="green", lwd=1)
            grid.curve(brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1.5, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1.5)
            iArrows(brass_sample_chr$pos1[j],max(ascat_sample_filter_chr$tot)+1.5, brass_sample_chr$pos2[j], 
                    max(ascat_sample_filter_chr$tot)+1.5,h.lwd=1, sh.lwd=1, sh.col="green", curve=0.00000005 , width=1, size=0.002)
          }else
          {
            if(brass_sample_chr$SVTYPE[j] == "INV")
            {
              segments(brass_sample_chr$pos1[j], 0, brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+2, col="dodgerblue", lwd=1)
              segments(brass_sample_chr$pos2[j], 0, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+2, col="dodgerblue", lwd=1)
              grid.curve(brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+2, brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+2)
              iArrows(brass_sample_chr$pos1[j],max(ascat_sample_filter_chr$tot)+2, brass_sample_chr$pos2[j], 
                      max(ascat_sample_filter_chr$tot)+2,h.lwd=1, sh.lwd=1, sh.col="dodgerblue", curve=0.00000005 , width=1, size=0.002)
            }else
            {
              if(brass_sample_chr$chrom1[j] ==i)
              {
                segments(brass_sample_chr$pos1[j], 0, brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1, col="black", lwd=1)
                if(brass_sample_chr$RGENUPS_Strand[j] == "+")
                {
                  segments(brass_sample_chr$pos1[j],  max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos1[j]+max(ascat_sample_filter_chr$end)/100, 
                           max(ascat_sample_filter_chr$tot)+2, col="black", lwd=1)
                  text(brass_sample_chr$pos1[j]+max(ascat_sample_filter_chr$end)/100, max(ascat_sample_filter_chr$tot)+2.3, brass_sample_chr$CHR2[j], cex=1.2)
                }else
                {
                  segments(brass_sample_chr$pos1[j], max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos1[j]-max(ascat_sample_filter_chr$end)/100, 
                           max(ascat_sample_filter_chr$tot)+2, col="black", lwd=1)
                  text(brass_sample_chr$pos1[j]-max(ascat_sample_filter_chr$end)/100, max(ascat_sample_filter_chr$tot)+2.3, brass_sample_chr$CHR2[j], cex=1.2)
                }}else
                {
                  segments(brass_sample_chr$pos2[j], 0, brass_sample_chr$pos2[j], 
                           max(ascat_sample_filter_chr$tot)+1, col="black", lwd=1)
                  if(brass_sample_chr$RGENUPS_Strand[j] == "+")
                  {
                    segments(brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos2[j]+max(ascat_sample_filter_chr$end)/100, 
                             max(ascat_sample_filter_chr$tot)+2, col="black", lwd=1)
                    text(brass_sample_chr$pos2[j]+max(ascat_sample_filter_chr$end)/100, max(ascat_sample_filter_chr$tot)+2.3,brass_sample_chr$CHROM[j], cex=1.2)
                  }else
                  {
                    segments(brass_sample_chr$pos2[j], max(ascat_sample_filter_chr$tot)+1, brass_sample_chr$pos2[j]-max(ascat_sample_filter_chr$end)/100, 
                             max(ascat_sample_filter_chr$tot)+2, col="black", lwd=1)
                    text(brass_sample_chr$pos2[j]-max(ascat_sample_filter_chr$end)/100, max(ascat_sample_filter_chr$tot)+2.3, brass_sample_chr$CHROM[j], cex=1.2)
                  }}
            }
          }
        }
      }}

    }}
  dev.off()
}


#graphics.off()
#plot(1,1)

