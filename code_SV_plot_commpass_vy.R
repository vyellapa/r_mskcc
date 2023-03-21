
setwd("~/Desktop/misc/")

cnv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-106917-T2-1-D1-2_subclones.txt",sep="\t",header=TRUE)
cnv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130719-T1-4-D1-2_subclones.txt",sep="\t",header=TRUE)
#cnv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130720-T1-2-D1-2_subclones.txt",sep="\t",header=TRUE)
#cnv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130718-T1-6-D1-2_subclones.txt",sep="\t",header=TRUE)
cnv=cnv[,c(1:14)]

svv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130719-T1-4-D1-2_vs_I-H-130719-N1-1-D1-2.annot.bedpe",sep="\t")
#svv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130720-T1-2-D1-2_vs_I-H-130720-N1-1-D1-2.annot.bedpe",sep="\t")
#svv=read.table("/Users/yellapav/Desktop/p220_wgs/misc/oncoplot/SV/I-H-130718-T1-6-D1-2_vs_I-H-130718-N1-1-D1-2.annot.bedpe",sep="\t")
#gene<- read.delim("ptcl_GENE_LIST.txt",sep="\t", stringsAsFactors = F)



colnames(svv)=c("chr1","start1","end1","chr2","start2","end2","id/name","brass_score","strand1","strand2","sample","svclass","bkdist","assembly_score","readpair names","readpair count","bal_trans","inv","occL","occH","copynumber_flag","range_blat","Brass Notation","non-template","micro-homology","assembled readnames","assembled read count","gene1","gene_id1","transcript_id1","strand1","end_phase1","region1","region_number1","total_region_count1","first/last1","gene2","gene_id2","transcript_id2","strand2","phase2","region2","region_number2","total_region_count2","first/last2","fusion_flag")
svv=svv[svv$sample!="I-H-130718-T1-6-D1-2,I-H-130718-N1-1-D1-2",]

svv$sample="I-H-130719-T1-4-D1-2"
svv=svv[svv$sample!="I-H-130720-T1-2-D1-2,I-H-130720-N1-1-D1-2",]


library(plotrix)
library(igraph)
library(quantsmooth)
###sample_list<- intersect(unique(sv$sample) , unique(cnv$IDA))

samp_regex="I-H"
for(bbg in grep(samp_regex,list.files("/Users/yellapav/Desktop/p220_2019/genome_plots/input",".*subclones.txt$",full.names=T), value=T)){
  brass = gsub("_subclones.txt","_brass.txt",bbg)
  
  svv=read.table(brass,sep="\t",header = T)
  colnames(svv)=tolower(colnames(svv))
  #colnames(svv)=c("chr1","start1","end1","chr2","start2","end2","id/name","brass_score","strand1","strand2","sample","svclass","bkdist","assembly_score","readpair names","readpair count","bal_trans","inv","occL","occH","copynumber_flag","range_blat","Brass Notation","non-template","micro-homology","assembled readnames","assembled read count","gene1","gene_id1","transcript_id1","strand1","end_phase1","region1","region_number1","total_region_count1","first/last1","gene2","gene_id2","transcript_id2","strand2","phase2","region2","region_number2","total_region_count2","first/last2","fusion_flag")
  svv$sample=substr(svv$sample,0,20)
  
  cnv=read.table(bbg,sep="\t",header=TRUE)
  cnv=cnv[,c(1:14)]
  print(c(bbg,brass))



chrom<- unique(cnv$chr)


svv$pos1<- as.numeric(as.character(svv$start1))
svv$pos2<- as.numeric(as.character(svv$start2))


###sample_list<- c("MMRF_1700_1_BM","MMRF_1700_2_BM" )
###setwd("/Volumes/GoogleDrive/My Drive/MSKCC/GRANT/SCOR/figure/new_jan/figure_2/cn_plots_fig_2/")
setwd("/Users/yellapav/Desktop/p220_2019/genome_plots/plots")
sample_list=unique(svv$sample)
for(w in (1:length(unique(svv$sample))))
{
  
  pdf(sprintf("%s_CN_COMMPASS_plots.pdf", sample_list[w]), height=4, width=12)
  ###ascat_sample2<- cnv[cnv$IDA == sample_list[w],]
  ###ascat_sample<- ascat_sample2[,c(1:6)]
  ###colnames(ascat_sample)<- c("sample","chrom","start","end","tot","min")
  cnv$diff<- cnv$endpos - cnv$startpos
  cnv<- cnv[order(cnv$chr, cnv$startpos),]
  kk<-sample_list[w]
  brass_sample2<- svv[svv$sample == kk,]
  brass_sample<- brass_sample2[,c("sample","chr1","start1","chr2","start2","svclass","strand1","strand2")]
  brass_sample$chr1=as.character(brass_sample$chr1)
  brass_sample$chr2=as.character(brass_sample$chr2)
  
  for(jj in (1:23))
  {  
    i<- as.character(chrom[jj])
    
    #gene_chr<- gene[gene$Chromosome.Name==i,]
    brass_sample_chr <- brass_sample[brass_sample$chr1 == i | brass_sample$chr2 ==i,]
    cnv_filter <- cnv[cnv$chr==i,]
    if(nrow(cnv_filter)>0){
      cnv_filter_chr<- cnv_filter
      cnv_filter_chr$ntot[is.na(cnv_filter_chr$ntot)]<-0
      #cnv_filter_chr$min[is.na(ascat_sample_filter_chr$min)]<-0
    iArrows <- igraph:::igraph.Arrows
    
    ### setwd("/Volumes/fm6/chapman_data/copy_number_plots/all_copy_number_plot_batt//")
    #header
    name<-sample_list[w]
    par(mar=c(5,7,2,5), xpd=F)
    par(mfrow=c(1,1))
    plot(1, type="n", xlab="", xlim=c(0, max(cnv_filter_chr$endpos)), 
         ylim=c(0, max(cnv_filter_chr$ntot)+3), 
         ylab="", las=2,bty ="n",
         xaxt="n", yaxt="n", main=paste(kk, "-","Chromosome", i, sep=" "))
    
    #axis and cytoband
    axis(2, line=0, at = c(0:(max(cnv_filter_chr$ntot) +2)),  col.axis = "black",
         labels=c(0:(max(cnv_filter_chr$ntot) +2)) ,lwd.ticks =1, cex.axis = 1, las=2)
    mtext("Copy Number", 2, +2.5, at = median(c(0:(max(cnv_filter_chr$ntot) +1))), cex=1)
    par(new=TRUE, xpd=T)
    paintCytobands(i, pos = c(0, -0.4), units = "bases", width = 0.3, cex.leg = 0.1,
                   bands = "major", legend=F, length.out =max(cnv_filter_chr$endpos)) 
    
    #grey bar on cytoband
    rect(par("usr")[1],par("usr")[3],par("usr")[2],0.5,col = "gray96", border = FALSE, lwd = 0.001)
    
    col<- rep(c("gray90","gray96"), max(cnv_filter_chr$ntot)+2)
    
    #paint alternate grey bars
    for(z in (1:(max(cnv_filter_chr$ntot)+2)))
    {
      
      rect(par("usr")[1],z-0.5,par("usr")[2],z+0.5,col = col[z], border = FALSE, lwd = 0.001)
    }
    
    #draw line for total and min
    segments(cnv_filter_chr$startpos, cnv_filter_chr$ntot, cnv_filter_chr$endpos,cnv_filter_chr$ntot,
             lwd=2.5, col="blue3")
    #segments(ascat_sample_filter_chr$start, ascat_sample_filter_chr$min, ascat_sample_filter_chr$end,ascat_sample_filter_chr$min,
    #         lwd=2.5, col="darkgoldenrod1", lty = 2)
    
    if(nrow(brass_sample_chr)>0)
    {  
      for(j in (1:nrow(brass_sample_chr)))
      {
        if(brass_sample_chr$svclass[j] == "deletion")
        {
          #drawline at breakpoint1 and 2
          segments(brass_sample_chr$start1[j], 0, brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1, col="dodgerblue", lwd=1)
          segments(brass_sample_chr$start2[j], 0, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1, col="dodgerblue", lwd=1)
          grid.curve(brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1)
          #connect the two lines
          iArrows(brass_sample_chr$start1[j],max(cnv_filter_chr$ntot)+1, brass_sample_chr$start2[j],
                  max(cnv_filter_chr$ntot)+1,h.lwd=1, sh.lwd=1, sh.col="dodgerblue", curve=0.00000005 , width=1, size=0.002)
        }else
        {
          if(brass_sample_chr$svclass[j] == "tandem-duplication")
          {
            segments(brass_sample_chr$start1[j], 0, brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1.5, col="firebrick2", lwd=1)
            segments(brass_sample_chr$start2[j], 0, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1.5, col="firebrick2", lwd=1)
            grid.curve(brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1.5, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1.5)
            iArrows(brass_sample_chr$start1[j],max(cnv_filter_chr$ntot)+1.5, brass_sample_chr$start2[j], 
                    max(cnv_filter_chr$ntot)+1.5,h.lwd=1, sh.lwd=1, sh.col="firebrick2", curve=0.00000005 , width=1, size=0.002)
          }else
          {
            
            if(brass_sample_chr$svclass[j] == "inversion")
            {
              segments(brass_sample_chr$start1[j], 0, brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+2, col="green", lwd=1)
              segments(brass_sample_chr$start2[j], 0, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+2, col="green", lwd=1)
              grid.curve(brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+2, brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+2)
              iArrows(brass_sample_chr$start1[j],max(cnv_filter_chr$ntot)+2, brass_sample_chr$start2[j], 
                      max(cnv_filter_chr$ntot)+2,h.lwd=1, sh.lwd=1, sh.col="green", curve=0.00000005 , width=1, size=0.002)
            }
            else
            {
              if(brass_sample_chr$chr1[j] ==i)
              {
                segments(brass_sample_chr$start1[j], 0, brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1, col="black", lwd=1)
                if(brass_sample_chr$strand1[j] == "+")
                {
                  segments(brass_sample_chr$start1[j],  max(cnv_filter_chr$ntot)+1, brass_sample_chr$start1[j]+max(cnv_filter_chr$endpos)/100, 
                           max(cnv_filter_chr$ntot)+2, col="black", lwd=1)
                  text(brass_sample_chr$start1[j]+max(cnv_filter_chr$endpos)/100, max(cnv_filter_chr$ntot)+sample(c(2,2.2,2.3,2.5,2.65,2.8,2.95),1), brass_sample_chr$chr2[j], cex=0.8)
                }else
                {
                  segments(brass_sample_chr$start1[j], max(cnv_filter_chr$ntot)+1, brass_sample_chr$start1[j]-max(cnv_filter_chr$endpos)/100, 
                           max(cnv_filter_chr$ntot)+2, col="black", lwd=1)
                  text(brass_sample_chr$start1[j]-max(cnv_filter_chr$endpos)/100, max(cnv_filter_chr$ntot)+sample(c(2,2.2,2.3,2.5,2.65,2.8,2.95),1), brass_sample_chr$chr2[j], cex=0.8)
                }}else
                {
                  segments(brass_sample_chr$start2[j], 0, brass_sample_chr$start2[j], 
                           max(cnv_filter_chr$ntot)+1, col="black", lwd=1)
                  if(brass_sample_chr$strand1[j] == "+")
                  {
                    segments(brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1, brass_sample_chr$start2[j]+max(cnv_filter_chr$endpos)/100, 
                             max(cnv_filter_chr$ntot)+2, col="black", lwd=1)
                    text(brass_sample_chr$start2[j]+max(cnv_filter_chr$endpos)/100, max(cnv_filter_chr$ntot)+sample(c(2.2,2.35,2.5,2.65,2.8,2.95,3.15),1),brass_sample_chr$chr1[j], cex=0.8)
                  }else
                  {
                    segments(brass_sample_chr$start2[j], max(cnv_filter_chr$ntot)+1, brass_sample_chr$start2[j]-max(cnv_filter_chr$endpos)/100, 
                             max(cnv_filter_chr$ntot)+2, col="black", lwd=1)
                    text(brass_sample_chr$start2[j]-max(cnv_filter_chr$endpos)/100, max(cnv_filter_chr$ntot)+sample(c(2.2,2.39,2.5,2.65,2.8,2.95,3.15),1), brass_sample_chr$chr1[j], cex=0.8)
                  }}
            }
          }
        }
      }}

    }}
  dev.off()
}
}


#graphics.off()
#plot(1,1)
head(svv)
