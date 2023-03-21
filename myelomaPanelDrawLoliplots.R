
#Library loading (Surprise!!)
options( java.parameters = "-Xmx4g" )
library(xlsx)
library(plotrix)
library(RMySQL)
library(GenomicRanges)
library(RColorBrewer)
library(trackViewer)
library(plyr)
setwd("/Users/yellapav/Desktop/lolli_plot")

#Read gene panel
genePanel = read.table("~/Desktop/lolli_plot/mm_genes",header=F,sep="\t", stringsAsFactors=FALSE)


#Function that returns a table with exon coordinates (w.o. UTR)
GetExonsFromName <- function(name,genome="hg19")
{
    con = dbConnect(MySQL(),user="genome",dbname=genome,host="genome-mysql.cse.ucsc.edu")
    gene = dbGetQuery(con,paste('select * from refGene where name2="',name,'"',sep=''))
    if(dim(gene)[1]>0)
    {
        exon.start = list()
        exon.end = list()
        refseq.list = list()
        chromosome.list = list()
        for(i in 1:dim(gene)[1])
        {
            chromosome = gene[i,3]
            cds.start = as.numeric(gene[i,7])
            cds.end = as.numeric(gene[i,8])
            gene.name = gene[i,13]
            gene.refseq = gene[i,2]
            
   
              exon.start.buffer = unlist(as.numeric(strsplit(gene[i,10],",",fixed=T)[[1]]))
              exon.end.buffer = unlist(as.numeric(strsplit(gene[i,11],",",fixed=T)[[1]]))

            exon.start.buffer = lapply(exon.start.buffer,function(start){return(max(cds.start,start))})
            exon.end.buffer = lapply(exon.end.buffer,function(start){return(max(cds.start,start))})
            exon.start.buffer = lapply(exon.start.buffer,function(end){return(min(cds.end,end))})
            exon.end.buffer = lapply(exon.end.buffer,function(end){return(min(cds.end,end))})
            exon.start = c(exon.start,exon.start.buffer)
            exon.end = c(exon.end,exon.end.buffer)
            chromosome.list = c(chromosome.list,rep(chromosome,length(exon.end.buffer)))
            refseq.list = c(refseq.list,rep(gene.refseq,length(exon.end.buffer)))
        }
        
        dbDisconnect(con)
        
        return(cbind(refseq.list,chromosome.list,exon.start,exon.end,name))	
    }
    else
    {
        dbDisconnect(con)
        cat("missing the following gene :",name,"\n")
        return(c("","","","",name))
    }
}


#Check gene size
ensemblExons = lapply(as.character(genePanel[,1]),GetExonsFromName)
geneSize = lapply(1:nrow(genePanel),function(iter){
#geneSize = lapply(1:2,function(iter){
  consolidatedIntervals = GRanges(seqnames=unlist(ensemblExons[[iter]][,2]),ranges=IRanges(start=as.numeric(ensemblExons[[iter]][,3]),end=as.numeric(ensemblExons[[iter]][,4])))
  buffer = as.matrix(coverage(consolidatedIntervals)[[1]])
  return(c(unique(unlist(ensemblExons[[iter]][,5])),length(which(buffer!=0))))
})
geneSize = do.call(rbind,geneSize)

#correct gene size in the gene panel
genePanel = merge(genePanel,geneSize,by=1,all.x=T,all.y=T)
write.table(genePanel,"genePanel.txt",row.names=F,sep="\t",quote=F)



#Protein domains
con = dbConnect(MySQL(),user="genome",dbname="hg19",host="genome-mysql.cse.ucsc.edu")
pfam = dbGetQuery(con,'select * from ucscGenePfam')
dbDisconnect(con)


#maf mutations

mafMuts = read.table("~/Desktop/MMRF/MMRF_Canonical_GLFilter_uniq.maf", sep = "\t", stringsAsFactors=FALSE, encoding= "utf-8", quote="", header=T)
mafMuts[,5] = paste("chr",MafMuts[,5],sep="")

mafMutsCoords = paste(MafMuts[,5],MafMuts[,6],sep=":")
mafMutsCoords = paste(mafMutsCoords,MafMuts[,7],sep="-")


drawMafPlot <- function(geneName,pfam,mafMuts,fileName){
  
  con = dbConnect(MySQL(),user="genome",dbname=genome,host="genome-mysql.cse.ucsc.edu")
  gene = dbGetQuery(con,paste('select * from refGene where name2="',geneName,'"',sep=''))
  dbDisconnect(con)

  chromosome = unique(gene[which(gene[,3]%in%c(paste("chr",1:22,sep=""),"chrX","chrY")),3])
  txStart = min(as.numeric(gene[,5]))
  txEnd = max(as.numeric(gene[,6]))
  geneRange = GRanges(chromosome,IRanges(start=txStart,end=txEnd))

  pfamBuffer = pfam[which(pfam[,2]==chromosome),]
  pfamRange = GRanges(chromosome,IRanges(start=pfamBuffer[,3],end=pfamBuffer[,4]))
  overlaps = findOverlaps(pfamRange,geneRange)
  overlaps =as.matrix(overlaps)
  pfamBuffer = pfamBuffer[overlaps[,1],]
  if(dim(pfamBuffer)[1]>0)
  {  
    pfamRange = GRanges(chromosome,IRanges(start=pfamBuffer[,3],end=pfamBuffer[,4],names=pfamBuffer[,5]))
    pfamColors = brewer.pal(n=length(unique(pfamBuffer[,5])),name="Dark2")
    pfamColors = lapply(pfamBuffer[,5],function(name){return(pfamColors[which(unique(pfamBuffer[,5])==name)])})
    pfamRange$fill = pfamColors
  }
  else
  {
    pfamRange = GRanges(NULL,IRanges(start=NULL,end=NULL))
  }
  
  mafPalette = brewer.pal(n=8,name="Set1")
  mafPalette = mafPalette[c(6,2,7,4,5,1,3,8)]
  mutBuffer = mafMuts[which(mafMuts[,"Chromosome"]==chromosome),]
  mafRange = GRanges(chromosome,IRanges(start=as.numeric(mutBuffer[,"Start_Position"]),end=as.numeric(mutBuffer[,"End_Position"])))
  overlaps = findOverlaps(mafRange,geneRange)
  overlaps =as.matrix(overlaps)
  mutBuffer = mutBuffer[overlaps[,1],]
  mutBuffer = mutBuffer[,which(colnames(mutBuffer)%in%c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele1","HGVSp_Short"))]
  mutBuffer = mutBuffer[mutBuffer$HGVSp_Short!="",]
  if(nrow(mutBuffer)<1) {print(paste(c("No non synonymous mutations in",geneName),collapse=""));break;}
  mutBuffer = count(mutBuffer,vars=colnames(mutBuffer))
  mutName = mutBuffer[,"HGVSp_Short"]
  colors = list(x=factor(mutBuffer[,"Variant_Classification"]))
  if(length(which(mutBuffer[,"Variant_Classification"]=="Intron"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Intron")] = mafPalette[1]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Missense_Mutation"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Missense_Mutation")] = mafPalette[2]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Silent"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Silent")] = mafPalette[3]
  if(length(which(mutBuffer[,"Variant_Classification"]=="3'UTR"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="3'UTR")] = mafPalette[4]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Splice_Site"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Splice_Site")] = mafPalette[5]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Nonstop_Mutation"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Nonstop_Mutation")] = mafPalette[6]
  if(length(which(mutBuffer[,"Variant_Classification"]=="5'UTR"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="5'UTR")] = mafPalette[4]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Nonsense_Mutation"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Nonsense_Mutation")] = mafPalette[6]
  if(length(which(mutBuffer[,"Variant_Classification"]=="5'Flank"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="5'Flank")] = mafPalette[1]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Translation_Start_Site"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Translation_Start_Site")] = mafPalette[5]
  if(length(which(mutBuffer[,"Variant_Classification"]=="In_Frame_Del"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="In_Frame_Del")] = mafPalette[7]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Frame_Shift_Del"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Frame_Shift_Del")] = mafPalette[8]
  if(length(which(mutBuffer[,"Variant_Classification"]=="Frame_Shift_Ins"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="Frame_Shift_Ins")] = mafPalette[8]
  if(length(which(mutBuffer[,"Variant_Classification"]=="In_Frame_Ins"))>0) colors[which(mutBuffer[,"Variant_Classification"]=="In_Frame_Ins")] = mafPalette[7]
  scores = mutBuffer[,"freq"]
  if(dim(mutBuffer)[1]>0)
  {
    mafRange = GRanges(chromosome,IRanges(start=as.numeric(mutBuffer[,"Start_Position"]),end=as.numeric(mutBuffer[,"Start_Position"]),names=mutName),color=unlist(colors),score=scores)
  }
  else
  {
    mafRange = GRanges(NULL,IRanges(start=NULL,end=NULL))
  }
  
  
  
  legend <- mafPalette
  names(legend) <- c("Intron", "Missense","Silent","UTR","Splice","Stop","In frame indel","Frameshift")
  
  if(length(pfamRange)==0){pfamRange = geneRange}

    pdf(fileName,width=15,height=10)
  if(max(scores)>10){
    lolliplot(GRangesList(Count=mafRange), pfamRange,ranges=geneRange,xaxis=F,legend=legend, cex=0.5) 
    #dandelion.plot(GRangesList(Count=mafRange), pfamRange,ranges=geneRange,xaxis=F,legend=legend) 
    }     
  else{
    print("Not available")
   #lolliplot(GRangesList(Count=mafRange), pfamRange,ranges=geneRange,xaxis=F,yaxis=F,legend=legend, cex=0.5) 
    #dandelion.plot(GRangesList(Count=mafRange), pfamRange,ranges=geneRange,xaxis=F,yaxis=F,legend=legend) 
  }
  dev.off()

}



Names = c(unlist(genePanel[,1]))


#geneNames = geneNames[which(!is.na(geneNames))]
#geneNames = geneNames[which(!geneNames%in%c("Size: ","Include ","Incl+maybe","Maybe","Lipkin genes","Immunotargets"))]

lapply(Names,function(gene){
  cat(gene,"\n")
  drawMafPlot(gene,pfam,MafMuts,paste(c("figs/",gene,".pdf"),collapse=""))
})

