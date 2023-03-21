args=(commandArgs(trailingOnly=TRUE))
snppos_file = args[1]
gccorrection_file = args[2]
gc_correction_dir = '/ifs/work/leukgen/home/gg10/soft/gc_correction/'


# this is done in a separate directory where the GC content files are stored
SNPpos<-read.table(file=snppos_file,header=TRUE,as.is=TRUE)
SNPpos = SNPpos[!is.na(SNPpos$chr),]

GCcorrect = NULL
for (chrindex in unique(SNPpos[,1])) {
  print(chrindex) 
  Position = SNPpos
  oldname<-names(Position)
  names(Position)<-c("Chr","Position")

  Position[,1]<-as.character(Position[,1])
  Position[,2]<-as.numeric(as.character(Position[,2]))
  Position<-Position[Position[,1]==chrindex,]
  temp<-Position

  base_tail<-substring(Position[,2],nchar(Position[,2]),nchar(Position[,2]))
  base_tail[as.numeric(base_tail)<6 & as.numeric(base_tail)>0]<-1
  base_tail[as.numeric(base_tail)>=6 | as.numeric(base_tail)==0]<-6
  posindex_want<-paste(substring(Position[,2],1,nchar(Position[,2])-1),base_tail,sep="")
  Position_tot<-cbind(Position,posindex_want)

  GCindexlist<-read.table(paste(gc_correction_dir, "out",chrindex,".txt",sep=""),skip=1)

  SNP_unique_index<-which(as.numeric(GCindexlist[,1])%in%as.numeric(unique(posindex_want)))
#  SNP_index<-cbind(Position=unique(posindex_want),Index=SNP_unique_index)
# adapted to make more robust...
  SNP_index<-cbind(Position=intersect(as.numeric(unique(posindex_want)),as.numeric(GCindexlist[,1])),Index=SNP_unique_index)
  SNP_index_merge<-merge(Position_tot,SNP_index,by.x="posindex_want",by.y="Position",all.x=TRUE)
  SNP_index_merge<-SNP_index_merge[,-1]
  SNP_index_merge$Position<-as.numeric(as.character(SNP_index_merge$Position))
  SNP_index_merge$Index<-as.numeric(as.character(SNP_index_merge$Index))

  SNP_index_merge<-SNP_index_merge[order(SNP_index_merge$Position),]

  for (i in c(0.000013,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064,0.0128,0.0256,0.0512,0.1024,0.2048,1,2,5,10))
  {
    cat("i=",i,"\n")
    temp_left<-(SNP_index_merge$Index-1000000*i/(2*5))
    temp_right<-(SNP_index_merge$Index+1000000*i/(2*5))
    start<-1
    end<-nrow(GCindexlist)
    temp_left[which(temp_left<0)]<-1
    temp_right[which(temp_right>end)]<-end
    #next lines added to solve NA problem.. -> if nothing around the position is found back, revert to the chromosomal average
    temp_left[which(is.na(temp_left))]<-1
    temp_right[which(is.na(temp_right))]<-end
    range<-cbind(temp_left,temp_right)
    GCcontent<-apply(range,1,function(x) mean(GCindexlist[(x[1]:x[2]),2]))
    Position<-cbind(Position,GCcontent)
  }
  colnames(Position)[3:ncol(Position)]<-c("Probe","200bp","400bp","800bp","1600bp","3200bp","6400bp","12800bp","25600bp","51200bp","102400bp","204800bp","1M","2M","5M","10M")
  
  GCcorrect = rbind(GCcorrect,Position)
}
write.table(GCcorrect,file=gccorrection_file, quote=F, row.names=TRUE, col.names=TRUE,sep="\t")
