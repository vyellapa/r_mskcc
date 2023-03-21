# edited species and amb names

rm(list=ls()) #clear all vars 
#library(DESeq)
library(DESeq2)
subtitle="Monica"

# == define load/save PATH 
prj_name="Haiming NHD13 vs WT (Project_4100)"

countPath="/Volumes/armstronglab/Noushin/Haiming/Project_4100/htseq_output_mm9"
prjPath="/Volumes/armstronglab/Noushin/Haiming/Project_4100/"

# == Default params
cond1="LSD"
cond2="Cntrl"
species="mouse"
MyLibType="paired-end";

count_suffix="intersection-nonempty.strand_no.a10.mm9"
#=================================================================================
#=================================================================================
#=================================================================================
dir.create(file.path(prjPath, "Results"), showWarnings = FALSE)
dir.create(file.path(prjPath, "RData"), showWarnings = FALSE)

savePath=file.path(prjPath, "Results");

dir.create(file.path(savePath, "Lists"), showWarnings = FALSE)
dir.create(file.path(savePath, "Reports"), showWarnings = FALSE)
dir.create(file.path(savePath, "Figs"), showWarnings = FALSE)

#=================================================================================
# == Load the htseq-count files
#=================================================================================
list.files(countPath)
myfiles <- list.files(countPath, pattern= '*.count'); #read all .count files 

#========================================================== 
# Load the Lookup table to specify Nan's samples per cond
#==========================================================
library("stringr")
lookup<- read.table(file.path(prjPath,"samples_info.csv"),header=1,sep=",") # row.names=1 means use the 1st col as row.names; 
myconds <- colnames(lookup)
lookup[,1]=paste(str_trim(lookup[,1]),".",count_suffix,".count",sep="");
lookup[,2]=paste(str_trim(lookup[,2]),".",count_suffix,".count",sep="");
# here I was not able to input .xls files so had to open the ".csv" file in vim and perform the following :%s/^M/\r/g where ^M is in fact (CTRL+V+M) to substitute all the ^M's in the file by line break
#==========================================================
# Build the count-read matrix
#==========================================================
condition=character();

for (i in 1:length(myfiles)) {
  # print(myfiles[i])
  A<- read.table(file.path(countPath,myfiles[i]), col.names=c("gene",myfiles[i])) # row.names=1 means use the 1st col as row.names; 
  if (exists("countTable")) {countTable <- merge(countTable,A, by= "gene")} else {countTable <- A}
  
  if (!(is.na(any(pmatch(lookup[,1],myfiles[i]))))) {
    condition=c(condition,myconds[1])} else {
      if (!(is.na(any(pmatch(lookup[,2],myfiles[i]))))) {
        condition=c(condition,myconds[2])} else {
          stop('there is somthing wrong with finding your sample in the loopup file!')     
        }   
    }
  print(c(myfiles[i],condition[i]))
  rm(A)
}

# == create .gct file for BEFORE normalization
gct_before_normal=countTable;
gct_before_normal$id <- countTable$gene;
gct_before_normal$gene <- NULL;
 

# == Add Gene Names 
source("/Volumes/armstronglab/Noushin/R-scripts/find_gene_names.R")
new_countTable = find_gene_names(gct_before_normal,species); # here the species is pre-set to "mouse"

# == save .GCT file
write.csv(new_countTable, file.path(savePath,"Lists",paste("GCT-before-normalization-",subtitle,".csv",sep="")))
#cat("The RAW (before normalization) .gct file has been saved at ",file.path(savePath,paste("GCT-before-normalization.csv")),"\n")
print(paste("The raw reads BEFORE normalization are saved at: ",file.path(savePath,"Data",paste("GCT-before-normalization-",subtitle,".csv",sep=""))))
row.names(countTable)= countTable$gene;
countTable$gene <- NULL;

remove <- c("__alignment_not_unique","__ambiguous","__no_feature","__not_aligned","__too_low_aQual");
notMapped <- countTable[which(rownames(countTable) %in% remove), ]
countTable <- countTable[-which(rownames(countTable) %in% remove), ]

#countTable <- countTable[6:nrow(countTable),]

#==========================================  
# === save data for next script

myDesign= data.frame(
  row.names=colnames(countTable),
  condition=condition,
  libType= rep(c(MyLibType),ncol(countTable))) 

save(countTable, myDesign, condition, file=file.path(prjPath,"RData","countTable.RData"))

#==========================================  
jpeg(file.path(savePath,"Figs",paste0("Fig-01-BarPlot_HTSeq_notMapped_",subtitle,".jpg")), width=1480, height=980)
#max(countTable[1:4,1])
yrange <- c(0,20e6)
#yrange <- c(0,8000)
x <- notMapped[1:4,1]
ticks <- pretty(yrange,5)
ylabels <- format(ticks, big.mark=",", scientific=FALSE)

mp <- barplot(x, ylim=yrange, axes=0)
#axis(2, at = ticks, labels = ylabels, las = 1, cex.axis=0.7) 
axis(1, at=c(0,mp), labels=c("",row.names(notMapped[1:4,])))
axis(2, ticks, ylabels) 
abline(h=seq(5e6,20e6, by=5e6), col="gray", lty="dotted")
title(main=prj_name, col.sub="blue", cex.lab=0.75)
dev.off()


print(paste("FINISHED! The countTable data is saved at ",file.path(savePath,"RData","countTable.RData")))
print(paste("The HTSeq report file is saved at: ",file.path(savePath,"Figs",paste0("BarPlot_HTSeq_notMapped_",subtitle,".jpg"))))

