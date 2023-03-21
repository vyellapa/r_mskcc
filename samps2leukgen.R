#Read tables staudt louis
mm=read.table("~/Desktop/MM_Sham/MM_forLeukform.txt",sep="\t",header=T,stringsAsFactors = FALSE)
samps=read.table("~/Desktop/MM_Sham/samples.tsv",sep=" ",header=F,stringsAsFactors = FALSE)
mm=mm[,c(1:9)]

mm$Institution=rep("MEMORIAL SLOAN KETTERING (MSK)",nrow(mm))
mm$Gender=rep("UNKNOWN",nrow(mm))
mm$specimen=mm$tiessuetype

mm$specimen[mm$specimen=="tumor"]="TUMOR"
mm$specimen[mm$typeofpatient=="Normal"]="NORMAL"
mm$specimen[mm$typeofpatient=="normal"]="NORMAL"
mm$diseaseAtCollection=mm$typeofpatient

mm$diseaseAtCollection[mm$diseaseAtCollection=="Normal"]="NOT APPLY (NA)"
mm$diseaseAtCollection[mm$diseaseAtCollection=="normal"]="NOT APPLY (NA)"
mm$diseaseAtCollection[mm$diseaseAtCollection=="Myeloma"]="MULTIPLE MYELOMA (MM)"
mm$diseaseAtCollection[mm$diseaseAtCollection=="myeloma"]="MULTIPLE MYELOMA (MM)"
mm$diseaseAtCollection[mm$diseaseAtCollection=="mgus"]="MONOCLONAL GAMMOPATHY OF UNDETERMINED SIGNIFICANCE (MGUS)"
mm$diseaseAtCollection[mm$diseaseAtCollection=="Smoldering"]="SMOLDERING MYELOMA (SMM)"
mm$extractionAnalyte=rep("RNA",nrow(mm))
mm$seqCenter=rep("OTHER",nrow(mm))
mm$seqTechnique=rep("RNA|WTA",nrow(mm))
mm$seqPlatform=rep("ILLUMINA-HISEQ-2500",nrow(mm))
mm$seqReadLength=rep("101",nrow(mm))
mm$seqReadType=rep("PAIR-END",nrow(mm))
mm$diseaseAtCollection[mm$diseaseAtCollection=="MGUS"]="MONOCLONAL GAMMOPATHY OF UNDETERMINED SIGNIFICANCE (MGUS)"
nrow(mm)
mm$SID=mm$Sample.ID..original.
mm$SID[mm$Sample.ID..original.==""]=mm$Sample.ID[mm$Sample.ID..original.==""]
dim(mm)
mm=mm[,c(1,2,3,4,5,20,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
mm$SID[mm$diseaseAtCollection=="cell lines"]=mm$Patient.ID[mm$diseaseAtCollection=="cell lines"]

mm=mm[mm$SID!="",]
mm$SID=gsub(" 100ng", "_100ng", mm$SID)
mm$SID=gsub(" 300ng", "_300ng", mm$SID)
mm$SID=gsub("NHMM ", "NHMM", mm$SID)
mm$SID=gsub("NORMAL ", "NORMAL_", mm$SID,perl=T)
mm$SID=gsub("BM ", "BM", mm$SID,perl=T)
mm$SID=gsub("normal ", "normal_", mm$SID,perl=T)


mm$SID=gsub('CD138[+] ', "CD138plus_", mm$SID)
mm$SID=gsub('CD138[+]', "CD138plus_", mm$SID)
mm$SID=gsub('CD138[-] ', "CD138_", mm$SID)
mm$SID=gsub('NHMM135 CD138[-]', "NHMM135", mm$SID)
mm$SID=gsub('CD138-NHMM33', "CD138_NHMM33", mm$SID)
mm$SID=gsub(' ', "", mm$SID)

mm$SID2=mm$SID


merged=merge(samps, mm, by.x="V1", by.y="SID", all=T)
merged=merged[,c(1,22,2:21)]

#merged$specimen[merged$typeofpatient=="cell lines"]="CELL_LINES"
merged$diseaseAtCollection[merged$typeofpatient=="cell lines"]="MULTIPLE MYELOMA (MM)"
merged$specimen[merged$specimen=="non-tumor"]="TUMOR"


write.table(merged,file="~/Desktop/MM_Sham/leukgen_form.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote=FALSE)

