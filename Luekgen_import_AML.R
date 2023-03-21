#dbDownload only possible on MKSGuest
library(SRAdb)
if(!file.exists('SRAdb.sqlite')) {  sqlfile <- getSRAdbFile() }

file.info('SRAmetadb.sqlite')

sra_con <- dbConnect(SQLite(),'SRAmetadb.sqlite')

sra_con

sra_tables <- dbListTables(sra_con)
sra_tables

dbListFields(sra_con, 'sra')
rs <- dbGetQuery(sra_con,"select * from sra")
rsSRR=rs[grep("^SRR.*", rs$run_accession, perl = TRUE),]

srrID=read.table("~/Desktop/AML_Leucegene/AML_SRA.tsv")
amlRs=rsSRR[rsSRR$run_accession %in% srrID$V1,]

amlRs$LK_EXTERNAL_ID=paste(amlRs$run_accession, unlist(lapply(strsplit(amlRs$experiment_title, ";"), function(x) {gsub(": ","__",x[[1]][1])})), sep="__")
amlRs$LK_GENDER=rep("UNKNOWN",nrow(amlRs))

amlRs$LK_SPECIES=rep("HUMAN",nrow(amlRs))
amlRs$LK_INSTITUTION=rep("OTHER",nrow(amlRs))
amlRs$LK_SPECIMEN_TYPE=rep("TUMOR",nrow(amlRs))

> amlRs$LK_ANALYTE=rep("RNA",nrow(amlRs))
> amlRs$LK_SEQ_TECHNIQUE=rep("RNA|WTA",nrow(amlRs))
amlRs$LK_SEQ_PLATFORM=rep("ILLUMINA-HISEQ-2000",nrow(amlRs))
amlRs$LK_READ_TYPE=rep("PAIR-END",nrow(amlRs))
amlRs$LK_READ_LEN=rep("100",nrow(amlRs))
amlRs$LK_SEQ_CEN_ID=amlRs$run_accession
amlRs$LK_Specimen_Notes=amlRs$sample_attribute
amlRs[grep("myeloid",amlRs$sample_attribute), ]$LK_DISEASE="ACUTE MYELOID LEUKEMIA (AML)"
amlRs[grep("lymphoblastic",amlRs$sample_attribute), ]$LK_DISEASE="ACUTE LYMPHOBLASTIC LEUKEMIA (ALL)"
amlRs[grep("cord",amlRs$sample_attribute), ]$LK_DISEASE="NOT APPLY (NA)"

amlRs[grep("cord",amlRs$sample_attribute), ]$LK_SPECIMEN_TYPE="NORMAL"

write.table(amlRs,file="~/Desktop/AML_Leucegene/SRAdb_AML_lkGen.tsv", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)



