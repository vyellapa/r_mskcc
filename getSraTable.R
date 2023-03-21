#dbDownload only possible on MKSGuest

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

write.table(rs,file="~/Desktop/AML_Leucegene/SRAdb_annotations.tsv", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
write.table(rsSRR,file="~/Desktop/AML_Leucegene/SRAdb_annotations_SRR.tsv", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
write.table(amlRs,file="~/Desktop/AML_Leucegene/SRAdb_AML.tsv", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)

#sqliteQuickSQL(sra_con, 'PRAGMA TABLE_INFO(study)')



#a <- dbGetQuery(sra_con,"select * from col_desc")
#b <- dbGetQuery(sra_con,"select * from experiment")
#c <- dbGetQuery(sra_con,"select * from metaInfo")
#d <- dbGetQuery(sra_con,"select * from run")
#e <- dbGetQuery(sra_con,"select * from sample")
#f <- dbGetQuery(sra_con,"select * from sra_ft")
#g <- dbGetQuery(sra_con,"select * from sra_ft_content")
#h <- dbGetQuery(sra_con,"select * from sra_ft_segdir")
#i <- dbGetQuery(sra_con,"select * from sra_ft_segments")
#j <- dbGetQuery(sra_con,"select * from study")
