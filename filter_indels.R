pon = read.table("/Users/yellapav/Desktop/record/MGUS_BRIAN/results/indels/normals/indel_pon.txt",sep = "\t",header = F) 
colnames(pon) = c("Freq","key")

for(i in list.files("~/Desktop/record/MGUS_BRIAN/results/indels/", pattern = "*.tsv.gz",full.names = T)){
  print(i)
  name = gsub(".indels.output.annot.tsv.gz","",basename(i))
  filter_file=sprintf("%s/%s.pindel_filter.tsv",dirname(i),name)
  filter_file1=sprintf("%s/%s.pindel_filter1.tsv",dirname(i),name)
  filter_file2=sprintf("%s/%s.pindel_filter2.tsv",dirname(i),name)
  filter_file_all=sprintf("%s/%s.pindel_filterAll.tsv",dirname(i),name)
  
  
  file = read.table(i, sep="\t",quote='~',header=T)
  file = file %>% mutate(key=paste(CHR,START,END,REF,ALT,sep=":")) %>% left_join(pon,by = c("key" = "key")) %>% 
    mutate(Freq = ifelse(is.na(Freq), 0, Freq))
  
  l = nchar(as.character(file$CONTEXT_5[1]))
  file$homo3=as.character(str_detect(substr(file$CONTEXT_3,1,6),"TTTTTT|AAAAAA|GGGGGG|CCCCCC"))
  file$homo5=as.character(str_detect(substr(file$CONTEXT_5,l-6,l),"TTTTTT|AAAAAA|GGGGGG|CCCCCC"))
  
  file = file %>% mutate(filter1 = if_else(Freq > 0,"Yes","No"), filter2=case_when(homo3=="TRUE" ~ 'Yes',
                                                                                   homo5=="TRUE" ~ "Yes",
                                                                                   TRUE ~ "No"))
  
  file = file[grep("P",file$PASSED_BY),] %>% dplyr::filter(pindel_READS_FORWARD>0 & pindel_READS_REVERSE>0)
  file1 = file %>% dplyr::filter(filter1=="No")
  file2 = file %>% dplyr::filter(filter2=="No")
  file3 = file %>% dplyr::filter(filter2=="No" & filter1=="No")
  
  write.table(file,file=filter_file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
  write.table(file1,file=filter_file1, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
  write.table(file2,file=filter_file2, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
  write.table(file3,file=filter_file_all, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
  
  
}

