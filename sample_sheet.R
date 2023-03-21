library(dplyr)
library(stringr)

map = read.table("~/Desktop/record/MGUS_BRIAN/EGAD00001001901/delimited_maps/Run_Sample_meta_info.map", sep =";", header=F)
map$V1 = gsub("gender=|phenotype=|subject_id=","",map$V1)
map$V2 = gsub("gender=|phenotype=|subject_id=","",map$V2)
map$V3 = gsub("gender=|phenotype=|subject_id=","",map$V3)
map = map %>% dplyr::select(-c(V4,V5)) %>% 
  dplyr::mutate(key = ifelse(V2=="Normal", paste0(V3,'-PB'),paste0(V3,'-BM')))

sample = read.table("~/Desktop/record/MGUS_BRIAN/EGAD00001001901/delimited_maps/Sample_File.map", sep ="\t", header=F) %>% 
  dplyr::select(-c(V2,V3))
sample$S = substr(sample$V1,12,16)



md5 = read.table("~/Desktop/record/MGUS_BRIAN/EGAD00001001901/md5_check.txt", sep =" ", header=F) %>% dplyr::select(-c(V1,V2))
temp = as.data.frame(matrix(unlist(strsplit(as.character(md5$V3),split = '/')),ncol = 2,byrow=T)) %>% dplyr::rename("EGAF" = "V1")
md5 = cbind(md5,temp) 

map1 = dplyr::left_join(md5,sample, by = c("EGAF" = "V4")) %>% dplyr::left_join(map, by = c("S" = "key"))
head(map1)
write.table(map1, "~/Desktop/record/MGUS_BRIAN/sample_sheet.txt",sep = "\t", eol = "\n", quote = F, col.names = F, row.names = F, append = F)