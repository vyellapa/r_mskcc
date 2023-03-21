setwd("/Users/yellapav/Desktop/p220_2019/weinhold/facets")
l=list()
files=list.files(path="/Users/yellapav/Desktop/p220_2019/weinhold/facets", pattern = "I.*_profile.txt$")
for(i in seq_len(length(files)))
{
  print(files[i])
  j=read.table(files[i],header=T,sep="\t")
  j$sample=gsub("_profile.txt","",files[i])
  l[[i]]=j
  
}

all=do.call("rbind",l)