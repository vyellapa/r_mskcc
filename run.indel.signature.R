setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/indel_signatures/vcfs/")

for(i in list.files(pattern=".*.txt")) {
  name=gsub(".txt", "",i)
  out.file=sprintf('%s.out.tsv',name)
  print(i)
    INDELS <- read.table(i, sep = "\t", header = T, stringsAsFactors = F) 
    out <- indelsClassification(mat = INDELS)
    table=as.data.frame(out[[1]])
    
    
    write.table(table, file=out.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    png(file = sprintf('%s.png',name), res=300, width=1200, height=1200)
    replayPlot(out[[2]])
    dev.off()
    
    }
