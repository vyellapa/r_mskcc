args=(commandArgs(TRUE))
tumor_file = args[1]
normal_file = args[2]

tumor = read.delim(file=tumor_file, header=F, sep="\t", comment.char='#')
normal = read.delim(file=normal_file, header=F, sep="\t", comment.char='#')
tumor$id = paste(tumor$V1, tumor$V2)
normal$id = paste(normal$V1, normal$V2)
common.ids = intersect(tumor$id, normal$id)
tumor = tumor[tumor$id %in% common.ids,]
normal = normal[normal$id %in% common.ids,]

write.table(t(c('#CHR','POS','Count_A','Count_C','Count_G','Count_T','Good_depth')), file=tumor_file, row.names=F, col.names=F, quote=F, sep="\t")
write.table(t(c('#CHR','POS','Count_A','Count_C','Count_G','Count_T','Good_depth')), file=normal_file, row.names=F, col.names=F, quote=F, sep="\t")
write.table(tumor[,1:7], file=tumor_file, row.names=F, col.names=F, quote=F, sep="\t", append=T)
write.table(normal[,1:7], file=normal_file, row.names=F, col.names=F, quote=F, sep="\t", append=T)
