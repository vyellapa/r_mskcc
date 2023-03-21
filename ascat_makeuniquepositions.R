args=(commandArgs(trailingOnly=TRUE))
posall.file = args[1]
pos.out = args[2]


posall = read.table(posall.file,sep="\t")

ctrans = 1:23
names(ctrans) = c(1:22,"X")

posg = ctrans[as.vector(posall[,1])]*1000000000+as.vector(posall[,2])

posgu = unique(posg)

posgus = sort(posgu)

btrans = c(1:22,"X")
names(btrans)=c(1:23)

pos = cbind(btrans[floor(posgus/1000000000)],posgus-1000000000*floor(posgus/1000000000))

write.table(pos,pos.out,sep="\t",quote=F,col.names=F,row.names=F)
