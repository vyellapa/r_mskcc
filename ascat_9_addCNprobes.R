SNPprobes = read.table("pos.snps.txt")

design = read.table("capturedesign.bed")
# these probes are all 120 bp

# join them into a smaller number of regions (if less than 100 bp apart):
ctrans = 1:24
names(ctrans) = c(1:22,"X","Y")
dpos = 1000000000 * ctrans[as.vector(design[,1])] + design[,2]
#diff = (c(dpos,24000000000)-c(0,dpos))[-1] There was a bug here. You need to take the absolute value of the distance!!!
diff = abs((c(dpos,24000000000)-c(0,dpos)))[-1]
join = diff<220 # 220 = 120 bp probes + 100 bp gap
startchr = as.vector(design[1,1])
startpos = design[1,2]
endpos = design[1,3]
joinedreg = NULL
for (i in 1:dim(design)[1]) {
  print(i)
  if(join[i]) {
    endpos = design[i+1,3]
  } 
  else {
    joinedreg = rbind(joinedreg,c(as.vector(startchr),startpos,endpos))
    if(i!=dim(design)[1]) {
      startchr = as.vector(design[i+1,1])
      startpos = design[i+1,2]
      endpos = design[i+1,3]
    }
  }
}

# this is just to check the distance between the SNPprobes:
#dpos = 1000000000 * ctrans[as.vector(SNPprobes[,1])] + SNPprobes[,2]
#diff = (c(dpos,24000000000)-c(0,dpos))[-1]
# 475 are <75bp - might need to change in next iteration

posSNP = 1000000000 * ctrans[as.vector(SNPprobes[,1])] + SNPprobes[,2]

finalpos = NULL
for (r in 1:dim(joinedreg)[1]) {
#for(r in 1:5642) {
  print(r)
  startpos = 1000000000 * ctrans[as.vector(joinedreg[r,1])] + as.numeric(joinedreg[r,2])
  endpos = 1000000000 * ctrans[as.vector(joinedreg[r,1])] + as.numeric(joinedreg[r,3])

  npos = startpos
  while(npos<=endpos) {
    if(sum(abs(posSNP-npos)<100)>0) {
      npos = posSNP[abs(posSNP-npos)<100][1]
    }
    else {
      # this is under an else condition because this doesn't capture all SNP probes, so we will add them later
      finalpos = c(finalpos,npos)
    } 
    npos = npos + 100
  }
}
finalpos=c(finalpos,posSNP)
finalpos=sort(finalpos)

btrans = c(1:22,"X")
names(btrans) = 1:23

finalloci = cbind(btrans[floor(finalpos/1000000000)],finalpos%%1000000000)

write.table(finalloci,"pos.both.txt", quote=F, col.names=F, row.names=F)
