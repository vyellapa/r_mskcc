Sys.setenv(GENCODE_DIR = "/Users/yellapav/Desktop/SV_paper/gtrack")
ge = track.gencode("/Users/yellapav/Desktop/SV_paper/gtrack/CDKN2C.gtf")

bedpe = jJ("/Users/yellapav/Desktop/SV_paper/ggenome/MMRF_1682_1_BM.bedpe")
gg = gG(juncs = bedpe)


seg = read.table("/Users/yellapav/Desktop/SV_paper/ggenome/MMRF_1682_1_BM.cnv",header = T,sep = "\t")
scna = seg2gr(seg)

amps = scna %Q% (seg.mean>-1 & num.mark > 50)
amp.score = as(coverage(amps), 'GRanges')
gt.amps = gTrack(amps,  col = 'blue', name = 'Amps')
gt.amp.score = gTrack(amp.score, y.field = 'score',col = 'red', name = 'Amp score', line = TRUE)


gt.amp.score = gTrack(amps, y.field = 'seg.mean',col = 'blue', name = 'Total CN', line = TRUE)


plot(c(ge,gt.amp.score,gg$gt),c("1:40736959-66820179","2:64428646-70411246","12:96650771-96990771","18:2503163-2803163","22:23133429-23473429"),links = bedpe$grl)

#rbind(r1,r2) %>% group_by(chrom1) %>% summarize(min = min(start1-10000), max = max(start1+10000)) %>% data.frame() %>% mutate(key = paste0(chrom1,":",min-10000,"-",max+10000)) %>% dplyr::select(key)
dev.copy2pdf(file='myfile.pdf',useDingbats=FALSE,family='sans')
