bedpe = jJ("/Users/yellapav/Desktop/SV_paper/ggenome/MMRF_1534_1_BM.bedpe")
gg = gG(juncs = bedpe)
plot(gg$gt, links = bedpe$grl)

r= read.table("/Users/yellapav/Desktop/SV_paper/ggenome/MMRF_1534_1_BM.bedpe", header= T , sep="\t")
r1 = dplyr::select(r,chrom1,start1)
r1$chrom1 = as.character(r1$chrom1)
r1$start1 = as.numeric(as.character(r1$start1))

r2 = dplyr::select(r,chrom2,start2) %>%  dplyr::rename(chrom1 = chrom2,start1 = start2)
r2$chrom1 = as.character(r2$chrom1)
r2$start1 = as.numeric(as.character(r2$start1))


rbind(r1,r2) %>% group_by(chrom1) %>% summarize(min = min(start1-10000), max = max(start1+10000)) %>% data.frame()
   chrom1       min       max
#1       1 117460620 206848496
#2      10  70942829  70963149
#3      11  67083598 114709844
#4      12 123523726 123543726
#5      13  48359648  56403171
#6      14 101337107 106174384
#7      15  67463331  94029624
#8      16  81015462  89797615
#9      17   4365556  75750521
#10     19   1694226  17051022
#11      2  88581777 172863614
#12     20   2437482  62835296
#13     21  34444649  35730703
#14     22  18032122  47818192
#15      5 139528711 157663600
#16      6  29794062 166894295
#17      7  73619125 127735201
#18      8  27571426 145889142
#19      9 125381605 126810349
#20      X  17625129 153308027


plot(gg$gt,c("1:117460620-206848496","2:88581777-172863614","5:139528711-157663600","6:29794062-166894295","7:73619125-127735201","8:27571426-145889142","9:125381605-126810349", "10:70942829-70963149","11:67083598-114709844","12:123523726-123543726","13:48359648-56403171","14:101337107-106174384","15:67463331-94029624","16:81015462-89797615","17:4365556-75750521","19:1694226-17051022","20:2437482-62835296","21:34444649-35730703","22:18032122-47818192"),links = bedpe$grl)


gencode = track.gencode(stack.gap = 1e5, cex.label = 0.8, height = 20)


