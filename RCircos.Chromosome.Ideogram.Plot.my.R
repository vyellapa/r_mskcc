RCircos.Chromosome.Ideogram.Plot.my <- function () 
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    right.side <- nrow(RCircos.Pos)/2
    outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
    inner.location <- RCircos.Par$chr.ideog.pos
    chroms <- unique(RCircos.Cyto$Chromosome)
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
            ]
        start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        end <- the.chr$Location[nrow(the.chr)]
        mid <- round((end - start + 1)/2, digits = 0) + start
        chr.color <- 'grey'
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
            RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
            RCircos.Pos[end:start, 2] * inner.location)
        polygon(pos.x, pos.y, border='grey', lwd=0.5)
        chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
        text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
            RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name, 
            srt = RCircos.Pos$degree[mid], col='grey', cex=0.6)
        lines(RCircos.Pos[start:end, ] * RCircos.Par$highlight.pos, 
            col = chr.color, lwd = 0.5)
    }
    for (a.band in 1:nrow(RCircos.Cyto)) {
        a.color <- RCircos.Cyto$BandColor[a.band]
        if (a.color == "white") {
            next
        }
        start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 
            1
        end <- RCircos.Cyto$Location[a.band]
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
            RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
            RCircos.Pos[end:start, 2] * inner.location)
        polygon(pos.x, pos.y, col = alpha(a.color,0.25), border = NA)
    }
}
