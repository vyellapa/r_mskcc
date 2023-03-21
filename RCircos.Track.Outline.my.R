RCircos.Track.Outline.my <- function (out.pos, in.pos, num.layers) 
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    subtrack.height <- (out.pos - in.pos)/num.layers
    chroms <- unique(RCircos.Cyto$Chromosome)
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
            ]
        start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        end <- the.chr$Location[nrow(the.chr)]
        polygon.x <- c(RCircos.Pos[start:end, 1] * out.pos, RCircos.Pos[end:start, 
            1] * in.pos)
        polygon.y <- c(RCircos.Pos[start:end, 2] * out.pos, RCircos.Pos[end:start, 
            2] * in.pos)
        polygon(polygon.x, polygon.y, col = NULL, lwd=0.3, border=RCircos.Par$grid.line.color)
        
        for (a.line in 1:(num.layers - 1)) {
            height <- out.pos - a.line * subtrack.height
            lines(RCircos.Pos[start:end, 1] * height, RCircos.Pos[start:end, 2] * height, col = RCircos.Par$grid.line.color, lwd=0.3)
        }
    }
}
