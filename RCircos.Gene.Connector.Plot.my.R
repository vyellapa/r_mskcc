RCircos.Gene.Connector.Plot.my <- function (genomic.data, track.num, side, in.pos = 1.32) 
{
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    gene.data <- RCircos.Get.Plot.Data(genomic.data, "plot")
    label.data <- RCircos.Get.Gene.Label.Locations(gene.data)
    connect.data <- data.frame(label.data$Location, label.data$Label.Position)
    locations <- RCircos.Track.Positions(side, track.num)

    out.pos <- locations[1] # 


    line.colors <- RCircos.Get.Plot.Colors(label.data, RCircos.Par$text.color)
    
    
    genomic.col <- ncol(connect.data) - 1
    label.col <- ncol(connect.data)
    chroms <- unique(connect.data[, 1])
    for (a.chr in 1:length(chroms)) {
        chr.row <- which(connect.data[, 1] == chroms[a.chr])
        total <- length(chr.row)
        for (a.point in 1:total) {
            top.loc <- out.pos
            bot.loc <- RCircos.Par$track.in.start  - sum(RCircos.Par$track.heights[1:length(RCircos.Par$track.heights)]) - sum(RCircos.Par$track.padding[1:length(RCircos.Par$track.padding)] ) - 0.02
  
            
            p1 <- connect.data[chr.row[a.point], genomic.col]
            p2 <- connect.data[chr.row[a.point], genomic.col] # p2 <- connect.data[chr.row[a.point], label.col]

            # lines(c(RCircos.Pos[p1, 1] * out.pos, RCircos.Pos[p1,1] * top.loc),
            #       c(RCircos.Pos[p1, 2] * out.pos, RCircos.Pos[p1, 2] * top.loc), col = 'red')
            
            # lines(c(RCircos.Pos[p2, 1] * bot.loc, RCircos.Pos[p2, 1] * in.pos),
            #       c(RCircos.Pos[p2, 2] * bot.loc, RCircos.Pos[p2, 2] * in.pos), col = 'green')

            
            lines(c(RCircos.Pos[p1, 1] * top.loc, RCircos.Pos[p2, 1] * bot.loc), # xs
                  c(RCircos.Pos[p1, 2] * top.loc, RCircos.Pos[p2, 2] * bot.loc), col =  alpha('black', 0.1), lwd=0.5) # ys
        }
    }
}
