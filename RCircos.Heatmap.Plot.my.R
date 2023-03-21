RCircos.Heatmap.Plot.my <- function (heatmap.data, data.col, track.num, side, plotTrack=TRUE, heatmap.ranges=NA, heatmap.color=NA) 
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()

    min.with <- 1000000
    heatmap.data$width <- heatmap.data$chromEnd - heatmap.data$chromStart
    heatmap.data <- heatmap.data[order(-heatmap.data$width),]  # make sure the narrowest plots are drawn as last
    narrow.cn <-  heatmap.data$width < min.with
    flank <- (min.with - heatmap.data$width[narrow.cn])/2
    heatmap.data$chromEnd[narrow.cn] <- heatmap.data$chromEnd[narrow.cn ] + flank
    heatmap.data$chromStart[narrow.cn ] <- heatmap.data$chromStart[narrow.cn ] - flank
    heatmap.data$chromStart[heatmap.data$chromStart<0] <- 0
    
    heatmap.data <- RCircos.Get.Plot.Data.nosort(heatmap.data, "plot")
    heatmap.data1 <- RCircos.Get.Plot.Data.nosort(data.frame(Chromosome=heatmap.data$Chromosome, chromStart=heatmap.data$chromStart, chromEnd=heatmap.data$chromStart), "plot")
    heatmap.data2 <- RCircos.Get.Plot.Data.nosort(data.frame(Chromosome=heatmap.data$Chromosome, chromStart=heatmap.data$chromEnd, chromEnd=heatmap.data$chromEnd), "plot")

    
    if ((length(heatmap.ranges)==1) && (is.na(heatmap.ranges))) {
        ColorLevel <- RCircos.Par$heatmap.ranges
    } else {
        ColorLevel <- heatmap.ranges
    }

    if ((length(heatmap.color)==1) && (is.na(heatmap.color))) {
        ColorRamp <- RCircos.Get.Heatmap.ColorScales(RCircos.Par$heatmap.color)
    } else {
        
    }
    
    


    
    columns <- 5:(ncol(heatmap.data) - 1)
    min.value <- min(as.matrix(heatmap.data[, columns]))
    max.value <- max(as.matrix(heatmap.data[, columns]))


    
    heatmap.locations1 <- as.numeric(heatmap.data1[, ncol(heatmap.data2)])
    heatmap.locations2 <- as.numeric(heatmap.data2[, ncol(heatmap.data2)])
    
    start <- heatmap.locations1 # -  RCircos.Par$heatmap.width/2
    end <- heatmap.locations2 # + RCircos.Par$heatmap.width/2
    data.chroms <- as.character(heatmap.data[, 1])
    chromosomes <- unique(data.chroms)
    cyto.chroms <- as.character(RCircos.Cyto$Chromosome)

    for (a.chr in 1:length(chromosomes)) {
        cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
        locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]) # chromosome locations
        chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]] # chromosome start
        chr.end <- max(locations) # chromosome end
        data.rows <- which(data.chroms == chromosomes[a.chr]) # points on this chromosome
        start[data.rows[start[data.rows] < chr.start]] <- chr.start # chromosome starts for each point
        end[data.rows[end[data.rows] > chr.end]] <- chr.end # chromosome end for each point
    }
    
    locations <- RCircos.Track.Positions.my(side, track.num)  # positions
    out.pos <- locations[1]
    in.pos <- locations[2]
    chroms <- unique(RCircos.Cyto$Chromosome)
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
            ]
        the.start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        the.end <- the.chr$Location[nrow(the.chr)]
        polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
            RCircos.Pos[the.end:the.start, 1] * in.pos)
        polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
            RCircos.Pos[the.end:the.start, 2] * in.pos)
        polygon(polygon.x, polygon.y, col = "white",  border = RCircos.Par$grid.line.color, lwd=0.3)
    }

    
    heatmap.value <- as.numeric(heatmap.data[, data.col])
    for (a.point in 1:length(heatmap.value)) {
        
        the.level <- which(ColorLevel <= heatmap.value[a.point])
        cell.color <- heatmap.color[max(the.level)] # establish the color
        
        the.start <- start[a.point]
        the.end <- end[a.point]
        if (is.na(the.start) |  is.na(the.end)) {
            browser()
        }
        
        polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, RCircos.Pos[the.end:the.start, 1] * in.pos)
        polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, RCircos.Pos[the.end:the.start, 2] * in.pos)
        polygon(polygon.x, polygon.y, col = cell.color, border = NA)
    }

}
