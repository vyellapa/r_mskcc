my.RCircos.Tile.Plot <- function (tile.data, track.num, side, tile.colors=NA) 
{
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    tile.data <- RCircos.Get.Plot.Data(tile.data, "plot")
    the.layer <- 1
    the.chr <- tile.data[1, 1]
    start <- tile.data[1, 2]
    end <- tile.data[1, 3]
    tile.layers <- rep(1, nrow(tile.data))
    if (nrow(tile.data)>1) {
        for (a.row in 2:nrow(tile.data)) {
            if (tile.data[a.row, 2] >= end) {
                the.layer <- 1
                start <- tile.data[a.row, 2]
                end <- tile.data[a.row, 3]
            }
            else if (tile.data[a.row, 1] != the.chr) {
                the.layer <- 1
                the.chr <- tile.data[a.row, 1]
                start <- tile.data[a.row, 2]
                end <- tile.data[a.row, 3]
            }
            else {
                the.layer <- the.layer + 1
                if (tile.data[a.row, 3] > end) {
                    end <- tile.data[a.row, 3]
                }
            }
                                        # tile.layers[a.row] <- the.layer
            tile.layers[a.row] <- 1
        }
    }
    locations <- RCircos.Track.Positions.my(side, track.num)
    out.pos <- locations[1]
    in.pos <- locations[2]
    layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers
    num.layers <- max(tile.layers)
    if (num.layers > RCircos.Par$max.layers) {
        if (side == "in") {
            in.pos <- out.pos - layer.height * num.layers
        }
        else {
            out.pos <- in.pos + layer.height * num.layers
        }
        cat(paste("Tiles plot will use more than one track.", 
            "Please select correct area for next track.\n"))
    }
    if (num.layers < RCircos.Par$max.layers) {
        layer.height <- RCircos.Par$track.height/num.layers
    }
    if (length(tile.colors)==1) {
    tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color)}
    RCircos.Track.Outline.my(out.pos, in.pos, num.layers)
    the.loc <- ncol(tile.data)
    for (a.row in 1:nrow(tile.data)) {
        tile.len <- tile.data[a.row, 3] - tile.data[a.row, 2]
        tile.range <- round(tile.len/RCircos.Par$base.per.unit/2, 
            digits = 0)
        start <- tile.data[a.row, the.loc] - tile.range
        end <- tile.data[a.row, the.loc] + tile.range
        layer.bot <- in.pos
        layer.top <- out.pos
        # layer.bot <- in.pos + layer.height * (tile.layers[a.row] - 1)
        # layer.top <- layer.bot + layer.height * 0.8
        # layer.top <- layer.bot + layer.height
        polygon.x <- c(RCircos.Pos[start:end, 1] * layer.top, 
            RCircos.Pos[end:start, 1] * layer.bot)
        polygon.y <- c(RCircos.Pos[start:end, 2] * layer.top, 
            RCircos.Pos[end:start, 2] * layer.bot)
        polygon(polygon.x, polygon.y, col = tile.colors[a.row], lwd=1, border=tile.colors[a.row])
    }
}
