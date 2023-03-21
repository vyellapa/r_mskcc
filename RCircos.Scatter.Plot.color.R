RCircos.Scatter.Plot.color <- function (scatter.data, data.col, track.num, side, by.fold = 0,  scatter.colors, draw.bg =TRUE, draw.scale=FALSE, no.sort=FALSE, data.ceiling=NA) 
{

    # scatter.data.original <- scatter.data

    # scatter.data.original.small <- scatter.data.original[c(1:10, 601:610),]
    # no.points.small <- nrow(scatter.data.original.small )
    # new.order.small <- sample(1:no.points.small, no.points.small)
    # scatter.data.original.small[new.order.small, ]
    # scatter.data.original.small$rate[scatter.data.original.small$col=='grey']

    # scatter.data <- scatter.data.original 

    no.points <- nrow(scatter.data)
    if (no.sort) {
        new.order <- 1:no.points
    } else {
        new.order <- sample(1:no.points, no.points)
    }

    # apply new ordering
    scatter.data <- scatter.data[new.order,]
    scatter.colors <- scatter.colors[new.order]
     
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    scatter.data <- RCircos.Get.Plot.Data.nosort(scatter.data, "plot")

    
    locations <- RCircos.Track.Positions.my(side, track.num)
    out.pos <- locations[1]
    in.pos <- locations[2]
    if (min(as.numeric(scatter.data[, data.col])) >= 0) {
        point.bottom <- in.pos
        if (is.na(data.ceiling)) {
            data.ceiling <- max(scatter.data[, data.col])
        }
    }
    else {
        point.bottom <- in.pos + (RCircos.Par$track.height/2)
        if (is.na(data.ceiling)) {
            data.ceiling <- 5
        }
    }
    sub.height <- out.pos - point.bottom

    if (draw.bg) {
        RCircos.Track.Outline.my(out.pos, in.pos, RCircos.Par$sub.tracks)
    }
    if (draw.scale) {
        text(RCircos.Pos[1, 1] * locations[1], RCircos.Pos[1, 2] * locations[1], round(data.ceiling), cex=0.5)
        text(RCircos.Pos[1, 1] * locations[2], RCircos.Pos[1, 2] * locations[2], '0', cex=0.5)
    }
    
    
    for (a.point in 1:nrow(scatter.data)) {
        the.point <- scatter.data[a.point, ncol(scatter.data)]
        color <- scatter.colors[a.point]
        if (scatter.data[a.point, data.col] > data.ceiling) {
            the.value <- data.ceiling
        }
        else if (scatter.data[a.point, data.col] < (-1 * data.ceiling)) {
            the.value <- data.ceiling * -1
        }
        else {
            the.value <- scatter.data[a.point, data.col]
        }

        if (by.fold > 0) {
            if (the.value >= by.fold) {
                color <- "red"
            }
            else if (the.value <= -by.fold) {
                color <- "blue"
            }
            else {
                color <- "black"
            }
        }
        height <- point.bottom + the.value/data.ceiling * sub.height
        points(RCircos.Pos[the.point, 1] * height,
               RCircos.Pos[the.point, 2] * height,
               col = color,
               pch = RCircos.Par$point.type, 
               cex = RCircos.Par$point.size)
    }
}
