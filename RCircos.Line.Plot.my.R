RCircos.Line.Plot.my <- function (line.data, data.col, track.num, side, lineCol, lineType) 
{
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()

    line.data1 <- RCircos.Get.Plot.Data.nosort(data.frame(chromosome=line.data$chromosome, start=line.data$start, end=line.data$start), "plot") # add a last column with integer locations
    line.data2 <- RCircos.Get.Plot.Data.nosort(data.frame(chromosome=line.data$chromosome, start=line.data$end, end=line.data$end), "plot")
    # get locations of the tracks
    locations <- RCircos.Track.Positions.my(side, track.num)
    
    out.pos <- locations[1] # posiitons of the track
    in.pos <- locations[2] # position of the track

    
    if (min(as.numeric(line.data[, data.col])) >= 0) {
        point.bottom <- in.pos
        
        data.ceiling <- max(line.data[, data.col]) # data between 0 and 5
    }
    else {
        point.bottom <- in.pos + (RCircos.Par$track.height/2)
        data.ceiling <- 3 # data between -5 and 5
    }
    sub.height <- out.pos - point.bottom
    # line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color)

    RCircos.Track.Outline.my(out.pos, in.pos, RCircos.Par$sub.tracks)

    for (a.point in 1:(nrow(line.data))) {
        point.one <- line.data1[a.point, ncol(line.data1)] # integer location of point  1
        point.two <- line.data2[a.point, ncol(line.data2)] # integer location of point 2

        
        # cut the values if needed
        # if (line.data[a.point, 1] != line.data[a.point + 1, 1]) {
        #     next
        # }
        if (line.data[a.point, data.col] > data.ceiling) {
            value.one <- data.ceiling
        }
        else if (line.data[a.point, data.col] < (-1 * data.ceiling)) {
            value.one <- data.ceiling * -1
        }
        else {
            value.one <- line.data[a.point, data.col]
        }

        if (line.data[a.point , data.col] > data.ceiling) {
            value.two <- data.ceiling
        }
        else if (line.data[a.point , data.col] < (-1 * data.ceiling)) {
            value.two <- data.ceiling * -1
        }
        else {
            value.two <- line.data[a.point, data.col]
        }
        
        height.one <- point.bottom + value.one/data.ceiling * sub.height # scale the y values
        height.two <- point.bottom + value.two/data.ceiling * sub.height # scale the y values

        # height <- out.pos - a.line * subtrack.height
        # lines(RCircos.Pos[start:end, 1] * height, RCircos.Pos[start:end, 2] * height, col = RCircos.Par$grid.line.color, lwd=0.3)

        lines(c(RCircos.Pos[point.one:point.two, 1] * height.one), # xs
              c(RCircos.Pos[point.one:point.two, 2] * height.one), # ys, RCircos.Pos[point.one, *] is always 1 anyway 
              col = lineCol[a.point], lty=lineType, lwd=1.5)
 
    }
}
