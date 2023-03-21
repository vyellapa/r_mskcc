RCircos.Track.Positions.my <- function (side, track.num) 
{
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
    side <- tolower(side)
    if (side == "in") {
         # out.pos <- RCircos.Par$track.in.start - (track.num - 
         #    1) * one.track
         # in.pos <- out.pos - RCircos.Par$track.height
        out.pos <- RCircos.Par$track.in.start 
        if (track.num>1) {
            out.pos <- RCircos.Par$track.in.start - sum( RCircos.Par$track.heights[1:(track.num-1)]) -          
                 sum(RCircos.Par$track.padding[1:(track.num - 1)])
        }
        
        in.pos <- out.pos - RCircos.Par $track.heights[track.num]       
    }
    else if (side == "out") {
        in.pos <- RCircos.Par$track.out.start + (track.num - 
            1) * one.track
        out.pos <- in.pos + RCircos.Par$track.height
    }
    else {
        stop("Incorrect track location. It must be \"in\" or \"out\".")
    }
    return(c(out.loc = out.pos, in.loc = in.pos))
}
