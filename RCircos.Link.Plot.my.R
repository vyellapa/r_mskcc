RCircos.Link.Plot.my <- function (link.data, track.num, by.chromosome = FALSE, link.colors=NA, lwd=0.5) 
{

    if (length(link.colors)==1) {
        link.colors <- rep('BurlyWood', nrow(link.data))
    }

    link.data$col <- link.colors

    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type = "link")
    locations <- RCircos.Track.Positions.my('in', track.num)
    start <- locations[['out.loc']]
    base.positions <- RCircos.Pos * start
    data.points <- matrix(rep(0, nrow(link.data) * 2), ncol = 2)
    for (a.link in 1:nrow(link.data)) {
        data.points[a.link, 1] <- RCircos.Data.Point(link.data[a.link, 
            1], link.data[a.link, 2])
        data.points[a.link, 2] <- RCircos.Data.Point(link.data[a.link, 
            4], link.data[a.link, 5])
        if (data.points[a.link, 1] == 0 || data.points[a.link, 
            2] == 0) {
            print("Error in chromosome locations ...")
            break
        }
    }
    # link.colors <- RCircos.Get.Link.Colors(link.data, by.chromosome)
    for (a.link in 1:nrow(data.points)) {
        point.one <- data.points[a.link, 1]
        point.two <- data.points[a.link, 2]
        if (point.one > point.two) {
            point.one <- data.points[a.link, 2]
            point.two <- data.points[a.link, 1]
        }
        P0 <- as.numeric(base.positions[point.one, ])
        P2 <- as.numeric(base.positions[point.two, ])
        links <- RCircos.Link.Line(P0, P2)
        # lines(links$pos.x, links$pos.y, type = "l", col = link.colors[a.link])
        lines(links$pos.x, links$pos.y, type = "l", col = link.data$col[a.link], lwd=lwd)
    }
}
