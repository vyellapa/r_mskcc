RCircos.Gene.Name.Plot.my <- function (gene.data, name.col, track.num, side, colors) 
{
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    gene.data <- RCircos.Get.Plot.Data.nosort(gene.data, "plot")
    gene.data <- RCircos.Get.Gene.Label.Locations(gene.data)
    side <- tolower(side)
    locations <- RCircos.Track.Positions(side, track.num)
    if (side == "in") {
        label.pos <- locations[1]
    }
    else {
        label.pos <- locations[2]
    }
    right.side <- nrow(RCircos.Pos)/2
    text.colors <- RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color)
    for (a.text in 1:nrow(gene.data)) {
        gene.name <- as.character(gene.data[a.text, name.col])
        the.point <- as.numeric(gene.data[a.text, ncol(gene.data)])
        rotation <- RCircos.Pos$degree[the.point]
        if (side == "in") {
            if (the.point <= right.side) {
                text.side <- 2
            }
            else {
                text.side <- 4
            }
        }
        else {
            if (the.point <= right.side) {
                text.side <- 4
            }
            else {
                text.side <- 2
            }
        }

         text(RCircos.Pos[the.point, 1] * label.pos, RCircos.Pos[the.point, 
            2] * label.pos, label = gene.name, pos = text.side, 
            cex = RCircos.Par$text.size, srt = rotation, offset = 0, 
            col = as.character(gene.data$color[a.text]))
    }
}

RCircos.Gene.Name.Plot.my2 <- function (gene.data, name.col, track.num, side) 
{
    side <- tolower(side)
    if (side != "in" && side != "out") {
        stop("Side must be either in or out.\n")
    }
    if (name.col < 4) {
        stop("Column number for gene names must be 4 or bigger.\n")
    }
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    gene.data <- RCircos.Get.Plot.Data(gene.data, "plot")
    gene.data <- RCircos.Get.Gene.Label.Locations(gene.data)
    right.side <- nrow(RCircos.Pos)/2
    locations <- RCircos.Track.Positions(side, track.num)
    the.points <- as.numeric(gene.data[, ncol(gene.data)])
    if (side == "in") {
        label.pos <- locations[1]
        text.side <- rep(4, nrow(gene.data))
        text.side[the.points <= right.side] <- 2
    }
    else {
        label.pos <- locations[2]
        text.side <- rep(2, nrow(gene.data))
        text.side[the.points <= right.side] <- 4
    }
    text.colors <- RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color)
    for (a.text in 1:nrow(gene.data)) {
        gene.name <- as.character(gene.data[a.text, name.col])
        rotation <- RCircos.Pos$degree[the.points[a.text]]
        text(RCircos.Pos[the.points[a.text], 1] * label.pos, 
            RCircos.Pos[the.points[a.text], 2] * label.pos, label = gene.name, 
            pos = text.side[a.text], cex = RCircos.Par$text.size, 
            srt = rotation, offset = 0, col = text.colors[a.text])
    }
}
