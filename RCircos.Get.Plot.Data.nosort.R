RCircos.Get.Plot.Data.nosort <- function (genomic.data, plot.type, validate=TRUE) 
{

    if (validate) {
    genomic.data <- RCircos.Validate.Genomic.Data.my(genomic.data, plot.type) }
    data.points <- rep(0, nrow(genomic.data))
    for (a.row in 1:nrow(genomic.data)) {
        if ((a.row %% 1000)==0) {
            cat(paste(a.row,sep=''))
        }
        chromosome <- as.character(genomic.data[a.row, 1])
        location <- round((genomic.data[a.row, 2] + genomic.data[a.row, 
            3])/2, digits = 0)
        data.points[a.row] <- RCircos.Data.Point(chromosome, location)
    }
    genomic.data["Location"] <- data.points
    # genomic.data <- genomic.data[order(genomic.data$Location), ]
    return(genomic.data)
}
