RCircos.Validate.Genomic.Data.my <- function (genomic.data, plot.type = c("plot", "link")) 
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    plot.type <- tolower(plot.type)
    if (plot.type == "plot") {
        chrom.col <- 1
    }
    else if (plot.type == "link") {
        chrom.col <- c(1, 4)
    }
    else {
        stop("Plot type must be \"plot\" or \"line\"")
    }
    for (a.col in 1:length(chrom.col)) {
        the.col <- chrom.col[a.col]
        genomic.data[, the.col] <- as.character(genomic.data[, 
            the.col])
        for (a.row in 1:nrow(genomic.data)) {
            if (length(grep("chr", genomic.data[a.row, the.col])) == 
                0) {
                genomic.data[a.row, the.col] <- paste("chr", 
                  genomic.data[a.row, the.col], sep = "")
            }
        }
        cyto.chroms <- unique(as.character(RCircos.Cyto$Chromosome))
        data.chroms <- unique(as.character(genomic.data[, the.col]))
        if (sum(data.chroms %in% cyto.chroms) < length(data.chroms)) {
            cat(paste("Some chromosomes are in genomic data only", 
                "and have been removed.\n\n"))
            all.chroms <- as.character(genomic.data[, the.col])
            genomic.data <- genomic.data[all.chroms %in% cyto.chroms, ]
        }
        data.chroms <- unique(as.character(genomic.data[, the.col]))
        if (min(genomic.data[, the.col + 1]) < 0) {
            stop("Error! chromStart position less than 0.")
        }
        if (min(genomic.data[, the.col + 2]) < 0) {
            stop("Error! chromEnd position less than 0.")
        }
        for (a.chr in 1:length(data.chroms)) {
            the.chr <- data.chroms[a.chr]
            in.data <- genomic.data[genomic.data[, the.col] == the.chr, ]
            cyto.data <- RCircos.Cyto[grep(the.chr, RCircos.Cyto$Chromosome), 
                ]

            bad.rows <- in.data[, the.col + 1] > max(cyto.data[, 3])
            in.data[bad.rows, the.col + 1] <- max(cyto.data[, 3])
            bad.rows <- in.data[, the.col + 2] > max(cyto.data[, 3])
            in.data[bad.rows, the.col + 2] <- max(cyto.data[, 3])

            genomic.data[genomic.data[, the.col] == the.chr, ] <- in.data 
            
            if (max(in.data[, the.col + 1]) > max(cyto.data[, 3]) | max(in.data[, the.col + 2]) > max(cyto.data[, 3])) {
                cat(paste(the.chr, max(in.data[, 2]), max(in.data[, 3]), "\n"))
                stop("Error! Location is outside of chromosome length.")
            }
        }
        for (a.row in 1:nrow(genomic.data)) {
            if (genomic.data[a.row, the.col + 1] > genomic.data[a.row, 
                the.col + 2]) {
                cat("chromStart greater than chromEnd.\n")
                stop(paste("Row:", a.row, genomic.data[a.row, 
                  2], genomic.data[a.row, 3]))
            }
        }
    }
    return(genomic.data)
}
