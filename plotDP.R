plot1D_DP <- function (density, polygon.data, pngFile = NA, dp.driver.genes = NULL, density.from = 0, 
    x.max = NA, y.max = NA, y = NULL, N = NULL, mutationCopyNumber = NULL, 
    no.chrs.bearing.mut = NULL, samplename = "", CALR = numeric(0), 
    cluster.locations = NULL, mutation.assignments = NULL) 
{
    if (!is.na(pngFile)) {
        png(filename = pngFile, , width = 1500, height = 1000)
    }
    xlabel = "Mutation Copy Number"
    if (is.null(mutationCopyNumber)) {
        print("No mutationCopyNumber. Using mutation burden")
        if (is.null(y) | is.null(N)) {
            print("When not supplying mutationCopyNumber, y (mutCount) and N (totalCount) are required")
            q(save = "no", status = 1)
        }
        mutationCopyNumber = y/N
        xlabel = "Mutation Burden"
    }
    if (!is.null(no.chrs.bearing.mut)) {
        mutationCopyNumber = mutationCopyNumber/no.chrs.bearing.mut
        xlabel = "Fraction of Tumour Cells"
    }
    xx = density[, 1]
    yy = density[, 2]
    if (is.na(y.max)) {
        y.max = ceiling(max(polygon.data))
    }
    hist(mutationCopyNumber[mutationCopyNumber <= x.max], breaks = seq(-0.1, 
        x.max, 0.025), col = "lightgrey", freq = FALSE, xlab = xlabel, 
        main = "", ylim = c(0, y.max), cex.axis = 1, cex.lab = 1)
    polygon(c(xx, rev(xx)), polygon.data, border = "plum4", col = cm.colors(1, 
        alpha = 0.3))
    lines(xx, yy, col = "plum4", lwd = 1)
    title('Density dist. of corrected VAF', cex.main = 1.2)
    if(!is.null(dp.driver.genes)) {
        abline(v=dp.driver.genes$AnnotateDPPosition, col='red')
        text(x=dp.driver.genes$AnnotateDPPosition, y=seq(from=max(yy), to=1, -1.1)[1:nrow(dp.driver.genes)], labels=dp.driver.genes$Gene, col='red')
    }
    if (!is.null(cluster.locations) & !is.null(mutation.assignments)) {
        assign.counts = table(mutation.assignments)
        for (cluster in unique(mutation.assignments)) {
            x = cluster.locations[cluster.locations[, 1] == cluster, 
                2]
            lines(x = c(x, x), y = c(0, y.max), col = "black", 
                lwd = 3)
            text(paste("Cluster", cluster, sep = " "), x = x + 
                0.01, y = (9/10) * y.max, adj = c(0, 0), cex = 2)
            text(paste(as.numeric(assign.counts[names(assign.counts) == 
                as.character(cluster)]), "mutations", sep = " "), 
                x = x + 0.01, y = (9/10) * y.max - 0.35, adj = c(0, 
                  0), cex = 2)
        }
    }
    if (length(CALR) > 0) {
        x.index = sapply(1:length(CALR), function(i) {
            which.min(abs(CALR[i] - xx))
        })
        CALR.yvals = (polygon.data[x.index] + polygon.data[2 * 
            length(xx) - x.index])/2
        points(CALR, CALR.yvals, pch = 20, col = "red", cex = 3)
    }
    if (!is.na(pngFile)) {
        dev.off()
    }
}
