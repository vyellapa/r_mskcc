cn_file = commandArgs(T)[1]
rgs_file = commandArgs(T)[2]
s = sub(".+/", "", cn_file)
s = sub("\\..+", "", s)

xlim = c(107350001, 106000000)
igh_annot = read.table("igh_components.ens.txt", header = T, sep = "\t", stringsAsFactors = F)
igh_type = ifelse(igh_annot[,7] == "IGHD", "IGHD", substr(igh_annot[,7], 4, 4)) 
library(RColorBrewer)
cols = brewer.pal(length(unique(igh_type)), "Set1")
names(cols) = c("V", "D", "J", "E", "M", "IGHD", "G", "A")

invert_dir = c("+" = "-", "-" = "+")

cn = read.table(cn_file, header = F, sep = "\t", stringsAsFactors = F)
if (file.info(rgs_file)$size == 0) {
    rgs = NULL
} else {
    rgs = read.table(rgs_file, header = F, sep = "\t", stringsAsFactors = F)
    rgs = rgs[rgs[,1] == "14" & rgs[,4] == "14", ]
}   

pdf(paste0(s, ".igh_region_rgs.pdf"), w = 10, h = 5)
par(mfrow = c(2, 1), mar = c(0, 4, 6, 2) + .1) 
plot(c(), axes = F, xlab = "", ylab = "", main = s, xlim = xlim, ylim = c(0,2))
rect(
    igh_annot[,4],
    0,  
    igh_annot[,5],
    .5, 
    border = cols[igh_type],
    lty = ifelse(igh_type %in% names(cols)[c(1, 3, 5, 7)], 1, 2)
)   
if (!is.null(rgs) && nrow(rgs) > 0) {
    segments(
        c(rowMeans(rgs[, 2:3]), rowMeans(rgs[, 5:6])),
        .5, 
        rep(rowMeans(rgs[, c(2:3, 5:6)]), 2), 
        1   
    )   
#     text(
#         rowMeans(rgs[, c(2:3, 5:6)]),
#         1.1,
#         labels = paste(invert_dir[rgs[,10]], invert_dir[rgs[,9]])
#     )   
}   
legend("topleft", ncol = 4, col = cols, lty = 1:2, legend = names(cols), bty = "n")

par(mar = c(5, 4, 0, 2) + .1) 
plot(c(), xlim = xlim, xlab = "Chromosome 14 position", ylab = "Read depth", main = "", ylim = quantile(cn[,4], c(0, 0.99)), las = 1)
segments(cn[,2], cn[,4], cn[,3], xpd = NA)

dev.off()

