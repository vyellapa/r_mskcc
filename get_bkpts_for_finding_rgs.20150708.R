library(changepoint)

BKPT_SLOP = 2000

input_bg = commandArgs(T)[1]
sample = sub("\\..+", "", input_bg)

d = read.table(input_bg, header = F, sep = "\t", stringsAsFactor = F)
d = d[order(d[,2]), ]
res = cpt.meanvar(d[,4], test.stat = "Poisson", method = "BinSeg", Q = 250, penalty = "Manual", pen.val = 1000)

pdf(paste0(sample, ".segs.pdf"))
plot(ts(d[,4]), main = sample, ylab = "Coverage")
abline(v = cpts(res), col = "grey")
segments(c(1, cpts(res)), param.est(res)[[1]], c(cpts(res), 4000), col = "blue", lwd = 2)
dev.off()

lambda = param.est(res)$lambda

out = data.frame(
    "14",
    d[cpts(res), 3] - BKPT_SLOP,
    d[cpts(res), 3] + BKPT_SLOP,
    round(lambda[-length(lambda)]),
    round(lambda[-1]),
    ifelse(lambda[-1] < lambda[-length(lambda)], "+", "-")
)   

write.table(
    out,
    paste0(sample, ".igh_region.bkpts.bed"),
    col.names = F,
    row.names = F,
    sep = "\t",
    quote = F 
)   

