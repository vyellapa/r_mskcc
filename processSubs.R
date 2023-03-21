processSubs <- function(subs) {

# new style
scatter.data  <-  calcIntermutDist(subs)
scatter.data$mutType <- paste(scatter.data$ref_base_pyrimidine_context,'>',scatter.data$mutant_base_pyrimidine_context ,sep='')

scatter.colors <- rep("", nrow(scatter.data))
scatter.colors[scatter.data$mutType=="C>A"] <- "royalblue"
scatter.colors[scatter.data$mutType=="C>G"] <- "black"
scatter.colors[scatter.data$mutType=="C>T"] <- "red"
scatter.colors[scatter.data$mutType=="T>A"] <- "grey"
scatter.colors[scatter.data$mutType=="T>C"] <- "green2"
scatter.colors[scatter.data$mutType=="T>G"] <- "hotpink"

result <- list()
result$scatter.colors <- scatter.colors
result$scatter.data <- scatter.data

result

# old style    
#subs$mutType <- paste(subs$ref_base_pyrimidine_context,'>',subs$mutant_base_pyrimidine_context ,sep='')
# inter-mutation distances computed within each mutation type
#subs.CtoA <- subset(subs, subset=mutType=="C>A")
#subs.CtoG <- subset(subs, subset=mutType=="C>G")
#subs.CtoT <- subset(subs, subset=mutType=="C>T")
#subs.TtoA <- subset(subs, subset=mutType=="T>A")
#subs.TtoC <- subset(subs, subset=mutType=="T>C")
#subs.TtoG <- subset(subs, subset=mutType=="T>G")

#subs.CtoA.processed <- calcIntermutDist(subs.CtoA)
#subs.CtoG.processed  <-  calcIntermutDist(subs.CtoG)
#subs.CtoT.processed  <-  calcIntermutDist(subs.CtoT)
#subs.TtoA.processed  <-  calcIntermutDist(subs.TtoA)
#subs.TtoC.processed  <-  calcIntermutDist(subs.TtoC)
#subs.TtoG.processed  <-  calcIntermutDist(subs.TtoG)

#scatter.data <- rbind(
#    subs.CtoA.processed,
#    subs.CtoG.processed,
#    subs.CtoT.processed,
#    subs.TtoA.processed,
#    subs.TtoC.processed,
#    subs.TtoG.processed     
#    )

#scatter.colors <- c(
#    rep("royalblue",nrow(subs.CtoA.processed)),
#    rep("black",nrow(subs.CtoG.processed)),
#    rep("red",nrow(subs.CtoT.processed)),
#    rep("grey",nrow(subs.TtoA.processed)),
#    rep("green2",nrow(subs.TtoC.processed)),
#    rep("hotpink",nrow(subs.TtoG.processed))
#    )



} 
