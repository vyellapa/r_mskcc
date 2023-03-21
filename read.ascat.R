# function for importing the result of ASCAT or ASCAT NGS.
# arguments
# FILE.CN - full path the the ascat output (extension .ascat.summary.csv or ascat_ngs.summary.csv)
read.ascat <- function(FILE.CN) {

	    cv.data <- read.table(FILE.CN,header=FALSE, sep=',') # ASCAT
            cd.data <- cv.data[,2:8]
            cat( paste( dim(cv.data)[1], ' copy-number segments \n'))
            names(cv.data) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
            cv.data$seg_no <- NULL
            cv.data$Chromosome <- as.character(cv.data$Chromosome)
            cv.data$Chromosome[cv.data$Chromosome=='23'] <- 'X'
            cv.data$Chromosome[cv.data$Chromosome=='24'] <- 'Y'
            
            if (nrow(cv.data)>0) {
                cv.data$Chromosome <- paste('chr', cv.data$Chromosome,sep='')
            }

            
            cv.data$major.copy.number.inTumour <- cv.data$total.copy.number.inTumour - cv.data$minor.copy.number.inTumour
            
            
            cv.data$major.copy.number.inTumour.temp <- pmax(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
            cv.data$minor.copy.number.inTumour.temp <- pmin(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
            
            cv.data$major.copy.number.inTumour <- cv.data$major.copy.number.inTumour.temp
            cv.data$minor.copy.number.inTumour <- cv.data$minor.copy.number.inTumour.temp

            cv.data
}
