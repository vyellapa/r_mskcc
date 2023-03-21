## Import Project 222 datafiles from caveman, pindel, delly, brass and cnvkit. 
## Subset data based on leukgen ID to include only 58 patients with RRMM.
## Store output as separate files in directory p222.results.RRMM
## v.1.3.2018

#run setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.seq'); source('/Users/rustade/Documents/MSKCC/RRMM/analysis.seq/p222.subset.R')

#Packages
library(openxlsx)

#Workdir
setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.seq')

#Import data
caveman.p222 <- read.csv('../datafiles/p222.results/caveman.csv', stringsAsFactors = FALSE)
pindel.p222 <- read.csv('../datafiles/p222.results/Pindel.csv', stringsAsFactors = FALSE)
delly.p222 <- read.csv('../datafiles/p222.results/delly.100bp.trans.oct30.csv', stringsAsFactors = FALSE)
brass.p222 <- read.csv('../datafiles/p222.results/brassR3.50bp.trans.oct30.csv', stringsAsFactors = FALSE)
cnvkit.p222 <- read.csv('../datafiles/p222.results/cnvkit.seg.csv', stringsAsFactors = FALSE)

#Import patient/sample key
key <- read.xlsx('../datafiles/FISH.xlsx', cols = c(1, 3))

#Uniform ID variable names
names(key) <- c('sampleID', 'sampleLeukid')
names(caveman.p222)[2] <- 'sampleLeukid'
caveman.p222$sampleID <- 'temp'
names(pindel.p222)[2] <- 'sampleLeukid'
pindel.p222$sampleID <- 'temp'
names(delly.p222)[13] <- 'sampleLeukid'
delly.p222$sampleID <- 'temp'
names(brass.p222)[9] <- 'sampleLeukid'
brass.p222$sampleID <- 'temp'
names(cnvkit.p222)[1] <- 'sampleLeukid'
cnvkit.p222$sampleID <- 'temp'

#Define function to subset files
subsetRRMM <- function(p222){
        output <- data.frame()
        for(i in 1:nrow(p222)) {
                ID <- p222$sampleLeukid[i]
                if(substr(ID, 1, nchar(ID)-5) %in% key$sampleLeukid){
                        p222$sampleID[i] <- as.character(key$sampleID[which(key$sampleLeukid == substr(ID, 1, nchar(ID)-5))])
                        if(nrow(output) >0){
                                output <- rbind(output, p222[i,])
                        } else {
                                output <- p222[i,]
                        }
                }
        } 
        return(output)
}

#Subset files
caveman.RRMM <- subsetRRMM(caveman.p222)
pindel.RRMM <- subsetRRMM(pindel.p222)
delly.RRMM <- subsetRRMM(delly.p222)
brass.RRMM <- subsetRRMM(brass.p222)
cnvkit.RRMM <- subsetRRMM(cnvkit.p222)

#Write output
write.csv(caveman.RRMM, '../datafiles/p222.results.RRMM/caveman.csv', row.names = FALSE)
write.csv(pindel.RRMM, '../datafiles/p222.results.RRMM/pindel.csv', row.names = FALSE)
write.csv(delly.RRMM, '../datafiles/p222.results.RRMM/delly.csv', row.names = FALSE)
write.csv(brass.RRMM, '../datafiles/p222.results.RRMM/brass.csv', row.names = FALSE)
write.csv(cnvkit.RRMM, '../datafiles/p222.results.RRMM/cnvkit.csv', row.names = FALSE)


