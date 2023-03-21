########
#R-Script to import all clinical datasets for project 222
#Run: setwd('/Users/rustade/Documents/MSKCC/Projects/222'); source('./All_clinical/import_all.R')

#Set working directory
setwd('/Users/rustade/Documents/MSKCC/Projects/222')

#Import required packages
library(openxlsx)
library(dplyr)

#Import clinical data excel spreadsheet 
invivoscribe <- read.xlsx('../invivoscribe_MRD/summary.xlsx', detectDates = TRUE)
NDMM <- read.xlsx('./NDMM/datafiles/NDMM_HOTB_Manual.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))
SMM <- read.xlsx('./SMM/datafiles/SMM_HOTB_Manual.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))
RRMM <- read.xlsx('./RRMM/datafiles/RRMM_HOTB_Manual.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))
NDMM$cohort <- 'NDMM'; names(NDMM) <- gsub('base_', '', names(NDMM)); names(NDMM) <- gsub('sample_', '', names(NDMM))
SMM$cohort <- 'SMM'; names(SMM) <- gsub('base_', '', names(SMM)); names(SMM) <- gsub('sample_', '', names(SMM))
RRMM$cohort <- 'RRMM'; names(RRMM) <- gsub('sample_', '', names(RRMM))

# Pre-processing variables
RRMM$disease_stage <- RRMM$context

# Create common summary
variables <- names(SMM); rem <- c("progressed", "progression_date"); variables <- variables[! variables %in% rem]

HOTB <- SMM[variables]
HOTB <- rbind(HOTB, 
              NDMM[variables],
              RRMM[variables]) 
names(HOTB)[1] <- 'sample_ID'

for(i in 1:nrow(HOTB)){
        if(HOTB$cohort[i] == HOTB$disease_stage[i]){
                HOTB$disease_stage[i] <- NA
        }
}

HOTB$capture <- NA
HOTB$tube_position <- NA
HOTB$tube_barcode <- NA
HOTB$invivoscribe_BMPC <- NA

for(i in 1:nrow(invivoscribe)){
        ID <- which(HOTB$sample_ID == invivoscribe$sample_ID[i])
        HOTB$capture[ID] <- invivoscribe$capture[i]
        HOTB$tube_position[ID] <- invivoscribe$tube_position[i]
        HOTB$tube_barcode[ID] <- invivoscribe$tube_barcode[i]
        HOTB$invivoscribe_BMPC[ID] <- invivoscribe$BMPC[i]
}

write.xlsx(HOTB, './All_clinical/analysis/HOTB_summary.xlsx', asTable = F)