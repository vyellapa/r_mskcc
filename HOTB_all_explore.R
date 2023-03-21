setwd('/Users/rustade/Documents/MSKCC/Projects')

library(openxlsx)
library(dplyr)
library(stringr)

HOTB_all <- read.xlsx('./HOTB_projects_all/HOTB_all_master.xlsx', detectDates = T)
HOTB_key <- read.xlsx('./HOTB_projects_all/HOTB_all_key.xlsx', detectDates = T)
invivoscribe <- read.xlsx('./invivoscribe_MRD/VDJ_data/summary.xlsx', detectDates = TRUE)
malin_adaptive <- read.xlsx('./invivoscribe_MRD/Adaptive_VDJ_clinical_Malin.xlsx', detectDates = TRUE)

HOTB_key$MRN <- str_pad(HOTB_key$MRN, 8, pad='0')

## Identify which patients are in invivoscribe file, but not in myType files. 

not_in_myType$MRN <- NA
for(i in 1:nrow(not_in_myType)){
        not_in_myType$MRN[i] <- HOTB_key$MRN[which(HOTB_key$sample_ID == not_in_myType$sample_ID[i])]
}

write.xlsx(not_in_myType, 'invivoscribe_missing_clinical.xlsx', asTable = F)

malin_SMM <- malin_adaptive[malin_adaptive$Sample_ID %in% SMM$ID,c('Sample_ID', 'Collection_date', 'BMPC', 'BMPC_comment')]
malin_SMM$BMPC <- round(as.numeric(malin_SMM$BMPC)*100,0)
malin_SMM$even_aspirate <- NA
for(i in 1:nrow(malin_SMM)){
        malin_SMM$even_aspirate[i] <- SMM$aspirate[which(SMM$ID == malin_SMM$Sample_ID[i])]
}
malin_SMM <- mutate(malin_SMM, aspirate_OK = ifelse(even_aspirate == BMPC, 'OK', 'NO'))

write.xlsx(malin_SMM, 'malin_SMM.xlsx', asTable = F)



RRMM_not_malin <- RRMM[!RRMM$ID %in% malin_RRMM$Sample_ID,]