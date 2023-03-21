#merge caveman annotation files from project 222

setwd('/Users/rustade/Documents/MSKCC/Projects/222/All_myType/caveman')
library(openxlsx)
library(plyr)
library(dplyr)
library(stringr)

#Import and format
ER <- read.xlsx('caveman_manual.xlsx')
names(ER)[13:15] <- paste(names(ER)[13:15], 'ER', sep = '_')
ER <- as.data.frame(sapply(ER, str_trim, side = 'both'), stringsAsFactors = F)
ER <- mutate(ER, MANUAL_OK_ER = ifelse(MANUAL_OK_ER %in% c('OK', 'NO'), MANUAL_OK_ER, 'UNCERTAIN'))

MH <- read.xlsx('P222_CAVEMAN_MH.xlsx')
MH <- MH[c('ID_VARIANT', "MANUAL_OK", "MANUAL_ANNOTATION", "MANUAL_COMMENT")]
names(MH)[2:4] <- paste(names(MH)[2:4], 'MH', sep = '_')
MH <- as.data.frame(sapply(MH, str_trim, side = 'both'), stringsAsFactors = F)
MH$MANUAL_OK_MH <- toupper(MH$MANUAL_OK_MH)
MH <- mutate(MH, MANUAL_OK_MH = ifelse(MANUAL_OK_MH %in% c('OK', 'NO'), MANUAL_OK_MH, 'UNCERTAIN'))

TA <- read.xlsx('p222_caveman_TA.xlsx')
TA <- TA[c('ID_VARIANT', "MANUAL_OK", "MANUAL_ANNOTATION", "MANUAL_COMMENT")]
names(TA)[2:4] <- paste(names(TA)[2:4], 'TA', sep = '_')
TA <- as.data.frame(sapply(TA, str_trim, side = 'both'), stringsAsFactors = F)
TA <- mutate(TA, MANUAL_OK_TA = ifelse(MANUAL_OK_TA %in% c('OK', 'NO'), MANUAL_OK_TA, 'UNCERTAIN'))

#Combine
combined <- merge(ER, MH, by = 'ID_VARIANT')
combined <- merge(combined, TA, by = 'ID_VARIANT')

combined$CONSENSUS <- NA
combined$AGREE <- 'AGREE'

#dput(names(combined)) prints names in R vector format

combined <- combined[c("CONSENSUS", "AGREE", "MANUAL_OK_MH", "MANUAL_ANNOTATION_MH", "MANUAL_COMMENT_MH", "MANUAL_OK_ER", "MANUAL_ANNOTATION_ER", 
                       "MANUAL_COMMENT_ER", "MANUAL_OK_TA", "MANUAL_ANNOTATION_TA", "MANUAL_COMMENT_TA", "CHR", "GENE", "cDNA_CHANGE", "PROTEIN_CHANGE", 
                       "TARGET_VAF", "REFERENCE_VAF", "Bolli_2014_gene", "Bolli_2017_oncpos_gene", 
                       "Lohr_2014_gene", "Walker_2015_gene", "Kortuem_2016_BCJ_gene", 
                       "Kortuem_2016_Blood_gene", "BIDIR", "DIRPROP", "EFFECT", "FILTER", 
                       "Bolli_Class", "Bolli_Frequency", "Bolli_Annotation", "Bolli_Positions", 
                       "MMRF_Class", "MMRF_Frequency", "MMRF_VAF", "MMRF_Q25", "MMRF_Q75", 
                       "MMRF_Positions", "HEME_EXACT", "COSMIC", "Normals_Frequency", 
                       "Normals_median_VAF", "MAX_MAF", "FIN_MAF", "NFE_MAF", "AFR_MAF", 
                       "AMR_MAF", "EAS_MAF", "SAS_MAF", "OTHER_MAF", "G1000_MAF", "TARGET_DEPTH", 
                       "REFERENCE_DEPTH", "START", "END", "REF", "ALT", "TARGET_NAME", 
                       "REFERENCE_NAME", "ID_VARIANT")]

#Checking
combined <- mutate(combined, AGREE = ifelse(MANUAL_OK_MH == 'UNCERTAIN' | MANUAL_OK_ER == 'UNCERTAIN' | MANUAL_OK_MH == 'UNCERTAIN', 'CHECK', AGREE))
combined <- mutate(combined, AGREE = ifelse(MANUAL_ANNOTATION_MH == MANUAL_ANNOTATION_ER & MANUAL_ANNOTATION_ER == MANUAL_ANNOTATION_TA, AGREE, 'CHECK'))
combined <- mutate(combined, AGREE = ifelse(is.na(MANUAL_ANNOTATION_MH) | is.na(MANUAL_ANNOTATION_ER) | is.na(MANUAL_ANNOTATION_TA), 'CHECK', AGREE))
combined <- mutate(combined, AGREE = ifelse(AGREE == 'CHECK' & (MANUAL_OK_MH == 'OK' | MANUAL_OK_ER == 'OK' | MANUAL_OK_MH == 'OK'), 'PRIORITY', AGREE))

# Output
#write.xlsx(combined, 'p222_caveman_merged.xlsx', asTable = F)