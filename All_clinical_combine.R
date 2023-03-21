#Script to make combined dataset for patient characteristics at the time of sampling
#SMM, NDMM and RRMM.
# v. 21.3.2018

library(plyr)
library(dplyr)
library(openxlsx)
library(stringr)
library(lubridate)

setwd('/Users/rustade/Documents/MSKCC/Projects/222/All_clinical/combined')

set.seed(353)
options(scipen=999)

#Import
NDMM <- read.xlsx('NDMM_dataline_manual.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))
SMM <- read.xlsx('SMM_dataline_combined.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))
RRMM <- read.xlsx('../../RRMM/datafiles/RRMM_HOTB_Manual.xlsx', detectDates = TRUE, na.strings = c('NA', 'ND'))

#Change variables

NDMM$disease_stage <- 'NDMM'
names(NDMM) <- gsub('base_', '', names(NDMM))
names(NDMM) <- gsub('sample_', '', names(NDMM))
SMM$disease_stage <- 'SMM'
names(SMM) <- gsub('base_', '', names(SMM))
names(SMM) <- gsub('sample_', '', names(SMM))
names(RRMM) <- gsub('sample_', '', names(RRMM))

variables <- names(SMM)
rem <- c("progressed", "progression_date")
variables <- variables[! variables %in% rem]

combined <- rbind(SMM[variables], 
                      NDMM[variables],
                      RRMM[variables]) 
names(combined)[1] <- 'sample_ID'

#changing variables
combined[-c(1:5)] <- sapply(combined[-c(1:5)], tolower)
combined$IFE_urine <- gsub('lambda and kappa', 'kappa and lambda', combined$IFE_urine)

combined <- mutate(combined, ISS = factor(ifelse(beta2m >= 5.5, 3, 
                                         ifelse(beta2m < 3.5 & alb >= 3.5, 1, 2)), levels = c(1:3)))

combined <- mutate(combined, FLC_involved_class = ifelse(IFE_light != 'negative', IFE_light,
                                                   ifelse(kappa>lambda & kappa > 1.94, 'kappa', 
                                                          ifelse(lambda>kappa & lambda > 2.63, 'lambda', 'none'))))
combined$IFE_heavy <- gsub('negative', 'none', combined$IFE_heavy)
combined <- mutate(combined, FLC_involved_value = as.numeric(ifelse(FLC_involved_class == 'kappa', kappa,
                                                         ifelse(FLC_involved_class == 'lambda', lambda, NA))))
combined <- mutate(combined, FLC_uninvolved_value = as.numeric(ifelse(FLC_involved_class == 'kappa', lambda,
                                                           ifelse(FLC_involved_class == 'lambda', kappa, NA))))
combined <- mutate(combined, FLC_uninvolved_value = ifelse(FLC_uninvolved_value == 0, 0.05, FLC_uninvolved_value)) #Set 0 values to detection limit.
combined <- mutate(combined, FLC_ratio = round(FLC_involved_value/FLC_uninvolved_value,0))

combined$aspirate <- as.numeric(combined$aspirate) #change strange formats
combined <- mutate(combined, biopsy = as.numeric(ifelse(grepl('<', biopsy), 
                                                        round(rnorm(1, mean=2.5),0), 
                                                        biopsy))) # impute normal bone marrows as random normal with mean 2.5 and sd=1

combined <- mutate(combined, IgA = round(as.numeric(ifelse(grepl('<', IgA), 0, IgA)),0)) #unmeasurable set to 0
combined <- mutate(combined, IgG = round(as.numeric(ifelse(grepl('<', IgG), 0, IgG)),0)) #unmeasurable set to 0
combined <- mutate(combined, IgM = round(as.numeric(ifelse(grepl('<', IgM), 0, IgM)),0)) #unmeasurable set to 0
combined <- mutate(combined, eGFR = round(as.numeric(ifelse(grepl('>', eGFR), 60, eGFR)),0)) #>60 set to 60
combined$SPEP <- as.numeric(combined$SPEP)
combined$UPEP <- as.numeric(combined$UPEP)
combined$ANC <- as.numeric(combined$ANC)
combined$CRP <- as.numeric(combined$CRP)
combined$WBC <- as.numeric(combined$WBC)
combined$ANC <- as.numeric(combined$ANC)
combined$kappa <- as.numeric(combined$kappa)
combined$lambda <- as.numeric(combined$lambda)
combined$plt <- as.numeric(combined$plt)
combined$alb <- as.numeric(combined$alb)
combined$beta2m <- as.numeric(combined$beta2m)
combined$ca <- as.numeric(combined$ca)
combined$crea <- as.numeric(combined$crea)
combined$hgb <- as.numeric(combined$hgb)
combined$LDH_value<- as.numeric(combined$LDH_value)
combined <- mutate(combined, gender = ifelse(gender == 'm', 'male',
                                             ifelse(gender == 'f', 'female', gender)))
combined$gender <- gsub("(^)([[:alpha:]])", "\\U\\2", combined$gender, perl=TRUE)
combined$FLC_involved_class <- gsub("(^)([[:alpha:]])", "\\U\\2", combined$FLC_involved_class, perl=TRUE)
combined <- mutate(combined, IFE_heavy = ifelse(IFE_heavy == 'igg', 'IgG',
                                                ifelse(IFE_heavy == 'iga', 'IgA',
                                                       ifelse(IFE_heavy == 'igd', 'IgD', 'None'))))

gsub("igg", "IgG", combined$IFE_heavy, perl=TRUE)

# Remove duplicates
combined <- combined[!duplicated(combined[c('MRN', 'sample_ID')]),]

## Output
write.xlsx(combined, 'HOTB_clinical_combined.xlsx')
