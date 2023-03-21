#Script to create clean output from dataline files and merge with previously collected data.
#For NDMM dataset
# Created 19.3.2018
# Amended 16.01.2019 to calculate age at diagnosis and remove PHI

library(plyr)
library(dplyr)
library(openxlsx)
library(stringr)
library(lubridate)

setwd('/Users/rustade/Documents/MSKCC/Projects/222/All_clinical/combined')

ptinfo <- read.xlsx('../from_dataline/NDMM_ptinfo.xlsx', detectDates = T)
labinfo <- read.xlsx('../from_dataline/NDMM_labinfo.xlsx', detectDates = T)
NDMM_in <- read.xlsx('../../clinical_data_collection/NDMM/datafiles/NDMM_HOTB_Manual.xlsx', detectDates = T)

# merge ptinfo and labinfo
names(ptinfo) <- c("MRN", "review_date", "DOB", 'gender', 'race', 'ethnicity', 'date_sample', 'treatment_date', 'last_visit', 'deceased_date')

labinfo <- join(labinfo, ptinfo, by = 'MRN', type = 'full')

## Create lab variables
labinfo$Base.IFE <- tolower(labinfo$Base.IFE)
labinfo <- mutate(labinfo, base_IFE_heavy = ifelse(grepl('^ig', Base.IFE), str_extract(Base.IFE, '^ig[agd]'), 'negative')) 
labinfo <- mutate(labinfo, base_IFE_light = ifelse(base_IFE_heavy != 'negative' & grepl('kappa|lambda', Base.IFE), str_extract(Base.IFE, 'kappa|lambda'), 
                                             ifelse(grepl('^no ', Base.IFE) | !grepl('kappa|lambda', Base.IFE), 'negative', str_extract(Base.IFE, 'kappa|lambda'))))

labinfo <- mutate(labinfo, base_IFE_heavy = ifelse(is.na(Base.IFE), NA, base_IFE_heavy))
labinfo <- mutate(labinfo, base_IFE_light = ifelse(is.na(Base.IFE), NA, base_IFE_light))


for(i in grep('SPEP', names(labinfo))[1]:grep('SPEP', names(labinfo))[length(grep('SPEP', names(labinfo)))]){
        labinfo[,i] <- ifelse(grepl('^[Nn]', labinfo[,i]), 0, labinfo[,i])
        labinfo[,i] <- as.numeric(labinfo[,i])
}
labinfo <- mutate(labinfo, base_UPEP = rowSums(labinfo[grep('Urine.SPEP', names(labinfo))[1]:grep('Urine.SPEP', names(labinfo))[length(grep('Urine.SPEP', names(labinfo)))]], na.rm = T))
for(i in 1:nrow(labinfo)){
        if(sum(is.na(labinfo[i,grep('Urine.SPEP', names(labinfo))])) == 4){
                labinfo$base_UPEP[i] <- NA
        }
}

labinfo <- mutate(labinfo, base_SPEP = rowSums(labinfo[grep('^SPEP', names(labinfo))[1]:grep('^SPEP', names(labinfo))[length(grep('^SPEP', names(labinfo)))]], na.rm = T))
for(i in 1:nrow(labinfo)){
        if(sum(is.na(labinfo[i,grep('^SPEP', names(labinfo))])) == 7){
                labinfo$base_SPEP[i] <- NA
        }
}

labinfo$IFE.Kappa.Urine <- tolower(labinfo$IFE.Kappa.Urine); labinfo$IFE.Kappa.Urine[is.na(labinfo$IFE.Kappa.Urine)] <- 'missing'
labinfo$IFE.Lambda.Urine <- tolower(labinfo$IFE.Lambda.Urine); labinfo$IFE.Lambda.Urine[is.na(labinfo$IFE.Lambda.Urine)] <- 'missing'

labinfo <- mutate(labinfo, base_IFE_urine = ifelse(IFE.Kappa.Urine == 'missing' & IFE.Lambda.Urine == 'missing', NA,
                                                   ifelse(grepl('sign', IFE.Kappa.Urine) & grepl('sign', IFE.Lambda.Urine), 'kappa and lambda',
                                                          ifelse(!grepl('sign', IFE.Kappa.Urine) & !grepl('sign', IFE.Lambda.Urine), 'negative',
                                                                 ifelse(grepl('sign', IFE.Kappa.Urine), 'kappa', 'lambda')))))

labinfo <- mutate(labinfo, eGFR = ifelse(grepl('BLACK', race), EGFR.African.American, `EGFR.Non-African.American`))
labinfo <- mutate(labinfo, base_kappa = round(as.numeric(Free.Kappa),1))
labinfo <- mutate(labinfo, base_lambda = round(as.numeric(Free.Lambda),1))
labinfo <- mutate(labinfo, CRP = round(ifelse(!is.na(as.numeric(`CRP-H`)), as.numeric(`CRP-H`), as.numeric(CRP)), 2))

labinfo <- mutate(labinfo, LDH_value = ifelse(LDH.Abnormality.Code == 'H', LDH, NA))
labinfo <- mutate(labinfo, LDH = ifelse(LDH.Abnormality.Code == 'H', 'elevated', 'normal'))

labinfo <- mutate(labinfo, plt = Plt)
labinfo <- mutate(labinfo, alb = Alb)
labinfo <- mutate(labinfo, beta2m = Beta2m)
labinfo <- mutate(labinfo, ca = Ca)
labinfo <- mutate(labinfo, crea = Crea)
labinfo <- mutate(labinfo, hgb = HGB)

#create new output file
labinfo <- labinfo[names(labinfo) %in% names(NDMM_in)]
labinfo <- labinfo[-which(names(labinfo) %in% c("date_sample", "treatment_date"))]
NDMM_in <- NDMM_in[-c(which(names(NDMM_in) %in% names(labinfo[-1])))]

NDMM_out <- join(NDMM_in, labinfo, type = 'full', by = 'MRN')

# Calculate age

NDMM_out$DOB <- as.POSIXct(strptime(NDMM_out$DOB, format = '%Y-%m-%d'))
NDMM_out$diagn_MM <- as.POSIXct(strptime(NDMM_out$diagn_MM, format = '%Y-%m-%d'))
NDMM_out$age  <-  year(as.period(interval(NDMM_out$DOB, NDMM_out$diagn_MM))) #age at time of diagnosis

# OS variable

NDMM_out$deceased <- ifelse(is.na(NDMM_out$deceased_date), 'No', 'Yes')

# Final variables without PHI
NDMM_out <- select(NDMM_out, -DOB, -MRD_sample, -aspirate_sample_comment)
NDMM_out <- filter(NDMM_out, !is.na(sample_ID))

write.xlsx(NDMM_out, 'p222_NDMM_withMRN.xlsx')
