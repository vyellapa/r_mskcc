########
#R-Script to import and clean clinical dataset for down stream analyses in R - RRMM seq project MSKCC. 
#Run: setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.clinical'); source('clinical.data.import.R')

#Set working directory
setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.clinical')


#Import required packages
library(openxlsx)
library(dplyr)
library(lubridate)

#Import clinical data excel spreadsheet 

df <- read.xlsx('../datafiles/RRMM_HOTB_Manual.xlsx', na.strings = c('ND', 'NA'), detectDates = TRUE)
df$sampleID <- as.factor(df$sampleID)

####Assigning class and label to each variable
#Charcter and numeric variables are imported in correct format and are thus not changed. 

#General patient information
df$recordReviewDate <- as.POSIXct(strptime(df$recordReviewDate, format = '%Y-%m-%d'))
df$Sex <- factor(df$Sex, levels = c(1:2),labels = c('Male', 'Female'))
df$DOB <- as.POSIXct(strptime(df$DOB, format = '%Y-%m-%d'))
df$diagnDateMM <- as.POSIXct(strptime(df$diagnDateMM, format = '%Y-%m-%d'))
df$sampleDate <- as.POSIXct(strptime(df$sampleDate, format = '%Y-%m-%d'))
df$deceasedDate <- as.POSIXct(strptime(df$deceasedDate, format = '%Y-%m-%d'))
df$lastVisit <- as.POSIXct(strptime(df$lastVisit, format = '%Y-%m-%d'))

#Disease status
df$monitoring <- factor(df$monitoring, levels = c(0:4), labels = c('Non-secretory', 'Oligo-secretory', 'Urine M-spike', 'Serum FLC', 'Serum M-spike'))
df$IFEHeavyBase <- factor(df$IFEHeavyBase, levels = c(0:3), labels = c('Negative', 'IgG', 'IgA', 'IgD'))
df$IFELightBase <- factor(df$IFELightBase, levels = c(0:2), labels = c('Negative', 'Lambda', 'Kappa'))
df$IFEUrineBase <- factor(df$IFEUrineBase, levels = c(0:3), labels = c('Negative', 'Lambda', 'Kappa', 'Lambda and Kappa'))
df$LDHBase <- factor(df$LDHBase, levels = c(0:1),labels = c('Normal', 'Elevated'))
df$lyticLesionBase <- as.logical(df$lyticLesionBase)
df$pathFractureBase <- as.logical(df$pathFractureBase)
df$PETPositiveBase <- as.logical(df$PETPositiveBase)
df$EMDBase <- as.logical(df$EMDBase)
df$CRABCaBase <- as.logical(df$CRABCaBase)
df$CRABRenalBase <- as.logical(df$CRABRenalBase)
df$CRABAnemiaBase <- as.logical(df$CRABAnemiaBase)
df$CRABBoneBase <- as.logical(df$CRABBoneBase)
df$CRABSubBase <- as.logical(df$CRABSubBase)
df$HRBiomarkerBMPC <- as.logical(df$HRBiomarkerBMPC)
df$HRBiomarkerFLC <- as.logical(df$HRBiomarkerFLC)
df$HRBiomarkerLesion <- as.logical(df$HRBiomarkerLesion)

df$IFEHeavySampl <- factor(df$IFEHeavySampl, levels = c(0:3), labels = c('Negative', 'IgG', 'IgA', 'IgD'))
df$IFELightSampl <- factor(df$IFELightSampl, levels = c(0:2), labels = c('Negative', 'Lambda', 'Kappa'))
df$IFEUrineSampl <- factor(df$IFEUrineSampl, levels = c(0:3), labels = c('Negative', 'Lambda', 'Kappa', 'Lambda and Kappa'))
df$LDHSampl <- factor(df$LDHSampl, levels = c(0:1),labels = c('Normal', 'Elevated'))
df$lyticLesionSampl <- as.logical(df$lyticLesionSampl)
df$pathFractureSampl <- as.logical(df$pathFractureSampl)
df$PETPositiveSampl <- as.logical(df$PETPositiveSampl)
df$EMDSampl <- as.logical(df$EMDSampl)

df$tempIFEBase <- addNA(df$IFELightBase)
levels(df$tempIFEBase) <- c(levels(df$tempIFEBase),999)
df$tempIFEBase[is.na(df$tempIFEBase)] <- 999

df <- mutate(df, involvedFLCBase = ifelse(df$tempIFEBase== 'Lambda', df$lambdaBase, 
                                          ifelse(df$tempIFEBase == 'Kappa', df$kappaBase, 
                                                 ifelse(df$lambdaBase>=10, df$lambdaBase,
                                                        ifelse(df$kappaBase>=10, df$kappaBase, NA)))))

df$tempIFESampl <- addNA(df$IFELightSampl)
levels(df$tempIFESampl) <- c(levels(df$tempIFESampl),999)
df$tempIFESampl[is.na(df$tempIFESampl)] <- 999

df <- mutate(df, involvedFLCSampl = ifelse(df$tempIFESampl== 'Lambda', df$lambdaSampl, 
                                          ifelse(df$tempIFESampl == 'Kappa', df$kappaSampl, 
                                                 ifelse(df$lambdaSampl>=10, df$lambdaSampl,
                                                        ifelse(df$kappaSampl>=10, df$kappaSampl, NA)))))

#Treatment and resistance data
df$ASCTSampl <- as.logical(df$ASCTSampl)
df$alloSCTSampl <- as.logical(df$alloSCTSampl)
df$thalSampl <- factor(df$thalSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$lenSampl <- factor(df$lenSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$pomSampl <- factor(df$pomSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$bortSampl <- factor(df$bortSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$carfSampl <- factor(df$carfSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$ixaSampl <- factor(df$ixaSampl, levels = c(0:3), labels = c('Not exposed', 'Not resistant', 'Resistant', 'Exposed, unevaluable'))
df$lastTreatmentSamplDate <- as.POSIXct(strptime(df$lastTreatmentSamplDate, format = '%Y-%m-%d'))
df$lastTreatmentSamplRegimen <- as.factor(df$lastTreatmentSamplRegimen)
df$lastTreatmentSamplResponse <- factor(df$lastTreatmentSamplResponse, levels = c(1:6), labels = c('CR', 'VGPR', 'PR', 'MR', 'SD', 'PD'))
df$contextSample <- factor(df$contextSample, levels = c(0:1), labels = c('Not relapse', 'Relapse'))
df$treatmentPostSampleDate <- as.POSIXct(strptime(df$treatmentPostSampleDate, format = '%Y-%m-%d'))
df$treatmentPostSampleRegimen <- as.factor(df$treatmentPostSampleRegimen)
df$treatmentPostSampleResponse <- factor(df$treatmentPostSampleResponse, levels = c(0:6), labels = c('Unevaluable', 'CR', 'VGPR', 'PR', 'MR', 'SD', 'PD'))
df$secondTreatmentPostSampleDate <- as.POSIXct(strptime(df$secondTreatmentPostSampleDate, format = '%Y-%m-%d'))
df$ASCTAfter <- as.logical(df$ASCTAfter)
df$alloSCTAfter <- as.logical(df$alloSCTAfter)
df$lastDrugs <- factor(df$lastDrugs, levels = c(1:6), labels = c('IMID maintenance', 'HDT-ASCT', 'VDT-PACE', 'IMID and PI', 'IMID', 'PI'))
df$onTreatmentSample <- factor(df$onTreatmentSample, levels = c(0:1), labels = c('Off treatment', 'On or shortly following treatment'))
df$firstDrugs <- factor(df$firstDrugs, levels = c(1:7), labels = c('IMID maintenance', 'HDT-ASCT', 'VDT-PACE', 'IMID and PI', 'IMID', 'PI', 'Other treatment'))
df$firstDrugsResponse <- factor(df$firstDrugsResponse, levels = c(0:6), labels = c('Unevaluable', 'CR', 'VGPR', 'PR', 'MR', 'SD', 'PD'))

####Merged variables
##ISS
df <- mutate(df, ISSBase = factor(ifelse(beta2mBase >= 5.5, 3, 
                                         ifelse(beta2mBase < 3.5 & albBase >= 3.5, 1, 2)), levels = c(1:3)))
df <- mutate(df, ISSSampl = factor(ifelse(beta2mSampl >= 5.5, 3, 
                                          ifelse(beta2mSampl < 3.5 & albSampl >= 3.5, 1, 2)), levels = c(1:3)))

#BMPC
df <- mutate(df, BMPCBase = ifelse(is.na(aspirateBase), biopsyBase, 
                                   ifelse(is.na(biopsyBase), aspirateBase, 
                                          ifelse(biopsyBase>aspirateBase, biopsyBase, aspirateBase))))
df <- mutate(df, BMPCSampl = ifelse(is.na(aspirateSampl), biopsySampl, 
                                    ifelse(is.na(biopsySampl), aspirateSampl, 
                                           ifelse(biopsySampl>aspirateSampl, biopsySampl, aspirateSampl))))

#Treatment and response categories
df <- mutate(df, treatmentsTotal = treatmentsSampl + treatmentsAfter)

df <- mutate(df, progressingSampl = factor(
        ifelse(as.numeric(onTreatmentSample) == 1 | as.numeric(contextSample) == 1, 6,
               ifelse(as.numeric(lastDrugs) == 5, 1, 
                      ifelse(as.numeric(lastDrugs) == 1, 2,
                             ifelse(as.numeric(lastDrugs) == 6, 3,
                                    ifelse(as.numeric(lastDrugs) == 4, 4,
                                           ifelse(as.numeric(lastDrugs) == 3, 5, NA)))))), 
        levels = c(1:6), labels = c('Relapse on IMID', 'Relapse on IMID maint.', 'Relapse on PI', 'Relapse on PI and IMID', 'Relapse on VDT-PACE', 'Not relapse or not on treatment')))
#Change 'progressingSampl' to 'refractorySampl' when patients are re-classified as appropriate.

df <- mutate(df, refractoryFirstDrugs = factor(
        ifelse(as.numeric(firstDrugsResponse) <=5, 6,
               ifelse(as.numeric(firstDrugs) %in% c(1,5), 1,
                      ifelse(as.numeric(firstDrugs) == 6, 2,
                             ifelse(as.numeric(firstDrugs) == 4, 3,
                                    ifelse(as.numeric(firstDrugs) %in% c(2:3), 4, 5))))), 
        levels = c(1:6), labels = c('IMID refractory', 'PI refractory', 'IMID and PI refractory', 'Intensive chemo refractory', 'Refractory other', 'Not refractory')))
#Intensive chemo = VDT-PACE or HDM-ASCT. Other= two patients on IbrutinibDex trial, one on AT-7519.
#Definition of refractory = primary SD or PD

df <- mutate(df, PDFirstDrugs = factor(
        ifelse(as.numeric(firstDrugsResponse) <=6, 6,
               ifelse(as.numeric(firstDrugs) %in% c(1,5), 1,
                      ifelse(as.numeric(firstDrugs) == 6, 2,
                             ifelse(as.numeric(firstDrugs) == 4, 3,
                                    ifelse(as.numeric(firstDrugs) %in% c(2:3), 4, 5))))), 
        levels = c(1:6), labels = c('IMID refractory', 'PI refractory', 'IMID and PI refractory', 'Intensive chemo refractory', 'Refractory other', 'Not refractory')))
#Definition of refractory = primary PD

df <- mutate(df, refractoryCombined = factor(
        ifelse(as.numeric(progressingSampl) %in% c(1,2) | as.numeric(refractoryFirstDrugs) == 1, 1,
               ifelse(as.numeric(progressingSampl) == 3 | as.numeric(refractoryFirstDrugs) == 2, 2,
                      ifelse(as.numeric(progressingSampl) == 4 | as.numeric(refractoryFirstDrugs) == 3, 3,
                             ifelse(as.numeric(progressingSampl) == 5 | as.numeric(refractoryFirstDrugs) == 4, 4, 5)))), 
        levels = c(1:5), labels = c('IMID refractory', 'PI refractory', 'IMID and PI refractory', 'Intensive chemo refractory', 'Not refractory')))

df <- mutate(df, IMIDSampl = factor(
        ifelse(as.numeric(thalSampl) >= 2, 2, 
               ifelse(as.numeric(lenSampl) >= 2, 2, 
                      ifelse(as.numeric(pomSampl) >= 2, 2, 1))), levels = c(1, 2), labels = c('Not Exposed', 'Exposed')))

df <- mutate(df, IMIDSamplResistant = factor(
        ifelse(as.numeric(thalSampl) == 3, 2, 
               ifelse(as.numeric(lenSampl) == 3, 2, 
                      ifelse(as.numeric(pomSampl) == 3, 2, 1))), levels = c(1, 2), labels = c('Not Resistant', 'Resistant')))

df <- mutate(df, PISampl = factor(
        ifelse(as.numeric(bortSampl) >= 2, 2, 
               ifelse(as.numeric(carfSampl) >= 2, 2, 
                      ifelse(as.numeric(ixaSampl) >= 2, 2, 1))), levels = c(1, 2), labels = c('Not Exposed', 'Exposed')))
df <- mutate(df, PISamplResistant = factor(
        ifelse(as.numeric(bortSampl) == 3, 2, 
               ifelse(as.numeric(carfSampl) == 3, 2, 
                      ifelse(as.numeric(ixaSampl) == 3, 2, 1))), levels = c(1, 2), labels = c('Not Resistant', 'Resistant')))

#Time and survival variables
df <- mutate(df, diagnYearMM = as.POSIXlt(diagnDateMM)$year + 1900)
df$diagnAge  <-  year(as.period(interval(df$DOB, df$diagnDateMM)))
df <- mutate(df, statusOS = factor(
        ifelse(is.na(deceasedDate), 1, 2), levels = c(1,2), labels = c('Alive', 'Dead')))

df <- mutate(df, dateOS = deceasedDate)
for(i in 1:length(df$dateOS)){
        if(is.na(df$dateOS[i])){df$dateOS[i] <- df$recordReviewDate[i]} 
        else{}}

df$timeOS <- month(as.period(interval(df$diagnDateMM, df$dateOS), unit = 'month'))
df$timeSampleOS <- month(as.period(interval(df$sampleDate, df$dateOS), unit = 'month'))
df$timeDiagnSample <- month(as.period(interval(df$diagnDateMM, df$sampleDate), unit = 'month'))
df$timeSampleTreatment <- month(as.period(interval(df$sampleDate, df$treatmentPostSampleDate), unit = 'month'))

df <- mutate(df, statusSamplePFS = factor(
        ifelse(!is.na(secondTreatmentPostSampleDate)|!is.na(deceasedDate), 2, 1), 
        levels = c(1,2), labels = c('Not progressed', 'Progressed')))

df <-  mutate(df, dateSamplePFS = secondTreatmentPostSampleDate)
for(i in 1:length(df$dateSamplePFS)){
        if(is.na(df$dateSamplePFS[i])){
                if(is.na(df$deceasedDate[i])){
                        df$dateSamplePFS[i] <- df$lastVisit[i]}
                else{df$dateSamplePFS[i] <- df$deceasedDate[i]
                }}else{}}
rm(i)

df$timeSamplePFS <- month(as.period(interval(df$treatmentPostSampleDate, df$dateSamplePFS), unit = 'month'))
