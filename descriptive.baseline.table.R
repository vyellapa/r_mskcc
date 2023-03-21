####Generating data frame with descriptives from RRMM cohort

#Run: setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical'); source('descriptive.baseline.table.R')

###Notes:
#Missing values do not count towards percentages.
#Lytic lesions variable defined identically to bone disease in CRAB.
#Missing R-ISS, for that I need to import cytogenetics.
#Include normal ranges. 
#BMPC is the highest of aspirate and biopsy

#Calculating no of missing data for variables
missingHeavyBase <- length(df$sampleID)-sum(c(sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgG', na.rm = TRUE),
                                              sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgA', na.rm = TRUE),
                                              sum(as.numeric(df$monitoring) %in% c(3,4), na.rm = TRUE),
                                              sum(as.numeric(df$monitoring) <=2, na.rm = TRUE)))
missingLightBase <- length(df$sampleID)-sum(c(sum(df$IFELightBase=='Kappa', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$kappaBase>=100, na.rm = TRUE),
                                              sum(df$IFELightBase=='Lambda', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$lambdaBase>=100, na.rm = TRUE),
                                              sum(as.numeric(df$monitoring) == 1, na.rm = TRUE)))
missingLytic <- sum(is.na(df$lyticLesionBase))
missingPET <- sum(is.na(df$PETPositiveBase))
missingFocal <- sum(is.na(df$focalLesionsBase))
missingEMD <- sum(is.na(df$EMDBase))
missingISS <- sum(is.na(df$ISSBase))

#Vector for column varnames
varnames <- c('Number of patients',
              'Sex, male',
              'Age at diagnosis, years, median (IQR)',
              'Heavy chain',
              'IgG',
              'IgA',
              'Light chain disease',
              'Non/oligo-secretory disease',
              'Light chain',
              'Kappa',
              'Lambda',
              'Non-secretory disease',
              'Laboratory values',
              'Bone marrow plasma cell percentage, median (IQR)',
              'Hemoglobin, g/dL, median (IQR)',
              'P-Calcium total, mg/dL, median (IQR)',
              'S-Creatinine, mg/dL, median (IQR)',
              'Neutrophil count, K/mcl, median (IQR)',
              'Platelet count, K/mcl, median (IQR)',
              'Imaging',
              'Lytic bone disease',
              'PET positive',
              'Focal lesions (>=1)',
              'Extra medullary disease',
              'International Staging System (ISS)',
              'I',
              'II',
              'III',
              'Revised ISS (R-ISS)'
              )

#Vector for column numbers
numbers <- c(length(df$sampleID),
             sum(df$Sex=='Male'),
             paste(median(df$diagnAge, na.rm = TRUE), ' (', IQR(df$diagnAge, na.rm = TRUE), ')', sep = ''),
             NA,
             sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgG', na.rm = TRUE),
             sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgA', na.rm = TRUE),
             sum(as.numeric(df$monitoring) %in% c(3,4), na.rm = TRUE),
             sum(as.numeric(df$monitoring) <=2, na.rm = TRUE),
             NA,
             sum(df$IFELightBase=='Kappa', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$kappaBase>=100, na.rm = TRUE),
             sum(df$IFELightBase=='Lambda', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$lambdaBase>=100, na.rm = TRUE),
             sum(as.numeric(df$monitoring)==1, na.rm = TRUE),
             NA,
             paste(round(median(df$BMPCBase, na.rm = TRUE)), ' (', round(IQR(df$BMPCBase, na.rm = TRUE)), ')', sep = ''),
             paste(round(median(df$hgbBase, na.rm = TRUE), digits = 1), ' (', round(IQR(df$hgbBase, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$caBase, na.rm = TRUE), digits = 1), ' (', round(IQR(df$caBase, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$creaBase, na.rm = TRUE), digits = 1), ' (', round(IQR(df$creaBase, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$ANCBase, na.rm = TRUE), digits = 1), ' (', round(IQR(df$ANCBase, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$plateletsBase, na.rm = TRUE), digits = 1), ' (', round(IQR(df$plateletsBase, na.rm = TRUE), digits = 1), ')', sep = ''),
             NA,
             sum(df$lyticLesionBase, na.rm = TRUE),
             sum(df$PETPositiveBase, na.rm = TRUE),
             sum(df$focalLesionsBase >= 1, na.rm = TRUE),
             sum(df$EMDBase, na.rm = TRUE),
             NA,
             sum(df$ISSBase == 1, na.rm = TRUE),
             sum(df$ISSBase == 2, na.rm = TRUE),
             sum(df$ISSBase == 3, na.rm = TRUE),
             NA
             )

#Vector for column percentages
percentages <- c(NA,
                 round((sum(df$Sex=='Male')/sum(!is.na(df$Sex)))*100),
                 NA,
                 NA,
                 round((sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgG', na.rm = TRUE)/(length(df$sampleID)-missingHeavyBase))*100),
                 round((sum(filter(df, as.numeric(monitoring)==5)$IFEHeavyBase=='IgA', na.rm = TRUE)/(length(df$sampleID)-missingHeavyBase))*100),
                 round((sum(as.numeric(df$monitoring) %in% c(3,4), na.rm = TRUE)/(length(df$sampleID)-missingHeavyBase))*100),
                 round((sum(as.numeric(df$monitoring) <=2, na.rm = TRUE)/(length(df$sampleID)-missingHeavyBase))*100),
                 NA,
                 round(((sum(df$IFELightBase=='Kappa', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$kappaBase>=100, na.rm = TRUE))/(length(df$sampleID)-missingLightBase))*100),
                 round(((sum(df$IFELightBase=='Lambda', na.rm = TRUE)+sum(filter(df, IFELightBase == 'Negative')$lambdaBase>=100, na.rm = TRUE))/(length(df$sampleID)-missingLightBase))*100),
                 round((sum(as.numeric(df$monitoring)==1, na.rm = TRUE)/(length(df$sampleID)-missingLightBase))*100),
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 round((sum(df$lyticLesionBase, na.rm = TRUE)/(length(df$sampleID)-missingLytic))*100),
                 round((sum(df$PETPositiveBase, na.rm = TRUE)/(length(df$sampleID)-missingPET))*100),
                 round((sum(df$focalLesionsBase >= 1, na.rm = TRUE)/(length(df$sampleID)-missingFocal))*100),
                 round((sum(df$EMDBase, na.rm = TRUE)/(length(df$sampleID)-missingEMD))*100),
                 NA,
                 round((sum(df$ISSBase == 1, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 round((sum(df$ISSBase == 2, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 round((sum(df$ISSBase == 3, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 NA
                 )

#Vector for column missing data
missing <- c(NA,
             NA,
             NA,
             missingHeavyBase,
             NA,
             NA,
             NA,
             NA,
             missingLightBase,
             NA,
             NA,
             NA,
             NA,
             sum(is.na(df$BMPCBase)),
             sum(is.na(df$hgbBase)),
             sum(is.na(df$caBase)),
             sum(is.na(df$creaBase)),
             sum(is.na(df$ANCBase)),
             sum(is.na(df$plateletsBase)),
             NA,
             missingLytic,
             missingPET,
             missingFocal,
             missingEMD,
             missingISS,
             NA,
             NA,
             NA,
             NA
             
)

descriptiveBaselineTable <- data.frame(variable = varnames, value = numbers, percentage = percentages, missing = missing, stringsAsFactors = FALSE)

colnames(descriptiveBaselineTable) <- c('Variable', 'no./value', '%', 'no. missing')