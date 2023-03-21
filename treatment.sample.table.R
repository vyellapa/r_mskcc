####Generating data frame with treatment information from the time of sampling, RRMM cohort

#Load data: 
#setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical')
#source('clinical.data.import.R')

#Run: source('treatment.sample.table.R')

###Notes:
#Percentages do not count missing values. Percentages of resistance are estimated from the population who were exposed to the drug.
#Define treatment line in legend.
#Define resistance in legend.

#Calculating no of missing data for variables
missingOnTreatment <- sum(is.na(df$onTreatmentSample))
missingIMID <- sum(is.na(df$IMIDSampl))
missingThal <- sum(is.na(df$thalSampl))
missingLen <- sum(is.na(df$lenSampl))
missingPom <- sum(is.na(df$pomSampl))
missingPI <- sum(is.na(df$PISampl))
missingBort <- sum(is.na(df$bortSampl))
missingCarf <- sum(is.na(df$carfSampl))
missingIxa <- sum(is.na(df$ixaSampl))
missingLastDrugs <- sum(is.na(df$lastDrugs))
missingFirstDrugs <- sum(is.na(df$firstDrugs))

#Vector for column varnames
varnames <- c('Treatment from diagnosis to sampling',
              'Prior lines of therapy, median (range)',
              'One or more IMIDs',
              'Thalidomide',
              'Lenalidomide',
              'Pomalidomide',
              'One or more PIs',
              'Bortezomib',
              'Carfilzomib',
              'Ixazomib',
              'HDM-ASCT',
              'alloSCT',
              'Resistance at any time before sampling\n (% based no. of exposed)',
              'Any IMID',
              'Thalidomide',
              'Lenalidomide',
              'Pomalidomide',
              'Any PI',
              'Bortezomib',
              'Carfilzomib',
              'Ixazomib',
              'Treatment status at time of sampling',
              'On or shortly following treatment',
              'Off treatment',
              'Last drugs received before sampling',
              'IMID maintenance',
              'HDT-ASCT',
              'VDT-PACE',
              'IMID and PI other than VDT-PACE',
              'IMID',
              'PI',
              'First drugs received after sampling',
              'IMID maintenance',
              'HDM-ASCT',
              'VDT-PACE',
              'IMID and PI other than VDT-PACE',
              'IMID',
              'PI',
              'Other treatment'
)

#Vector for column numbers
numbers <- c(NA,
             paste(median(df$treatmentsSampl, na.rm = TRUE), ' (', min(df$treatmentsSampl, na.rm = TRUE), '-', max(df$treatmentsSampl, na.rm = TRUE), ')', sep = ''),
             sum(as.numeric(df$IMIDSampl) == 2, na.rm = TRUE),
             sum(as.numeric(df$thalSampl) >= 2, na.rm = TRUE),
             sum(as.numeric(df$lenSampl) >= 2, na.rm = TRUE),
             sum(as.numeric(df$pomSampl) >= 2, na.rm = TRUE),
             sum(as.numeric(df$PISampl) == 2, na.rm = TRUE),
             sum(as.numeric(df$bortSampl) >= 2, na.rm = TRUE),
             sum(as.numeric(df$carfSampl) >= 2, na.rm = TRUE),
             sum(as.numeric(df$ixaSampl) >= 2, na.rm = TRUE),
             sum(df$ASCTSampl, na.rm = TRUE),
             sum(df$alloSCTSampl, na.rm = TRUE),
             NA,
             sum(as.numeric(df$IMIDSamplResistant) == 2, na.rm = TRUE),
             sum(as.numeric(df$thalSampl) == 3, na.rm = TRUE),
             sum(as.numeric(df$lenSampl) == 3, na.rm = TRUE),
             sum(as.numeric(df$pomSampl) == 3, na.rm = TRUE),
             sum(as.numeric(df$PISamplResistant) == 2, na.rm = TRUE),
             sum(as.numeric(df$bortSampl) == 3, na.rm = TRUE),
             sum(as.numeric(df$carfSampl) == 3, na.rm = TRUE),
             sum(as.numeric(df$ixaSampl) == 3, na.rm = TRUE),
             NA,
             sum(as.numeric(df$onTreatmentSample) == 2, na.rm = TRUE),
             sum(as.numeric(df$onTreatmentSample) == 1, na.rm = TRUE),
             NA,
             sum(as.numeric(df$lastDrugs) == 1, na.rm = TRUE),
             sum(as.numeric(df$lastDrugs) == 2, na.rm = TRUE),
             sum(as.numeric(df$lastDrugs) == 3, na.rm = TRUE),
             sum(as.numeric(df$lastDrugs) == 4, na.rm = TRUE),
             sum(as.numeric(df$lastDrugs) == 5, na.rm = TRUE),
             sum(as.numeric(df$lastDrugs) == 6, na.rm = TRUE),
             NA,
             sum(as.numeric(df$firstDrugs) == 1, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 2, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 3, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 4, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 5, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 6, na.rm = TRUE),
             sum(as.numeric(df$firstDrugs) == 7, na.rm = TRUE)
)

#Vector for column percentages
percentages <- c(NA,
                 NA,
                 round((sum(as.numeric(df$IMIDSampl) == 2, na.rm = TRUE)/(length(df$sampleID)-missingIMID))*100),
                 round((sum(as.numeric(df$thalSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingThal))*100),
                 round((sum(as.numeric(df$lenSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingLen))*100),
                 round((sum(as.numeric(df$pomSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingPom))*100),
                 round((sum(as.numeric(df$PISampl) == 2, na.rm = TRUE)/(length(df$sampleID)-missingPI))*100),
                 round((sum(as.numeric(df$bortSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingBort))*100),
                 round((sum(as.numeric(df$carfSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingCarf))*100),
                 round((sum(as.numeric(df$ixaSampl) >= 2, na.rm = TRUE)/(length(df$sampleID)-missingIxa))*100),
                 round((sum(df$ASCTSampl, na.rm = TRUE)/(length(df$sampleID)-sum(is.na(df$ASCTSampl))))*100),
                 round((sum(df$alloSCTSampl, na.rm = TRUE)/(length(df$sampleID)-sum(is.na(df$alloSCTSampl))))*100),
                 NA,
                 round((sum(as.numeric(df$IMIDSamplResistant) == 2, na.rm = TRUE)/(length(df$sampleID)-missingIMID))*100),
                 round((sum(as.numeric(df$thalSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$thalSampl) >= 2, na.rm = TRUE)))*100),
                 round((sum(as.numeric(df$lenSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$lenSampl) >= 2, na.rm = TRUE)))*100),
                 round((sum(as.numeric(df$pomSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$pomSampl) >= 2, na.rm = TRUE)))*100),
                 round((sum(as.numeric(df$PISamplResistant) == 2, na.rm = TRUE)/(length(df$sampleID)-missingPI))*100),
                 round((sum(as.numeric(df$bortSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$bortSampl) >= 2, na.rm = TRUE)))*100),
                 round((sum(as.numeric(df$carfSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$carfSampl) >= 2, na.rm = TRUE)))*100),
                 round((sum(as.numeric(df$ixaSampl) == 3, na.rm = TRUE)/(sum(as.numeric(df$ixaSampl) >= 2, na.rm = TRUE)))*100),
                 NA,
                 round((sum(as.numeric(df$onTreatmentSample) == 2, na.rm = TRUE)/(length(df$sampleID)-missingOnTreatment))*100),
                 round((sum(as.numeric(df$onTreatmentSample) == 1, na.rm = TRUE)/(length(df$sampleID)-missingOnTreatment))*100),
                 NA,
                 round((sum(as.numeric(df$lastDrugs) == 1, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 round((sum(as.numeric(df$lastDrugs) == 2, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 round((sum(as.numeric(df$lastDrugs) == 3, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 round((sum(as.numeric(df$lastDrugs) == 4, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 round((sum(as.numeric(df$lastDrugs) == 5, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 round((sum(as.numeric(df$lastDrugs) == 6, na.rm = TRUE)/(length(df$sampleID)-missingLastDrugs))*100),
                 NA,
                 round((sum(as.numeric(df$firstDrugs) == 1, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 2, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 3, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 4, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 5, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 6, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100),
                 round((sum(as.numeric(df$firstDrugs) == 7, na.rm = TRUE)/(length(df$sampleID)-missingFirstDrugs))*100)
)

#Vector for column missing data
missing <- c(NA,
             sum(is.na(df$treatmentsSampl)),
             missingIMID,
             missingThal,
             missingLen,
             missingPom,
             missingPI,
             missingBort,
             missingCarf,
             missingIxa,
             sum(is.na(df$ASCTSampl)),
             sum(is.na(df$alloSCTSampl)),
             NA,
             missingIMID,
             missingThal,
             missingLen,
             missingPom,
             missingPI,
             missingBort,
             missingCarf,
             missingIxa,
             missingOnTreatment,
             NA,
             NA,
             missingLastDrugs,
             NA,
             NA,
             NA,
             NA,
             NA,
             NA,
             missingFirstDrugs,
             NA,
             NA,
             NA,
             NA,
             NA,
             NA,
             NA

)

treatmentSampleTable <- data.frame(variable = varnames, value = numbers, percentage = percentages, missing = missing, stringsAsFactors = FALSE)

colnames(treatmentSampleTable) <- c('Variable', 'no./value', '%', 'no. missing')
