####Generating data frame with descriptives from RRMM cohort

#Run: setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical'); source('clinical.data.import.R')

###Notes:
#Like Samplline.
#

#Calculating no of missing data for variables
missingContext <- sum(is.na(df$contextSample))
missingLytic <- sum(is.na(df$lyticLesionSampl))
missingPET <- sum(is.na(df$PETPositiveSampl))
missingFocal <- sum(is.na(df$focalLesionsSampl))
missingEMD <- sum(is.na(df$EMDSampl))
missingISS <- sum(is.na(df$ISSSampl))

#Vector for column varnames
varnames <- c('Time from diagnosis to sampling, years, median (IQR)',
              'Disease status at time of sampling',
              'Relapse',
              'Re-staging during treatment line',
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
numbers <- c(paste(round(median(df$timeDiagnSample/12, na.rm = TRUE), digits = 1), ' (', round(IQR(df$timeDiagnSample/12, na.rm = TRUE), digits = 1), ')', sep = ''),
             NA,
             sum(as.numeric(df$contextSample) == 2, na.rm = TRUE),
             sum(as.numeric(df$contextSample) == 1, na.rm = TRUE),
             NA,
             paste(round(median(df$BMPCSampl, na.rm = TRUE)), ' (', round(IQR(df$BMPCSampl, na.rm = TRUE)), ')', sep = ''),
             paste(round(median(df$hgbSampl, na.rm = TRUE), digits = 1), ' (', round(IQR(df$hgbSampl, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$caSampl, na.rm = TRUE), digits = 1), ' (', round(IQR(df$caSampl, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$creaSampl, na.rm = TRUE), digits = 1), ' (', round(IQR(df$creaSampl, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$ANCSampl, na.rm = TRUE), digits = 1), ' (', round(IQR(df$ANCSampl, na.rm = TRUE), digits = 1), ')', sep = ''),
             paste(round(median(df$plateletsSampl, na.rm = TRUE), digits = 1), ' (', round(IQR(df$plateletsSampl, na.rm = TRUE), digits = 1), ')', sep = ''),
             NA,
             sum(df$lyticLesionSampl, na.rm = TRUE),
             sum(df$PETPositiveSampl, na.rm = TRUE),
             sum(df$focalLesionsSampl >= 1, na.rm = TRUE),
             sum(df$EMDSampl, na.rm = TRUE),
             NA,
             sum(df$ISSSampl == 1, na.rm = TRUE),
             sum(df$ISSSampl == 2, na.rm = TRUE),
             sum(df$ISSSampl == 3, na.rm = TRUE),
             NA
)

#Vector for column percentages
percentages <- c(NA,
                 NA,
                 round((sum(as.numeric(df$contextSample) == 2, na.rm = TRUE)/(length(df$sampleID)-missingContext))*100),
                 round((sum(as.numeric(df$contextSample) == 1, na.rm = TRUE)/(length(df$sampleID)-missingContext))*100),
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 round((sum(df$lyticLesionSampl, na.rm = TRUE)/(length(df$sampleID)-missingLytic))*100),
                 round((sum(df$PETPositiveSampl, na.rm = TRUE)/(length(df$sampleID)-missingPET))*100),
                 round((sum(df$focalLesionsSampl >= 1, na.rm = TRUE)/(length(df$sampleID)-missingFocal))*100),
                 round((sum(df$EMDSampl, na.rm = TRUE)/(length(df$sampleID)-missingEMD))*100),
                 NA,
                 round((sum(df$ISSSampl == 1, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 round((sum(df$ISSSampl == 2, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 round((sum(df$ISSSampl == 3, na.rm = TRUE)/(length(df$sampleID)-missingISS))*100),
                 NA
)

#Vector for column missing data
missing <- c(sum(is.na(df$timeDiagnSample)),
             missingContext,
             NA,
             NA,
             NA,
             sum(is.na(df$BMPCSampl)),
             sum(is.na(df$hgbSampl)),
             sum(is.na(df$caSampl)),
             sum(is.na(df$creaSampl)),
             sum(is.na(df$ANCSampl)),
             sum(is.na(df$plateletsSampl)),
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

descriptiveSampleTable <- data.frame(variable = varnames, value = numbers, percentage = percentages, missing = missing, stringsAsFactors = FALSE)

colnames(descriptiveSampleTable) <- c('Variable', 'no./value', '%', 'no. missing')