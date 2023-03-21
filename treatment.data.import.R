########
#R-Script to import treatment data for down stream analyses in R - RRMM seq project MSKCC. 
#Run: setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical'); source('treatment.data.import.R')

#Set working directory
setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical')

#Import required packages
library(openxlsx)
library(dplyr)
library(lubridate)

#Import treatment spreadsheet
treatment <- read.xlsx('../datafiles/treatment.data.xlsx', na.strings = c('ND', 'NA'), detectDates = TRUE)

####Assigning class and label to each variable
treatment$sampleID <- factor(treatment$sampleID)
treatment$treatnr <- as.integer(treatment$treatnr)
treatment$start.date <- as.POSIXct(strptime(treatment$start.date, format = '%Y-%m-%d'))
treatment$stop.date <- as.POSIXct(strptime(treatment$stop.date, format = '%Y-%m-%d'))
treatment$start.context <- factor(treatment$start.context)
treatment$drugs <- factor(treatment$drugs)
treatment$response <- factor(treatment$response)
treatment$stop.reason <- factor(treatment$stop.reason)


#compute variables
#Need to make a variable with the number of days since diagnosis at the start and end of each treatment line.

#####Figure script
library(ggplot2)


