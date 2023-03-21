# Script to combine CNV calls for p222
# Specific to the agreed column order and file format 
# Version 3.8.18

setwd('/Users/rustade/Documents/MSKCC/Projects/222/All_myType/CNV')
library(openxlsx)
library(plyr)
library(dplyr)

# Importe all CNV files
cnames <- c('SAMPLE_LEUKID', "Del_1p", "Gain_1q", "Del_6q", "Del_8p", "Gain_11q", "Del_12p", "Del_13q", "Del_16q", "Del_17p", "HRD", "Gender", 'Comment')

ER <- read.xlsx('p222_CNV_annotated_ER.xlsx'); ER <- ER[,-c(2,3)]
names(ER) <- cnames; names(ER)[-1] <- paste(names(ER)[-1], 'ER', sep = '_')
ER$Gender_ER <- gsub('FEMALE', 'F', ER$Gender_ER, fixed = T); ER$Gender_ER <- gsub('MALE', 'M', ER$Gender_ER, fixed = T)

Teja <- read.table('p222_cna_annotations_Teja.txt', sep = '\t', header = T); Teja$Comment <- NA
names(Teja) <- cnames; names(Teja)[-1] <- paste(names(Teja)[-1], 'Teja', sep = '_')
Teja[,-1] <- mutate_all(Teja[,-1], .funs=toupper); Teja$HRD_Teja <- gsub('HYPERDIPLOIDY', 'YES', Teja$HRD_Teja, fixed = T)

TA <- read.xlsx('CNV_project222_TA.xlsx'); TA <- TA[,-c(2,3,16)]
names(TA) <- cnames; names(TA)[-1] <- paste(names(TA)[-1], 'TA', sep = '_')

MH <- read.xlsx('CNV_project222_MH030818.xlsx'); MH <- MH[,-c(2,3)]
names(MH) <- c('SAMPLE_LEUKID', "Del_1p", "Gain_1q", "Gain_11q", "Del_12p", "Del_13q", "Del_16q", "Del_17p", "HRD", "Gender", 'Comment')
names(MH)[-1] <- paste(names(MH)[-1], 'MH', sep = '_')
MH$Gender_MH <- gsub('WOMAN', 'F', MH$Gender_MH, fixed = T); MH$Gender_MH <- gsub('MAN', 'M', MH$Gender_MH, fixed = T)

# Merge
all <- join(Teja, ER, by = 'SAMPLE_LEUKID'); all <- join(all, TA, by = 'SAMPLE_LEUKID'); all <- join(all, MH, by = 'SAMPLE_LEUKID')

# Sort variables
all$CONSENSUS <- NA
all$AGREE <- 'AGREE'
all$CHECK <- NA
all <- all[, order(names(all), decreasing = T)]

cnames <- c("SAMPLE_LEUKID", "CONSENSUS", "CHECK", "AGREE", "HRD_Teja", "HRD_TA", "HRD_MH", "HRD_ER", 
            "Gender_Teja", "Gender_TA", "Gender_MH", "Gender_ER", "Gain_1q_Teja", 
            "Gain_1q_TA", "Gain_1q_MH", "Gain_1q_ER", "Gain_11q_Teja", "Gain_11q_TA", 
            "Gain_11q_MH", "Gain_11q_ER", "Del_8p_Teja", "Del_8p_TA", "Del_8p_ER", 
            "Del_6q_Teja", "Del_6q_TA", "Del_6q_ER", "Del_1p_Teja", "Del_1p_TA", 
            "Del_1p_MH", "Del_1p_ER", "Del_17p_Teja", "Del_17p_TA", "Del_17p_MH", 
            "Del_17p_ER", "Del_16q_Teja", "Del_16q_TA", "Del_16q_MH", "Del_16q_ER", 
            "Del_13q_Teja", "Del_13q_TA", "Del_13q_MH", "Del_13q_ER", "Del_12p_Teja", 
            "Del_12p_TA", "Del_12p_MH", "Del_12p_ER", "Comment_Teja", 
            "Comment_TA", "Comment_MH", "Comment_ER")

all <- all[, cnames]


for(i in 1:nrow(all)){
        if((length(unique(as.character(all[i,grep('Del_1p', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_1p'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_1p')
                }}
        if((length(unique(as.character(all[i,grep('Gain_1q', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Gain_1q'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Gain_1q')
                }}
        if((length(unique(as.character(all[i,grep('Del_6q', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_6q'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_6q')
                }}
        if((length(unique(as.character(all[i,grep('Del_8p', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_8p'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_8p')
                }}
        if((length(unique(as.character(all[i,grep('Gain_11q', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Gain_11q'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Gain_11q')
                }}
        if((length(unique(as.character(all[i,grep('Del_12p', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_12p'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_12p')
                }}
        if((length(unique(as.character(all[i,grep('Del_13q', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_13q'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_13q')
                }}
        if((length(unique(as.character(all[i,grep('Del_16q', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_16q'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_16q')
                }}
        if((length(unique(as.character(all[i,grep('Del_17p', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Del_17p'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Del_17p')
                }}
        if((length(unique(as.character(all[i,grep('HRD', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'HRD'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'HRD')
                }}
        if((length(unique(as.character(all[i,grep('Gender', names(all), value = T)])))) > 1){
                all$AGREE[i] <- 'CHECK'
                if(is.na(all$CHECK[i])){
                        all$CHECK[i] <- 'Gender'
                } else {
                        all$CHECK[i] <- paste(all$CHECK[i], 'Gender')
                }}
}

# Export
write.xlsx(all, 'p222_CNV_merged.xlsx', asTable = F)
