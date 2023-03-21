# Dependencies
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(stringi)))

# Import signatures function definition
main_combine <- function(clinical_path, survival_path, sv_path, cnv_path, tpm_path){
    ## Import data and combine key variables to form main analysis dataset
    ## Generate summary variables
    
    ## Clinical data - legend
    # ftrttrpl	Triplet as first treatment regimen flag
    # D_PT_therclass	Therapy classification
    # line1sct	SCT during first line of therapy flag
    # mainttrt	Treatment following first SCT
    # maintday	Duration of post-SCT maintenance treatment
    clinical <- read.csv(clinical_path, stringsAsFactors = F)
    clinical <- select(clinical, PUBLIC_ID, D_PT_gender, D_PT_age, D_PT_iss, D_PT_therclass, ftrttrpl, line1sct, mainttrt, maintday)
    names(clinical)[1] <- 'public_id'
    
    clinical$D_PT_therclass <- sub('combined ', '', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('bortezomib/carfilzomib-based', 'Bort/Car-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('bortezomib/IMIDs-based', 'Bort/IMID-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('bortezomib/IMIDs/carfilzomib-based', 'Bort/Car/IMID-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('IMIDs/carfilzomib-based', 'Car/IMID-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('Bortezomib-based', 'Bort-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('Carfilzomib-based', 'Car-based', clinical$D_PT_therclass)
    clinical$D_PT_therclass <- sub('IMIDs-based', 'IMID-based', clinical$D_PT_therclass)
    
    clinical <- mutate(clinical, Maintenance = ifelse(mainttrt == "", 0, 1))
    names(clinical)[which(names(clinical) == 'ftrttrpl')] <- 'Triplet_induction'
    names(clinical)[which(names(clinical) == 'line1sct')] <- 'ASCT_upfront'
    names(clinical)[which(names(clinical) == 'D_PT_iss')] <- 'ISS'
    
    
    ## Survival - legend
    # deathdy	Date of Death
    # lstalive	Last known alive
    # ttpfs	Time to PFS
    # pfscdy	PFS censored date
    # censpfs	Censor flag: progression-free survival

    survival <- read.csv(survival_path, stringsAsFactors = F)
    survival <- select(survival, public_id, deathdy, lstalive, maxline, ttpfs, censpfs, pfscdy)
    
    survival <- mutate(survival, 
                       timepfs = ifelse(is.na(ttpfs), pfscdy, ttpfs),
                       timeos = ifelse(is.na(deathdy), lstalive, deathdy),
                       censos = ifelse(is.na(deathdy), 0, 1))
    survival <- mutate(survival,
                       censos = ifelse(is.na(timeos), NA, censos))
    
    ## SVs - legend
    #t(11;14) -- CCND1
    #t(6;14) -- CCND3
    #t(4;14) -- WHSC1/MMSET
    #t(14;16) -- MAF (APOBEC associated)
    #t(14;20) -- MAFB (APOBEC associated)
    SV <- read.table(sv_path, sep = '\t', header = T, stringsAsFactors = F)
    SV <- SV[grep('1_BM', SV$Study_Visit_iD),]
    SV <- SV[c(1,grep('_CALL$', names(SV)))]
    names(SV)[1] <- 'sample'
    names(SV)[which(names(SV) == 'SeqWGS_WHSC1_CALL')] <- 'MMSET'
    names(SV)[which(names(SV) == 'SeqWGS_MYC_CALL')] <- 'MYC'
    names(SV)[which(names(SV) == 'SeqWGS_CCND1_CALL')] <- 'CCND1'
    SV <- mutate(SV, public_id = stri_sub(sample, 0,9))
    SV <- mutate(SV, 
                 MAF_MAFA_MAFB = factor(ifelse((SeqWGS_MAF_CALL == 1) | (SeqWGS_MAFA_CALL == 1) | (SeqWGS_MAFB_CALL == 1), 1, 0)))
    
    # CNVs - legend
    CNV <- read.table(cnv_path, sep = '\t', header = T, stringsAsFactors = F)
    CNV <- CNV[grep('1_BM', CNV$Study_Visit_ID),]
    names(CNV)[1] <- 'sample'
    CNV <- mutate(CNV, public_id = stri_sub(sample, 0,9))
    names(CNV)[which(names(CNV) == 'SeqWGS_Cp_Hyperdiploid_Call')] <- 'HRD'
    names(CNV)[which(names(CNV) == 'SeqWGS_Cp_17p13_50percent')] <- 'Del17p'
    names(CNV)[which(names(CNV) == 'SeqWGS_Cp_1q21_50percent')] <- 'Gain1q'
    
    # GENE EXPRESSION - TPM
    
    TPM <- read.table(tpm_path, sep = '\t', header = T, stringsAsFactors = F)
    TPM <- TPM[c(1,grep('1_BM', names(TPM)))]
    TPM_GENES <- filter(TPM, GENE_ID %in% c('ENSG00000138185', 'ENSG00000188389', 'ENSG00000111291')) # CD39/ENTPD1, PDL1/PDCD1 and GPRC5D
    TPM_GENES <- melt(TPM_GENES, variable.name = 'sample', value.name = 'counts', id = 'GENE_ID')
    TPM_GENES <- dcast(TPM_GENES, sample ~ GENE_ID, value.var = 'counts')
    names(TPM_GENES)[which(names(TPM_GENES) == 'ENSG00000138185')] <- 'CD39'
    names(TPM_GENES)[which(names(TPM_GENES) == 'ENSG00000188389')] <- 'PDL1'
    names(TPM_GENES)[which(names(TPM_GENES) == 'ENSG00000111291')] <- 'GPRC5D'
    TPM_GENES$sample <- as.character(TPM_GENES$sample)
    TPM_GENES <- mutate(TPM_GENES, public_id = stri_sub(sample, 0, 9))
    
    # Combine datasets
    combined <- left_join(clinical, survival, by = 'public_id')
    combined <- left_join(combined, TPM_GENES, by = 'public_id')
    combined <- left_join(combined, select(SV, public_id, MAF_MAFA_MAFB, MMSET, MYC, CCND1), by = 'public_id')
    combined <- left_join(combined, select(CNV, public_id, HRD, Del17p, Gain1q), by = 'public_id')
    
    # Generate variables
    combined <- within(combined, 
                            CD39_3 <- as.integer(cut(CD39, 
                                                     c(0, 2, 10, max(CD39, na.rm = T)), 
                                                     include.lowest=TRUE)))

    combined <- within(combined, 
                            CD39_2 <- as.integer(cut(CD39, 
                                                     c(0, 2, max(CD39, na.rm = T)), 
                                                     include.lowest=TRUE)))
    combined <- within(combined, 
                            PDL1_2 <- as.integer(cut(PDL1, 
                                                     c(0, 1, max(PDL1, na.rm = T)), 
                                                     include.lowest=TRUE)))
    combined <- within(combined, 
                        GPRC5D_med <- as.integer(cut(GPRC5D, 
                                                 c(0, median(GPRC5D, na.rm = T), max(GPRC5D, na.rm = T)),
                                                 include.lowest=TRUE)))
    combined <- within(combined, 
                        GPRC5D_Q <- as.integer(cut(GPRC5D, 
                                                 quantile(GPRC5D, probs=0:4/4, na.rm = T), 
                                                 include.lowest=TRUE)))
    
    combined <- mutate(combined,
                       D_PT_therclass = factor(D_PT_therclass,
                                            levels = c('Bort-based', 'Car-based', 'IMID-based', 'Bort/Car-based', 
                                                       'Bort/IMID-based', 'Car/IMID-based', 'Bort/Car/IMID-based')),
                       ISS = factor(ISS),
                       D_PT_gender = factor(D_PT_gender),
                       MAF_MAFA_MAFB = factor(MAF_MAFA_MAFB),
                       MMSET = factor(MMSET),
                       MYC = factor(MYC),
                       CCND1 = factor(CCND1),
                       HRD = factor(HRD), 
                       Del17p = factor(Del17p), 
                       Gain1q = factor(Gain1q),
                       Triplet_induction = factor(Triplet_induction),
                       ASCT_upfront = factor(ASCT_upfront),
                       Maintenance = factor(Maintenance),
                       CD39_3 = factor(CD39_3),
                       CD39_2 = factor(CD39_2),
                       PDL1_2 = factor(PDL1_2),
                       GPRC5D_med = factor(GPRC5D_med),
                       GPRC5D_Q = factor(GPRC5D_Q)
                       ) 
    
    return(combined)
}



