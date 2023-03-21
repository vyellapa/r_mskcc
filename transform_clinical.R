# Dependencies
set.seed(353)
suppressMessages(library(dplyr))
suppressWarnings(library(stringi))

# Function definition
transform_clinical <- function(clinical){
    # Remove "treated" subgroup
    clinical <- filter(clinical, disease_stage != 'Treated')
    
    # Remove duplicated samples
    clinical <- filter(clinical, !tube_barcode %in% c("8019250148", "8019250208"))
    
    # 'imputing' BM aspirates
    ## Bone marrow aspirate differential was missing for 6 samples. 
    ## In 5 samples, the pathologist commented on likely or markedly hemodiluted sample. 
    ## For these samples, I set aspirate to random normal with mean 2.5 and SD 1. 
    ## For the last sample (159888), there was no comment, but aspirate differential was missing. 
    ## Flow cytometry of that sample showed 9 % plasma cells, which I use in the analysis. 

    clinical[clinical$sample_ID=='159888',]$aspirate <- 9
    clinical <- mutate(clinical, aspirate = (ifelse(is.na(aspirate), 
                                                            round(rnorm(1, mean=2.5),0), 
                                                            aspirate)))

    # Changing patients with no light chain class to lambda or kappa based on LC restriction.
    clinical[clinical$sample_ID=='124132',]$FLC_involved_class <- 'Lambda'
    clinical[clinical$sample_ID=='69056',]$FLC_involved_class <- 'Kappa'
    clinical[clinical$sample_ID=='32867',]$FLC_involved_class <- 'Lambda'
    # Sample '23842' is non-secretory, no light chain restriction in BMPC.
    
    # Generate new variables
    clinical <- mutate(clinical, lambda = as.factor(ifelse(FLC_involved_class == 'Lambda', 1, 0)))
    clinical <- mutate(clinical, flow_PC = round(as.numeric(flow_PC), 1))
    clinical <- mutate(clinical, capture_summary = factor(ifelse(All_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_num = as.numeric(ifelse(capture_summary == 'Success', 1, 0)))
    clinical <- mutate(clinical, myTYPE = factor(ifelse(is.na(ANY_myTYPE), 'ND', 
                                             ifelse(ANY_myTYPE == 1, 'pos', 'neg')),
                                             levels = c('pos', 'neg', 'ND')))
    clinical$DOB <- as.POSIXct(strptime(clinical$DOB, format = '%Y-%m-%d'))
    clinical$date_sample <- as.POSIXct(strptime(clinical$date_sample, format = '%Y-%m-%d'))
    clinical$age  <-  year(as.period(interval(clinical$DOB, clinical$date_sample))) #age at time of sample
    
    clinical <- mutate(clinical, 
                   SMM_FACTOR = ifelse(disease_stage == 'SMM', 1, 0),
                   NDMM_FACTOR = ifelse(disease_stage == 'NDMM', 1, 0),
                   RRMM_FACTOR = ifelse(disease_stage == 'RRMM', 1, 0)
                  )
    clinical <- mutate(clinical, 
                   KV_simple = stri_split_regex(KV_gene, '[1-9]', simplify = T)[,1],
                   KJ_simple = stri_split_regex(KJ_gene, '[1-9]', simplify = T)[,1],
                  )
    clinical <- mutate(clinical, 
                   KV_simple = gsub('IGK', '', KV_simple),
                   KJ_simple = gsub('IGK', '', KJ_simple),
                  )
    clinical <- mutate(clinical, 
                   IGK_REARR = ifelse(!is.na(KV_gene), paste(KV_gene, KJ_gene, sep = ':'), NA),
                   IGK_simple = ifelse(!is.na(KV_gene), paste(KV_simple, KJ_simple, sep = ':'), NA)
                  )
    # Capture rate with different assay usage - neither of these groups include leader.
    clinical <- mutate(clinical, capture_no_FR1 = factor(ifelse(FR2_call == 'C' | FR3_call == 'C' | IGK_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_no_FR2 = factor(ifelse(FR1_call == 'C' | FR3_call == 'C' | IGK_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_no_FR3 = factor(ifelse(FR1_call == 'C' | FR2_call == 'C' | IGK_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_no_IGK = factor(ifelse(FR1_call == 'C' | FR2_call == 'C' | FR3_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_no_leader = factor(ifelse(FR1_call == 'C' | FR2_call == 'C' | FR3_call == 'C' | IGK_call == 'C', 
                                                                   'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, capture_IGK_FR3 = factor(ifelse(FR3_call == 'C' | IGK_call == 'C', 'Success', 'Failure'),
                                                  levels = c('Success', 'Failure'),
                                                  ordered = T))
    clinical <- mutate(clinical, no_assays_success = factor(0, levels = c(0:4), ordered = T))
    
    clinical <- mutate(clinical,
                  FR1 = ifelse(FR1_call == 'C', 1, 0),
                  FR2 = ifelse(FR2_call == 'C', 1, 0),
                  FR3 = ifelse(FR3_call == 'C', 1, 0),
                  IGK = ifelse(IGK_call == 'C', 1, 0))
    
    for(i in 1:nrow(clinical)){
            inc <- 0
            if(clinical$FR1_call[i] == 'C'){
                    inc = inc + 1
            }
            if(clinical$FR2_call[i] == 'C'){
                    inc = inc + 1
            }
            if(clinical$FR3_call[i] == 'C'){
                    inc = inc + 1
            }
            if(clinical$IGK_call[i] == 'C'){
                    inc = inc + 1
            }
            clinical$no_assays_success[i] = inc
    }    
    return(clinical)
}