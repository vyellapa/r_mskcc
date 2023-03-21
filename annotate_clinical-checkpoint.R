# Dependencies
suppressMessages(library(dplyr))

# Function call

annotate_clinical <-function(clinical,
                             invivoscribe,
                             leader,
                             calls,
                             flow, 
                             key, 
                             somatic,
                             cytogenetics) {
    # Merging clinical data with invivoscribe results
    ## Drop duplicated sample runs - selecting the ones with capture success if different
    invivoscribe$sample_ID <- as.character(invivoscribe$sample_ID)
    invivoscribe <- filter(invivoscribe, !tube_barcode %in% c('8019250217', '8019250150')) 

    ## Drop samples/patients not in invivoscribe datafile
    clinical <- clinical[clinical$sample_ID %in% invivoscribe$sample_ID,] 

    ## Merge with tube barcodes
    clinical <- left_join(clinical, 
                     invivoscribe[c("sample_ID", "tube_barcode")], 
                     by = 'sample_ID')

    ## Merge with FR1-3 and IGK calls
    clinical <- left_join(clinical, 
                     calls, 
                     by = 'tube_barcode')
    
    ## Merge with Leader calls and update 'All_call' column with leader data
    names(leader)[1] <- 'tube_barcode'
    clinical <- left_join(clinical, 
                     leader[c('tube_barcode', 'Leader')], 
                     by = 'tube_barcode')
    clinical <- mutate(clinical, 
                       All_call = ifelse(Leader %in% c('C'), 'C', All_call))
    
    # Merge with flow data
    clinical <- left_join(clinical, 
                     flow[c("sample_ID", "flow_PC", "CD56_pos", 
                                      "CD20_pos", "CD117_pos")], 
                     by = 'sample_ID')
    
    # Merge with cell preparation data
    key$sample_ID <- as.character(key$sample_ID)
    key <- key[!duplicated(key$MRN),]
    clinical <- left_join(clinical, 
                     key[c('sample_ID', 'cell_type')], 
                     by = 'sample_ID')
    
    # Merge with mutation data
    names(clinical)[which(names(clinical) == 'LEUKGEN_ID')] <- 'LEUKID'
    ## Define helper function
    p <- function(v) {
        Reduce(f=paste, x = v)
    }

    mytype_done <- data.frame(LEUKID = unique(cytogenetics$LEUKID),
                             ANY_myTYPE = 0,
                             stringsAsFactors = F)
    somatic <- somatic %>% 
                    group_by(LEUKID) %>%
                    summarise(mutations = p(GENE))
    cytogenetics <- cytogenetics %>%
                        filter(SEQ_result == 1) %>%
                        group_by(LEUKID) %>% 
                        summarise(cytogenetics = p(Cytogenetic_abnormality))
    clinical <- left_join(clinical, 
                     somatic, 
                     by = 'LEUKID')
    clinical <- left_join(clinical, 
                     cytogenetics, 
                     by = 'LEUKID')
    clinical <- left_join(clinical, 
                     mytype_done, 
                     by = 'LEUKID')
    clinical <- mutate(clinical, ANY_myTYPE = as.factor(ifelse(is.na(cytogenetics) & is.na(mutations), ANY_myTYPE, 1)))
    
    return(clinical)

}