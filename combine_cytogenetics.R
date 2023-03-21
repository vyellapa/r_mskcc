# Dependencies
suppressMessages(library(stringi))
suppressMessages(library(dplyr))

# Function definition
# Combining cytogenetics

combine_cytogenetics <- function(cnv, sv){
    names(cnv)[2] <- 'Cytogenetic_abnormality'
    names(sv)[2] <- 'Cytogenetic_abnormality'
    combined <- rbind(cnv, sv)
    fish_seq_vars <- c("Gain(1q)", "Del(1p)", "Del(12p)", 
                     "Del(13q)", "Del(16q)", "Del(17p)", unique(grep('t(', combined$Cytogenetic_abnormality, value = T, fixed = T)))
    combined <- filter(combined, Cytogenetic_abnormality %in% fish_seq_vars)
    
    combined$FISH_result[is.na(combined$FISH_result)] <- 3

    combined <- mutate(combined, concordance = 
                       factor(ifelse(SEQ_result == 1 & FISH_result == 1, 1,
                                     ifelse(SEQ_result == 1 & FISH_result == 0, 2,
                                            ifelse(SEQ_result == 0 & FISH_result == 1, 3,
                                                   ifelse(SEQ_result == 0 & FISH_result == 0, 4,
                                                          ifelse(SEQ_result == 0 & FISH_result == 3, 5, 6))))),
                              levels = c(1:6),
                              labels = c('Both positive', 'SEQ only', 'FISH only', 
                                         'Both negative', 'SEQ neg - FISH NA', 'SEQ pos - FISH NA')))

    cyto_levels = c("Del(1p)", "Gain(1q)", "Del(12p)", "Del(13q)", "Del(16q)",
                    "Del(17p)", "t(4;14)", "t(6;14)", "t(11;14)", "t(14;16)", "t(14;20)")
    
    combined$Cytogenetic_abnormality <- factor(combined$Cytogenetic_abnormality,
                                               levels = cyto_levels)

    return(combined)
    }


    
