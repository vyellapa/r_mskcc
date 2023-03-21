# Dependencies
suppressMessages(library(dplyr))
suppressWarnings(suppressMessages(library(stringi)))
library(reshape2)

# Function definitions
cnv_names <- function(cnv, origin){
    # Origin = FISH or SEQ
    #set names vector
    ns <- names(cnv)
    
    # Set sample ID
    ns[grep('TARGET_NAME|ID|LEUK', ns)] <- 'LEUKID'
    
    # Define helper function 
    helperfun <- function(n, type){
        return(paste(origin,
                     '_', type, '(',
            str_extract(n, "\\d+\\w"),
                     ')',
                     sep = ''))}
    
    # Apply helper function to names vector
    ns[grep('Del|del', ns)] <- sapply(ns[grep('Del|del', ns)], 
                                      FUN = helperfun, 
                                      type = 'Del')
    ns[grep('gain|Gain|amp|Amp', ns)] <- sapply(ns[grep('gain|Gain|amp|Amp', ns)], 
                                      FUN = helperfun, 
                                      type = 'Gain')
    return(ns)
}


cnv_addfrac <- function(fishseq, fracs){
    fracs[-1] <- suppressWarnings(as.data.frame(sapply(fracs[-1], FUN = as.numeric)))
    fracs[-1] <- suppressWarnings(as.data.frame(sapply(fracs[-1], FUN = round, digits = 3)))
    names(fracs) <- cnv_names(cnv_frac, 'FRAC')
    fraclong <- melt(fracs, 
                 id.vars = 'LEUKID',
                 variable.name = 'Cytogenetic_abnormality',
                 value.name = 'FISH_frac'
    )
    fraclong <- transform(fraclong, Cytogenetic_abnormality = stri_split(Cytogenetic_abnormality, fixed = '_', simplify = T)[,2])
    fraclong <- filter(fraclong, !is.na(FISH_frac) & FISH_frac > 0)
    fishseq <- suppressWarnings(left_join(fishseq, fraclong, by = c('LEUKID', 'Cytogenetic_abnormality')))
    
    
    # Re-apply lost factor levels
    cyto_levels = c("Del(1p)", "Gain(1q)", "Del(12p)", "Del(13q)", "Del(16q)",
                    "Del(17p)", "t(4;14)", "t(6;14)", "t(11;14)", "t(14;16)", "t(14;20)")
    
    fishseq$Cytogenetic_abnormality <- factor(fishseq$Cytogenetic_abnormality,
                                               levels = cyto_levels)
    return(fishseq)
    
}