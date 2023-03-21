# Dependencies
suppressMessages(library(stringi))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
library(reshape2)

# changing cnv names
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

cnv_tolong <- function(seq, fish) {
    seq[-1] <- as.data.frame(sapply(seq[-1], FUN =  gsub, pattern = 'YES', replacement = 1), stringsAsFactors = F)
    seq[-1] <- as.data.frame(sapply(seq[-1], FUN =  gsub, pattern = 'NO', replacement = 0), stringsAsFactors = F)
    seq[-1] <- as.data.frame(sapply(seq[-1], FUN = as.numeric))
    fish[-1] <- suppressWarnings(as.data.frame(sapply(fish[-1], FUN = as.numeric)))
    
    names(seq) <- cnv_names(seq, 'SEQ')
    names(fish) <- cnv_names(fish, 'FISH')
    overlap <- c("Gain(1q)", "Del(1p)", "Del(12p)", 
                 "Del(13q)", "Del(16q)", "Del(17p)")
    
    fish$FISH_NONE <- NA

    long <- left_join(seq, fish, by = 'LEUKID')
    
    long <- melt(long, 
                         id.vars = c('LEUKID', grep('FISH', names(long), value = T)),
                         variable.name = 'SEQ_variable',
                         value.name = 'SEQ_result'
    )
    
    long <- melt(long, 
                         id.vars = c('LEUKID', grep('SEQ', names(long), value = T)),
                         variable.name = 'FISH_variable',
                         value.name = 'FISH_result'
    )
    
    long <- transform(long, SEQ_variable = stri_split(SEQ_variable, fixed = '_', simplify = T)[,2])
    long <- transform(long, FISH_variable = stri_split(FISH_variable, fixed = '_', simplify = T)[,2])
    
    long <- filter(long, SEQ_variable == FISH_variable | (FISH_variable == 'NONE' & !SEQ_variable %in% overlap))

    long <- long[-4]
    names(long)[2] <- 'CNV'
    return(long)
}