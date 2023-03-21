# Dependencies
suppressMessages(library(stringi))
suppressMessages(library(dplyr))
library(reshape2)


# changing sv names
sv_names <- function(sv){
    
    #set names vector
    ns <- names(sv)
    
    # Set sample ID
    ns[grep('TARGET_NAME|ID|LEUK', ns)] <- 'LEUKID'
    
    # Define helper function 
    helperfun <- function(n, type){
        chrs <- c(as.numeric(stri_extract_first_regex(n, '\\d+')),
                  as.numeric(stri_extract_last_regex(n, '\\d+'))
                  )
        chrs = sort(chrs)
        return(paste(type,
                     '_t(',
                     chrs[1],
                     ';',
                     chrs[2],
                     ')',
                     sep = ''))}
    
    # Apply helper function to names vector
    ns[grep('FISH', ns)] <- sapply(ns[grep('FISH', ns)], 
                                      FUN = helperfun, 
                                      type = 'FISH')
    ns[grep('SEQ', ns)] <- sapply(ns[grep('SEQ', ns)], 
                                      FUN = helperfun, 
                                      type = 'SEQ')
    return(ns)
}

sv_tolong <- function(sv) {
    sv[-1] <- suppressWarnings(as.data.frame(sapply(sv[-1], FUN = as.numeric)))
    
    names(sv) <- sv_names(sv)
    
    long <- melt(sv,
                 id.vars = c('LEUKID', grep('FISH', names(sv), value = T)),
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
    
    long <- filter(long, SEQ_variable == FISH_variable)

    long <- long[-4]
    names(long)[2] <- 'SV'
    return(long)
}