# Notes:
## Main results table, capture results by patient cohort

# Function definition
make_table_2 <- function(clinical){
    # Subgroup definitions
    myTYPE_success <- filter(clinical, ANY_myTYPE == 1)
    SMM <- clinical[clinical$disease_stage == 'SMM',]
    NDMM <- clinical[clinical$disease_stage == 'NDMM',]
    RRMM <- clinical[clinical$disease_stage == 'RRMM',]
    
    # Defining table data frame    
    table_2 <- data.frame(
        group = c(
                'Clonality overall',                    
                'myTYPE positive',                 
                'Clonality in myTYPE pos'
        ),
        All = c(
                paste(nrow(clinical[clinical$capture_summary == 'Success',]), ' (',
                      round(nrow(clinical[clinical$capture_summary == 'Success',])/nrow(clinical)*100), 
                      ' %)', sep = ''),
                paste(nrow(myTYPE_success), '/', nrow(clinical[!is.na(clinical$ANY_myTYPE),]), ' (',
                      round(nrow(myTYPE_success)/nrow(clinical[!is.na(clinical$ANY_myTYPE),])*100), 
                      ' %)', sep =''),
                paste(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success',])
                            /nrow(myTYPE_success)*100), ' %)', sep = '')
        ),
        SMM = c(
                paste(nrow(SMM[SMM$capture_summary == 'Success',]), ' (',
                      round(nrow(SMM[SMM$capture_summary == 'Success',])/nrow(SMM)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$disease_stage == 'SMM',]), 
                      '/', nrow(SMM[!is.na(SMM$ANY_myTYPE),]),
                      ' (', round(nrow(myTYPE_success[myTYPE_success$disease_stage == 'SMM',])/
                                          nrow(SMM[!is.na(SMM$ANY_myTYPE),])*100), 
                      ' %)', sep =''),
                paste(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                  myTYPE_success$disease_stage == 'SMM',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                        myTYPE_success$disease_stage == 'SMM',])
                            /nrow(myTYPE_success[myTYPE_success$disease_stage == 'SMM',])*100), ' %)', 
                      sep = '')
        ),
        NDMM = c(
                paste(nrow(NDMM[NDMM$capture_summary == 'Success',]), ' (',
                      round(nrow(NDMM[NDMM$capture_summary == 'Success',])/nrow(NDMM)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$disease_stage == 'NDMM',]), 
                      '/', nrow(NDMM[!is.na(NDMM$ANY_myTYPE),]),
                      ' (', round(nrow(myTYPE_success[myTYPE_success$disease_stage == 'NDMM',])/
                                          nrow(NDMM[!is.na(NDMM$ANY_myTYPE),])*100), 
                      ' %)', sep =''),
                paste(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                  myTYPE_success$disease_stage == 'NDMM',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                        myTYPE_success$disease_stage == 'NDMM',])/
                                    nrow(myTYPE_success[myTYPE_success$disease_stage == 'NDMM',])*100), 
                      ' %)', sep = '')
        ),
        RRMM = c(
                paste(nrow(RRMM[RRMM$capture_summary == 'Success',]), ' (',
                      round(nrow(RRMM[RRMM$capture_summary == 'Success',])/nrow(RRMM)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$disease_stage == 'RRMM',]), 
                      '/', nrow(RRMM[!is.na(RRMM$ANY_myTYPE),]),
                      ' (', round(nrow(myTYPE_success[myTYPE_success$disease_stage == 'RRMM',])/
                                          nrow(RRMM[!is.na(RRMM$ANY_myTYPE),])*100), 
                      ' %)', sep =''),
                paste(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                  myTYPE_success$disease_stage == 'RRMM',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$capture_summary == 'Success' &
                                                        myTYPE_success$disease_stage == 'RRMM',])/
                                    nrow(myTYPE_success[myTYPE_success$disease_stage == 'RRMM',])*100), 
                      ' %)', sep = '')
        ))

    # Set final column names
    names(table_2) <- c(
        'variable',
        paste('All (n = ', nrow(clinical), ')', sep = ''),
        paste('Smoldering multiple myeloma (n = ', nrow(SMM), ')', sep = ''),
        paste('Newly diagnosed multiple myelona (n = ', nrow(NDMM), ')', sep = ''),
        paste('Relapsing/refractory multiple myeloma (n = ', nrow(RRMM), ')', sep = '')
    )
    return(table_2)
}
