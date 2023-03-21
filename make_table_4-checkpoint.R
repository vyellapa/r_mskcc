make_table_4 <- function(clinical){
    # Subgroup definitions
    myTYPE_success <- filter(clinical, ANY_myTYPE == 1)
    
    # Define table df
    table_4 <- data.frame(
        Assays = c(
                'All',
                'FR1',
                'FR2',
                'FR3',
                'IGK'
                ),
        Whole_cohort = c(
                paste(count(clinical, vars = capture_summary)$n[1], ' (',
                      round(count(clinical, vars = capture_summary)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(nrow(clinical[clinical$FR1_call == 'C',]), ' (',
                      round(nrow(clinical[clinical$FR1_call == 'C',])/nrow(clinical)*100), ' %)', sep = ''),
                paste(nrow(clinical[clinical$FR2_call == 'C',]), ' (',
                      round(nrow(clinical[clinical$FR2_call == 'C',])/nrow(clinical)*100), ' %)', sep = ''),
                paste(nrow(clinical[clinical$FR3_call == 'C',]), ' (',
                      round(nrow(clinical[clinical$FR3_call == 'C',])/nrow(clinical)*100), ' %)', sep = ''),
                paste(nrow(clinical[clinical$IGK_call == 'C',]), ' (',
                      round(nrow(clinical[clinical$IGK_call == 'C',])/nrow(clinical)*100), ' %)', sep = '')
                ),
        myTYPE_positive = c(
                paste(count(myTYPE_success, vars = capture_summary)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_summary)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$FR1_call == 'C',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$FR1_call == 'C',])/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$FR2_call == 'C',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$FR2_call == 'C',])/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$FR3_call == 'C',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$FR3_call == 'C',])/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(nrow(myTYPE_success[myTYPE_success$IGK_call == 'C',]), ' (',
                      round(nrow(myTYPE_success[myTYPE_success$IGK_call == 'C',])/nrow(myTYPE_success)*100), ' %)', sep = '')
                ),
        Lambda = c(
                paste(count(filter(clinical, lambda ==1), vars = capture_summary)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_summary)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR1_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR1_call == 'C',])/count(clinical, lambda)$n[2]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR2_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR2_call == 'C',])/count(clinical, lambda)$n[2]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR3_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$FR3_call == 'C',])/count(clinical, lambda)$n[2]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$IGK_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==1)[filter(clinical, lambda ==1)$IGK_call == 'C',])/count(clinical, lambda)$n[2]*100), ' %)', sep = '')
                ),
        kappa = c(
                paste(count(filter(clinical, lambda ==0), vars = capture_summary)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_summary)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR1_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR1_call == 'C',])/count(clinical, lambda)$n[1]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR2_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR2_call == 'C',])/count(clinical, lambda)$n[1]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR3_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$FR3_call == 'C',])/count(clinical, lambda)$n[1]*100), ' %)', sep = ''),
                paste(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$IGK_call == 'C',]), ' (',
                      round(nrow(filter(clinical, lambda ==0)[filter(clinical, lambda ==0)$IGK_call == 'C',])/count(clinical, lambda)$n[1]*100), ' %)', sep = '')
                )
    )
    
    # Set column names
    names(table_4) <- c(
        'Assay usage',
        paste('All (n = ', nrow(clinical), ')', sep = ''),
        paste('myTYPE positive (n = ', nrow(myTYPE_success), ')', sep = ''),
        paste('Lambda (n = ', count(clinical, lambda)$n[2], ')', sep = ''),
        paste('Kappa (n = ', count(clinical, lambda)$n[1], ')', sep = '')
    )
    return(table_4)
}