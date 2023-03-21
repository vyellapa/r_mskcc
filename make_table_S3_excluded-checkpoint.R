# Function definition
make_table_S3 <- function(clinical){
    # Subgroup definitions
    myTYPE_success <- filter(clinical, ANY_myTYPE == 1)
    
    # Define table df
    table_S3 <- data.frame(
        Assays = c(
                'All',
                'Leader excluded',
                'FR1 excluded',
                'FR2 excluded',
                'FR3 excluded',
                'IGK excluded',
                'IGK and FR3'
                ),
        Whole_cohort = c(
                paste(count(clinical, vars = capture_summary)$n[1], ' (',
                      round(count(clinical, vars = capture_summary)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_no_leader)$n[1], ' (',
                      round(count(clinical, vars = capture_no_leader)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_no_FR1)$n[1], ' (',
                      round(count(clinical, vars = capture_no_FR1)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_no_FR2)$n[1], ' (',
                      round(count(clinical, vars = capture_no_FR2)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_no_FR3)$n[1], ' (',
                      round(count(clinical, vars = capture_no_FR3)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_no_IGK)$n[1], ' (',
                      round(count(clinical, vars = capture_no_IGK)$n[1]/nrow(clinical)*100), ' %)', sep = ''),
                paste(count(clinical, vars = capture_IGK_FR3)$n[1], ' (',
                      round(count(clinical, vars = capture_IGK_FR3)$n[1]/nrow(clinical)*100), ' %)', sep = '')
                ),
        myTYPE_positive = c(
                paste(count(myTYPE_success, vars = capture_summary)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_summary)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_no_leader)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_no_leader)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_no_FR1)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_no_FR1)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_no_FR2)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_no_FR2)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_no_FR3)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_no_FR3)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_no_IGK)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_no_IGK)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = ''),
                paste(count(myTYPE_success, vars = capture_IGK_FR3)$n[1], ' (',
                      round(count(myTYPE_success, vars = capture_IGK_FR3)$n[1]/nrow(myTYPE_success)*100), ' %)', sep = '')
                ),
        Lambda = c(
                paste(count(filter(clinical, lambda ==1), vars = capture_summary)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_summary)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_no_leader)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_no_leader)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_no_FR1)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_no_FR1)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_no_FR2)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_no_FR2)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_no_FR3)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_no_FR3)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_no_IGK)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_no_IGK)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==1), vars = capture_IGK_FR3)$n[1], ' (',
                      round(count(filter(clinical, lambda ==1), vars = capture_IGK_FR3)$n[1]/nrow(filter(clinical, lambda ==1))*100), ' %)', sep = '')
                ),
        Kappa = c(
                paste(count(filter(clinical, lambda ==0), vars = capture_summary)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_summary)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_no_leader)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_no_leader)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_no_FR1)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_no_FR1)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_no_FR2)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_no_FR2)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_no_FR3)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_no_FR3)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_no_IGK)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_no_IGK)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = ''),
                paste(count(filter(clinical, lambda ==0), vars = capture_IGK_FR3)$n[1], ' (',
                      round(count(filter(clinical, lambda ==0), vars = capture_IGK_FR3)$n[1]/nrow(filter(clinical, lambda ==0))*100), ' %)', sep = '')
                ))
    # Set column names
    names(table_S3) <- c(
            'Assay usage',
            paste('All (n = ', nrow(clinical), ')', sep = ''),
            paste('myTYPE positive (n = ', nrow(myTYPE_success), ')', sep = ''),
            paste('Lambda (n = ', count(clinical, lambda)$n[2], ')', sep = ''),
            paste('Kappa (n = ', count(clinical, lambda)$n[1], ')', sep = '')
    )
    
    return(table_S3)
}