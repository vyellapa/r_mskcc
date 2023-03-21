# Function definition
make_table_S4 <- function(clinical){
    # Define table df
    table_S4 <- data.frame(
        kappa_rearrangement = c(
            'Intron-Kde',
            'V-Kde',
            'V-J',
            'V-undetermined'),
        lambda_mm = data.frame(table(filter(clinical, lambda == 1)$IGK_simple))[2],
        kappa_mm = data.frame(table(filter(clinical, lambda == 0)$IGK_simple))[2]
        )

        # Set column names
        names(table_S4) <- c(
                'Kappa rearrangement',
                paste('Lambda myeloma (n = ', count(filter(clinical, IGK == 1), lambda)$n[2], ')', sep = ''),
                paste('Kappa myeloma (n = ', count(filter(clinical, IGK == 1), lambda)$n[1], ')', sep = '')
    )
    
    return(table_S4)
}