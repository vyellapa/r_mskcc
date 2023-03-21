# Function definition
make_table_3 <- function(clinical) {
    # Aspirate
    aspirate_data <- clinical[c("aspirate", "capture_num")]
    aspirate_data <- aspirate_data[complete.cases(aspirate_data),]
    aspirate_logistic <- glm(aspirate_data$capture_num ~ I(aspirate_data$aspirate/10), 
                             family = 'binomial')
    # Flow
    flow_data <- clinical[c("flow_PC", "capture_num")]
    flow_data <- flow_data[complete.cases(flow_data),]
    flow_logistic <- glm(flow_data$capture_num ~ I(flow_data$flow_PC/10), 
                         family = 'binomial')
    # Biopsy
    biopsy_data <- clinical[c("biopsy", "capture_num")]
    biopsy_data <- biopsy_data[complete.cases(biopsy_data),]
    biopsy_logistic <- glm(biopsy_data$capture_num ~ I(biopsy_data$biopsy/10),
                           family = 'binomial')
    # myTYPE
    myTYPE_data <- clinical[c("ANY_myTYPE", "capture_num")]
    myTYPE_data <- myTYPE_data[complete.cases(myTYPE_data),]
    myTYPE_logistic <- glm(myTYPE_data$capture_num ~ myTYPE_data$ANY_myTYPE, 
                           family = 'binomial')
    # Ig light
    light_data <- clinical[c("capture_num","lambda")]
    light_data <- light_data[complete.cases(light_data),]
    light_logistic <- glm(light_data$capture_num ~ light_data$lambda, 
                          family = 'binomial')
    # Combined
    combined_data <- clinical[c('aspirate', 'ANY_myTYPE', 'capture_num', 'capture_summary', 'lambda')]
    combined_data <- combined_data[complete.cases(combined_data),]
    combined_logistic <- glm(combined_data$capture_num ~ I(combined_data$aspirate/10) + combined_data$ANY_myTYPE + combined_data$lambda,
                             family = 'binomial')
    combined_data <- cbind(combined_data, fitted = combined_logistic$fitted)
    combined_data <- mutate(combined_data, myTYPE = factor(ifelse(ANY_myTYPE == 1, 'Yes', 'No'), 
                                                           levels = c('Yes', 'No')))
    # Table definition
    table_3 <- data.frame(
            Predictor = c('Univariate models',
                          'BMPC aspirate (10 %)',
                          'BMPC biopsy (10 %)',
                          'BMPC flow (10 %)',
                          'Lambda light chain',
                          'myTYPE positive',
                          'Multivariate model',
                          'BMPC aspirate (10 %)',
                          'myTYPE positive',
                          'Lambda light chain'
            ),
            n = c(NA,
                  nrow(aspirate_data),
                  nrow(biopsy_data),
                  nrow(flow_data),
                  nrow(light_data),
                  nrow(myTYPE_data),
                  nrow(combined_data),
                  NA,
                  NA,
                  NA
            ),
            OR = c(NA,
                   round(exp(summary(aspirate_logistic)$coeff[2]), 2),
                   round(exp(summary(biopsy_logistic)$coeff[2]), 2),
                   round(exp(summary(flow_logistic)$coeff[2]), 2),
                   round(exp(summary(light_logistic)$coeff[2]), 2),
                   round(exp(summary(myTYPE_logistic)$coeff[2]), 2),
                   NA,
                   round(exp(summary(combined_logistic)$coeff[2]), 2),
                   round(exp(summary(combined_logistic)$coeff[3]), 2),
                   round(exp(summary(combined_logistic)$coeff[4]), 2)

            ),
            CI = c(NA,
                   paste(round(exp(confint(aspirate_logistic)[2,1]),2),
                         round(exp(confint(aspirate_logistic)[2,2]),2), sep = "-"),
                   paste(round(exp(confint(biopsy_logistic)[2,1]),2),
                         round(exp(confint(biopsy_logistic)[2,2]),2), sep = "-"),
                   paste(round(exp(confint(flow_logistic)[2,1]),2), 
                         round(exp(confint(flow_logistic)[2,2]),2), sep = "-"),
                   paste(round(exp(confint(light_logistic)[2,1]),2), 
                         round(exp(confint(light_logistic)[2,2]),2), sep = "-"),
                   paste(round(exp(confint(myTYPE_logistic)[2,1]),2),
                         round(exp(confint(myTYPE_logistic)[2,2]),2), sep = "-"),
                   NA,
                   paste(round(exp(confint(combined_logistic)[2,1]),2),
                         round(exp(confint(combined_logistic)[2,2]),2), sep = "-"),
                   paste(round(exp(confint(combined_logistic)[3,1]),2),
                         round(exp(confint(combined_logistic)[3,2]),2), sep = "-"),
                   paste(round(exp(confint(combined_logistic)[4,1]),2),
                         round(exp(confint(combined_logistic)[4,2]),2), sep = "-")
            ),
            p = c(NA,
                  round(summary(aspirate_logistic)$coeff[8],3),
                  round(summary(biopsy_logistic)$coeff[8],3),
                  round(summary(flow_logistic)$coeff[8],3),
                  round(summary(light_logistic)$coeff[8],3),
                  round(summary(myTYPE_logistic)$coeff[8],3),
                  NA,
                  round(summary(combined_logistic)$coeff[14],3),
                  round(summary(combined_logistic)$coeff[15],3),
                  round(summary(combined_logistic)$coeff[16],3)
            ))
    # Final column names
    names(table_3) <- c("Predictor", "n", "Odds ratio", "95 % CI", "p-value")
    return(table_3)
}