make_table_S2 <- function(clinical){
    
    # SPEP
    SPEP_data <- clinical[c("SPEP", "capture_num")];SPEP_data <- SPEP_data[complete.cases(SPEP_data),]
    SPEP_logistic <- glm(SPEP_data$capture_num ~ SPEP_data$SPEP, family = 'binomial')
    
    # FLC
    FLC_data <- clinical[c("FLC_involved_value", "capture_num")];FLC_data <- FLC_data[complete.cases(FLC_data),]
    FLC_logistic <- glm(FLC_data$capture_num ~ I(FLC_data$FLC_involved_value/100), family = 'binomial')
    
    # Gender
    gender_logistic <- glm(clinical$capture_num ~ clinical$gender, family = 'binomial')
    
    # Age
    age_logistic <- glm(clinical$capture_num ~ clinical$age, family = 'binomial')
    
    # Ig Heavy
    heavy_data <- clinical[c("capture_num","IFE_heavy")];heavy_data <- heavy_data[complete.cases(heavy_data),]
    heavy_data <- mutate(heavy_data, IgA = ifelse(IFE_heavy == 'IgA', 1, 0))
    heavy_data <- mutate(heavy_data, IgG = ifelse(IFE_heavy == 'IgG', 1, 0))
    heavy_data <- mutate(heavy_data, None = ifelse(IFE_heavy == 'None', 1, 0))
    IgG_logistic <- glm(heavy_data$capture_num ~ heavy_data$IgG, family = 'binomial')
    IgA_logistic <- glm(heavy_data$capture_num ~ heavy_data$IgA, family = 'binomial')
    None_logistic <- glm(heavy_data$capture_num ~ heavy_data$None, family = 'binomial')

    #Table
    table_S2 <- data.frame(
        Predictor = c('IgG',
                      'IgA',
                      'No heavy chain',
                      'Serum M-spike (1 g/dl)',
                      'Involved s-FLC (100 mg/dl)',
                      'Gender (male)',
                      'Age at sample (1 year)'
                      ),
        n = c(nrow(heavy_data),
              nrow(heavy_data),
              nrow(heavy_data),
              nrow(SPEP_data),
              nrow(FLC_data),
              nrow(clinical),
              nrow(clinical)
        ),
        OR = c(round(exp(summary(IgG_logistic)$coeff[2]), 2),
               round(exp(summary(IgA_logistic)$coeff[2]), 2),
               round(exp(summary(None_logistic)$coeff[2]), 2),
               round(exp(summary(SPEP_logistic)$coeff[2]), 2),
               round(exp(summary(FLC_logistic)$coeff[2]), 2),
               round(exp(summary(gender_logistic)$coeff[2]), 2),
               round(exp(summary(age_logistic)$coeff[2]), 2)
        ),
        CI = c(paste(round(exp(confint(IgG_logistic)[2,1]),2), 
                     round(exp(confint(IgG_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(IgA_logistic)[2,1]),2), 
                     round(exp(confint(IgA_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(None_logistic)[2,1]),2), 
                     round(exp(confint(None_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(SPEP_logistic)[2,1]),2), 
                     round(exp(confint(SPEP_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(FLC_logistic)[2,1]),2), 
                     round(exp(confint(FLC_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(gender_logistic)[2,1]),2),
                     round(exp(confint(gender_logistic)[2,2]),2), sep = "-"),
               paste(round(exp(confint(age_logistic)[2,1]),2), 
                     round(exp(confint(age_logistic)[2,2]),2), sep = "-")
        ),
        p = c(round(summary(IgG_logistic)$coeff[8],2),
              round(summary(IgA_logistic)$coeff[8],2),
              round(summary(None_logistic)$coeff[8],2),
              round(summary(SPEP_logistic)$coeff[8],2),
              round(summary(FLC_logistic)$coeff[8],2),
              round(summary(gender_logistic)$coeff[8],2),
              round(summary(age_logistic)$coeff[8],2)
        ))

    names(table_S2) <- c("Predictor", 
                                 "n", 
                                 "Odds ratio", 
                                 "95 % CI", 
                                 "p-value"
                        )
    return(table_S2)
}