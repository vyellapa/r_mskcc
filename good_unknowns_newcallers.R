suppressMessages(library(dplyr))

good_unknowns <- function(combined) {
    # Formatting
    combined$VAF_IF_CORRECT <- as.numeric(combined$VAF_IF_CORRECT)
    combined$TARGET_DEPTH<- as.numeric(combined$TARGET_DEPTH)

    # Thresholds for good unknowns
    snps <- filter(combined, MANUAL_ANNOTATION == 'SNP')
    artifacts <- filter(combined, MANUAL_ANNOTATION == 'ARTIFACT')

    snps_95 <- quantile(snps$VAF_IF_CORRECT, probs = c(0.025, 0.975))
    artifacts_95 <- quantile(artifacts$VAF_IF_CORRECT, probs = c(0.025, 0.975))

    # Define good unknowns
    combined <- mutate(combined, 
                   ANNOTATION_CORRECTED = ifelse(
                           MANUAL_ANNOTATION != 'UNKNOWN',
                           MANUAL_ANNOTATION,
                           ifelse(
                                   VAF_IF_CORRECT > artifacts_95[2] & VAF_IF_CORRECT < snps_95[1] & 
                                           FLAGS_ALL == 'PASS' & TARGET_DEPTH >= 100, 
                                   'GOOD_UNKNOWN',
                                   MANUAL_ANNOTATION)))
    
    return(combined)
}