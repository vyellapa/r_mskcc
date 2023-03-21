suppressMessages(library(dplyr))
suppressMessages(library(stringi))

combine <- function(snv, indel) {
    snv$VT <- "SUB"
    snv$TARGET_VAFCORRECT_VAF <- NA
    if (!"TARGET_VAFCORRECT_VAF" %in% names(indel)) {
        indel$TARGET_VAFCORRECT_VAF <- NA
    }
    indel$VT <- "IND"
    common_vars <- intersect(names(snv), names(indel))
    snv <- snv[common_vars]
    indel <- indel[common_vars]
    combined <- rbind(snv, indel)
    combined <- combined[!is.na(combined$MANUAL_ANNOTATION),]
    combined <- mutate(combined, VAF_IF_CORRECT = ifelse(
        !is.na(TARGET_VAFCORRECT_VAF),
        TARGET_VAFCORRECT_VAF,
        TARGET_VAF 
    ))
    names(combined)[which(names(combined) == 'TARGET_NAME')] <- "LEUKID"
    combined <- mutate(combined, LEUKID = stri_sub(LEUKID, 1, -6))
    return(combined)
}