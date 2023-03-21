# Dependencies
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))

# Purity function definition

estimate_purity <- function(somatic_variants, CNV_seg){
        samples <- unique(somatic_variants$LEUKID)
        somatic_variants$CN <- 2
    
        CNV_seg <- mutate(CNV_seg, CN = ifelse(seg.mean > 0, 3, 1))
        
        output <- data.frame(matrix(ncol = ncol(somatic_variants), nrow = 0))
        names(output) <- names(somatic_variants)
        
        # Determine copy number state
        
        for (i in 1:length(samples)){
            pt_var <- filter(somatic_variants, LEUKID == samples[i])
            pt_cnv <- filter(CNV_seg, LEUKID == samples[i])
            if (nrow(pt_cnv) == 0){
                    output <- rbind(output, pt_var)
                    next
            }
            cnv_gr <- with(pt_cnv, GRanges(chrom, IRanges(start = start, end = end)))
            var_gr <- with(pt_var, GRanges(CHR, IRanges(start = START, end = START)))
            ol <- suppressWarnings(findOverlaps(query = var_gr, subject = cnv_gr, type = "within"))
            if(length(ol) > 0){
                pt_var$CN[queryHits(ol)] <- pt_cnv$CN[subjectHits(ol)]
                output <- rbind(output, pt_var)
                next
            }
            output <- rbind(output, pt_var)
        }
        
        # Correcting VAF for CN by multiplying VAF and CN state (1, 2 or 3)
        # I need to take into account amplification of mut vs wt.  
        
        output <- mutate(output, VAF_CN_CORRECTED = TARGET_VAF * CN)
        
        # Taking the highest VAF variant from each patient as estimated purity
        purity_by_pt <- output %>% 
                group_by(LEUKID) %>%
                summarise(PURITY = max(VAF_CN_CORRECTED))
        
        output <- left_join(output, purity_by_pt[c('LEUKID', 'PURITY')], by = 'LEUKID')

        return(output)
}
