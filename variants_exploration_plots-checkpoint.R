suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))



good_unknowns_plot <- function(combined, plotname) {
    # Good unknowns for plotting
    snps <- filter(combined, MANUAL_ANNOTATION == 'SNP')
    artifacts <- filter(combined, MANUAL_ANNOTATION == 'ARTIFACT')

    snps_95 <- quantile(snps$VAF_IF_CORRECT, probs = c(0.025, 0.975))
    artifacts_95 <- quantile(artifacts$VAF_IF_CORRECT, probs = c(0.025, 0.975))

    unknown <- filter(combined, MANUAL_ANNOTATION == 'UNKNOWN')
    unknown <- mutate(unknown, GOOD_UNKNOWN = ifelse(
            FILTER == 'PASS' & TARGET_DEPTH >= 100, 'GOOD_UNKNOWN', 'UNKNOWN'))
    
    # Color palette
    annpal <- c("#CC33FF", "#377EB8", "#984EA3", "#4DAF4A")
    names(annpal) <- c('GOOD_UNKNOWN', 'UNKNOWN', 'SNP', 'ARTIFACT')
    colScale <- scale_color_manual(name = 'Variant annotation category', values = annpal)

    # Plotting
    plot1 <- ggplot()+
                geom_density(data = snps, aes(VAF_IF_CORRECT, col = MANUAL_ANNOTATION), size = 1.5)+
                geom_density(data = artifacts, aes(VAF_IF_CORRECT, col = MANUAL_ANNOTATION), size = 1.5)+
                geom_jitter(data = unknown, aes(VAF_IF_CORRECT, y = 6, col = GOOD_UNKNOWN),
                            height = 4, size = 1.5)+
                geom_vline(xintercept=snps_95[1], color = "#984EA3", size = 1, linetype = 2)+
                geom_vline(xintercept=snps_95[2], color = "#984EA3", size = 1, linetype = 2)+
                geom_vline(xintercept=artifacts_95[1], color = "#4DAF4A", size = 1, linetype = 2)+
                geom_vline(xintercept=artifacts_95[2], color = "#4DAF4A", size = 1, linetype = 2)+
                colScale+
                #scale_color_manual(values=rev(c("#E41A1C", "#FF7F00", "#377EB8", "#984EA3", "#4DAF4A")))+
                scale_x_continuous(breaks=seq(0,1,by=0.1), 
                                limits = c(0,1))+
                labs(
                    x = 'Target VAF',
                    y = 'Density',
                    col = 'Variant annotation category'
                )
                theme_bw()

    ggsave(paste(plotname, '_good_unknowns.png', sep = ''), plot1)
    
    return(plot1)
}
    

all_variants_densityplot <- function(combined, plotname) {

    combined$ANNOTATION_CORRECTED <- factor(combined$ANNOTATION_CORRECTED, 
                    levels = c('ONCOGENIC', 'LIKELY', 'GOOD_UNKNOWN', 'UNKNOWN', 'SNP', 'ARTIFACT'), 
                    ordered = T)

    # Color palette
    annpal <- c("#E41A1C", "#FF7F00", "#CC33FF", "#377EB8", "#984EA3", "#4DAF4A")
    names(annpal) <- levels(combined$ANNOTATION_CORRECTED)
    colScale <- scale_color_manual(name = 'Variant annotation category', values = annpal)

    # Plotting
    plot1 <- ggplot()+
                geom_density(data = filter(combined, ANNOTATION_CORRECTED != "GOOD_UNKNOWN"), 
                            aes(VAF_IF_CORRECT, col = ANNOTATION_CORRECTED), 
                            size = 1.5)+
                geom_density(data = filter(combined, ANNOTATION_CORRECTED == "GOOD_UNKNOWN"), 
                            aes(VAF_IF_CORRECT, col = ANNOTATION_CORRECTED), 
                            size = 1.5)+
                scale_x_continuous(
                            breaks=seq(0,1,by=0.1), 
                            limits = c(0,1))+
                colScale+
                labs(
                    x = 'Target VAF',
                    y = 'Density'
                )+
                theme_bw()

    ggsave(paste(plotname, '_annot_density.png', sep = ''), plot1)
    
    return(plot1)
}