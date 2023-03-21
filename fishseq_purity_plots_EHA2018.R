suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))

fishseq_purity_EHA <- function(fishseq, purity, plotname) {
    # Annotate fishseq with purity
    toplot <- left_join(fishseq, select(purity, LEUKID, PURITY), by = "LEUKID")
    
    # Remove sorted CD138 cells, low purity and aberrations with no positives
    CD138 <- c("I-H-130564-T1-1", "I-H-130582-T1-1", "I-H-130587-T1-1", 
               "I-H-130603-T1-1", "I-H-130642-T1-1", "I-H-130672-T1-1")
    
    toplot <- filter(toplot, PURITY > 0.2, !LEUKID %in% CD138, Cytogenetic_abnormality!='t(14;20)')

    # Color palette
    cytopal <- c("#0571B0", "#92C5DE", "#CD5C5C", "#A9A9A9", "#DCDCDC", '#228B22')
    names(cytopal) <- c('Both positive', 'SEQ only', 'FISH only',
                        'Both negative', 'SEQ neg - FISH NA', 'SEQ pos - FISH NA')

    # Plotting
    plot1 <- ggplot(toplot, aes(Cytogenetic_abnormality, reorder(LEUKID, 10/as.numeric(concordance))))+
        geom_tile(aes(fill = concordance), color = 'white', size = 0.3)+
        scale_fill_manual(values = cytopal)+
        theme(legend.title = element_text(size = 14, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title=element_text(size=14,face="bold"),
              axis.text.x = element_text(size = 12),
              legend.position = 'bottom'
              ) +
        scale_x_discrete(position = "top")+
        labs(fill = "FISH and SEQ status",
             y = 'LEUKGEN_ID',
             x = 'Cytogenetic abnormality')
    
    ggsave(paste(plotname, '_fishseq_purity.pdf', sep = ''), plot1, height = 4, width = 10)
    
    return(plot1)
}


