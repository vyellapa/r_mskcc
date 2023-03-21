suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))

fishseq_purity <- function(fishseq, purity, plotname) {
    # Annotate fishseq with purity
    
    toplot <- left_join(fishseq, select(purity, LEUKID, PURITY), by = "LEUKID")
    
    toplot <- filter(toplot, PURITY > 0.2)

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
    
    ggsave(paste(plotname, '_fishseq_purity.png', sep = ''), plot1, height = 6, width = 10)
    
    return(plot1)
}


