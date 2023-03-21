suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))

fishseq_all <- function(fishseq, plotname) {
    
    # Color palette
    cytopal <- c("#0571B0", "#92C5DE", "#CD5C5C", "#A9A9A9", "#DCDCDC", '#228B22')
    names(cytopal) <- c('Both positive', 'SEQ only', 'FISH only',
                        'Both negative', 'SEQ neg - FISH NA', 'SEQ pos - FISH NA')

    # Plotting
    plot1 <- ggplot(fishseq, aes(Cytogenetic_abnormality, reorder(LEUKID, 10/as.numeric(concordance))))+
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

    ggsave(paste(plotname, '_fishseq_all.png', sep = ''), plot1, height = 25, width = 10)
    
    return(plot1)
}


