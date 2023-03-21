# Dependencies
library(ggplot2)
suppressMessages(library(ggpubr))

# Function definition
make_figure_2 <- function(kappa){
    fig2a <- ggplot(filter(kappa, 
                       Rank < 7, 
                       IGK_call != 'I'), 
                aes(x = as.character(Rank), 
                    y = total_reads_percent, 
                    group = tube_barcode, 
                    col = clonality))+
            geom_point()+
            geom_line()+
            facet_wrap(~lambda)+
            labs(x = 'Clonal fraction rank',
                 y = 'Clonal fraction (%)',
                col = 'IGK assay result')+
            theme_bw()+
            theme(strip.background = element_rect(fill="white"),
                  text = element_text(size=15),
                  legend.position="top")
    
    fig2b <- ggplot(filter(kappa, 
                       Rank == 1, 
                       IGK_call == "C", 
                       KV_simple == 'V'), 
                aes(lambda, I(100-Avg_identity_percent)))+
                geom_point()+
                geom_hline(aes(yintercept = 2), 
                           col= 'red', 
                           size = 1.5)+
                geom_boxplot(fill=NA, 
                             size = 1,
                            outliers.shape = NA)+
                labs(x = 'Light chain restriction',
                     y = 'Somatic hypermutation of kappa variable region (%)')+
            theme_bw()+
            theme(text = element_text(size=15))

    figure_2 = ggarrange(fig2a, 
                     fig2b,
                     widths = c(1.8,1),
                     ncol = 2,
                     nrow = 1,
                     labels = c('A', 'B'))
    return(figure_2)
}