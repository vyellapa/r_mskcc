# Dependencies
library(ggplot2)
suppressMessages(library(ggpubr))

make_figure_S1 <- function(FR1, FR2, FR3, Leader, IGK){
    corr_eqn <- function(x,y, digits = 3) {
        xy <- data.frame(x, y)
        corr_coef <- round(cor(xy$x, xy$y)^2, digits = digits)
        paste("R^2 == ", corr_coef)
    }


    # Make plot
    
    labels_S1a <- data.frame(x = 7, 
                     y = 9, 
                     label = corr_eqn(FR1$clonal_input_percent, FR1$clonal_reads_percent))
    
    fig_S1a <- ggplot(FR1, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F)+
        geom_point(size = 2)+
        labs(title = 'IGH FR1',
             x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1a, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE)+
        scale_x_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        scale_y_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    
    labels_S1b <- data.frame(x = 7, 
                     y = 9, 
                     label = corr_eqn(FR2$clonal_input_percent, FR2$clonal_reads_percent))
    
    fig_S1b <- ggplot(FR2, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F)+
        geom_point(size = 2)+
        labs(title = 'IGH FR2',
             x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1b, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE)+
        scale_x_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        scale_y_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    
    labels_S1c <- data.frame(x = 8, 
                     y = 14, 
                     label = corr_eqn(FR3$clonal_input_percent, FR3$clonal_reads_percent))
    
    fig_S1c <- ggplot(FR3, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F)+
        geom_point(size = 2)+
        labs(title = 'IGH FR3',
             x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1c, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE)+
        scale_x_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        scale_y_continuous(breaks=seq(0,15,3), 
                           limits = c(0,15))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    
    labels_S1d <- data.frame(x = 8, 
                     y = 9, 
                     label = corr_eqn(Leader$clonal_input_percent, Leader$clonal_reads_percent))
    
    fig_S1d <- ggplot(Leader, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F)+
        geom_point(size = 2)+
        labs(title = 'IGHV Leader',
             x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1d, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE)+
        scale_x_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        scale_y_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    
    labels_S1e <- data.frame(x = 15, 
                     y = 18, 
                     label = corr_eqn(IGK$clonal_input_percent, IGK$clonal_reads_percent))
    
    fig_S1e <- ggplot(IGK, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F)+
        geom_point(size = 2)+
        labs(title = 'IGK',
             x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1e, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE)+
        scale_x_continuous(breaks=seq(0,20,4), 
                           limits = c(0,20))+
        scale_y_continuous(breaks=seq(0,20,4), 
                           limits = c(0,20))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
        
    figure_S1 = ggarrange(fig_S1a,
                         fig_S1b,
                         fig_S1c,
                         fig_S1d,
                         fig_S1e,
                         ncol = 2,
                         nrow = 3,
                         labels = c('A', 'B', 'C', 'D', 'E'))
                    
    return(figure_S1)
}