# Dependencies
library(ggplot2)

make_figure_S1_FR1 <- function(FR1){
    corr_eqn <- function(x,y, digits = 3) {
        xy <- data.frame(x, y)
        corr_coef <- round(cor(xy$x, xy$y)^2, digits = digits)
        paste("R^2 == ", corr_coef)
    }

    # Make plot
    
    labels_S1a <- data.frame(x = 6, 
                     y = 9, 
                     label = corr_eqn(FR1$clonal_input_percent, FR1$clonal_reads_percent))
    
    fig_S1a <- ggplot(FR1, aes(clonal_input_percent, clonal_reads_percent))+
        geom_smooth(aes(x = clonal_input_percent, y = clonal_reads_percent), method = 'lm', se = F, size = 1.5)+
        geom_point(size = 2)+
        labs(x = 'Clonal DNA input (%)',
            y = 'Clonal reads (%)')+
        geom_text(data = labels_S1a, 
                  aes(x = x, 
                      y = y,
                      label = label),
                  parse = TRUE,
                  size = 8)+
        scale_x_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        scale_y_continuous(breaks=seq(0,10,2), 
                           limits = c(0,10))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              text=element_text(size=22))
                    
    return(fig_S1a)
}