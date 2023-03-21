# Dependencies
library(ggplot2)
suppressMessages(library(ggpubr))

# Function definition
make_figure_3 <- function(clinical){
    fig3a <- ggplot(data = clinical, aes(biopsy, aspirate))+
                geom_point(size=1.5)+
                geom_smooth(method='lm', size=1.5, col = 'black')+
                labs(y='Aspirate smear (%)',
                     x='Core biopsy (%)')+
                scale_y_continuous(breaks = seq(0,100,20))+
                scale_x_continuous(breaks = seq(0,100,20))+
                geom_line(data = data.frame(x=c(0:100), 
                                            y=c(0:100)),
                         aes(x, y),
                         size = 1,
                         linetype = 5)+
                theme_bw()+
                theme(text=element_text(size=15))
    
    fig3b <- ggplot(data = clinical, aes(aspirate, flow_PC))+
                geom_point(size=1.5)+
                geom_smooth(method='lm', size=1.5, col = 'black')+
                labs(x='Aspirate smear (%)',
                     y='Flow cytometry (%)')+
                scale_y_continuous(breaks = seq(0,100,20))+
                scale_x_continuous(breaks = seq(0,100,20))+
                geom_line(data = data.frame(x=c(0:100), 
                                        y=c(0:100)),
                     aes(x, y),
                     size = 1,
                     linetype = 5)+
                theme_bw()+
                theme(text=element_text(size=15))
    
    figure_3 = ggarrange(fig3a, 
                         fig3b, 
                         ncol = 2,
                         nrow = 1,
                         labels = c('A', 'B'))
    return(figure_3)
}