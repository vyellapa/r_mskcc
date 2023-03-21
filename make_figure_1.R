# Dependencies
library(ggplot2)
suppressMessages(library(ggpubr))

make_figure_1 <- function(clinical){
    # Build logistic regression model - identical code as for Table 3.
    combined_data <- clinical[c('aspirate', 'ANY_myTYPE', 'capture_num', 'capture_summary', 'lambda')]
    combined_data <- combined_data[complete.cases(combined_data),]
    combined_logistic <- glm(combined_data$capture_num ~ I(combined_data$aspirate/10) + combined_data$ANY_myTYPE + combined_data$lambda,
                             family = 'binomial')
    combined_data <- cbind(combined_data, fitted = combined_logistic$fitted)
    combined_data <- mutate(combined_data, myTYPE = factor(ifelse(ANY_myTYPE == 1, 'Yes', 'No'), 
                                                           levels = c('Yes', 'No')))
    # List of labels for figure facet_wrap
    lambda_names <- list(
        '0' = paste('Kappa light chain (n = ', 
                    count(combined_data, 
                          lambda)$n[1], 
                    ')', 
                    sep = ''),
        '1' = paste('Lambda light chain (n = ', 
                    count(combined_data, 
                          lambda)$n[2], 
                    ')', 
                    sep = '')
    )
    
    # Labeller to set facet_wrap labels
    regression_labeller <- function(variable, value){
            return(lambda_names[value])
    }
    
    # Make plot
    
    fig_1a <- ggplot(data = clinical, 
                     aes(reorder(capture_summary,
                                 desc(capture_summary)),
                         aspirate))+
                     geom_dotplot(aes(fill = capture_summary), 
                         col = NA, 
                         binaxis = 'y',
                         stackdir = 'center', 
                         dotsize = 0.75, 
                         stackratio = 1.5)+
                    geom_boxplot(fill=NA, 
                             size = 0.75,
                             outlier.size = -1)+
                    labs(y='BMPC (%)',
                         x='Clonality detection')+
                    scale_y_continuous(breaks = seq(0,100,20))+
                    scale_fill_manual(values = c("#E41A1C", "#0066ff"))+
                    theme_bw()+
                    theme(text=element_text(size=15))+
                    guides(fill=FALSE)
    
    fig_1b <- ggplot(data = combined_data, 
                     aes(reorder(myTYPE, desc(myTYPE)), 
                         aspirate))+
                     geom_dotplot(aes(fill = myTYPE), 
                         col = NA, 
                         binaxis = 'y',
                         stackdir = 'center', 
                         dotsize = 0.75, 
                         stackratio = 1.5)+
                    geom_boxplot(fill=NA, 
                             size = 0.75,
                             outlier.size = -1)+
                    labs(y='BMPC (%)',
                         x='Somatic aberration \ndetected by myTYPE')+
                    scale_y_continuous(breaks = seq(0,100,20))+
                    scale_fill_manual(values = c("#E41A1C", "#0066ff"))+
                    theme_bw()+
                    guides(fill=FALSE)+
                    theme(text=element_text(size=15))
    
    fig_1c <- ggplot(combined_data, 
                       aes(aspirate, 
                           fitted, 
                           col = myTYPE))+
                    geom_point(data = filter(combined_data, 
                                             myTYPE == 'No'), 
                               aes(aspirate, 
                                   capture_num), 
                               alpha = 0.7,
                               position=position_jitter(height=0.02, width=0))+
                    geom_point(data = filter(combined_data, 
                                             myTYPE == 'Yes'), 
                               aes(aspirate, capture_num), 
                               alpha = 0.7, 
                               position=position_jitter(height=0.02, 
                                                        width=0))+
                    geom_line(size = 1.5)+
                    scale_colour_manual(values = c("#0066ff", "#E41A1C"))+
                    scale_y_continuous(breaks = seq(0,1,0.2), 
                                       sec.axis = sec_axis(trans = ~., breaks = c(0,1), 
                                                           labels = c('Failure', 'Success'), 
                                                           name = 'Clonality detection results \n(dots)'))+
                    labs(x = 'BMPC (%)',
                         y = 'Probability of clonality detection \n(regression lines)',
                         col = 'Somatic aberration \ndetected by myTYPE')+
                    facet_wrap(~lambda, 
                               labeller = regression_labeller)+
                    theme_bw()+
                    theme(strip.background =element_rect(fill="white"),
                          legend.position="bottom",
                          text=element_text(size=15))
    
    fig_1ab = ggarrange(fig_1a, 
                 fig_1b, 
                 ncol = 1,
                 nrow = 2,
                 labels = c('A', 'B'))
    
    figure_1 = ggarrange(fig_1ab, 
                     fig_1c, 
                     ncol = 2,
                     nrow = 1,
                     labels = c(NA, 'C'),
                     widths = c(1, 2))
                    
    return(figure_1)
}