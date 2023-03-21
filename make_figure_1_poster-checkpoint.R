# Dependencies
library(ggplot2)

make_figure_1_poster <- function(clinical){
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
        '0' = paste('Kappa (n = ', 
                    count(combined_data, 
                          lambda)$n[1], 
                    ')', 
                    sep = ''),
        '1' = paste('Lambda (n = ', 
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
    figure_1 <- ggplot(combined_data, 
                       aes(aspirate, 
                           fitted, 
                           col = myTYPE))+
                    geom_point(data = filter(combined_data, 
                                             myTYPE == 'No'), 
                               aes(aspirate, 
                                   capture_num), 
                               alpha = 0.8, 
                               size = 2,
                               position=position_jitter(height=0.03, width=0))+
                    geom_point(data = filter(combined_data, 
                                             myTYPE == 'Yes'), 
                               aes(aspirate, capture_num), 
                               alpha = 0.8,
                               size = 2,
                               position=position_jitter(height=0.03, 
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
                          text=element_text(size=20))
    return(figure_1)
}