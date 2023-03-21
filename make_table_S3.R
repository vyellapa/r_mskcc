make_table_S3 <- function(clinical){
    table_S3 <- data.frame(
            Marker = c(
                    'CD20',
                    'CD56',
                    'CD117'
                    ),
            Positive = c(
                    paste(round(count(clinical[clinical$capture_summary == 'Success',], CD20_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD20_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Success',], CD20_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD20_pos)), ")", sep = ""),
                    paste(round(count(clinical[clinical$capture_summary == 'Success',], CD56_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD56_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Success',], CD56_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD56_pos)), ")", sep = ""),
                    paste(round(count(clinical[clinical$capture_summary == 'Success',], CD117_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD117_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Success',], CD117_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Success',]$CD117_pos)), ")", sep = "")
            ),
            Negative = c(
                    paste(round(count(clinical[clinical$capture_summary == 'Failure',], CD20_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD20_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Failure',], CD20_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD20_pos)), ")", sep = ""),
                    paste(round(count(clinical[clinical$capture_summary == 'Failure',], CD56_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD56_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Failure',], CD56_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD56_pos)), ")", sep = ""),
                    paste(round(count(clinical[clinical$capture_summary == 'Failure',], CD117_pos)$n[2]/sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD117_pos))*100), '% (',
                          count(clinical[clinical$capture_summary == 'Failure',], CD117_pos)$n[2], "/", sum(!is.na(clinical[clinical$capture_summary == 'Failure',]$CD117_pos)), ")", sep = "")
            ),
            p = c(
                    round(chisq.test(clinical[!is.na(clinical$CD20_pos),]$CD20_pos, clinical[!is.na(clinical$CD20_pos),]$capture_summary)$p.value,2),
                    round(chisq.test(clinical[!is.na(clinical$CD56_pos),]$CD56_pos, clinical[!is.na(clinical$CD56_pos),]$capture_summary)$p.value,2),
                    round(chisq.test(clinical[!is.na(clinical$CD117_pos),]$CD117_pos, clinical[!is.na(clinical$CD117_pos),]$capture_summary)$p.value,2)
            ))

    names(table_S3) <- c(
            'Cell surface marker',
            'Clonality detected',
            'Clonality not detected',
            'Chi^2 p-value'
    )
    return(table_S3)
}