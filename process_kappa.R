# Dependencies
suppressMessages(library(dplyr))

# Function call

process_kappa <- function(kappa_path, kappablast_path, clinical){
    ## Import Kappa
    kappa <- read.table(kappa_path, header = F, sep = ',', stringsAsFactors = F, blank.lines.skip = T, strip.white = T, na.strings = '')
    kappa <- kappa[!is.na(kappa$V1),]
    names(kappa) <- c("Rank", "Sequence", "Length", "Merge_count", "V_gene", "J_gene", 
    "total_reads_percent", "Cumulative_percent", "V_coverage")
    kappa <- kappa[-which(kappa$Rank == 'Rank'),]
    kappa$tube_barcode <- NA
    kappa$total_count <- NA
    for(i in seq(1,nrow(kappa), 12)){
        tube = kappa[i,2]
        count = kappa[(i+1), 2]
        kappa$tube_barcode[i:(i+11)] <- tube
        kappa$total_count[i:(i+11)] <- count   
    }
    kappa <- kappa[-which(kappa$Rank %in% c('Sample', 'Total count')),]
    rownames(kappa) <- NULL
    kappa$total_count <- gsub(',', '.', kappa$total_count)
    kappa$total_reads_percent <- as.numeric(kappa$total_reads_percent)
    kappa$Rank <- as.numeric(kappa$Rank)
    
    # Merge Kappa    
    kappa$tube_barcode <- as.numeric(kappa$tube_barcode)
    kappa <- filter(kappa, tube_barcode %in% clinical$tube_barcode)
    kappa <- left_join(kappa, 
                       select(clinical, capture_summary, IGK_call, lambda, tube_barcode, KV_simple, KJ_simple), 
                       by = 'tube_barcode')
        
    ## Import kappa blast
    kappa_blast <- read.csv(kappablast_path, header = T, stringsAsFactors = F)
    kappa_blast$tube_barcode <- as.numeric(kappa_blast$tube_barcode)

    ## Add kappa blast
    kappa <- left_join(kappa, select(kappa_blast, Productive, tube_barcode, Avg_identity_percent), by = "tube_barcode")
    kappa$Avg_identity_percent <- as.numeric(kappa$Avg_identity_percent)
    kappa <- mutate(kappa, IGK_SHM = ifelse(Avg_identity_percent >98, 1, 0))
    
    # Mutate kappa
    kappa <- mutate(kappa, lambda = ifelse(lambda == 1, 'Lambda', 'Kappa'),
                           clonality = ifelse(IGK_call == "C", 'Clonality detected', 'Negative'))
                    
    return(kappa)
}