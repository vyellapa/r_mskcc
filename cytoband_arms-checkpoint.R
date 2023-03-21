# Dependencies
suppressPackageStartupMessages(library(dplyr))
suppressMessages(library(stringi))
options(scipen=999)

# Define function
cytoband_arms <- function(path){
    cyto <- read.csv(path, stringsAsFactors = F, sep = '\t', header = F)

    # Cytoband format
    cyto <- cyto[-5]
    names(cyto) <- c('chrom', 'start', 'end', 'band')
    cyto <- filter(cyto, !chrom %in% c('chrX', 'chrY'))
    cyto$chrom <- stri_sub(cyto$chrom, 4)
    cyto$arm <- stri_sub(cyto$band, 0, 1)

    arms <- data.frame(cyto %>% 
            group_by(chrom, arm) %>%
            summarise(start = min(start), 
                      end = max(end)))

    arms <- mutate(arms, arm_size = end - start)
    return(arms)
}