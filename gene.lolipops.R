#### Function to create gene mutation lolipop plots
## Takes Caveman output in MAF format.
## v.1.17.2018

#### Dependencies
library(dplyr)
library(maftools)

setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.seq')

## Import dependancy function from script ' toMAF' v. 1.10.2018
createMAF <- function(caveman) {
        MAFdf <- select(caveman, GENE, CHR, START, END, REF, ALT, EFFECT, VT, sampleLeukid, PROTEIN_CHANGE, TARGET_VAF)
        names(MAFdf) <- c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'AAChange', 't_vaf')
        MAFdf$Variant_Type <- 'SNP'
        MAFdf$Variant_Classification <-  gsub('non_synonymous_codon', 'Missense_Mutation', MAFdf$Variant_Classification)
        MAFdf$Variant_Classification <-  gsub('splice_site_variant', 'Splice_Site', MAFdf$Variant_Classification)
        MAFdf$Variant_Classification <-  gsub('initiator_codon_change', 'Start_Codon_SNP', MAFdf$Variant_Classification)
        MAFdf$Variant_Classification <-  gsub('stop_gained', 'Nonsense_Mutation', MAFdf$Variant_Classification)
        MAF <- read.maf(MAFdf, verbose = FALSE)
        return(MAF)
}

## Function for basic lollipop plot
loll <- function(gene.string, maf){
        lollipopPlot(maf = maf, 
                     gene = gene.string, 
                     labelPos = 'all', 
                     cBioPortal = TRUE, 
                     showDomainLabel = FALSE)
}

## Function to make and export lollipop plot for a list of genes
makeLoll <- function(gene.list, data){
        input.maf <- createMAF(data)
        for(i in 1:length(gene.list)){
                tryCatch({
                        g <- loll(gene.list[i], input.maf)
                        ggsave(paste('../datafiles/RRMM.annot.work/figures/', 
                                     as.character(gene.list[i]), 
                                     '.png', 
                                     sep = ''))},
                        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
}

## Run function
makeLoll(unique(cm.a.2$GENE), cm.a.2) 
