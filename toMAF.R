#### Function to convert Caveman output to MAF format
## Takes a data frame containing the original columns of Caveman output 'Project 222'
## Subsets most important columns for MAF and transforms to appropriate data format for Bioconductor package 'maftools'
## Store output as MAF object
## v.1.17.2018

#### About MAF
## Mandatory fields for MAF (Mutation Annotation Format): 
## Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
## Recommended optional fields: non MAF specific fields containing vaf and amino acid change information.
## Complete specification: https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

#### About maftools
## Link https://www.bioconductor.org/packages/3.7/bioc/vignettes/maftools/inst/doc/maftools.html

#### Dependencies
library(dplyr)
library(maftools)

#### Define function
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
