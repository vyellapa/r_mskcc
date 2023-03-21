library(openxlsx)
library(dplyr)
library(reshape2)

pts <- read.xlsx("./patients.xlsx")
pts <- pts[-4]

snvs <- read.csv("./p222_all_snv_indel_newcallers.csv", stringsAsFactors = F)
cyto <- read.csv("./p222_all_fishseq.csv", stringsAsFactors = F)

key <- read.xlsx("./HOTB_clinical_combined.xlsx")


myTYPE_data <- key %>%
        select(MRN, LEUKGEN_ID) %>%
        rename(LEUKID = "LEUKGEN_ID") %>%
        filter(!is.na(LEUKID) & LEUKID != 0) %>%
        left_join(cyto %>%
                          filter(SEQ_result == 1) %>%
                          select(LEUKID, Cytogenetic_abnormality) %>%
                          dcast(LEUKID~Cytogenetic_abnormality, length, fill = 0)) %>%
        left_join(snvs %>%
                          mutate(MUTATIONS = paste(VAG_GENE, VAG_PROTEIN_CHANGE, VAG_cDNA_CHANGE, TARGET_VAF, MANUAL_ANNOTATION, sep = '_')) %>%
                          select(LEUKID, MUTATIONS) %>%
                          group_by(LEUKID) %>%
                          summarise(MUTATIONS = paste(MUTATIONS, collapse = "|")))


out <- left_join(pts, myTYPE_data)        

write.xlsx(out, "./mytype_p222_selected.xlsx")