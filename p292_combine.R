library(openxlsx)
library(dplyr)

tra <- read.xlsx('./p292.SNParray.myTYPE.unique.Feb26.xlsx', sheet = 6)
names(tra)[1] <- 'Pathology.ID'

tra <- tra %>%
                select(Pathology.ID, translocation) %>%
                filter(translocation != 't(8;14)')

clin <- read.xlsx('./p292.clinical.data.xlsx')

comb <- left_join(clin, tra)

write.xlsx(comb, "./p292_clinical_molecular.xlsx")