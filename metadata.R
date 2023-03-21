library(openxlsx)
library(dplyr)

clin <- read.xlsx("./HOTB_clinical_combined.xlsx", detectDates = T)
key <- read.xlsx("./HOTB_all_key_2.xlsx", detectDates = T)

out <- clin %>%
                left_join(key %>%
                                  select(sample_ID, cell_type) %>%
                                  mutate(sample_ID = as.character(sample_ID))) %>%
                filter(!LEUKGEN_ID %in% c(NA, "")) %>%
                select(LEUKGEN_ID, date_sample, disease_stage, cell_type)

out <- out %>%
        mutate(disease_stage = ifelse(disease_stage == 'NDMM', 'Newly diagnosed multiple myeloma',
                                      ifelse(disease_stage == 'RRMM', 'Relapsed/refractory multiple myeloma',
                                             ifelse(disease_stage == 'SMM', 'Smoldering multiple myeloma', 'Multiple myeloma during active treatment')))) %>%
        mutate(cell_type = ifelse(cell_type == 'MNC', 'Bone marrow mononuclear cells', 'Bone marrow CD138 positive cells'))

write.xlsx(out, "p222_sample_meta.xlsx")