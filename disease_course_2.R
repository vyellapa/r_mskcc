library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)

patient <- read.csv('../180510_clinical.csv', stringsAsFactors = F)

clinical <- read.csv('../180510_clinical_baseline.csv', stringsAsFactors = F)
clinical <- select(clinical, Patient, Diagn_date, Death_date)

patient <- left_join(patient, clinical, by = 'Patient')

patient$Data <- factor(patient$Data)
patient$Type_nr <- factor(patient$Type_nr)
patient$Datatype <- factor(patient$Datatype)
patient$Patient <- factor(patient$Patient)
patient$Date_start <- mdy(patient$Date_start)
patient$Date_end <- mdy(patient$Date_end)
patient$Diagn_date <- mdy(patient$Diagn_date)
patient$Death_date <- mdy(patient$Death_date)

patient$Month_start <- month(as.period(interval(patient$Diagn_date, patient$Date_start), unit = 'month'))
patient$Month_end <- month(as.period(interval(patient$Diagn_date, patient$Date_end), unit = 'month'))
patient$Month_death <- month(as.period(interval(patient$Diagn_date, patient$Death_date), unit = 'month'))
patient <- mutate(patient, Month_interval = Month_end - Month_start)


ggplot(data=filter(patient, Datatype=='Treatment'), 
       aes(x=Month_start, y=Patient))+
        geom_segment(aes(x=Month_start, xend=Month_end, y=Patient, yend=Patient, color=Data), 
                     size=5)+
        geom_point(data = filter(patient, Datatype=='Progression'),
                   aes(x = Month_start, y = Patient, shape = 'Progression'),
                   size = 2.5)+
        geom_point(data = filter(patient, Datatype=='Progression' & Type_nr == 1),
                   aes(x = Month_death, y = Patient, shape = 'Death'),
                   size = 4)+
        theme(legend.position = "bottom")+
        labs(y = 'Patient',
             x = 'Time from diagnosis of Multiple Myeloma (Months)',
             col = 'Treatment')

ggsave('Rapid_autopsy_patients.png',
       width = 15,
       height = 6)