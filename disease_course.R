library(dplyr)
library(ggplot2)
library(lubridate)

patient <- read.csv('180510_clinical.csv', stringsAsFactors = F)

patient$Datatype <- factor(patient$Datatype)
patient$Patient <- factor(patient$Patient)
patient$Date_start <- mdy(patient$Date_start)
patient$Date_end <- mdy(patient$Date_end)

patient <- mutate(patient, treatstart = min(Date_start))

patient$Month_start <- month(as.period(interval(patient$treatstart, patient$Date_start), unit = 'month'))
patient$Month_end <- month(as.period(interval(patient$treatstart, patient$Date_end), unit = 'month'))
patient <- mutate(patient, Month_interval = Month_end - Month_start)

ggplot()+
        geom_step(data = filter(patient, Datatype == 'Treatment'),
                  aes(Month_start, desc(Type_nr)))+
        geom_text(data = filter(patient, Datatype == 'Treatment'),
                  aes(Month_start, desc(Type_nr), label = Data),
                  hjust = -0.1,
                  vjust = -1)+
        geom_point(data = filter(patient, Datatype == 'Progression'),
                   aes(Month_start, desc(Type_nr)),
                   size = 3,
                   col = 'red')+
        geom_label(data = filter(patient, Datatype == 'Progression'),
                   aes(Month_start, desc(Type_nr), label = Data),
                  hjust = 0.4,
                  vjust = 2)+
        theme_bw()
        
