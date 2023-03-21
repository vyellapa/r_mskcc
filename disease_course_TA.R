setwd("C:/Users/admin/Documents/R")
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)

patient <- read.csv('180510_clinical.csv', stringsAsFactors = F)

patient$Data <- factor(patient$Data)
patient$Type_nr <- factor(patient$Type_nr)
patient$Datatype <- factor(patient$Datatype)
patient$Patient <- factor(patient$Patient)
patient$Date_start <- mdy(patient$Date_start)
patient$Date_end <- mdy(patient$Date_end)

patient <- mutate(patient, treatstart = min(Date_start))

patient$Month_start <- month(as.period(interval(patient$treatstart, patient$Date_start), unit = 'month'))
patient$Month_end <- month(as.period(interval(patient$treatstart, patient$Date_end), unit = 'month'))
patient <- mutate(patient, Month_interval = Month_end - Month_start)
patient <- mutate(patient, y=5)
patient <- mutate(patient, y=ifelse(patient$Datatype=='Treatment', 5, NA))


ggplot(data=filter(patient, Datatype=='Treatment'), 
       aes(x=Month_start, y=y))+
  geom_segment(aes(x=Month_start, xend=Month_end, y=y, yend=y, color=Data), size=5) + 
  theme(legend.position = "bottom")












###################################################################################################
#Don't need this(However would be useful for plotting two different dataframes together):

start <- patient[c("Patient", "Datatype", "Type_nr", "Date_start", "Data", "Date_end", 
                   "Best_response", "Imaging", "BM", "treatstart", "Month_start", "Month_interval", "y")]

end <- patient[c("Patient", "Datatype", "Type_nr", "Date_start", "Data", "Date_end", 
                 "Best_response", "Imaging", "BM", "treatstart", 
                 "Month_end", "Month_interval", "y")]

start <- filter(start, Datatype=='Treatment')
end <- filter(end, Datatype=='Treatment')

start <- rename(start, month=Month_start)
end <- rename(end, month=Month_end)


#grouped by Data
ggplot(data = start, aes(x = month, y = y), group = Data, shape = Data) +
  geom_segment(data = merge(start, end, by = 'Data'), 
               aes(x = month.x, xend = month.y, y = y.x, yend = y.y, color=Data), size=5) +
  theme(legend.position="bottom")

#grouped by Type_nr
ggplot(data = start, aes(x = month, y = y), group = Type_nr, shape = Type_nr) +
  geom_segment(data = merge(start, end, by = 'Type_nr'), 
               aes(x = month.x, xend = month.y, y = y.x, yend = y.y, color=Type_nr), size=5) +
  theme(legend.position="bottom")
