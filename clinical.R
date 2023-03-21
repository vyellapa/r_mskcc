library(dplyr)
library(tidyverse)

treat = read.table("~/Desktop/p220_2019/misc/clinical_unix.txt",header=T,sep="\t")
treat = treat %>% dplyr::filter(Date_end!="")


cc = c("Patient","Datatype","Type_nr","Date_start","Data","Update","Date_end","Best_response","Imaging","BM")
cc = c("Patient","Datatype","Date_start","Data","Date_end")
treat=treat[,cc] %>% unique() 
colnames(treat) = c("Patient","Datatype","start","value","end")
treat = treat %>% dplyr::filter(Patient == "RA17-7" | Patient == "RA16-14" | Patient == "RA16-4" | Patient == "RA16-1" )
write.table(treat,"~/Desktop/p220_2019/misc/curated.txt",row.names = F,col.names = TRUE,append = F, quote = F,sep = "\t" )

new_t = read.table("~/Desktop/p220_2019/misc/curated_death.txt",header=T,sep="\t")
#new_t %>% filter(Datatype!="Death")
t1 = new_t %>% filter(Patient == "RA16-1") %>% 
  dplyr::mutate(std=(as.Date(as.character(start), format="%m/%d/%Y") - as.Date(as.character("3/23/06"), format="%m/%d/%Y")), 
                Patient = "I106917",
                endd=(as.Date(as.character(end), format="%m/%d/%Y") - as.Date(as.character("3/23/06"), format="%m/%d/%Y")))

t2 = new_t %>% filter(Patient == "RA16-4") %>% 
  dplyr::mutate(std=(as.Date(as.character(start), format="%m/%d/%Y") - as.Date(as.character("3/4/13"), format="%m/%d/%Y")),
                Patient = "I130718",
                endd=(as.Date(as.character(end), format="%m/%d/%Y") - as.Date(as.character("3/4/13"), format="%m/%d/%Y")))

t3 = new_t %>% filter(Patient == "RA16-14") %>% 
  dplyr::mutate(std=(as.Date(as.character(start), format="%m/%d/%Y") - as.Date(as.character("4/10/07"), format="%m/%d/%Y")),
                Patient = "I130719",
                endd=(as.Date(as.character(end), format="%m/%d/%Y") - as.Date(as.character("4/10/07"), format="%m/%d/%Y")))

t4 = new_t %>% filter(Patient == "RA17-7") %>% 
  dplyr::mutate(std=(as.Date(as.character(start), format="%m/%d/%Y") - as.Date(as.character("8/6/15"), format="%m/%d/%Y")),
                Patient = "I130720",
                endd=(as.Date(as.character(end), format="%m/%d/%Y") - as.Date(as.character("8/6/15"), format="%m/%d/%Y")))


new_tt=rbind(t1,t2)
new_tt=rbind(new_tt,t3)
new_tt=rbind(new_tt,t4)
new_tt = new_tt %>% dplyr::mutate(start=(as.numeric(std))/30, end=(as.numeric(endd))/30) %>% dplyr::mutate(start=ifelse(start<0,0,start))

#new_t$date_diff <- (as.Date(as.character(new_t$end), format="%m/%d/%Y")-
#  as.Date(as.character(new_t$start), format="%m/%d/%Y"))





m=melt(treat, id.vars = c("Patient","value"))
colnames(m) = c("Patient","task","variable","value")
m$value=as.Date(m$value, format = "%m/%d/%Y")


ggplot(m, aes(value, task)) + facet_wrap(Patient~.,nrow = 4,scales = "free_y")+
  geom_line(size = 10) +
  labs(x = '', y = '', title = 'Gantt chart using ggplot2') +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(colour="black", linetype = "dashed"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels = "%Y %b", limits = c(start_date, NA), date_breaks = '6 months')

new_tt = new_tt %>% dplyr::select(-c("Datatype","std","endd"))
m=melt(new_tt, id.vars = c("Patient","value")) 
colnames(m) = c("Patient","task","variable","value")
m = m %>% dplyr::arrange(value)
m=m %>% dplyr::mutate(task=ifelse(task=="Death","Disease Course",as.character(task)))
m$task <- factor(m$task, levels=rev(c(grep("Disease",unique(m$task),invert=TRUE,value=TRUE),"Disease Course")))


ggplot(m, aes(value, task)) + facet_wrap(Patient~.,nrow = 4,scales = "free_y")+
  geom_line(size = 6, color="dodgerblue") +
  labs(x = 'Time in Months', y = '', title = 'Clinical Course') +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(colour="black", linetype = "dashed"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0),
        panel.spacing = unit(0.01, "lines"),
        strip.background =element_rect(fill="#FFFFFF"))

m$color="a"
m[grep("MEL",m$task),]$color="b"
m[grep("PAC",m$task),]$color="c"

ggplot(m, aes(value, task,color=factor(color))) + facet_grid(Patient~.,scales = "free_y",space="free")+
  geom_line(size = 6.5 ) + scale_colour_manual(values = c("dodgerblue", "gray", "#fe0000")) +
  labs(x = 'Time in Months', y = '', title = 'Clinical Course') +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(colour="black", linetype = "dashed"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0),
        panel.spacing = unit(0.01, "lines"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 39),
        strip.background =element_rect(fill="#FFFFFF"))

