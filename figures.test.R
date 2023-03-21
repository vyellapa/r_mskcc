ggplot(df, aes(x=factor(1), fill=refractoryFirstDrugs))+
  geom_bar(width = 1)+
  coord_polar("y")+
  facet_wrap(~progressingSampl,nrow = 3)

ggplot(df, aes(x=factor(1), fill=refractoryFirstDrugs))+
  geom_bar(width = 1)
  
ggplot(df, 
       aes(
 x = timeSampleTreatment,
  fill = contextSample))+
  geom_bar()+
  scale_color_brewer(
    palette = 'Accent')+
  coord_flip()



####Identify patients who started a new treatment line within few months of sample
drop <- grep('cont', df$treatmentPostSampleRegimen, fixed = TRUE)
newline <- select(df, sampleID, timeSampleTreatment, treatmentPostSampleRegimen, contextSample, onTreatmentSample)[-drop,]
sameline <- select(df, sampleID, timeSampleTreatment, lastTreatmentSamplResponse, treatmentPostSampleRegimen, firstDrugsResponse, contextSample, onTreatmentSample)[drop,]