#Compatison of baseline and sample data
#Run: source('comparison.baseline.sample.R')
setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.clinical')

source('clinical.data.import.R')

library(gplots)
library(RColorBrewer)

#Make clean dataset
base <- select(df,
               sampleID,
               BMPCBase,
               caBase, 
               creaBase, 
               hgbBase, 
               WBCBase, 
               ANCBase,
               plateletsBase,
               beta2mBase, 
               albBase
               )

sampl <- select(df,
                sampleID,
               BMPCSampl,
               caSampl, 
               creaSampl, 
               hgbSampl, 
               WBCSampl, 
               ANCSampl,
               plateletsSampl,
               beta2mSampl, 
               albSampl
)

names <- c('sampleID',
           'BMPC',
            'Ca', 
            'Crea', 
            'Hgb', 
            'WBC', 
            'ANC',
            'Plate',
            'Beta2m', 
            'Alb'
)

colnames(base) <- names
colnames(sampl) <- names

#Make matrix of change from diagnosis to sample
change <- cbind(base[1],round(sampl[-1]/base[-1],2))
change$na_count <- apply(change, 1, function(x) sum(is.na(x)))
change <- filter(change, na_count <=4)
rownames(change) <- change$sampleID
change <- change[,-c(1,11)]
change <- change-1
change <- as.matrix(change)

#Make heatmap
#heatmap.2(change,
#        col=brewer.pal(11, 'RdYlBu'),
#       trace = 'none',
#        density.info = 'density',
#        key.title = 'Change from baseline',
#        key.xlab = 'Fold change',
#        scale = 'row',
#        ylab = 'SampleID'
#        )

#Combined clean dataset for figures
#base <- mutate(base, status = 'base')
#sampl <- mutate(sampl, status = 'sampl')
#combined <- rbind(base, sampl)
