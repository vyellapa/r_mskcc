##Consort diagram

#Load data:
#setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.clinical')
#source('clinical.data.import.R')

#Run: source('context.resistance.figure.R')

library(diagram)

par(mfrow=c(1,1), ps = 12, font = 1, mar=c(0,0,0,0))
##initialize new grphics device
openplotmat()
##number of elements per row
elpos<-coordinates (c(1,1,1,3))
##draw arrows from each row to next row
treearrow(from=elpos[1,],to=elpos[2,],lwd=3)  
treearrow(from=elpos[1,],to=elpos[3,],lwd=3)  
treearrow(from=elpos[3,],to=elpos[4:6,],lwd=3)  

##create a generic 3-lined label for each textbox
labels = vector(length=6)
labels[1] = paste('Pre-treated MM\n(n=', length(df$sampleID), ')', sep = '')
labels[2] = paste('Sample at progression\n(n=', sum(as.numeric(df$contextSample) == 2), ')', sep = '')
labels[3] = paste('Progression on treatment\n(n=', sum(as.numeric(df$progressingSampl) <= 5), ')', sep = '')
labels[4] = paste('IMID\n Len maint (n=', sum(as.numeric(df$progressingSampl) == 2), ')\nOther (n=', sum(as.numeric(df$progressingSampl) == 1), ')', sep = '')
labels[5] = paste('IMID and PI\n(n=', sum(as.numeric(df$progressingSampl) == 4), ')', sep = '')
labels[6] = paste('PI\n(n=', sum(as.numeric(df$progressingSampl) == 3), ')', sep = '')

##plot text boxes
for ( i in 1:3) textround (elpos[i,],radx=0.14,rady=0.05,lab=labels[i], shadow.size = 0)
for ( i in 4:6) textround (elpos[i,],radx=0.09,rady=0.07,lab=labels[i], shadow.size = 0)
