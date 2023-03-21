library(dplyr)
library(stringr)
prior1 = 0.5
p01 = p10 = 0.1

##################################################################################################################
##### Please change file path here if trying it on test data ############################################
sham = read.table("~/Downloads/tiny_example.sham", header = F, sep = "\t", colClasses = c("integer","character"))
##################################################################################################################

# Calculate the frequency of 1's and coverage at each base
mn = min(sham$V1)
sham = sham %>% 
  dplyr::mutate(last.pos = V1 + str_length(V2)) 
mx = max(sham$last.pos)
priors1 = as.numeric(rep(prior1,(mx-mn+1)))
priors0 = 1 - priors1


frequency = as.numeric(rep("na",(mx-mn+1)))
priorl = vector(mode = "list", length = dim(sham)[1])
freql = vector(mode = "list", length = dim(sham)[1])
#For each for of sham
for(i in 1:dim(sham)[1]) {
  freq.temp = frequency
  chars = unlist(str_split(sham[i,2],"")) #fetch read string into a vector
  start.pos = sham[i,1]
  
  for(j in seq_along(chars)) { #iterate over each position
    if(chars[j]=="1") {
      freq.temp[start.pos+j]=1
    } else {
      freq.temp[start.pos+j]=0
      #priors1.temp[start.pos+j]=(1-priors1.temp[start.pos+j])*p01}
    } 
  }
  
 freql[[i]] = freq.temp
}


freq1 = as.data.frame(t(matrix(unlist(freql), byrow = T, ncol=length(frequency))))
freq2 = freq1
freq2$frequency1=prior1
freq2$coverage = 0
f = as.matrix(freq1)
f[f == 0 ] <- 1
freq3 = as.data.frame(f)
for(i in 1:dim(freq1)[1]) {
  freq2[i,]$frequency1 = sum(as.vector(freq1[i,]),na.rm=TRUE)
  freq2[i,]$coverage = sum(as.vector(freq3[i,]),na.rm=TRUE)
}



## Assuming binomial distribution, 
## use the sucessses (1's observed)  and trials (covearge from reads) at each particular site to caluclate the posteriors 
## dbinom package in R is used for this purpose
## https://www.rdocumentation.org/packages/stats/versions/3.3/topics/Binomial

freq2 %>% 
  dplyr::mutate(posterior = dbinom(frequency1, coverage, 0.1)) %>%
  dplyr::mutate(posterior=ifelse(coverage==0,prior1,posterior)) %>% 
  dplyr::mutate(posterior2=1-posterior)
                                                                                                          

