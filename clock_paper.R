library(readr) 
library(tidyr)
library(splitstackshape)
library(dplyr)
library(readxl)
library(magrittr)
library(lmerTest)
library(knitr)
library(ggplot2)
library(reshape2)
library(MASS)
library(lme4)
library(pander)
library(survival)
library(MCMCpack)
library(extraDistr)
options("scipen"=7, "digits"=4)
library("seqinr")
library("Rsamtools")
library("GenomicRanges")
library("pheatmap")
library("RColorBrewer")
maxiter = 2000
options("scipen"=50, "digits"=4)
library(stringr)

 
####################################################
### signature contribution for each DP cluster
#################################################### 
setwd("/Users/yellapav/Desktop/record/clock_paper/data")
sig<- read.delim("all_patients_DP_contribution.txt", sep="\t", stringsAsFactors = F)
sig$sampleID<- rownames(sig)
sig_2<- cbind.data.frame(sig$sampleID, sig$Signature.Subs.05  * sig$mutations) #### file with signature 5 contribution for each cluster
colnames(sig_2)<-c("branches","burden")

### the columns 4th-6th may not be corrected in absolute number -  Kevin re-ran DP on Sanger data generating slightely different solution 
tree2<- read.delim("summary_fro_timing.txt", sep="\t", header=T, stringsAsFactors = F) 
tree<- tree2

#### V0D58I has 2 clusters that overlap >90%. I considered them as trunk collapsing these 2 clusters in one. 
### @Im assuming this clusters 2 and 1 are now the trunk?

tree$tree_structure[tree$patient == "V0D58I"]<-"2>1"
sig_2$burden[sig_2$branches == "V0D58I_2"]<- 1447.15 + 693.27
sig_2<- sig_2[-which(sig_2$branches == "V0D58I_3"),]


out <- do.call("rbind",strsplit(tree$tree_structure,"") )
out_df<- as.data.frame(out)
out_df[,1:4]<- apply(out_df[,1:4], 2, function(x) paste(tree$patient, x, sep="_"))
tree$clonal_five<- NA
tree$burden<- NA

for(i in (1:nrow(tree)))
{
  vec<-  as.data.frame(strsplit(tree$tree_structure[i],">"))
  vec2<- as.numeric(as.character(vec[,1]))
  names_branch<- paste(tree$patient[i], vec2, sep="_")
  sig_bran<- sig_2[which(sig_2$branches %in% names_branch),] #### file sig_2 with signature 5 contribution for each cluster
  tree$burden[i]<-  sum(sig_bran$burden)
  tree$clonal_five[i] <- sig_bran[sig_bran$branches == names_branch[1],2]
}

final<- tree
head(final) ###### final contains branched samples


### create final matrix for clock analysis follwoing Peter Campbell code (Mitchel, et al. Cell 2018)
clin<- read.delim("Supplementary_Table_1.txt", sep="\t")
colnames(clin)[1]<-"sample_ID"
def<- merge(final, clin[,1:3], by="sample_ID")
# def$burden<- rowSums(def[,12:15])
def2<- def[,c("sample_ID","patient","patient","tree_structure","Num.clonal","Num.branch.mutations","burden","clonal_five","Age","code_btranch")]

colnames(def2)<- c("Paz","Sample","LRISample","Name","Num.clonal.all","Num.branch.mutations","Num.mutations","Num.clonal" ,"Age","code_btranch" )

branch.df3 <- def2[def2$code_btranch =="branch",]  #### only multiple sample (chapman samples without dirichlet)
branch.df2<- def2[complete.cases(def2$Age),] #### only multiple sample (chapman samples without dirichlet)



##### start the clock process
ggplot(data = branch.df2, aes(x=Age, y=Num.mutations, color=LRISample)) +
  geom_point(shape=16) + xlim(c(0,80)) + ylim(c(0,10000)) +
  xlab("Age (years)") + ylab("Number of mutations per subclone") ##  ylab("Number of clonal mutations")

### 2 options:
###   1) runnig all samples == only way is to use linear mixed effects models costraining the origin
###   2) runnig only branched samples == having multiple samples you can run linear mixed effects models with unconstrained origin and with variable slope

#### run only patients with multiple samples

branch.df<- branch.df3 ###### only sanger with dirichlet data

### random effect (Age | Sample) and fixed effect (Age + ) - on the slope which is unconstrain
muts.per.year.lmer.with.pt.intercepts <- lmer(Num.mutations ~ Age + (Age | Sample), data=branch.df, REML=FALSE)
### variable slope without costraining to the origin
muts.per.year.lmer.with.pop.intercept <- lmer(Num.mutations ~ Age -1 + (Age | Sample), data=branch.df, REML=FALSE)
### costrained to pass throught zero (the origin)
muts.per.year.lmer <- lmer(Num.mutations ~ Age - 1 + (Age - 1 | Sample ), data=branch.df, REML=FALSE)




