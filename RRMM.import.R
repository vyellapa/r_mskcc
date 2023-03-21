## Import Project 222 subset "RRMM" datafiles from caveman, pindel, delly, brass and cnvkit. 
## v.1.3.2018

#run setwd('/Users/rustade/Documents/MSKCC/RRMM/analysis.seq'); source('/Users/rustade/Documents/MSKCC/RRMM/analysis.seq/RRMM.import.R')

#Workdir
setwd('/Users/rustade/Documents/MSKCC/Projects/222/RRMM/analysis.seq')

#Import data
caveman <- read.csv('../datafiles/p222.results.RRMM/caveman.csv', stringsAsFactors = FALSE)
pindel <- read.csv('../datafiles/p222.results.RRMM/pindel.csv', stringsAsFactors = FALSE)
delly <- read.csv('../datafiles/p222.results.RRMM/delly.csv', stringsAsFactors = FALSE)
brass <- read.csv('../datafiles/p222.results.RRMM/brass.csv', stringsAsFactors = FALSE)
cnvkit <- read.csv('../datafiles/p222.results.RRMM/cnvkit.csv', stringsAsFactors = FALSE)

