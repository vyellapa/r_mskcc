##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpBattenberg.
#
# cgpBattenberg is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

args=(commandArgs(TRUE))
lib_path<-toString(args[1])
imputeInfoFile<-toString(args[2])
inputFile<-toString(args[3])
normalInputFile<-toString(args[4])
outputDir<-toString(args[5])
minCount = as.numeric(args[6])
samplename = toString(args[7])
thougenloc<-toString(args[8])

source("/ifs/work/leukgen/home/gg10/soft/complete_genomics/GetBAFLogR.R")
source(paste(lib_path,"ascat.R",sep="/"))

getBAFsAndLogRs(imputeInfoFile,inputFile,normalInputFile,paste(outputDir,"normalBAF.tab",sep="",collapse=""),paste(outputDir,"mutantBAF.tab",sep="",collapse=""),paste(outputDir,"normalLogR.tab",sep="",collapse=""),paste(outputDir,"mutantLogR.tab",sep="",collapse=""),minCount=minCount,samplename=samplename,thougenloc=thougenloc)

# This is a workaround such that BB doesnt have to be rerun across a lot of samples. The problem is that when ASCAT reads in a file it changes two things: If the samplename starts with a number it places an X in front of the name and underscores are replaced by dots. ascat.loadData was extended to take the real samplenames and place those in the header of the read in data files. However that changes the other underscores and dots thing as well and that changes filenames. We therefore place the dots in the samplenames such that ASCAT will produce the file with dots in it. A complete fix will be added at a later date.
samplename = gsub("\\-", "\\.", samplename)
ascat.bc = ascat.loadData(paste(outputDir,"mutantLogR.tab",sep=""),paste(outputDir,"mutantBAF.tab",sep=""),paste(outputDir,"normalLogR.tab",sep=""),paste(outputDir,"normalBAF.tab",sep=""), samplenames=c(samplename))
ascat.plotRawData(ascat.bc, parentDir=outputDir)

q(save="no")
