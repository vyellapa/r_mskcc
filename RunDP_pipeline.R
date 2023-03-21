#'
#' DPClust 1D/nD pipeline v0.1.1
#' 
#' sd11@sanger.ac.uk
#'
# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/RunDP_pipeline.R /ifs/work/leukgen/home/gg10/soft 1 10000 1000 /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/nd_dp_clustering_input_current/nd_dp_input /ifs/res/papaemme/users/gg10/renal/sdhb_mutant/nd_dp_clustering_input_current/sample2purity.txt nd_dp F 1 2 NA

setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/dirichlet_step1")
outdir = getwd()
args = commandArgs(TRUE)
libdir = toString(args[1]) # Directory where the pipeline is installed
run = as.integer(args[2]) # The sample to be run. Integer that selects from a list of unique samplenames
no.iters = as.integer(args[3]) # Number of iters
no.iters.burn.in = as.double(args[4]) # Number of iters used for burnin
datpath = toString(args[5]) # Where are input files stored
purity_file = toString(args[6]) # A file containing samplenames and purity
analysis_type = toString(args[7]) # String that represents the analysis that is to be run
parallel = as.logical(args[8]) # Supply true or false whether to run parts of the method in parallel
no.of.threads = as.integer(args[9]) # Integer that determines how many threads to use when running parts in parallel
mut.assignment.type = as.integer(args[10]) # Integer that determines which mutation assignment method is to be used in 1d/nd cases
num_muts_sample = as.integer(args[11]) # Integer that determines how many mutations to sample

# Optional arguments
if (length(args) >= 12) {
  bin.size = as.double(args[12]) # Binning of mutations not suitable for 1D/nD.
  if (length(args) >= 13) {
    blockid = as.integer(args[13]) # Select a seed
    if (length(args) >= 14) {
      no.of.blocks = as.integer(args[14]) # Not used by 1D/nD clustering.
    } else {
      no.of.blocks = 1
    }
  } else { 
    blockid = 1
    no.of.blocks = 1
  }
} else {
  bin.size = NA
  blockid = 1
  no.of.blocks = 1
}

# Hard coded for now
is.male = T
is.vcf = F

# Set the seed. If it's a block run we use the blockid to select a different seed
seeds = c(123,321,213,231,456,654,465,645,789,987,978,798)
set.seed(seeds[blockid])

# Check whether a supported analysis_type was supplied
supported_commands = c('1d_dp', 'nd_dp', 'replot_1d', 'replot_nd', 'sample_muts', 'reassign_muts_1d')
if (!(analysis_type %in% supported_commands)) {
  print(paste("Type of analysis", analysis_type, "unknown."))
  print(paste(c("Specify either ", supported_commands)), sep=" ")
  q(save="no", status=1)
}

# Check whether the mut.assignment.type is supported
supported_mut.assignment.methods = c(1,2)
if (!(mut.assignment.type %in% supported_mut.assignment.methods)) {
  print(paste("Type of mutation assignment method", mut.assignment.type, "unknown."))
  print(paste(c("Specify either ", supported_mut.assignment.methods)), sep=" ")
  q(save="no", status=1)
}

# Source the required files
setwd(libdir)
source("RunDP.R")
setwd(outdir)

library(DPClust)

# Parse the input file and obtain the required data for this run
sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
if(analysis_type=='1d_dp') {
	samplename = unique(sample2purity$sample)[run]
	datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
	subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
	cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
} else if (analysis_type=='nd_dp') {
	samplename = unique(sample2purity$sample)[run]
	datafiles = sample2purity$datafile
	subsamples = sample2purity$sample
	cellularity = sample2purity$cellularity
}
print("")
print(paste("Running:", samplename, sep=" "))
print(paste("Working dir:", outdir, sep=" "))
print(paste("Analysis type:", analysis_type, sep=" "))
print("Datafiles:")
print(datafiles)
print("")

# Set the name of the output directory
outdir = paste(outdir, "/", samplename, "_DPoutput_", no.iters,"iters_",no.iters.burn.in,"burnin", sep="")

# Create the output directory
if (!file.exists(outdir)) { dir.create(outdir) }

if (file.exists(paste(outdir, "/dataset.RData", sep=""))) {
  # Wait a certain number of seconds before loading - this is required for when starting multiple processes/threads on this file
  if (!is.na(blockid) & blockid != "NA") {
    Sys.sleep(blockid*3)
  }
	load(paste(outdir, "/dataset.RData", sep=""))
} else {
  list_of_datafiles = paste(datpath, datafiles, sep="/")
  
	dataset = load.data(list_of_datafiles, 
                      cellularity=cellularity, 
                      Chromosome="chr", 
                      position="end",
                      WT.count="WT.count", 
                      mut.count="mut.count", 
                      subclonal.CN="subclonal.CN", 
                      no.chrs.bearing.mut="no.chrs.bearing.mut", 
                      mutation.copy.number="mutation.copy.number", 
                      subclonal.fraction="subclonal.fraction", 
  		                is.male=is.male,
                      is.vcf=is.vcf,
  		                ref.genome.version="hg19")

  print(num_muts_sample)
  print(class(num_muts_sample))
  if (!is.na(num_muts_sample) & num_muts_sample!="NA") {
    dataset = sample_mutations(dataset, num_muts_sample)
  }
}

# Save the dataset
save(file=paste(outdir, "/dataset.RData", sep=""), dataset)

RunDP(analysis_type=analysis_type, 
      dataset=dataset, 
      samplename=samplename, 
      subsamples=subsamples, 
      no.iters=no.iters, 
      no.iters.burn.in=no.iters.burn.in, 
      outdir=outdir, 
      conc_param=0.01, 
      cluster_conc=5, 
      resort.mutations=T, 
      parallel=parallel, 
      blockid=blockid, 
      no.of.blocks=no.of.blocks, 
      remove.node.frequency=12,
      remove.branch.frequency=51,
      mut.assignment.type=mut.assignment.type,
      annotation=vector(mode="character",length=nrow(dataset$mutCount)),
      init.alpha=0.01, 
      shrinkage.threshold=0.1,
      bin.size=bin.size,
      muts.sampled=!is.na(num_muts_sample))
