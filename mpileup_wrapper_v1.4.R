# R wrapper script for mpileup -> nucleotide frequencies
# mpileup_wrapper_v1.4.R
# version: 1.4
# Komal S Rathi
# Memorial Sloan Kettering Cancer Center
# Created ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified ##------ Mon Feb 15 17:12:17 2016 ------##
# Usage:
# Rscript mpileup_wrapper_v1.4.R --help
# Changes in this version:
# takes position range as well as list of positions
# take input directory or file with list of bam files
# writes pileup output for samples with 0 depth
# writes pileup output for positions with 0 coverage
# takes commandline parameters
# computes counts for overlapping regions
# allows subsetting output to get only depth
# adds multi threading capability

timestamp()
suppressMessages(library(Rsamtools, quietly = TRUE))
suppressMessages(library(plyr, quietly = TRUE))
suppressMessages(library(reshape2, quietly = TRUE))
suppressMessages(library(optparse, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))
suppressMessages(library(tools, quietly = TRUE))
suppressMessages(library(doMC, quietly = TRUE))

option_list = list(make_option(c("-i", "--input"), type = "character",
                               help = "input directory, sample list file or single bam file",
                               metavar = "path", action = "store"),
                   make_option(c("-p", "--pos"), type = "character",
                               help = "positions file", metavar = "path"),
                   make_option(c("-f", "--fasta"), type = "character",
                               help = "reference fasta", metavar = "path"),
                   make_option(c("-s", "--subset"), type = "logical",
                               default = FALSE,
                               help = "get a subset of the output? [default %default]",
                               metavar = "logical"),
                   make_option(c("-o","--outdir"), type = "character",
                               help = "output directory path",
                               metavar = "path"),
                   make_option(c("-d", "--depth"), type = "integer",
                               help = "max depth [default %default]",
                               default = 1000, metavar = "integer"),
                   make_option(c("-bq", "--minbq"), type = "integer",
                               help = "min base quality [default %default]",
                               default = 0, metavar = "integer"),
                   make_option(c("-mq", "--minmq"), type = "integer",
                               help = "min mapq [default %default]",
                               default = 0, metavar = "integer"),
                   make_option(c("-ins", "--insertions"), type = "logical",
                               default = FALSE,
                               help = "include insertions? [default %default]",
                               metavar = "logical"),
                   make_option(c("-del", "--deletions"), type = "logical",
                               default = FALSE,
                               help = "include deletions? [default %default]",
                               metavar = "logical"),
                   make_option(c("-t","--threads"), type = "integer",
                               help = "threads to use [default %default]",
                               default = 1, metavar = "integer"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# get argument values
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Provide input file (list of bams) or input directory", call. = FALSE)
}

if (is.null(opt$pos)){
    print_help(opt_parser)
    stop("Provide positions file", call. = FALSE)
}

if (is.null(opt$fasta)){
    print_help(opt_parser)
    stop("Provide reference fasta", call. = FALSE)
}

if(is.null(opt$outdir)){
    print_help(opt_parser)
    stop("Provide output directory", call. = FALSE)
}

# get values of commandline args
positions_file <- opt$pos
input <- opt$input
fasta_file <- opt$fasta
max_depth <- opt$depth
min_base_quality <- opt$minbq
min_mapq <- opt$minmq
include_insertions <- opt$insertions
include_deletions <- opt$deletions
getsubset <- opt$subset
outdir <- opt$outdir
threads <- opt$threads

# register threads for ddply to run in parallel
registerDoMC(threads)

# read positions file
positions_file <- read.delim(positions_file, header = F)

# read fasta reference
fasta_file <- FaFile(file = fasta_file)

# if input is directory, process all bams in that directory
# else if input is file, process all bams listed in the text file
# you can also give path to a single bam file as input
if(dir.exists(input)){
    files <- list.files(path = input,
                        pattern = "*.sorted.bam$",
                        full.names = T,
                        recursive = T)
} else if(file.exists(input) && file_ext(input) == "txt"){
    files <- read.table(input, comment.char = "")
    files <- as.vector(files$V1)
} else if(file.exists(input) && file_ext(input) == "bam"){
    print("opening bam")
    files <- input
    print("done bam")
}

print ("mofo")
# format positions file and get reference base
if(ncol(positions_file) == 3 && nrow(positions_file) >= 1){
    # when you have chr, pos and tag list with ref and alt info
    print("positions file 1")
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))
    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
    names(positions_file) <- c('CHR','POS','TAG','REF')

} else if(ncol(positions_file) == 4 && nrow(positions_file) == 1){
    # range file only one range
    print("positions file 2")
    tags <- positions_file$V4
    positions_file <- data.frame(V1 = positions_file$V1,
                                 V2 = seq(from = positions_file$V2,
                                          to = positions_file$V3),
                                 V3 = as.character(positions_file$V4))
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))

    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
    names(positions_file) <- c('CHR','POS','TAG','REF')

} else if(ncol(positions_file) == 4 && nrow(positions_file) > 1){
    print("positions file 3")
    # range file with multiple ranges
    mylist <- apply(positions_file,1,
                    function(x) data.frame(V1 = as.numeric(x[1]),
                                           V2 = seq(from = as.numeric(x[2]),
                                                    to = as.numeric(x[3])),
                                           V3 = as.character(x[4])))
    positions_file <- do.call(rbind,
                              lapply(mylist,
                                     data.frame,
                                     stringsAsFactors=FALSE))
    tags <- positions_file$V3
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))
    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
    names(positions_file) <- c('CHR','POS','TAG','REF')

} else if(ncol(positions_file) == 5 && nrow(positions_file) >= 1){
    print("positions file 4")
    # when you have chr, pos and tag list without ref and alt info
    names(positions_file) <- c('CHR','POS','TAG','REF','ALT')
}

# function to compute pileup
compute.pileup <- function(x){
    param <- ScanBamParam(which = GRanges(x$CHR,
                          IRanges(start = as.numeric(x$POS),
                                  end = as.numeric(x$POS))))

    p_param <- PileupParam(distinguish_strand = TRUE,
                           distinguish_nucleotides = TRUE,
                           max_depth = 1000,
                           include_deletions = F,
                           include_insertions = F,
                           min_base_quality = 0,
                           min_mapq = 0,
                           min_nucleotide_depth = 0,
                           min_minor_allele_depth = 0)

    res <- pileup(bf, scanBamParam = param, pileupParam = p_param)

    # no pileup information so do the following:
    if(nrow(res) == 0){
        res <- data.frame(seqnames = x$CHR,
                          pos = x$POS,
                          which_label = paste(x$CHR, ":", x$POS, "-", x$POS,
                                              sep = ''))
        res$strand <- rep(c('-','+'), length.out = nrow(x))
        res$nucleotide <- rep(c('A','T','G','C'), length.out = nrow(x))
        res$count <- 0
        res <- res[,c(1,2,4,5,6,3)]
    }

    res <- merge(res, x,
                 by.x = c('seqnames', 'pos'),
                 by.y = c('CHR', 'POS'),
                 all.y = T)

    # factorize labels
    res$strand <- factor(res$strand, levels = c('-','+'))
    res$nucleotide <- factor(res$nucleotide, levels = c('A','T','G','C'))
    res$seqnames <- factor(res$seqnames, levels = unique(res$seqnames))

    # if the input does not contain REF ALT
    if(ncol(res) == 8){
        # long to wide format
        results <- dcast(res, seqnames+pos~nucleotide+strand,
                         value.var = "count", fill = 0, drop = FALSE)

        # add reference base
        results$REF <- x$REF
        results <- results[,c(1,2,ncol(results),3:(ncol(results)-1))]

        # remove columns that have NA in the name
        if(length(grep("NA$|^NA",colnames(results)))>0){
            results <- results[, -grep("NA$|^NA", colnames(results))]
        }

        # calculate depth - total, forward and reverse
        results$D <- apply(results[, 4:ncol(results)], 1, sum)
        print(str(results))
        results$D_forward <- apply(results[, grep('[+]', colnames(results))],
                               1, sum)
        results$D_reverse <- apply(results[, grep('[-]', colnames(results))],
                               1, sum)
        results <- results[, c(1, 2, 3,
                           grep('D', colnames(results)),
                           grep('[+]', colnames(results)),
                           grep('[-]', colnames(results)))]
    } else if(ncol(res) == 9){
        res$ALT <- as.character(res$ALT)
        res$REF <- as.character(res$REF)

        # long to wide format
        results <- dcast(res, seqnames+pos+REF+ALT~nucleotide+strand,
                         value.var = "count", fill = 0, drop = FALSE)

        print(str(results))
        # remove columns that have NA in the name
        if(length(grep("NA$|^NA",colnames(results)))>0){
            results <- results[, -grep("NA$|^NA", colnames(results))]
        }

        # calculate depth - total, forward and reverse
        results$D <- apply(results[, 5:ncol(results)], 1, sum)
        results$D_forward <- apply(results[, grep('[+]', colnames(results))],
                               1, sum)
        results$D_reverse <- apply(results[, grep('[-]', colnames(results))],
                               1, sum)
        results <- results[, c(1, 2, 3, 4,
                           grep('D', colnames(results)),
                           grep('[+]', colnames(results)),
                           grep('[-]', colnames(results)))]

    }
    colnames(results) <- sub('[+]', 'forward', colnames(results))
    colnames(results) <- sub('[-]', 'reverse', colnames(results))
    colnames(results)[1:2] <- c('CHR', 'POS')

    return(results)
}

# get pileup for each file
for(i in files){
    print("here we go")
    print(paste("Processing...", basename(i), sep = ''))
    bamfile <- i
    bf <- BamFile(bamfile)

    # for each tag call process
    t <- ddply(.data = positions_file, .variables = 'TAG',
               .fun = compute.pileup, .parallel = TRUE)
    # print(str(t

    # write output
    outfile <- sub('[.]bam', '.out', sub('.*/', '', i))
    outfile <- paste(outdir,'/',outfile,sep='')
    write.table(x = t, file = outfile, quote = F, row.names = F, sep = '\t')
}

print('Total time taken...')
time <- proc.time()
print(paste(time[[1]], 'secs', sep = ' '))
timestamp()
