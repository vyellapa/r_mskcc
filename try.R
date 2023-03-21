suppressMessages(library(optparse, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))


option_list = list(make_option(c("-i", "--input"), type = "character",
                               help = "input file",
                               metavar = "path", action = "store"),
                   make_option(c("-o", "--output"), type = "character",
                               help = "output directory", metavar = "path"),
                   make_option(c("-d","--diplogR"), type = "double",
                               help = "alternate diplogR to use",
                               default = 0.1, metavar = "integer"));



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt$d=1.01
cat (opt$i,opt$o, opt$d,"\n")
