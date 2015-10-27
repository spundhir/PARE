#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--bedFile"), help="input BED file or stdin"),
    make_option(c("-l", "--list"), help="input is a list", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$bedFile)) {
	cat("\nProgram: bedStat.R (R script to compute statistics on BED file)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

#if(is.na(file.info("stdin")$size)) { print("Hello"); }
if(identical(opt$bedFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$bedFile)
}

if(!is.null(opt$list)){
    cat(sprintf("Median=%0.0f\n", median(data$V1)))
    cat(sprintf("Mean=%0.0f\n", mean(data$V1)))
    cat(sprintf("Sum=%0.0f\n", sum(data$V1)))
    cat(sprintf("Max=%0.0f\n", max(data$V1)))
    cat(sprintf("Min=%0.0f\n", min(data$V1)))
    cat(sprintf("Count=%0.0f\n", length(data$V1)))
}else{
    cat(sprintf("Median=%0.0f\n", median(data$V3-data$V2)))
    cat(sprintf("Mean=%0.0f\n", mean(data$V3-data$V2)))
    cat(sprintf("Max=%0.0f\n", max(data$V3-data$V2)))
    cat(sprintf("Min=%0.0f\n", min(data$V3-data$V2)))
    cat(sprintf("Count=%0.0f\n", length(data$V3-data$V2)))
}
