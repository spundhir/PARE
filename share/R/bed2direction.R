#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing results from bed2direction script"),
	make_option(c("-j", "--sessionFile"), help="input session file containing SVM model"),
	make_option(c("-o", "--outFile"), help="output file containing predictions")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: bed2direction.R (R script to predict NFR directionality)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(randomForest))

load(opt$sessionFile)
test <- read.table(opt$inFile)
res <- predict(model, test[,c(7:14)])
test[,(ncol(test)+1)] <- res
write.table(test, opt$outFile, sep="\t", row.names=F, col.names=F, quote=F)
