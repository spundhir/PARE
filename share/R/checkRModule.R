#!/usr/bin/env Rscript

# PARESuite: a computational pipeline to Predict Active Regulatory Elements
# Copyright (C) 2015  Sachin Pundhir (pundhir@binf.ku.dk)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="package name"),
    make_option(c("-l", "--list"), help="instead list all installed packages in this output file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input)) & (is.null(opt$list))) {
    cat("\nProgram: checkRModule.R (check if a R module is installed)\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

if((!is.null(opt$input))) {
    is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

    cat(is.installed(opt$input))
} else {
    write.table(installed.packages()[,c(1,2,3,16)], file=opt$list, col.names=T, row.names=F, quote=F)
}
