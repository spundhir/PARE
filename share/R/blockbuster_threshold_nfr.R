#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(gridExtra))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--histoneFile"), help="input file containing tag count at histone enriched regions"),
	make_option(c("-j", "--bkgFile"), help="input file containing tag count at random regions"),
	make_option(c("-k", "--tfFile"), default="NA", help="input file containing tag count at TF bound regions"),
	make_option(c("-m", "--mode"), default=1, help="cut off mode (default: 1=histone percentile; 2=bkg 3rd quantile; 3=intersection; 4=bkg 95th quantile)"),
	make_option(c("-p", "--percentile"), default=0.0005, help="histone percentile (default: %default)"),
	make_option(c("-o", "--outFile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$histoneFile) | is.null(opt$bkgFile) | is.null(opt$outFile)) {
	cat("\nProgram: blockbuster_threshold_nfr.R (R script to determine optimal block threshold)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

#histone <- read.table("histone.reads.countStat")
#bkg <- read.table("bkg.reads.countStat")
#tf <- read.table("tf.reads.countStat")
histone <- read.table(opt$histoneFile)
bkg <- read.table(opt$bkgFile)
histone$class <- "histone"
bkg$class <- "bkg"

if(file.exists(opt$tfFile)) {
    median <- vector()
    histone_percentile <- vector()
    bkg_quartile_3rd <- vector()
    bkg.intersection.point <- vector()
    bkg_percentile_95 <- vector()
    for (col in 1:(ncol(histone)-1)) {
        tf <- read.table(opt$tfFile)
        tf$class <- "tf"

        data <- rbind(histone, bkg, tf)

        # determine intersection point
        lower.limit <- min(log(data[,col]+1))
        upper.limit <- max(log(data[,col]+1))

        data.subset <- subset(data, class=="histone")
        histone.density <- density(log(data.subset[,col]), from=lower.limit, to=upper.limit, n=2^10)
        data.subset <- subset(data, class=="bkg")
        bkg.density <- density(log(data.subset[,col]), from=lower.limit, to=upper.limit, n=2^10)
        data.subset <- subset(data, class=="tf")
        tf.density <- density(log(data.subset[,col]), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y
        density.difference.tf <- histone.density$y - tf.density$y

        bkg.intersection.point.vec <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point <- bkg.intersection.point.vec[1]
        tf.intersection.point.vec <- round(histone.density$x[which(diff(density.difference.tf > 0) != 0) + 1], digits=2)
        tf.intersection.point <- tf.intersection.point.vec[1]

        tf_percentile_95 <- round(log(quantile(tf[,col], probs=0.95, type=8)), digits=2)

        tf_quartile_3rd <- round(log(as.list(summary(tf[,col])[5])[[1]]), digits=2)
        bkg_quartile_3rd <- round(log(as.list(summary(bkg[,col])[5])[[1]]), digits=2)

        median[col] <- sprintf("bkg: %0.2f, tf: %0.2f, histone: %0.2f", median(log(bkg[,col])), median(log(tf[,col])), median(log(histone[,col])))
    }

    pdf(opt$outFile)
    for (col in 1:(ncol(histone)-1)) {
        data.subset <- data[,c(col,ncol(data))]
        colnames(data.subset) <- c("value", "class")
        p1 <- ggplot(data.subset, aes(x=log(value+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("peak median; %s", median)) + xlim(c(0,8)) +
        geom_vline(xintercept=tf_quartile_3rd, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_quartile_3rd+0.4, label=tf_quartile_3rd), y=0.10) +
        geom_vline(xintercept=bkg_quartile_3rd, colour="red", linetype = "longdash") +
        geom_text(aes(x=bkg_quartile_3rd+0.4, label=bkg_quartile_3rd), y=0.20) +
        geom_vline(xintercept=tf.intersection.point, color="blue") +
        geom_text(aes(x=tf.intersection.point+0.4, label=tf.intersection.point), y=0.30) +
        geom_vline(xintercept=bkg.intersection.point, color="red") +
        geom_text(aes(x=bkg.intersection.point+0.4, label=bkg.intersection.point), y=0.40) +
        geom_vline(xintercept=tf_percentile_95, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_percentile_95+0.4, label=tf_percentile_95), y=0.0)

        grid.arrange(p1, top=sprintf("block threshold: %s", threshold_bkg))
    }
    dev.off()

    # set threshold (0: based on 3rd quartile value, 1: based on intersect point value)
    if(as.numeric(opt$mode)==0) {
        threshold_tf <- (round(exp(mean(c(tf_quartile_3rd)))))
        threshold_bkg <- (round(exp(mean(c(bkg_quartile_3rd)))))
    } else if(as.numeric(opt$mode)==1) {
        threshold_tf <- (round(exp(mean(c(tf_percentile_95)))))
        threshold_bkg <- 0
    } else {
        threshold_tf <- (round(exp(mean(c(tf.intersection.point)))))
        threshold_bkg <- (round(exp(mean(c(bkg.intersection.point)))))
    }

    if(threshold_tf > threshold_bkg) {
        cat(threshold_tf)
    } else {
        cat(threshold_bkg)
    }
} else {
    median <- vector()
    histone_percentile <- vector()
    bkg_quartile_3rd <- vector()
    bkg.intersection.point <- vector()
    bkg_percentile_95 <- vector()
    for (col in 1:(ncol(histone)-1)) {
        data <- rbind(histone, bkg)

        # determine intersection point
        lower.limit <- min(log(data[,col]+1))
        upper.limit <- max(log(data[,col]+1))

        # histone percentile
        histone_percentile[col] <- round(log(quantile(histone[,col], probs=as.numeric(as.numeric(opt$percentile)), type=8)), digits=2)
        histone_percentile[col] <- max(histone_percentile[col], 0)

        # bkg quartile
        bkg_quartile_3rd[col] <- round(log(as.list(summary(bkg[,col])[5])[[1]]), digits=2)
        bkg_quartile_3rd[col] <- max(bkg_quartile_3rd[col], 0)

        # intersection point
        data.subset <- subset(data, class=="histone") 
        histone.density <- density(log(data.subset[,col]), from=lower.limit, to=upper.limit, n=2^10)
        data.subset <- subset(data, class=="bkg") 
        bkg.density <- density(log(data.subset[,col]), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y

        bkg.intersection.point.vec <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point[col] <- bkg.intersection.point.vec[1]

        # bkg percentile
        bkg_percentile_95[col] <- round(log(quantile(bkg[,col], probs=0.95, type=8)), digits=2)
        bkg_percentile_95[col] <- max(bkg_percentile_95[col], 0)

        median[col] <- sprintf("bkg median: %0.2f; histone median: %0.2f", median((bkg[,col])), median((histone[,col])))
    }

    # set threshold (0: based on 3rd quartile value, 1: based on intersect point value)
    if(as.numeric(opt$mode)==1) {
        threshold_bkg <- (round(exp(mean(c(histone_percentile)))))
    } else if(as.numeric(opt$mode)==2) {
        threshold_bkg <- (round(exp(mean(c(bkg_quartile_3rd)))))
    } else if(as.numeric(opt$mode)==3) {
        threshold_bkg <- (round(exp(mean(c(bkg.intersection.point)))))
    } else {
        threshold_bkg <- (round(exp(mean(c(bkg_percentile_95)))))
    }

    if(threshold_bkg < 5) {
        tmp <- (round(exp(mean(c(bkg_quartile_3rd)))))
        if(!is.na(tmp)) {
            threshold_bkg <- tmp
        }
    }

    if(threshold_bkg < 5) {
        tmp <- (round(exp(mean(c(bkg.intersection.point)))))
        if(!is.na(tmp)) {
            threshold_bkg <- tmp
        }
    }

    if(threshold_bkg < 5) {
        tmp <- (round(exp(mean(c(bkg_percentile_95)))))
        if(!is.na(tmp)) {
            threshold_bkg <- tmp
        }
    }

    if(threshold_bkg < 5) {
        threshold_bkg <- 5
    }

    cat(threshold_bkg)
    cat("\n")
    #print(threshold_bkg)

    pdf(opt$outFile)
    for (col in 1:(ncol(histone)-1)) {
        data.subset <- data[,c(col,ncol(data))]
        colnames(data.subset) <- c("value", "class")
        p1 <- ggplot(data.subset, aes(x=log(value+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("threshold: %0.0f; %s", threshold_bkg, median[col])) + xlim(c(0,8)) +
        #geom_vline(xintercept=bkg_quartile_3rd, colour="green", linetype = "longdash") +
        #geom_text(aes(x=bkg_quartile_3rd+0.4, label=bkg_quartile_3rd), y=0.10) +
        #geom_vline(xintercept=bkg.intersection.point, color="#9970ab", linetype = "longdash") +
        #geom_text(aes(x=bkg.intersection.point+0.4, label=bkg.intersection.point), y=0.20) +
        #geom_vline(xintercept=bkg_percentile_95, colour="blue", linetype = "longdash") +
        #geom_text(aes(x=bkg_percentile_95+0.4, label=bkg_percentile_95), y=0.0) + 
        #geom_vline(xintercept=histone_percentile, colour="red", linetype = "longdash") +
        #geom_text(aes(x=histone_percentile+0.4, label=histone_percentile), y=0.0)
        geom_vline(xintercept=log(threshold_bkg), colour="red", linetype = "longdash")
        #geom_text(aes(x=log(threshold_bkg)+0.4, label=threshold_bkg), y=0.0)

        grid.arrange(p1, top=sprintf("block threshold: %s", threshold_bkg))
    }
    dev.off()
}
