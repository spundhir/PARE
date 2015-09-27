suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(gridExtra))

# load arguments
ARGV <- commandArgs(TRUE)

# for old blockbuster_threshold_nfr output files
if(is.null(ARGV[4])) {
    cebpa <- read.table(ARGV[1])
    bkg <- read.table(ARGV[2])
    cebpa$class <- "cebpa"
    bkg$class <- "bkg"
    data <- rbind(cebpa, bkg)
    p1 <- ggplot(data, aes(x=log(V1+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 1 (upstream)") + xlim(c(0,8))
    p2 <- ggplot(data, aes(x=log(V2+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 2 (upstream)") + xlim(c(0,8))
    p3 <- ggplot(data, aes(x=log(V3+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 1 (peak)") + xlim(c(0,8))
    p4 <- ggplot(data, aes(x=log(V4+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 2 (peak)") + xlim(c(0,8))
    p5 <- ggplot(data, aes(x=log(V5+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 1 (downstream)") + xlim(c(0,8))
    p6 <- ggplot(data, aes(x=log(V6+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") + ggtitle("replicate 2 (downstream)") + xlim(c(0,8))
    pdf(ARGV[3])
    grid.arrange(p1, p2, p3, p4, p5, p6)
    dev.off()
} else {
    #histone <- read.table("histone.reads.countStat")
    #bkg <- read.table("bkg.reads.countStat")
    #tf <- read.table("tf.reads.countStat")
    histone <- read.table(ARGV[1])
    bkg <- read.table(ARGV[2])
    histone$class <- "histone"
    bkg$class <- "bkg"

    if(file.exists(ARGV[3])) {
        tf <- read.table(ARGV[3])
        tf$class <- "tf"
    }

    if(file.exists(ARGV[3])) {
        data <- rbind(histone, bkg, tf)

        pdf(ARGV[4])
        ## replicate 1
        # determine intersection point
        lower.limit <- min(log(data$V1+1))
        upper.limit <- max(log(data$V1+1))

        histone.density <- density(log(subset(data, class=="histone")$V1), from=lower.limit, to=upper.limit, n=2^10)
        bkg.density <- density(log(subset(data, class=="bkg")$V1), from=lower.limit, to=upper.limit, n=2^10)
        tf.density <- density(log(subset(data, class=="tf")$V1), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y
        density.difference.tf <- histone.density$y - tf.density$y

        bkg.intersection.point.rep1 <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point.rep1 <- bkg.intersection.point.rep1[1]
        tf.intersection.point.rep1 <- round(histone.density$x[which(diff(density.difference.tf > 0) != 0) + 1], digits=2)
        tf.intersection.point.rep1 <- tf.intersection.point.rep1[1]

        tf_percentile_95_rep1 <- round(log(quantile(tf$V1, probs=0.95, type=8)), digits=2)

        tf_quartile_3rd_rep1 <- round(log(as.list(summary(tf$V1)[5])[[1]]), digits=2)
        bkg_quartile_3rd_rep1 <- round(log(as.list(summary(bkg$V1)[5])[[1]]), digits=2)

        median_rep1 <- sprintf("bkg: %0.2f, tf: %0.2f, histone: %0.2f", median(log(bkg$V1)), median(log(tf$V1)), median(log(histone$V1)))

        p1 <- ggplot(data, aes(x=log(V1+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("replicate 1 (peak median; %s)", median_rep1)) + xlim(c(0,8)) +
        geom_vline(xintercept=tf_quartile_3rd_rep1, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_quartile_3rd_rep1+0.4, label=tf_quartile_3rd_rep1), y=0.10) +
        geom_vline(xintercept=bkg_quartile_3rd_rep1, colour="red", linetype = "longdash") +
        geom_text(aes(x=bkg_quartile_3rd_rep1+0.4, label=bkg_quartile_3rd_rep1), y=0.20) +
        geom_vline(xintercept=tf.intersection.point.rep1, color="blue") +
        geom_text(aes(x=tf.intersection.point.rep1+0.4, label=tf.intersection.point.rep1), y=0.30) +
        geom_vline(xintercept=bkg.intersection.point.rep1, color="red") +
        geom_text(aes(x=bkg.intersection.point.rep1+0.4, label=bkg.intersection.point.rep1), y=0.40) +
        geom_vline(xintercept=tf_percentile_95_rep1, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_percentile_95_rep1+0.4, label=tf_percentile_95_rep1), y=0.0)

        ## replicate 2
        # determine intersection point
        lower.limit <- min(log(data$V2+1))
        upper.limit <- max(log(data$V2+1))

        histone.density <- density(log(subset(data, class=="histone")$V2), from=lower.limit, to=upper.limit, n=2^10)
        bkg.density <- density(log(subset(data, class=="bkg")$V2), from=lower.limit, to=upper.limit, n=2^10)
        tf.density <- density(log(subset(data, class=="tf")$V2), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y
        density.difference.tf <- histone.density$y - tf.density$y

        bkg.intersection.point.rep2 <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point.rep2 <- bkg.intersection.point.rep2[1]
        tf.intersection.point.rep2 <- round(histone.density$x[which(diff(density.difference.tf > 0) != 0) + 1], digits=2)
        tf.intersection.point.rep2 <- tf.intersection.point.rep2[1]

        tf_percentile_95_rep2 <- round(log(quantile(tf$V2, probs=0.95, type=8)), digits=2)
        
        tf_quartile_3rd_rep2 <- round(log(as.list(summary(tf$V2)[5])[[1]]), digits=2)
        bkg_quartile_3rd_rep2 <- round(log(as.list(summary(bkg$V2)[5])[[1]]), digits=2)

        median_rep2 <- sprintf("bkg: %0.2f, tf: %0.2f, histone: %0.2f", median(log(bkg$V2)), median(log(tf$V2)), median(log(histone$V2)))

        p2 <- ggplot(data, aes(x=log(V2+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("replicate 2 (peak median; %s)", median_rep2)) + xlim(c(0,8)) +
        geom_vline(xintercept=tf_quartile_3rd_rep2, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_quartile_3rd_rep2+0.4, label=tf_quartile_3rd_rep2), y=0.10) +
        geom_vline(xintercept=bkg_quartile_3rd_rep2, colour="red", linetype = "longdash") +
        geom_text(aes(x=bkg_quartile_3rd_rep2+0.4, label=bkg_quartile_3rd_rep2), y=0.20) +
        geom_vline(xintercept=tf.intersection.point.rep2, color="blue") +
        geom_text(aes(x=tf.intersection.point.rep2+0.4, label=tf.intersection.point.rep2), y=0.30) +
        geom_vline(xintercept=bkg.intersection.point.rep2, color="red") +
        geom_text(aes(x=bkg.intersection.point.rep2+0.4, label=bkg.intersection.point.rep2), y=0.40) +
        geom_vline(xintercept=tf_percentile_95_rep2, colour="blue", linetype = "longdash") +
        geom_text(aes(x=tf_percentile_95_rep2+0.4, label=tf_percentile_95_rep2), y=0.0)

        grid.arrange(p1, p2)
        dev.off()

        # set threshold (0: based on 3rd quartile value, 1: based on intersect point value)
        if(ARGV[5]==0) {
            threshold_tf <- (round(exp(mean(c(tf_quartile_3rd_rep1, tf_quartile_3rd_rep2)))))
            threshold_bkg <- (round(exp(mean(c(bkg_quartile_3rd_rep1, bkg_quartile_3rd_rep2)))))
        } else if(ARGV[5]==1) {
            threshold_tf <- (round(exp(mean(c(tf_percentile_95_rep1, tf_percentile_95_rep2)))))
            threshold_bkg <- 0
        } else {
            threshold_tf <- (round(exp(mean(c(tf.intersection.point.rep1, tf.intersection.point.rep2)))))
            threshold_bkg <- (round(exp(mean(c(bkg.intersection.point.rep1, bkg.intersection.point.rep2)))))
        }

        if(threshold_tf > threshold_bkg) {
            cat(threshold_tf)
        } else {
            cat(threshold_bkg)
        }
    } else {
        data <- rbind(histone, bkg)

        pdf(ARGV[4])
        ## replicate 1
        # determine intersection point
        lower.limit <- min(log(data$V1+1))
        upper.limit <- max(log(data$V1+1))

        histone.density <- density(log(subset(data, class=="histone")$V1), from=lower.limit, to=upper.limit, n=2^10)
        bkg.density <- density(log(subset(data, class=="bkg")$V1), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y

        bkg.intersection.point.rep1 <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point.rep1 <- bkg.intersection.point.rep1[1]

        #histone_percentile_rep1 <- round(log(min(histone$V1)), digits=2)
        histone_percentile_rep1 <- round(log(quantile(histone$V1, probs=as.numeric(ARGV[6]), type=8)), digits=2)
        histone_percentile_rep1 <- max(histone_percentile_rep1, 0)
        bkg_percentile_95_rep1 <- round(log(quantile(bkg$V1, probs=0.90, type=8)), digits=2)
        bkg_percentile_95_rep1 <- max(bkg_percentile_95_rep1, 0)
        bkg_quartile_3rd_rep1 <- round(log(as.list(summary(bkg$V1)[5])[[1]]), digits=2)
        bkg_quartile_3rd_rep1 <- max(bkg_quartile_3rd_rep1, 0)

        median_rep1 <- sprintf("bkg: %0.2f, histone: %0.2f", median(log(bkg$V1)), median(log(histone$V1)))

        p1 <- ggplot(data, aes(x=log(V1+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("replicate 1 (peak median; %s)", median_rep1)) + xlim(c(0,8)) +
        geom_vline(xintercept=bkg_quartile_3rd_rep1, colour="green", linetype = "longdash") +
        geom_text(aes(x=bkg_quartile_3rd_rep1+0.4, label=bkg_quartile_3rd_rep1), y=0.10) +
        geom_vline(xintercept=bkg.intersection.point.rep1, color="#9970ab", linetype = "longdash") +
        geom_text(aes(x=bkg.intersection.point.rep1+0.4, label=bkg.intersection.point.rep1), y=0.20) +
        geom_vline(xintercept=bkg_percentile_95_rep1, colour="red", linetype = "longdash") +
        geom_text(aes(x=bkg_percentile_95_rep1+0.4, label=bkg_percentile_95_rep1), y=0.0) + 
        geom_vline(xintercept=histone_percentile_rep1, colour="blue", linetype = "longdash") +
        geom_text(aes(x=histone_percentile_rep1+0.4, label=histone_percentile_rep1), y=0.0)

        ## replicate 2
        # determine intersection point
        lower.limit <- min(log(data$V2+1))
        upper.limit <- max(log(data$V2+1))

        histone.density <- density(log(subset(data, class=="histone")$V2), from=lower.limit, to=upper.limit, n=2^10)
        bkg.density <- density(log(subset(data, class=="bkg")$V2), from=lower.limit, to=upper.limit, n=2^10)

        density.difference.bkg <- histone.density$y - bkg.density$y

        bkg.intersection.point.rep2 <- round(histone.density$x[which(diff(density.difference.bkg > 0) != 0) + 1], digits=2)
        bkg.intersection.point.rep2 <- bkg.intersection.point.rep2[1]

        #histone_percentile_rep2 <- round(log(min(histone$V2)), digits=2)
        histone_percentile_rep2 <- round(log(quantile(histone$V2, probs=as.numeric(ARGV[6]), type=8)), digits=2)
        histone_percentile_rep2 <- max(histone_percentile_rep2, 0)
        bkg_percentile_95_rep2 <- round(log(quantile(bkg$V2, probs=0.90, type=8)), digits=2)
        bkg_percentile_95_rep2 <- max(bkg_percentile_95_rep2, 0)

        bkg_quartile_3rd_rep2 <- round(log(as.list(summary(bkg$V2)[5])[[1]]), digits=2)
        bkg_quartile_3rd_rep2 <- max(bkg_quartile_3rd_rep2, 0)

        median_rep2 <- sprintf("bkg: %0.2f, histone: %0.2f", median(log(bkg$V2)), median(log(histone$V2)))

        p2 <- ggplot(data, aes(x=log(V2+1), fill=class)) + geom_density(alpha=.3) + xlab("read count (log)") +
        ggtitle(sprintf("replicate 2 (peak median; %s)", median_rep2)) + xlim(c(0,8)) +
        geom_vline(xintercept=bkg_quartile_3rd_rep2, colour="green", linetype = "longdash") +
        geom_text(aes(x=bkg_quartile_3rd_rep2+0.4, label=bkg_quartile_3rd_rep2), y=0.10) +
        geom_vline(xintercept=bkg.intersection.point.rep2, color="#9970ab", linetype = "longdash") +
        geom_text(aes(x=bkg.intersection.point.rep2+0.4, label=bkg.intersection.point.rep2), y=0.20) +
        geom_vline(xintercept=bkg_percentile_95_rep2, colour="red", linetype = "longdash") +
        geom_text(aes(x=bkg_percentile_95_rep2+0.4, label=bkg_percentile_95_rep2), y=0.0) +
        geom_vline(xintercept=histone_percentile_rep2, colour="blue", linetype = "longdash") +
        geom_text(aes(x=histone_percentile_rep2+0.4, label=histone_percentile_rep2), y=0.0)

        grid.arrange(p1, p2)
        dev.off()

        # set threshold (0: based on 3rd quartile value, 1: based on intersect point value)
        if(ARGV[5]==1) {
            #threshold_bkg <- (round(exp(mean(c(bkg_percentile_95_rep1, bkg_percentile_95_rep2)))))
            #threshold_bkg <- max((round(exp(mean(c(histone_percentile_rep1, histone_percentile_rep2))))),10)
            threshold_bkg <- (round(exp(mean(c(histone_percentile_rep1, histone_percentile_rep2)))))
        } else if(ARGV[5]==2) {
            threshold_bkg <- (round(exp(mean(c(bkg_quartile_3rd_rep1, bkg_quartile_3rd_rep2)))))
        } else if(ARGV[5]==3) {
            threshold_bkg <- (round(exp(mean(c(bkg.intersection.point.rep1, bkg.intersection.point.rep2)))))
        } else {
            threshold_bkg <- (round(exp(mean(c(bkg_percentile_95_rep1, bkg_percentile_95_rep2)))))
        }

        if(threshold_bkg < 5) {
            tmp <- (round(exp(mean(c(bkg_quartile_3rd_rep1, bkg_quartile_3rd_rep2)))))
            if(!is.na(tmp)) {
                threshold_bkg <- tmp
            }
        }

        if(threshold_bkg < 5) {
            tmp <- (round(exp(mean(c(bkg.intersection.point.rep1, bkg.intersection.point.rep2)))))
            if(!is.na(tmp)) {
                threshold_bkg <- tmp
            }
        }

        if(threshold_bkg < 5) {
            tmp <- (round(exp(mean(c(bkg_percentile_95_rep1, bkg_percentile_95_rep2)))))
            if(!is.na(tmp)) {
                threshold_bkg <- tmp
            }
        }

        if(threshold_bkg < 5) {
            threshold_bkg <- 5
        }

        cat(threshold_bkg)
    }
}
