#!/usr/bin/env Rscript
# segment copy number data and generate plot

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("CGHcall"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--centromereFile", default = NULL, type = "character", action = "store", help ="centromere file"),
                make_option("--alpha", default = 0.000001, type = "double", action = "store", help ="alpha"),
                make_option("--smoothRegion", default = 10, type = "double", action = "store", help ="smooth region"),
                make_option("--outlierSDscale", default = 2.5, type = "double", action = "store", help ="outlier SD scale"),
                make_option("--undoSD", default = 2, type = "double", action = "store", help ="undo SD"),
                make_option("--prefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"),
		make_option("--perchromplots", default = FALSE, type="logical", action = "store", help = "should be make per chrom plots"))

parser <- OptionParser(usage = "%prog [options] inDir", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need copy number file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$prefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
} else {
    cnFile <- arguments$args;
}

chroms <- c(1:22, "X")

cn <- lapply(cnFile, read.table, header=T, as.is=T)
cn <- lapply(cn, function(x) {x$log2 <- scale(x$log2, scale=F); x})
cn <- do.call("rbind", cn)
#cn <- read.table(cnFile, header=T, as.is=T)
cn[,1] <- gsub("chr", "", cn[,1])
cn <- cn[order(factor(cn[,1], levels=chroms), as.numeric(cn[,2])),,drop=F]
keep <- which(cn[,1] %in% chroms)
chroms <- unique(cn[keep,1])
if (length(rm) > 0) { cn <- cn[keep,]}
cn[which(cn[,1]=="X"),1] <- 23
cn[,1] <- as.numeric(cn[,1])
cn <- cn[order(cn[,1], cn[,2], cn[,3]),]
cn <- cbind(name = paste(cn[,1], cn[,2], cn[,3], sep="_"), cn[,c(1:3,6)])
cn <- cn[which(!duplicated(cn$name)),]
cgh <- make_cghRaw(cn)
normalized <- normalize(cgh, smoothOutliers=T, trim=0.025, smooth.region=opt$smoothRegion, outlier.SD.scale=opt$outlierSDscale)
segmented <- segmentData(normalized, relSDlong=3, undo.splits="sdundo", undo.SD=opt$undoSD, alpha=opt$alpha, trim=0.025)

fn <- paste(opt$prefix, '.segment.Rdata', sep = '')
save(segmented, file = fn)
Data <- cbind(fData(segmented), copynumber(segmented), segmented(segmented))
colnames(Data)[5] <- "log2_ratio_seg"
write.table(Data, file = paste(opt$prefix, ".seg.txt", sep=""), col.names=NA, quote=F, sep="\t")

f <- factor(paste(Data$Chromosome, Data$log2_ratio_seg))
sData <- split(Data, f)
collapsedData <- do.call('rbind', lapply(sData, function(x) {
                            c(Chromosome = x[1,"Chromosome"], Start = x[1, "Start"], End = x[nrow(x), "End"], nBins = nrow(x), log2Ratio = x[1, "log2_ratio_seg"]) }))
oo <- order(collapsedData[, "Chromosome"], collapsedData[, "Start"])
collapsedData <- collapsedData[oo, ]
write.table(collapsedData, file = paste(opt$prefix, ".collapsed_seg.txt", sep = ""), row.names = F, quote = F, sep = "\t")

ylim <- c(min(as.numeric(Data[,4])), max(as.numeric(Data[,4])))
ylim[2] <- ylim[2]+0.5

rlechr <- rle(Data$Chr)

pdf(paste(opt$prefix,".seg_plot.pdf", sep=""), height=5, width=18)
plot(as.numeric(Data[,4]), pch=20, xlab='Position', ylab="Copy number", xaxt='n', ylim=ylim)
points(as.numeric(Data[,5]), pch = 20, col = 'blue')
abline(v=cumsum(rlechr$lengths), col="red", lty=3)
abline(h=0, col="darkgrey", lty=3)
text(cumsum(rlechr$lengths)-(rlechr$lengths/2), ylim[2]-0.25, rlechr$values)

if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1]==j)[1],3]
        index <- which(Data$Chromosome==j & Data$Start > pos)[1]
        if (!is.na(index)) {
            abline(v=index, col="darkgrey", lty=3)
        }
    }
}
dev.off()

png(paste(opt$prefix,".seg_plot.png", sep=""), type = 'cairo-png', height=400, width=2000)
plot(as.numeric(Data[,4]), pch=20, xlab='Position', ylab="Copy number", xaxt='n', ylim=ylim)
points(as.numeric(Data[,5]), pch = 20, col = 'blue')
abline(v=cumsum(rle(Data$Chr)$lengths), col="red", lty=3)
abline(h=0, col="darkgrey", lty=3)
text(cumsum(rlechr$lengths)-(rlechr$lengths/2),	ylim[2]-0.25, rlechr$values)

if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1]==j)[1],3]
        index <- which(Data$Chromosome==j & Data$Start > pos)[1]
        if (!is.na(index)) {
            abline(v=index, col="darkgrey", lty=3)
        }
    }
}
dev.off()

if(opt$perchromplots) {
if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
}

for (chr in chroms) {
    chrData <- subset(Data, Chromosome == chr)
    if (nrow(chrData) > 0) {
        ylim <- c(min(as.numeric(chrData[,4])), max(as.numeric(chrData[,4])))
        ylim[2] <- ylim[2]+0.5
        png(paste(opt$prefix,".seg_plot.", chr, ".png", sep=""), type = 'cairo-png', height=500, width=500)
        plot(as.numeric(chrData[,4]), pch=20, xlab='Position', ylab="Copy number", ylim=ylim, main = paste('Chromosome', chr))
        points(as.numeric(chrData[,5]), pch = 20, col = 'blue')

        if (!is.null(opt$centromereFile)) {
            for (j in unique(cen[,1])) {
                pos <- cen[which(cen[,1]==j)[1],3]
                index <- which(chrData$Chromosome==j & chrData$Start > pos)[1]
                if (!is.na(index)) {
                    abline(v=index, col="darkgrey", lty=3)
                }
                text(cumsum(rle(chrData$Chromosome)$lengths)-((rle(chrData$Chromosome)$lengths)/2), ylim[2]-0.25)
            }
        }
        dev.off()
    }
}
}

