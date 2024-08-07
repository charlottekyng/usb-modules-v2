#!/usr/bin/env Rscript
# cluster vcf and output dendrogram

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("gplots"));

options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
	make_option("--genome", default = 'b37', help = "genome build [default %default]"),
	make_option("--outPrefix", default = NULL, help = "output prefix [default %default]"),
	make_option("--clusterType", default = "hetDiffHom", help = "hetDiffHom or hetSameHom, should het be treated the same as hom?"))

parser <- OptionParser(usage = "%prog vcf.files", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

vcfFile <- arguments$args[1]


vcf <- readVcf(vcfFile, opt$genome)

if (ncol(vcf)<3) { stop("Not enough samples to do clustering")}
gt <- geno(vcf)$GT
X <- matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
X[gt == "0/0"] <- 0
X[gt == "0/1"] <- 1
if (opt$clusterType=="hetDiffHom") {
	X[gt == "1/1"] <- 2 
} else { X[gt == "1/1"] <- 1 }
X[!gt %in% c("0/0", "0/1", "1/1")] <- NA
#plot(hclust(dist(t(X), method = 'manhattan')))

gt <- matrix(as.integer(factor(gt)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))

fn <- paste(opt$outPrefix, ".clust.png", sep = '')
png(fn, height = 600, width = 3000)
null <- plot(hclust(dist(t(gt)), method = 'ward.D2'))
dev.off()

fn <- paste(opt$outPrefix, ".heatmap.ward.pdf", sep = '')
pdf(fn, height = max(10, 5+(ncol(gt)/5)), width = max(10, 5+(ncol(gt)/5)))
null <- heatmap.2(as.matrix(dist(t(gt))), hclustfun=function(x){hclust(x, method='ward.D2')}, 
	scale = 'none', trace = 'none', keysize = 0.3, margins = c(10,10))
dev.off()

fn <- paste(opt$outPrefix, ".heatmap.completelink.pdf", sep = '')
pdf(fn, height = max(10, 5+(ncol(gt)/5)), width = max(10, 5+(ncol(gt)/5)))
null <- heatmap.2(as.matrix(dist(t(gt))),
        scale = 'none', trace = 'none', keysize = 0.3, margins = c(10,10))
dev.off()
