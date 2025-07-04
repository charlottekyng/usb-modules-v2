#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

#---------------
# initialization
#---------------

# load base libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("crayon"))
suppressPackageStartupMessages(library("foreach"))
# suppressPackageStartupMessages(library("Cairo"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("facets"))

#--------------
# parse options
#--------------

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--summaryType", default = c("GL_ASCNA", "GL_LRR", "cnlr.median", "tcn.em", "lcn.em", "cf.em")),
	make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need facets output files\n")
	print_help(parser);
	stop();
} else if (is.null(opt$outFile)) {
	cat("Need output prefix\n")
	print_help(parser);
	stop();
} else if (is.null(opt$genesFile)) {
        cat("Need genes files\n")
        print_help(parser);
        stop();
} else if (!all(opt$summaryType %in% c("GL_ASCNA", "GL_LRR", "tcn.em", "lcn.em", "cnlr.median", "cf.em"))) {
	cat("summaryType can only be one or more of GL_ASCNA, GL_LRR, tcn.em, lcn.em, cnlr.median, cf.em\n")
	print_help(parser);
	stop();
} else {
	facetsFiles <- arguments$args
}

genes <- read.delim(opt$genesFile, as.is=T, check.names=F)
genes$chrom <- gsub("chr", "", genes$chrom)

chroms <- unique(genes$chrom)
chrom_levels <- c(1:22, "X", "Y", "M", "MT")
chrom_levels <- c(chrom_levels, setdiff(chroms, chrom_levels))

genes$chrom <- factor(genes$chrom, levels=chrom_levels)

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)

cat ("Start reading Rdata files\n")
mm <- lapply(facetsFiles, function(f) {

	# Load only the required objects into a local environment
	cat("Reading file ", f, "\n")
	e <- new.env()
	load(f, envir = e)

	tab <- e$fit$cncf
	tab$chrom[which(tab$chrom==23)] <- "X"
	tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
	mcols(tabGR) <- tab %>% select(num.mark,cnlr.median:mafR.clust,cf.em:lcn.em)

	fo <- findOverlaps(tabGR, genesGR)

	df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))
	if ("GL_ASCNA" %in% opt$summaryType) {

		ploidy <- median(unlist(apply(cbind(df$tcn.em, df$num.mark),1,function(x){rep(x[1], x[2])})))

		df$GL_ASCNA <- 0
		df$GL_ASCNA[df$tcn.em < ploidy] <- -1
		df$GL_ASCNA[df$tcn.em == 0] <- -2
		df$GL_ASCNA[df$tcn.em > ploidy] <- 1
		df$GL_ASCNA[df$tcn.em >= ploidy + 4] <- 2
	}
	if ("GL_LRR" %in% opt$summaryType) {

        lrr <- sort(e$out$jointseg$cnlr)
        lrr <- lrr[round(0.25 * length(lrr)):round(0.75 * length(lrr))]

		df$GL_LRR <- 0
		df$GL_LRR[df$cnlr.median < median(lrr)-(2.5*sd(lrr))] <- -1
		df$GL_LRR[df$cnlr.median < median(lrr)-(7*sd(lrr))] <- -2
		df$GL_LRR[df$cnlr.median > median(lrr)+(2*sd(lrr))] <- 1
		df$GL_LRR[df$cnlr.median > median(lrr)+(6*sd(lrr))] <- 2
	}
	
	df %<>% select(hgnc, GL_ASCNA, GL_LRR, tcn.em, lcn.em, cnlr.median, cf.em) %>% ungroup()
	return(df)
	rm(lrr, tab, ploidy, fo, tabGR, e, df)
    gc()

	cat("Finished reading file ", f, "\n")


})
names(mm) <- facetsFiles
for (f in facetsFiles) {
	n <- sub('\\..*', '', sub('.*/', '', f))
	colnames(mm[[f]])[2:ncol(mm[[f]])] <- paste(n,  colnames(mm[[f]])[2:ncol(mm[[f]])], sep="_")
}
cat ("Finished reading RData files\n")

save.image(paste(opt$outFile, ".RData", sep=""))

cat ("Start bind_cols...\n")
mm <- lapply(mm, function(x){
	x[match(genes$hgnc, x$hgnc),-1]
})
mm <- cbind(genes, bind_cols(mm))

cat ("Finished bind_cols...\n")
save.image(paste(opt$outFile, ".RData", sep=""))

seg_sample <- seg_chr <- seg_band <- seg_start <- seg_end <- seg_cnlr <- seg_genes <- seg_type <- seg_GLtype <- seg_cf.em <- NA

for (i in grep("_GL_ASCNA$|_GL_LRR$", colnames(mm))) {
	for(chr in intersect(c(1:22,"X"), unique(mm$chrom))) {
		tt <- mm[which(mm$chrom==chr),c(1:5,i,grep(sub("(.+)_GL_ASCNA$|(.+)_GL_LRR$", "\\1\\2_cf.em", colnames(mm)[i]), colnames(mm))), drop=F]
		tt[which(is.na(tt[,6])),6] <- -1000
		rr <- rle(tt[,6]); 
		if (rr$values[1]== -1000 & length(rr$values)>1) {
			rr$values[1] <- rr$values[2]
		}
		if (rr$values[length(rr$values)]== -1000 & length(rr$values)>1) {
			rr$values[length(rr$values)] <- rr$values[length(rr$values)-1]
		}
		if (length(rr$values)>=3) {
			for ( idx in which(rr$values== -1000)) {
				if (rr$values[idx-1]== rr$values[idx+1]) { rr$values[idx] <- rr$values[idx-1]}
				else {rr$values[idx] <- 0}
			}
		}
		mm[which(mm$chrom==chr),i] <- as.vector(unlist(apply(cbind(rr$value,rr$length), 1, function(x){rep(x[1],x[2])})))

		tt <- mm[which(mm$chrom==chr),c(1:5,i,grep(sub("(.+)_GL_ASCNA$|(.+)_GL_LRR$", "\\1\\2_cf.em", colnames(mm)[i]), colnames(mm))), drop=F]
		rr <- rle(tt[,6]); 
		if (length(rr$length)>1) {
			cs <- cumsum(rr$lengths)
			start <- c(1,cs[1:(length(cs)-1)]+1)
			end <- cs
		} else {start <- 1; end <- rr$lengths[1] }

		for (idx in which(rr$values %in% c(-2,2))) {
			if (rr$values[idx] %in% c(-2,2)) {
				seg_sample <- c(seg_sample, sub("(.+)_GL_ASCNA$|(.+)_GL_LRR$", "\\1\\2", colnames(mm)[i]))
				seg_chr <- c(seg_chr, chr)
				seg_band <- c(seg_band, paste(tt[start[idx],"band"], tt[end[idx],"band"], sep="-"))
				seg_start <- c(seg_start, tt[start[idx],"start"])
				seg_end <- c(seg_end, tt[end[idx],"end"])
#				seg_genes <- c(seg_genes, toString(mm[start[idx]:end[idx],"hgnc"]))
				seg_genes <- c(seg_genes, toString(tt[start[idx]:end[idx],"hgnc"]))
				seg_type <- c(seg_type, rr$values[idx])
				seg_GLtype <- c(seg_GLtype, sub(".+_GL_(ASCNA)$|.+_GL_(LRR)$", "\\1\\2",colnames(mm)[i]))
				seg_cf.em <- c(seg_cf.em, tt[start[idx],grep("cf.em", colnames(tt))])
			}
		}		

	}
}



seg_type[which(seg_type ==  2)] <- "amp"
seg_type[which(seg_type == -2)] <- "del"
mm1 <- data.frame(cbind(seg_sample, seg_chr, seg_band, seg_start, seg_end, seg_genes, seg_type, seg_GLtype, seg_cf.em))
# First raw is always empty. Check just in case and drop it if 'seg_sample' is NA
mm1 <- mm1[-which(is.na(mm1$seg_sample)), ]
write.table(mm1, file=paste(opt$outFile, ".ampdel.txt", sep=""), sep="\t", row.names=F, na="", quote=F)

cat ("Done writing amp/del files\n")
lapply(opt$summaryType, function(c){
	mm2 <- cbind(mm[,1:5], mm[,grep(c, colnames(mm))])
	colnames(mm2) <- gsub(paste("_", c, sep=""), "", colnames(mm2))
	write.table(mm2, file=paste(opt$outFile, ".", c, ".txt", sep=""), sep="\t", row.names=F, na="", quote=F)
})

cat ("All done\n")
