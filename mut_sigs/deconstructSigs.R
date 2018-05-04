suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--chr_col", default = "CHROM", type= 'character', help = "column name for chr"),
	make_option("--pos_col", default = "POS", type='character', help = "column name for pos"),
	make_option("--ref_col", default = "REF", type='character', help = "column name for ref"),
	make_option("--alt_col", default = "ALT", type='character', help = "column name for alt"),
	make_option("--tri.count.method", default = "exome2genome", help = "tri.count.method for deconstructSigs"),
	make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutations file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    muts_file <- arguments$args[1];
}

plotdir="plots"
if(!file.exists(plotdir)){dir.create(plotdir)}


allmuts <- read.delim(muts_file, as.is=T)

allmuts <- allmuts[grep("interrogation", allmuts$FILTER,invert=T),]

allmuts$chr <- paste("chr", allmuts$CHROM, sep="")
sigs <- mut.to.sigs.input(allmuts, "TUMOR_SAMPLE", "chr","POS", "REF", "ALT")

ws <- lapply(rownames(sigs), function(sample) {
	whichSignatures(tumor.ref = sigs, sample.id=sample,
		signatures.ref = signatures.cosmic, contexts.needed = T,
		tri.counts.method = "exome2genome")
})
names(ws) <- rownames(sigs)
lapply(names(ws), function(sample) {

	pdf(paste(plotdir, "/", outPrefix, "_", sample, ".pdf", sep=""))
	plotSignatures(ws[[sample]])
	dev.off()
})

signatures <- do.call("rbind", lapply(ws, function(w) { w$weights}))

mutrate <- as.numeric(table(allmuts$TUMOR_SAMPLE))
mutrate_snv <- as.numeric(table(allmuts$TUMOR_SAMPLE[which(allmuts$REF %in% c("A", "C", "G", "T") & allmuts$ALT %in% c("A", "C", "G", "T"))]))
mutrate_indel <- mutrate-mutrate_snv

mutrate <- data.frame(
	TOTAL = mutrate,
	SNV = mutrate_snv,
	INDEL = mutrate_indel,
	SNV_prop =mutrate_snv/mutrate,
	INDEL_prop = mutrate_indel/mutrate)
rownames(mutrate) <- sort(unique(allmuts$TUMOR_SAMPLE))

signatures <- signatures[match(rownames(mutrate), rownames(signatures)),]

save(signatures, sigs, file=paste(outPrefix, ".RData", sep=""))



