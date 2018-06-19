suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("parallel"));
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--chr_col", default = "CHROM", type= 'character', help = "column name for chr"),
	make_option("--pos_col", default = "POS", type='character', help = "column name for pos"),
	make_option("--ref_col", default = "REF", type='character', help = "column name for ref"),
	make_option("--alt_col", default = "ALT", type='character', help = "column name for alt"),
	make_option("--tri.count.method", default = "exome2genome", help = "tri.count.method for deconstructSigs"),
	make_option("--num_iter", default = NA, type='integer', help = "number of re-sampling with replacement"),
	make_option("--num_cores", default = 1, type='integer', help = "number of cores to use"),
	make_option("--seed", default = 1237, type='integer', help = "seed for randomization"),
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

#plotdir="plots"
#if(!file.exists(plotdir)){dir.create(plotdir)}


allmuts <- read.delim(muts_file, as.is=T)

allmuts <- allmuts[grep("interrogation", allmuts$FILTER,invert=T),]

allmuts$chr <- paste("chr", allmuts$CHROM, sep="")

if (is.na(opt$num_iter)){
	muts_for_input <- allmuts
} else {
	muts_for_input <- do.call("rbind", lapply(unique(allmuts$TUMOR_SAMPLE), function(s){
		x <- subset(allmuts, TUMOR_SAMPLE==s)
		set.seed(opt$seed)
		y <- x[sample(1:nrow(x), opt$num_iter*nrow(x), replace=T),]
		y$TUMOR_SAMPLE = paste(unique(y$TUMOR_SAMPLE), rep(1:opt$num_iter, nrow(x)), sep="_")
		y
	}))
}

sigs <- mut.to.sigs.input(muts_for_input, "TUMOR_SAMPLE", "chr","POS", "REF", "ALT")

if (opt$num_cores>1) {	
	cl <- makeCluster(opt$num_cores, "SOCK")
	ws <- parLapply(cl, rownames(sigs), function(sample,sigs) {
		library(deconstructSigs)
		whichSignatures(tumor.ref = sigs, sample.id=sample,
                signatures.ref = signatures.cosmic, contexts.needed = T,
                tri.counts.method = "exome2genome")
	}, sigs)
else {
	ws <- lapply(rownames(sigs), function(sample) {
		whichSignatures(tumor.ref = sigs, sample.id=sample,
		signatures.ref = signatures.cosmic, contexts.needed = T,
		tri.counts.method = "exome2genome")
	})

}
names(ws) <- rownames(sigs)
if (!is.na(opt$num_iter)){
	if(opt$num_iter>1){
	summarise_whichSignatures <- function(x, signatures=signatures.cosmic, sampleName) {
		weights_mat <- do.call("rbind", lapply(x, function(y){y$weights}))
		weights <- matrix(colMeans(weights_mat),nrow=1)
		weightsSD <- apply(weights_mat, 2, sd)
		colnames(weights) <- colnames(weights_mat)
		rownames(weights) <- sampleName
		unknown <- 1 - sum(weights)
		product <- weights %*% as.matrix(signatures)
		tumor_mat <- do.call("rbind", lapply(x, function(y){y$tumor}))
		tumor <- matrix(colMeans(tumor_mat), nrow=1)
		colnames(tumor) <- colnames(tumor_mat)
		rownames(tumor) <- sampleName
		diff <- tumor - product
		out <- list(weights, tumor, product, diff, unknown, weightsSD)
		names(out) <- c("weights", "tumor", "product", "diff", "unknown", "weightsSD")
		return(out)
	}
	ws2 <- lapply(unique(allmuts$TUMOR_SAMPLE), function(s){
		summarise_whichSignatures(ws[grep(paste(s, "_", sep=""), names(ws))], sampleName=s)
	})
	names(ws2) <- unique(allmuts$TUMOR_SAMPLE)
	ws <- ws2
}
}

lapply(names(ws), function(sample) {

	pdf(paste(opt$outPrefix, "_", sample, ".pdf", sep=""))
	plotSignatures(ws[[sample]])
	dev.off()
})

signatures <- do.call("rbind", lapply(ws, function(w) { w$weights}))

if (!is.na(opt$num_iter)){
	signaturesSD <- do.call("rbind", lapply(ws, function(w) { w$weightsSD}))
	colnames(signaturesSD) <- paste(colnames(signaturesSD), "SD", sep="")
	signatures <- cbind(signatures, signaturesSD)
	
}
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

write.table(signatures, file=paste(opt$outPrefix, ".txt", sep="")
save(signatures, sigs, file=paste(opt$outPrefix, ".RData", sep=""))



