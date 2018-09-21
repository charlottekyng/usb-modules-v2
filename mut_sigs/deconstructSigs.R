cat ("Running deconstructSigs.R\n\n")

suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("parallel"));
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--sample_col", default = "TUMOR_SAMPLE", type='character', help = "column name for tumor sample"),
	make_option("--chr_col", default = "CHROM", type= 'character', help = "column name for chr"),
	make_option("--pos_col", default = "POS", type='character', help = "column name for pos"),
	make_option("--ref_col", default = "REF", type='character', help = "column name for ref"),
	make_option("--alt_col", default = "ALT", type='character', help = "column name for alt"),
	make_option("--tri.count.method", default = "exome2genome", help = "tri.count.method for deconstructSigs"),
	make_option("--num_iter", default = NA, type='integer', help = "number of re-sampling with replacement (at least 10, otherwise NA)"),
	make_option("--num_cores", default = 1, type='integer', help = "number of cores to use"),
	make_option("--min_muts_to_include", default = 20, type='integer', help = "minimum number of mutations required to derive signature"),
	make_option("--seed", default = 1237, type='integer', help = "seed for randomization"),
	make_option("--tumorSample", default = NULL, type='character', help = "tumor samples to run"),
	make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [mutation_summary_file]", option_list = optList);

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

cat ("Reading mutation summary file\n")
allmuts <- read.delim(muts_file, as.is=T)

if (!is.null(opt$tumorSample)) {
	samples <- opt$tumorSample[which(opt$tumorSample %in% allmuts[,opt$sample_col])]
	missing_samples <- opt$tumorSample[which(!opt$tumorSample %in% allmuts[,opt$sample_col])]
	if (length(samples)>0) {
		cat ("These samples have been found in the mutations file", samples, "\n")
	}
	if (length(missing_samples)>0) { 
		cat ("At least one sample appears to have no mutations:", missing_samples, "\n")
	}
	allmuts <- allmuts[which(allmuts[,opt$sample_col] %in% opt$tumorSample),,drop=F]
	sample_levels <- opt$tumorSample
} else { cat ("Using all samples in the mutations file\n") 
	sample_levels <- sort(unique(allmuts[,opt$sample_col]))}

cat ("Only using SNVs for deconstructSigs\n")
sn <- c("A", "C", "G", "T")
pointmuts <- allmuts[which(allmuts[,opt$ref_col] %in% sn & allmuts[,opt$alt_col] %in% sn),,drop=F]


if(nrow(pointmuts)>0) {

	if (length(grep("chr", pointmuts[1,opt$chr_col]))==0) {
		cat("chr_col colum appears not to have chr - appending chr\n")
		pointmuts$chr <- paste("chr", pointmuts[,opt$chr_col], sep="")
	} else { pointmuts$chr <- pointmuts[,opt$chr_col] }
	pointmuts$chr <- gsub("chrMT", "chrM", pointmuts$chr)

	if (is.na(opt$num_iter)){
		cat("No bootstrapping to be performed\n")
		muts_for_input <- pointmuts
	} else {
		if (opt$num_iter<10) {
			cat ("Number of iterations too small (at least 10, otherwise what is the point?). Not going to perform bootstrapping.\n")
			muts_for_input <- pointmuts
		} else {
			muts_for_input <- do.call("rbind", lapply(unique(pointmuts[,opt$sample_col]), function(s){
				x <- pointmuts[which(pointmuts[,opt$sample_col] ==s),,drop=F]
				set.seed(opt$seed)
				y <- x[sample(1:nrow(x), opt$num_iter*nrow(x), replace=T),]
				y$TUMOR_SAMPLE = paste(unique(y$TUMOR_SAMPLE), rep(1:opt$num_iter, nrow(x)), sep="_")
				y
			}))
			cat("Finished bootstrapping mutations\n")
		}
	}

	sigs <- mut.to.sigs.input(muts_for_input, opt$sample_col, "chr", opt$pos_col, opt$ref_col, opt$alt_col)

	rs <- rowSums(sigs)

	toofewmuts <- rownames(sigs)[which(rs<opt$min_muts_to_include)]
	if(length(toofewmuts)>0) {
		cat ("These samples had <", opt$min_muts_to_include, "point mutations:", toofewmuts, " Removing these samples\n")
		sigs <- sigs[which(!rownames(sigs) %in% toofewmuts),,drop=F]
	}

	if(nrow(sigs) > 0) {
	
		cat ("Starting whichSignatures:", date(), "\n")
		if (opt$num_cores>1) {
			cl <- makeCluster(opt$num_cores, "SOCK")
			ws <- parLapply(cl, rownames(sigs), function(sample,sigs, ...) {
				library(deconstructSigs)
				whichSignatures(tumor.ref = sigs, sample.id=sample,
      	          signatures.ref = signatures.cosmic, contexts.needed = T, ...)
			}, sigs, tri.counts.method = opt$tri.count.method)
			stopCluster(cl)
		} else {
			ws <- lapply(rownames(sigs), function(sample) {
				whichSignatures(tumor.ref = sigs, sample.id=sample,
				signatures.ref = signatures.cosmic, contexts.needed = T,
				tri.counts.method = opt$tri.count.method)
			})
		}
		names(ws) <- rownames(sigs)
		cat("Finished whichSignatures", date(), "\n")

		if (!is.na(opt$num_iter)){
			if(opt$num_iter>=10){
				cat("Summarising bootstrapped signatures\n")
				summarise_whichSignatures <- function(x, signatures=signatures.cosmic, sampleName) {
					weights_mat <- do.call("rbind", lapply(x, function(y){y$weights}))
					tumor_mat <- do.call("rbind", lapply(x, function(y){y$tumor}))
					prod_mat <- do.call("rbind", lapply(x, function(y){ y$product}))
					diff_mat <- tumor_mat-prod_mat

					error <- sqrt(rowSums(diff_mat^2))
					corr <- unlist(lapply(1:nrow(prod_mat), function(n){cor(prod_mat[n,], tumor_mat[n,], method="spearman")}))

					weights <- matrix(colMeans(weights_mat),nrow=1)
					weightsSD <- apply(weights_mat, 2, sd)
					colnames(weights) <- colnames(weights_mat)
					rownames(weights) <- sampleName
					unknown <- 1 - sum(weights)

					tumor <- matrix(colMeans(tumor_mat), nrow=1)
					colnames(tumor) <- colnames(tumor_mat)
					rownames(tumor) <- sampleName

					product <- weights %*% as.matrix(signatures)
					diff <- tumor - product
					out <- list(weights, tumor, product, diff, unknown, weightsSD, weights_mat, tumor_mat, prod_mat, diff_mat, error, corr)
					names(out) <- c("weights", "tumor", "product", "diff", "unknown", "weightsSD",
						"weights_mat", "tumor_mat", "prod_mat", "diff_mat", "error", "corr")
					return(out)
				}
				ws2 <- lapply(unique(allmuts[,opt$sample_col]), function(s){
					summarise_whichSignatures(ws[grep(paste(s, "_", sep=""), names(ws))], sampleName=s)
				})
				names(ws2) <- unique(allmuts[,opt$sample_col])
				ws <- ws2
			}
		}
	
		plot_signature_bars_sanger_style <- function(x, main="", ylim=c(0,20), beside=T, border=NA, 
			cex.axis=1.5, cex.lab=1, las=2, ylab="Proportion of mutations",
			col=rep(c("cyan2", "black", "red", "grey", "chartreuse3", "lightpink1"), each=16),...) {
			barplot(x*100, beside=T, border=border, main=main, ylim=ylim, cex.axis=cex.axis, cex.lab=cex.lab, ylab=ylab,
			col=col, las=las, ...)
		}

		lapply(names(ws), function(sample) {
			pdf(paste(opt$outPrefix, ".pdf", sep=""), height=3, width=6)
			plot_signature_bars_sanger_style(ws[[sample]]$tumor)
			dev.off()
		})

		signatures <- do.call("rbind", lapply(ws, function(w) { w$weights }))

		if (!is.na(opt$num_iter)){
			if(opt$num_iter>=10){
				signaturesSD <- do.call("rbind", lapply(ws, function(w) { w$weightsSD }))
				colnames(signaturesSD) <- paste(colnames(signaturesSD), "SD", sep="")
				corr <- do.call("rbind", lapply(ws, function(w){ c(mean(w$corr), sd(w$corr))}))
				colnames(corr) <- c("Spearman_expected_vs_observed", "Spearman_SD")
				error <- do.call("rbind", lapply(ws, function(w){ c(mean(w$error), sd(w$corr))}))
				colnames(error) <- c("Error_expected_vs_observed", "Error_SD")
				signatures <- cbind(signatures, signaturesSD, corr, error)
			}	
		}
	}

	cat ("Computing mutational burden\n")
	sn <- c("A", "C", "G", "T")
	subs <- allmuts[,opt$ref_col] %in% sn & allmuts[,opt$alt_col] %in% sn
	mutcounts <- table(factor(allmuts[,opt$sample_col], levels=sample_levels), factor(subs, levels=c("TRUE", "FALSE")))

	mutrate <- data.frame(
		TOTAL = rowSums(mutcounts),
		SNV = mutcounts[,"TRUE"],
		INDEL = mutcounts[,"FALSE"])
	rownames(mutrate) <- sample_levels

	mutrate$SNV_prop = as.numeric(mutrate$SNV)/as.numeric(mutrate$TOTAL)
	mutrate$INDEL_prop = as.numeric(mutrate$INDEL)/as.numeric(mutrate$TOTAL)

	if(nrow(sigs) > 0) {
		signatures <- signatures[match(rownames(mutrate), rownames(signatures)),,drop=F]
		signatures <- cbind(mutrate, signatures)
	} else {
		signatures <- mutrate
	}
} else {
	mutrate <- data.frame(
		TOTAL = rep(0, length(sample_levels)),
		SNV = rep(0, length(sample_levels)),
		INDEL = rep(0, length(sample_levels)))
	rownames(mutrate) <- sample_levels

	mutrate$SNV_prop = rep(NA, length(sample_levels))
	mutrate$INDEL_prop = rep(NA, length(sample_levels))
	signatures <- mutrate
}
cat ("Writing and saving results\n")

write.table(signatures, file=paste(opt$outPrefix, ".txt", sep=""), sep="\t", col.names=NA, quote=F, na="")
if(exists("ws")) {
	save(signatures, ws, file=paste(opt$outPrefix, ".RData", sep=""))
} else {
	save(signatures, file=paste(opt$outPrefix, ".RData", sep=""))
}



