suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(rtracklayer));

optList <- list(
	make_option('--inputRSEMFile', action='store', default = 'all.genes.expected_count.results', help = 'input RSEM file to be normalized'),
	make_option('--gtf', action='store', default = NULL, help = 'GTF annotation file if gene subset is required'),
	make_option('--geneBiotype', action='store', default = "protein_coding", help = 'gene biotype/s to include'),
	make_option('--outputFile', action='store', default = NULL, help = 'output file'),
	make_option('--normalizationMethod', action='store', default = 'uq', help = 'which normalization method to use'),
	make_option('--threshold_for_uq', action='store', type= 'integer', default = 1000, help = 'the threshold for UQ normalization'))

parser <- OptionParser(usage = "%prog", option_list = optList);
opt <- parse_args(parser, positional_arguments = F);
#opt <- arguments$options;
#print(opt)
#print(arguments)

if (is.null(opt$inputRSEMFile)) {
    cat("Need input RSEM file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outputFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
}

if (!is.null(opt$gtf)){
	gtf <- import(opt$gtf)
	if (!is.null(opt$geneBiotype)){
		gtf <- gtf[which(gtf$gene_biotype %in% opt$geneBiotype),]
	} else { cat("No geneBiotype provided, using all genes in the GTF.\n") }
} else {
	cat ("No GTF provided, using all genes\n")
}
rsem <- read.delim(opt$inputRSEMFile, as.is=T, row.names=1)
rsem <- rsem[which(rownames(rsem) %in% gtf$gene_id),]


## This performs quantile-normalization
## The default is upper-quartile normalization to a fixed value of 1000,
## as described in the TCGA LIHC paper
quantile_normalize_rsem <- function(obj, quantile_val=0.75, fixed_val=1000) {
	norm_factor = fixed_val/apply(obj,2,quantile,quantile_val)
	t(t(obj)*norm_factor)
}


if (opt$normalizationMethod=='uq') {
	if (!is.null(opt$threshold_for_uq)) {
		rsem_norm <- quantile_normalize_rsem(rsem, fixed_val=opt$threshold_for_uq)
	} else {
		cat ('No other normalization implemented at the moment...')
	}
} else { cat ('No other normalization implemented at the moment...')}

write.table(rsem_norm, file=opt$outputFile, sep="\t", col.names=NA, quote=F, na="")


