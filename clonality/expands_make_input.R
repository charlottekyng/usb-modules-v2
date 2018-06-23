suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--type", default = NULL, help = "mutations or cna")
        
)

parser <- OptionParser(usage = "%prog [options] mutations/CNA file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input file\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

dat <- read.delim(files[1], as.is=T)

if (opt$type == "mutations") {
	out <- dat[,c("CHROM", "POS", "TUMOR.FA")]
	colnames(out) <- c("chr", "startpos", "AF_Tumor")
	out$PN_B <- 0
} else if (opt$type == "cna") {
	out <- dat[,c("chrom", "start", "end", "expected_cn")]
	colnames(out) <- c("chr", "startpos", "endpos", "CN_Estimate")
}

write.table(out, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

