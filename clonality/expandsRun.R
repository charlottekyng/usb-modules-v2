suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("expands"));

optList <- list(
        make_option("--mutations", default = NULL, help = "mutations file"),
        make_option("--segs", default = NULL, help = "segments file"),
	make_option("--numcores", default = 4, help = "number of cores")
        make_option("--outPrefix", default = NULL, help = "output prefix")        
)

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = F);
opt <- arguments$options;

if (is.null(opt$outPrefix)) {
    cat("Need output file prefix\n");
    print_help(parser);
    stop();
} else if (is.null(opt$mutations)) {
    cat("Need output mutations prefix\n");
    print_help(parser);
    stop();
} else if (is.null(opt$segs)) {
    cat("Need output segs prefix\n");
    print_help(parser);
    stop();
}

expands <- runExPANdS(SNV=opt$mutations, CBS=opt$segs, nc=opt$numcores, snvF=opt$outPrefix)
save(expands, file=paste(opt$outPrefix, ".RData", sep="")
