suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("dplyr"));


optList <- list(
	make_option("--outFile", default = NULL, type="character", action = "store", help = "output file"),
	make_option("--summaryType", default = "median", type="character", action="store", help = "summary Type: median, mean, max, min")
)

parser <- OptionParser(usage = "%prog [options] [star files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need star ReadsPerGene files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else if (!opt$summaryType %in% c("mean", "median", "max", "min")) {
	cat("summaryType has to be one of mean median max or min\n")
	print_help(parser)
	stop()
} else {
    starFiles <- arguments$args
}

dat <- lapply(starFiles, read.delim, as.is=T, header=F)
if (opt$summaryType=="median") { func=median
} else if (opt$summaryType=="mean") { func=mean 
} else if (opt$summaryType=="max") { func=max 
} else if (opt$summaryType=="min") { func=min 
}

res <- do.call("cbind", lapply(2:4, function(x) {
	apply(bind_cols(lapply(dat, function(y) { y[,x,drop=F]})),1,func,na.rm=T)}))
		
rownames(res) <- dat[[1]][,1]

write.table(res, file=opt$outFile, row.names=T, col.names=F, quote=F, na="", sep="\t")












