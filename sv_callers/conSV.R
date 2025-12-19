suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("data.table"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("conSV"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

print(commandArgs(trailingOnly = TRUE))

optList <- list(
	make_option("--min_callers", type = "integer", dest = "min_callers",
	            help = "Minimum number of callers supporting an SV"),
	make_option("--slope", type = "integer", dest = "slope",
	            help = "Distance allowed overlap in bp"),
	make_option("--sv_callers", type = "character", dest = "sv_callers",
	            help = "Dot-separated list of SV callers"),
	make_option("--output", type = "character", dest = "output",
	            help = "Output file"),
	make_option("--input", type = "character", dest = "input", 
	            help = "input samples (separated by dot)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

opt <- parse_args(parser, positional_arguments = TRUE)$options;

if (is.null(opt$input)) {
    cat("Need input files\n");
    print_help(parser);
    stop();
} else if (is.null(opt$min_callers)) {
    cat("Need min-callers\n");
    print_help(parser);
    stop();
} else if (is.null(opt$slope)) {
    cat("Need slope\n");
    print_help(parser);
    stop();
} else if (is.null(opt$sv_callers)) {
    cat("Need sv-callers\n");
    print_help(parser);
    stop();
} else if (is.null(opt$output)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} 

sv_callers <- strsplit(opt$sv_callers, "\\.")[[1]]
sv_list <- list()

for (caller in sv_callers) {
  path <- paste0("conSV/variantExtr/", opt$input, ".", caller, ".varExtr.vcf")
  if (!file.exists(path)) next

  dt <- fread(path, sep = "\t", header = TRUE, skip = "#CHROM")
  dt <- dt %>% rename(CHROM = "#CHROM")

  formatted <- switch(
    caller,
    brass  = format_brass(dt),
    delly  = format_delly(dt),
    gridss = format_gridss(dt),
    manta  = format_manta(dt),
    svaba  = format_svaba(dt)
  )

  if (nrow(formatted) == 0) next

  sv_list[[caller]] <- format_caller(formatted) %>%
    mutate(
      SVTYPE = mapply(rename_SVclass, SVTYPE, CHROM, CHROM2, STRAND1, STRAND2),
      sv_method = caller
    ) %>%
    select(CHROM, START1, END1, CHROM2, START2, END2, SVTYPE, STRAND1, STRAND2,sv_method)
}

if (length(sv_list) == 0) {
  message("No consensus SVs found for ", opt$input)
  quit(status = 0)
}

# Concatenate
sv_all <- bind_rows(sv_list)

# Generate consensus
consensus_sv <- make_sv_consensus(sv_all, max_dist = opt$slope)

# Save
write.table(
consensus_sv  %>% filter(n_caller_support >= opt$min_callers),
opt$output, 
sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)