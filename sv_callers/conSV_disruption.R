suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("data.table"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("readxl"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));
suppressPackageStartupMessages(library("biomaRt"));
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"));

generateGRangesFromBedpe <- function(bedpe){	
	
	# For deletions (DEL): single range from START1 to END2 for true DEL < 10Mbs, otherwise two ranges at START1-END1 and START2-END2
	gr_del <-c( 
	with(subset(bedpe, CONS_SVTYPE == "deletion" & START2-START1 < 10000000),
	GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END2), svclass = CONS_SVTYPE)),
	with(subset(bedpe, CONS_SVTYPE == "deletion" & START2-START1 >= 10000000), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), svclass = CONS_SVTYPE)),
	with(subset(bedpe, CONS_SVTYPE == "deletion" & START2-START1 >= 10000000), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), svclass = CONS_SVTYPE))
	)
	# For inversions (INV): breakpoints at START1-END1 and START2-END2
	gr_inv <- c(
	with(subset(bedpe, CONS_SVTYPE == "inversion"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), svclass = CONS_SVTYPE)),
	with(subset(bedpe, CONS_SVTYPE == "inversion"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), svclass = CONS_SVTYPE)))
	# For tandem-duplication (DUP): single range from START1 to END1
	gr_dup <- with(
	subset(bedpe, CONS_SVTYPE == "tandem-duplication"),
	GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), svclass = CONS_SVTYPE)
	)
	# For translocations (TRA)
	gr_tra <- c(
	with(subset(bedpe, CONS_SVTYPE == "translocation"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), svclass = CONS_SVTYPE)),
	with(subset(bedpe, CONS_SVTYPE == "translocation"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), svclass = CONS_SVTYPE))
	)
	# For insertion (INS)
	gr_ins <- c(
	with(subset(bedpe, CONS_SVTYPE == "insertion"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), svclass = CONS_SVTYPE)),
	with(subset(bedpe, CONS_SVTYPE == "insertion"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), svclass = CONS_SVTYPE))
	)

	# Combine
	gr_all <- c(gr_del, gr_inv,gr_dup,gr_tra, gr_ins)
	return(gr_all)
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--oncokb_file", default = "conSV/cancerGeneList_edit_oncokb.xlsx",  type = "character", dest = "oncokb_file",
	            help = "Oncokb file"),
	make_option("--oncokb_minresources", default=2, type = "integer", dest = "oncokb_minresources",
	            help = "Oncokb min resources threshold"),
	make_option("--output", type = "character", dest = "output",
	            help = "Output file"),
	make_option("--input", type = "character", dest = "input", 
	            help = "Input file"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

opt <- parse_args(parser, positional_arguments = TRUE)$options;

if (is.null(opt$input)) {
    cat("Need input files\n");
    print_help(parser);
    stop();
} else if (is.null(opt$oncokb_file)) {
	cat("Need oncokb file\n");
	print_help(parser);
	stop();
} else if (is.null(opt$oncokb_minresources)) {
	cat("Need oncokb min resources\n");
	print_help(parser);
	stop();	
} else if (is.null(opt$oncokb_file)) {
	cat("Need oncokb file\n");
	print_help(parser);
	stop();
} else if (is.null(opt$output)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} 



disruptions_df <- data.frame(
  chrom = character(),
  start = integer(),
  end = integer(),
  width = integer(),
  sv_type = character(),
  cancer_genes = character(),
  cancer_promoters = character(),
  genes = character(),
  promoters = character()
)

# Generate GRanges from bedpe
onco <- read_excel(opt$oncokb_file) %>%  rename(
    "# of occurrence within resources (Column J-P)" = "resources","Hugo Symbol" = "Hugo_Symbol", "Gene Aliases" ="Gene_Aliases", "Gene Type" = "Gene_Type"
)

# filter for genes with at least 3 resources
onco_list <- onco %>%
  filter(resources >= opt$oncokb_minresources) %>%
  mutate(
    AllSymbol = ifelse(is.na(Gene_Aliases), 
                       Hugo_Symbol, 
                       paste(Hugo_Symbol, Gene_Aliases, sep = ","))
  ) %>%
  pull(AllSymbol) %>% 
  strsplit(",") %>% 
  unlist() %>% 
  unique()

# Annotate with genes and promoters
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
# promoter defined as 2000bp upstream and 50bp downstream of TSS
promoters <- promoters(genes, upstream = 2000, downstream = 50)

# Read bedpe
if (!file.exists(opt$input)) next
bedpe <- fread(opt$input, header = TRUE, stringsAsFactors = FALSE)
if (length(bedpe) == 0) next

# Generate GRanges from bedpe function defined earlier
gr_all = generateGRangesFromBedpe(bedpe)

#find overlaps
hits_gene <- as.data.frame(findOverlaps(gr_all, genes, ignore.strand=TRUE))
hits_promoters <- as.data.frame(findOverlaps(gr_all, promoters, ignore.strand=TRUE))
#annotate hits
hits_gene$GENE_SYMBOL <- biomaRt::select(org.Hs.eg.db, genes[hits_gene$subjectHits]$gene_id, "SYMBOL")$SYMBOL
hits_promoters$PROMOTER_SYMBOL <- biomaRt::select(org.Hs.eg.db, promoters[hits_promoters$subjectHits]$gene_id, "SYMBOL")$SYMBOL

if ((length(hits_gene) == 0)&&(length(hits_promoters) == 0)) next

#generate collapsed hits
collapsed_genes <- hits_gene %>%
	group_by(queryHits) %>%
	summarise(GENE_SYMBOL = paste(unique(GENE_SYMBOL), collapse = ","))
collapsed_promoters <- hits_promoters %>%
	group_by(queryHits) %>%
	summarise(PROMOTER_SYMBOL = paste(unique(PROMOTER_SYMBOL), collapse = ","))

#annotate gr
gr_all$GENE_SYMBOL <- ""
gr_all$PROMOTER_SYMBOL <- ""
gr_all$GENE_SYMBOL[collapsed_genes$queryHits] <- collapsed_genes$GENE_SYMBOL
gr_all$PROMOTER_SYMBOL[collapsed_promoters$queryHits] <- collapsed_promoters$PROMOTER_SYMBOL


#annotate w oncogenes
disrupted_genes = as.data.frame(gr_all)  %>% 
	mutate(
		onco_gene_hit = sapply(GENE_SYMBOL, function(x) {
			genes <- unlist(strsplit(x, ","))
			paste(intersect(genes, onco_list), collapse = ",")
		}),
		onco_promoter_hit = sapply(PROMOTER_SYMBOL, function(x) {
			promoters <- unlist(strsplit(x, ","))
			paste(intersect(promoters, onco_list), collapse = ",")
		}))

#rename columns to match disruptions_df
colnames(disrupted_genes) = c("chrom","start","end","width","strand","sv_type","genes","promoters","cancer_genes","cancer_promoters")


#bind results
disruptions_df = rbind(disruptions_df, disrupted_genes %>%  dplyr::select(colnames(disruptions_df)))

# Save
write.table(
disruptions_df,
opt$output, 
sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

