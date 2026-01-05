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

# Genes to consider for fusion events
landscape = c(
"KRAS", "PIK3CA", "FBXW7" ,"GNAS", "ERBB2", "MAP2K4", "ELF3", "ERBB3", "EPHA6",
"TP53", "ATM", "RECQL4", "CDKN2A", "CDKN2B", 
"BRCA1","BRCA2","PALB2","BARD1",
"APC", "RNF43", "SOX9","CTNND1", "CTNNB1", "AXIN1", 
"SMAD4", "ACVR2A", "TGFBR2", "ACVR2B", "TGFBR1", 
"ARID2","KDM6A", "ARID1A","PBRM1", 
"LOXHD1", "SLITRK5", "CNTN4","ACVR1B","PTEN","EPHA3","BRAF",
"MLH1","MSH2","MSH3","MSH6","PMS2","POLD1", "POLD2","EPCAM",
"FHIT", "MTAP", "PRKN",  "PTPRD", "CDK12", "LRP1B", "SMYD3",
"MYC", "GATA6",
"ALK","EGFR","FGFR1","FGFR2","FGFR3","FGR","MET","NTRK1","NTRK2","NTRK3","PDGFRA","PIK3CA","PKN1","PRKCA","PRKCB","RAF1","RET", "ROS1" 
)
# Function to check if fusion is valid
valid_fusions <- function(x, landscape) {
				x <- trimws(x)                      # remove surrounding spaces
				all(x != "") &&                     # no empty names
				length(unique(x)) > 1  &&           # check 2 different genes
				any(x %in% landscape)				# at least one gene in landscape
				}
				
# Function to write empty output and quit in case of no fusions
write_empty_output <- function(outfile) {
  write.table(
    data.frame(),
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  quit(save = "no", status = 0)
}

# Function to generate GRanges from bedpe
generateGRangesFromBedpe <- function(bedpe){	
	# For deletions (DEL)
	gr_del <-c( 
	with(subset(bedpe, CONS_SVTYPE == "deletion"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), strand=STRAND1, svclass = CONS_SVTYPE,cluster = cluster)),
	with(subset(bedpe, CONS_SVTYPE == "deletion"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), strand=STRAND2, svclass = CONS_SVTYPE,cluster = cluster))
	)
	# For inversions (INV): 
	gr_inv <- c(
	with(subset(bedpe, CONS_SVTYPE == "inversion"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), strand=STRAND1, svclass = CONS_SVTYPE,cluster = cluster)),
	with(subset(bedpe, CONS_SVTYPE == "inversion"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), strand=STRAND2, svclass = CONS_SVTYPE,cluster = cluster)))
	# For translocations (TRA)
	gr_tra <- c(
	with(subset(bedpe, CONS_SVTYPE == "translocation"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), strand=STRAND1, svclass = CONS_SVTYPE,cluster = cluster)),
	with(subset(bedpe, CONS_SVTYPE == "translocation"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), strand=STRAND2, svclass = CONS_SVTYPE,cluster = cluster))
	)
	# For insertion (INS)
	gr_ins <- c(
	with(subset(bedpe, CONS_SVTYPE == "insertion"), GRanges(seqnames = CHROM, ranges = IRanges(start = START1, end = END1), strand=STRAND1, svclass = CONS_SVTYPE,cluster = cluster)),
	with(subset(bedpe, CONS_SVTYPE == "insertion"), GRanges(seqnames = CHROM2, ranges = IRanges(start = START2, end = END2), strand=STRAND2, svclass = CONS_SVTYPE,cluster = cluster))
	)

	# Combine
	gr_all <- c(gr_del,gr_ins, gr_inv,gr_tra)
	return(gr_all)
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

# Command line options: just input and output file names
optList <- list(
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
} else if (is.null(opt$output)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} 

# Read bedpe
if (!file.exists(opt$input)) write_empty_output(opt$output)
bedpe <- fread(opt$input, header = TRUE, stringsAsFactors = FALSE)
if (length(bedpe) == 0) write_empty_output(opt$output)

### create df
fusion_df <- data.frame(
  chrom = character(),
  start = integer(),
  end = integer(),
  gene1=character(),
  strand1=character(), 
  chrom2 = character(),
  start2 = integer(),
  end2 = integer(),
  gene2= character(),
  strand2 =character(), 
  sv_type = character()
)

# Annotate with genes and promoters
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Generate GRanges from bedpe function defined earlier
gr_all = generateGRangesFromBedpe(bedpe)

#find overlaps
hits_gene <- as.data.frame(findOverlaps(gr_all, genes, ignore.strand=TRUE))
#annotate hits
hits_gene$GENE_SYMBOL <- biomaRt::select(
	org.Hs.eg.db,
	keys = genes[hits_gene$subjectHits]$gene_id,
	columns = "SYMBOL"
)$SYMBOL

if (length(hits_gene) == 0) write_empty_output(opt$output)

#generate collapsed hits
collapsed_genes <- hits_gene %>%
	group_by(queryHits) %>%
	summarise(GENE_SYMBOL = paste(unique(GENE_SYMBOL), collapse = ","))

#annotate gr
gr_all$GENE_SYMBOL <- ""
gr_all$GENE_SYMBOL[collapsed_genes$queryHits] <- collapsed_genes$GENE_SYMBOL

#split by cluster and check valid fusions
gr_list <- split(gr_all, gr_all$cluster)

# Filter valid fusions
gr_valid_by_cluster <- lapply(gr_list, function(gr) {
  keep <- sapply(
    strsplit(gr$GENE_SYMBOL, ","),
    valid_fusions,
    landscape = landscape
  )
  gr[keep]
})

# keep clusters with exactly 2 breakpoints
gr_valid_by_cluster <- gr_valid_by_cluster[
	lengths(gr_valid_by_cluster) == 2
]

if (!length(gr_valid_by_cluster)) write_empty_output(opt$output)

# Build fusion calls
for (cl in names(gr_valid_by_cluster)) {
	
	gr_cl <- gr_valid_by_cluster[[cl]]

	fusion_genes <- data.frame(
		cluster = cl,
		chrom1 = as.character(seqnames(gr_cl[1])),
		start1 = start(gr_cl[1]),
		end1 = end(gr_cl[1]),
		gene1 = gr_cl$GENE_SYMBOL[1],
		strand1 = as.character(strand(gr_cl[1])),
		chrom2 = as.character(seqnames(gr_cl[2])),
		start2 = start(gr_cl[2]),
		end2 = end(gr_cl[2]),
		gene2 = gr_cl$GENE_SYMBOL[2],
		strand2 = as.character(strand(gr_cl[2])),
		sv_type = gr_cl$svclass[1],
		stringsAsFactors = FALSE
	)

	fusion_df <- rbind(fusion_df, fusion_genes)
}

# Save
write.table(
fusion_df,
opt$output, 
sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

