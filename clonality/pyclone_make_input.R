suppressPackageStartupMessages(library("optparse"));

optList <- list(
	make_option("--outFile", default = NULL, help = "output file"),
	make_option("--totalOrAllelic", default = "allelic", help ="total or allelic copy number")
)

parser <- OptionParser(usage = "%prog [options] mutation_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input mutation table files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

make_pyclone_input <- function(muts, sample_col="TUMOR_SAMPLE", 
	mutation_id_col=NULL, gene_col="GENE", aa_col="HGVS_P", cdna_col="HGVS_C",
	chrom_col="CHROM", pos_col="POS", ref_col="REF", alt_col="ALT",
	ref_counts_col=NULL, var_counts_col=NULL, ad_col="TUMOR.AD", MAF_col="TUMOR.FA", dp_col="TUMOR.DP", 
	major_cn_col=NULL, minor_cn_col="facetsLCN_EM", total_cn_col="facetsTCN_EM", 
	total_or_allelic="allelic") {

	empty_table <- 	matrix(nrow=0, ncol=6)
	colnames(empty_table) <- c("mutation_id", "ref_counts", "var_counts", "normal_cn", "major_cn", "minor_cn")

	if (total_or_allelic =="allelic") { 
		if (!minor_cn_col %in% colnames(muts)){stop ("minor_cn_col not found")}
		muts <- muts[which(muts[,minor_cn_col] %in% seq(0,1000)),,drop=F] 
	}
	if(nrow(muts)==0) { return(empty_table)}

	if (!is.null(ref_counts_col) & !is.null(var_counts_col)) {
		if (ref_counts_col %in% colnames(muts) & var_counts_col %in% colnames(muts)) {
			cat("Using ref_counts_col and var_counts_col\n")
			ref_counts = muts[,ref_counts_col] 
			var_counts = muts[,var_counts_col]
		} else { cat("ref_counts_col and var_counts_col not valid. Trying something else\n") }
	} else if (!is.null(ad_col)) {
		if (ad_col %in% colnames(muts)) { 
			cat("Using ad_col\n")	
			tt <- lapply(muts[,ad_col], strsplit, split=",")
			ref_counts <- unlist(lapply(tt, function(x){x[[1]][1]}))
			var_counts <- unlist(lapply(tt, function(x){x[[1]][2]}))
		} else { cat ("ad_col not in column names. Trying something else.") }
	} else if (!is.null(MAF_col) & !is.null(dp_col)) {
		if (MAF_col %in% colnames(muts) & dp_col %in% colnames(muts)) {
			cat("Using MAF_col and dp_col \n")
			muts[,MAF_col] <- gsub("%", "", muts[,MAF_col]); 
			if (any(as.numeric(muts[,MAF_col])>1)) { muts[,MAF_col] <- as.numeric(muts[,MAF_col])/100 }
			var_counts <- round(as.numeric(muts[,MAF_col])*as.numeric(muts[,dp_col])) 
			ref_counts <- as.numeric(muts[,dp_col])-var_counts 
		} else { cat ("MAF_col and dp_col not in column names. Trying something else.") }
	} else { stop ("ref_counts_col and var_counts_col, ad_col, or MAF_col and dp_col required.\n") }
		
	if (is.null(major_cn_col) & total_or_allelic =="allelic") { 
		if (!total_cn_col %in% colnames(muts) | !minor_cn_col %in% colnames(muts)){
			stop("total_cn_col or minor_cn_col not found")}
		major_cn <- as.numeric(muts[,total_cn_col])-as.numeric(muts[,minor_cn_col]) 
	}
	if (is.null(major_cn_col) & total_or_allelic =="total") { 
		if (!total_cn_col %in% colnames(muts)) { stop ("total_cn_col not found")}
		major_cn <- muts[,total_cn_col] 
	}
	if (total_or_allelic =="allelic") { 
		if (!minor_cn_col %in% colnames(muts)) { stop ("minor_cn_col not found")}
		minor_cn <- muts[,minor_cn_col]
	}
	if (total_or_allelic =="total") { minor_cn <- 0 }

	if (!is.null(mutation_id_col)) { rownames(muts) <- muts[,mutation_id_col] 
	} else {
		id_fields <- c(chrom_col, pos_col, ref_col, alt_col, gene_col, aa_col, cdna_col)
		id_fields <- id_fields[which(id_fields %in% colnames(muts))]
		cat ("mutation_id_col not provided. Constructing ids from: \n")
		cat (toString(id_fields)); cat("\n")
		x <- muts[,id_fields]
		rownames(muts) <- gsub(" ", "", gsub(", ", "_", apply(x, 1, toString)), fixed=T)
	}

	res <- cbind(rownames(muts), ref_counts, var_counts, 2, major_cn, minor_cn)
	colnames(res) <- c("mutation_id", "ref_counts", "var_counts", "normal_cn", "major_cn", "minor_cn")
	res <- res[which(as.numeric(res[,"major_cn"]) > 0),,drop=F]
	if(nrow(muts)==0) { return(empty_table)}

	res
}

dat <- read.delim(files[1], as.is=T)
write.table(make_pyclone_input(dat), file=opt$outFile, sep="\t", row.names=F, quote=F, na="")
	
	
	
	
	
	
