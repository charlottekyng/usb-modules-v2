#################################################################
#1) formatting SV callers different outputs in standard BEDPE
#################################################################

#'Title: Format Brass
#' @description 
#'This function extract required fields from output vcf 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param brass Data frame with BRASS output
#'@export
format_brass = function(brass){
	brass2 = brass  %>% mutate(
	#SVTYPE=str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
	CHROM2=str_extract(ALT,"chr[0-9XYM]+"),
	POS2=str_extract(ALT,"(?<=:)\\d+"),
	MATEID=str_extract(INFO,"(?<=MATEID=)[^;]+"),
	SVTYPE=str_extract(INFO,"(?<=SVCLASS=)[^;]+"),
	BKDIST=str_extract(INFO,"(?<=BKDIST=)[^;]+"),
	BALS=ifelse(str_detect(INFO, "BALS"),str_extract(INFO,"(?<=BALS=)[^;]+"),""))

	valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

	brass2 <- brass2 %>%
		filter(CHROM %in% valid_chroms & CHROM2 %in% valid_chroms)

	return(brass2)

}

#'Title: Format Delly
#' @description 
#'This function extract required fields from output vcf 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param delly Data frame with delly output
#'@export
format_delly = function(delly){
	delly2 = delly  %>%  mutate(
    SVTYPE = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
    RDRATIO = str_extract(INFO, "(?<=RDRATIO=)[^;]+"),
    SVLEN = str_extract(INFO, "(?<=SVLEN=)[^;]+"),
    INSLEN = str_extract(INFO, "(?<=INSLEN=)[^;]+"),
    HOMLEN = str_extract(INFO, "(?<=HOMLEN=)[^;]+"),
    CHROM2 = if_else(str_detect(INFO, "CHR2="), str_extract(INFO, "(?<=CHR2=)[^;]+"), CHROM),
    END = str_extract(INFO, "(?<=\\bEND=)[0-9]+"),
    POS2 = if_else(str_detect(INFO, "POS2="), str_extract(INFO, "(?<=POS2=)[^;]+"), END))


	valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

	delly2 <- delly2 %>%
		filter(CHROM %in% valid_chroms & CHROM2 %in% valid_chroms)

	return(delly2)
}

#'Title: Format GRIDSS2
#' @description 
#'This function extract required fields from output vcf 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param gridss Data frame with GRIDSS2 output
#'@export
format_gridss = function(gridss){
	gridss2 = gridss %>% 
	mutate(
		SVTYPE = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
		EVENT = str_extract(INFO, "(?<=EVENT=)[^;]+"),
		MATEID = str_extract(INFO, "(?<=MATEID=)[^;]+"),
		CHROM2=ifelse(str_detect(ALT,"chr"),str_extract(ALT,"chr[0-9XYM]+"),CHROM),
		POS2 = ifelse(str_detect(ALT, "chr"),
					as.numeric(str_extract(ALT,"(?<=:)\\d+")),
					ifelse(nchar(ALT)>nchar(REF),POS + (nchar(ALT)-nchar(REF)),POS + (nchar(REF)-nchar(ALT)))),
		LOCAL_LINKED_BY=str_extract(INFO,"(?<=LOCAL_LINKED_BY=)[^;]+"),
		REMOTE_LINKED_BY=str_extract(INFO,"(?<=REMOTE_LINKED_BY=)[^;]+"))
	
	valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

	gridss2 <- gridss2 %>%
		filter(CHROM %in% valid_chroms & CHROM2 %in% valid_chroms)

	return(gridss2)
}

#'Title: Format Manta
#' @description 
#'This function extract required fields from output vcf 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param manta Data frame with manta output
#'@export
format_manta = function(manta){
	manta2 = manta %>% 
	mutate(
		CHROM2=ifelse(str_detect(ALT,"chr"),str_extract(ALT,"chr[0-9XYM]+"),CHROM),
		POS2 = ifelse(str_detect(ALT, "chr"),
			as.numeric(str_extract(ALT,"(?<=:)\\d+")),
			as.numeric(str_extract(INFO,"(?<=END=)[^;]+"))),
		SVLEN = str_extract(INFO, "(?<=SVLEN=)[^;]+"),
		SVTYPE = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
		EVENT = str_extract(INFO, "(?<=EVENT=)[^;]+"),
		HOMLEN =str_extract(INFO, "(?<=HOMLEN=)[^;]+"),
		SVINSLEN =str_extract(INFO, "(?<=SVINSLEN=)[^;]+"),
		MATEID = str_extract(INFO, "(?<=MATEID=)[^;]+"))

	valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

	manta2 <- manta2 %>%
		filter(CHROM %in% valid_chroms & CHROM2 %in% valid_chroms)

	return(manta2)
}

#'Title: Format svaba
#' @description 
#'This function extract required fields from output vcf 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param svaba Data frame with svaba output
#'@export
format_svaba = function(svaba){
	svaba2 = svaba %>% 
		mutate(
		CHROM2=ifelse(str_detect(ALT,"chr"),str_extract(ALT,"chr[0-9XYM]+"),CHROM),
		POS2 = ifelse(str_detect(ALT, "chr"),
			as.numeric(str_extract(ALT,"(?<=:)\\d+")),
			ifelse(nchar(ALT)>nchar(REF),POS + (nchar(ALT)-nchar(REF)),POS + (nchar(REF)-nchar(ALT)))),
		SCTG= str_extract(INFO, "(?<=SCTG=)[^;]+"),
		SVLEN = str_extract(INFO, "(?<=SPAN=)[^;]+"),
		SVTYPE = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
		MATEID = str_extract(INFO, "(?<=MATEID=)[^;]+")) 

	valid_chroms <- paste0("chr", c(1:22, "X", "Y"))

	svaba2 <- svaba2 %>%
		filter(CHROM %in% valid_chroms & CHROM2 %in% valid_chroms)
	
	return(svaba2)
}


#'Title: Format output in bedpe
#' @description 
#'This function format each caller output in standard bedpe
#'strands and positions, in particular
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param df from formatted function
#'@export
format_caller <- function(df) {
	df2 = df %>%   
	mutate(parsed = Map(parse_strands_from_alt, ALT, REF)) %>%
	unnest_wider(parsed, names_sep = "") %>%
	rename(STRAND1 = parsedstrand1, STRAND2 = parsedstrand2) %>% 
	mutate(START1 = as.numeric(POS)-1, END1 = as.numeric(POS), START2 = as.numeric(POS2)-1, END2 = as.numeric(POS2))
	# the output bed file should have 0-based starts and 1-based ends
	return(df2)
}

#################################################################
#2) parsing STRAND information from ALT field
#################################################################

#'Title: Parse strand from ALT/SVTYPE field
#' @description 
#'This function infer strand from ALT field 
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param alt is the ALT field which can have shorthand notation (e.g. <DEL>), square brackets or genomic sequence
#'@param ref  REF field which is the reference
#'@export

parse_strands_from_alt <- function(alt, ref) {

	# Handle symbolic alleles first
	if (alt == "<DEL>") {
	return(c(strand1 = "+", strand2 = "+"))
	} else if (alt == "<INV>") { # only ++ because BND should have been 'normalised' such as pos1 is smaller than pos2
	return(c(strand1 = "+", strand2 = "-"))
	} else if (alt == "<DUP>") {
	return(c(strand1 = "-", strand2 = "-"))
	} else if (alt == "<DUP:TANDEM>") {
	return(c(strand1 = "-", strand2 = "-"))
	} else if (alt == "<INS>") { #not sure about this
	return(c(strand1 = "+", strand2 = "+"))
	}
	# Handle BND alleles
	# According to VCF spec:
	# - A[chr:pos[  => this base is joined after the other one (strand1 = +, strand2 = +)
	# - A]chr:pos]  => (strand1 = +, strand2 = -)
	# - ]chr:pos]A  => (strand1 = -, strand2 = -)
	# - [chr:pos[A  => (strand1 = -, strand2 = +)
	# these have been confirmed from BRASS, which returns both bed and vcf
	if (grepl("\\[|\\]", alt)) {
		if (grepl("^[A-Za-z]+\\[", alt)) { # deletion
		strand1 <- "+"
		strand2 <- "+"
		} else if (grepl("^[A-Za-z]+\\]", alt)) { # inversion
		strand1 <- "+"
		strand2 <- "-"
		} else if (grepl("^\\]", alt)) { #tandem duplication
		strand1 <- "-"
		strand2 <- "-"
		} else if (grepl("^\\[", alt)) { # inversion
		strand1 <- "-"
		strand2 <- "+"
		}
		return(c(strand1 = strand1, strand2 = strand2))
	}else{ # these are the simple REF/ALT alleles
		strand1 <- "+"
		strand2 <- "+"
		return(c(strand1 = strand1, strand2 = strand2))
	}
	
}


#################################################################
#3) homogenise SVTYPE field
#################################################################

#'Title: Format and harmonize SVtype field
#' @description 
#'This function homogenize SVTYPE field
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param SVTYPE from formatted function
#'@param chr1 first chromosome coordinate
#'@param chr2 second chromosome coordinate
#'@param strand1 of first BND
#'@param strand2 of second BND
#'@export
rename_SVclass <- function(SVTYPE, chr1, chr2, ref, alt, strand1 = NA, strand2 = NA) {

	if (chr1 != chr2) {
	return("translocation")
	}

	if (SVTYPE %in% c("DUP", "DEL", "INV","INS", "TRA")) {
		return(dplyr::case_match(
		SVTYPE,
		"DUP" ~ "tandem-duplication",
		"DUP:TANDEM" ~ "tandem-duplication",
		"DEL" ~ "deletion",
		"INV" ~ "inversion",
		"INS" ~ "insertion",
		"TRA" ~ "translocation"
		))
	}

	  if (grepl("^[A-Za-z]$", ref) & grepl("^[\\.A-Za-z\\.]+$", alt)) return("insertion")
    if (grepl("^[A-Za-z\\.]+$", ref) & grepl("^[\\.A-Za-z\\.]$", alt)) return("deletion")

    if (strand1 == "+" & strand2 == "+") return("deletion") 
    if (strand1 == "-" & strand2 == "-") return("tandem-duplication")
    if (strand1 == "+" & strand2 == "-") return("inversion")
    if (strand1 == "-" & strand2 == "+") return("inversion")

}

###################################################################
#4) create a consenus from different callers based on range overlap
###################################################################


#'Title: Make consensus from different callers
#' @description 
#'This function to create consensus from different callers
#' @import data.table
#' @import dplyr
#' @import stringr
#' @import tidyr
#'@param df concatenate bedpe format SV call results
#'@param max_dist amx distance to consider to events (default 200bp)
#'@export
make_sv_consensus <- function(df, max_dist = 200) {
  required_cols <- c("CHROM", "START1", "END1", "CHROM2", "START2", "END2", "SVTYPE", "STRAND1", "STRAND2", "sv_method")
  if (!all(required_cols %in% names(df))) {
    stop(paste("Missing required columns:", paste(setdiff(required_cols, names(df)), collapse = ", ")))
  }

  df <- df %>%
    arrange(CHROM, START1, CHROM2, START2)

  cluster <- integer(nrow(df))
  cluster_id <- 1

  for (i in seq_len(nrow(df))) {
    if (cluster[i] != 0) next

    cluster[i] <- cluster_id

    close_idx <- which(
      df$CHROM == df$CHROM[i] &
      df$CHROM2 == df$CHROM2[i] &
      abs(df$START1 - df$START1[i]) <= max_dist &
      abs(df$START2 - df$START2[i]) <= max_dist
    )

    cluster[close_idx] <- cluster_id
    cluster_id <- cluster_id + 1
  }

  consensus <- df %>%
    mutate(cluster = cluster) %>%
    group_by(cluster) %>%
    summarise(
      CHROM   = first(CHROM),
      START1  = ceiling(median(START1)),
      END1    = ceiling(median(END1)),
      CHROM2  = first(CHROM2),
      START2  = ceiling(median(START2)),
      END2    = ceiling(median(END2)),
      STRAND1 = names(sort(table(STRAND1), decreasing = TRUE))[1],
      STRAND2 = names(sort(table(STRAND2), decreasing = TRUE))[1],
	  CONS_SVTYPE  = names(sort(table(SVTYPE), decreasing = TRUE))[1],
	  SVTYPE_counts = paste0(names(table(SVTYPE)), "(", table(SVTYPE), ")", collapse = ", "),
      sv_method = paste(unique(sv_method), collapse = ","),
      n_caller_support = str_count(sv_method, ",") + 1,
	  n_BND_support = n(),
      .groups = "drop"
    ) %>%
    arrange(CHROM, START1)

  return(consensus)
}
