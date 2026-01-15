suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("data.table"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("gridExtra"));
suppressPackageStartupMessages(library("cowplot"));
suppressPackageStartupMessages(library("ShatterSeek"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

print(commandArgs(trailingOnly = TRUE))

optList <- list(
	make_option("--sample", type = "character", dest = "sample",
	            help = "Sample name"),
	make_option("--CN", type = "character", dest = "CN",
	            help = "Facet's cncf file"),
	make_option("--SV", type = "character", dest = "SV",
	            help = "Consensus SV file"),
	make_option("--outFile", type = "character", dest = "outFile",
	            help = "Output file"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

opt <- parse_args(parser, positional_arguments = TRUE)$options;

if (is.null(opt$CN)) {
    cat("Need CN file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$SV)) {
    cat("Need SV file\n");
    print_help(parser);
    stop();
}

#load CN data (Facets)
CN_df <- fread(opt$CN, sep='\t', header = TRUE)
CN_df = CN_df  %>% mutate(chrom = ifelse(chrom==23, "X", as.character(chrom)))
CN_data <- CNVsegs(chrom=CN_df$chrom,
				   start=CN_df$start,
				   end=CN_df$end,
				   total_cn=CN_df$tcn.em)

#load SV data (consensus SVs)
SV_df <- fread(opt$SV,sep='\t', header = TRUE)
#format chr name for ShatterSeek
SV_df <- SV_df %>% mutate(CHROM=gsub("chr","",as.character(CHROM))) %>% mutate(CHROM2=gsub("chr","",as.character(CHROM2)))
SV_df = SV_df  %>% filter(CHROM != "Y" & CHROM2 != "Y")
# add sv type as requested by shatterseek
SV_df <- SV_df  %>% mutate(CONS_SVTYPE=case_when(
	CHROM!=CHROM2 ~ "TRA",
	STRAND1=="+" & STRAND2=="-" ~ "DEL",
	STRAND1=="-" & STRAND2=="+" ~ "DUP",
	STRAND1=="+" & STRAND2=="+" ~ "h2hINV",
	STRAND1=="-" & STRAND2=="-" ~ "t2tINV"	
))
SV_data <- SVs(chrom1=as.character(SV_df$CHROM), 
			   pos1=as.numeric(SV_df$START1), 
			   chrom2=as.character(SV_df$CHROM2), 
			   pos2=as.numeric(SV_df$END2), 
			   strand1=as.character(SV_df$STRAND1),
			   strand2=as.character(SV_df$STRAND2),
			   SVtype=as.character(SV_df$CONS_SVTYPE))

#run ShatterSeek
chromothripsis <- shatterseek(
                SV.sample=SV_data,
                seg.sample=CN_data,
                genome="hg38")

#output results
summary = chromothripsis@chromSummary
#compute intra and inter chromosomial SVs
summary = summary  %>% mutate(intraSV = rowSums(select(.,number_DEL,number_h2hINV, number_t2tINV,number_DUP),na.rm = TRUE ))  
summary = summary  %>% mutate(interSV = clusterSize_including_TRA - intraSV)  

#criteria to define High or Low confidence chromothripsis events
summary =summary %>%
  mutate(
    is_high_1 =
      (intraSV >= 6 & #6 interleaved intra-chromosomal SVs
      max_number_oscillating_CN_segments_2_states >= 7 & #7 contiguousoscillating CN segments between 2 states
      pval_fragment_joins < 0.05 & #significant fragment joins (p-value < 0.05)
      (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05)), #either significant chromosome breakpoint enrichment (p-value < 0.05) or significant exponential distribution of breakpooints test (p-value < 0.05)

	is_high_2 =		
	  (intraSV >=3 & #3 interleaved intra-chromosomal SVs
	  interSV >=4 & # 4 inter-chromosomal SVs
	  max_number_oscillating_CN_segments_2_states >= 7 & #7 contiguous oscillating CN segments between 2 states
	  pval_fragment_joins < 0.05), #significant fragment joins test (p-value < 0.05)

    is_low =
      intraSV >= 6  & #6 interleaved intra-chromosomal SVs
      (max_number_oscillating_CN_segments_2_states >= 4 & max_number_oscillating_CN_segments_2_states <= 6) & #4-6 contiguous oscillating CN segments between 2 states
      pval_fragment_joins < 0.05 & #fragment joins test (p-value < 0.05)
      (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05), #either chromosome breakpoint enrichment (p-value < 0.05) or exponential distribution of breakpooints test (p-value < 0.05)

    ## cascade with priority
    cth = case_when(
      is_high_1 | is_high_2 ~ "High",
      !is_high_1 & !is_high_2 & is_low ~ "Low",
      TRUE ~ ""
    )
  )  

#format results table
res = summary %>% select(
	chrom, start, end, #cluster coordinates
	intraSV, interSV, clusterSize_including_TRA, # intra inter and tot cluster SVs
	max_number_oscillating_CN_segments_2_states, 
	pval_fragment_joins, 
	pval_exp_cluster,
	chr_breakpoint_enrichment,
	is_high_1, is_high_2, is_low, cth,
	inter_other_chroms, inter_other_chroms_coords_all) # also keep track of interchromosomial breakpoints

#write results
write.table(
	res,
	opt$outFile,
	sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE
)

#generate plots for chromothripsis events
chrom_to_plot = res %>% filter(cth=="Low" | cth=="High") %>% pull(chrom)

for(chr in chrom_to_plot){
plot_cth = plot_chromothripsis(chromothripsis, chr = chr)
plot_cth2 = arrangeGrob(plot_cth[[1]], plot_cth[[2]], plot_cth[[3]], plot_cth[[4]], nrow=4, ncol=1, heights = c(.2,.4,.4,.4))
pdf(paste0("shatterseek/",opt$sample,"/chr_",chr,"_cth.pdf"),  width = 5, height = 5)
print(plot_grid(plot_cth2))
dev.off()
}