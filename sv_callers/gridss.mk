#template for gridss
MUT_CALLER = gridss

# Run gridss on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/gridss.$(NOW)
PHONY += gridss

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

gridss : gridss_vcfs

gridss_dir:
	mkdir -p gridss

gridss_vcfs : gridss_dir
	$(foreach pair,$(SAMPLE_PAIRS),gridss/$(tumor.$(pair)).somatic.vcf)

# PoN generation rules
define gridss-pon

#create folder for PoN
gridss_pon_dir: gridss_dir
	mkdir -p gridss/pondir

# Main call for normal pairs needed for PoN
gridss/pondir/%.normal.vcf : bam/%.bam gridss_pon_dir
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(findstring hg38,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed,) \
	$$(if $$(findstring b37,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed $$(BED_DIR)/ENCFF001TDO.bed,) \
	-r $$(REF_FASTA) \
	-o $$@ \
	$$<")

#creating the PoN
gridss/pondir/gridss_pon_breakpoint.bedpe gridss/pondir/gridss_pon_single_breakend.bed: gridss/pondir/*.normal.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
		$$(GRIDSS) gridss.GeneratePonBedpe \
		$$(ls -1 gridss/pondir/*.normal.vcf | awk ' { print "INPUT=" $0 }' | head -$n) \
		O=gridss/pondir/gridss_pon_breakpoint.bedpe \
		SBO=gridss/pondir/gridss_pon_single_breakend.bed \
		REFERENCE_SEQUENCE=$ref_genome")

endef

# Somatic analysis depends on PoN files
define gridss-tumor-normal

# Main call for tumor-normal pairs
gridss/$1.vcf : bam/$2.bam bam/$1.bam gridss/pondir/gridss_pon_breakpoint.bedpe gridss/pondir/gridss_pon_single_breakend.bed
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(or $$(findstring hg38,$$(REF)), $$(findstring b37,$$(REF))),\
	    -b $$(if $$(findstring hg38,$$(REF)),$$(BED_DIR)/ENCFF356LFX.bed,\ 
		   $$(BED_DIR)/ENCFF001TDO.bed),) \
	-r $$(REF_FASTA) \
	-o $$@ \
	$$^")

# Somatic filtering
gridss/$1.somatic.vcf : gridss/$1.vcf gridss/pondir/gridss_pon_breakpoint.bedpe gridss/pondir/gridss_pon_single_breakend.bed
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss_somatic_filter \
  --pondir gridss/pondir/ \
  --input $$< \
  --output $$@ \
  --fulloutput gridss/$1.full.vcf.gz \
  --scriptdir $(dirname $(which gridss_somatic_filter)) \
  -n 1 \
  -t 2")

endef

# Apply the PoN and tumor-normal rules
$(foreach sample,$(PANEL_OF_NORMAL_SAMPLES),$(eval $(call gridss-pon,$(sample))))
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call gridss-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))