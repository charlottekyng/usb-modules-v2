include usb-modules-v2/Makefile.inc

LOGDIR ?= log/sv_pon.$(NOW)

PHONY += sv_pon

.DELETE_ON_ERROR:
.PHONY : $(PHONY)

#setup reference to run just once for all samples
gridss_ref:
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	-s setupreference \
	-r $$(REF_FASTA)"

sv_pon : gridss/pondir/gridss_pon_breakpoint.bedpe gridss/pondir/gridss_pon_single_breakend.bed

define gridss-pon
# Main call for normal pairs needed for PoN
gridss/pondir/%.normal.vcf : bam/%.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(findstring hg38,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed,) \
	$$(if $$(findstring b37,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed $$(BED_DIR)/ENCFF001TDO.bed,) \
	-r $$(REF_FASTA) \
	-o $$@ \
	$$<")"

#creating the main PoN
gridss/pondir/gridss_pon_breakpoint.bedpe gridss/pondir/gridss_pon_single_breakend.bed: gridss/pondir/*.normal.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
		$$(GRIDSS) GeneratePonBedpe \
		-INPUT $$< \
		-OUTPUT_BEDPE gridss/pondir/gridss_pon_breakpoint.bedpe \
		-OUTPUT_BED gridss/pondir/gridss_pon_single_breakend.bed \
		-REFERENCE_SEQUENCE $$(REF_FASTA)"

endef

$(foreach normal,$(PANEL_OF_NORMAL_SAMPLES), $(eval $(call gridss-pon,$(normal))))