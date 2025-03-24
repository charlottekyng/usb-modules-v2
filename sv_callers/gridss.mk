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


define gridss-tumor-normal

# Main call
gridss/$1.vcf : bam/$2.bam bam/$1.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(or $$(findstring hg38,$$(REF)), $$(findstring b37,$$(REF))),\
	    -b $$(if $$(findstring hg38,$$(REF)),$$(BED_DIR)/ENCFF356LFX.bed,\
		   $$(BED_DIR)/ENCFF001TDO.bed),) \
	-r $$(REF_FASTA) \
	-o $$@ \
	$$^")


#creating the PoN
mkdir -p pondir
java -Xmx8g \
	-cp ~/dev/gridss/target/gridss-2.10.2-gridss-jar-with-dependencies.jar \
	gridss.GeneratePonBedpe \
	$(ls -1 *.vcf.gz | awk ' { print "INPUT=" $0 }' | head -$n) \
	O=pondir/gridss_pon_breakpoint.bedpe \
	SBO=pondir/gridss_pon_single_breakend.bed \
	REFERENCE_SEQUENCE=$ref_genome

# Somatic filtering
gridss_somatic_filter \
  --pondir refdata/hg19/dbs/gridss/pon3792v1/ \
  --input all_calls.vcf \
  --output high_confidence_somatic.vcf.gz \
  --fulloutput high_and_low_confidence_somatic.vcf.gz \
  --scriptdir $(dirname $(which gridss_somatic_filter)) \
  -n 1 \
  -t 2

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call gridss-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))