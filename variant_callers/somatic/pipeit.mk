MUT_CALLER := pipeit

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

FILTER_VARIANTS = false

LOGDIR ?= log/pipeit.$(NOW)

PHONY += pipeit
pipeit : pipeit_vcfs pipeit_tables

VARIANT_TYPES ?= pipeit_snps pipeit_indels
pipeit_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,$(type))))
pipeit_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define pipeit-vcf
vcf/$1_$2.pipeit.vcf : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call RUN,4,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),,"\
	$$(PIPEIT) -t $$< -n $$(word 3,$$^) -o $$(@D) -e $(TARGETS_FILE_INTERVALS) -u $$(PRIMER_TRIM_BED) \
	-o $$@ -s $$(MIN_TUMOR_AD) -r $$(MIN_NORMAL_DEPTH) -m $$(MIN_TUMOR_DEPTH) -f $$(MIN_TN_AD_RATIO) \
	-j $$(TVC_SOMATIC_JSON) -a true")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))
