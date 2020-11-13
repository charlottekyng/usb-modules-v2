MUT_CALLER := tvc
FILTER_VARIANTS := false

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pipeit.$(NOW)

PHONY += pipeit
pipeit : pipeit_vcfs pipeit_tables

VARIANT_TYPES ?= pipeit
pipeit_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,$(type))))
pipeit_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define pipeit-vcf
vcf/$1_$2.pipeit.vcf : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call RUN,4,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	$$(SINGULARITY_RUN) -B $$(dir $$(TARGETS_FILE_INTERVALS)) -B $$(dir $$(PRIMER_TRIM_BED))  \
	$$(PIPEIT_IMG) -t ./$$< -n ./$$(word 3,$$^) -e $$(TARGETS_FILE_INTERVALS) \
	-o $1_$2 -s $$(MIN_TUMOR_AD) -r $$(MIN_NORMAL_DEPTH) -m $$(MIN_TUMOR_DEPTH) -f $$(MIN_TN_AD_RATIO) \
	$$(ifdef $$(PIPEIT_JSON),-j ./$$(PIPEIT_JSON),,) -a true && ln PipeIT/results/$1_$2/$1_$2.PipeIT.vcf $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call pipeit-vcf,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
