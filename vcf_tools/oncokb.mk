include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

# Usage: consensus or single caller

TUMOR_TYPE ?=
TYPE ?= snvs
CONSENSUS ?= consensus
CALLERS ?= strelka2_$(TYPE) mutect2
CALLERS := $(strip $(CALLERS))
CALLER_STRING := $(subst $(space),.,$(CALLERS))
# Generate the suffix
CONSENSUS_VCFS := $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING))
FIRST_VCF := $(firstword $(CONSENSUS_VCFS))
VCF_SUFFIX := $(shell echo $(FIRST_VCF) |cut -d. -f$$(($(words $(CALLERS))+3))- | sed s/.vcf//)

LOGDIR ?= log/oncokb.$(NOW)

PHONY += all

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

all: alltables/allTN.consensus.$(CALLER_STRING).$(VCF_SUFFIX).$(TUMOR_TYPE).oncokb.maf

#oncokb_vcfs : $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING)) 
#oncokb_tables : $(call MAKE_TABLE_FILE_LIST,consensus.$(CALLER_STRING))

# convert VCF to MAF
define convert
maf/$1_$2.consensus.$(CALLER_STRING).$(VCF_SUFFIX).maf : vcf/$1_$2.consensus.$(CALLER_STRING).$(VCF_SUFFIX).vcf
endef
# Loop through each sample pair and caller to preprocess each VCF
$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call convert,$(tumor.$(pair)),$(normal.$(pair)))) \
)

# oncokb MAF annotator
define oncokb_annotator
maf/$1_$2.consensus.$(CALLER_STRING).$(VCF_SUFFIX).$(TUMOR_TYPE).oncokb.maf : maf/$1_$2.consensus.$(CALLER_STRING).$(VCF_SUFFIX).maf
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE)," 	sleep 5 && python $$(ONCOKB_MAF_ANNO) -i $$< -o $$@ $$(if $$(findstring hg38,$$(REF)), -r GRCh38) -t $$(TUMOR_TYPE) -b $$(ONCOKB_TOKEN) -d && $$(RM) $$<")
endef
# Loop through each sample pair and generate consensus VCFs
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call oncokb_annotator,$(tumor.$(pair)),$(normal.$(pair)))) \
)

# concatenate all MAF in a alltable
alltables/allTN.consensus.$(CALLER_STRING).$(VCF_SUFFIX).$(TUMOR_TYPE).oncokb.maf  : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).consensus.$(CALLER_STRING).$(VCF_SUFFIX).$(TUMOR_TYPE).oncokb.maf)

include usb-modules-v2/vcf_tools/vcftools.mk
