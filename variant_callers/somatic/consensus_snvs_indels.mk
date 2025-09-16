include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

# Usage: make TYPE=snvs or make TYPE=indels
TYPE ?= snvs
CALLERS ?= strelka2_$(TYPE) mutect2
CALLERS := $(strip $(CALLERS))
CALLER_STRING := $(subst $(space),.,$(CALLERS))
MIN_CALLERS ?= 2

LOGDIR ?= log/consensus_$(TYPE).$(NOW)

PHONY += all clean

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

#clean intermediate files
all: consensus_tables consensus_vcfs
	rm -f vcf/*.norm2.vcf.gz vcf/*.norm2.vcf.gz.csi

consensus_vcfs : $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING)) 
consensus_tables : $(call MAKE_TABLE_FILE_LIST,consensus.$(CALLER_STRING))

# Generate the suffix
CONSENSUS_VCFS := $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING))
FIRST_VCF := $(firstword $(CONSENSUS_VCFS))
VCF_SUFFIX := $(shell echo $(FIRST_VCF) |cut -d. -f$$(($(words $(CALLERS))+3))- | sed s/.vcf//)

# Preprocess step to normalize and annotate each caller's VCF
define preprocess
vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.norm2.vcf.gz : vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.vcf
vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.norm2.vcf.gz.csi : vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.norm2.vcf.gz
endef
# Loop through each sample pair and caller to preprocess each VCF
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach caller,$(CALLERS), \
		$(eval $(call preprocess,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
	) \
)

# Generate consensus VCF for each sample pair and all callers with bcftools isec
define consensus
vcf/$1_$2.consensus.$(CALLER_STRING).vcf : $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).som_ad_ft.nft.hotspot.pass.norm2.vcf.gz.csi)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -n+$$(MIN_CALLERS) -w1 $$(basename $$^) -o $$@")

vcf/$1_$2.consensus.$(CALLER_STRING).$(VCF_SUFFIX).vcf : vcf/$1_$2.consensus.$(CALLER_STRING).vcf
endef
# Loop through each sample pair and generate consensus VCFs
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call consensus,$(tumor.$(pair)),$(normal.$(pair)))) \
)

include usb-modules-v2/vcf_tools/vcftools.mk
