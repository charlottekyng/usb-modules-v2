include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

# Usage: make TYPE=snvs or make TYPE=indels
TYPE ?= snvs
CALLERS ?= strelka2_$(TYPE) mutect2
CALLERS := $(strip $(CALLERS))
CALLER_STRING := $(subst $(space),.,$(CALLERS))
MIN_CALLERS ?= 2

LOGDIR ?= log/consensus_$(TYPE).$(NOW)

PHONY += all consensus_tables

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

VCF_LONGEST_SUFFIX = $(shell ls vcf/*.vcf |grep -v consensus| awk '{ print length, $$0 }' | sort -nr | head -1 | cut -d. -f3- | sed s/.vcf//)

# Main target to create consensus VCF
all: consensus_tables $(foreach pair,$(SAMPLE_PAIRS),vcf/$(tumor.$(pair))_$(normal.$(pair)).consensus.$(CALLER_STRING).$(VCF_LONGEST_SUFFIX).vcf)

consensus_tables : $(call MAKE_TABLE_FILE_LIST,consensus.$(CALLER_STRING))


# Preprocess step to normalize and annotate each caller's VCF
define preprocess

vcf/$1_$2.$3.norm.tagged.vcf.gz: vcf/$1_$2.$3.$(VCF_LONGEST_SUFFIX).vcf

	$$(INIT) module load $$(BCFTOOLS_MODULE); \
	bcftools norm -f $(REF_FASTA) -m -both $$< -Oz -o $$@
	bcftools index $$@
endef
# Loop through each sample pair and caller to preprocess each VCF
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach caller,$(CALLERS), \
		$(eval $(call preprocess,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
	) \
)

# Generate consensus VCF for each sample pair and all callers with bcftools isec
define consensus
vcf/$1_$2.consensus.$(CALLER_STRING).$(VCF_LONGEST_SUFFIX).vcf: $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).norm.tagged.vcf.gz)
	$$(INIT) module load $$(BCFTOOLS_MODULE); \
	bcftools isec -n+$(MIN_CALLERS) -w1 $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).norm.tagged.vcf.gz) -o $$@
	rm -f $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).norm.tagged.vcf.gz)
	rm -f $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).norm.tagged.vcf.gz.csi)
endef


# Loop through each sample pair and generate consensus VCFs
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call consensus,$(tumor.$(pair)),$(normal.$(pair)))) \
)

include usb-modules-v2/vcf_tools/vcftools.mk