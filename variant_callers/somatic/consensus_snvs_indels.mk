include usb-modules-v2/Makefile.inc

# Usage: make TYPE=snvs or make TYPE=indels
TYPE ?= snvs
CALLERS ?= strelka2_$(TYPE) mutect2
CALLERS := $(strip $(CALLERS))
CALLER_STRING := $(subst $(space),.,$(CALLERS))
MIN_CALLERS ?= 2

LOGDIR ?= log/consensus_$(TYPE).$(NOW)

PHONY += all

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

$(info TYPE = $(TYPE))
$(info CALLERS = $(CALLERS))
$(info CALLER_STRING = $(CALLER_STRING))
$(info SAMPLE_PAIRS = $(SAMPLE_PAIRS))

# Main target to create consensus VCF
all: $(foreach pair,$(SAMPLE_PAIRS), \
	vcf/$(tumor.$(pair))_$(normal.$(pair)).consensus_$(TYPE).$(CALLER_STRING).som_ad_ft.nft.hotspot.pass.vcf)


# Preprocess step to normalize and annotate each caller's VCF
define preprocess
vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.norm.tagged.vcf.gz: vcf/$1_$2.$3.som_ad_ft.nft.hotspot.pass.vcf
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
vcf/$1_$2.consensus_$(TYPE).$(CALLER_STRING).som_ad_ft.nft.hotspot.pass.vcf: $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).som_ad_ft.nft.hotspot.pass.norm.tagged.vcf.gz)
	$$(INIT) module load $$(BCFTOOLS_MODULE); \
	bcftools isec -n+$(MIN_CALLERS) -w1 $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).som_ad_ft.nft.hotspot.pass.norm.tagged.vcf.gz) -o $$@
	rm -f $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).som_ad_ft.nft.hotspot.pass.norm.tagged.vcf.gz)
	rm -f $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).som_ad_ft.nft.hotspot.pass.norm.tagged.vcf.gz.csi)
endef


# Loop through each sample pair and generate consensus VCFs
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call consensus,$(tumor.$(pair)),$(normal.$(pair)))) \
)