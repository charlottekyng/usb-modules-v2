ifndef SOMATIC_VARIANT_CALLER_INC



### Filters

SOMATIC_FILTERS = $(if $(TARGETS_FILE_INTERVALS),target_ft) som_ad_ft
ifeq ($(shell test $(words $(NORMAL_SAMPLES)) -gt 2; echo $$?),0)
	SOMATIC_FILTERS += $(if $(findstring true,$(USE_NFT)),nft)
endif
SOMATIC_FILTERS += $(if $(findstring GRCm38,$(REF)),,hotspot) pass

#SOMATIC_FILTER_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_FILTERS,$1)) pass)
SOMATIC_FILTER_SUFFIX = $(subst $( ),.,$(strip $(SOMATIC_FILTERS)))


### Annotation
SOMATIC_ANNS = $(if $(findstring true,$(USE_SUFAM)),sufam) \
	$(if $(findstring GRCm38,$(REF)),eff dbsnp,\
	eff dbsnp cosmic exac_nontcga clinvar gene_ann nsfp)

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
	SOMATIC_ANNS += amplicons
endif

ifeq ($(findstring true,$(ANN_FACETS)),true)
	SOMATIC_ANNS += facets
endif

SOMATIC_ANN_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_ANNS,$1)))



SOMATIC_SUFFIX = $(subst $( ),.,$(strip $1 $(SOMATIC_FILTER_SUFFIX) $(SOMATIC_ANN_SUFFIX)))

SOMATIC_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATIC_SUFFIX,$(type)))

SOMATIC_VCFS = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1),\
	vcf/$(pair).$(suff).vcf))


SOMATIC_EFF_TYPES ?= nonsynonymous_synonymous_hotspot nonsynonymous_synonymous_hotspot_lincRNA_upstream
SOMATIC_TABLE_EXTENSIONS = tab $(foreach eff,$(SOMATIC_EFF_TYPES),tab.$(eff))
SOMATIC_TABLE_SUFFIXES = $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1),\
	 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))

SOMATIC_TABLES = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),\
	tables/$(pair).$(suff).txt))
SOMATIC_TABLES += $(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),alltables/allTN$(PROJECT_PREFIX).$(suff).txt)


SUFAM_ANN_SUFFIX = $(subst $( ),.,$(strip $(if $(findstring GRCm38,$(REF)),,hotspot) $(call SOMATIC_ANNS,$1)))
SUFAM_SUFFIX = $(subst $( ),.,$(strip $1 $(SUFAM_ANN_SUFFIX)))
SUFAM_VCF_SUFFIXES = $(foreach type,$1,$(call SUFAM_SUFFIX,$(type)))
SUFAM_TABLE_SUFFIXES = $(foreach suff,$(call SUFAM_VCF_SUFFIXES,$1),\
	 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))

endif
SOMATIC_VARIANT_CALLER_INC = true
