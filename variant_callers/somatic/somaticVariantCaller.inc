ifndef SOMATIC_VARIANT_CALLER_INC



### Filters

ifeq ($(findstring true,$(FILTER_VARIANTS)),true)
VCF_FILTERS = $(if $(TARGETS_FILE_INTERVALS),target_ft) som_ad_ft
ifeq ($(shell test $(words $(NORMAL_SAMPLES)) -gt 2; echo $$?),0)
	VCF_FILTERS += $(if $(findstring true,$(USE_NFT)),nft)
endif
VCF_FILTERS += $(if $(findstring GRCm38,$(REF)),,hotspot) pass
else
VCF_FILTERS =
endif
VCF_FILTER_SUFFIX = $(subst $( ),.,$(strip $(VCF_FILTERS)))


### Annotation
ifeq ($(findstring true,$(ANNOTATE_VARIANTS)),true)

VCF_ANNS = $(if $(findstring true,$(USE_SUFAM)),sufam) \
	$(if $(findstring GRCm38,$(REF)),eff dbsnp,\
	eff dbsnp cosmic exac_nontcga clinvar gene_ann nsfp) \
	$(if $(findstring true,$(ANN_NFT)),nft)

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
	VCF_ANNS += amplicons
endif

ifeq ($(findstring true,$(ANN_FACETS)),true)
	VCF_ANNS += facets
endif
else 
VCF_ANNS =
endif
VCF_ANNS_SUFFIX = $(subst $( ),.,$(strip $(VCF_ANNS)))

### 


VCF_SUFFIX = $(subst $( ),.,$(strip $1 $(VCF_FILTER_SUFFIX) $(VCF_ANNS_SUFFIX)))

VCF_SUFFIXES = $(foreach type,$1,$(call VCF_SUFFIX,$(type)))

MAKE_VCF_FILE_LIST = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call VCF_SUFFIXES,$1),\
	vcf/$(pair).$(suff).vcf))

ifeq ($(findstring false,$(INCLUDE_LNCRNA_IN_SUMMARY)),false)
TABLE_EFF_TYPES ?= nonsynonymous_synonymous_hotspot
else
TABLE_EFF_TYPES ?= nonsynonymous_synonymous_hotspot_lincRNA_upstream
endif

TABLE_EXTENSIONS = tab $(foreach eff,$(TABLE_EFF_TYPES),tab.$(eff))
TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
	 $(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

MAKE_TABLE_FILE_LIST = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call TABLE_SUFFIXES,$1),\
	tables/$(pair).$(suff).txt))
MAKE_TABLE_FILE_LIST += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/allTN$(PROJECT_PREFIX).$(suff).txt)

DOWNSTREAM_EFF_TYPES ?= $(word 1,$(TABLE_EFF_TYPES))
DOWMSTREAM_VCF_TABLE_SUFFIX ?= $(call VCF_SUFFIXES,$1).tab.$(DOWNSTREAM_EFF_TYPES)

SUFAM_ANN_SUFFIX = $(subst $( ),.,$(strip $(if $(findstring GRCm38,$(REF)),,hotspot) $(call SOMATIC_ANNS,$1)))
SUFAM_SUFFIX = $(subst $( ),.,$(strip $1 $(SUFAM_ANN_SUFFIX)))
SUFAM_VCF_SUFFIXES = $(foreach type,$1,$(call SUFAM_SUFFIX,$(type)))
SUFAM_TABLE_SUFFIXES = $(foreach suff,$(call SUFAM_VCF_SUFFIXES,$1),\
	 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))

endif
SOMATIC_VARIANT_CALLER_INC = true
